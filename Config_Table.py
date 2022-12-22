import numpy as np
import Config_with_SNSPDs_and_QuadRF as Config
import Config_with_SNSPDs_and_QuadRF_Sprint as Config_Sprint # opx configuration for sprint experiments

# ------------- AOM Double-Pass calibration ---------------------------------------
from scipy.interpolate import griddata
c = np.load('UtilityResources\\calibrationFunction_griddata_cubic.npz', allow_pickle = True)['data']
c = dict(enumerate(c.flatten(), 1))[1] # Converting load data to dictionary
foo = lambda f,op : float(griddata(c['calibration_data'][0], c['calibration_data'][1],([f],[op]), method='cubic'))

calibrationData = {'ch_1':
                       {'calibration_func': foo,'affected_phases_names':['Pulse_1'], 'calibrate_RF_Power': c['calibrate_RF_Power'], 'calibrate_Optical_Power': 0.2},
                    'ch_2':
                       {'calibration_func': foo,'affected_phases_names':['Pulse_1'], 'calibrate_RF_Power': c['calibrate_RF_Power'], 'calibrate_Optical_Power': 0.2}
                   }

# calibrationData['ch_1']['calibration_func'] = lambda f, op: 30 + 0.82 * ((f*1e-6-113)*(f*1e-6-92.3)/(100))**2 if f <= 113e6 else 31.08
calibrationData = None
# Calibration data is created by placing power meter in front of the back exit of TOP1 (after double pass!) and running CalibrateLockAOM.py
# ---------------------------------------------------

Initial_Values = {
    # 'Operation_Mode': 'Magnetic_fountain',
    # 'Operation_Mode': 'Imaging',
    # 'Operation_Mode': 'PrePGC_Fountain',
    # 'Operation_Mode': 'OD_FS',
    # 'Operation_Mode': 'Depump',
    # 'Operation_Mode': 'Transit_Exp',
    # 'Operation_Mode': 'Spectrum_Exp',
    # 'Operation_Mode': 'CRUS_Exp',
    'Operation_Mode': 'SPRINT_Exp',
    # 'Operation_Mode': 'Continuous',
    'Imaging_Phase': 'Pulse_1',
    'Triggering_Phase': -1,  # Don't change this. Triggering phase should be defined within each operation mode (see below)
    # 'Operation_Mode': 'Double_PGC_with_Microwave',
    # 'Operation_Mode': 'Double_PGC_with_OD',
    'OD_Free_Space': False,

    # Operation frequencies:
    'AOM_TOP_2_freq': 100e6,
    'Flash_freq': 121.6625e6,
    'AOM_Repump_freq': 78.4735e6,
    'AOM_Depump_freq': 133.325e6,
    'Repump_PGC_freq': 78.4735e6,  # The final frequency of the Repump in the PGC sequence

    # MOT parameters:
    'MOT_freq': 113e6,
    'MOT_AOM_freq': 110e6,
    'MOT_duration': 2000,          # [msec]
    # 'MOT_duration': 300,          # [msec]
    'MOT_rep': -1,
    'AntiHelmholtz_delay': 0.1,    # [msec]
    'Post_MOT_delay': 1,  # [msec] All lights are off (except depump) between MOT and PGC prep stages

    # PGC parameters:
    'PGC_final_freq': 93e6, #98e6 - until 13.11,
    'PGC_initial_freq': -1,  # By default, taken from previous phase (see @Operation_Modes['DefaulValues'])
    'PGC_duration': 5,# 4 until 30.10.22
    # 'PGC_prep_duration': 4,        # [msec]
    'PGC_prep_duration': 5,  # 4 until 30.10.22        # [msec]
    'PGC_final_amp': 0.1, # 0.244 - till 13.11,         # Relative AOM amplitude between 0 to 1 - (0.12 yields 3.5 mW)
    'PGC_final_amp_0': 1,         # Relative AOM amplitude between 0 to 1 - (0.12 yields 3.5 mW)
    'PGC_final_amp_minus': 1,         # Relative AOM amplitude between 0 to 1 - (0.12 yields 3.5 mW)
    'PGC_final_amp_plus': 1,         # Relative AOM amplitude between 0 to 1 - (0.26 yields 3.6 mW)
    'PGC_initial_amp': -1,         # By default, taken from previous phase (see @Operation_Modes['DefaulValues'])
    'PGC_initial_amp_0': -1,         # By default, taken from previous phase (see @Operation_Modes['DefaulValues'])
    'PGC_initial_amp_minus': -1,         # By default, taken from previous phase (see @Operation_Modes['DefaulValues'])
    'PGC_initial_amp_plus': -1,        # By default, taken from previous phase (see @Operation_Modes['DefaulValues'])
    'PGC_initial_Delta_freq': 0,
    'PGC_final_Delta_freq': 0,

    # Pre PGC Fountain parameters:
    'Fountain_initial_freq': -1,  # By default, should be taken from previous phase
    'Fountain_final_freq': -1,  # By default, should be taken from next (?) phase
    'Fountain_duration': 0,        # [msec]
    'Fountain_prep_duration': 1,   # [msec], Can't be zero!!!
    'Fountain_initial_amp_0': 1,  # Relative amplitude between 0 to 1;
    'Fountain_initial_amp_minus': 1,  # Relative amplitude between 0 to 1;
    'Fountain_initial_amp_plus': 1,
    'Fountain_final_amp_0': 1,        # Relative AOM amplitude between 0 to 1
    'Fountain_final_amp_minus': 1,        # Relative AOM amplitude between 0 to 1
    'Fountain_final_amp_plus': 1,     # Relative AOM amplitude between 0 to 1
    'Fountain_final_Delta_freq': 0.35e6,  # 0.38e6 - until 30.10.22
    'Fountain_initial_Delta_freq': -1,         # By default, should be taken from previous phase

    # Free-Fall parameters:
    # The time take the atoms to reach the Toroids and fall back = 2 * sqrt(2 * (Toroid height[m] - MOT height[m]) / g[m/sec^2]) * 1e3[msec]
    # 'FreeFall_duration': 1.4 * np.sqrt(2 * 0.01 / 9.8) * 1e3,  # We use 1.4 because there is a limit for the pulse duration = 63ms
    'FreeFall_duration': 90,  # We use 1.4 because there is a limit for the pulse duration = 63ms
    'Coil_timing': 30,             # Zeeman coils turn on time from start of free fall[msec]
    'AOM_Off_Detuning': 80e6,      # [Hz]; Frequency for essentially turnning light coming from AONs off

    # Imaging parameters:
#   'Snap_Time': 1,              # [msec]
    'Trigger_delay': -0.1,         # [msec]
    'N_Snaps': 1,
    'Buffer_Cycles': -1, # By default, any change to N_Snaps will set this to 3. Don't change this value here

    # Depump measure:
    'Depump_pulse_duration': 0,  # [msec]
    'Depump_pulses_spacing': 0.2,  # [msec]
    'Depump_Start': 0.2,  # [msec]

    # OD free-space measure:
    'OD_FS_pulse_duration': 0,  # [msec]
    'OD_FS_pulses_spacing': 0.2,  # [msec]
    'OD_FS_Start': 0,  # [msec]
    'OD_FS_sleep': 0.2,  # [msec]

    # Transits / In-fiber OD:
    'M_delay': 35,  # [msec]
    'OD_delay': 0,  # [msec]
    'M_window': int(Config.readout_pulse_len),  # [nsec]
    'OD_duration_pulse1': 0,  # [msec]
    'OD_sleep': 2,  # [msec]
    'OD_duration_pulse2': 0,  # [msec]
    'M_time': 20,  # [msec]
    'M_off_time': 0,  # [msec]

    # Measuring pulses parameters:
    'PrePulse_duration': 1,        # [msec]
    'Shutter_open_time': 0,  # [msec]

    'Pulse_1_amp_i': 1,              # relative amplitude 0 to 1 (change in db is calculated by script)
    'Pulse_1_amp_f': 1,              # relative amplitude 0 to 1 (change in db is calculated by script)
    'Pulse_1_CH1_Freq_i': -1,        # [Hz]
    'Pulse_1_CH1_Freq_f': -1,        # [Hz]
    'Pulse_1_CH_2_3_Freq': -1,     # [Hz], usually should be equal to MOT_AOM_freq
    'Pulse_1_CH4_Freq': -1,         # [Hz], usually should be equal to Repump_PGC_freq

    'Pulse_1_duration': 0.2,       # [msec]
    'Pulse_1_decay_duration': 0,   # [msec] # by default, pulse 1 is sharp (no decay)
    'Pulse_1_final_amp_0': 1,         # Relative AOM amplitude between 0 to 1 - (0.12 yields 3.5 mW)
    'Pulse_1_final_amp_minus': 1,         # Relative AOM amplitude between 0 to 1 - (0.12 yields 3.5 mW)
    'Pulse_1_final_amp_plus': 1,         # Relative AOM amplitude between 0 to 1 - (0.26 yields 3.6 mW)
    'Pulse_1_initial_amp': 0,         # By default, taken from previous phase (see @Operation_Modes['DefaulValues'])
    'Pulse_1_initial_amp_0': 0,         # By default, taken from previous phase (see @Operation_Modes['DefaulValues'])
    'Pulse_1_initial_amp_minus': 0,         # By default, taken from previous phase (see @Operation_Modes['DefaulValues'])
    'Pulse_1_initial_amp_plus': 0,        # By default, taken from previous phase (see @Operation_Modes['DefaulValues'])
    'Pulse_1_Repump_amp': 1,  # Relative AOM amplitude between 0 to 1

    'InterPulses_duration': 0,     # [msec]; Time between pulse 1 & 2.

    'Pulse_2_CH1_Freq': -1,
    'Pulse_2_duration': 0,        # [msec]; If 0 -> no 2nd pulse

    'PostPulses_freq': -1,

    'opticalPowerCalibrationFunction': calibrationData # Used for calibration of the power of the lock AOM (TOP1)
}
## @Phases_Names: names of different phases. Note: order is important, as the index of this array defines self.triggerOnPhaseIndex in quadRFMOTController
# Phases_Names = ['MOT', 'PGC', 'Fountain', 'Free_Fall', 'Pulse_1', 'Inter_Pulses', 'Pulse_2', 'Post_Pulse']
Phases_Names = ['MOT', 'Fountain', 'PGC', 'Free_Fall', 'Pulse_1', 'Inter_Pulses', 'Pulse_2', 'Post_Pulse']

## @Operation_Modes defines the differences between different opeartion modes, compared to Initial_Values.
## e.g., 'Video' is exactly the same as @Initial_Values, except for N_Snaps, Null_Cycles etc.

Operation_Modes = {
                    'Default_Values': {
                                        'Triggering_Phase': 'Pulse_1',
                                        'Pulse_1_CH1_Freq_i': Initial_Values['MOT_freq'],
                                        # 'Pulse_1_CH1_Freq_i': Initial_Values['Flash_freq'],
                                        'Pulse_1_CH1_Freq_f': Initial_Values['MOT_freq'],
                                        # 'Pulse_1_CH1_Freq_f': Initial_Values['Flash_freq'],
                                        'Pulse_1_CH_2_3_Freq': Initial_Values['MOT_AOM_freq'],   # [Hz], usually should be equal to MOT_AOM_freq
                                        'Pulse_1_CH4_Freq': Initial_Values['Repump_PGC_freq'],  # [Hz], usually should be equal to Repump_PGC_freq
                                        'Pulse_2_CH1_Freq': Initial_Values['MOT_freq'],
                                        'Pulse_2_CH_2_3_Freq': Initial_Values['MOT_AOM_freq'],
                                        'Pulse_2_CH4_Freq': Initial_Values['AOM_Repump_freq'],
                                        # Fountain
                                        'Fountain_initial_Delta_freq': 0,         # By default, should be taken from previous phase
                                        'Fountain_initial_freq': Initial_Values['MOT_freq'],  # By default, should be taken from previous phase
                                        'Fountain_final_freq': Initial_Values['MOT_freq'],  # By default, should be taken from next (?) phase
                                        'PGC_initial_amp': 1,         # By default, should be taken from previous phase
                                        'PGC_initial_amp_0': 1,         # By default, should be taken from previous phase
                                        'PGC_initial_amp_minus': 1,         # By default, should be taken from previous phase
                                        'PGC_initial_amp_plus': 1,   # By default, should be taken from previous phase
                                        'PGC_initial_freq': Initial_Values['MOT_freq'],                # By default, should be taken from previous phase
                                        'PostPulses_freq': Initial_Values['MOT_freq'],
                                      },
                    'Continuous': {  # Note, this mode is different than others. It is, well, continuous.
                                    #'Triggering_Phase': 'Pulse_1',
                                    'CH1_freq': '113MHz',
                                    'CH2_freq': '100MHz',
                                    'CH3_freq': '110MHz',
                                    'CH4_freq': '78.4735MHz', #repump
                                    # 'CH4_freq': '133.325MHz', - depump

                                    # In order to turn a channel off, set amp to '0x0'
                                    'CH1_amp': '15.25dbm',
                                    'CH2_amp': '2dbm',
                                    'CH3_amp': '1.5dbm',
                                    'CH4_amp': '5.75dbm'
                                 },
                    'Imaging': {'Triggering_Phase': 'Pulse_1',
                                'Pulse_1_CH1_Freq_f': Initial_Values['MOT_freq'],
                                'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'],
                                'N_Snaps': 1,
                                'Buffer_Cycles': 0,
                                'Imaging_Phase': 'Pulse_1'
                                },
                    'OD_FS': {'Triggering_Phase': 'Pulse_1',
                              'Pulse_1_CH1_Freq_f': Initial_Values['MOT_freq'],
                              'Pulse_1_CH4_Freq': 190e6, # Means Rempump is off (on = 78e6)
                              'Pulse_1_Repump_amp': 0.0000001,
                              'N_Snaps': 1,
                              'Buffer_Cycles': 0,
                              'Imaging_Phase': 'Pulse_1',
                              'OD_FS_pulse_duration': 0.2,  # [msec]
                              'Pulse_1_duration': Initial_Values['OD_FS_Start'] + Initial_Values['OD_FS_pulses_spacing'] + 2 * 0.2 + 2 * Initial_Values['OD_FS_sleep'],  # [msec]
                              ## If with fountain:
                              # 'Pre_PGC_Fountain_duration': 1,  # [msec]
                              # 'PGC_duration': 5,
                              # 'PGC_prep_duration': 5,        # [msec]
                              # 'Fountain_duration': 0,  # [msec]
                              # 'Fountain_prep_duration': 1,  # [msec], Can't be zero!!!
                              # 'Fountain_initial_amp_0': 1,  # Relative amplitude between 0 to 1;
                              # 'Fountain_initial_amp_minus': 1,  # Relative amplitude between 0 to 1;
                              # 'Fountain_initial_amp_plus': 1,  # Relative amplitude between 0 to 1;
                              # 'Fountain_final_amp_0': 1,  # Relative amplitude between 0 to 1;
                              # 'Fountain_final_amp_minus': 1,  # Relative amplitude between 0 to 1;
                              # 'Fountain_final_amp_plus': 1,  # Relative amplitude between 0 to 1;
                              # 'OD_Free_Space': True, # ZA : what is this?
                              },
                    'Depump': {'Triggering_Phase': 'Pulse_1',
                               'Pulse_1_CH1_Freq_f': Initial_Values['MOT_freq'],
                               'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'] + 30e6, # Means off
                               'Pulse_1_Repump_amp': 1, # Turned off by detuning
                               'N_Snaps': 1,
                               'Buffer_Cycles': 0,
                               'Imaging_Phase': 'Pulse_1',
                               'Depump_pulse_duration': 0.2,  # [msec]
                               'Pulse_1_duration': Initial_Values['Depump_Start'] + Initial_Values['Depump_pulses_spacing'] + 2 * 0.2,  # [msec]
                               'Fountain_duration': 1,        # [msec]
                               'Fountain_prep_duration': 1,  # [msec], Can't be zero!!!
                               'Fountain_initial_amp_0': 1,  # Relative amplitude between 0 to 1;
                               'Fountain_initial_amp_minus': 1,  # Relative amplitude between 0 to 1;
                               'Fountain_initial_amp_plus': 1,  # Relative amplitude between 0 to 1;
                               'Fountain_final_amp_0': 1,  # Relative amplitude between 0 to 1;
                               'Fountain_final_amp_minus': 1,  # Relative amplitude between 0 to 1;
                               'Fountain_final_amp_plus': 1  # Relative amplitude between 0 to 1;
                               },
                    'Transit_Exp': {'Triggering_Phase': 'Free_Fall',
                                    'Pulse_1_CH1_Freq_f': Initial_Values['MOT_freq'],
                                    'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'] + 30e6,
                                    'N_Snaps': 1,
                                    'Buffer_Cycles': 0,
                                    'Imaging_Phase': 'Pulse_1',
                                    'PrePulse_duration': 10,  # [msec]
                                    'Shutter_open_time': 4,  # [msec]
                                    'Pulse_1_Repump_amp': 0.000001,
                                    'Pulse_1_duration': 60,  # [msec]
                                    'M_time': 60,  # [msec]
                                    'M_off_time': 1,  # [msec]
                                    ## If with fountain:
                                    'Fountain_duration': 0.5,  # [msec]
                                    'Fountain_prep_duration': 0.5,  # [msec], Can't be zero!!!
                                    'Fountain_final_Delta_freq': 0.45e6,  # 0.38e6 - until 30.10.22
                                    'PGC_duration': 5.1  # [msec] EXTREMELY IMOPRTANT for OPX-QuadRF sync
                                    },
                    'Spectrum_Exp': {'Triggering_Phase': 'Free_Fall',
                                    'Pulse_1_CH1_Freq_f': Initial_Values['MOT_freq'],
                                    'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'] + 30e6,
                                    'N_Snaps': 1,
                                    'Buffer_Cycles': 0,
                                    'Imaging_Phase': 'Pulse_1',
                                    'PrePulse_duration': 20,  # [msec]
                                    'Shutter_open_time': 3,  # [msec]
                                    'Pulse_1_Repump_amp': 0.000001,
                                    'Pulse_1_duration': len(Config.Spectrum_Exp_Gaussian_samples) * 1000 * 21 * 4 / 1e6,  # [msec]
                                    ## If with fountain:
                                    'Fountain_duration': 0.5,  # [msec]
                                    'Fountain_prep_duration': 0.5,  # [msec], Can't be zero!!!
                                    'M_window': int(Config.readout_pulse_spectrum_len), # [nsec]
                                    'M_time': len(Config.Spectrum_Exp_Gaussian_samples) * 1000 * 21 * 4 / 1e6,  # Pulse_length[nsec] * 1000 repetitions * (Bandwidth[MHz] * frequency steps[MHz]) * 4 / 1e6[nsec/msec] - [msec]
                                     },
                    'CRUS_Exp': {'Triggering_Phase': 'Free_Fall',
                                 'Pulse_1_CH1_Freq_f': Initial_Values['MOT_freq'],
                                 'Pulse_1_amp_f': 1, # decides the power distribution between on-res & detuned pulses
                                 'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'] + 30e6,
                                 'N_Snaps': 1,
                                 'Buffer_Cycles': 0,
                                 'Imaging_Phase': 'Pulse_1',
                                 'PrePulse_duration': 10,  # [msec]
                                 'Shutter_open_time': 3,  # [msec]
                                 'Pulse_1_Repump_amp': 0.000001,
                                 'Pulse_1_duration': (Config.readout_CRUS_pulse_len) / 1e6,  # [msec]
                                 ## If with fountain:
                                 'Fountain_duration': 0.5,  # [msec]
                                 'Fountain_prep_duration': 0.5,  # [msec], Can't be zero!!!
                                 'M_window': int(Config.readout_CRUS_pulse_len), # [nsec]
                                 # 'M_off_time': 10, # [msec]
                                 'M_off_time': 10, # [msec]
                                 'M_time': 3*(Config.readout_CRUS_pulse_len) / 1e6,  # [msec]
                                 'Fountain_final_Delta_freq': 0.45e6,  # 0.38e6 - until 30.10.22
                                 'PGC_duration': 5.1 #[msec] EXTREMELY IMOPRTANT for OPX-QuadRF sync
                                 },
                    'SPRINT_Exp': {'Triggering_Phase': 'Free_Fall',
                                   'Fountain_final_Delta_freq': 0.45e6,  # 0.38e6 - until 30.10.22
                                   'Pulse_1_CH1_Freq_f': Initial_Values['MOT_freq'],
                                    'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'] + 30e6,
                                    'N_Snaps': 1,
                                    'Buffer_Cycles': 0,
                                    'Imaging_Phase': 'Pulse_1',
                                    'PrePulse_duration': 10,  # [msec]
                                    'Shutter_open_time': 3,  # [msec]
                                    'Pulse_1_Repump_amp': 0.000001,
                                    'Pulse_1_duration': 3*int(Config_Sprint.readout_pulse_sprint_len_N)/ 1e6,  # [msec]
                                    ## If with fountain:
                                    'Fountain_duration': 0.5,  # [msec]
                                    'Fountain_prep_duration': 0.5,  # [msec], Can't be zero!!!
                                    'M_window': int(Config_Sprint.readout_pulse_sprint_len_N), # [nsec]
                                    'M_time': 3*int(Config_Sprint.readout_pulse_sprint_len_N) / 1e6,  # Pulse_length[nsec] * 1000 repetitions * (Bandwidth[MHz] * frequency steps[MHz]) * 4 / 1e6[nsec/msec] - [msec]
                                    'M_off_time': 5,  # [msec] - should be at least 5 ms, to sync quadrf and OPX
                                   },
                    'PrePGC_Fountain': {'Triggering_Phase': 'Pulse_1',
                                        'Fountain_final_Delta_freq': 0.45e6,
                                        'Pulse_1_CH1_Freq_f': Initial_Values['Flash_freq'],
                                        'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'],
                                        # Fountain
                                        'Fountain_duration': 0.5,        # [msec]
                                        'Fountain_prep_duration': 0.5,  # [msec], Can't be zero!!!
                                        # Imaging
                                        'Imaging_Phase': 'Pulse_1',
                                        'PrePulse_duration': 1,  # [msec]
                                        },
                    'Fountain': {
                        # 'PGC_duration': 20,
                        # 'PGC_prep_duration': 7,  # [msec]
                        # 'PGC_final_amp_0': 0.03,  # Relative AOM amplitude between 0 to 1 - (0.12 yields 3.5 mW)
                        # 'PGC_final_amp_plus': 0.06,  # Relative AOM amplitude between 0 to 1 - (0.26 yields 3.6 mW)
                        'Triggering_Phase': 'Pulse_1',  # This is the name of the phase on which trigger from OPX is sent; this is the first thing one should check if there's a mismatch between OPX & QuadRF
                        'Pulse_1_CH1_Freq_f': Initial_Values['MOT_freq'],
                        'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'],
                        # 'MOT_duration': 3000,  # [msec]
                        'N_Snaps': 1,
                        'Buffer_Cycles': 0,
                        'Fountain_duration': 1.2,  # [msec]
                        'Fountain_prep_duration': 1,  # [msec]
                        # 'Fountain_final_amp_0': 0.03,  # Relative AOM amplitude between 0 to 1
                        # 'Fountain_final_amp_plus': 0.06,  # Relative AOM amplitude between 0 to 1
                        # 'PGC_final_Delta_freq': 0,
                        # 'Fountain_Delta_freq': 0.35e6,
                        'Imaging_Phase': 'Pulse_1'
                    },
                    'Fountain_detuning': {
                        'Triggering_Phase': 'Pulse_1',  # This is the name of the phase on which trigger from OPX is sent; this is the first thing one should check if there's a mismatch between OPX & QuadRF
                        'Pulse_1_CH1_Freq_f': 118e6,
                        # 'Pulse_1_CH1_Freq_f': Initial_Values['MOT_freq'],
                        'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'],
                        'PGC_final_freq': 93e6,
                        'N_Snaps': 1,
                        'Fountain_duration': 1.2,  # [msec]
                        'Fountain_prep_duration': 1,  # [msec]
                        'Buffer_Cycles': 0,
                        'PGC_duration': 5,
                        'PGC_prep_duration': 5,  # [msec]
                        'PGC_final_amp_0': 0.3,  # Relative AOM amplitude between 0 to 1 - (0.12 yields 3.5 mW)
                        'PGC_final_amp_plus': 0.3,  # Relative AOM amplitude between 0 to 1 - (0.26 yields 3.6 mW)
                        'Fountain_initial_amp_0': 0.3,   # Relative amplitude between 0 to 1; By default, should be taken from previous phase
                        'Fountain_initial_amp_plus': 0.3,  # Relative amplitude between 0 to 1; By default, should be taken from previous phase
                        'PGC_final_Delta_freq': 0,
                        'Imaging_Phase': 'Pulse_1'
                    },
                    'Fountain_detuning_with_OD': {
                        'Triggering_Phase': 'Pulse_1',  # This is the name of the phase on which trigger from OPX is sent; this is the first thing one should check if there's a mismatch between OPX & QuadRF
                        'Pulse_1_CH1_Freq_f': 118e6,
                        # 'Pulse_1_CH1_Freq_f': Initial_Values['MOT_freq'],
                        'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'],
                        'OD_ON': True,
                        'PGC_final_freq': 93e6,
                        'N_Snaps': 1,
                        'Buffer_Cycles': 0,
                        'PGC_duration': 4.1,
                        'PGC_prep_duration': 4,  # [msec]
                        'PGC_final_amp_0': 1,  # Relative AOM amplitude between 0 to 1 - (0.12 yields 3.5 mW)
                        'PGC_final_amp_plus': 1,  # Relative AOM amplitude between 0 to 1 - (0.26 yields 3.6 mW)
                        'PGC_final_Delta_freq': 0.38e6,
                        'OD_duration': 0.2,
                        'Wait_duration': 2,
                        'Imaging_Phase': 'Pulse_1'
                    },
                    'Magnetic_fountain': {
                        'Triggering_Phase': 'Pulse_1',  # This is the name of the phase on which trigger from OPX is sent; this is the first thing one should check if there's a mismatch between OPX & QuadRF
                        # 'Pulse_1_CH1_Freq_f': Initial_Values['MOT_freq'],
                        'Pulse_1_CH1_Freq_f': 113e6,
                        'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'],
                        'PGC_final_freq': 100e6,
                        'N_Snaps': 1,
                        'Buffer_Cycles': 0,
                        'PGC_duration': 4.5,
                        'PGC_prep_duration': 2,  # [msec]
                        'PGC_final_amp_0': 0.02,  # Relative AOM amplitude between 0 to 1 - (0.12 yields 3.5 mW)
                        'PGC_final_amp_plus': 0.06,  # Relative AOM amplitude between 0 to 1 - (0.26 yields 3.6 mW)
                        'PGC_final_Delta_freq': 0,
                        'Zeeman_delay': 0.6,
                        'AntiHelmholtz_delay': 0.1,
                        'Post_MOT_delay': 0.6,
                        'Imaging_Phase': 'Pulse_1'
                    },
                    'Video': {'N_Snaps': 20,
                              'Buffer_Cycles': 0,
                              },
                    'Double_PGC': {
                                   'Triggering_Phase': 'Pulse_2',
                                   'Pulse_1_decay_duration': 2,  # [msec]
                                   'Pulse_1_duration': 3,  # [msec]
                                   'InterPulses_duration': 1,  # [msec]
                                   'Pulse_2_duration': 0.2,  # [msec]
                                   'Pulse_1_amp_i': 1,
                                   'Pulse_1_amp_f': 0.12,
                                   'PrePulse_duration': 10,
                                   'Imaging_Phase': 'Pulse_2'
                                   },
                    # 'Double_PGC_with_OD': {
                    #                       'Triggering_Phase': 'PGC',
                    #                       'Pulse_1_decay_duration': 2,  # [msec]
                    #                       'Pulse_1_duration': 3,  # [msec]
                    #                       'InterPulses_duration': 1,  # [msec]
                    #                       'Pulse_2_duration': 1.5,  # [msec]
                    #                       'Pulse_1_amp_i': 1,
                    #                       'Pulse_1_amp_f': 0.12,
                    #                       'Pulse_2_CH1_Freq': Initial_Values['Flash_freq_det'],
                    #                       'PrePulse_duration': 10,        # [msec] #
                    #                       },
                    # 'Double_PGC_with_Microwave': {
                    #                            # Pulses Duration
                    #                            'PrePulse_duration': 10,  # [msec] #
                    #                            'Pulse_1_decay_duration': 2,  # [msec]
                    #                            'Pulse_1_duration': 3,  # [msec]
                    #                            'InterPulses_duration': 0,  # [msec]
                    #                            'Pulse_2_duration': 1,  # [msec]
                    #                             # Pulses amplitudes
                    #                            'Pulse_1_amp_i': 1,
                    #                            'Pulse_1_amp_f': 0.12,
                    #                             # Pulses Frequency
                    #                            'Pulse_2_CH1_Freq': Initial_Values['Flash_freq_det'],
                    #                            'Pulse_2_CH_2_3_Freq':  Initial_Values['MOT_AOM_freq'] + Initial_Values['AOM_Off_Detuning'], # turn Ch. 2 & 3 off
                    #                            'Pulse_2_CH4_Freq': Initial_Values['AOM_Depump_freq']
                    #                        },
                    'Microwave': {
                                'Pulse_1_duration': 3,  # [msec]
                                'Pulse_1_CH_2_3_Freq': Initial_Values['MOT_AOM_freq'] + Initial_Values['AOM_Off_Detuning'], # turn Ch. 2 & 3 off
                                'Pulse_1_CH4_Freq': Initial_Values['AOM_Depump_freq'] # This is the depump pulse, so depump freq
                    },
                    # 'OD_FS': {
                    #         'Triggering_Phase': 'Pulse_1',
                    #         'OD_Free_Space': True,
                    #         'Pulse_1_duration': 0.2,  # [msec]
                    #         'Pulse_1_CH_2_3_Freq': Initial_Values['MOT_AOM_freq'] + Initial_Values['AOM_Off_Detuning'],  # turn Ch. 2 & 3 off
                    #         'Pulse_1_CH4_Freq': Initial_Values['AOM_Depump_freq'],  # This is the depump pulse, so depump freq
                    #         'Imaging_Phase': 'Pulse_1'
                    # },
                   }
def updateOperationMode(OM = None, values = None):
    values = Initial_Values if values is None else values
    if OM is None: OM = Initial_Values['Operation_Mode']
    if OM not in Operation_Modes:
        print('\033[91m' + f'Warning: {OM} is not a known operation mode.' + '\033[0m')
    for key in Operation_Modes[OM]:
        values[key] = Operation_Modes[OM][key]
    return(values)

Initial_Values = updateOperationMode(OM = 'Default_Values') # First, make sure default values are there
Initial_Values = updateOperationMode()                      # Then, update the values pertaining to the specific operation mode


## ------------------------------------------------------------------------------------------------------------------------  ##
##                      Frequency Scan Values.
## ------------------------------------------------------------------------------------------------------------------------  ##
calibrationData = {'ch_2':
                       {'calibration_func': foo, 'calibrate_RF_Power': c['calibrate_RF_Power'], 'calibrate_Optical_Power': 0.2}
                   }
Frequency_Scan_Config_Values = {
    'Start_Delay_time': -1, # [ms], Time to wait before scan begins
    'Saw_Tooth_Scan': True, # If TrueL jump. If false, reverse frequency direction (tiangle)
    'Frequency_Range': (93e6, 133e6),# Order matters. If lower freq is first, we scan up. Otherwise: scan down.
    'Start_Frequency': None, # if None, start at first frequency. Otherwise, start scan closest to @Start_Frequency as possible
    'Total_Scan_Duration': Initial_Values['Pulse_1_duration'], # [ms] This should be the same as the TOTAL length od the measurement
    'N_Different_Frequencies': 21, # This defines the resolution of the scan
    'N_Scan_Repetitions':4, # Number of time we go through all frequencies during @Total_Scan_Duration
    'opticalPowerCalibrationFunction': calibrationData  # Used for calibration of the power of the lock AOM (TOP1)
}
Frequency_Scan_Config_Values['Start_Delay_time'] = 0 if Initial_Values['Triggering_Phase'] == 'Pulse_1' else Initial_Values['PrePulse_duration']
## ------------------------------------------------------------------------------------------------------------------------  ##
##                      Other values, not operation modes related.
## ------------------------------------------------------------------------------------------------------------------------  ##

# @Values_Factor hold factoring for values - each key has an array:
# So, as Prep_duration is written in ms, we have to factor it into an int value, but in 4ns, for the OPX.
# These are only factored in function "updateValue" in OPX control
# Values_Factor = {
#     'MOT_duration': [int, 1e6 / Config.MOT_pulse_len],
#     'AntiHelmholtz_delay': [int, 1e6 / 4],  # [msec]
#     'PGC_duration': [int, 1e6 / 4],  # [msec]
#     'Post_MOT_delay': [int, 1e6 / 4],  # [msec]
#     'PGC_prep_duration': [int, 1e6 / 4],  # [msec]
#     'PGC_final_amp': [int, None, lambda x: np.log(int(1 / (1 - x)))],
#     'PGC_final_amp_0': [int, None, lambda x: np.log(int(1 / (1 - x)))],
#     'Fountain_duration': [int, 1e6 / 4],  # [msec]
#     'Fountain_prep_duration': [int, 1e6 / 4],  # [msec]
#
#     # The frequency in the frame of reference where the atoms have just enough velocity to reach the toroids [Hz].
#     # Free-Fall parameters:
#     # The time take the atoms to reach the Toroids and fall back = 2 * sqrt(2 * (Toroid height[m] - MOT height[m]) / g[m/sec^2]) * 1e3[msec]
#     'FreeFall_duration': [int, 1e6 / 4],
#     # We use 1.4 because there is a limit for the pulse duration. = 63ms
#     'Coil_timing': [int, 1e6 / 4],  # Zeeman coils turn on time from start of free fall[msec]
#     # Imaging parameters:
#     'Snapshot_Intervals': [int, 1e6 / 4],  # [msec]
#     'PrePulse_duration': [int, 1e6 / 4],  # [msec]
#     'Trigger_delay': [int, 1e6 / 4],  # [msec]
#     'Pulse_1_duration': [int, 1e6 / 4],  # [msec]
#     'Pulse_1_decay_duration': [int, 1e6 / 4],  # [msec]
#     'InterPulses_duration': [int, 1e6 / 4],
#     'Pulse_2_duration': [int, 1e6 / 4],  # [msec]
#     # OD parameters:
#     'OD_duration': [int, 1e6 / 4],  # [msec]
#     'Wait_duration': [int, 1e6 / 4],  # [msec]
#     'OD_Start': [int, 1e6 / 4],  # [msec]
# }


IOParametersMapping = {  # These are chans. in OPX, and should all be int(s). The values have no real significance, except for in Parameters update in OPX control
  "MOT_switch": 1,
  "Linear_PGC_switch": 2,
  "Transit_Exp_switch": 3,
  "Spectrum_Exp_switch": 4,
  "CRUS_Exp_switch": 4,
  "AntiHelmholtz_Delay_switch": 5,
  "Max_Probe_counts_switch": 6,
  "AntiHelmholtz_delay": 10,
  "MOT_duration": 11,
  "Post_MOT_delay": 12,
  "PGC_duration": 20,
  "PGC_prep_duration": 21,
  "PGC_initial_amp_0": 22,
  "PGC_initial_amp_minus": 23,
  "PGC_initial_amp_plus": 24,
#  "PGC_final_amp": -1,
  "PGC_final_amp_0": 25,
  "PGC_final_amp_minus": 26,
  "PGC_final_amp_plus": 27,
  "PGC_pulse_duration_0": 28,
  "PGC_pulse_duration_minus": 29,
  "PGC_pulse_duration_plus": 30,
  "Fountain_duration": 31,
  "Fountain_aom_chirp_rate": 32,
  "Trigger_delay": 41,
  "PrePulse_duration": 42,
  "Pulse_1_duration": 43,
  "Pulse_1_decay_duration": 44,
  "N_Snaps": 45,
  "Buffer_Cycles": 46,
  "OD_FS_Start": 51,
  "OD_FS_pulse_duration": 52,
  "OD_FS_pulses_spacing": 53,
  "OD_delay": 54,
  "OD_SNSPDs_duration": 55,
  "Depump_Start": 56,
  "Depump_pulse_duration": 57,
  "Depump_pulses_spacing": 58,
  "Shutter_opening_time": 59
}

# For each key, the channels affected by the change of value. This way, only channels which we actually change are updated on the QuadRF

Key_to_Channel = {
    'Operation_Mode': [],
    'Triggering_Phase':[1,4],
    'MOT_freq': [1],
    'PGC_final_freq': [1],
    'Fountain_initial_freq':[1],
    'Fountain_final_freq':[1],
    'Flash_freq': [1],
    'AOM_Repump_freq': [4],
    'MOT_AOM_freq': [2],
    'MOT_duration': [1,4],        # [msec]
    'MOT_rep': [1,4],
    'PGC_duration': [1,4],
    'Post_MOT_delay': [1,4],
    'PGC_prep_duration': [1,4],      # [msec]
    'PGC_final_amp':[1],
    'PGC_final_amp_0': [1],       # Relative AOM amplitude between 0 to 1
    'PGC_final_amp_plus': [1],       # Relative AOM amplitude between 0 to 1
    'Repump_PGC_freq': [4],     # The final frequency of the Repump in the PGC sequence
    'Fountain_duration': [1,4],    # [msec]
    'Fountain_prep_duration': [1,4],          # [msec]
    'Fountain_initial_amp_0': [],  # (also minus) Relative AOM amplitude between 0 to 1
    'Fountain_initial_amp_plus': [],
    'Fountain_initial_amp_minus': [],
    'Fountain_final_amp_0': [],      # (also minus) Relative AOM amplitude between 0 to 1
    'Fountain_final_amp_plus': [],   # Relative AOM amplitude between 0 to 1
    'PGC_final_Delta_freq': [],        # The frequency in the frame of reference where the atoms have just enough velocity to reach the toroids [Hz].
    'PGC_initial_Delta_freq': [],        # The frequency in the frame of reference where the atoms have just enough velocity to reach the toroids [Hz].
    'Fountain_Delta_freq': [],        # The frequency in the frame of reference where the atoms have just enough velocity to reach the toroids [Hz].
    'Fountain_final_Delta_freq': [],
    # The frequency in the frame of reference where the atoms have just enough velocity to reach the toroids [Hz].

    # Free-Fall parameters:
    # The time take the atoms to reach the Toroids and fall back = 2 * sqrt(2 * (Toroid height[m] - MOT height[m]) / g[m/sec^2]) * 1e3[msec]
    'FreeFall_duration': [1,4],  # We use 1.4 because there is a limit for the pulse duration. = 63ms
    'OD_FS_sleep':[1,4],
    # Imaging parameters:
    'AOM_Off_Detuning': [4],      # [Hz]; Frequency for essentially turnning light coming from AONs off
    'Snapshot_Intervals': [1,4],   # [msec]
    'PrePulse_duration': [1,4],            # [msec]
    'Trigger_delay': [1,4],          # [msec]
    'Pulse_1_duration': [1,4],       # [msec]
    'InterPulses_duration': [1, 4], # [msec]; Time between pulse 1 & 2.
    'Pulse_2_duration': [1, 4],      # [msec]; If 0 -> no snd pulse
    'Zeeman_delay': [1,4],          # [msec]
    'N_Snaps': [1,4],
    # 'Buffer_Cycles': [1,2,3,4]
    'Pulse_1_decay_duration': [1,  4],
    'Pulse_1_amp_i': [1, 4],
    'Pulse_1_amp_f': [1, 4],
    'Pulse_1_CH1_Freq_i': [1],
    'Pulse_1_CH1_Freq_f': [1],
    'Pulse_2_CH1_Freq': [1],
    'Pulse_2_CH_2_3_Freq': [2, 3],
    'Pulse_2_CH4_Freq': [4],
    'PostPulses_freq': [1],
    # Continuous params
    'CH1_freq': [1],
    'CH2_freq': [2],
    'CH3_freq': [3],
    'CH4_freq': [4],
    'CH1_amp': [1],
    'CH2_amp': [2],
    'CH3_amp': [3],
    'CH4_amp': [4]

}
# verifyKeysForCaseSensitivity goes over all pairs of keys an verifies there isn't a pair which is the same up-to cases.
# i.e. it will generate a warning if we have 'MOT_Duration' in one dictionary and 'MOT_duration' in another.
import itertools
allDictionaries = [Initial_Values, IOParametersMapping, Key_to_Channel]#, Values_Factor]

def verifyKeysForCaseSensitivity():
    keys = []
    for dic in allDictionaries:
        keys = keys + list(dic.keys())
    # Now keys is one list of all the keys in @allDictionaries
    for pair in itertools.combinations(keys, 2):
        if pair[0].lower() == pair[1].lower() and pair[0] != pair[1]:
            print('\033[91m' + f'Warning: {pair[0]} is not the same as {pair[1]}' + '\033[0m')

verifyKeysForCaseSensitivity()


def updateMOT_rep():
    Initial_Values['MOT_rep'] = int(np.ceil((Initial_Values['MOT_duration'] * 1e6) / Config.MOT_pulse_len))
    return(Initial_Values['MOT_rep'])

updateMOT_rep()

