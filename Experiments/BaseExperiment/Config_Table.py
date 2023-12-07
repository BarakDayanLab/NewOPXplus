from enum import Enum
import numpy as np
import Experiments.BaseExperiment.Config_Experiment as Config
import Experiments.SPRINT.Config_with_SNSPDs_and_QuadRF_Sprint as Config_Sprint  # opx configuration for sprint experiments
import Experiments.QRAM.Config_with_SNSPDs_and_QuadRF_QRAM as Config_QRAM  # opx configuration for sprint experiments
from UtilityResources.AOMCalibration import calibrate, calibration_data
from Utilities.Utils import Utils


Operation_Modes_Enumeration = [
    'Magnetic_fountain',
    'Imaging',
    'PrePGC_Fountain',
    'OD_FS',
    'Depump',
    'Teansit_Exp',
    'Spectrum_Exp',
    'CRUS_Exp',
    'SPRINT_2-3_Exp',
    'SPRINT_Exp',
    'QRAM_Exp',
    'Transits_Exp',
    'Continuous',
    'Double_PGC_with_Microwave',
    'Double_PGC_with_OD'
]


#-------------------------------------------
# Initial Values
#-------------------------------------------

# TODO 1: Merge between the Initial_Values and the Operation_Mode('default_values')
Initial_Values = {
    # TODO - eventually remove this! Now QuadRFMOTController uses Operation_Mode to check if it's "Continuous"
    #'Operation_Mode': 'PrePGC_Fountain',
    'Operation_Mode': 'dont care!',

    'Experiment_Name': 'Default',

    'Imaging_Phase': 'Pulse_1',
    'Triggering_Phase': -1,  # Don't change this. Triggering phase should be defined within each operation mode (see below)
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
    'PGC_final_freq': 91e6, #98e6 - until 13.11,
    'PGC_initial_freq': -1,  # By default, taken from previous phase (see @Operation_Modes['DefaulValues'])
    'PGC_duration': 2.5,# 4 until 30.10.22
    # 'PGC_prep_duration': 4,        # [msec]
    'PGC_prep_duration': 2.5,  # 4 until 30.10.22        # [msec]
    'PGC_final_amp': 0.05, # 0.244 - till 13.11,         # Relative AOM amplitude between 0 to 1 - (0.12 yields 3.5 mW)
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
    'FreeFall_duration': 84,  # We use 1.4 because there is a limit for the pulse duration = 63ms
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
    # TODO Q: Why are we using both "M_window" and "readout_pulse_len" ?
    'M_window': int(Config.readout_pulse_len),  # [nsec]
    'OD_duration_pulse1': 0,  # [msec]
    'OD_sleep': 2,  # [msec]
    'OD_duration_pulse2': 0,  # [msec]
    'M_time': 20,  # [msec]
    'M_off_time': 0,  # [msec]

    # Measuring pulses parameters:
    'PrePulse_duration': 1,  # [msec]
    'Shutter_open_time': 0,  # [msec]
    'PrePulse_Repump_amp': 1,  # relative

    'Pulse_1_amp_i': 1,              # relative amplitude 0 to 1 (change in db is calculated by script)
    'Pulse_1_amp_f': 1,              # relative amplitude 0 to 1 (change in db is calculated by script)
    'Pulse_1_CH1_Freq_i': -1,        # [Hz]
    'Pulse_1_CH1_Freq_f': -1,        # [Hz]
    'Pulse_1_CH_2_3_Freq': -1,     # [Hz], usually should be equal to MOT_AOM_freq
    'Pulse_1_CH4_Freq': -1,         # [Hz], usually should be equal to Repump_PGC_freq

    'Pulse_1_duration': 5,       # [msec]
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

    'OPX_Quad_Misalignment_Delay': 4000,  # = 4us [ns]

    'opticalPowerCalibrationFunction': None  # calibrationData  # Used for calibration of the power of the lock AOM (TOP1)
}

#-------------------------------------------
# Phase Names
#
# Names of different phases in the experiment.
# The order of the phases is important - the index of each phase defines self.triggerOnPhaseIndex in quadRFMOTController
# Index starts with 0
# TODO: Can we change the below into Enumeration? Need to also change the code in QuadRFController
#-------------------------------------------
Phases_Names = ['MOT', 'Fountain', 'PGC', 'Free_Fall', 'Pulse_1', 'Inter_Pulses', 'Pulse_2', 'Post_Pulse']
class Phases(Enum):
    MOT = 0
    FOUNTAIN = 1
    PGC = 2
    FREE_FALL = 3,
    PULSE_1 = 4,
    INTER_PULSES = 5,
    PULSE_2 = 6,
    POST_PULSE = 7
# Usage:
# >>> x = Phases.PGC
# >>> Phases.PGC.name
# >>> Phases.PGC.value

#-------------------------------------------
# Operation Modes
#
# Operation_Modes defines the differences between different operations, compared to Initial_Values.
# e.g., 'Video' is exactly the same as @Initial_Values, except for N_Snaps, Null_Cycles etc.
#-------------------------------------------

# TODO: eventually this goes away! As all parameters where put separately in relevan experiment.
Operation_Modes__DEPRECATED = {
    'Default_Values': {
        'MOT_rep': int(np.ceil((Initial_Values['MOT_duration'] * 1e6) / Config.MOT_pulse_len)),
        # TODO: use Phases.PULSE_1 instead
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
        'PrePulse_CH2_freq': 133.325e6,  # Hz

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
    'Continuous': {
        # Note, this mode is different than others. It is, well, continuous.
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
    'Imaging': {
        'Triggering_Phase': 'Pulse_1',
        'Pulse_1_CH1_Freq_f': Initial_Values['MOT_freq'],
        'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'],
        'N_Snaps': 1,
        'Buffer_Cycles': 0,
        'Imaging_Phase': 'Pulse_1'
    },
    'OD_FS': {
        'Triggering_Phase': 'Pulse_1',
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
    'Depump': {
        'Triggering_Phase': 'Pulse_1',
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
    'Transit_Exp': {
        'Triggering_Phase': 'Free_Fall',
        'Fountain_final_Delta_freq': 0.45e6,  # 0.38e6 - until 30.10.22
        # 'Fountain_final_Delta_freq': 0,  # 0.38e6 - until 30.10.22
        'Pulse_1_CH1_Freq_f': Initial_Values['MOT_freq'],
        'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'],
        # 'Pulse_1_Repump_amp': 0.000001,
        'N_Snaps': 1,
        'Buffer_Cycles': 0,
        'Imaging_Phase': 'Pulse_1',
        'PrePulse_duration': 10,  # [msec]
        'Shutter_open_time': 5,  # [msec]
        'Pulse_1_duration': (Config.readout_CRUS_pulse_len) / 1e6,  # [msec]
        'M_time': (Config.readout_CRUS_pulse_len) / 1e6,  # [msec]
        'M_off_time': 5,  # [msec]
        ## If with fountain:
        'Fountain_duration': 0.5,  # [msec]
        'Fountain_prep_duration': 0.5,  # [msec], Can't be zero!!!
        'M_window': int(Config.readout_CRUS_pulse_len),  # [nsec]
        # 'PGC_duration': 5  # [msec] EXTREMELY IMOPRTANT for OPX-QuadRF sync
    },
    'Spectrum_Exp': {
        'Triggering_Phase': 'Free_Fall',
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
    'CRUS_Exp': {
        'Triggering_Phase': 'Free_Fall',
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
    'SPRINT_Exp':  {
        'Triggering_Phase': 'Free_Fall',
        'MOT_rep': int(
            np.ceil((Initial_Values['MOT_duration'] * 1e6) / Config_Sprint.MOT_pulse_len)),
        'Fountain_final_Delta_freq': 0.45e6,  # 0.38e6 - until 30.10.22
        'PrePulse_Repump_amp': 0.000001,  # relative
        'PrePulse_CH2_freq': 133.325e6, # Hz
        'Pulse_1_CH1_Freq_f': Initial_Values['MOT_freq'],
        'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'] + 30e6,
        'Pulse_1_Repump_amp': 0.000001,
        'N_Snaps': 1,
        'Buffer_Cycles': 0,
        'Imaging_Phase': 'Pulse_1',
        'PrePulse_duration': 13,  # [msec]
        'Shutter_open_time': 3.5,  # [msec]
        'Pulse_1_duration': int(max(Config_Sprint.readout_pulse_sprint_len_N,
                                    Config_Sprint.readout_pulse_sprint_len_S)) / 1e6,  # [msec]
        ## If with fountain:
        'Fountain_duration': 0.5,  # [msec]
        'Fountain_prep_duration': 0.5,  # [msec], Can't be zero!!!
        'M_window': int(max(Config_Sprint.readout_pulse_sprint_len_N,
                            Config_Sprint.readout_pulse_sprint_len_S)), # [nsec]
        'M_time': int(max(Config_Sprint.readout_pulse_sprint_len_N,
                          Config_Sprint.readout_pulse_sprint_len_S)) / 1e6,  # Pulse_length[nsec] * 1000 repetitions * (Bandwidth[MHz] * frequency steps[MHz]) * 4 / 1e6[nsec/msec] - [msec]
        'M_off_time': 5,  # [msec] - should be at least 5 ms, to sync quadrf and OPX
    },
    'QRAM_Exp':  {
        'Triggering_Phase': 'Free_Fall',
        'MOT_rep': int(np.ceil((Initial_Values['MOT_duration'] * 1e6) / Config_QRAM.MOT_pulse_len)),
        'Fountain_final_Delta_freq': 0.45e6,  # 0.38e6 - until 30.10.22
        'PrePulse_Repump_amp': 0.000001,  # relative
        'PrePulse_CH2_freq': 133.325e6,  # Hz  (QuadRF - Depump)
        'Pulse_1_CH1_Freq_f': Initial_Values['MOT_freq'],  # (QuadRF - AOM - Top F2)
        'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'] + 30e6,  # (QuadRF - Repump)
        'Pulse_1_Repump_amp': 0.000001,
        'N_Snaps': 1,
        'Buffer_Cycles': 0,
        'Imaging_Phase': 'Pulse_1',
        # 'PrePulse_duration': 4,  # [msec]
        # 'PrePulse_duration': 14,  # [msec]
        'PrePulse_duration': 12,  # [msec]
        'Shutter_open_time': 1,  # [msec]
        'Pulse_1_duration': int(max(Config_QRAM.readout_pulse_sprint_len_N,
                                      Config_QRAM.readout_pulse_sprint_len_S)) / 1e6,  # [msec]
        ## If with fountain:
        'Fountain_duration': 0.5,  # [msec]
        'Fountain_prep_duration': 0.5,  # [msec], Can't be zero!!!
        'M_window': int(max(Config_QRAM.readout_pulse_sprint_len_N,
                              Config_QRAM.readout_pulse_sprint_len_S)), # [nsec]
        'M_time': int(max(Config_QRAM.readout_pulse_sprint_len_N,
                            Config_QRAM.readout_pulse_sprint_len_S)) / 1e6,  # Pulse_length[nsec] * 1000 repetitions * (Bandwidth[MHz] * frequency steps[MHz]) * 4 / 1e6[nsec/msec] - [msec]
        'M_off_time': 1,  # [msec] - should be at least 5 ms, to sync quadrf and OPX
    },
    'Transits_Exp':  {
        'Triggering_Phase': 'Free_Fall',
        'MOT_rep': int(np.ceil((Initial_Values['MOT_duration'] * 1e6) / Config_QRAM.MOT_pulse_len)),
        'Fountain_final_Delta_freq': 0.45e6,  # 0.38e6 - until 30.10.22
        'PrePulse_Repump_amp': 1,  # relative
        'PrePulse_CH2_freq': 133.325e6, # Hz
        'Pulse_1_CH1_Freq_f': Initial_Values['MOT_freq'],
        'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'],
        'Pulse_1_Repump_amp': 1,
        'N_Snaps': 1,
        'Buffer_Cycles': 0,
        'Imaging_Phase': 'Pulse_1',
        # 'PrePulse_duration': 4,  # [msec]
        # 'PrePulse_duration': 14,  # [msec]
        'PrePulse_duration': 12,  # [msec]
        'Shutter_open_time': 3.5,  # [msec]
        'Pulse_1_duration': int(max(Config_QRAM.readout_pulse_sprint_len_N,
                                      Config_QRAM.readout_pulse_sprint_len_S)) / 1e6,  # [msec]
        ## If with fountain:
        'Fountain_duration': 0.5,  # [msec]
        'Fountain_prep_duration': 0.5,  # [msec], Can't be zero!!!
        'M_window': int(max(Config_QRAM.readout_pulse_sprint_len_N,
                              Config_QRAM.readout_pulse_sprint_len_S)), # [nsec]
        'M_time': int(max(Config_QRAM.readout_pulse_sprint_len_N,
                            Config_QRAM.readout_pulse_sprint_len_S)) / 1e6,  # Pulse_length[nsec] * 1000 repetitions * (Bandwidth[MHz] * frequency steps[MHz]) * 4 / 1e6[nsec/msec] - [msec]
        'M_off_time': 1,  # [msec] - should be at least 5 ms, to sync quadrf and OPX
    },
    'SPRINT_2-3_Exp':  {
        'Triggering_Phase': 'Free_Fall',
        'Fountain_final_Delta_freq': 0.45e6,  # 0.38e6 - until 30.10.22
        'PrePulse_Repump_amp': 0.000001,  # relative
        'PrePulse_CH2_freq': 133.325e6,  # Hz
        'Pulse_1_CH1_Freq_f': Initial_Values['MOT_freq'],
        'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'] + 30e6,
        # 'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'],
        'Pulse_1_Repump_amp': 0.000001,
        'N_Snaps': 1,
        'Buffer_Cycles': 0,
        'Imaging_Phase': 'Pulse_1',
        # 'PrePulse_duration': 3,  # [msec]
        'PrePulse_duration': 13,  # [msec]
        'Shutter_open_time': 5,  # [msec]
        'Pulse_1_duration': 8,  # [msec]
        ## If with fountain:
        'Fountain_duration': 0.5,  # [msec]
        'Fountain_prep_duration': 0.5,  # [msec], Can't be zero!!!
        'M_window': int(8e6), # [nsec]
        'M_time': 8,  # Pulse_length[nsec] * 1000 repetitions * (Bandwidth[MHz] * frequency steps[MHz]) * 4 / 1e6[nsec/msec] - [msec]
        'M_off_time': 5,  # [msec] - should be at least 5 ms, to sync quadrf and OPX
    },
    'PrePGC_Fountain': {
        'Triggering_Phase': 'Pulse_1',
        'Fountain_final_Delta_freq': 0.45e6,
        'PrePulse_CH2_freq': 133.325e6,  # Hz #Ziv Added for Cooling optimization
        'Pulse_1_CH1_Freq_f': Initial_Values['Flash_freq'],
        'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'],
        # Fountain
        'Fountain_duration': 0.5,        # [msec]
        'Fountain_prep_duration': 0.5,  # [msec], Can't be zero!!!
        # Imaging
        'Imaging_Phase': 'Pulse_1',
        'PrePulse_duration': 1,  # [msec]
        'Pulse_1_duration': 0.2,       # [msec]
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
    'Video': {
        'N_Snaps': 20,
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

# TODO: Remove all the below - the one liner using merge_multiple_jsons
# What it does:
# 1. Starts with taking Initial Values
# 2. Adds Default Values on top of it
# 3. Adds the values of the operation mode
#
# def updateOperationMode(OM = None, values = None):
#     values = Initial_Values if values is None else values
#     if OM is None:
#         OM = Initial_Values['Operation_Mode']
#     if OM not in Operation_Modes:
#         print('\033[91m' + f'Warning: {OM} is not a known operation mode.' + '\033[0m')
#     for key in Operation_Modes[OM]:
#         values[key] = Operation_Modes[OM][key]
#     return values
#
# Initial_Values = updateOperationMode(OM = 'Default_Values') # First, make sure default values are there
# Initial_Values = updateOperationMode()                      # Then, update the values pertaining to the specific operation mode

# Merge all configurations - (a) initial_values; (b) default_values; (c) the operation mode values
#Initial_Values = Utils.merge_multiple_jsons([Initial_Values, Operation_Modes['Default_Values'], Operation_Modes[Initial_Values['Operation_Mode']]])
# TODO: Eventually, we simply need to merge them (manually) to one json, and get rid of Operation_Modes__DEPRECATED['Default_Values']
Default_Values = Utils.merge_multiple_jsons([Initial_Values, Operation_Modes__DEPRECATED['Default_Values']])

#-------------------------------------------------------------------------------------------------------------------
# Frequency Scan Values
#-------------------------------------------------------------------------------------------------------------------
calibrationData = {
    'ch_2': {
               'calibration_func': calibrate,
               'calibrate_RF_Power': calibration_data['calibrate_RF_Power'],
               'calibrate_Optical_Power': 0.2
    }
}

start_delay_time = 0 if Initial_Values['Triggering_Phase'] == 'Pulse_1' else Initial_Values['PrePulse_duration']
Frequency_Scan_Config_Values = {
    'Start_Delay_time': start_delay_time,  # [ms], Time to wait before scan begins
    'Saw_Tooth_Scan': True,  # If TrueL jump. If false, reverse frequency direction (tiangle)
    'Frequency_Range': (93e6, 133e6),  # Order matters. If lower freq is first, we scan up. Otherwise: scan down.
    'Start_Frequency': None,  # if None, start at first frequency. Otherwise, start scan closest to @Start_Frequency as possible
    'Total_Scan_Duration': Initial_Values['Pulse_1_duration'],  # [ms] This should be the same as the TOTAL length of the measurement
    'N_Different_Frequencies': 21,  # This defines the resolution of the scan
    'N_Scan_Repetitions': 4,  # Number of time we go through all frequencies during @Total_Scan_Duration
    'opticalPowerCalibrationFunction': calibrationData  # Used for calibration of the power of the lock AOM (TOP1)
}

#-------------------------------------------
# IOParametersMapping
#
# These are channels in OPX, and should all be integers.
# The values have no real significance, except for in Parameters update in OPX control
#-------------------------------------------
# TODO: Remove this once we see everything works fine with IOP
IOParametersMapping_ZZZ = {
  "MOT_switch": 1,
  "Linear_PGC_switch": 2,
  "Transit_Exp_switch": 3,
  "Spectrum_Exp_switch": 4,
  "SPRINT_Exp_switch": 4,
  "CRUS_Exp_switch": 4,
  "QRAM_Exp_switch": 4,

  "Experiment_Switch": 4,

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
  "Fountain_prep_time": 32,
  "Fountain_final_Delta_freq": 39,
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

#----------------------------------------
# Key_to_Channel Map
#
# This map indicates for each key, what channels are affected by a value change.
# This way, only channels which we actually change are updated on the QuadRF
# (e.g. send MOT_duration to both channels 1 & 4)
#----------------------------------------
Key_to_Channel = {
    'Operation_Mode': [],
    'Triggering_Phase': [1, 4],
    'MOT_freq': [1],
    'PGC_final_freq': [1],
    'Fountain_initial_freq': [1],
    'Fountain_final_freq': [1],
    'Flash_freq': [1],
    'AOM_Repump_freq': [4],
    'MOT_AOM_freq': [2],
    'MOT_duration': [1, 4],        # [msec]
    'MOT_rep': [1, 4],
    'PGC_duration': [1, 4],
    'Post_MOT_delay': [1, 4],
    'PGC_prep_duration': [1, 4],      # [msec]
    'PGC_final_amp': [1],
    'PGC_final_amp_0': [1],       # Relative AOM amplitude between 0 to 1
    'PGC_final_amp_plus': [1],       # Relative AOM amplitude between 0 to 1
    'Repump_PGC_freq': [4],     # The final frequency of the Repump in the PGC sequence
    'Fountain_duration': [1, 4],    # [msec]
    'Fountain_prep_duration': [1, 4],          # [msec]
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
    'FreeFall_duration': [1, 4],  # We use 1.4 because there is a limit for the pulse duration. = 63ms
    'OD_FS_sleep': [1, 4],
    # Imaging parameters:
    'AOM_Off_Detuning': [4],      # [Hz]; Frequency for essentially turnning light coming from AONs off
    'Snapshot_Intervals': [1, 4],   # [msec]
    'PrePulse_duration': [1, 4],            # [msec]
    'Trigger_delay': [1, 4],          # [msec]
    'Pulse_1_duration': [1, 4],       # [msec]
    'InterPulses_duration': [1, 4], # [msec]; Time between pulse 1 & 2.
    'Pulse_2_duration': [1, 4],      # [msec]; If 0 -> no snd pulse
    'Zeeman_delay': [1, 4],          # [msec]
    'N_Snaps': [1, 4],
    # 'Buffer_Cycles': [1,2,3,4]
    'Pulse_1_decay_duration': [1, 4],
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

# Run a sanity check on keys in all dictionaries
Utils.verify_keys_for_case_sensitivity([Initial_Values, Key_to_Channel])


