from enum import Enum
import numpy as np
from Experiments.BaseExperiment.Phases import Phases
import Experiments.BaseExperiment.Config_Experiment as Config
from UtilityResources.AOMCalibration import calibrate, calibration_data
from Utilities.Utils import Utils


# -------------------------------------------
# Initial Values
# -------------------------------------------

Initial_Values = {
    'Experiment_Name': 'Default',

    'Imaging_Phase': Phases.PULSE_1,
    'Triggering_Phase': -1,  # Don't change this. It is advised that each experiment will define its own trigger phase
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
    #'Snap_Time': 1,  # [msec]
    'Trigger_delay': -0.1,  # [msec]
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
    'PrePulse_CH1_freq': -1,  # [Hz]

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

Default_Values = {
    'MOT_rep': int(np.ceil((Initial_Values['MOT_duration'] * 1e6) / Config.MOT_pulse_len)),
    'Triggering_Phase': Phases.PULSE_1,
    'Pulse_1_CH1_Freq_i': Initial_Values['MOT_freq'],
    # 'Pulse_1_CH1_Freq_i': Initial_Values['Flash_freq'],
    'Pulse_1_CH1_Freq_f': Initial_Values['MOT_freq'],
    # 'Pulse_1_CH1_Freq_f': Initial_Values['Flash_freq'],
    'Pulse_1_CH_2_3_Freq': Initial_Values['MOT_AOM_freq'],  # [Hz], usually should be equal to MOT_AOM_freq
    'Pulse_1_CH4_Freq': Initial_Values['Repump_PGC_freq'],  # [Hz], usually should be equal to Repump_PGC_freq
    'Pulse_2_CH1_Freq': Initial_Values['MOT_freq'],
    'Pulse_2_CH_2_3_Freq': Initial_Values['MOT_AOM_freq'],
    'Pulse_2_CH4_Freq': Initial_Values['AOM_Repump_freq'],

    'PrePulse_CH1_freq': Initial_Values['MOT_freq'] - 10e6,
    'PrePulse_CH2_freq': 133.325e6,  # Hz

    # Fountain
    'Fountain_initial_Delta_freq': 0,  # By default, should be taken from previous phase
    'Fountain_initial_freq': Initial_Values['MOT_freq'],  # By default, should be taken from previous phase
    'Fountain_final_freq': Initial_Values['MOT_freq'],  # By default, should be taken from next (?) phase
    'PGC_initial_amp': 1,  # By default, should be taken from previous phase
    'PGC_initial_amp_0': 1,  # By default, should be taken from previous phase
    'PGC_initial_amp_minus': 1,  # By default, should be taken from previous phase
    'PGC_initial_amp_plus': 1,  # By default, should be taken from previous phase
    'PGC_initial_freq': Initial_Values['MOT_freq'],  # By default, should be taken from previous phase
    'PostPulses_freq': Initial_Values['MOT_freq'],
}

# -------------------------------------------------------------------------------------------------------------------
# Frequency Scan Values
# -------------------------------------------------------------------------------------------------------------------
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
    'PrePulse_CH1_freq': [1],  # [Hz]
    'Trigger_delay': [1, 4],  # [msec]
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


