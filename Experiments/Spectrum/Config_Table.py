import numpy as np
from Experiments.Enums.Phases import Phases
from Experiments.BaseExperiment.Config_Table import Initial_Values
import Experiments.Spectrum.Config_Experiment as Config


# These values will be added/override the ones in Initial_Values (e.g. Common Experiment Values)
Experiment_Values = {
    'Experiment_Name': 'Spectrum',
    'Phases_Order': [Phases.MOT, Phases.FOUNTAIN, Phases.PGC, Phases.FREE_FALL, Phases.PULSE_1, Phases.INTER_PULSES, Phases.PULSE_2, Phases.POST_PULSE],
    'Triggering_Phase': Phases.FREE_FALL,
    'MOT_rep': int(np.ceil((Initial_Values['MOT_duration'] * 1e6) / Config.MOT_pulse_len)),
    'Fountain_final_Delta_freq': 0.45e6,  # 0.38e6 - until 30.10.22
    'PrePulse_Repump_amp': 1,  # relative
    'PrePulse_CH2_freq': 133.325e6, # Hz
    'Pulse_1_CH1_Freq_f': Initial_Values['MOT_freq'],
    'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'],
    'Pulse_1_Repump_amp': 1,
    'N_Snaps': 1,
    'Buffer_Cycles': 0,
    'Imaging_Phase': Phases.PULSE_1,
    'PrePulse_duration': 10,  # [msec]
    'Shutter_open_time': 3.5,  # [msec]
    'Pulse_1_duration': int(max(Config.readout_pulse_sprint_len_N, Config.readout_pulse_sprint_len_S)) / 1e6,  # [msec]
    'Fountain_duration': 0.5,  # [msec]  # Make duration 0 to disable fountain! # TODO: should have a dedicated flag
    'Fountain_prep_duration': 0.5,  # [msec], Can't be zero!!!
    'M_window': int(max(Config.readout_pulse_sprint_len_N, Config.readout_pulse_sprint_len_S)), # [nsec]
    'M_time': int(max(Config.readout_pulse_sprint_len_N, Config.readout_pulse_sprint_len_S)) / 1e6,  # Pulse_length[nsec] * 1000 repetitions * (Bandwidth[MHz] * frequency steps[MHz]) * 4 / 1e6[nsec/msec] - [msec]
    'M_off_time': 1,  # [msec] - should be at least 5 ms, to sync quadrf and OPX

    'OPX_Quad_Misalignment_Delay': 4000  # = 4us [ns]
}
