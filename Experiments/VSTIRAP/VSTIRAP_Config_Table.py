import numpy as np
from Experiments.Enums.Phases import Phases
from Experiments.BaseExperiment.Config_Table import Initial_Values
import Experiments.VSTIRAP.VSTIRAP_Config_Experiment as Config



# These values will be added/override the ones in Initial_Values (e.g. Common Experiment Values)
Experiment_Values = {
    'Experiment_Name': 'VSTIRAP',
    'Triggering_Phase': Phases.FREE_FALL,
    'Phases_Order': [Phases.MOT, Phases.FOUNTAIN, Phases.PGC, Phases.FREE_FALL, Phases.PULSE_1, Phases.INTER_PULSES, Phases.PULSE_2, Phases.POST_PULSE],
    'MOT_rep': int(np.ceil((Initial_Values['MOT_duration'] * 1e6) / Config.MOT_pulse_len)),
    # 'Fountain_final_Delta_freq': 0.33e6,  # 0.37e6 - until 26.02.24
      'Fountain_final_Delta_freq': 0.46e6,
    # 'Fountain_final_Delta_freq': 0.35e6,  # 0.315e6 - until 28.02.24
    'PrePulse_Repump_amp': 0.000001,  # relative  (QuadRF)
    'PrePulse_CH2_freq': 133.325e6,  # Hz  (QuadRF - Depump)
    'Pulse_1_CH1_Freq_f': Initial_Values['MOT_freq'],  # (QuadRF - AOM - Top F2)
    'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'] + 30e6,  # (QuadRF - Repump)
    'Pulse_1_Repump_amp': 1e-6,  # This value is translated to dB and is added to the amplitude of the channel meaning: Ch_Amp = Ch_Amp + 10 * np.log10('Pulse_1_Repump_amp')
    'Pulse_1_Depump_amp': 1e-6,  # This value is translated to dB and is added to the amplitude of the channel meaning: Ch_Amp = Ch_Amp + 10 * np.log10('Pulse_1_Depump_amp')
    'N_Snaps': 1,
    'Buffer_Cycles': 0,
    'Imaging_Phase': Phases.PULSE_1,
    # 'PrePulse_duration': 4,  # [msec]
    'PrePulse_duration': 20,  # [msec]
    # 'PrePulse_duration': 20,  # [msec]
    'Shutter_open_time': 1,  # [msec]
    'Pulse_1_duration': int(max(Config.readout_pulse_sprint_len_N, Config.readout_pulse_sprint_len_S)) / 1e6,  # [msec]
    ## If with fountain:
    'Fountain_duration': 0.5,  # [msec]
    'Fountain_prep_duration': 0.5,  # [msec], Can't be zero!!!
    'M_window': int(max(Config.readout_pulse_sprint_len_N, Config.readout_pulse_sprint_len_S)),  # [nsec]
    'ignored_marginals': 5e5,  # [nsec] data at the beginning and at the end of a window to "throw away" because of shutters noise
    'M_time': int(max(Config.readout_pulse_sprint_len_N, Config.readout_pulse_sprint_len_S)) / 1e6,
    # Pulse_length[nsec] * 1000 repbetitions * (Bandwidth[MHz] * frequency steps[MHz]) * 4 / 1e6[nsec/msec] - [msec]
    'M_off_time': 1.5,  # [msec] - should be at least 5 ms, to sync quadrf and OPX

    # TOP F1 pulser freq - for transits 2-3' experiment using the pulser as repump
    'CH3_continuous_freq': '130MHz',  # QuadRF # This is used in QuadRF as ch3 that goes to TopF1 pulser
    'CH3_continuous_amp': '31dbm',  # QuadRF

    # TODO: we should have a different value here - and fix opx code to use this value
    'OPX_Quad_Misalignment_Delay': 4000,  # = 4us [ns]


    'VSTIRAP_pulser_N_amp': 0.0,
    'VSTIRAP_pulser_S_amp': 0.0,
    'VSTIRAP_beam_amp': 0.015
    # 'VSTIRAP_beam_amp': 0.015
}
