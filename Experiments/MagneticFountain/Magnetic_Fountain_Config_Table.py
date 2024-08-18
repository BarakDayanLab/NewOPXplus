from Experiments.Enums.Phases import Phases
from Experiments.BaseExperiment.Config_Table import Initial_Values


# These values will be added/override the ones in Initial_Values (e.g. Common Experiment Values)
Experiment_Values = {
    'Experiment_Name': 'MagneticFountain',
    'Triggering_Phase': Phases.PULSE_1,
    'Phases_Order': [Phases.MOT, Phases.FOUNTAIN, Phases.PGC, Phases.FREE_FALL, Phases.PULSE_1, Phases.INTER_PULSES, Phases.PULSE_2, Phases.POST_PULSE],
    # 'Fountain_final_Delta_freq': 0.30e6, #17.03.24
    # 'Fountain_final_Delta_freq': 0.33e6, #0.37e6
    'Fountain_final_Delta_freq': 0.46e6,
    'PrePulse_CH2_freq': 133.325e6, # for depump during mot, change back to 70 mhz for killing depump # Hz #Ziv Added for Cooling optimization
    'Pulse_1_CH1_Freq_f': Initial_Values['Flash_freq'],
    'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq_Off'],
    # Fountain
    'Fountain_duration': 0.5,  # [msec]
    'Fountain_prep_duration': 0.5,  # [msec], Can't be zero!!!
    # Imaging
    'Imaging_Phase': Phases.PULSE_1,
    'PrePulse_duration': 1,  # [msec]
    'Pulse_1_duration': 0.2,  # [msec]

    # 'OPX_Quad_Misalignment_Delay': 4000  # = 4us [ns]
    'OPX_Quad_Misalignment_Delay': 4000 - 140000,  # = 4us [ns]
}
