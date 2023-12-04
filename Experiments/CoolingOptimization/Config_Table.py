from Experiments.BaseExperiment.Config_Table import Initial_Values


# These values will be added/override the ones in Initial_Values (e.g. Common Experiment Values)
Experiment_Values = {
    'Experiment_Name': 'CoolingOptimization',

    'Triggering_Phase': 'Pulse_1',
    'Fountain_final_Delta_freq': 0.45e6,
    'PrePulse_CH2_freq': 133.325e6,  # Hz #Ziv Added for Cooling optimization
    'Pulse_1_CH1_Freq_f': Initial_Values['Flash_freq'],
    'Pulse_1_CH4_Freq': Initial_Values['AOM_Repump_freq'],
    # Fountain
    'Fountain_duration': 0.5,  # [msec]
    'Fountain_prep_duration': 0.5,  # [msec], Can't be zero!!!
    # Imaging
    'Imaging_Phase': 'Pulse_1',
    'PrePulse_duration': 1,  # [msec]
    'Pulse_1_duration': 0.2,  # [msec]
}