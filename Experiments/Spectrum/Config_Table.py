from Experiments.BaseExperiment.Config_Table import Initial_Values
import Experiments.Spectrum.Config_Experiment as Config


# These values will be added/override the ones in Initial_Values (e.g. Common Experiment Values)
Experiment_Values = {
    'Experiment_Name': 'Spectrum',

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
    'M_window': int(Config.readout_pulse_spectrum_len),  # [nsec]
    'M_time': len(Config.Spectrum_Exp_Gaussian_samples) * 1000 * 21 * 4 / 1e6
}
