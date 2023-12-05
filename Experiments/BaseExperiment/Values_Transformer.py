from Utilities.BDLogger import BDLogger
from Experiments.BaseExperiment import Config_Experiment as Config  # Attempt to load the default config (may be overriden later)

#----------------------------------------
# Values_Factor Map
#
# This map holds factoring for values - each key has an array:
# - One element - just setting it
# - Two elements - multiplying by the factor and casting to the required type
# - Three elements - calling a transformation function and casting
#
# E.g. "Prep_duration" is written in ms, we have to factor it into an int value, but in 4ns, for the OPX.
# Both the updateValue uses factors and also the convenience functions


class Values_Transformer:

    MODES = ['SET_VALUE', 'FACTOR_AND_CAST', 'TRANSFORM']
    
    def __init__(self):
        self.logger = BDLogger()
        # Create internal dictionary from Values Factors
        self._entries_dict = {}
        for k, v in Values_Factor.items():
            self._entries_dict[k.upper()] = v
        pass

    def knows(self, key):
        """
        Returns True/False - whether this transformer deals with this key
        """
        key = key.upper()
        return key in self._entries_dict

    def mode(self, key):
        """
        Returns the mode that the transformer has registered for this key
        """
        key = key.upper()
        if not self.knows(key):
            self.logger.error(f'Cannot find {key} in Values_Factor map')
        return self.MODES[len(self._entries_dict[key])-1]

    def factor_and_cast(self, key, value):
        # TODO
        pass

    def transform(self, key, value):
        # TODO: get the string - find the function and invoke it
        pass

Values_Factor = {
    # TODO: we believe we do not need these in Values_Factor - can we remove them?
    # "Transit_Exp_switch": [bool, None, Transit_Exp_switch],
    # "Spectrum_Exp_switch": [bool, None, Spectrum_Exp_switch],
    # "QRAM_Exp_switch": [bool, None, QRAM_Exp_switch],

    "Experiment_Switch": [bool],

    # TODO: do we need these? Just for updating values thru the console? Are we switching modes during experiment?
    "Transit_Exp_switch": [bool],
    "Spectrum_Exp_switch": [bool],
    "QRAM_Exp_switch": [bool],

    # TODO: this is the only setting that relies on Config - can we replace it with something?
    'MOT_duration': [int, 1e6 / Config.MOT_pulse_len],
    'AntiHelmholtz_delay': [int, 1e6 / 4],  # [msec]
    'Post_MOT_delay': [int, 1e6 / 4],  # [msec]
    
    # TODO: is the below needed?
    # 'N_Snaps': [int, None, update_N_Snaps],  # TODO: Dor said it's not in use
    
    'Buffer_Cycles': [int, 1],
    ## PGC parameters ##
    "PGC_duration": [int, 1e6 / 4],  # [msec]
    "PGC_prep_duration": [int, 1e6 / 4],

    ## Fountain parameters ##
    "Pre_PGC_Fountain_duration": [int, 1e6 / 4],
    "Fountain_duration": [int, 1e6 / 4],
    
    # TODO - do we need these for SPECTRUM experiment? OR - we can stay w/o the function calls
    # "Fountain_prep_duration": [int, None, update_fountain_prep_time],
    # "Fountain_final_Delta_freq": [float, None, update_fountain_Delta_freq],
    # "Fountain_aom_chirp_rate": [float, None, update_fountain_Delta_freq],
    
    "Fountain_prep_duration": [int, 1e6 / 4],
    "Fountain_final_Delta_freq": [float, (1e3) / (1e6 / 4) * 4],
    # self.fountain_aom_chirp_rate = int(df * 1e3 / (self.fountain_prep_duration * 4))  # mHz/nsec
    "Fountain_aom_chirp_rate": [float, (1e3) / (1e6 / 4) * 4],
    # self.fountain_aom_chirp_rate = int(df * 1e3 / (self.fountain_prep_duration * 4))  # mHz/nsec

    # The frequency in the frame of reference where the atoms have just enough velocity to reach the toroids [Hz].
    # Free-Fall parameters:
    # The time take the atoms to reach the Toroids and fall back = 2 * sqrt(2 * (Toroid height[m] - MOT height[m]) / g[m/sec^2]) * 1e3[msec]
    'FreeFall_duration': [int, 1e6 / 4],
    # We use 1.4 because there is a limit for the pulse duration. = 63ms
    'Coil_timing': [int, 1e6 / 4],  # Zeeman coils turn on time from start of free fall[msec]
    
    ## Imaging parameters ##
    'Snapshot_Intervals': [int, 1e6 / 4],  # [msec]
    'PrePulse_duration': [int, 1e6 / 4],  # [msec]
    'Trigger_delay': [int, 1e6 / 4],  # [msec]
    'Pulse_1_duration': [int, 1e6 / 4],  # [msec]
    'Pulse_1_decay_duration': [int, 1e6 / 4],  # [msec]
    'InterPulses_duration': [int, 1e6 / 4],
    'Pulse_2_duration': [int, 1e6 / 4],  # [msec]
    
    ## OD parameters ##
    'OD_FS_Start': [int, 1e6 / 4],  # [msec]
    'OD_FS_pulse_duration': [int, 1e6 / 4],  # [msec]
    'OD_FS_pulses_spacing': [int, 1e6 / 4],  # [msec]
    'Depump_Start': [int, 1e6 / 4],  # [msec]
    'Depump_pulse_duration': [int, 1e6 / 4],  # [msec]
    'Depump_pulses_spacing': [int, 1e6 / 4],  # [msec]
    'OD_duration': [int, 1e6 / 4],  # [msec]
    'Wait_duration': [int, 1e6 / 4],  # [msec]
    'Shutter_open_time': [int, 1e6 / 4],  # [msec]
    
    # SNSPDs measurement parameters:
    'M_delay': [int, 1e6 / 4],  # [msec]
    'OD_delay': [int, 1e6 / 4],  # [msec]
    'OD_duration_pulse1': [int, 1e6 / 4],  # [msec]
    'OD_sleep': [int, 1e6 / 4],  # [msec]
    'OD_duration_pulse2': [int, 1e6 / 4],  # [msec]
    'M_time': [int, 1e6 / 4]  # [msec]'
}


