import inspect
from enum import Enum
# TODO: Document all parameters in this file

# --------------------------------------------------------------------------------------------------------
# IOParameters
#
# IO_Parameters are used to bridge between Python and OPX code. Each Enum is a parameter that can be controlled
# and updated from the Python code (or console) and the OPX can relate to it and act upon the change.
#
# The int values have no real significance (as long as the OPX code also used these ENUMs and not hard-coded values)
# --------------------------------------------------------------------------------------------------------

class IOParameters(Enum):
    MOT_SWITCH_ON = 1
    MOT_SWITCH_OFF = 2  # TODO: Do we want to leave this or make a change in OPX code?
    EXPERIMENT_SWITCH = 4  # Turn experiment on/off
    ANTIHELMHOLTZ_DELAY_SWITCH = 5
    MAX_PROBE_COUNTS_SWITCH = 6
    ANTIHELMHOLTZ_DELAY = 10
    MOT_DURATION = 11
    POST_MOT_DELAY = 12
    MOT_REPETITION = 13
    PGC_DURATION = 20
    PGC_PREP_DURATION = 21
    PGC_INITIAL_AMP_0 = 22
    PGC_INITIAL_AMP_MINUS = 23
    PGC_INITIAL_AMP_PLUS = 24
    #  PGC_FINAL_AMP = -1  # TODO: IS THIS NEEDED?
    PGC_FINAL_AMP_0 = 25   # TODO: IS THIS THE BEST NAME FOR IT?
    PGC_FINAL_AMP_MINUS = 26
    PGC_FINAL_AMP_PLUS = 27
    PGC_PULSE_DURATION_0 = 28
    PGC_PULSE_DURATION_MINUS = 29
    PGC_PULSE_DURATION_PLUS = 30
    FOUNTAIN_DURATION = 31
    FOUNTAIN_PREP_TIME = 32
    FOUNTAIN_PULSE_DURATION_0 = 36
    FOUNTAIN_PULSE_DURATION_MINUS = 37
    FOUNTAIN_PULSE_DURATION_PLUS = 38
    FOUNTAIN_FINAL_DELTA_FREQ = 39
    TRIGGER_DELAY = 41
    PREPULSE_DURATION = 42
    PULSE_1_DURATION = 43
    PULSE_1_DECAY_DURATION = 44
    N_SNAPS = 45
    BUFFER_CYCLES = 46
    PREPULSE_CH1_FREQ = 49
    OD_FREQUENCY = 50
    OD_FS_START = 51
    OD_FS_PULSE_DURATION = 52
    OD_FS_PULSES_SPACING = 53
    OD_DELAY = 54
    OD_SNSPDS_DURATION = 55
    OD_CONTINUOUS_ATTENUATION = 56

    DEPUMP_START = 86
    DEPUMP_PULSE_DURATION = 87
    DEPUMP_PULSES_SPACING = 88
    SHUTTER_OPENING_TIME = 89

    MW_SPEC_FREQ = 60
    MW_SPEC_MW_PULSE = 61
    MW_SPEC_OD_PULSE = 62
    MW_SPEC_REPETITION_TIMES = 63
    MW_SPEC_DELTA_FREQ = 64

    PUSHBEAM_DURATION = 70
    PUSHBEAM_AMPLITUDE = 71
    PUSHBEAM_FREQUENCY = 72

    @classmethod
    def has(cls, key):
        """
        Given a string key, check if it appears in the list of Enums defined above as an IO Parameter
        Key is case-insensitive - it is always converted to upper-case.
        """
        return key.upper() in cls.__members__

    @classmethod
    def value_from_string(cls, key):
        """
        Given a string key, return the integer value assigned to the specific IO Parameter.
        """
        return cls.__members__[key.upper()].value

    # def __init__(self):
    #     self._dict = {}
    #     for i in inspect.getmembers(self):
    #         # Ignore private/protected functions and methods that do not start with a underscore
    #         if not i[0].startswith('_') and not inspect.ismethod(i[1]):
    #             self._dict[i[1]] = i[0]
    #
    # def __val__(self, name):
    #     pass
    #
    # def __str__(self, key):
    #     return self._dict[key]
    #
    # def __repr__(self, key):
    #     return self._dict[key]
