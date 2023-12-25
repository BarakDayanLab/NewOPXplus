from enum import Enum
from enum import auto

# -----------------------------------------------------------------------
# Termination Reason Enumeration
#
# - SUCCESS - Experiment completed as planned
# - USER - The user terminated (e.g. pressed ESC, or other option that the code presents)
# - ERROR - Some internal fatal error occurred that prevents the proper run of the experiment
# -----------------------------------------------------------------------


class TerminationReason(Enum):
    USER = auto()
    ERROR = auto()
    SUCCESS = auto()
    PLAYBACK_END = auto()
