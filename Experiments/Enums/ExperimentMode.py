from enum import Enum
from enum import auto

# -----------------------------------------------------------------------
# Experiment Mode Enumeration
#
# - LIVE - Experiment runs in lab with all devices connected
# - OFFLINE - Experiment does not connect to lab devices
# - PLAYBACK - Pre-recorded data is replayed
# -----------------------------------------------------------------------


class ExperimentMode(Enum):
    LIVE = auto()
    OFFLINE = auto()
    PLAYBACK = auto()
