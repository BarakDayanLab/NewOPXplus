import os
from enum import Enum, auto
import playsound


class SOUNDS(Enum):
    DEFAULT: auto()
    INDICATOR_1: auto()
    INDICATOR_2: auto()
    WARN: auto()
    FATAL: auto()


class BDSound:

    def __init__(self):

        # TODO: make this dependant on the user - are we on what machine?
        self.media_folder = r'C:\Windows\Media'
        pass

    def play(self, sound_enum):

        name = 'ding'  # default sound
        if sound_enum == SOUNDS.INDICATOR_1:
            name = 'ding'
        elif sound_enum == SOUNDS.INDICATOR_2:
            name = 'chord'
        elif sound_enum == SOUNDS.WARN:
            name = 'chimes'
        elif sound_enum == SOUNDS.FATAL:
            name = 'Alarm05'
        else:
            print(f'Unknown sound type: {sound_enum}')

        sound_file = os.path.join(self.media_folder, name + '.wav')
        playsound(sound_file)
