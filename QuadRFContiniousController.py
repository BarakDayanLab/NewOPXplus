from quadRFMOTController import QuadRFMOTController
import numpy as np

QuadRFMOTController(initialValues= {'Operation_Mode': 'Continuous',
                                    'CH3_freq': '110MHz',
                                    'CH3_amp': '0.2dbm'}, updateChannels=[3], debugging=True)    # In order to turn a channel off, set amp to '0x0'
