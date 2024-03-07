import numpy as np


controller = 'con1'

Measure_Time = 100  # ms

# Parameters:0

# Pulse_durations
readout_pulse_len = int(50 * 1e3)
north_const_pulse_len = 500
south_const_pulse_len = 500
analyzer_const_pulse_len = 500
MOT_pulse_len = int(50 * 1e3)
PGC_pulse_len = 500
Probe_pulse_len = 500
Fountain_pulse_len = 500
Depump_pulse_len = 500
OD_pulse_len = int(50 * 1e3)
Repump_pulse_len = 500
Trigger_pulse_len = 1
AntiHelmholtz_pulse_len = 60

# Intermediate frequencies
IF_TOP1_MOT = 113e6
IF_TOP1_PGC = 103e6
IF_TOP1_Flash = 121.6625e6
IF_TOP2 = 90e6
IF_AOM_MOT = 110e6
# IF_AOM_OD = 133.325e6
IF_AOM_OD = 136.325e6
IF_AOM_OD = 93e6
IF_AOM_Depump = 133.325e6
# IF_AOM_Depump = 139.325e6
IF_AOM_Repump = 78.4735e6
IF_AOM_N = 129.2368e6
IF_AOM_S = 129.2368e6
IF_AOM_LO = 129.2368e6
IF_AOMs_MZ = 110e6
IF_AOM_Anc = 184e6
IF_AOM_Spectrum = 133.325e6 / 2

IF_Divert = 20e6
# IF_AOM_N = 127.1e6
# IF_AOM_S = 90e6
IF_AOM_SigmaPlus = 114.58e6
IF_AOM_SigmaMinus = 114.58e6
IF_AOM_Pi = 75.34e6
IF_AOM_Analyzer = np.abs(IF_AOM_N - IF_AOM_S) * 2
IF_CRUS_pulser = 125e6

SQUARE_wf = [0.45] * 10000

config = {

    'version': 1,

    'controllers': {
        controller: {
            'type': 'opx1',
            'analog_outputs': {
                1: {'offset': +0.0},  # AOM Lock TOP 1
                2: {'offset': +0.0},  # AOM main TOP 2
                3: {'offset': +0.0},  # AOM 0
                4: {'offset': +0.0},  # AOM + / AOM 2-2'
                5: {'offset': +0.0},  # AOM - / AOM 2-3'
                6: {'offset': +0.0},  # AOM 1-2' - Repump
                7: {'offset': 0.0},  # MW_I

                # changed for MW_SPEC
                # 7: {'offset': +0.0},  # Pulse_EOM
                8: {'offset': +0.0},  # AOM_N

                9: {'offset': +0.0},  # AOM_S
                10: {'offset': +0.0},  # AOM_Analyzer
            },
            'digital_outputs': {
                1: {},  # Switch AOM + / AOM 2-2'
                2: {},  # Switch AOM - / AOM 2-3'
                3: {},  # AntiHelmholtz Coils
                4: {},  # Helmholtz Coils
                5: {},  # Camtrigger
                6: {},  # APDs
                7: {},  # Trigger STIRAP
                8: {},  # Trigger FS
                9: {},  # trigger
                10: {},  # EOM pulses
            },
            'analog_inputs': {
                1: {'offset': +0.0},  # DET10
                2: {'offset': +0.0},  # Summing amp

            },
            'digital_inputs': {
                1: {'polarity': 'RISING', "threshold": 0.5, "deadtime": 4},
                2: {'polarity': 'RISING', "threshold": 0.5, "deadtime": 4},
                3: {'polarity': 'RISING', "threshold": 0.5, "deadtime": 4},
                4: {'polarity': 'RISING', "threshold": 0.5, "deadtime": 4},
                5: {'polarity': 'RISING', "threshold": 0.5, "deadtime": 4},
                6: {'polarity': 'RISING', "threshold": 0.5, "deadtime": 4},
                7: {'polarity': 'RISING', "threshold": 0.5, "deadtime": 4},
                8: {'polarity': 'RISING', "threshold": 0.5, "deadtime": 4},
                9: {'polarity': 'RISING', "threshold": 0.5, "deadtime": 4},
                10: {'polarity': 'RISING', "threshold": 0.5, "deadtime": 4},
            },
        }
    },

    'elements': {
        # "AOM_TOP_1": {
        #     "thread": "a",
        #     "singleInput": {
        #         "port": (controller, 1)
        #     },
        #     'operations': {
        #         'MOT': "MOT_lock",
        #     },
        #     'intermediate_frequency': IF_TOP1_MOT,
        #     # 'hold_offset': {'duration': 10},
        # },
        #
        "AOM_2-2'": {
            "thread": "a",
            "singleInput": {
                "port": (controller, 1)
            },
            'operations': {
                'Depump': "Depump_pulse",
            },
            'intermediate_frequency': IF_AOM_Depump,
            # 'hold_offset': {'duration': 10},
        },

        "Pulser_CRUS": {
            "singleInput": {
                "port": (controller, 8)
            },
            # 'digitalInputs': {
            #     "Shutter_Switch": {
            #         "port": (controller, 5),
            #         "delay": 0,
            #         "buffer": 0,
            #     },
            # },
            'operations': {
                'CRUS_pulse': "CRUS_pulses",
            },
            'intermediate_frequency': IF_CRUS_pulser,
        },

        "digital_detectors_N": {
            "singleInput": {
                "port": (controller, 5)
            },
            'digitalInputs': {
                "switch1": {
                    "port": (controller, 2),
                    "delay": 0,
                    "buffer": 0,
                },
                "switch2": {
                    "port": (controller, 5),
                    "delay": 0,
                    "buffer": 0,
                },
            },
            "digitalOutputs": {
                "out1": (controller, 1),
                "out2": (controller, 2),
                "out3": (controller, 3),
                "out4": (controller, 4),
                "out5": (controller, 5),
                "out6": (controller, 6),
                "out7": (controller, 7),
                "out8": (controller, 8),
                "out9": (controller, 9),
                "out10": (controller, 10),
            },
            'outputs': {
                'out1': (controller, 1)
            },
            'operations': {
                'readout': "digital_readout",
                'OD_measure': "OD_pulse"
            },
            'time_of_flight': 1200 + 28,
            'smearing': 0,
            'intermediate_frequency': IF_AOM_OD,
        },

        "digital_detectors_S": {
            "singleInput": {
                "port": (controller, 1)
            },
            'digitalInputs': {
                "switch1": {
                    "port": (controller, 2),
                    "delay": 0,
                    "buffer": 0,
                },
                "switch2": {
                    "port": (controller, 5),
                    "delay": 0,
                    "buffer": 0,
                },
            },
            "digitalOutputs": {
                "out1": (controller, 1),
                "out2": (controller, 2),
                "out3": (controller, 3),
                "out4": (controller, 4),
                "out5": (controller, 5),
                "out6": (controller, 6),
                "out7": (controller, 7),
                "out8": (controller, 8),
                "out9": (controller, 9),
                "out10": (controller, 10),
            },
            'outputs': {
                'out1': (controller, 1)
            },
            'operations': {
                'readout': "digital_readout",
                'OD_measure': "OD_pulse"
            },
            'time_of_flight': 1200 + 28,
            'smearing': 0,
            'intermediate_frequency': IF_TOP2,
        },

        "AntiHelmholtz_Coils": {
            'digitalInputs': {
                # "AntiHelmholtz": {
                #     "port": (controller, 3),
                #     "delay": 0,
                #     "buffer": 0,
                # },
                "Helmholtz": {
                    "port": (controller, 4),
                    "delay": 0,
                    "buffer": 0,
                },
            },
            'operations': {
                'AntiHelmholtz_MOT': "AntiHelmholtz_on",
            },
        },

        "PULSER_E": {
            "singleInput": {
                "port": (controller, 6),
            },
            'digitalInputs': {
                "AWG_Switch": {
                    "port": (controller, 10),
                    # "delay": 176, # OPX control EOM
                    "delay": 160,  # OPX control EOM double pass
                    # "delay": 400, # AWG control EOM
                    "buffer": 0,
                },
            },
            'operations': {
                'Const_open': "Pulser_ON",
                'Const_early_open': "Pulser_ON_early",
                'Square_Pulse': "square_pulse",
            },
            # 'intermediate_frequency': IF_AOM_LO,
            'intermediate_frequency': IF_AOMs_MZ,
            # 'intermediate_frequency': IF_AOM_Anc,
        },

        "PULSER_L": {
            "singleInput": {
                "port": (controller, 7),
            },
            'digitalInputs': {
                "AWG_Switch": {
                    "port": (controller, 10),
                    # "delay": 176, # OPX control EOM
                    "delay": 160,  # OPX control EOM double pass
                    # "delay": 400, # AWG control EOM
                    "buffer": 0,
                },
            },
            'operations': {
                'Const_open': "Pulser_ON",
                'Const_late_open': "Pulser_ON_late",
                'Square_Pulse': "square_pulse",
            },
            # 'intermediate_frequency': IF_AOM_LO,
            'intermediate_frequency': IF_AOMs_MZ,
            # 'intermediate_frequency': IF_AOM_Anc,
        },

        "PULSER_LO": {
            "singleInput": {
                "port": (controller, 7),
            },
            'operations': {
                'Const_open': "MOT_lock",
                'Square_Pulse': "square_pulse",
            },
            # 'intermediate_frequency': IF_AOM_LO,
            'intermediate_frequency': IF_AOMs_MZ,
        },

        "AOM_Spectrum": {
            'singleInput': {
                "port": (controller, 8)
            },
            'operations': {
                'Spectrum_pulse': "Frequency_Sweep",
            },
            'intermediate_frequency': IF_AOM_Spectrum,
        },

        "PULSER_N": {
            "singleInput": {
                "port": (controller, 9),
            },
            'digitalInputs': {
                "FS_North": {
                    "port": (controller, 8),
                    "delay": 500,
                    "buffer": 0,
                },
            },
            'operations': {
                'Const_open': "MOT_lock",
                'Const_open_triggered': "MOT_lock_ON",
            },
            'intermediate_frequency': IF_AOM_N,
        },

        "PULSER_S": {
            "singleInput": {
                "port": (controller, 10),
            },
            'digitalInputs': {
                "APD_Switch": {
                    "port": (controller, 6),
                    "delay": 0,
                    "buffer": 0,
                },
                "SouthtoNorth_Shutter": {
                    "port": (controller, 9),
                    "delay": 0,
                    "buffer": 0,
                },
            },
            'operations': {
                'Const_open': "MOT_lock",
                'Const_open_triggered': "MOT_lock_ON",
            },
            'intermediate_frequency': IF_AOM_S,
        },

    },

    "pulses": {
        "digital_readout": {
            'length': 100,
            'operation': 'control',
            'waveforms': {
                'single': 'zero_wf'
            },
            'digital_marker': 'ON'
        },
        "OD_pulse": {
            'operation': 'control',
            'length': OD_pulse_len,
            'waveforms': {
                'single': 'const_wf'
            },
            'digital_marker': 'ON'
        },
        "Depump_pulse": {
            'operation': 'control',
            'length': MOT_pulse_len,
            'waveforms': {
                'single': 'const_wf'
            },
        },
        "MOT_lock": {
            'operation': 'control',
            'length': MOT_pulse_len,
            'waveforms': {
                'single': 'const_wf'
            },
            # 'digital_marker': 'ON',
        },
        "Pulser_ON": {
            'operation': 'control',
            'length': MOT_pulse_len,
            'waveforms': {
                'single': 'const_wf'
            },
            # 'digital_marker': 'ON',
            'digital_marker': 'Trig_EOM_MZ'
        },
        "Pulser_ON_early": {
            'operation': 'control',
            'length': MOT_pulse_len,
            'waveforms': {
                'single': 'const_early_wf'
            },
            'digital_marker': 'ON',
        },
        "Pulser_ON_late": {
            'operation': 'control',
            'length': MOT_pulse_len,
            'waveforms': {
                'single': 'const_late_wf'
            },
            'digital_marker': 'ON',
        },
        "MOT_lock_ON": {
            'operation': 'control',
            'length': MOT_pulse_len,
            'waveforms': {
                'single': 'const_wf'
            },
            'digital_marker': 'ON',
            # 'digital_marker': 'ON_pulse_switch'
        },
        "square_pulse": {
            'operation': 'control',
            'length': len(SQUARE_wf),
            'waveforms': {
                'single': 'square_wf'
            }
        },
        "Frequency_Sweep": {
            'operation': 'control',
            'length': MOT_pulse_len,
            'waveforms': {
                'single': 'const_wf'
            },
        },
        "AntiHelmholtz_on": {
            'operation': 'control',
            'length': MOT_pulse_len,
            'digital_marker': 'ON'
        },
        "CRUS_pulses": {
            'operation': 'control',
            'length': OD_pulse_len,
            'waveforms': {
                'single': 'const_wf'
            },
        },
    },

    'integration_weights': {
        'integW': {
            'cosine': [(1.0, int(readout_pulse_len))],
            'sine': [(0.0, int(readout_pulse_len))]
        },
        'Det_opt': {
            'cosine': [1.0],
            'sine': [0.0]
        }
    },

    "waveforms": {
        'zero_wf': {
            'type': 'constant',
            'sample': 0.0
        },
        'MOT_wf': {
            'type': 'constant',
            'sample': 0.3
        },
        'const_wf': {
            'type': 'constant',
            # 'sample': 0.1
            'sample': 0.35
        },
        'const_late_wf': {
            'type': 'constant',
            # 'sample': 0.1
            'sample': 0.4
            # 'sample': 0.29
        },
        'const_early_wf': {
            'type': 'constant',
            # 'sample': 0.1
            'sample': 0.49
        },
        'square_wf': {
            'type': 'arbitrary',
            # 'sample': 0.1
            'samples': SQUARE_wf
        },
    },

    "digital_waveforms": {
        "ON": {
            "samples": [(1, 0)]
        },
        "Trig": {
            "samples": [(1, 20), (0, 0)]
        },
        "Trig_EOM_MZ": {
            "samples": [(0, 250), (1, 250)] * int(MOT_pulse_len / 500) + [(0, 0)]
            # "samples": [(0, 270), (1, 230)] * MZ_balancing_seq_rep + [(0, 0)]
        },
        "trig_wf0": {
            "samples": [(1, 0)]
        }
    },
}


streams = {
    "Detector_1_Avg_Counts": {
        "number": 1,
        "name": "Detector_1_Avg_Counts",
        "type": "int",
        "binary": "I",  # Unsigned int
        "playback": "Det1_avg_counts.npz",
        "save_raw": True
    },
    "Detector_2_Avg_Counts": {
        "number": 2,
        "name": "Detector_2_Avg_Counts",
        "type": "int",
        "binary": "I",  # Unsigned int
        "playback": "Det2_avg_counts.npz",
        "save_raw": True
    },
    "Detector_3_Avg_Counts": {
        "number": 1,
        "name": "Detector_3_Avg_Counts",
        "type": "int",
        "binary": "I",  # Unsigned int
        "playback": "Det3_avg_counts.npz",
        "save_raw": True
    },
    "Detector_4_Avg_Counts": {
        "number": 1,
        "name": "Detector_4_Avg_Counts",
        "type": "int",
        "binary": "I",  # Unsigned int
        "playback": "Det4_avg_counts.npz",
        "save_raw": True
    },
    "Detector_5_Avg_Counts": {
        "number": 1,
        "name": "Detector_5_Avg_Counts",
        "type": "int",
        "binary": "I",  # Unsigned int
        "playback": "Det5_avg_counts.npz",
        "save_raw": True
    },
    "Detector_6_Avg_Counts": {
        "number": 1,
        "name": "Detector_6_Avg_Counts",
        "type": "int",
        "binary": "I",  # Unsigned int
        "playback": "Det6_avg_counts.npz",
        "save_raw": True
    },
    "Detector_7_Avg_Counts": {
        "number": 1,
        "name": "Detector_7_Avg_Counts",
        "type": "int",
        "binary": "I",  # Unsigned int
        "playback": "Det7_avg_counts.npz",
        "save_raw": True
    },
    "Detector_8_Avg_Counts": {
        "number": 1,
        "name": "Detector_8_Avg_Counts",
        "type": "int",
        "binary": "I",  # Unsigned int
        "playback": "Det8_avg_counts.npz",
        "save_raw": True
    },
    "Detector_9_Avg_Counts": {
        "number": 1,
        "name": "Detector_9_Avg_Counts",
        "type": "int",
        "binary": "I",  # Unsigned int
        "playback": "Det9_avg_counts.npz",
        "save_raw": True
    },
    "Detector_10_Avg_Counts": {
        "number": 1,
        "name": "Detector_10_Avg_Counts",
        "type": "int",
        "binary": "I",  # Unsigned int
        "playback": "Det10_avg_counts.npz",
        "save_raw": True
    }
}