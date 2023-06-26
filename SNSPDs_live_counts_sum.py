                                                                                                                     # %%
from quadRFMOTController import QuadRFMOTController
from qm.QuantumMachinesManager import QuantumMachinesManager
from qm.qua import *
from qm import SimulationConfig
from qm import LoopbackInterface
from scipy import signal
from scipy.fftpack import fft, fftshift
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import tkinter as tk
import numpy as np
import time
import math
import json
import csv

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(np.abs(array) - value)).argmin()
    return idx


def r_squared(y, vec, apd):
    yhat = y(vec)  # fi
    ybar = np.sum(apd) / len(apd)  # y_mean
    ssreg = np.sum((yhat - ybar) ** 2)  # regression sum of squares
    sstot = np.sum((apd - ybar) ** 2)  # total sum of squares
    r_squared = ssreg / sstot
    return r_squared

def calc_cmat(correction_vars):
    """
    calculating the correction matrix required for IQ mixer using the variable \theta   k
    :param correction_vars: array of correction variables
    theta = correction_vars[0]
    theta = correction_vars[1]
    :return: correction matrix
    """

    theta = correction_vars[0]
    k = correction_vars[1]
    R = [[np.sin(theta), np.cos(theta)], [np.cos(theta), np.sin(theta)]]
    c = [[k, 0], [0, 1 / k]]
    return np.matmul(c, R).flatten().tolist()


QMm = QuantumMachinesManager(host='132.77.54.230', port='80')

controller = 'con1'

# Parameters:0

# Pulse_durations
readout_pulse_len = int(50 * 1e6)
north_const_pulse_len = 500
south_const_pulse_len = 500
analyzer_const_pulse_len = 500
MOT_pulse_len = int(50 * 1e6)
PGC_pulse_len = 500
Probe_pulse_len =500
Fountain_pulse_len = 500
Depump_pulse_len = 500
OD_pulse_len = int(50 * 1e6)
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

IF_Divert = 20e6
#IF_AOM_N = 127.1e6
#IF_AOM_S = 90e6
IF_AOM_SigmaPlus = 114.58e6
IF_AOM_SigmaMinus = 114.58e6
IF_AOM_Pi = 75.34e6
IF_AOM_Analyzer = np.abs(IF_AOM_N - IF_AOM_S) * 2
IF_CRUS_pulser = 125e6

SQUARE_wf = [0.45] * 10000

Super_Sprint_Config = {

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
                10: {'offset': +0.0}, # AOM_Analyzer
            },
            'digital_outputs': {
                1: {},  # Switch AOM + / AOM 2-2'
                2: {},  # Switch AOM - / AOM 2-3'
                3: {},  # AntiHelmholtz Coils
                5: {},  # Camtrigger
                7: {},  # Trigger STIRAP
                8: {},  # Trigger FS
                9: {},  # trigger
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
                "AntiHelmholtz": {
                    "port": (controller, 3),
                    "delay": 0,
                    "buffer": 0,
                },
            },
            'operations': {
                'AntiHelmholtz_MOT': "AntiHelmholtz_on",
            },
        },

        "PULSER_E/L": {
            "singleInput": {
                "port": (controller, 6),
            },
            'operations': {
                'Const_open': "MOT_lock",
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
                "FS_North": {
                    "port": (controller, 5),
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
                'single':  'zero_wf'
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
            }
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
            'sample': 0.45
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
        "trig_wf0": {
            "samples": [(1, 0)]
        }
    },
}

qm_ss = QMm.open_qm(Super_Sprint_Config)
QMm.clear_all_job_results()

with program() as dig:
    # QuadRFMOTController(initialValues={'Operation_Mode': 'Continuous', 'CH1_freq': '113MHz', 'CH1_amp': '10.074226808222061dbm'},
    # QuadRFMOTController(initialValues={'Operation_Mode': 'Continuous', 'CH1_freq': '113MHz', 'CH1_amp': '16.95dbm'},
    #                     updateChannels=[1], debugging=False, continuous=False)  # updates values on QuadRF (uploads table) #
    QuadRFMOTController(initialValues={'Operation_Mode': 'Continuous', 'CH3_freq': '90MHz', 'CH3_amp': '31dbm'},
                        updateChannels=[3], debugging=False, continuous=False)  # updates values on QuadRF (uploads table) #

    counts1 = declare(int)
    counts2 = declare(int)
    counts3 = declare(int)
    counts4 = declare(int)
    counts5 = declare(int)
    counts6 = declare(int)
    counts7 = declare(int)
    counts8 = declare(int)
    counts9 = declare(int)
    counts10 = declare(int)

    counts_st_N = declare_stream()
    counts_st_S = declare_stream()
    counts_st1 = declare_stream()
    counts_st2 = declare_stream()
    counts_st3 = declare_stream()
    counts_st4 = declare_stream()
    counts_st5 = declare_stream()
    counts_st6 = declare_stream()
    counts_st7 = declare_stream()
    counts_st8 = declare_stream()
    counts_st9 = declare_stream()
    counts_st10 = declare_stream()

    n = declare(int)

    m_window = OD_pulse_len  # [nsec]
    # diff = declare(int)
    # g2 = declare(fixed, size=m_window)
    # g2_idx = declare(int)
    # g2_st = declare_stream()
    Measuring_time = 500 * 1e6  # [nsec]
    rep = int(Measuring_time / m_window)
    # with infinite_loop_():
    #     play("AntiHelmholtz_MOT", "AntiHelmholtz_Coils")
        # play("Depump", "AOM_2-2'")
        # play("MOT", "AOM_TOP_1")
    with infinite_loop_():
        with for_(n, 0, n < rep, n+1):

            # play("Const_open_triggered", "PULSER_N")
            play("Const_open", "PULSER_N")
            play("Const_open", "PULSER_S")
            play("Const_open", "PULSER_E/L")
            # play("Square_Pulse", "PULSER_LO")
            # play("Const_open"*amp(0.7), "PULSER_LO")
            # play("AntiHelmholtz_MOT", "AntiHelmholtz_Coils")
            # play("CRUS_pulse", "Pulser_CRUS")

            # measure("OD_measure", "digital_detectors_S", None,
            measure("OD_measure", "digital_detectors_N", None,
                    counting.digital(counts1, m_window, element_outputs="out1"),
                    counting.digital(counts2, m_window, element_outputs="out2"),
                    counting.digital(counts3, m_window, element_outputs="out3"),
                    counting.digital(counts4, m_window, element_outputs="out4"),
                    counting.digital(counts5, m_window, element_outputs="out5"),
                    counting.digital(counts6, m_window, element_outputs="out6"),
                    counting.digital(counts7, m_window, element_outputs="out7"),
                    counting.digital(counts8, m_window, element_outputs="out8"),
                    counting.digital(counts9, m_window, element_outputs="out9"),
                    counting.digital(counts10, m_window, element_outputs="out10"),
                    )

            ## Save Data: ##
            save(counts1, counts_st1)
            save(counts2, counts_st2)
            save(counts3, counts_st3)
            save(counts4, counts_st4)
            save(counts5, counts_st5)
            save(counts6, counts_st6)
            save(counts7, counts_st7)
            save(counts8, counts_st8)
            save(counts9, counts_st9)
            save(counts10, counts_st10)

    with stream_processing():
        counts_st1.buffer(rep).save("avg_counts_1")
        counts_st2.buffer(rep).save("avg_counts_2")
        counts_st3.buffer(rep).save("avg_counts_3")
        counts_st4.buffer(rep).save("avg_counts_4")
        counts_st5.buffer(rep).save("avg_counts_5")
        counts_st6.buffer(rep).save("avg_counts_6")
        counts_st7.buffer(rep).save("avg_counts_7")
        counts_st8.buffer(rep).save("avg_counts_8")
        counts_st9.buffer(rep).save("avg_counts_9")
        counts_st10.buffer(rep).save("avg_counts_10")


job = qm_ss.execute(dig)
avg_count1_handle = job.result_handles.get('avg_counts_1')
avg_count2_handle = job.result_handles.get('avg_counts_2')
avg_count3_handle = job.result_handles.get('avg_counts_3')
avg_count4_handle = job.result_handles.get('avg_counts_4')
avg_count5_handle = job.result_handles.get('avg_counts_5')
avg_count6_handle = job.result_handles.get('avg_counts_6')
avg_count7_handle = job.result_handles.get('avg_counts_7')
avg_count8_handle = job.result_handles.get('avg_counts_8')
avg_count9_handle = job.result_handles.get('avg_counts_9')
avg_count10_handle = job.result_handles.get('avg_counts_10')

avg_count1_handle.wait_for_values(1)
avg_count2_handle.wait_for_values(1)
avg_count3_handle.wait_for_values(1)
avg_count4_handle.wait_for_values(1)
avg_count5_handle.wait_for_values(1)
avg_count6_handle.wait_for_values(1)
avg_count7_handle.wait_for_values(1)
avg_count8_handle.wait_for_values(1)
avg_count9_handle.wait_for_values(1)
avg_count10_handle.wait_for_values(1)

south_vals = []
north_vals = []
SPCMs_vals = []

fig = plt.figure()
font = font_manager.FontProperties(family='Comic Sans MS', weight='bold', style='normal', size=16)

while avg_count1_handle.is_processing():

    avg_counts_res1 = avg_count1_handle.fetch_all()
    avg_counts_res2 = avg_count2_handle.fetch_all()
    avg_counts_res3 = avg_count3_handle.fetch_all()
    avg_counts_res4 = avg_count4_handle.fetch_all()
    avg_counts_res5 = avg_count5_handle.fetch_all()
    avg_counts_res6 = avg_count6_handle.fetch_all()
    avg_counts_res7 = avg_count7_handle.fetch_all()
    avg_counts_res8 = avg_count8_handle.fetch_all()
    avg_counts_res9 = avg_count9_handle.fetch_all()
    avg_counts_res10 = avg_count10_handle.fetch_all()

    print(str(sum(avg_counts_res1 + avg_counts_res2 + avg_counts_res3 + avg_counts_res4 + avg_counts_res8)) + ' , ' + str(sum(avg_counts_res6 + avg_counts_res7 + avg_counts_res5)))

    # plot:
    # north_vals.append(sum(avg_counts_res1 + avg_counts_res2 + avg_counts_res3 + avg_counts_res4))
    # south_vals.append(sum(avg_counts_res5 + avg_counts_res6 + avg_counts_res7 + avg_counts_res8))
    #Ziv
    south_vals.append(sum(avg_counts_res6 + avg_counts_res7 + avg_counts_res5))
    north_vals.append(sum(avg_counts_res1 + avg_counts_res2 + avg_counts_res3 + avg_counts_res4 + avg_counts_res8))
    SPCMs_vals.append(sum(avg_counts_res9 + avg_counts_res10))
    ##
    N_counts = "{:,}".format(north_vals[-1]*2)
    S_counts = "{:,}".format(south_vals[-1]*2)
    SPCMs_counts = "{:,}".format(SPCMs_vals[-1]*2)

    plt.clf()
    plt.plot(north_vals[-100:], label='North Counts: ' + N_counts + ' Hz')
    plt.plot(south_vals[-100:], label='South Counts: ' + S_counts + ' Hz')
    plt.plot(SPCMs_vals[-100:], label='SPCMs Counts: ' + SPCMs_counts + ' Hz')
    plt.title("counts")
    plt.legend(loc='upper left', prop=font)
    plt.show()

    plt.pause(0.5)


print('finished')
