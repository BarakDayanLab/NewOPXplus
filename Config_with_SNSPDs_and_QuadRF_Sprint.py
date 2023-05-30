import numpy as np
from scipy import signal
import math
import json
import sys
import logging
from logging import StreamHandler, Formatter, INFO, WARN, ERROR
from matplotlib import pyplot as plt
import math

def gaussian(x, mu, sigma):
    return 1 / sigma / np.sqrt(2 * np.pi) * np.exp(-1.0 / 2.0 * ((x - mu) / sigma) ** 2)


def gaussian2(x, mu, sigma):
    return 0.2 * np.exp(-1.0 / 2.0 * ((x - mu) / sigma) ** 2)


def gauss_adaptive(amplitude, length):
    t = np.linspace(-length / 2, length / 2, length)
    sigma = length / 6  # defining the width of the pulse to match the length desired
    gauss_wave = amplitude * np.exp(-(t ** 2) / (2 * sigma ** 2))
    return [float(x) for x in gauss_wave]


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


controller = 'con1'

# Parameters:

# time tags vector size

# parameters of sizes
vec_size = 10000
num_of_detectors = 6
opx_max_per_window = vec_size*num_of_detectors/2
dets_number = [1, 2, 3, 6, 7, 8]

# Pulse_durations
readout_pulse_len = 10000000
MOT_pulse_len = 10e6
Short_pulse_len = 40
PGC_pulse_len = 40
Fountain_pulse_len = 40
FreeFall_pulse_len = 40
Probe_pulse_len = 5e6
Depump_pulse_len = 40
Measuring_pulse_len = 40
OD_pulse_len = 10e6
Repump_pulse_len = 10e6
north_const_pulse_len = 500
south_const_pulse_len = 500
analyzer_const_pulse_len = 500
Detection_pulse_len = 1000
Flash_pulse_len = 10000
spectrum_pulse = 240

# Intermediate frequencies
IF_TOP1_MOT = 113e6
IF_TOP1_PGC = 93e6
IF_TOP1_Flash = 121.6625e6
IF_TOP2 = 90e6
IF_AOM_MOT = 110e6
IF_AOM_MOT_OFF = 80e6
IF_AOM_OD = 92675000 # (226 - 266.65 / 2) * 1e6
# IF_AOM_OD = 133.325e6
IF_AOM_Depump = 133.325e6
IF_AOM_Repump = 78.4735e6
# IF_AOM_N = 89.2368e6
IF_AOM_N = 129.2368e6
# IF_AOM_S = 89.2368e6
IF_AOM_S = 129.2368e6
# IF_AOM_LO = 89.2368e6
IF_AOM_LO = 129.2368e6
IF_AOM_SigmaPlus = 114.58e6
IF_AOM_SigmaMinus = 114.58e6
IF_AOM_Pi = 75.34e6
IF_CRUS_pulser = 125e6

IF_Divert = 20e6
IF_AOM_Analyzer = np.abs(IF_AOM_N - IF_AOM_S) * 2

# Waveforms

Linear = 0.45 * np.linspace(1, 0, 1000)
tau = -(1000 - 1) / np.log(0.01)
Exponential = 0.49 * signal.exponential(1000, 0, tau, False)

T_f = 1e4
phase = math.pi / 2
p = np.poly1d([-(6 * math.pi) / (2 * T_f ** 2), (6 * math.pi) / (2 * T_f ** 3)])
t = np.linspace(-T_f, T_f, int(T_f))
PolyCosine = signal.sweep_poly(t, p)
PolySine = signal.sweep_poly(t, p, 90)

trig_samples = [0, 0, 0, 0] + [0.3, 0.3, 0.3, 0.3] * 16 + [0, 0, 0, 0]

Det_Gaussian_samples = ((signal.gaussian(200, std=(24 / 30)) * 0.8 - 0.4).tolist() + [-0.4] * 200) * 4

# SPRINT parameters
# num_of_photons_det_pulses = 1.5 # alpha^2
num_of_photons_det_pulses = 1.65 # alpha^2
# num_of_photons_det_pulses = 2 # alpha^2
num_of_photons_sprint_pulses = 0.13 # alpha^2
# num_of_photons_sprint_pulses = 0 # alpha^2  # For only det pulses sequence

num_of_det_pulses_S = 3
# num_of_det_pulses_S = 2
num_of_sprint_pulses_S = 1
# num_of_sprint_pulses_S = 0  # For only det pulses sequence
num_of_det_pulses_N = 3
# num_of_det_pulses_N = 2
num_of_sprint_pulses_N = 3
# num_of_sprint_pulses_N = 0  # For only det pulses sequence

# parameters for window len
# efficiency = 0.5 # the efficiency of the system
# efficiency = 0.29 # the efficiency of the system
efficiency = 0.168 # the efficiency of the system with sagnac
num_of_photons_per_sequence_S = num_of_photons_det_pulses * num_of_det_pulses_S + num_of_photons_sprint_pulses * num_of_sprint_pulses_S
num_of_photons_per_sequence_N = num_of_photons_det_pulses * num_of_det_pulses_N + num_of_photons_sprint_pulses * num_of_sprint_pulses_N

def Sprint_Exp_Gaussian_samples(sprint_pulse_len=110,det_pulse_len = 30, det_pulses_amp = [0.4]*6,sprint_pulses_amp = [0.4]*4
                                ,num_init_zeros = 10,num_between_zeros = 10, num_fin_zeros = 0):
    Sprint_Exp_Gaussian_samples = [0] * num_init_zeros
    Sprint_Exp_Gaussian_samples += [0] * 12  # Due to undesired reflections (Echos)
    for n in det_pulses_amp:
        # Sprint_Exp_Gaussian_samples += [n] * 50 + [0] * num_between_zeros
        Sprint_Exp_Gaussian_samples += (signal.gaussian(det_pulse_len, std=(det_pulse_len * 0.5 / 2.355)) * n).tolist() + [0] * num_between_zeros
    # Sprint_Exp_Gaussian_samples += (signal.gaussian(prep_pulse_len, std=(prep_pulse_len / 2.355)) * prep_pulse_amp).tolist() + [0] * num_between_zeros
    Sprint_Exp_Gaussian_samples += [0] * 12  # Due to undesired reflections (Echos)
    for m in sprint_pulses_amp:
        # Sprint_Exp_Gaussian_samples += [m] * 110 + [0] * num_between_zeros
        Sprint_Exp_Gaussian_samples += (signal.gaussian((sprint_pulse_len-4), std=((sprint_pulse_len-4) / 2.355)) * m).tolist() + [0] * (num_between_zeros + 4)
    Sprint_Exp_Gaussian_samples += [0] * num_fin_zeros
    # return Sprint_Exp_Gaussian_samples[:-num_between_zeros]
    return Sprint_Exp_Gaussian_samples


def QRAM_Exp_Gaussian_samples(sprint_pulse_len=110, det_pulse_len=30, det_pulses_amp=[0.4]*6, sprint_pulses_amp=[0.4]*4,
                              num_init_zeros=10, num_between_zeros=10, num_fin_zeros=0):
    qram_exp_gaussian_samples = [0] * num_init_zeros
    for n in det_pulses_amp:
        # Sprint_Exp_Gaussian_samples += [n] * 50 + [0] * num_between_zeros
        qram_exp_gaussian_samples += (signal.gaussian(det_pulse_len, std=(det_pulse_len / 2.355)) * n).tolist() + [0] * num_between_zeros
    # Sprint_Exp_Gaussian_samples += (signal.gaussian(prep_pulse_len, std=(prep_pulse_len / 2.355)) * prep_pulse_amp).tolist() + [0] * num_between_zeros
    for m in sprint_pulses_amp:
        # Sprint_Exp_Gaussian_samples += [m] * 110 + [0] * num_between_zeros
        qram_exp_gaussian_samples += (signal.gaussian(sprint_pulse_len, std=(sprint_pulse_len / 2.355)) * m).tolist() + [0] * num_between_zeros
    qram_exp_gaussian_samples += [0] * num_fin_zeros
    return qram_exp_gaussian_samples[:-num_between_zeros]

def QRAM_Exp_samples(delta=240, pulse_len=10000000):
    QRAM_exp_samples = []
    for n in range(int(pulse_len//(3*delta))):
        QRAM_exp_samples += [0.4]*delta + [0]*2*delta
    return QRAM_exp_samples

det_pulse_len = 40
# det_pulse_len = 80
# num_init_zeros_S = 30
# num_init_zeros_S = 10  # For only det pulses sequence
num_init_zeros_S = 0  # For QRAM seq
# num_fin_zeros_S = 20
# num_fin_zeros_S = 10  # For only det pulses sequence
num_fin_zeros_S = 0  # For only det pulses sequence
# num_between_zeros = 10
num_between_zeros = 20
prep_pulse_amp_S = 0.4
prep_pulse_len = 50

# pulses structure
# det_pulse_amp_S = [0.45, 0, 0.45, 0, 0.45, 0]
det_pulse_amp_S = [0.45, 0, 0.45, 0, 0.45, 0, 0.45, 0] # QRAM seq
# det_pulse_amp_S = [0.45, 0, 0.45, 0]  # SPRINT with 2-3'
# det_pulse_amp_S = [0.45, 0, 0, 0, 0, 0]

# sprint_pulse_amp_S = [0.075, 0]
# sprint_pulse_amp_S = [0.075, 0, 0, 0]
# sprint_pulse_amp_S = [0.075, 0.075, 0.075, 0.075]
sprint_pulse_amp_S = [0, 0, 0, 0]

# sprint_pulse_amp_S = [0, 0, 0]
# sprint_pulse_amp_S = []  # For only det pulses sequence
sprint_pulse_len = 100
Sprint_Exp_Gaussian_samples_S = Sprint_Exp_Gaussian_samples(sprint_pulse_len=sprint_pulse_len,
                                                            det_pulse_len=det_pulse_len,
                                                            det_pulses_amp=det_pulse_amp_S,
                                                            sprint_pulses_amp=sprint_pulse_amp_S, num_init_zeros=num_init_zeros_S,
                                                            num_between_zeros=num_between_zeros, num_fin_zeros=num_fin_zeros_S)

# num_init_zeros_N = 30
# num_init_zeros_N = 10  # For only det pulses sequence
num_init_zeros_N = 0  # For QRAM seq
# num_fin_zeros_N = 20
# num_fin_zeros_N = 10  # For only det pulses sequence
num_fin_zeros_N = 0  # For only det pulses sequence
# det_pulse_amp_N = [0, 0.45, 0, 0.45, 0, 0.45]
det_pulse_amp_N = [0, 0.45, 0, 0.45, 0, 0.45, 0, 0.45]
# det_pulse_amp_N = [0, 0.45, 0, 0.45]  # SPRINT with 2-3'
# det_pulse_amp_N = [0.45, 0, 0, 0, 0, 0]
# sprint_pulse_amp_N = [0, 0.085]

# SPRINT N PULSES
#for q-ram z-measurment - for 1/2 power from S=0.075 sprint pulses
# sprint_pulse_amp_N = [0, 0.065, 0, 0.065]
# for sprint
# sprint_pulse_amp_N = [0, 0.085, 0.085, 0.085]
# sprint_pulse_amp_N = [0, 0.085, 0.085, 0]
# sprint_pulse_amp_N = [0.085, 0.085, 0.085, 0]
sprint_pulse_amp_N = [0, 0, 0, 0]

# sprint_pulse_amp_N = [0, 0, 0]
# sprint_pulse_amp_N = []  # For only det pulses sequence

Sprint_Exp_Gaussian_samples_N = Sprint_Exp_Gaussian_samples(sprint_pulse_len=sprint_pulse_len,
                                                            det_pulse_len=det_pulse_len,
                                                            det_pulses_amp=det_pulse_amp_N,
                                                            sprint_pulses_amp=sprint_pulse_amp_N, num_init_zeros=num_init_zeros_N,
                                                            num_between_zeros=num_between_zeros, num_fin_zeros=num_fin_zeros_N)

# readout_pulse_sprint_len_N = math.ceil(((opx_max_per_window/4)/(efficiency*1e6*num_of_photons_per_sequence_N))*len(Sprint_Exp_Gaussian_samples_N))*1e6# [ns] length of the measurment window for North, the 4's are for division in 4
# readout_pulse_sprint_len_N = math.ceil(((opx_max_per_window/1.5)/(efficiency*1e6*num_of_photons_per_sequence_N))*len(Sprint_Exp_Gaussian_samples_N))*1e6# [ns] length of the measurment window for North, the 4's are for division in 4
readout_pulse_sprint_len_N = 10*1e6# [ns] length of the measurment window for North, the 4's are for division in 4
# readout_pulse_sprint_len_S = math.ceil(((opx_max_per_window/4)/(efficiency*1e6*num_of_photons_per_sequence_S))*len(Sprint_Exp_Gaussian_samples_S))*1e6# [ns] length of the measurment window for South, the 4's are for division in 4
# readout_pulse_sprint_len_S = math.ceil(((opx_max_per_window/1.5)/(efficiency*1e6*num_of_photons_per_sequence_S))*len(Sprint_Exp_Gaussian_samples_S))*1e6# [ns] length of the measurment window for South, the 4's are for division in 4
readout_pulse_sprint_len_S = 10*1e6# [ns] length of the measurment window for South, the 4's are for division in 4

SPRINT_Exp_TOP2_samples = [0.45]*int(max(readout_pulse_sprint_len_N, readout_pulse_sprint_len_S))
QRAM_Exp_TOP2_samples = QRAM_Exp_samples(delta=240, pulse_len=24000)

CRUS_probe_samples = [0] * 340 + ([0.4] * (256 * 2 + 10)) + [0] * 162 # twice the 512ns period of the AWG, with scope triggered
CRUS_pulser_samples = [0] * 128 + ([0.4] * (256 + 128)) + [0] * 128 * 4 # twice the 512ns period of the AWG, with scope triggered
readout_CRUS_pulse_len = int((len(CRUS_pulser_samples) / 2) * 40000)
# Det_Gaussian_samples = ((signal.gaussian(100, std=(24 / 2.355)) * 0.4).tolist() + [0] * 100) * 4
SWAP_Gaussian_samples = [0] * 24 + (signal.gaussian(400, std=(100 / 2.355)) * 0.3).tolist() + [0] * 24
DC_cal_samples = [0] * 20 + ([0.3] * 80) + [0] * 20
EOM_pulse_seq_samples = Det_Gaussian_samples + SWAP_Gaussian_samples + DC_cal_samples + [0] * 24

delay = 1
# Det_square_samples = ([0] * 16 + ([0.3] * 64) + [0] * 16) * 4
Det_square_samples = (([0.4] * 100) + [0] * 200) * 4
Det_square_single_samples = ([0.4] * 20) + [0] * 100
SWAP_square_samples = ([0] * 20 + ([0.3] * 408) + [0] * 20)
square_cal_samples = ([0] * 20 + ([0.3] * 88) + [0] * 20)
AOMs_pulse_seq_samples = [0] * delay + Det_square_samples + SWAP_square_samples + square_cal_samples + [0] * (
            24 - delay)

square_samples = ([0] * 400 + [0.49] * 408 + [0] * 400)
x = np.linspace(0, len(square_samples), len(square_samples))
Depump_pulse_samples = np.convolve(square_samples, gaussian(x, 600, 20))[900 // 4 * 4:1500 // 4 * 4]
OD_pulse_samples = np.convolve(square_samples, gaussian(x, 600, 20))[900 // 4 * 4:1500 // 4 * 4]

square_samples2 = ([0] * 400 + [0.8] * 400 + [0] * 400)
xx = np.linspace(0, len(square_samples2), len(square_samples2))
Gaussian_pulse_samples = gauss_adaptive(0.2, 1000)

## Attenuators (global) for AOMs 0, + & - of MOT sequence
# Factor 0.0-1.0
AOM_0_Attenuation = 1 # 0.85 seems to be about right, for 0.8 the slopes are exactly the same
AOM_Plus_Attenuation = 0.42
# AOM_Plus_Attenuation = 0.38
AOM_Minus_Attenuation = 1 # 0.85 seems to be about right, for 0.8 the slopes are exactly the same

## For Homodyne ##
LO_pulse_samples = ([0.4] * 20) + [0] * 100

## MW Spectroscopy parameters

with open('conf_args.json', 'r') as fp:
    confargs = json.load(fp)
dc_i = confargs['offsets']['I']
dc_q = confargs['offsets']['Q']
i_port = confargs['ports']['I']
q_port = confargs['ports']['Q']
lo_freq = confargs['lo_frequency']
im_freq = confargs['im_frequency']
wf_samples = [confargs['wf_samples']['wf1'], confargs['wf_samples']['wf2']]
correction_matrix = calc_cmat([confargs['correction_vars']['theta'], confargs['correction_vars']['k']])
pulse_time = confargs['pulse_time']

config = {

    'version': 1,

    'controllers': {
        controller: {
            'type': 'opx1',
            'analog_outputs': {
                1: {'offset': +0.0},  # TOP2
                2: {'offset': +0.0},  # MOT AOM 0
                3: {'offset': +0.0},  # MOT AOM -
                4: {'offset': +0.0},  # MOT AOM +
                5: {'offset': +0.0},  # AOM 2-3' (OD)
                6: {'offset': +0.0},  # AOM S
                7: {'offset': dc_i},  # MW_I
                8: {'offset': dc_q},  # MW_Q
                9: {'offset': +0.0},  # AOM N
                10: {'offset': +0.0}, # AOM S
            },

            'digital_outputs': {
                1: {},  # Trigger TOP1
                2: {},  # Trigger AOM -/0
                3: {},  # Switch ON/OFF anti-Helmholtz coils (for MOT)
                4: {},  # Switch ON/OFF Helmholtz coils (for Super-SPRINT) and Shutters to SNSPDs
                5: {},  # Pulses for fast switch OFF
                6: {},  # Pulses for fast switch ON
                7: {},  # Measurement Trigger
                8: {},  # MW Trigger
                9: {},  # Blue/Red Trigger
                10: {},  # Probe Trigger
            },

            'analog_inputs': {
                1: {'offset': +0.0},  # FLR
                2: {'offset': +0.0},  # HOMODYNE and not OD Free-Space
                # 2: {'offset': +0.197359321899414038085},  # Summing amp / Homodyne
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

        "Cooling_Sequence": {
            'digitalInputs': {
                "TOP1_Trigger": {
                    "port": (controller, 1),
                    "delay": 0,
                    "buffer": 0,
                },
            },
            'operations': {
                'C_Seq': "C_seq_trigger",
            },
        },

        "MOT_AOM_0": {
            'singleInput': {
                "port": (controller, 2)
            },
            'operations': {
                'MOT': "MOT_lock",
                'Linear': "Linear_pulse",
                'Const': "Const_pulse",
            },
            'intermediate_frequency': IF_AOM_MOT,
        },

        "MOT_AOM_-": {
            'singleInput': {
                "port": (controller, 3)
            },
            'operations': {
                'MOT': "MOT_lock",
                'Linear': "Linear_pulse",
                'Const': "Const_pulse",
            },
            'intermediate_frequency': IF_AOM_MOT,
        },

        "MOT_AOM_+": {
            'singleInput': {
                "port": (controller, 4)
            },
            'operations': {
                'MOT': "MOT_lock",
                'Linear': "Linear_pulse",
                'Const': "Const_pulse",
            },
            'intermediate_frequency': IF_AOM_MOT,
        },

        "AOM_2-3'_for_interference": {
            "singleInput": {
                "port": (controller, 5)
            },
            'digitalInputs': {
                "OD_Switch": {
                    "port": (controller, 2),
                    "delay": 0,
                    "buffer": 0,
                },
            },
            'operations': {
                'OD_FS': "OD_FS_pulse",
            },
            'intermediate_frequency': IF_AOM_OD,
        },

        "AOM_2-2/3'": {
            "singleInput": {
                "port": (controller, 5)
            },
            'digitalInputs': {
                "OD_Switch": {
                    "port": (controller, 2),
                    "delay": 0,
                    "buffer": 0,
                },
                "Shutter_Switch": {
                    "port": (controller, 5),
                    "delay": 0,
                    "buffer": 0,
                },
            },
            'operations': {
                'OD': "OD_pulse",
                'OD_FS': "OD_FS_pulse",
                'Depump': "Depump_pulse",
            },
            'intermediate_frequency': IF_AOM_OD,
        },

        "AOM_2-2'": {
            'singleInput': {
                "port": (controller, 5)
            },
            'operations': {
                'Depump': "Depump_pulse",
            },
            'intermediate_frequency': IF_AOM_Depump,
        },

        "Measurement": {
            'digitalInputs': {
                "Measure_trigger": {
                    "port": (controller, 7),
                    "delay": 0,
                    "buffer": 0,
                },
            },
            'operations': {
                'Measure': "Measure_trigger",
            },
        },

        "detector_with_dig_out": {
            ## fake port ##
            "singleInput": {
                "port": (controller, 1)
            },
            'digitalInputs': {
                "AWG_Switch": {
                    "port": (controller, 10),
                    "delay": 0,
                    "buffer": 0,
                },
            },

            ###############
            "digitalOutputs": {
                "out": (controller, dets_number[0]),
            },
            'outputs': {
                  'out1': (controller, 1)
            },
            'operations': {
                'readout': "digital_readout",
                'readout_SPRINT': "digital_readout_sprint",
                'readout_QRAM': "digital_readout_QRAM",
            },
            'time_of_flight': 36,
            'smearing': 0,
            'intermediate_frequency': IF_TOP2,
        },
        # all other detectors are added at the end of the script as "detectors_no_dig_out_{i}"

        "FLR_detection": {
            # open fake:#
            "singleInput": {
                "port": (controller, 1)
            },
            #############
            "outputs": {
                'out1': (controller, 1)
            },
            'operations': {
                'Detection': "Det_pulse",
            },
            'time_of_flight': 36,
            'smearing': 0,
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

        "Zeeman_Coils": {
            'digitalInputs': {
                "Helmholtz": {
                    "port": (controller, 9),
                    "delay": 0,
                    "buffer": 0,
                },
                # "Shutter_Switch": {
                #     "port": (controller, 5),
                #     "delay": 0,
                #     "buffer": 0,
                # },
            },
            'operations': {
                'ZeemanSplit': "Zeeman_on",
                'ZeemanOFF': "Zeeman_off",
            },
        },

        "AOM_ANALYZER_N": {
            "singleInput": {
                "port": (controller, 6),
            },
            'operations': {
                'Const_open': "MOT_lock",
                'Detection_pulses': "Square_detection_pulses",
                'Homodyne_Pulse': "Homodyne_Pulse",
                'Sprint_experiment_pulses_N': "Gaussian_Sprint_pulse_N"
            },
            'intermediate_frequency': IF_AOM_N,
        },

        "AOM_ANALYZER_S": {
            "singleInput": {
                "port": (controller, 7),
            },
            'operations': {
                'Const_open': "MOT_lock",
                'Detection_pulses': "Square_detection_pulses",
                'Homodyne_Pulse': "Homodyne_Pulse",
                'Sprint_experiment_pulses_S': "Gaussian_Sprint_pulse_N"
            },
            'intermediate_frequency': IF_AOM_N,
        },
        # for bell-experiment
        # "PULSER_ANCILLA": {
        #     "singleInput": {
        #         "port": (controller, 8),
        #     },
        #     'operations': {
        #         'Const_open': "MOT_lock",
        #         'Detection_pulses': "Square_detection_pulses",
        #         'Homodyne_Pulse': "Homodyne_Pulse",
        #         'Sprint_experiment_pulses_N': "Gaussian_Sprint_pulse_N"
        #     },
        #     'intermediate_frequency': IF_AOM_N,
        # },
        # for bell-experiment
        "PULSER_LO": {
            "singleInput": {
                "port": (controller, 8),
            },
            'operations': {
                'Const_open': "MOT_lock",
            },
            'intermediate_frequency': IF_AOM_LO,
        },

        "PULSER_N": {
            "singleInput": {
                "port": (controller, 9),
            },
            'operations': {
                'Const_open': "MOT_lock",
                'Detection_pulses': "Square_detection_pulses",
                'Homodyne_Pulse': "Homodyne_Pulse",
                'Sprint_experiment_pulses_N': "Gaussian_Sprint_pulse_N"
            },
            'intermediate_frequency': IF_AOM_N,
            # 'time_of_flight': 0,
            # 'smearing': 0,
        },

        "PULSER_S": {
            "singleInput": {
                "port": (controller, 10),
            },
            'digitalInputs': {
                "Shutter_Switch": {
                    "port": (controller, 5),
                    "delay": 0,
                    "buffer": 0,
                },
            },
            'operations': {
                'Const_open': "MOT_lock",
                'Detection_pulses': "Square_detection_pulses",
                'Homodyne_Pulse': "Homodyne_Pulse",
                'Sprint_experiment_pulses_S': "Gaussian_Sprint_pulse_S"
            },
            'intermediate_frequency': IF_AOM_S,
        },
    },

    "pulses": {

        "Probe_lock": {
            'operation': 'control',
            'length': Probe_pulse_len,
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
            'digital_marker': 'ON'
        },

        "MOT_lock_ON": {
            'operation': 'control',
            'length': MOT_pulse_len,
            'waveforms': {
                'single': 'const_wf'
            },
            'digital_marker': 'ON_pulse_switch'
        },

        "MOT_lock_OFF": {
            'operation': 'control',
            'length': MOT_pulse_len,
            'waveforms': {
                'single': 'const_wf'
            },
            'digital_marker': 'OFF_pulse_switch'
        },

        "PGC_lock": {
            'operation': 'control',
            'length': PGC_pulse_len,
            'waveforms': {
                'single': 'const_wf'
            }
        },

        "CRUS_pulses": {
            'operation': 'control',
            'length': len(CRUS_pulser_samples),
            'waveforms': {
                'single': 'CRUS_pulser_wf'
            },
        },

        "CRUS_probe_pulses": {
            'operation': 'control',
            'length': len(CRUS_probe_samples),
            'waveforms': {
                'single': 'CRUS_probe_wf'
            },
            'digital_marker': 'ON'
        },

        "Linear_pulse": {
            'operation': 'control',
            'length': len(Linear),
            'waveforms': {
                'single': 'linear_wf'
            }
        },

        "Exponential_pulse": {
            'operation': 'control',
            'length': len(Exponential),
            'waveforms': {
                'single': 'exp_wf'
            }
        },

        "Const_pulse": {
            'operation': 'control',
            'length': Short_pulse_len,
            'waveforms': {
                'single': 'const_wf'
            },
        },

        "Const_pulse_2dBm": {
            'operation': 'control',
            'length': Short_pulse_len,
            'waveforms': {
                'single': 'const_wf_2dBm'
            },
        },

        "OD_pulse": {
            'operation': 'measurement',
            'length': OD_pulse_len,
            'waveforms': {
                'single': 'const_wf'
            },
            'integration_weights': {
                'Detection_opt': 'integW'
            },
            'digital_marker': 'ON'
        },

        "OD_FS_pulse": {
            'operation': 'control',
            'length': OD_pulse_len,
            'waveforms': {
                'single': 'const_wf'
            },
            'digital_marker': 'ON'
        },

        "Depump_pulse": {
            'operation': 'control',
            'length': Depump_pulse_len,
            'waveforms': {
                'single': 'const_wf'
            },
        },

        "Square_detection_pulses": {
            'operation': 'control',
            'length': len(Det_square_samples),
            'waveforms': {
                'single': 'Detection_Square_wf'
            },
        },

        "Gaussian_detection_pulses": {
            'operation': 'control',
            'length': len(Det_Gaussian_samples),
            'waveforms': {
                'single': 'Detection_Gaussian_wf'
            },
        },

        "Gaussian_detection_pulses_trigged": {
            'operation': 'control',
            'length': len(Det_Gaussian_samples),
            'waveforms': {
                'single': 'Detection_Gaussian_wf'
            },
            'digital_marker': 'ON'
        },

        "Gaussian_Sprint_pulse_S": {
            'operation': 'control',
            'length': len(Sprint_Exp_Gaussian_samples_S),
            'waveforms': {
                'single': 'Sprint_Gaussian_wf_S'
            },
            'digital_marker': 'ON'
        },

        "Gaussian_Sprint_pulse_N": {
            'operation': 'control',
            'length': len(Sprint_Exp_Gaussian_samples_N),
            'waveforms': {
                'single': 'Sprint_Gaussian_wf_N'
            }
        },

        "DC_cal_pulse": {
            'operation': 'control',
            'length': len(DC_cal_samples),
            'waveforms': {
                'single': 'DC_cal_wf'
            }
        },

        "Homodyne_Det_pulse": {
            'operation': 'measurement',
            'length': len(LO_pulse_samples),
            'waveforms': {
                'single': 'LO_pulse_wf'
            },
            'integration_weights': {
                'Detection_opt': 'Homodyne_Det_opt'
            },
            'digital_marker': 'ON'
        },

        "Homodyne_Pulse": {
            'operation': 'control',
            'length': len(LO_pulse_samples),
            'waveforms': {
                'single': 'LO_pulse_wf'
            },
        },

        "Det_pulse": {
            'operation': 'measurement',
            'length': MOT_pulse_len,
            'waveforms': {
                'single': 'zero_wf'
            },
            'integration_weights': {
                'Detection_opt': 'integW'
            },
            # 'digital_marker': 'ON'
        },

        "C_seq_trigger": {
            'operation': 'control',
            'length': MOT_pulse_len,
            'digital_marker': 'Trig'
        },

        "Measure_trigger": {
            'operation': 'control',
            'length': Measuring_pulse_len,
            'digital_marker': 'ON'
        },

        "Repump_trigger": {
            'operation': 'control',
            'length': Repump_pulse_len,
            'digital_marker': 'ON'
        },

        "AntiHelmholtz_on": {
            'operation': 'control',
            'length': MOT_pulse_len,
            'digital_marker': 'ON'
        },

        "Zeeman_on": {
            'operation': 'control',
            'length': FreeFall_pulse_len,
            'digital_marker': 'ON'
        },

        "Zeeman_off": {
            'operation': 'control',
            'length': FreeFall_pulse_len,
            'digital_marker': 'OFF'
        },

        "Shutter_on": {
            'operation': 'control',
            'length': OD_pulse_len,
            'digital_marker': 'ON'
        },

        "Snapshot_Flash": {
            'operation': 'control',
            'length': Flash_pulse_len,
            'digital_marker': 'Trig'
        },

        # "PolySine_pulse": {
        #     'operation': 'control',
        #     'length': len(PolySine),
        #     'waveforms': {
        #         'single': 'PolySine_wf'
        #     },
        #     'digital_marker': 'Trig'
        # },

        # "PolyCosine_pulse": {
        #     'operation': 'control',
        #     'length': len(PolyCosine),
        #     'waveforms': {
        #         'single': 'PolyCosine_wf'
        #     },
        #     'digital_marker': 'Trig'
        # },

        "digital_readout": {
            'length': readout_pulse_len,
            'operation': 'control',
            'waveforms': {
                'single': 'zero_wf'
            },
            'digital_marker': 'ON'
        },

        "digital_readout_sprint": {
            'length': int(max(readout_pulse_sprint_len_N, readout_pulse_sprint_len_S)),
            'operation': 'control',
            'waveforms': {
                'single': 'const_wf_2dBm'
            },
            'digital_marker': 'ON'
        },

        "digital_readout_QRAM": {
            'length': len(QRAM_Exp_TOP2_samples),
            'operation': 'control',
            'waveforms': {
                'single': 'qram_wf'
            },
            'digital_marker': 'ON'
        },

        "digital_readout_CRUS": {
            'length': readout_CRUS_pulse_len,
            'operation': 'control',
            'waveforms': {
                'single': 'zero_wf'
            },
            'digital_marker': 'ON'
        },

        'pulse1_in': {
            'operation': 'control',
            'length': pulse_time,
            'waveforms': {
                'I': 'wf1',
                'Q': 'wf2'
            },
        },

        'pulse_gaus': {
            'operation': 'control',
            'length': len(Gaussian_pulse_samples),
            'waveforms': {
                'I': 'wf_gaus',
                'Q': 'wf_0'
            }
        },

        'pulse1_in': {
            'operation': 'control',
            'length': pulse_time,
            'waveforms': {
                'I': 'wf1',
                'Q': 'wf2'
            },
        },

        'pulse_gaus': {
            'operation': 'control',
            'length': len(Gaussian_pulse_samples),
            'waveforms': {
                'I': 'wf_gaus',
                'Q': 'wf_0'
            }
        }

    },

    'integration_weights': {
        'integW': {
            'cosine': [1.0],
            'sine': [0.0]
        },
        'Homodyne_Det_opt': {
            'cosine': [1.0] * int(len(LO_pulse_samples) / 4),
            'sine': [0.0] * int(len(LO_pulse_samples) / 4)
        },
        'Det_opt': {
            'cosine': [1.0],
            'sine': [0.0]
        }
    },

    "waveforms": {
        'const_wf': {
            'type': 'constant',
            'sample': 0.45
            # 'sample': 0.3
        },
        'const_wf_2dBm': {
            'type': 'constant',
            'sample': 0.4
            # 'sample': 0.3
        },
        'zero_wf': {
            'type': 'constant',
            'sample': 0.0
        },
        'linear_wf': {
            'type': 'arbitrary',
            'samples': Linear
        },
        'exp_wf': {
            'type': 'arbitrary',
            'samples': Exponential
        },
        'Detection_Gaussian_wf': {
            'type': 'arbitrary',
            'samples': Det_Gaussian_samples
        },
        'Detection_Square_wf': {
            'type': 'arbitrary',
            'samples': Det_square_samples
        },
        'Sprint_Gaussian_wf_S': {
            'type': 'arbitrary',
            'samples': Sprint_Exp_Gaussian_samples_S
        },
        'Sprint_Gaussian_wf_N': {
            'type': 'arbitrary',
            'samples': Sprint_Exp_Gaussian_samples_N
        },
        'qram_wf': {
            'type': 'arbitrary',
            'samples': QRAM_Exp_TOP2_samples
        },
        'CRUS_pulser_wf': {
            'type': 'arbitrary',
            'samples': CRUS_pulser_samples
        },
        'CRUS_probe_wf': {
            'type': 'arbitrary',
            'samples': CRUS_probe_samples
        },
        'DC_cal_wf': {
            'type': 'arbitrary',
            'samples': DC_cal_samples
        },
        'LO_pulse_wf': {
            'type': 'arbitrary',
            'samples': LO_pulse_samples
        },
        # 'PolySine_wf': {
        #     'type': 'arbitrary',
        #     'samples': PolySine
        # },
        # 'PolyCosine_wf': {
        #     'type': 'arbitrary',
        #     'samples': PolyCosine
        # },
        'Depump_wf': {
            'type': 'arbitrary',
            'samples': Depump_pulse_samples
        },
        'wf1': {
            'type': 'constant',
            'sample': wf_samples[0]
        },
        'wf2': {
            'type': 'constant',
            'sample': wf_samples[1]
        },
        'wf_gaus': {
            'type': 'arbitrary',
            'samples': Gaussian_pulse_samples
        },
        'wf_0': {
            'type': 'constant',
            'sample': 0
        }
    },

    "digital_waveforms": {
        "ON": {
            "samples": [(1, 0)]
        },
        "OFF": {
            "samples": [(0, 0)]
        },
        "ON_pulse_switch": {
            "samples": [(1, 100000), (0, 0)]
        },
        "OFF_pulse_switch": {
            "samples": [(0, 80), (1, 100000), (0, 0)]
        },
        "Trig": {
            "samples": [(1, 20000), (0, 0)]
        },
        "Trig_ns": {
            "samples": [(1, 40), (0, 0)]
        }
    },

    "mixers": {
        "my_mixer": [
            {"intermediate_frequency": im_freq, "lo_frequency": lo_freq, "correction": correction_matrix}
        ]
    }

}
# add Dig detectors elements without output
config['elements'].update({f"detectors_no_dig_out_{i}":
                     {"singleInput":
                          {"port": (controller, 1)},
                      "digitalOutputs":
                          {"out": (controller, j)},
                      'outputs': {
                          'out1': (controller, 1)
                      },
                      'operations': {
                          'readout': "digital_readout",
                          'readout_SPRINT': "digital_readout_sprint",
                          'readout_QRAM': "digital_readout_QRAM",
                      },
                      'time_of_flight': 36,
                      'smearing': 0,
                      } for i,j in enumerate(dets_number[1:])})

