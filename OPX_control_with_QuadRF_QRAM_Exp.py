# from Config import config
import Config_with_SNSPDs_and_QuadRF_QRAM as Config
from Config_Table import Initial_Values, Phases_Names  # , Values_Factor
from quadRFMOTController import QuadRFMOTController
from quadRFFrequencyScannerController import QuadRFFrequencyScannerController

from qm.QuantumMachinesManager import QuantumMachinesManager
from qm.qua import *
from qm import generate_qua_script
from qm.qua import lib
from scipy import signal
from qm import SimulationConfig
# import matplotlib
# matplotlib.use("Qt5agg")
import matplotlib.pyplot as plt
import numpy as np
import time, json
import sys
import os
import math
import itertools
import operator
from pynput import keyboard
from playsound import playsound
import logging
from logging import StreamHandler, Formatter, INFO, WARN, ERROR
import pymsgbox
import Config_Table
from UtilityResources.HMP4040Control import HMP4040Visa

# --------- Camera functionality ------------##
try:
    from UtilityResources import MvCameraController
    # from OPXcontrol.OPX_control_New_v1 import OPX
    from mvIMPACT import acquire
except:
    print('Warning! could not import MvCamera module; in MvCamera.py')
###  ------------------------------------ ###

## Logging for MW spectroscopy ##
logger = logging.getLogger("MWSpectroscopy")
logger.setLevel(INFO)
handler = StreamHandler(sys.stdout)
handler.setFormatter(Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
logger.addHandler(handler)


###            Support functions          ###
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


def assign_variables_to_element(element, *variables):
    """
    Forces the given variables to be used by the given element thread. Useful as a workaround for when the compiler
    wrongly assigns variables which can cause gaps.
    To be used at the beginning of a program, will add a 16ns wait to the given element. Use an `align()` if needed.

    Example::

        >>> with program() as program_name:
        >>>     a = declare(int)
        >>>     b = declare(fixed)
        >>>     assign_variables_to_element('resonator', a, b)
        >>>     align()
        >>>     ...

    """
    _exp = variables[0]
    for variable in variables[1:]:
        _exp += variable
    wait(4 + 0 * _exp, element)


###  ------------------------------------ ###

all_elements = ["Cooling_Sequence", "MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+", "AntiHelmholtz_Coils", "Measurement"]


def MOT(mot_repetitions):
    """
    The MOT function is used to play the MOT. To that end, we send RF signal to AOM_0, AOM_+, AOM_- for the duration of the MOT.

    Parameters
    ----------
    :param mot_repetitions: derived from the MOT duration, and is calculated in the Experiment parameters section
    """
    FLR = declare(fixed)
    align("Cooling_Sequence", "MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+", "AntiHelmholtz_Coils", "Zeeman_Coils",
          "AOM_2-2/3'", "AOM_2-3'_for_interference", "FLR_detection", "Measurement") # , "Dig_detectors") # , "PULSER_N", "PULSER_S")

    ## MOT build-up ##
    n = declare(int)
    m = declare(int)
    play("Detection" * amp(FLR * 0), "FLR_detection", duration=4)  # we dont know why this works, but Yoav from QM made us write this line to solve an alignment problem we had in the next 2 for loops
    with for_(n, 1, n <= mot_repetitions, n + 1):
        play("MOT_with_EOM" * amp(Config.AOM_0_Attenuation), "MOT_AOM_0")
        play("MOT" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-")
        play("MOT" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+")
        # play("Const_open", "PULSER_N")
        # play("OD_FS" * amp(0.5), "AOM_2-2/3'")
        play("AntiHelmholtz_MOT", "AntiHelmholtz_Coils")
    with for_(m, 1, m <= (mot_repetitions - 1), m + 1):
        measure("Detection", "FLR_detection", None, integration.full("Detection_opt", FLR, "out1"))
        # play("OD_FS" * amp(0.1), "AOM_2-3'_for_interference")

    align("Cooling_Sequence", "MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+", "AntiHelmholtz_Coils", "Zeeman_Coils",
          "AOM_2-2/3'", "FLR_detection", "Measurement") # , "Dig_detectors") #, "PULSER_N", "PULSER_S")

    return FLR

def Pulse_const(total_pulse_duration):
    """
    This pulse id devided to two parts:
        1) Preparation sequence where different parameters changes such as the RF amplitudes to the 0 - + AOMs (scan each separately)
        2) Part where we keep the different values constant for the 0 - + AOMs

    Parameters
    ----------
    :param total_pulse_duration: duration of the whole pulse.
    :param prep_duration: duration of the preparation part of the pulse
    :param initial_amp_0: the initial amplitude of the RF signal to AOM 0 which in turn control the MOT Beams intensity.
    :param initial_amp_minus: the initial amplitude of the RF signal to AOM + which in turn control the MOT Beams intensity.
    :param initial_amp_plus: the initial amplitude of the RF signal to AOM - which in turn control the MOT Beams intensity.
    :param final_amp_0: the final amplitude of the RF signal to AOM 0 which in turn control the MOT Beams intensity.
    :param final_amp_plus: the final amplitude of the RF signal to AOM + which in turn control the MOT Beams intensity.
    :param final_amp_minus: the final amplitude of the RF signal to AOM - which in turn control the MOT Beams intensity.
    :param zero_pulse_duration: derived from both the exact preparation duration and the final relative amplitude of the RF signal to the AOM 0 which in turn control the MOT Beams intensity.
    :param plus_pulse_duration: derived from both the exact preparation duration and the final relative amplitude of the RF signal to the AOM + which in turn control the MOT Beams intensity.
    :param minus_pulse_duration: derived from both the exact preparation duration and the final relative amplitude of the RF signal to the AOM - which in turn control the MOT Beams intensity.
    """

    ## Playing the pulses to the AOMs for the constant part. (Qua) ##
    align("MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+")
    play("Const" * amp(Config.AOM_0_Attenuation), "MOT_AOM_0", duration=total_pulse_duration)
    play("Const" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-", duration=total_pulse_duration)
    play("Const" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+", duration=total_pulse_duration)


def Pulse_with_prep(total_pulse_duration, prep_duration, zero_pulse_duration, plus_pulse_duration,
                    minus_pulse_duration):
    """
    This pulse id devided to two parts:
        1) Preparation sequence where different parameters changes such as the RF amplitudes to the 0 - + AOMs (scan each separately)
        2) Part where we keep the different values constant for the 0 - + AOMs

    Parameters
    ----------
    :param total_pulse_duration: duration of the whole pulse.
    :param prep_duration: duration of the preparation part of the pulse
    :param initial_amp_0: the initial amplitude of the RF signal to AOM 0 which in turn control the MOT Beams intensity.
    :param initial_amp_minus: the initial amplitude of the RF signal to AOM + which in turn control the MOT Beams intensity.
    :param initial_amp_plus: the initial amplitude of the RF signal to AOM - which in turn control the MOT Beams intensity.
    :param final_amp_0: the final amplitude of the RF signal to AOM 0 which in turn control the MOT Beams intensity.
    :param final_amp_plus: the final amplitude of the RF signal to AOM + which in turn control the MOT Beams intensity.
    :param final_amp_minus: the final amplitude of the RF signal to AOM - which in turn control the MOT Beams intensity.
    :param zero_pulse_duration: derived from both the exact preparation duration and the final relative amplitude of the RF signal to the AOM 0 which in turn control the MOT Beams intensity.
    :param plus_pulse_duration: derived from both the exact preparation duration and the final relative amplitude of the RF signal to the AOM + which in turn control the MOT Beams intensity.
    :param minus_pulse_duration: derived from both the exact preparation duration and the final relative amplitude of the RF signal to the AOM - which in turn control the MOT Beams intensity.
    """

    ## Playing the pulses to the AOMs for the preparation sequence. (Qua) ##
    align("MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+")
    # align("MOT_AOM_-", "MOT_AOM_+")
    with if_(prep_duration > 0):
        play("Linear" * amp(Config.AOM_0_Attenuation), "MOT_AOM_0", duration=zero_pulse_duration,
             truncate=prep_duration)
        play("Linear" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-", duration=minus_pulse_duration,
             truncate=prep_duration)
        play("Linear" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+", duration=plus_pulse_duration,
             truncate=prep_duration)

    ## Playing the pulses to the AOMs for the constant part. (Qua) ##
    align("MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+")
    # align("MOT_AOM_-", "MOT_AOM_+")
    with if_(total_pulse_duration > prep_duration):
        play("Const" * amp(Config.AOM_0_Attenuation), "MOT_AOM_0",
             duration=(total_pulse_duration - prep_duration))
        play("Const" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-",
             duration=(total_pulse_duration - prep_duration))
        play("Const" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+",
             duration=(total_pulse_duration - prep_duration))


def Pulse_with_prep_with_chirp(total_pulse_duration, prep_duration, zero_pulse_duration, plus_pulse_duration,
                               minus_pulse_duration, aom_chirp_rate, delta_f):
    """
    This pulse id devided to two parts:
        1) Preparation sequence where different parameters changes such as amplitudes and frequencies of the 0 - + AOMs (scan each separately)
        2) Part where we keep the different values constant for the 0 - + AOMs

    Parameters
    ----------
    :param total_pulse_duration: duration of the whole pulse.
    :param prep_duration: duration of the preparation part of the pulse
    :param initial_amp_0: the initial amplitude of the RF signal to AOM 0 which in turn control the MOT Beams intensity.
    :param initial_amp_minus: the initial amplitude of the RF signal to AOM + which in turn control the MOT Beams intensity.
    :param initial_amp_plus: the initial amplitude of the RF signal to AOM - which in turn control the MOT Beams intensity.
    :param final_amp_0: the final amplitude of the RF signal to AOM 0 which in turn control the MOT Beams intensity.
    :param final_amp_plus: the final amplitude of the RF signal to AOM + which in turn control the MOT Beams intensity.
    :param final_amp_minus: the final amplitude of the RF signal to AOM - which in turn control the MOT Beams intensity.
    :param zero_pulse_duration: derived from both the exact preparation duration and the final relative amplitude of the RF signal to the AOM 0 which in turn control the MOT Beams intensity.
    :param plus_pulse_duration: derived from both the exact preparation duration and the final relative amplitude of the RF signal to the AOM + which in turn control the MOT Beams intensity.
    :param minus_pulse_duration: derived from both the exact preparation duration and the final relative amplitude of the RF signal to the AOM - which in turn control the MOT Beams intensity.
    :param aom_chirp_rate: the slope at which the pulse reach it's final frequency. The final frequency of the RF signals to AOMs + and - which in turn control the relative frame of reference for the atoms.
    """

    ## Playing the pulses to the AOMs for the preparation sequence. (Qua) ##
    align("MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+")
    with if_(prep_duration > 0):
        with if_(zero_pulse_duration == prep_duration):
            play("Const" * amp(Config.AOM_0_Attenuation), "MOT_AOM_0", duration=prep_duration)
            play("Const" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-", duration=prep_duration,
                 chirp=(aom_chirp_rate, "mHz/nsec"))
            play("Const" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+", duration=prep_duration,
                 chirp=(-aom_chirp_rate, "mHz/nsec"))
        with else_():
            play("Linear" * amp(Config.AOM_0_Attenuation), "MOT_AOM_0", duration=zero_pulse_duration,
                 truncate=prep_duration)
            play("Linear" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-", duration=minus_pulse_duration,
                 truncate=prep_duration, chirp=(aom_chirp_rate, "mHz/nsec"))
            play("Linear" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+", duration=plus_pulse_duration,
                 truncate=prep_duration, chirp=(-aom_chirp_rate, "mHz/nsec"))

    ## Playing the pulses to the AOMs for the constant part. (Qua) ##
    align("MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+")
    update_frequency("MOT_AOM_-", Config.IF_AOM_MOT + delta_f)
    update_frequency("MOT_AOM_+", Config.IF_AOM_MOT - delta_f)
    with if_(total_pulse_duration > prep_duration):
        play("Const" * amp(Config.AOM_0_Attenuation), "MOT_AOM_0",
             duration=(total_pulse_duration - prep_duration))
        play("Const" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-",
             duration=(total_pulse_duration - prep_duration))
        play("Const" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+",
             duration=(total_pulse_duration - prep_duration))


def FreeFall(freefall_duration, coils_timing):
    """
    The FreeFall function is used to stop the MOT beams and start "free-fall" section of the experiment.
    The function uses the freefall_duration and coil_timings for the purpose of turning on the Zeeman coils needed for the super-sprint experiment.

    Parameters
    ----------
    :param freefall_duration: The time needed for the atoms to reach the toroids and return
    :param coils_timing: The time to turn the Zeeman coils - "Z_bias"
    """

    # update_frequency("MOT_AOM_0", Config.IF_AOM_MOT)
    update_frequency("MOT_AOM_-", Config.IF_AOM_MOT)
    update_frequency("MOT_AOM_+", Config.IF_AOM_MOT)

    ## Aligning all the different elements used during the freefall time of the experiment ##
    align("Cooling_Sequence", "MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+", "Zeeman_Coils", "AOM_2-2/3'",
          "Measurement", "AOM_Early", "AOM_Late", "PULSER_ANCILLA", "PULSER_N", "PULSER_S") # , "Dig_detectors")

    ## Zeeman Coils turn-on sequence ##
    wait(coils_timing, "Zeeman_Coils")
    # play("ZeemanOFF", "Zeeman_Coils", duration=(freefall_duration - coils_timing))
    play("ZeemanSplit", "Zeeman_Coils", duration=(freefall_duration - coils_timing))


def Measure(measuring_duration):
    """
    The Measure function is an all purpose tool for any measurement is needed.

    Parameters
    ----------
    :param measuring_duration: Given in [msec] and plays a trigger to the measuring component.
    :return:
    """

    play('Measure', "Measurement", duration=measuring_duration)


def Depump_Measure(Depump_pulse_duration, spacing_duration):
    """
    The Measure function is an all purpose tool for any measurement is needed.

    Parameters
    ----------
    :param Depump_pulse_duration: The duration of the OD pulses in [msec] and plays a trigger to the measuring component.
    :param spacing_duration: The spacing duration between pulses in [msec].
    :return:
    """
    update_frequency("AOM_2-2/3'", Config.IF_AOM_Depump)
    play("Depump", "AOM_2-2/3'", duration=Depump_pulse_duration)
    wait(spacing_duration, "AOM_2-2/3'")
    play("Depump", "AOM_2-2/3'", duration=Depump_pulse_duration)
    update_frequency("AOM_2-2/3'", Config.IF_AOM_OD)


def OD_Measure(OD_pulse_duration, spacing_duration, OD_sleep):
    """
    The Measure function is an all purpose tool for any measurement is needed.

    Parameters
    ----------
    :param OD_pulse_duration: The duration of the OD pulses in [msec] and plays a trigger to the measuring component.
    :param spacing_duration: The spacing duration between pulses in [msec].
    :return:
    """
    # update_frequency("AOM_2-2/3'", Config.IF_AOM_OD)
    play("OD_FS", "AOM_2-2/3'", duration=OD_pulse_duration)
    update_frequency("AOM_2-2/3'", Config.IF_AOM_Depump)
    wait(OD_sleep, "AOM_2-2/3'")
    play("Depump", "AOM_2-2/3'", duration=spacing_duration)
    update_frequency("AOM_2-2/3'", Config.IF_AOM_OD)
    wait(OD_sleep, "AOM_2-2/3'")
    play("OD_FS", "AOM_2-2/3'", duration=OD_pulse_duration)


def MZ_balancing(m_time, m_window, shutter_open_time, phase_rep, points_for_sum,
                 counts_st_B, counts_st_D, phase_correction_st, phase_for_min_st, tt_st_B, tt_st_D):

    counts_B = declare(int)
    # tt_vec_B2 = declare(int, size=Config.vec_size)
    counts_B_sum = declare(int, value=0)
    counts_D = declare(int)
    # tt_vec_D2 = declare(int, size=Config.vec_size)
    counts_D_sum = declare(int, value=0)

    assign_variables_to_element("AOM_Early", counts_B_sum)
    assign_variables_to_element("AOM_Late", counts_D_sum)
    t = declare(int)
    i = declare(fixed)
    j = declare(int)
    k = declare(int)

    phase_scan = declare(fixed, value=[1, 0.25])
    # phase_scan = declare(fixed, value=[1.0, 1.0])
    phase_correction = declare(fixed, value=(1/phase_rep))
    phase_correction_per_scan = declare(fixed)
    phase_correction_value = declare(fixed)
    phase_correction_min = declare(fixed, value=0.0)
    phase_correction_start = declare(fixed, value=0.0)
    min_counts_D = declare(int, value=0)
    max_counts_B = declare(int, value=0)

    assign(max_counts_B, 0)
    assign(min_counts_D, 0)
    assign(phase_correction_min, 0.0)
    assign(phase_correction_start, 0.0)

    align("AOM_Early", "AOM_Late", "PULSER_N", "PULSER_S")
    play("Const_open_triggered" * amp(0), "PULSER_N", duration=shutter_open_time)
    play("Const_open" * amp(0), "PULSER_S", duration=shutter_open_time)

    align("AOM_Early", "AOM_Late", "PULSER_N", "PULSER_S")
    # play("Const_open_triggered", "PULSER_N", duration=m_time)
    # play("Const_open", "PULSER_S", duration=m_time)

    with for_each_(i, phase_scan):
        assign(phase_correction_per_scan, phase_correction * i)
        assign(phase_correction_start, phase_correction_min - i * 0.5)
        assign(max_counts_B, 0)
        reset_frame("AOM_Early", "AOM_Late")
        frame_rotation_2pi(phase_correction_start, "AOM_Early")
        with for_(j, 0, j < phase_rep, j + 1):
            assign(counts_B_sum, 0)
            assign(counts_D_sum, 0)
            with for_(k, 1, k <= points_for_sum, k + 1):
                play("MZ_balancing_pulses", "PULSER_N")
                # play("MZ_balancing_pulses", "PULSER_S")
                measure("MZ_balancing_pulses", "AOM_Early", None,
                        counting.digital(counts_B, m_window, element_outputs="OutBright1"))
                        # time_tagging.digital(tt_vec_B2, m_window, element_output="OutBright2", targetLen=counts_B2))

                measure("MZ_balancing_pulses", "AOM_Late", None,
                        counting.digital(counts_D, m_window, element_outputs="OutDark2"))
                        # time_tagging.digital(tt_vec_D2, m_window, element_output="OutDark2", targetLen=counts_D2))

                assign(counts_B_sum, counts_B_sum + counts_B)
                assign(counts_D_sum, counts_D_sum + counts_D)

            assign(phase_correction_value, phase_correction_start + Cast.mul_fixed_by_int(phase_correction_per_scan, j))
            save(phase_correction_value, phase_correction_st)
            save(counts_B_sum, counts_st_B)
            save(counts_D_sum, counts_st_D)

            # with if_(t == 1):
            #     assign(min_counts_D, (counts_D1 + counts_D2_sum))
            # with elif_((counts_D1 + counts_D2_sum) < min_counts_D):
            #     assign(min_counts_D, (counts_D1 + counts_D2_sum))
            #     assign(phase_correction_min, phase_correction_value)
            # assign(min_counts_D, Util.cond((counts_D1 + counts_D2) > min_counts_D,
            #                                 min_counts_D, (counts_D1 + counts_D2)))

            with if_((counts_B_sum - counts_D_sum) > max_counts_B):
                assign(max_counts_B, (counts_B_sum - counts_D_sum))
                assign(phase_correction_min, phase_correction_value)

            save(phase_correction_min, phase_for_min_st)
            frame_rotation_2pi(phase_correction_per_scan, "AOM_Early")

    reset_frame("AOM_Early", "AOM_Late")
    frame_rotation_2pi(phase_correction_min, "AOM_Early")


# def MZ_balancing(m_time, m_window, rep, points_for_sum,
#                  counts_st_B1, counts_st_B2, counts_st_D1, counts_st_D2, phase_correction_st, phase_for_min_st):
#
#     counts_B1 = declare(int)
#     counts_B2 = declare(int)
#     counts_D1 = declare(int)
#     counts_D2 = declare(int)
#
#     t = declare(int)
#
#     last_extinction_ratio = declare(fixed, value=0)
#     phase_correction = declare(fixed, value=(2/rep))
#     phase_correction_value = declare(fixed)
#     phase_correction_min = declare(fixed, value=(2/rep))
#     min_counts_D = declare(int, value=0)
#     max_counts_B = declare(int, value=0)
#
#     assign(max_counts_B, 0)
#     assign(min_counts_D, 0)
#     assign(counts_B1, 0)
#     assign(counts_D1, 0)
#     reset_frame("AOM_Early")
#     align("AOM_Early", "AOM_Late", "PULSER_N", "PULSER_S")
#     play("Const_open", "PULSER_S", duration=m_time)
#     play("Const_open", "PULSER_N", duration=m_time)
#
#     wait(287, "AOM_Early", "AOM_Late")
#     # with for_(t, 0, t < (m_time * 4), t + m_window):
#     with for_(t, 1, t <= rep, t + 1):
#
#         measure("MZ_balancing_pulses", "AOM_Early", None,
#                 # counting.digital(counts_B1, int(len(Config.QRAM_MZ_balance_pulse_Early)), element_outputs="OutBright1"),
#                 counting.digital(counts_B2, int(len(Config.QRAM_MZ_balance_pulse_Early)), element_outputs="OutBright2"))
#         measure("MZ_balancing_pulses", "AOM_Late", None,
#                 # counting.digital(counts_D1, int(len(Config.QRAM_MZ_balance_pulse_Late)), element_outputs="OutDark1"),
#                 counting.digital(counts_D2, int(len(Config.QRAM_MZ_balance_pulse_Late)), element_outputs="OutDark2"))
#
#         assign(phase_correction_value, Cast.mul_fixed_by_int(phase_correction, t))
#         save(phase_correction_value, phase_correction_st)
#
#         save(counts_B1, counts_st_B1)
#         save(counts_B2, counts_st_B2)
#         save(counts_D1, counts_st_D1)
#         save(counts_D2, counts_st_D2)
#
#         # with if_((t > 1) | ((counts_D1 + counts_D2) < min_counts_D)):
#         #     assign(min_counts_D, (counts_D1 + counts_D2))
#         #     assign(phase_correction_min, phase_correction_value)
#         #     # assign(min_counts_D, Util.cond((counts_D1 + counts_D2) > min_counts_D,
#         #     #                                 min_counts_D, (counts_D1 + counts_D2)))
#         # with else_():
#         #     assign(min_counts_D, (counts_D1 + counts_D2))
#         with if_((counts_B1 + counts_B2) > max_counts_B):
#             assign(max_counts_B, (counts_B1 + counts_B2))
#             assign(phase_correction_min, phase_correction_value)
#
#         save(phase_correction_min, phase_for_min_st)
#
#         frame_rotation_2pi(phase_correction, "AOM_Early")
#         # frame_rotation_2pi(0.1, "AOM_Early")
#         #
#         # assign(last_extinction_ratio, Cast.to_fixed((counts_D1 + counts_D2)/(counts_D1 + counts_D2
#         #                                                                                + counts_B1 + counts_B2)))
#         # assign(sign, 1)
#     reset_frame("AOM_Early")
#     frame_rotation_2pi(phase_correction_min, "AOM_Early")

def MZ_balancing_check(m_time, m_window, rep,
                       counts_st_B_balanced, counts_st_D_balanced):

    counts_B = declare(int)
    counts_D = declare(int)

    t = declare(int)

    align("AOM_Early", "AOM_Late", "PULSER_N", "PULSER_S")
    with for_(t, 0, t < rep, t + 1):
        play("MZ_balancing_pulses", "PULSER_N")
        # play("MZ_balancing_pulses", "PULSER_S")
        measure("MZ_balancing_pulses", "AOM_Early", None,
                counting.digital(counts_B, m_window, element_outputs="OutBright2"))
        measure("MZ_balancing_pulses", "AOM_Late", None,
                counting.digital(counts_D, m_window, element_outputs="OutDark2"))
        save(counts_B, counts_st_B_balanced)
        save(counts_D, counts_st_D_balanced)


def QRAM_Exp(m_off_time, m_time, m_window, shutter_open_time,
             ON_counts_st1, ON_counts_st2, ON_counts_st3,
             ON_counts_st4, ON_counts_st5, ON_counts_st6,
             ON_counts_st7, ON_counts_st8,
             tt_st_1, tt_st_2, tt_st_3, tt_st_4, tt_st_5, tt_st_6, tt_st_7, tt_st_8, rep_st,
             balancing_check_window, length,
             rep_MZ_check, counts_st_B_balanced, counts_st_D_balanced):
    """
     Generates train of 8 pulses (to be configured from config if each of thenm is north ore south)
     for SPRINT experiments

       Parameters
       ----------
       :param M_delay: The time delay from end of PGC until the first 2-3' pulse (OD + trigger).
       :param rep: The number of measuring window repetitions, derived from OD measuring duration (M_time/M_window).
       :param m_time: The duration of the entire measurement time.
       :param m_window: The duration of each measuring window - fpr each window there is 28 nsec "deadtime".
       :param counts_st1: The stream array for the number of photons counted at each measuring window for detector 1.
       :param counts_st2: The stream array for the number of photons counted at each measuring window for detector 2.
       :param counts_st3: The stream array for the number of photons counted at each measuring window for detector 3.
       :param counts_st4: The stream array for the number of photons counted at each measuring window for detector 4.
       :param counts_st5: The stream array for the number of photons counted at each measuring window for detector 5.
       :param counts_st6: The stream array for the number of photons counted at each measuring window for detector 6.
       :param counts_st7: The stream array for the number of photons counted at each measuring window for detector 7.
       :param counts_st8: The stream array for the number of photons counted at each measuring window for detector 8.
       """

    vec_size = Config.vec_size

    counts1 = declare(int)
    counts2 = declare(int)
    counts3 = declare(int)
    counts4 = declare(int)
    counts5 = declare(int)
    counts6 = declare(int)
    counts7 = declare(int)
    counts8 = declare(int)

    tt_vec1 = declare(int, size=vec_size)
    tt_vec2 = declare(int, size=vec_size)
    tt_vec3 = declare(int, size=vec_size)
    tt_vec4 = declare(int, size=vec_size)
    tt_vec5 = declare(int, size=vec_size)
    tt_vec6 = declare(int, size=vec_size)
    tt_vec7 = declare(int, size=vec_size)
    tt_vec8 = declare(int, size=vec_size)

    assign_variables_to_element("Dig_detectors", counts1)

    n = declare(int, value=0)
    t = declare(int)
    m = declare(int)

    # assign_variables_to_element("Dig_detectors", tt_vec1[0], counts1, m_window)
    align("AOM_Early", "AOM_Late", "PULSER_ANCILLA", "PULSER_N", "PULSER_S", "Dig_detectors") #, "AOM_2-2'")
    # play("Depump", "AOM_2-2'", duration=shutter_open_time)
    # play("Const_open_triggered" * amp(0), "PULSER_N", duration=shutter_open_time)
    # play("Const_open" * amp(0), "PULSER_S", duration=shutter_open_time)
    # align("AOM_Early", "AOM_Late", "PULSER_ANCILLA", "PULSER_N", "PULSER_S", "Dig_detectors")

    with for_(t, 0, t < (m_time + m_off_time) * 4, t + int(len(Config.QRAM_Exp_Gaussian_samples_S))): #assaf comment debbuging
        play("QRAM_experiment_pulses_Ancilla", "PULSER_ANCILLA")
        play("QRAM_experiment_pulses_S", "PULSER_S")
        play("QRAM_experiment_pulses_N", "PULSER_N")
        play("QRAM_experiment_pulses_Early", "AOM_Early")
        play("QRAM_experiment_pulses_Late", "AOM_Late")

    # wait(297, "Dig_detectors")
    wait(300, "Dig_detectors")
    # with for_(n, 0, n < m_time * 4, n + m_window):
    measure("readout_QRAM", "Dig_detectors", None,
            time_tagging.digital(tt_vec1, m_window, element_output="out1", targetLen=counts1),
            time_tagging.digital(tt_vec2, m_window, element_output="out2", targetLen=counts2),
            time_tagging.digital(tt_vec3, m_window, element_output="out3", targetLen=counts3),
            time_tagging.digital(tt_vec4, m_window, element_output="out4", targetLen=counts4),
            time_tagging.digital(tt_vec5, m_window, element_output="out5", targetLen=counts5),
            time_tagging.digital(tt_vec6, m_window, element_output="out6", targetLen=counts6),
            time_tagging.digital(tt_vec7, m_window, element_output="out7", targetLen=counts7),
            time_tagging.digital(tt_vec8, m_window, element_output="out8", targetLen=counts8),
            )

    MZ_balancing_check(balancing_check_window, length, rep_MZ_check, counts_st_B_balanced, counts_st_D_balanced)

    ## Save Data: ##

    # Number of Photons (NOP) Count stream for each detector: ##
    save(counts1, ON_counts_st1)
    save(counts2, ON_counts_st2)
    save(counts3, ON_counts_st3)
    save(counts4, ON_counts_st4)
    save(counts5, ON_counts_st5)
    save(counts6, ON_counts_st6)
    save(counts7, ON_counts_st7)
    save(counts8, ON_counts_st8)

    save(counts1, tt_st_1)
    save(counts2, tt_st_2)
    save(counts3, tt_st_3)
    save(counts4, tt_st_4)
    save(counts5, tt_st_5)
    save(counts6, tt_st_6)
    save(counts7, tt_st_7)
    save(counts8, tt_st_8)

    with for_(m, 0, m < vec_size, m + 1):
        wait(1000)
        save(tt_vec1[m], tt_st_1)
        save(tt_vec2[m], tt_st_2)
        save(tt_vec3[m], tt_st_3)
        save(tt_vec4[m], tt_st_4)
        save(tt_vec5[m], tt_st_5)
        save(tt_vec6[m], tt_st_6)
        save(tt_vec7[m], tt_st_7)
        save(tt_vec8[m], tt_st_8)
        save(n, rep_st)


def opx_control(obj, qm):
    with program() as opx_control_prog:
        ## declaring program variables: ##
        i = declare(int)
        Trigger_Phase = declare(int, value=Phases_Names.index(obj.Exp_Values['Triggering_Phase']))
        Imaging_Phase = declare(int, value=Phases_Names.index(obj.Exp_Values['Imaging_Phase']))
        x = declare(int)

        # Boolean variables:
        AntiHelmholtz_ON = declare(bool, value=True)
        QRAM_Exp_ON = declare(bool, value=False)

        # MOT variables
        MOT_Repetitions = declare(int, value=obj.Exp_Values['MOT_rep'])
        post_MOT_delay = declare(int, value=int(obj.Exp_Values['Post_MOT_delay'] * 1e6 / 4))

        # AntiHelmholtz delay after MOT:
        antihelmholtz_delay = declare(int, value=int(obj.Exp_Values['AntiHelmholtz_delay'] * 1e6 / 4))

        # PGC variables:
        pgc_duration = declare(int, value=int(obj.pgc_duration))

        # Fountain variables:
        fountain_duration = declare(int, value=int(obj.fountain_duration))
        # TODO: pre_PGC_fountain_duration = declare(int, value=int(obj.pre_PGC_fountain_duration))
        fountain_pulse_duration_0 = declare(int, value=int(obj.fountain_pulse_duration_0))  # The relative duration to reach the desired amplitude
        fountain_pulse_duration_minus = declare(int, value=int(obj.fountain_pulse_duration_minus))  # The relative duration to reach the desired amplitude
        fountain_pulse_duration_plus = declare(int, value=int(obj.fountain_pulse_duration_plus))  # The relative duration to reach the desired amplitude
        fountain_aom_chirp_rate = declare(int, value=int(obj.fountain_aom_chirp_rate))  # mHz/nsec
        fountain_delta_f = declare(int, value=int(obj.Exp_Values['Fountain_final_Delta_freq']))

        # Free fall duration:
        FreeFall_duration = declare(int, value=int(obj.Exp_Values['FreeFall_duration'] * 1e6 / 4))
        coils_timing = declare(int, value=int(obj.Exp_Values['Coil_timing'] * 1e6 / 4))

        # Imaging variables:
        N_Snaps = declare(int, value=obj.Exp_Values['N_Snaps'])
        Buffer_Cycles = declare(int, value=obj.Exp_Values['Buffer_Cycles'])

        # General measurements variables:
        PrePulse_duration = declare(int, value=int(obj.Exp_Values['PrePulse_duration'] * 1e6 / 4))
        Balancing_check_window = declare(int, value=int(obj.Balancing_check_window * 1e6 / 4))
        # Pulse_1_parameters:
        Pulse_1_duration = declare(int, value=int(obj.Exp_Values['Pulse_1_duration'] * 1e6 / 4))
        Pulse_1_decay_time = declare(int, value=int(obj.Exp_Values['Pulse_1_decay_duration'] * 1e6 / 4))
        pulse_1_duration_0 = declare(int, value=int((obj.Exp_Values['Pulse_1_initial_amp_0'] * obj.Exp_Values[
            'Pulse_1_decay_duration'] / (obj.Exp_Values['Pulse_1_initial_amp_0'] - obj.Exp_Values[
            'Pulse_1_final_amp_0'])) * 1e6 / 4))  # The relative duration to reach the desired amplitude
        pulse_1_duration_minus = declare(int, value=int((obj.Exp_Values['Pulse_1_initial_amp_minus'] * obj.Exp_Values[
            'Pulse_1_decay_duration'] / (obj.Exp_Values['Pulse_1_initial_amp_minus'] - obj.Exp_Values[
            'Pulse_1_final_amp_minus'])) * 1e6 / 4))  # The relative duration to reach the desired amplitude
        pulse_1_duration_plus = declare(int, value=int((obj.Exp_Values['Pulse_1_initial_amp_plus'] * obj.Exp_Values[
            'Pulse_1_decay_duration'] / (obj.Exp_Values['Pulse_1_initial_amp_plus'] - obj.Exp_Values[
            'Pulse_1_final_amp_plus'])) * 1e6 / 4))  # The relative duration to reach the desired amplitude

        # Depump measurement variables:
        Depump_pulse_duration = declare(int, value=int(obj.Depump_pulse_duration * 1e6 / 4))
        Depump_pulses_spacing = declare(int, value=int(obj.Depump_pulses_spacing * 1e6 / 4))
        Depump_start = declare(int, value=int(obj.Depump_Start * 1e6 / 4))

        # OD measurement variables:
        ## In fiber:
        M_off_time = declare(int, value=int(obj.M_off_time / 4))
        shutter_open_time = declare(int, value=int(obj.Shutter_open_time * 1e6 / 4))

        # Stream processing:
        counts_st_B = declare_stream()
        counts_st_B_balanced = declare_stream()
        tt_st_B = declare_stream()
        counts_st_D = declare_stream()
        counts_st_D_balanced = declare_stream()
        tt_st_D = declare_stream()
        phase_correction_st = declare_stream()
        phase_for_min_st = declare_stream()
        ON_counts_st1 = declare_stream()
        ON_counts_st2 = declare_stream()
        ON_counts_st3 = declare_stream()
        ON_counts_st4 = declare_stream()
        ON_counts_st5 = declare_stream()
        ON_counts_st6 = declare_stream()
        ON_counts_st7 = declare_stream()
        ON_counts_st8 = declare_stream()
        tt_st_1 = declare_stream()
        tt_st_2 = declare_stream()
        tt_st_3 = declare_stream()
        tt_st_4 = declare_stream()
        tt_st_5 = declare_stream()
        tt_st_6 = declare_stream()
        tt_st_7 = declare_stream()
        tt_st_8 = declare_stream()
        rep_st = declare_stream()
        AntiHelmholtz_ON_st = declare_stream()
        FLR_st = declare_stream()
        assign(IO1, 0)
        assign(IO2, 0)

        wait(100)

        with infinite_loop_():
            assign(i, IO1)

            ##########################
            ## Cooling Sequence ##
            ##########################

            # MOT sequence:

            FLR = MOT(MOT_Repetitions)
            play("AntiHelmholtz_MOT", "AntiHelmholtz_Coils", duration=antihelmholtz_delay)

            # Delay before fountain:
            wait(post_MOT_delay, "Cooling_Sequence")
            align(*all_elements)

            # Fountain sequence:
            wait(fountain_duration, "Cooling_Sequence")
            Pulse_with_prep_with_chirp(fountain_duration, obj.fountain_prep_duration,
                                       fountain_pulse_duration_0, fountain_pulse_duration_minus,
                                       fountain_pulse_duration_plus, fountain_aom_chirp_rate,
                                       fountain_delta_f)

            # PGC sequence:
            wait(pgc_duration, "Cooling_Sequence")
            Pulse_const(pgc_duration)
            align(*all_elements)

            # FreeFall sequence:
            with if_(QRAM_Exp_ON):
                assign(x, (30678780 - 3106 + 656000 * 2 + 4) // 4) # TODO -  added 38688900 to fix new delay due to wait(1000) in saving sprint data with vector size 10000, should be fixed as well
            with else_():
                assign(x, - 3106)
            # FreeFall(FreeFall_duration, coils_timing)
            FreeFall(FreeFall_duration - x, coils_timing)

            ##########################
            ## Measurement Sequence ##
            ##########################

            with if_(Trigger_Phase == 3):  # when trigger on PrePulse
                ## Trigger QuadRF Sequence #####################
                play("C_Seq", "Cooling_Sequence", duration=2500)
                ################################################

            align("Cooling_Sequence", "AOM_Early", "AOM_Late")
            # wait(2500, "AOM_2-2'")
            with if_(QRAM_Exp_ON):
                wait(PrePulse_duration, "Cooling_Sequence")
                # play("Depump", "AOM_2-2'", duration=(PrePulse_duration - shutter_open_time))
                MZ_balancing((PrePulse_duration - shutter_open_time - Balancing_check_window), len(Config.QRAM_MZ_balance_pulse_Late),
                             shutter_open_time, obj.phase_rep_MZ, obj.points_for_sum,
                             counts_st_B, counts_st_D, phase_correction_st, phase_for_min_st, tt_st_B, tt_st_D)
                MZ_balancing_check(Balancing_check_window, len(Config.QRAM_MZ_balance_pulse_Late),
                                   obj.rep_MZ_check, counts_st_B_balanced, counts_st_D_balanced)
            with else_():
                # play("Depump", "AOM_2-2'", duration=PrePulse_duration)
                wait(PrePulse_duration, "Cooling_Sequence")
            align(*all_elements, "AOM_2-2/3'", "AOM_Early", "AOM_Late", "PULSER_ANCILLA", "PULSER_N", "PULSER_S", "Dig_detectors")

            with if_(Trigger_Phase == 4):  # when trigger on pulse 1
                ## Trigger QuadRF Sequence #####################
                play("C_Seq", "Cooling_Sequence", duration=2500)
                ################################################
            ## For taking an image:
            with if_((Pulse_1_duration > 0) & ~QRAM_Exp_ON):
                align(*all_elements)
                Pulse_with_prep(Pulse_1_duration, Pulse_1_decay_time, pulse_1_duration_0,
                                pulse_1_duration_minus, pulse_1_duration_plus)
                Measure(Pulse_1_duration)  # This triggers camera (Control 7)
                align(*all_elements)

            with if_((Pulse_1_duration > 0) & QRAM_Exp_ON):
                align("Dig_detectors", "AOM_Early", "AOM_Late", "PULSER_ANCILLA", "PULSER_N", "PULSER_S")
                QRAM_Exp(M_off_time, Pulse_1_duration, obj.M_window, shutter_open_time,
                         ON_counts_st1, ON_counts_st2, ON_counts_st3,
                         ON_counts_st4, ON_counts_st5, ON_counts_st6,
                         ON_counts_st7, ON_counts_st8,
                         tt_st_1, tt_st_2, tt_st_3, tt_st_4, tt_st_5, tt_st_6, tt_st_7, tt_st_8, rep_st,
                         Balancing_check_window, len(Config.QRAM_MZ_balance_pulse_Late),
                         obj.rep_MZ_check, counts_st_B_balanced, counts_st_D_balanced)
                align("Dig_detectors", "AOM_Early", "AOM_Late", "PULSER_ANCILLA", "PULSER_N", "PULSER_S")

            save(AntiHelmholtz_ON, AntiHelmholtz_ON_st)
            with if_(AntiHelmholtz_ON):
                save(FLR, FLR_st)

            assign(N_Snaps, 1)
            assign(Buffer_Cycles, 0)
            assign(i, IO1)

            ## PARAMETERS UPDATE ##
            with if_(i > 0):
                pause()
            with while_(i > 0):
                ## Boolean variables control: ##
                with if_(i == 1):
                    update_frequency("MOT_AOM_0", Config.IF_AOM_MOT)
                with if_(i == 2):
                    update_frequency("MOT_AOM_0", Config.IF_AOM_MOT_OFF)
                with if_(i == 4):
                    assign(QRAM_Exp_ON, IO2)
                # ## AntiHelmholtz control ##
                # with if_(i == 10):
                #     assign(antihelmholtz_delay, IO2)
                # ## MOT variables control ##
                # with if_(i == 11):  # Live control over the
                #     assign(MOT_Repetitions, IO2)
                # with if_(i == 12):  # Live control over the
                #     assign(post_MOT_delay, IO2)

                # ## PGC variables control ##
                # with if_(i == 20):  # Live control over the PGC duration
                #     assign(pgc_duration, IO2)
                #
                # ## Fountain variables control: ##
                # with if_(i == 31):  # Live control over the fountain duration
                #     assign(fountain_duration, IO2)
                # with if_(i == 36):  # Live control over the final amplitude of the fountain AOM 0
                #     assign(fountain_pulse_duration_0, IO2)
                # with if_(i == 37):  # Live control over the final amplitude of the fountain AOM -
                #     assign(fountain_pulse_duration_minus, IO2)
                # with if_(i == 38):  # Live control over the final amplitude of the fountain AOM +
                #     assign(fountain_pulse_duration_plus, IO2)
                # with if_(i == 39):  # Live control over the fountain final frequency of the MOT + and - AOMs
                #     assign(fountain_aom_chirp_rate, IO2)
                #
                # ## Measurement variables control: ##
                with if_(i == 42):
                    assign(PrePulse_duration, IO2)
                with if_(i == 43):
                    assign(Pulse_1_duration, IO2)
                # with if_(i == 44):
                #     assign(Pulse_1_decay_time, IO2)
                # with if_(i == 45):  # The number of frames(snaps) to take for a video with growing snap time intervals
                #     assign(N_Snaps, IO2)
                # with if_(i == 46):
                #     assign(Buffer_Cycles, IO2)
                #
                # ## OD and N_atoms measuring variable control ##
                # with if_(i == 56):  # Live control of the Depump measurement start time
                #     assign(Depump_start, IO2)
                # with if_(i == 57):  # Live control of the Depump measurement pulses duration
                #     assign(Depump_pulse_duration, IO2)
                # with if_(i == 58):  # Live control of the Depump measurement wait duration between the 2 pulses
                #     assign(Depump_pulses_spacing, IO2)
                # with if_(i == 59):  # Live control of the delay due to shutter opening time.
                #     assign(shutter_open_time, IO2)

                pause()
                assign(i, IO1)

        with stream_processing():
            # counts_st_B.buffer(obj.total_phase_rep_MZ).save('Bright_Port_Counts')
            # counts_st_B_balanced.save('Bright_Port_Counts_check')
            counts_st_B.buffer(obj.total_phase_rep_MZ).zip(counts_st_B_balanced.buffer(obj.rep_MZ_check*2)).save('Bright_Port_Counts')
            # counts_st_D.buffer(obj.total_phase_rep_MZ).save('Dark_Port_Counts')
            # counts_st_D_balanced.save('Dark_Port_Counts_check')
            counts_st_D.buffer(obj.total_phase_rep_MZ).zip(counts_st_D_balanced.buffer(obj.rep_MZ_check*2)).save('Dark_Port_Counts')
            phase_correction_st.buffer(obj.total_phase_rep_MZ).zip(phase_for_min_st.buffer(obj.total_phase_rep_MZ)).save('Phase_Correction_array')
            ON_counts_st1.buffer(obj.rep).save('Det1_Counts')
            ON_counts_st2.buffer(obj.rep).save('Det2_Counts')
            ON_counts_st3.buffer(obj.rep).save('Det3_Counts')
            ON_counts_st4.buffer(obj.rep).save('Det4_Counts')
            ON_counts_st5.buffer(obj.rep).save('Det5_Counts')
            ON_counts_st6.buffer(obj.rep).save('Det6_Counts')
            ON_counts_st7.buffer(obj.rep).save('Det7_Counts')
            ON_counts_st8.buffer(obj.rep).save('Det8_Counts')
            (tt_st_1 + rep_st).buffer(obj.vec_size * obj.rep + 1).save('Det1_Probe_TT')
            (tt_st_2 + rep_st).buffer(obj.vec_size * obj.rep + 1).save('Det2_Probe_TT')
            (tt_st_3 + rep_st).buffer(obj.vec_size * obj.rep + 1).save('Det3_Probe_TT')
            (tt_st_4 + rep_st).buffer(obj.vec_size * obj.rep + 1).save('Det4_Probe_TT')
            (tt_st_5 + rep_st).buffer(obj.vec_size * obj.rep + 1).save('Det5_Probe_TT')
            (tt_st_6 + rep_st).buffer(obj.vec_size * obj.rep + 1).save('Det6_Probe_TT')
            (tt_st_7 + rep_st).buffer(obj.vec_size * obj.rep + 1).save('Det7_Probe_TT')
            (tt_st_8 + rep_st).buffer(obj.vec_size * obj.rep + 1).save('Det8_Probe_TT')
            FLR_st.save('FLR_measure')
            AntiHelmholtz_ON_st.save("antihelmholtz_on")

    # sourceFile = open('serialized_code.py', 'w')
    # print(generate_qua_script(opx_control_prog, Config.config), file=sourceFile)
    # sourceFile.close()

    job = qm.execute(opx_control_prog, flags=['auto-element-thread'])


    return job


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def printError(s=''): print(f"{bcolors.WARNING}{s}{bcolors.ENDC}")


def printGreen(s=''): print(f"{bcolors.OKGREEN}{s}{bcolors.ENDC}")


class OPX:
    def __init__(self, config=Config.config):

        ##########################
        # EXPERIMENT PARAMETERS: #
        ##########################

        # ---------- Handle QuadRF ------------
        self.Exp_Values = Initial_Values  # Initialize experiment values to be as in Config_Table.py
        self.QuadRFControllers = []
        # Note: So as not to connect again and again to QuadRF each time we update table, we now save the MOGDevic (actual QuadRF device) connected,
        # we hold this connection until update is finished, the we close the connection.
        # we do still hold the QuadRFController objects, for access to the table (read only!) when the experiment is running.
        qrfContr = QuadRFMOTController(initialValues=self.Exp_Values, updateChannels=(1, 2, 4), topticaLockWhenUpdating=False,
                                       debugging=False, continuous=False)
        self.QuadRFControllers.append(qrfContr)  # updates values on QuadRF (uploads table)
        self.QuadRFControllers.append(QuadRFMOTController(MOGdevice=qrfContr.dev, initialValues={'Operation_Mode': 'Continuous', 'CH3_freq': '90MHz', 'CH3_amp': '31dbm'},
                                                          updateChannels=[3], debugging=False, continuous=False))  # updates values on QuadRF (uploads table)
        #self.QuadRFControllers.append(QuadRFFrequencyScannerController(MOGdevice = qrfContr.dev, channel=2, debugging=False))  # updates values on QuadRF (uploads table)

        self.Update_QuadRF_channels = set({})  # Only update these channels on QuadRF when UpdateParameters method is called [note: this is a python set]
        qrfContr.disconnectQuadRF()
        # ---------- Finish handle QuadRF ------------

        # Free fall variables:
        self.FreeFall_duration = int(self.Exp_Values['FreeFall_duration'] * 1e6 / 4)
        self.Coils_timing = int(self.Exp_Values['Coil_timing'] * 1e6 / 4)
        if (self.Exp_Values['FreeFall_duration'] - self.Exp_Values['Coil_timing']) > 60:
            raise ValueError("FreeFall_duration - Coils_timing can't be larger than 60 ms")

        # PGC variables:
        self.pgc_duration = int(self.Exp_Values['PGC_duration'] * 1e6 / 4)
        self.pgc_prep_duration = int(self.Exp_Values['PGC_prep_duration'] * 1e6 / 4)
        if self.Exp_Values['PGC_initial_amp_0'] == self.Exp_Values['PGC_final_amp_0']:
            self.pgc_pulse_duration_0 = self.pgc_prep_duration
            self.pgc_pulse_duration_minus = self.pgc_prep_duration
            self.pgc_pulse_duration_plus = self.pgc_prep_duration
        else:
            self.pgc_pulse_duration_0 = int((self.Exp_Values['PGC_initial_amp_0'] * self.Exp_Values[
                'PGC_prep_duration'] / (self.Exp_Values['PGC_initial_amp_0'] - self.Exp_Values[
                'PGC_final_amp_0'])) * 1e6 / 4)  # The relative duration to reach the desired amplitude
            self.pgc_pulse_duration_minus = int((self.Exp_Values['PGC_initial_amp_minus'] * self.Exp_Values[
                'PGC_prep_duration'] / (self.Exp_Values['PGC_initial_amp_minus'] - self.Exp_Values[
                'PGC_final_amp_minus'])) * 1e6 / 4)  # The relative duration to reach the desired amplitude
            self.pgc_pulse_duration_plus = int((self.Exp_Values['PGC_initial_amp_plus'] * self.Exp_Values[
                'PGC_prep_duration'] / (self.Exp_Values['PGC_initial_amp_plus'] - self.Exp_Values[
                'PGC_final_amp_plus'])) * 1e6 / 4)  # The relative duration to reach the desired amplitude
            if self.pgc_pulse_duration_0 > int(60 * 1e6 / 4):  # longer then 60ms
                printError(
                    'The values for PGC_initial_amp_0 and PGC_final_amp_0 are too close or PGC_prep_duration is too long, might cause an ERROR!!!')
            if self.pgc_pulse_duration_minus > int(60 * 1e6 / 4):  # longer then 60ms
                printError(
                    'The values for PGC_initial_amp_minus and PGC_final_amp_minus are too close or PGC_prep_duration is too long, might cause an ERROR!!!')
            if self.pgc_pulse_duration_plus > int(60 * 1e6 / 4):  # longer then 60ms
                printError(
                    'The values for PGC_initial_amp_plus and PGC_final_amp_plus are too close or PGC_prep_duration is too long, might cause an ERROR!!!')
        self.pgc_initial_amp_0 = self.Exp_Values['PGC_initial_amp_0']
        self.pgc_initial_amp_minus = self.Exp_Values['PGC_initial_amp_minus']
        self.pgc_initial_amp_plus = self.Exp_Values['PGC_initial_amp_plus']
        self.pgc_final_amp_0 = self.Exp_Values['PGC_final_amp_0']
        self.pgc_final_amp_minus = self.Exp_Values['PGC_final_amp_minus']
        self.pgc_final_amp_plus = self.Exp_Values['PGC_final_amp_plus']
        # self.pgc_aom_chirp_rate = int(self.Exp_Values['PGC_final_Delta_freq'] * 1e3 / (self.Exp_Values[
        #                                                                                    'PGC_prep_duration'] * 1e6))  # [mHz/nsec], If needed pgc preparation duration must be constant!!!

        # Fountain variables:
        self.fountain_duration = int(self.Exp_Values['Fountain_duration'] * 1e6 / 4)
        self.fountain_prep_duration = int(self.Exp_Values['Fountain_prep_duration'] * 1e6 / 4)
        if self.Exp_Values['Fountain_initial_amp_0'] == self.Exp_Values['Fountain_final_amp_0']:
            self.fountain_pulse_duration_0 = self.fountain_prep_duration
            self.fountain_pulse_duration_minus = self.fountain_prep_duration
            self.fountain_pulse_duration_plus = self.fountain_prep_duration
        else:
            self.fountain_pulse_duration_0 = int(
                self.Exp_Values['Fountain_initial_amp_0'] * self.Exp_Values['Fountain_prep_duration'] * 1e6 / 4 / (
                        self.Exp_Values['Fountain_initial_amp_0'] - self.Exp_Values[
                    'Fountain_final_amp_0']))  # The relative duration to reach the desired amplitude
            self.fountain_pulse_duration_minus = int(
                self.Exp_Values['Fountain_initial_amp_minus'] * self.Exp_Values['Fountain_prep_duration'] * 1e6 / 4 / (
                        self.Exp_Values['Fountain_initial_amp_minus'] - self.Exp_Values[
                    'Fountain_final_amp_minus']))  # The relative duration to reach the desired amplitude
            self.fountain_pulse_duration_plus = int(
                self.Exp_Values['Fountain_initial_amp_plus'] * self.Exp_Values['Fountain_prep_duration'] * 1e6 / 4 / (
                        self.Exp_Values['Fountain_initial_amp_plus'] - self.Exp_Values[
                    'Fountain_final_amp_plus']))  # The relative duration to reach the desired amplitude
            if self.fountain_pulse_duration_0 > int(60 * 1e6 / 4):  # longer then 60ms
                printError(
                    'The values for Fountain_initial_amp_0 and Fountain_final_amp_0 are too close or Fountain_prep_duration is too long, might cause an ERROR!!!')
            if self.fountain_pulse_duration_minus > int(60 * 1e6 / 4):  # longer then 60ms
                printError(
                    'The values for Fountain_initial_amp_minus and Fountain_final_amp_minus are too close or Fountain_prep_duration is too long, might cause an ERROR!!!')
            if self.fountain_pulse_duration_plus > int(60 * 1e6 / 4):  # longer then 60ms
                printError(
                    'The values for Fountain_initial_amp_plus and Fountain_final_amp_plus are too close or Fountain_prep_duration is too long, might cause an ERROR!!!')
        self.fountain_initial_amp_0 = self.Exp_Values['Fountain_initial_amp_0']
        self.fountain_initial_amp_minus = self.Exp_Values['Fountain_initial_amp_minus']
        self.fountain_initial_amp_plus = self.Exp_Values['Fountain_initial_amp_plus']
        self.fountain_final_amp_0 = self.Exp_Values['Fountain_final_amp_0']
        self.fountain_final_amp_minus = self.Exp_Values['Fountain_final_amp_minus']
        self.fountain_final_amp_plus = self.Exp_Values['Fountain_final_amp_plus']
        self.fountain_aom_chirp_rate = int(self.Exp_Values['Fountain_final_Delta_freq'] * 1e3 / (
                    self.Exp_Values['Fountain_prep_duration'] * 1e6))  # mHz/nsec

        # OD and Depump measurement parameters:
        self.Depump_pulse_duration = self.Exp_Values['Depump_pulse_duration']  # [msec]
        self.Depump_pulses_spacing = self.Exp_Values['Depump_pulses_spacing']  # [msec]
        self.Depump_Start = self.Exp_Values['Depump_Start']  # [msec]
        ## OD Free space:
        self.OD_FS_pulse_duration = self.Exp_Values['OD_FS_pulse_duration']  # [msec]
        self.OD_FS_pulses_spacing = self.Exp_Values['OD_FS_pulses_spacing']  # [msec]
        self.OD_FS_Start = self.Exp_Values['OD_FS_Start']  # [msec]
        self.OD_FS_sleep = self.Exp_Values['OD_FS_sleep']
        ## OD In-fiber/Transits:
        self.Transit_switch = False
        self.prepulse_duration = self.Exp_Values['PrePulse_duration']
        self.OD_delay = self.Exp_Values['OD_delay']  # [msec]
        self.M_window = self.Exp_Values['M_window']  # [nsec]
        self.OD_duration_pulse1 = self.Exp_Values['OD_duration_pulse1']  # [msec]
        self.OD_sleep = self.Exp_Values['OD_sleep']  # [msec]
        self.OD_duration_pulse2 = self.Exp_Values['OD_duration_pulse2']  # [msec]
        self.M_time = int(self.Exp_Values['M_time'] * 1e6)  # [nsec]
        self.Shutter_open_time = self.Exp_Values['Shutter_open_time']  # [msec]
        self.M_off_time = int(self.Exp_Values['M_off_time'] * 1e6)  # [nsec]
        # self.rep = int(self.M_time / (self.M_window + 28 + 170))
        self.rep = int(self.M_time / self.M_window)
        self.vec_size = Config.vec_size

        # MZ balancing:
        self.Balancing_check_window = 2 # [msec]
        self.rep_MZ_fast_scan = int(0.5 * (self.prepulse_duration - self.Shutter_open_time - self.Balancing_check_window)
                                    * 1e6 / len(Config.QRAM_MZ_balance_pulse_Late))
        self.rep_MZ_slow_scan = int(0.5 * (self.prepulse_duration - self.Shutter_open_time - self.Balancing_check_window)
                                    * 1e6 / len(Config.QRAM_MZ_balance_pulse_Late))
        self.points_for_sum = 4
        if self.rep_MZ_fast_scan == 0:
            self.phase_rep_MZ = 1
        else:
            self.phase_rep_MZ = int(self.rep_MZ_fast_scan/self.points_for_sum)
        self.total_phase_rep_MZ = int(2 * self.rep_MZ_fast_scan/self.points_for_sum)
        self.rep_MZ_check = int(self.Balancing_check_window * 1e6 / len(Config.QRAM_MZ_balance_pulse_North))

        # MW spectroscopy parameters:
        self.MW_start_frequency = int(100e6)  # [Hz]
        self.Pulse_Length_MW = 400  # [usec]
        self.Pulse_Length_OD = 20  # [usec]

        # Main Experiment:
        self.TOP2_pulse_len = int(Config.Probe_pulse_len / 4)  # [nsec]
        self.Calibration_time = 10  # [msec]

        # run daily experiment
        self.Stop_run_daily_experiment = False


        try:
            self.qmm = QuantumMachinesManager(host='132.77.54.230', port='80')
            self.qmm.clear_all_job_results()
            self.qmm.reset_data_processing()
            self.qmm.close_all_quantum_machines()
            self.qm = self.qmm.open_qm(config)
            self.job = opx_control(self, self.qm)

        except KeyboardInterrupt:
            self.job.halt()
            self.qmm.reset_data_processing()

        self.io1_list = []
        self.io2_list = []

        self.keyPress = None  # initing key listener

    def __del__(self):
        # self.job.halt()
        self.qm.close()
        self.qmm.close()

    def updateValue(self, key, value, update_parameters=False):
        if key in self.Exp_Values:
            self.Exp_Values[key] = value
            # Special cases
            if key == 'Operation_Mode':
                if value not in Config_Table.Operation_Modes: print(
                    '\033[91m' + f'Warning: {key} is not a known operation mode!' + '\033[0m')
                for k in Config_Table.Operation_Modes[value]: self.updateValue(k,
                                                                               Config_Table.Operation_Modes[value][k])
            elif key == 'MOT_duration':  # This is a patch. It should work, but still.
                MOT_rep = Config_Table.updateMOT_rep()
                self.updateValue('MOT_rep', MOT_rep)
            elif key == 'OD_FS_pulse_duration':
                self.updateValue('Pulse_1_duration',
                                 self.OD_FS_Start + self.OD_FS_pulses_spacing + 2 * value + 2 * self.OD_FS_sleep)
                self.OD_FS_pulse_duration = value
            elif key == 'OD_FS_pulses_spacing':
                self.updateValue('Pulse_1_duration',
                                 self.OD_FS_Start + value + 2 * self.OD_FS_pulse_duration + 2 * self.OD_FS_sleep)
                self.OD_FS_pulses_spacing = value
            elif key == 'OD_FS_Start':
                self.updateValue('Pulse_1_duration',
                                 value + self.OD_FS_pulses_spacing + 2 * self.OD_FS_pulse_duration + 2 * self.OD_FS_sleep)
                self.OD_FS_Start = value
            elif key == 'Depump_pulse_duration':
                self.updateValue('Pulse_1_duration', self.Depump_Start + self.Depump_pulses_spacing + 2 * value)
                self.Depump_pulse_duration = value
            elif key == 'Depump_pulses_spacing':
                self.updateValue('Pulse_1_duration', self.Depump_Start + value + 2 * self.Depump_pulse_duration)
                self.Depump_pulses_spacing = value
            elif key == 'Depump_Start':
                self.updateValue('Pulse_1_duration',
                                 value + self.Depump_pulses_spacing + 2 * self.Depump_pulse_duration)
                self.Depump_Start = value
            # Keep track of which channels are updated
            if key in Config_Table.Key_to_Channel:
                for ch in Config_Table.Key_to_Channel[key]: self.Update_QuadRF_channels.add(ch)
            else:
                print(
                    '\033[93m' + 'Warning: ' + key + ' is not in Key_to_Channel. What channel does this key belong to? [change in Config_Table.py]' + '\033[0m')

        # See whether OPX should also be updated
        if key in Config_Table.IOParametersMapping:
            io1 = Config_Table.IOParametersMapping[key]
            io2 = value
            if key in self.Values_Factor:
                valueType = self.Values_Factor[key][0]
                if len(self.Values_Factor[key]) == 2:
                    valueFactor = self.Values_Factor[key][1]
                    io2 = valueType(value * valueFactor)
                    self.update_io_parameter(io1, io2)
                elif len(self.Values_Factor[key]) == 3:  # If a function is attached to value change in Config_Table
                    valueFunction = self.Values_Factor[key][2]
                    valueFunction(self, value)  # Call the predefined function.
                    # Note: the function should call self.update_io_parameter(io1, io2),otherwise nothing happens
                    # io2 = valueType((valueFunction(self, value)))
            else:
                printError(str(key) + ' not in Values_Factor!')
            # self.update_io_parameter(io1, io2) # TODO: note this here was commented out

        if update_parameters: self.update_parameters()

    def update_io_parameter(self, io1, io2):
        self.io1_list.append(io1)
        self.io2_list.append(io2)

    def update_parameters(self):
        if len(self.Update_QuadRF_channels) != 0:
            quadController = QuadRFMOTController(initialValues=self.Exp_Values,
                                                 updateChannels=self.Update_QuadRF_channels, armChannels=True,
                                                 debugging=False,
                                                 continuous=False)  # updates values on QuadRF (uploads table)
            # Note at this point QuadRF is not yet armed, meanning it's not waiting for trigger, until OPX is updated
            self.Update_QuadRF_channels.clear()  # After only updating changed channels, reset this variable.

        # Don't update if there is nothing to update
        if self.io1_list.__len__() == 0:
            return

        # Sets the last value of both of them to be 0, such that the last loop execution will cause the OPX to continue
        # self.update_io_parameter(6, 3)  # Update filler to do 2 more cycles to let the MOT stabilize
        self.io1_list.append(0)
        self.io2_list.append(0)

        # Sending of 1st parameter pair
        self.qm.set_io_values(self.io1_list[0], self.io2_list[0])
        # IOs are sent one after the other, meaning that even if IO1 > 0, doesn't mean that IO2 has been updated.
        # This is solved by adding an initial pause/resume:
        while not self.job.is_paused():
            pass
        self.job.resume()

        # Go over the rest of the parameters:
        for i in range(1, self.io1_list.__len__(), 1):
            while not self.job.is_paused():
                pass
            self.qm.set_io_values(self.io1_list[i], self.io2_list[i])
            self.job.resume()

        # Reset the lists
        self.io1_list = []
        self.io2_list = []

        # Save Config_Table to a file with timestamp:
        self.saveConfigTable()

        # quadController.plotTables()

    def saveConfigTable(self, path='.\\Config_Table Archive\\'):
        timeStamp = time.strftime("%d-%m-%Y %H_%M_%S", time.localtime())
        helmCoilsState = 'Working on it'  # self.getHelmholtzCoilsState()
        expValuesToSave = self.Exp_Values
        if 'opticalPowerCalibrationFunction' in expValuesToSave:
            expValuesToSave[
                'opticalPowerCalibrationFunction'] = 'Experiment values contained power calibration function. It is, however, impossible to save a function as JSON. Thus, this line is added as a warning'
        saveData = {'Time Stamp': timeStamp, 'Exp_Values': self.Exp_Values, 'Helmholtz Coils': helmCoilsState}
        with open(f'{path}{timeStamp}.json', 'w') as file:
            json.dump(saveData, file, indent=4)

    def getHelmholtzCoilsState(self):
        printGreen('Getting Helmholtz coils state...')
        res = {}
        try:
            HCoilsController = HMP4040Visa()
            for ch in [1, 2, 3]:
                HCoilsController.setOutput(ch=ch)
                state = HCoilsController.getOutputState()
                current = HCoilsController.getCurrent()
                voltage = HCoilsController.getVoltage()
                res['Channel %s' % str(ch)] = {'Current': current, 'Voltage': voltage, 'state': state}
            res['Output state'] = HCoilsController.getGeneralOutputState()
            HCoilsController = None  # closing controller
        except Exception:
            res = 'Error getting Helmholtz coils data'
            printError(res)
        return res

    # ## Boolean functions: ##
    def update_pgc_final_amplitude(self, final_amplitude):
        x = np.log(int(1 / (1 - final_amplitude)))
        print('Called update_pgc_final_amplitude', final_amplitude)
        # self.update_lin_pgc_final_amplitude(final_amplitude)
        # self.update_exp_pgc_final_amplitude(final_amplitude)
        return x

    def MOT_switch(self, with_atoms):
        if with_atoms:
            self.update_io_parameter(1, 0)
        else:
            self.update_io_parameter(2, 0)

    def Linear_PGC_switch(self, Bool):
        self.update_io_parameter(2, Bool)

    def Transit_Exp_switch(self, Bool):
        self.Transit_switch = Bool
        self.update_io_parameter(3, Bool)

    def CRUS_Exp_switch(self, Bool):
        self.update_io_parameter(4, Bool)

    def Spectrum_Exp_switch(self, Bool):
        self.update_io_parameter(4, Bool)

    def QRAM_Exp_switch(self, Bool):
        self.update_io_parameter(4, Bool)

    def AntiHelmholtz_Delay_switch(self, Bool):
        self.update_io_parameter(5, Bool)

    def Max_Probe_counts_switch(self, Bool):
        self.update_io_parameter(6, Bool)

    # ## Antihelmholtz delay variable update functions: ##

    def update_AntiHelmholtz_delay(self, new_delay):  # In units of [msec]
        self.update_io_parameter(10, int(new_delay * 1e6 / 4))

    ## update N_snaps: ##

    ## PGC variable update functions: ##
    def update_N_Snaps(self, N_Snaps):  # In units of [msec]
        if self.Exp_Values['Buffer_Cycles'] < 0:
            self.update_io_parameter(46, 3)  # Update 3 buffer cycles
        self.update_io_parameter(45, N_Snaps)  # Update N_snaps

    def update_PGC_prep_time(self, prep_time):  # In units of [msec]
        self.pgc_prep_duration = int(prep_time * 1e6 / 4)
        self.update_io_parameter(21, int(prep_time * 1e6 / 4))

    ## Fountain variable update functions: ##

    def update_fountain_prep_time(self, prep_time):  # In units of [msec]
        self.fountain_prep_duration = int(prep_time * 1e6 / 4)
        self.update_io_parameter(32, int(prep_time * 1e6 / 4))

    def update_fountain_Delta_freq(self, df):
        self.fountain_aom_chirp_rate = int(df * 1e3 / (self.fountain_prep_duration * 4))  # mHz/nsec
        print(self.fountain_aom_chirp_rate)
        self.update_io_parameter(39, self.fountain_aom_chirp_rate)

    ## Measuring variable update functions: ##

    def MeasureOD(self, OD_Time):
        self.update_io_parameter(41, int(OD_Time * 1e6 / 4))

    def MeasureOD_SNSPDs_delay(self, OD_Delay):
        self.update_io_parameter(44, int(OD_Delay * 1e6 / 4))

    def MeasureOD_SNSPDs_duration(self, duration):
        self.update_io_parameter(45, int(duration / self.M_window * 1e6 / 4))

    def Get_Max_Probe_counts(self, repetitions):

        self.Max_Probe_counts_switch(True)
        self.update_parameters()

        Probe_N_handle = self.job.result_handles.get("North_Probe")
        Probe_S_handle = self.job.result_handles.get("South_Probe")

        Probe_N_handle.wait_for_values(1)
        Probe_S_handle.wait_for_values(1)

        Probe_N_res = Probe_N_handle.fetch_all()
        Probe_S_res = Probe_S_handle.fetch_all()

        Probe_counts_North = [np.sum(Probe_N_res.tolist())]
        Probe_counts_South = [np.sum(Probe_S_res.tolist())]
        print(Probe_counts_South)

        for n in range(repetitions - 1):
            while np.sum(Probe_S_res.tolist()) == Probe_counts_South[-1]:
                Probe_N_res = Probe_N_handle.fetch_all()
                Probe_S_res = Probe_S_handle.fetch_all()
            Probe_counts_North.append(np.sum(Probe_N_res.tolist()))
            Probe_counts_South.append(np.sum(Probe_S_res.tolist()))

        print(r'$Max Probe_N = %.1f$' % ((np.average(Probe_counts_North) * 1000) / self.M_time,) + '[MPhotons/sec]')
        print(r'$Max Probe_S = %.1f$' % ((np.average(Probe_counts_South) * 1000) / self.M_time,) + '[MPhotons/sec]')

        self.Max_Probe_counts_switch(False)
        self.update_parameters()

        return [(np.average(Probe_counts_North) * 1000) / self.M_time,
                (np.average(Probe_counts_South) * 1000) / self.M_time]

    def get_avg_num_of_photons_in_det_pulse(self, det_pulse_len, sprint_sequence_delay, num_of_det_pulses, num_of_sprint_sequences):
        self.avg_num_of_photons_in_det_pulse = np.zeros(num_of_det_pulses)
        for i in range(num_of_det_pulses):
            detection_puls_ind =\
            list(np.arange((sprint_sequence_delay + i * det_pulse_len),(sprint_sequence_delay + (i + 1) * det_pulse_len)))
            self.avg_num_of_photons_in_det_pulse[i] = \
                np.sum(self.tt_S_SPRINT_events[detection_puls_ind]) / num_of_sprint_sequences

    def get_handles_from_OPX_Server(self, Num_Of_dets):
        '''
         gets handles of timetags and counts from OPX
        :return:
        '''
        Counts_handle = []
        tt_handle = []

        for i in Num_Of_dets:
            Counts_handle.append(self.job.result_handles.get("Det"+str(i)+"_Counts"))
            tt_handle.append(self.job.result_handles.get("Det"+str(i)+"_Probe_TT"))
        self.MZ_BP_counts_handle = self.job.result_handles.get("Bright_Port_Counts")
        self.MZ_DP_counts_handle = self.job.result_handles.get("Dark_Port_Counts")
        self.Phase_Correction_values_handle = self.job.result_handles.get("Phase_Correction_array")
        # self.Phase_Correction_handle = self.job.result_handles.get("Phase_Correction_for_min")
        FLR_handle = self.job.result_handles.get("FLR_measure")
        return Counts_handle, tt_handle, FLR_handle

    def get_tt_from_handles(self, Num_Of_dets, Counts_handle, tt_handle, FLR_handle):
        self.counts_res = []
        self.tt_res = []

        # handles wait for values
        for i in range(len(Num_Of_dets)):
            tt_handle[i].wait_for_values(1)
            Counts_handle[i].wait_for_values(1)
        FLR_handle.wait_for_values(1)
        self.MZ_BP_counts_handle.wait_for_values(1)
        self.MZ_DP_counts_handle.wait_for_values(1)
        self.Phase_Correction_values_handle.wait_for_values(1)
        # self.Phase_Correction_handle.wait_for_values(1)

        # add the tt and counts to python vars
        for i in range(len(Num_Of_dets)):
            self.counts_res.append(Counts_handle[i].fetch_all())
            self.tt_res.append(tt_handle[i].fetch_all())

        self.MZ_BP_counts_res = self.MZ_BP_counts_handle.fetch_all()
        self.MZ_DP_counts_res = self.MZ_DP_counts_handle.fetch_all()
        self.Phase_Correction_values = self.Phase_Correction_values_handle.fetch_all()
        # self.Phase_Correction = self.Phase_Correction_handle.fetch_all()
        self.FLR_res = -FLR_handle.fetch_all()

        # clear tt's vecs from last data
        self.tt_measure = []
        self.tt_BP_measure = []
        self.tt_DP_measure = []
        self.tt_N_measure = []
        self.tt_S_measure = []
        self.tt_FS_measure = []

        # remove zero-padding from tt-res and append into tt_measure
        for i in range(len(Num_Of_dets)):  # for different detectors
            # self.tt_measure.append([self.tt_res[i][(index * Config.vec_size): (index * Config.vec_size + counts)].tolist() for index, counts in
            #                         enumerate(self.counts_res[i])])
            self.tt_measure.append(self.tt_res[i][1:(self.tt_res[i][0])].tolist())
            self.tt_measure[i] = [elm + Config.detector_delays[i] for elm in self.tt_measure[i]
                                  if ((elm % self.M_window != 0) & (elm != 9999984) &
                                      ((elm + Config.detector_delays[i]) <= self.M_window))]  # Due to an unresolved bug in the OPD there are "ghost" readings of timetags equal to the maximum time of measuring window.
            self.tt_measure[i].sort()
        # unify detectors and windows within detectors and create vector of tt's for each direction (bright port, dark port, north, south and from FS) and sort them
        self.tt_BP_measure = sorted(sum(self.tt_measure[:2], [])) # unify detectors 1-3 and windows within detectors
        self.tt_DP_measure = sorted(sum(self.tt_measure[2:4], [])) # unify detectors 1-3 and windows within detectors
        self.tt_N_measure = sorted(sum(self.tt_measure[7:], [])) # unify detectors 1-3 and windows within detectors
        self.tt_S_measure = sorted(sum(self.tt_measure[4:5], [])) # unify detectors 6-8 and windows within detectors
        self.tt_FS_measure = sorted(sum(self.tt_measure[5:7], [])) # unify detectors 6-8 and windows within detectors
        # plt.plot(self.tt_S_measure)
        self.Phase_Correction_vec = self.Phase_Correction_values['value_0']
        self.Phase_Correction_min_vec = self.Phase_Correction_values['value_1']
        self.Phase_Correction_value = self.Phase_Correction_values['value_1'][-1]


    def get_pulses_bins(self,det_pulse_len,sprint_pulse_len,num_of_det_pulses,num_of_sprint_pulses,
                        sprint_sequence_delay,num_of_sprint_sequences,num_init_zeros,num_fin_zeros,num_between_zeros):
        '''
        generating bin vwctor with edges at the efges of pulses (such that the histogram would give the number of
        clicks within pulse).
        :param det_pulse_len:
        :param sprint_pulse_len:
        :param num_of_det_pulses:
        :param num_of_sprint_pulses:
        :param sprint_sequence_delay:
        :return:
        '''
        pulses_bins = [np.maximum(sprint_sequence_delay,0)]
        for i in range(num_of_sprint_sequences):
            pulses_bins += (pulses_bins[-1]+np.arange(det_pulse_len,(num_of_det_pulses+1)*det_pulse_len,det_pulse_len)).tolist()
            pulses_bins += (pulses_bins[-1]+np.arange(sprint_pulse_len,(num_of_sprint_pulses+1)*sprint_pulse_len,sprint_pulse_len)).tolist()
            pulses_bins[-1] += num_fin_zeros+num_init_zeros-num_between_zeros
        return pulses_bins


    def divide_to_reflection_trans(self,det_pulse_len,sprint_pulse_len,sprint_sequence_delay_S,sprint_sequence_delay_N,
                                   num_of_det_pulses,
                                   num_of_sprint_pulses,num_of_sprint_sequences):
        '''
        dividing south and north tt's vectros into reflection and transmission vectors, by:
        1. converting them into counts vector in same time spacing as the initial pulse.
        2. taking the relevant pulses from the sequence to each vector (reflection/transmission)
        :return:
        '''
        # create histogram with self.M_window bins
        self.pulses_bins_S = self.get_pulses_bins(det_pulse_len,sprint_pulse_len,num_of_det_pulses,num_of_sprint_pulses
                             ,sprint_sequence_delay_S,num_of_sprint_sequences,Config.num_init_zeros_S
                                                  ,Config.num_fin_zeros_S,Config.num_between_zeros)
        self.pulses_bins_N = self.get_pulses_bins(det_pulse_len,sprint_pulse_len,num_of_det_pulses,num_of_sprint_pulses
                             ,sprint_sequence_delay_N,num_of_sprint_sequences,Config.num_init_zeros_N
                                                  ,Config.num_fin_zeros_N,Config.num_between_zeros)
        self.tt_histogram_N,_ = np.histogram(self.tt_N_measure, self.pulses_bins_N)
        self.tt_histogram_S,_ = np.histogram(self.tt_S_measure, self.pulses_bins_S)

        ### build reflection and transmission pulses indices, by setting non zero at indices elements ###
        # build transmission South indices

        # transmission and reflection would take alternately reflection and transmission pulses
        # transmission
        tt_histogram_transmission = np.zeros(np.size(self.tt_histogram_S))
        tt_histogram_reflection = np.zeros(np.size(self.tt_histogram_S))
        for i in range(len(self.tt_histogram_S)):
            if i % (num_of_det_pulses + num_of_sprint_pulses) % 2 == 0:
                tt_histogram_transmission[i] = self.tt_histogram_S[i]
                tt_histogram_reflection[i] = self.tt_histogram_N[i]
            else:
                tt_histogram_transmission[i] = self.tt_histogram_N[i]
                tt_histogram_reflection[i] = self.tt_histogram_S[i]

        return tt_histogram_transmission, tt_histogram_reflection

    def divide_tt_to_reflection_trans(self, sprint_pulse_len, num_of_det_pulses):
        '''
        A function designed to count the number of photons reflected or transmitted for each sequence, such that,
        for the detection pulses the number of photons will be accumulated for each sequence and for the SPRINT
        pulses there will be number of reflected or transmitted photons for each SPRINT pulse.
        :param sprint_pulse_len: the length in [ns] of the SPRINT pulses in the sequence.
        :param num_of_det_pulses: the number of detection pulses in the sequence.
        :return:
        '''
        self.num_of_det_reflections_per_seq_S = np.zeros(self.number_of_QRAM_sequences)
        self.num_of_det_reflections_per_seq_N = np.zeros(self.number_of_QRAM_sequences)
        self.num_of_det_transmissions_per_seq_S = np.zeros(self.number_of_QRAM_sequences)
        self.num_of_det_transmissions_per_seq_N = np.zeros(self.number_of_QRAM_sequences)

        self.num_of_SPRINT_reflections_per_seq_S = np.zeros([self.number_of_QRAM_sequences, self.number_of_SPRINT_pulses_per_seq])
        self.num_of_SPRINT_reflections_per_seq_N = np.zeros([self.number_of_QRAM_sequences, self.number_of_SPRINT_pulses_per_seq])
        self.num_of_SPRINT_transmissions_per_seq_S = np.zeros([self.number_of_QRAM_sequences, self.number_of_SPRINT_pulses_per_seq])
        self.num_of_SPRINT_transmissions_per_seq_N = np.zeros([self.number_of_QRAM_sequences, self.number_of_SPRINT_pulses_per_seq])

        # tt_small_perturb = []
        for element in self.tt_N_measure + self.tt_BP_measure + self.tt_DP_measure:
            tt_inseq = element % self.QRAM_sequence_len
            if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                self.num_of_det_reflections_per_seq_S[(element-1)//self.QRAM_sequence_len] += self.filter_S[tt_inseq]
                self.num_of_det_transmissions_per_seq_N[(element-1)//self.QRAM_sequence_len] += self.filter_N[tt_inseq]
            # else:  # The part of the SPRINT pulses in the sequence
            #     SPRINT_pulse_num = (tt_inseq - self.end_of_det_pulse_in_seq) // (sprint_pulse_len + Config.num_between_zeros)
            #     if SPRINT_pulse_num < self.number_of_SPRINT_pulses_per_seq:
            #         self.num_of_SPRINT_reflections_per_seq_S[(element-1) // self.QRAM_sequence_len][SPRINT_pulse_num] += 1
            #         self.num_of_SPRINT_transmissions_per_seq_N[(element-1) // self.QRAM_sequence_len][SPRINT_pulse_num] += 1
            #         tt_small_perturb+=[element]

        for element in self.tt_S_measure + self.tt_FS_measure:
            tt_inseq = element % self.QRAM_sequence_len
            if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                self.num_of_det_reflections_per_seq_N[(element-1)//self.QRAM_sequence_len] += self.filter_N[tt_inseq]
                self.num_of_det_transmissions_per_seq_S[(element-1)//self.QRAM_sequence_len] += self.filter_S[tt_inseq]
            # else:  # The part of the SPRINT pulses in the sequence
            #     SPRINT_pulse_num = (tt_inseq - self.end_of_det_pulse_in_seq) // (sprint_pulse_len + Config.num_between_zeros)
            #     if SPRINT_pulse_num < self.number_of_SPRINT_pulses_per_seq:
            #         self.num_of_SPRINT_reflections_per_seq_N[(element-1) // self.QRAM_sequence_len][SPRINT_pulse_num] += 1
            #         self.num_of_SPRINT_transmissions_per_seq_S[(element-1) // self.QRAM_sequence_len][SPRINT_pulse_num] += 1

    def divide_BP_and_DP_counts(self, num_of_seq_per_count=50):
        '''
        A function designed to count the number of photons reflected or transmitted for each sequence, such that,
        for the detection pulses the number of photons will be accumulated for each sequence and for the SPRINT
        pulses there will be number of reflected or transmitted photons for each SPRINT pulse.
        :param sprint_pulse_len: the length in [ns] of the SPRINT pulses in the sequence.
        :param num_of_det_pulses: the number of detection pulses in the sequence.
        :return:
        '''

        # self.num_of_BP_counts_per_n_sequences = np.zeros(self.number_of_QRAM_sequences)
        # self.num_of_DP_counts_per_n_sequences = np.zeros(self.number_of_QRAM_sequences)
        self.num_of_BP_counts_per_n_sequences = np.zeros(self.number_of_QRAM_sequences//num_of_seq_per_count)
        self.num_of_DP_counts_per_n_sequences = np.zeros(self.number_of_QRAM_sequences//num_of_seq_per_count)

        for element in self.tt_BP_measure:
            tt_inseq = element % self.QRAM_sequence_len
            if tt_inseq > self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                self.num_of_BP_counts_per_n_sequences[(element-1)//(self.QRAM_sequence_len * num_of_seq_per_count)] += \
                    np.ceil(Config.QRAM_Exp_Square_samples_Late[tt_inseq])

        for element in self.tt_DP_measure:
            tt_inseq = element % self.QRAM_sequence_len
            if tt_inseq > self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                self.num_of_DP_counts_per_n_sequences[(element-1)//(self.QRAM_sequence_len * num_of_seq_per_count)] += \
                    np.ceil(Config.QRAM_Exp_Square_samples_Late[tt_inseq])

    def plot_folded_tt_histogram(self):

        plt.figure()
        plt.plot(sum(np.reshape(self.tt_histogram_N,[int(len(self.tt_histogram_N)/9),9])),label='tt_hist_N')
        plt.plot(sum(np.reshape(self.tt_histogram_S,[int(len(self.tt_histogram_S)/9),9])),label='tt_hist_S')
        plt.legend()

    def find_transits_and_sprint_events(self, detection_condition, num_of_det_pulses, num_of_sprint_pulses):
        '''
        find transits of atoms by searching events with detection_condition[0] photons before
        the sprint sequence and detection_condition[1] after.
        :param detection_condition:
        :param num_of_det_pulses:
        :param num_of_sprint_pulses:
        :return:
        '''
        transit_sequences_init_tt = []
        sprints_events = []
        num_of_pulses = num_of_det_pulses + num_of_sprint_pulses

        num_of_detected_atom = 0
        num_of_sprints = 0
        detection_pulse_range = np.zeros([self.number_of_QRAM_sequences, num_of_det_pulses], dtype=int)
        sprint_pulse_range = np.zeros([self.number_of_QRAM_sequences, num_of_sprint_pulses], dtype=int)

        #define ranges for tt_histogram (a vector that counts number of photons in each pulse)
        for i in range(self.number_of_QRAM_sequences):
            detection_pulse_range[i] =\
                list(range(i*num_of_pulses, i*num_of_pulses+num_of_det_pulses))
            # sprint pulse range starts from last detection pulse andfo num_of_sprint_pulses
            sprint_pulse_range[i] = \
                list(range(int(detection_pulse_range[i][-1])+1, int(detection_pulse_range[i][-1])+1+num_of_sprint_pulses))
        # find transits and sprint events
        for j in range(self.number_of_QRAM_sequences - 2):
            if \
                sum(self.tt_histogram_reflection[detection_pulse_range[j]]) >= detection_condition[0] and \
                        sum(self.tt_histogram_reflection[detection_pulse_range[j + 1]]) >= detection_condition[1]:
                    num_of_detected_atom += 1
                    transit_sequences_init_tt.append(j * self.QRAM_sequence_len)
                    if self.tt_histogram_reflection[sprint_pulse_range[j][0]]:
                        num_of_sprints += 1
                        sprints_events.append(self.tt_histogram_reflection[sprint_pulse_range[j]])
        sprints_data = [sprints_events, num_of_sprints]
        atom_detect_data = [transit_sequences_init_tt, num_of_detected_atom]
        return atom_detect_data, sprints_data

    def find_transits_and_sprint_events_changed(self, cond=[2,2,2], minimum_number_of_seq_detected=2):
        '''
        Find transits of atoms by searching events that satisfy the number of reflected photons per sequence required at
        each cond place with minimum number of conditions needed to be satisfied, defined by the minimum_number_of_seq_detected.
        For example:
        given cond=[2,1,2] and minimum_number_of_seq_detected=2, if in 3 consecutive sequences we get either 2-1-0 or
        0-2-1 or 2-0-2, the condition is satisfied and it is defined as a transit.
        :param cond: The condition that need to be met for number of reflections per detection pulses in sequence.
        :param minimum_number_of_seq_detected: The number of terms needed to be satisfied per cond vector.
        :return:
        '''

        current_transit = []
        self.all_transits_seq_indx = []  # Array of the sequence indexes of all recognized transits per cycle. The length of it will be the number of all transits at the current cycle.
        self.reflection_SPRINT_data_per_transit = []  # Array of vectors with data on the number of reflections per SPRINT pulse in sequence.
        self.transmission_SPRINT_data_per_transit = []  # Array of vectors with data on the number of transmissions per SPRINT pulse in sequence.

        for i in range(len(self.num_of_det_reflections_per_seq) - len(cond) + 1):
            cond_check = (self.num_of_det_reflections_per_seq[i:(i + len(cond))] >= cond).astype(int)
            if sum(cond_check) >= minimum_number_of_seq_detected:
                current_transit = np.unique(
                    current_transit + [*range(i + np.where(cond_check != 0)[0][0], (i + len(cond)))]).tolist()
            elif len(current_transit) > 1:
                current_transit = current_transit[:np.where(self.num_of_det_reflections_per_seq[current_transit] >= min(cond))[0][-1] + 1]
                if self.all_transits_seq_indx:
                    if bool(set(current_transit) & set(self.all_transits_seq_indx[-1])):
                        current_transit = self.all_transits_seq_indx[-1] + current_transit[1:]
                        self.all_transits_seq_indx = self.all_transits_seq_indx[:-1]
                        # self.reflection_SPRINT_data_per_transit = self.reflection_SPRINT_data_per_transit[:-1]
                        # self.transmission_SPRINT_data_per_transit = self.transmission_SPRINT_data_per_transit[:-1]
                self.all_transits_seq_indx.append(current_transit)
                # self.reflection_SPRINT_data_per_transit.append([self.num_of_SPRINT_reflections_per_seq[elem].tolist()
                #                                                 for elem in current_transit[:-1]])
                # self.transmission_SPRINT_data_per_transit.append([self.num_of_SPRINT_transmissions_per_seq[elem].tolist()
                #                                                   for elem in current_transit[:-1]])
                current_transit = []
        if len(current_transit) > 1:
            current_transit = current_transit[
                              :np.where(self.num_of_det_reflections_per_seq[current_transit] >= min(cond))[0][-1] + 1]
            if self.all_transits_seq_indx:
                if bool(set(current_transit) & set(self.all_transits_seq_indx[-1])):
                    current_transit = self.all_transits_seq_indx[-1] + current_transit[1:]
                    self.all_transits_seq_indx = self.all_transits_seq_indx[:-1]
                    # self.reflection_SPRINT_data_per_transit = self.reflection_SPRINT_data_per_transit[:-1]
                    # self.transmission_SPRINT_data_per_transit = self.transmission_SPRINT_data_per_transit[:-1]
            self.all_transits_seq_indx.append(current_transit)
            # self.reflection_SPRINT_data_per_transit.append([self.num_of_SPRINT_reflections_per_seq[elem].tolist()
            #                                                 for elem in current_transit[:-1]])
            # self.transmission_SPRINT_data_per_transit.append([self.num_of_SPRINT_transmissions_per_seq[elem].tolist()
            #                                                   for elem in current_transit[:-1]])

    def get_pulses_location_in_seq(self, delay, seq=Config.QRAM_Exp_Gaussian_samples_S, smearing = int(Config.num_between_zeros/2)):
        '''
        A function that uses the original sequence samples that the OPX uses, in order to obtain the location of the
        pulses in the sequence and build a filter. The user may add smearing which is the value that is added before and
        after each pulse in the sequence to match the filter to the performance of the physical system (AOMs).
        :param delay: Between the actual sequence pulse location to the location of the folded data
        :param seq: The sequence of pulses from which the filter is generated.
        :param smearing: The value that is added to the filter before and after the each pulse in the sequence.
        :return:
        '''
        seq_filter = (np.array(seq) > 0).astype(int)
        seq_filter = np.roll(seq_filter, delay)
        seq_indx = np.where(seq_filter > 0)[0]
        seq_filter_with_smearing = seq_filter
        pulses_loc = []
        if seq_indx.any():
            start_indx = seq_indx[0]
            for i in range(1, len(seq_indx)):
                if seq_indx[i] - seq_indx[i-1] > 1:
                    pulses_loc.append((start_indx - int(smearing), seq_indx[i-1] + int(smearing)))
                    for j in range(start_indx - int(smearing), seq_indx[i-1] + int(smearing)):
                        seq_filter_with_smearing[j] = 1
                    start_indx = seq_indx[i]
            pulses_loc.append((start_indx - int(smearing), seq_indx[-1] + int(smearing)))
            for j in range(start_indx - int(smearing), seq_indx[-1] + int(smearing)):
                seq_filter_with_smearing[j] = 1
        return pulses_loc,  seq_filter_with_smearing

    def get_avg_num_of_photons_in_seq_pulses(self, seq, pulse_loc, tt_measure):
        avg_num_of_photons_in_seq_pulses = []
        try:
            real_number_of_seq = math.ceil(max(tt_measure)/len(Config.QRAM_Exp_Gaussian_samples_S))
            # print('Real number of seq = %d' %real_number_of_seq)
        except:
            real_number_of_seq = self.number_of_QRAM_sequences
            # print('Max number of seq')
        for t in pulse_loc:
            avg_num_of_photons_in_seq_pulses.append((sum(seq[t[0]:t[1]]) + seq[t[1]]) / (real_number_of_seq * 0.167))  # Sagnac configuiration efficiency 16.7%
        return avg_num_of_photons_in_seq_pulses

    def get_max_value_in_seq_pulses(self, seq, pulse_loc):
        max_value_in_seq_pulses = []
        for t in pulse_loc:
            if t[1] < (len(seq) + 1):
                max_value_in_seq_pulses.append(max(seq[t[0]:t[1]]))
            else:
                max_value_in_seq_pulses.append(max(seq[t[0]:]))
        return max_value_in_seq_pulses

    def num_of_photons_txt_box_loc(self, pulse_loc):
        if pulse_loc:
            box_loc = (np.sum(pulse_loc, axis=1)/2).astype(int)
        else:
            box_loc = np.array([])
        return box_loc

    def fold_tt_histogram(self, exp_sequence_len):

        self.folded_tt_S = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_N = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_BP = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_DP = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_FS = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_S_directional = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_N_directional = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_BP_timebins = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_DP_timebins = np.zeros(exp_sequence_len, dtype=int)

        for x in [elem for elem in self.tt_S_measure]:
            self.folded_tt_S[x % exp_sequence_len] += 1
        for x in [elem for elem in self.tt_N_measure]:
            self.folded_tt_N[x % exp_sequence_len] += 1
        for x in [elem for elem in self.tt_BP_measure]:
            self.folded_tt_BP[x % exp_sequence_len] += 1
        for x in [elem for elem in self.tt_DP_measure]:
            self.folded_tt_DP[x % exp_sequence_len] += 1
        for x in [elem for elem in self.tt_FS_measure]:
            self.folded_tt_FS[x % exp_sequence_len] += 1

        self.folded_tt_S_directional = (np.array(self.folded_tt_S) + np.array(self.folded_tt_FS))
        # self.folded_tt_N_directional = self.folded_tt_N
        self.folded_tt_N_directional[:self.end_of_det_pulse_in_seq] = \
            np.array(self.folded_tt_N_directional[:self.end_of_det_pulse_in_seq]) \
            + np.array(self.folded_tt_BP[:self.end_of_det_pulse_in_seq])\
            + np.array(self.folded_tt_DP[:self.end_of_det_pulse_in_seq])

        self.folded_tt_BP_timebins[self.end_of_det_pulse_in_seq:] = self.folded_tt_BP[self.end_of_det_pulse_in_seq:]
        self.folded_tt_DP_timebins[self.end_of_det_pulse_in_seq:] = self.folded_tt_DP[self.end_of_det_pulse_in_seq:]
        if self.pulses_location_in_seq_A or (Config.sprint_pulse_amp_N[0] > 0):
            self.folded_tt_N_directional[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]] = \
                (np.array(self.folded_tt_N_directional[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]])
                 + (Config.sprint_pulse_amp_Early[1]
                 * np.array(self.folded_tt_BP[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]]))
                 + (Config.sprint_pulse_amp_Early[1]
                 * np.array(self.folded_tt_DP[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]]))
                 )
            self.folded_tt_BP_timebins[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]] = \
                (np.array(self.folded_tt_BP_timebins[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]])
                 - (Config.sprint_pulse_amp_Early[1]
                 * np.array(self.folded_tt_BP[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]]))
                 )
            self.folded_tt_DP_timebins[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]] = \
                (np.array(self.folded_tt_DP_timebins[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]])
                 - (Config.sprint_pulse_amp_Early[1]
                 * np.array(self.folded_tt_DP[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]]))
                 )


    def plt_pulses_roll(self, delay_det_N, delay_det_S, delay_rf_N, delay_rf_S):
        plt.figure()
        plt.plot(np.roll(self.folded_tt_S_delayed, delay_det_S), label='folded_S')
        plt.plot(np.roll(self.folded_tt_N_delayed, delay_det_N), label='folded_N')
        plt.plot(np.roll(self.S_pulses_loc_delayed, delay_rf_S), label='sent_S')
        plt.plot(np.roll(self.N_pulses_loc_delayed, delay_rf_N), label='sent_N')
        plt.legend()
        plt.title([delay_det_N, delay_det_S, delay_rf_N, delay_rf_S])

    def save_tt_to_batch(self, Num_Of_dets, N):
        for i in range(len(Num_Of_dets)):
            self.tt_measure_batch[i] = self.tt_measure_batch[i][-(N - 1):] + [self.tt_measure[i]]
        self.tt_S_measure_batch = self.tt_S_measure_batch[-(N - 1):] + [self.tt_S_measure]
        self.tt_N_measure_batch = self.tt_N_measure_batch[-(N - 1):] + [self.tt_N_measure]
        self.tt_BP_measure_batch = self.tt_BP_measure_batch[-(N - 1):] + [self.tt_BP_measure]
        self.tt_DP_measure_batch = self.tt_DP_measure_batch[-(N - 1):] + [self.tt_DP_measure]
        self.tt_FS_measure_batch = self.tt_FS_measure_batch[-(N - 1):] + [self.tt_FS_measure]
        self.MZ_BP_counts_balancing_batch = self.MZ_BP_counts_balancing_batch[-(N - 1):] + [self.MZ_BP_counts_res['value_0']]
        self.MZ_BP_counts_balancing_check_batch = self.MZ_BP_counts_balancing_check_batch[-(N - 1):] + [self.MZ_BP_counts_res['value_1']]
        self.MZ_DP_counts_balancing_batch = self.MZ_DP_counts_balancing_batch[-(N - 1):] + [self.MZ_DP_counts_res['value_0']]
        self.MZ_DP_counts_balancing_check_batch = self.MZ_DP_counts_balancing_check_batch[-(N - 1):] + [self.MZ_DP_counts_res['value_1']]
        self.Phase_Correction_vec_batch = self.Phase_Correction_vec_batch[-(N - 1):] + [self.Phase_Correction_vec]
        self.Phase_Correction_min_vec_batch = self.Phase_Correction_min_vec_batch[-(N - 1):] + [self.Phase_Correction_min_vec]


    def plot_sprint_figures(self, fig, ax, Num_Of_dets):

        ax[0].clear()
        ax[1].clear()
        ax[2].clear()
        ax[3].clear()
        ax[4].clear()
        ax[5].clear()
        ax[6].clear()
        # ax[7].clear()

        #
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        if self.acquisition_flag:
            props_thresholds = dict(boxstyle='round', edgecolor='green', linewidth=2, facecolor='green', alpha=0.5)
        else:
            playsound('C:/Windows/Media/cricket-1.wav')
            props_thresholds = dict(boxstyle='round', edgecolor='red', linewidth=2, facecolor='red', alpha=0.5)

        avg_BP = np.average(self.num_of_BP_counts_per_n_sequences)
        # avg_BP_before = np.average(self.MZ_BP_counts_res['value_1'][:self.rep_MZ_check])
        # avg_BP_after = np.average(self.MZ_BP_counts_res['value_1'][self.rep_MZ_check:])
        avg_DP = np.average(self.num_of_DP_counts_per_n_sequences)
        # avg_DP_before = np.average(self.MZ_DP_counts_res['value_1'][:self.rep_MZ_check])
        # avg_DP_after = np.average(self.MZ_DP_counts_res['value_1'][self.rep_MZ_check:])

        textstr_thresholds = '# %d - ' % self.Counter + 'Reflections: %d, ' % self.sum_for_threshold + \
                             'Efficiency: %.2f, ' % self.lockingEfficiency + \
                             'Flr: %.2f, ' % (1000 * np.average(self.FLR_res.tolist())) + \
                             'Lock Error: %.3f' % self.lock_err
        textstr_total_reflections = 'Total reflections per cycle "N" = %d \n' % (sum(self.num_of_det_reflections_per_seq_N),)\
                                    + 'Total reflections per cycle "S" = %d \n' % (sum(self.num_of_det_reflections_per_seq_S),) \
                                    +'Average reflections per cycle = %.2f' % (sum(self.num_of_det_reflections_per_seq_accumulated/self.Counter),)
        # textstr_avg_reflections = r'Average reflections per cycle = %.2f' % (sum(self.num_of_det_reflections_per_seq_accumulated/self.Counter),)
        # textstr_total_SPRINT_reflections = 'Total reflections from SPRINT pulses during transits = %d' % (self.total_number_of_reflections_from_SPRINT_pulses_in_transits,)
        textstr_BP_DP = 'Average Bright counts = %.2f \n' % (avg_BP,) + \
                        'Average Dark counts = %.2f \n' % (avg_DP,) + \
                        'Infidelity = %.2f' % (avg_DP/(avg_DP + avg_BP),)
        textstr_BP_DP_BA = 'Infidelity before = %.2f \n' % (self.Infidelity_before,) + \
                           'Infidelity after= %.2f' % (self.Infidelity_after,)

        # for n in range(len(Num_Of_dets)):
        #     ax[0].plot(self.single_det_folded[n], label='detectors %d' % (n+1))
        # # ax[0].plot(self.folded_tt_N, label='"N" detectors')
        # # ax[0].plot(self.folded_tt_S, label='"S" detectors')
        # # ax[0].plot((self.filter_S) * max(self.folded_tt_N + self.folded_tt_S), '--', color='orange', label='Filter "S"')
        # ax[0].set_title('binned timetags from all detectors folded (Live)', fontweight="bold")
        # ax[0].legend(loc='upper right')
        # # ax[0].text(0.4, 1.4, textstr_thresholds, transform=ax[0].transAxes, fontsize=28,
        # #            verticalalignment='top', bbox=props_thresholds)

        # for m in range(len(Num_Of_dets)):
        #     ax[1].plot(self.single_det_folded_accumulated[m], label='detectors %d' % (m+1))
        # # ax[1].plot(self.folded_tt_N_batch, label='"N" detectors')
        # # ax[1].plot(self.folded_tt_S_batch, label='"S" detectors')
        # # ax[1].plot((self.filter_N) * max(self.folded_tt_N_batch + self.folded_tt_S_batch), '--b', label='Filter "N"')
        # ax[1].set_title('binned timetags from all detectors folded (Averaged)', fontweight="bold")
        # ax[1].legend(loc='upper right')

        # ax[0].plot(self.folded_tt_N, label='"N" detectors')
        # ax[0].plot(self.folded_tt_S, label='"S" detectors')
        # ax[0].plot(self.folded_tt_BP, label='"BP" detectors')
        # ax[0].plot(self.folded_tt_DP, label='"DP" detectors')
        # ax[0].plot(self.folded_tt_FS, label='"FS" detectors')
        ax[0].plot(self.folded_tt_N_directional, label='"N" detectors')
        ax[0].plot(self.folded_tt_S_directional, label='"S" detectors')
        ax[0].plot(self.folded_tt_BP_timebins, label='"BP" detectors')
        ax[0].plot(self.folded_tt_DP_timebins, label='"DP" detectors')
        ax[0].plot((self.filter_S) * max(self.folded_tt_N_directional + self.folded_tt_S_directional),
                   '--', color='orange', label='Filter "S"')
        for i in range(len(self.Num_of_photons_txt_box_y_loc_live)):
            ax[0].text(self.Num_of_photons_txt_box_x_loc.tolist()[i], self.Num_of_photons_txt_box_y_loc_live[i],
                       '%.2f' % self.avg_num_of_photons_per_pulse_live[i],
                       horizontalalignment='center', fontsize=12, fontweight='bold', family=['Comic Sans MS'])
        for j in range(len(self.Num_of_photons_txt_box_y_loc_live_MZ)):
            ax[0].text(self.Num_of_photons_txt_box_x_loc_for_MZ_ports.tolist()[j], self.Num_of_photons_txt_box_y_loc_live_MZ[j][0],
                       '%.2f' % self.avg_num_of_photons_per_pulse_live_MZ[j][0],
                       horizontalalignment='center', fontsize=12, fontweight='bold', family=['Comic Sans MS'], color='#2ca02c')
            ax[0].text(self.Num_of_photons_txt_box_x_loc_for_MZ_ports.tolist()[j], self.Num_of_photons_txt_box_y_loc_live_MZ[j][1],
                       '%.2f' % self.avg_num_of_photons_per_pulse_live_MZ[j][1],
                       horizontalalignment='center', fontsize=12, fontweight='bold', family=['Comic Sans MS'], color='#d62728')
        # ax[0].set_ylim(0, 400)
        ax[0].set_title('binned timetags from all detectors folded (Live)', fontweight="bold")
        ax[0].legend(loc='upper right')
        ax[0].text(0.05, 1.4, textstr_thresholds, transform=ax[0].transAxes, fontsize=28,
                   verticalalignment='top', bbox=props_thresholds)

        # ax[1].plot(self.folded_tt_BP_batch, label='"BP" detectors')
        # ax[1].plot(self.folded_tt_DP_batch, label='"DP" detectors')
        ax[1].plot(self.folded_tt_N_directional_batch, label='"N" detectors')
        ax[1].plot(self.folded_tt_S_directional_batch, label='"S" detectors')
        ax[1].plot(self.folded_tt_BP_timebins_batch, label='"BP" detectors')
        ax[1].plot(self.folded_tt_DP_timebins_batch, label='"DP" detectors')
        ax[1].plot((self.filter_N) * max(self.folded_tt_N_directional_batch + self.folded_tt_S_directional_batch),
                   '--b', label='Filter "N"')
        for i in range(len(self.Num_of_photons_txt_box_y_loc)):
            ax[1].text(self.Num_of_photons_txt_box_x_loc.tolist()[i], self.Num_of_photons_txt_box_y_loc[i],
                     '%.2f' % self.avg_num_of_photons_per_pulse[i],
                     horizontalalignment='center', fontsize=12, fontweight='bold', family=['Comic Sans MS'])
        for j in range(len(self.Num_of_photons_txt_box_y_loc_MZ)):
            ax[1].text(self.Num_of_photons_txt_box_x_loc_for_MZ_ports.tolist()[j], self.Num_of_photons_txt_box_y_loc_MZ[j][0],
                       '%.2f' % self.avg_num_of_photons_per_pulse_MZ[j][0],
                       horizontalalignment='center', fontsize=12, fontweight='bold', family=['Comic Sans MS'], color='#2ca02c')
            ax[1].text(self.Num_of_photons_txt_box_x_loc_for_MZ_ports.tolist()[j], self.Num_of_photons_txt_box_y_loc_MZ[j][1],
                       '%.2f' % self.avg_num_of_photons_per_pulse_MZ[j][1],
                       horizontalalignment='center', fontsize=12, fontweight='bold', family=['Comic Sans MS'], color='#d62728')
        ax[1].set_title('binned timetags from all detectors folded (Averaged)', fontweight="bold")
        ax[1].legend(loc='upper right')

        ax[2].plot(self.MZ_BP_counts_res['value_0'], label='MZ Bright port')
        ax[2].plot(self.MZ_DP_counts_res['value_0'], label='MZ Dark port')
        ax[2].plot(self.MZ_BP_counts_res['value_0'] - self.MZ_DP_counts_res['value_0'], label='Dif ports')
        ax[6].tick_params(axis="y", labelcolor='#8c564b')
        ax[6].plot(self.Phase_Correction_vec, label='Phase correction values', color='#8c564b')
        ax[6].plot(self.Phase_Correction_min_vec, label='Phase correction values', color='#9467bd')
        ax[2].set_ylim(0, 1.1 * np.max([self.MZ_BP_counts_res['value_0'], self.MZ_DP_counts_res['value_0']]))
        ax[2].set_title('MZ outputs while locking', fontweight="bold")
        ax[2].legend(loc='upper right')
        ax[6].legend(loc='upper left')
        print(self.Phase_Correction_value)

        max_reflect_avg = max(self.num_of_det_reflections_per_seq_accumulated/self.Counter)
        max_reflect = max(self.num_of_det_reflections_per_seq)
        ax[3].plot(self.num_of_det_reflections_per_seq_accumulated/self.Counter, label='Num of reflections per sequence')
        ax[3].plot(self.num_of_det_reflections_per_seq*0.5*max_reflect_avg/max_reflect*0.3, label='Num of reflections per sequence (Live)')
        ax[3].set_title('Num of reflections per sequence', fontweight="bold")
        ax[3].legend(loc='upper right')
        # ax[3].text(0.1, 0.9*max(self.num_of_det_transmissions_per_seq_accumulated/self.Counter), textstr_avg_reflections, fontsize=14,
        ax[3].text(0.1, 0.9*max_reflect_avg, textstr_total_reflections, fontsize=14,
                   verticalalignment='top', bbox=props)

        ax[4].plot(self.MZ_BP_counts_res['value_1'], label='MZ BP counts before and after')
        ax[4].plot(self.MZ_DP_counts_res['value_1'], label='MZ DP port counts before and after')
        ax[4].axvline(len(self.MZ_DP_counts_res['value_1'])/2, linestyle='--', color='red')
        ax[4].set_title('MZ outputs around experiment', fontweight="bold")
        ax[4].legend(loc='upper right')
        ax[4].text(0.05, 0.6, textstr_BP_DP_BA, transform=ax[4].transAxes, fontsize=14,
                   verticalalignment='top', bbox=props)

        # ax[5].plot(self.num_of_BP_counts_per_n_sequences, label='MZ BP counts per %d seq' % 50)
        # ax[5].plot(self.num_of_DP_counts_per_n_sequences, label='MZ DP counts per %d seq' % 50)
        # ax[5].set_title('MZ outputs during experiment', fontweight="bold")
        # ax[5].legend(loc='upper right')
        # ax[5].text(0.05, 0.6, textstr_BP_DP, transform=ax[5].transAxes, fontsize=14,
        #            verticalalignment='top', bbox=props)

        if self.number_of_transits_live:
            textstr_transit_counts = r'$N_{Transits} = %s $' % (self.number_of_transits_live,) + r'$[Counts]$'
        else:
            textstr_transit_counts = r'$N_{Transits} = %s $' % (0,) + r'$[Counts]$'
        textstr_transit_event_counter = r'$N_{Transits Total} = %s $' % (self.number_of_transits_total,) + '[Counts]\n'\
                                        + textstr_transit_counts
                                        # + textstr_total_SPRINT_reflections

        # ax[4].plot(range(self.number_of_QRAM_sequences), self.seq_transit_events_live, label='Current transit events')
        # ax[4].set(xlabel='Sequence [#]', ylabel='Counts [Photons]')
        # ax[4].set_title('Transits per sequence (live)', fontweight="bold")
        # ax[4].text(0.05, 0.95, textstr_transit_counts, transform=ax[4].transAxes, fontsize=14,
        #            verticalalignment='top', bbox=props)
        ax[5].plot(range(self.number_of_QRAM_sequences), self.seq_transit_events_batched, label='Transit Events Accumulated')
        ax[5].plot(range(self.number_of_QRAM_sequences), self.seq_transit_events_live, label='Transit Events Live')
        ax[5].set(xlabel='Sequence [#]', ylabel='Counts [Photons]')
        ax[5].set_title('Transits per sequence', fontweight="bold")
        ax[5].legend(loc='upper right')
        ax[5].text(0.05, 0.95, textstr_transit_event_counter, transform=ax[5].transAxes, fontsize=14,
                   verticalalignment='top', bbox=props)

        plt.show(block=False)
        plt.pause(1.0)
        ###########################################################################################################

    def init_params_for_save_sprint(self, qram_sequence_len, Num_Of_dets):
        # define empty variables
        self.QRAM_sequence_len=qram_sequence_len
        self.number_of_QRAM_sequences = math.ceil(self.M_window / self.QRAM_sequence_len)
        self.number_of_SPRINT_pulses_per_seq = len(Config.sprint_pulse_amp_S)

        self.tt_measure = []
        self.tt_S_measure = []
        self.tt_measure_batch = [[]] * len(Num_Of_dets)
        self.tt_S_measure_batch = []
        self.tt_N_measure_batch = []
        self.tt_BP_measure_batch = []
        self.tt_DP_measure_batch = []
        self.tt_FS_measure_batch = []
        self.transit_sequences_batch = []
        self.all_transits_seq_indx_batch = []
        self.reflection_SPRINT_data_per_transit_batch = []
        self.transmission_SPRINT_data_per_transit_batch = []
        self.MZ_BP_counts_balancing_batch = []
        self.MZ_BP_counts_balancing_check_batch = []
        self.MZ_DP_counts_balancing_batch = []
        self.MZ_DP_counts_balancing_check_batch = []
        self.Phase_Correction_vec_batch = []
        self.Phase_Correction_min_vec_batch = []

        self.folded_transmission = np.zeros(len(Config.QRAM_Exp_Gaussian_samples_S))
        self.folded_reflection = np.zeros(len(Config.QRAM_Exp_Gaussian_samples_S))

        self.tt_S_binning = np.zeros(self.number_of_QRAM_sequences + 1)
        self.seq_transit_events_live = np.zeros(self.number_of_QRAM_sequences)
        self.seq_transit_events_batched = np.zeros(self.number_of_QRAM_sequences)
        self.tt_S_SPRINT_events = np.zeros(self.QRAM_sequence_len)
        self.tt_S_SPRINT_events_batch = np.zeros(self.QRAM_sequence_len)
        self.single_det_folded = np.zeros((len(Num_Of_dets), self.QRAM_sequence_len))
        self.single_det_folded_accumulated = np.zeros((len(Num_Of_dets), self.QRAM_sequence_len))
        self.num_of_det_reflections_per_seq_accumulated = np.zeros(self.number_of_QRAM_sequences)
        self.num_of_det_transmissions_per_seq_accumulated = np.zeros(self.number_of_QRAM_sequences)

    def Save_SNSPDs_QRAM_Measurement_with_tt(self, N, qram_sequence_len, preComment, lock_err_threshold, transit_condition,
                                             max_probe_counts, filter_delay, reflection_threshold, reflection_threshold_time,
                                             photons_per_det_pulse_threshold, FLR_threshold, exp_flag, with_atoms,
                                             MZ_inidelity_threshold):
        """
        Function for analyzing,saving and displaying data from sprint experiment.
        :param N: Number of maximum experiments (free throws) saved and displayed.
                 program we are looking for transits of atoms next to the toroid and record them.
        :param qram_sequence_len: the number of sprint sequences (detection and sprint pulses combination)
        :param preComment: The comment added at the start of the experiment, usually consisting of unique experiment
                           parameters.
        :param transit_condition:
        :param lock_err_threshold: The maximal error in locking resonator to rb line (in ms, depends on the scan width/vpp)
        :param max_probe_counts: The maximum counts probe counts at each direction measured when cavity isn't locked.
        :param filter_delay: Delay of the filter window compared to original location (can be taken by comparing the
                             time of the 1st detection pulse pick location to the original sequence)
        :param reflection_threshold:
        :param reflection_threshold_time:
        :return:
        """
        if not preComment:
            preComment = pymsgbox.prompt('Add comment to measurement: ', default='', timeout=int(30e3))

        # set constant parameters for the function
        Num_Of_dets = [1, 2, 3, 4, 5, 6, 7, 8]
        delay_in_detection_N = 30 # choose the correct delay in samples to the first detection pulse # TODO: 40?
        delay_in_detection_S = 20 # choose the correct delay in samples to the first detection pulse # TODO: 40?
        det_pulse_len = Config.det_pulse_len+Config.num_between_zeros
        sprint_pulse_len = Config.sprint_pulse_len+Config.num_between_zeros


        # initialize parameters - set zeros vectors and empty lists
        self.init_params_for_save_sprint(qram_sequence_len, Num_Of_dets)

        # get pulses location south and north
        self.pulses_location_in_seq, self.filter_gen = self.get_pulses_location_in_seq(0,
                                                                                       Config.QRAM_Exp_Gaussian_samples_General,
                                                                                       smearing=5)  # smearing=int(Config.num_between_zeros/2))
        self.pulses_location_in_seq_S, self.filter_S = self.get_pulses_location_in_seq(filter_delay[0],
                                                                                       Config.QRAM_Exp_Gaussian_samples_S,
                                                                                       smearing=5)  # smearing=int(Config.num_between_zeros/2))
        self.pulses_location_in_seq_N, self.filter_N = self.get_pulses_location_in_seq(filter_delay[1],
                                                                                       Config.QRAM_Exp_Gaussian_samples_N,
                                                                                       smearing=5)  # smearing=int(Config.num_between_zeros/2))
        self.pulses_location_in_seq_A, self.filter_A = self.get_pulses_location_in_seq(filter_delay[2],
                                                                                       Config.QRAM_Exp_Gaussian_samples_Ancilla,
                                                                                       smearing=5)  # smearing=int(Config.num_between_zeros/2))
        self.Num_of_photons_txt_box_x_loc = np.concatenate((self.num_of_photons_txt_box_loc(self.pulses_location_in_seq_S),
                                                           self.num_of_photons_txt_box_loc(self.pulses_location_in_seq_N),
                                                           self.num_of_photons_txt_box_loc(self.pulses_location_in_seq_A)))
        self.Num_of_photons_txt_box_x_loc_for_MZ_ports = self.num_of_photons_txt_box_loc(self.pulses_location_in_seq[-2:])
        self.num_of_detection_pulses = len(Config.det_pulse_amp_N)
        self.end_of_det_pulse_in_seq = self.pulses_location_in_seq[self.num_of_detection_pulses - 1][1]

        ## Listen for keyboard
        listener = keyboard.Listener(on_press=self.on_key_press)
        listener.start()  # start to listen on a separate thread
        self.keyPress = None
        print('\033[94m' + 'Press ESC to stop measurement.' + '\033[0m')  # print blue

        ####     get tt and counts from OPX to python   #####
        Counts_handle, tt_handle, FLR_handle = self.get_handles_from_OPX_Server(Num_Of_dets)

        ## take data only if
        start = True
        # take threshold from npz ( error from resonator lock PID)
        if exp_flag:
            try:
                self.lock_err = np.abs(np.load(
                    'U:\Lab_2023\Experiment_results\QRAM\Locking_PID_Error\locking_err.npy', allow_pickle=True))  # the error of locking the resontor to Rb line
            except:
                print('error in loading file')
        else:
            self.lock_err = lock_err_threshold / 2

        self.sum_for_threshold = reflection_threshold
        cycle = 0
        # Place holders for results # TODO: ask dor - is it works as we expect?
        while ((self.lock_err > lock_err_threshold) or (self.sum_for_threshold > reflection_threshold) and exp_flag) or start:
            if self.keyPress == 'ESC':
                print('\033[94m' + 'ESC pressed. Stopping measurement.' + '\033[0m')  # print blue
                self.updateValue("Sprint_Exp_switch", False)
                self.update_parameters()
                self.Stop_run_daily_experiment=True
                # Other actions can be added here
                break
            if start:
                if cycle > 2:
                    start = False
                else:
                    cycle += 1
            else:
                print('Above Threshold')
            self.get_tt_from_handles(Num_Of_dets, Counts_handle, tt_handle, FLR_handle)
            self.divide_tt_to_reflection_trans(sprint_pulse_len, self.num_of_detection_pulses)
            self.divide_BP_and_DP_counts(50)
            self.num_of_det_reflections_per_seq = self.num_of_det_reflections_per_seq_S \
                                                  + self.num_of_det_reflections_per_seq_N
            self.num_of_det_transmissions_per_seq = self.num_of_det_transmissions_per_seq_S \
                                                  + self.num_of_det_transmissions_per_seq_N
            # self.num_of_SPRINT_reflections_per_seq = self.num_of_SPRINT_reflections_per_seq_S \
            #                                          + self.num_of_SPRINT_reflections_per_seq_N
            # self.num_of_SPRINT_transmissions_per_seq = self.num_of_SPRINT_transmissions_per_seq_S \
            #                                            + self.num_of_SPRINT_transmissions_per_seq_N
            self.sum_for_threshold = sum(self.num_of_det_reflections_per_seq[-int(reflection_threshold_time//len(Config.QRAM_Exp_Gaussian_samples_S)):])  # summing over the reflection from detection pulses of each sequence corresponding the the reflection_threshold_time
            print(self.lock_err, self.lock_err > lock_err_threshold, self.sum_for_threshold)
            if exp_flag:
                try:
                    self.lock_err = np.abs(np.load(
                        'U:\Lab_2023\Experiment_results\QRAM\Locking_PID_Error\locking_err.npy'))  # the error of locking the resontor to Rb line
                except:
                    pass
                    print('error in loading file')
        ####    end get tt and counts from OPX to python   #####

        self.num_of_det_reflections_per_seq_accumulated += self.num_of_det_reflections_per_seq_S \
                                                           + self.num_of_det_reflections_per_seq_N
        self.num_of_det_transmissions_per_seq_accumulated += self.num_of_det_transmissions_per_seq_S \
                                                             + self.num_of_det_transmissions_per_seq_N

        # divide south and north into reflection and transmission
        self.tt_histogram_transmission, self.tt_histogram_reflection = \
        self.divide_to_reflection_trans(det_pulse_len=det_pulse_len, sprint_pulse_len=sprint_pulse_len,
                                        sprint_sequence_delay_S=delay_in_detection_S, sprint_sequence_delay_N=delay_in_detection_N,
                                        num_of_det_pulses=len(Config.det_pulse_amp_S),
                                        num_of_sprint_pulses=len(Config.sprint_pulse_amp_S),
                                        num_of_sprint_sequences=self.number_of_QRAM_sequences)

        # TODO: needed?
        self.tt_S_SPRINT_events = np.zeros(qram_sequence_len)
        self.tt_S_SPRINT_events_batch = np.zeros(qram_sequence_len)
        self.tt_Single_det_SPRINT_events = np.zeros((len(Num_Of_dets), qram_sequence_len))
        self.tt_Single_det_SPRINT_events_batch = np.zeros((len(Num_Of_dets), qram_sequence_len))
        self.tt_N_det_SPRINT_events_batch = np.zeros(qram_sequence_len)
        self.tt_S_det_SPRINT_events_batch = np.zeros(qram_sequence_len)

        # fold reflections and transmission
        self.fold_tt_histogram(exp_sequence_len=self.QRAM_sequence_len)

        # fold data from different detectors: # TODO: is needed? @Dor
        for i in range(len(Num_Of_dets)):
            # for x in [elem for elem in self.tt_S_measure if elem < self.M_window]: - for debugging assaf
            for x in [elem for elem in self.tt_measure[i]]:
                self.single_det_folded[i][x % self.QRAM_sequence_len] += 1
                self.single_det_folded_accumulated[i][x % self.QRAM_sequence_len] += 1

        # Batch folded tt "N", "S", BP, DP, FS
        self.folded_tt_S_batch = self.folded_tt_S
        self.folded_tt_N_batch = self.folded_tt_N
        self.folded_tt_BP_batch = self.folded_tt_BP
        self.folded_tt_DP_batch = self.folded_tt_DP
        self.folded_tt_FS_batch = self.folded_tt_FS
        self.folded_tt_S_directional_batch = self.folded_tt_S_directional
        self.folded_tt_N_directional_batch = self.folded_tt_N_directional
        self.folded_tt_BP_timebins_batch = self.folded_tt_BP_timebins
        self.folded_tt_DP_timebins_batch = self.folded_tt_DP_timebins

        # get the average number of photons in detection pulse
        self.avg_num_of_photons_per_pulse_S_live = self.get_avg_num_of_photons_in_seq_pulses(
            self.folded_tt_S_directional, self.pulses_location_in_seq_S, self.tt_FS_measure)
        self.avg_num_of_photons_per_pulse_N_live = self.get_avg_num_of_photons_in_seq_pulses(
            self.folded_tt_N_directional, self.pulses_location_in_seq_N, self.tt_BP_measure + self.tt_DP_measure)
        self.avg_num_of_photons_per_pulse_A_live = self.get_avg_num_of_photons_in_seq_pulses(
            (np.array(self.folded_tt_S_directional) + np.array(self.folded_tt_BP_timebins)
            + np.array(self.folded_tt_DP_timebins)).tolist(), self.pulses_location_in_seq_A,
            self.tt_BP_measure + self.tt_DP_measure)
        self.avg_num_of_photons_per_pulse_BP_live = self.get_avg_num_of_photons_in_seq_pulses(
            self.folded_tt_BP_timebins, self.pulses_location_in_seq[-2:], self.tt_BP_measure)
        self.avg_num_of_photons_per_pulse_DP_live = self.get_avg_num_of_photons_in_seq_pulses(
            self.folded_tt_DP_timebins, self.pulses_location_in_seq[-2:], self.tt_DP_measure)
        self.avg_num_of_photons_per_pulse_live = self.avg_num_of_photons_per_pulse_S_live + \
                                                 self.avg_num_of_photons_per_pulse_N_live + \
                                                 self.avg_num_of_photons_per_pulse_A_live
        self.avg_num_of_photons_per_pulse_live_MZ = [[x]+[y] for x, y in zip(self.avg_num_of_photons_per_pulse_BP_live,
                                                                             self.avg_num_of_photons_per_pulse_DP_live)]
            #                                         self.avg_num_of_photons_per_pulse_BP_live + \
            #                                         self.avg_num_of_photons_per_pulse_DP_live
        self.avg_num_of_photons_per_pulse = self.avg_num_of_photons_per_pulse_live
        self.avg_num_of_photons_per_pulse_MZ = self.avg_num_of_photons_per_pulse_live_MZ

        # get box location on the y-axis:
        self.max_value_per_pulse_S_live = self.get_max_value_in_seq_pulses(self.folded_tt_S_directional, self.pulses_location_in_seq_S)
        self.max_value_per_pulse_N_live = self.get_max_value_in_seq_pulses(self.folded_tt_N_directional, self.pulses_location_in_seq_N)
        self.max_value_per_pulse_A_live = self.get_max_value_in_seq_pulses(
            (np.array(self.folded_tt_S_directional) + np.array(self.folded_tt_BP_timebins) +
             np.array(self.folded_tt_DP_timebins)).tolist(), self.pulses_location_in_seq_A)
        self.max_value_per_pulse_BP_live = self.get_max_value_in_seq_pulses(
            self.folded_tt_BP_timebins, self.pulses_location_in_seq[-2:])
        self.max_value_per_pulse_DP_live = self.get_max_value_in_seq_pulses(
            self.folded_tt_DP_timebins, self.pulses_location_in_seq[-2:])
        self.Num_of_photons_txt_box_y_loc_live = self.max_value_per_pulse_S_live + \
                                                 self.max_value_per_pulse_N_live + \
                                                 self.max_value_per_pulse_A_live

        self.Num_of_photons_txt_box_y_loc_live_MZ = [[x]+[y] for x, y in zip(self.max_value_per_pulse_BP_live,
                                                                             self.max_value_per_pulse_DP_live)]
        # self.Num_of_photons_txt_box_y_loc_live_MZ = self.max_value_per_pulse_BP_live + self.max_value_per_pulse_DP_live
        self.Num_of_photons_txt_box_y_loc = self.Num_of_photons_txt_box_y_loc_live
        self.Num_of_photons_txt_box_y_loc_MZ = self.Num_of_photons_txt_box_y_loc_live_MZ

        avg_BP_before = np.average(self.MZ_BP_counts_res['value_1'][:self.rep_MZ_check])
        avg_BP_after = np.average(self.MZ_BP_counts_res['value_1'][self.rep_MZ_check:])
        avg_DP_before = np.average(self.MZ_DP_counts_res['value_1'][:self.rep_MZ_check])
        avg_DP_after = np.average(self.MZ_DP_counts_res['value_1'][self.rep_MZ_check:])
        self.Infidelity_before = avg_DP_before / (avg_DP_before + avg_BP_before)
        self.Infidelity_after = avg_DP_after / (avg_DP_after + avg_BP_after)

        ## record time
        timest = time.strftime("%Y%m%d-%H%M%S")
        datest = time.strftime("%Y%m%d")
        FLR_measurement = []
        Exp_timestr_batch = []

        FLR_measurement = FLR_measurement[-(N - 1):] + [self.FLR_res.tolist()]
        Exp_timestr_batch = Exp_timestr_batch[-(N - 1):] + [timest]

        ### Dor version ###
        self.find_transits_and_sprint_events_changed(cond=transit_condition, minimum_number_of_seq_detected=2)
        self.seq_transit_events_live[[elem for vec in self.all_transits_seq_indx for elem in vec]] += 1
        self.seq_transit_events_batched[[elem for vec in self.all_transits_seq_indx for elem in vec]] += 1
        self.all_transits_seq_indx_batch = self.all_transits_seq_indx_batch[-(N-1):] + [self.all_transits_seq_indx]
        # self.reflection_SPRINT_data_per_transit_batch = self.reflection_SPRINT_data_per_transit_batch[-(N-1):]\
        #                                                 + [self.reflection_SPRINT_data_per_transit]
        # self.transmission_SPRINT_data_per_transit_batch = self.transmission_SPRINT_data_per_transit_batch[-(N-1):]\
        #                                                   + [self.transmission_SPRINT_data_per_transit]
        self.number_of_transits_live = len(self.all_transits_seq_indx)
        self.number_of_transits_total = len([vec for lst in self.all_transits_seq_indx_batch for vec in lst])
        # self.total_number_of_reflections_from_SPRINT_pulses_in_transits = \
        #     np.sum(np.sum([vec for lst in self.reflection_SPRINT_data_per_transit_batch for vec in lst]))
        ###################

        self.save_tt_to_batch(Num_Of_dets, N)
        self.Counter = 1  # Total number of successful cycles
        self.repitions = 1  # Total number of cycles
        self.acquisition_flag = True
        #
        # create figures template
        fig = plt.figure()
        ax1 = plt.subplot2grid((3, 4), (0, 0), colspan=2, rowspan=1)
        ax2 = plt.subplot2grid((3, 4), (0, 2), colspan=2, rowspan=1)
        ax3 = plt.subplot2grid((3, 4), (1, 0), colspan=2, rowspan=1)
        ax4 = plt.subplot2grid((3, 4), (1, 2), colspan=2, rowspan=1)
        ax5 = plt.subplot2grid((3, 4), (2, 0), colspan=2, rowspan=1)
        # ax6 = plt.subplot2grid((3, 4), (2, 1), colspan=0, rowspan=1)
        ax6 = plt.subplot2grid((3, 4), (2, 2), colspan=2, rowspan=1)
        ax7 = ax3.twinx()
        #

        while True:
            if self.keyPress == 'ESC':
                print('\033[94m' + 'ESC pressed. Stopping measurement.' + '\033[0m')  # print blue
                self.updateValue("QRAM_Exp_switch", False)
                self.Stop_run_daily_experiment = True
                self.update_parameters()
                # Other actions can be added here
                break
            self.lockingEfficiency = self.Counter / self.repitions
            print(timest, self.Counter, 'Eff: %.2f' % self.lockingEfficiency, 'Flr: %.2f' % (1000 * np.average(self.FLR_res.tolist())))
            self.repitions += 1
            ######################################## PLOT!!! ###########################################################
            ax = [ax1, ax2, ax3, ax4, ax5, ax6, ax7] # , ax8]
            self.plot_sprint_figures(fig, ax, Num_Of_dets)
            ############################################################################################################
            count = 1

            while True:
                # record time:
                # timest = time.strftime("%Y%m%d,%H%M%S") # TODO: is it needed? already writen above..
                timest = time.strftime("%H%M%S")
                datest = time.strftime("%Y%m%d")
                self.get_tt_from_handles(Num_Of_dets, Counts_handle, tt_handle, FLR_handle)

                # Check if new tt's arrived:
                lenS = min(len(self.tt_S_measure), len(self.tt_S_measure_batch[-1]))
                lenN = min(len(self.tt_N_measure), len(self.tt_N_measure_batch[-1]))
                lenDP = min(len(self.tt_DP_measure), len(self.tt_DP_measure_batch[-1]))
                lenBP = min(len(self.tt_BP_measure), len(self.tt_BP_measure_batch[-1]))
                lenFS = min(len(self.tt_FS_measure), len(self.tt_FS_measure_batch[-1]))
                # Check if the number of same values in the new and last vector are less then 1/2 of the total number of values.
                is_new_tts_S = sum(np.array(self.tt_S_measure[:lenS]) == np.array(self.tt_S_measure_batch[-1][:lenS])) < lenS/2
                is_new_tts_N = sum(np.array(self.tt_N_measure[:lenN]) == np.array(self.tt_N_measure_batch[-1][:lenN])) < lenN/2
                is_new_tts_BP = sum(np.array(self.tt_BP_measure[:lenBP]) == np.array(self.tt_BP_measure_batch[-1][:lenBP])) < lenBP/2
                is_new_tts_DP = sum(np.array(self.tt_DP_measure[:lenDP]) == np.array(self.tt_DP_measure_batch[-1][:lenDP])) < lenDP/2
                is_new_tts_FS = sum(np.array(self.tt_FS_measure[:lenFS]) == np.array(self.tt_FS_measure_batch[-1][:lenFS])) < lenFS/2

                # # Check if the number of same values in the new and last vector are less then 1/2 of the total number of values.
                # is_new_tts_S = self.tt_S_measure[0] != self.tt_S_measure_batch[-1][0]
                # is_new_tts_N = self.tt_N_measure[0] != self.tt_N_measure_batch[-1][0]
                # is_new_tts_BP = self.tt_BP_measure[0] != self.tt_BP_measure_batch[-1][0]
                # is_new_tts_DP = self.tt_DP_measure[0] != self.tt_DP_measure_batch[-1][0]
                # is_new_tts_FS = self.tt_FS_measure[0] != self.tt_FS_measure_batch[-1][0]

                print(count)
                count += 1
                if is_new_tts_N or is_new_tts_S or is_new_tts_BP or is_new_tts_DP or is_new_tts_FS:
                    break
                if self.keyPress == 'ESC':
                    print('\033[94m' + 'ESC pressed. Stopping measurement.' + '\033[0m')  # print blue
                    self.updateValue("QRAM_Exp_switch", False)
                    self.update_parameters()
                    # Other actions can be added here
                    break

            # assaf - if x=self.M_window the index is out of range so i added 1
            try:
                self.lock_err = np.abs(np.load(
                    'U:\Lab_2023\Experiment_results\QRAM\Locking_PID_Error\locking_err.npy'))  # the error of locking the resontor to Rb line
            except:
                pass
                print('error in loading file')

            self.divide_tt_to_reflection_trans(sprint_pulse_len, self.num_of_detection_pulses)
            self.divide_BP_and_DP_counts(50)

            self.num_of_det_reflections_per_seq = self.num_of_det_reflections_per_seq_S \
                                                  + self.num_of_det_reflections_per_seq_N
            self.num_of_det_transmissions_per_seq = self.num_of_det_transmissions_per_seq_S \
                                                  + self.num_of_det_transmissions_per_seq_N
            # self.num_of_SPRINT_reflections_per_seq = self.num_of_SPRINT_reflections_per_seq_N \
            #                                          + self.num_of_SPRINT_reflections_per_seq_S
            # self.num_of_SPRINT_transmissions_per_seq = self.num_of_SPRINT_transmissions_per_seq_N \
            #                                            + self.num_of_SPRINT_transmissions_per_seq_S
            self.sum_for_threshold = sum(self.num_of_det_reflections_per_seq[-int(reflection_threshold_time//len(Config.QRAM_Exp_Gaussian_samples_S)):])  # summing over the reflection from detection pulses of each sequence corresponding the the reflection_threshold_time

            # fold reflections and transmission
            self.single_det_folded = np.zeros((len(Num_Of_dets), self.QRAM_sequence_len))
            self.fold_tt_histogram(exp_sequence_len=self.QRAM_sequence_len)

            # get the average number of photons in detection pulse
            self.avg_num_of_photons_per_pulse_S_live = self.get_avg_num_of_photons_in_seq_pulses(
                self.folded_tt_S_directional, self.pulses_location_in_seq_S, self.tt_FS_measure)
            self.avg_num_of_photons_per_pulse_N_live = self.get_avg_num_of_photons_in_seq_pulses(
                self.folded_tt_N_directional, self.pulses_location_in_seq_N, self.tt_BP_measure + self.tt_DP_measure)
            self.avg_num_of_photons_per_pulse_A_live = self.get_avg_num_of_photons_in_seq_pulses(
                (np.array(self.folded_tt_S_directional) + np.array(self.folded_tt_BP_timebins)
                 + np.array(self.folded_tt_DP_timebins)).tolist(), self.pulses_location_in_seq_A,
                self.tt_BP_measure + self.tt_DP_measure)
            self.avg_num_of_photons_per_pulse_BP_live = self.get_avg_num_of_photons_in_seq_pulses(
                self.folded_tt_BP_timebins, self.pulses_location_in_seq[-2:], self.tt_BP_measure)
            self.avg_num_of_photons_per_pulse_DP_live = self.get_avg_num_of_photons_in_seq_pulses(
                self.folded_tt_DP_timebins, self.pulses_location_in_seq[-2:], self.tt_DP_measure)
            self.avg_num_of_photons_per_pulse_live = self.avg_num_of_photons_per_pulse_S_live + \
                                                     self.avg_num_of_photons_per_pulse_N_live + \
                                                     self.avg_num_of_photons_per_pulse_A_live
            self.avg_num_of_photons_per_pulse_live_MZ = [[x] + [y] for x, y in
                                                         zip(self.avg_num_of_photons_per_pulse_BP_live,
                                                             self.avg_num_of_photons_per_pulse_DP_live)]

            # self.avg_num_of_photons_per_pulse_live_MZ = self.avg_num_of_photons_per_pulse_BP_live + \
            #                                             self.avg_num_of_photons_per_pulse_DP_live

            # get box location on the y-axis:
            self.max_value_per_pulse_S_live = self.get_max_value_in_seq_pulses(self.folded_tt_S_directional,
                                                                               self.pulses_location_in_seq_S)
            self.max_value_per_pulse_N_live = self.get_max_value_in_seq_pulses(self.folded_tt_N_directional,
                                                                               self.pulses_location_in_seq_N)
            self.max_value_per_pulse_A_live = self.get_max_value_in_seq_pulses(
                (np.array(self.folded_tt_S_directional) + np.array(self.folded_tt_BP_timebins) +
                 np.array(self.folded_tt_DP_timebins)).tolist(), self.pulses_location_in_seq_A)
            self.max_value_per_pulse_BP_live = self.get_max_value_in_seq_pulses(
                self.folded_tt_BP_timebins, self.pulses_location_in_seq[-2:])
            self.max_value_per_pulse_DP_live = self.get_max_value_in_seq_pulses(
                self.folded_tt_DP_timebins, self.pulses_location_in_seq[-2:])
            self.Num_of_photons_txt_box_y_loc_live = self.max_value_per_pulse_S_live + \
                                                     self.max_value_per_pulse_N_live + \
                                                     self.max_value_per_pulse_A_live
            self.Num_of_photons_txt_box_y_loc_live_MZ = [[x] + [y] for x, y in zip(self.max_value_per_pulse_BP_live,
                                                                                   self.max_value_per_pulse_DP_live)]
            # self.Num_of_photons_txt_box_y_loc_live_MZ = self.max_value_per_pulse_BP_live + self.max_value_per_pulse_DP_live

            avg_BP_before = np.average(self.MZ_BP_counts_res['value_1'][:self.rep_MZ_check])
            avg_BP_after = np.average(self.MZ_BP_counts_res['value_1'][self.rep_MZ_check:])
            avg_DP_before = np.average(self.MZ_DP_counts_res['value_1'][:self.rep_MZ_check])
            avg_DP_after = np.average(self.MZ_DP_counts_res['value_1'][self.rep_MZ_check:])
            self.Infidelity_before = avg_DP_before / (avg_DP_before + avg_BP_before)
            self.Infidelity_after = avg_DP_after / (avg_DP_after + avg_BP_after)

            if exp_flag:
                if (self.lock_err > lock_err_threshold) or \
                        ((1000 * np.average(self.FLR_res.tolist()) < FLR_threshold) and with_atoms) or \
                        (np.average(experiment.avg_num_of_photons_per_pulse_live) > photons_per_det_pulse_threshold) or \
                        (self.Infidelity_before > MZ_inidelity_threshold) or \
                        (self.Infidelity_after > MZ_inidelity_threshold):
                    self.acquisition_flag = False
                else:
                    self.acquisition_flag = True

            if (self.sum_for_threshold < reflection_threshold or not exp_flag) and self.acquisition_flag:
                print('Sum of reflections: %d' % self.sum_for_threshold)
                self.num_of_det_reflections_per_seq_accumulated += self.num_of_det_reflections_per_seq_S \
                                                                   + self.num_of_det_reflections_per_seq_N
                self.num_of_det_transmissions_per_seq_accumulated += self.num_of_det_transmissions_per_seq_S \
                                                                   + self.num_of_det_transmissions_per_seq_N

                self.seq_transit_events_live = np.zeros(self.number_of_QRAM_sequences)

                ### Find transits and build histogram:  ###
                self.find_transits_and_sprint_events_changed(cond=transit_condition, minimum_number_of_seq_detected=2)
                self.seq_transit_events_live[[vec for elem in self.all_transits_seq_indx for vec in elem]] += 1
                self.seq_transit_events_batched[[vec for elem in self.all_transits_seq_indx for vec in elem]] += 1
                self.all_transits_seq_indx_batch = self.all_transits_seq_indx_batch[-(N - 1):]\
                                                   + [self.all_transits_seq_indx]
                # self.reflection_SPRINT_data_per_transit_batch = self.reflection_SPRINT_data_per_transit_batch[-(N - 1):] \
                #                                                 + [self.reflection_SPRINT_data_per_transit]
                # self.transmission_SPRINT_data_per_transit_batch = self.transmission_SPRINT_data_per_transit_batch[-(N - 1):] \
                #                                                   + [self.transmission_SPRINT_data_per_transit]
                self.number_of_transits_live = len(self.all_transits_seq_indx)
                self.number_of_transits_total = len([vec for lst in self.all_transits_seq_indx_batch for vec in lst])
                # self.total_number_of_reflections_from_SPRINT_pulses_in_transits = \
                #     np.sum(np.sum([vec for lst in self.reflection_SPRINT_data_per_transit_batch for vec in lst]))
                ###########################################

                # Batch folded tt "N" and "S"
                self.folded_tt_S_batch = (self.folded_tt_S_batch * (self.Counter - 1) + self.folded_tt_S) / self.Counter
                self.folded_tt_N_batch = (self.folded_tt_N_batch * (self.Counter - 1) + self.folded_tt_N) / self.Counter
                self.folded_tt_BP_batch = (self.folded_tt_BP_batch * (self.Counter - 1) + self.folded_tt_BP) / self.Counter
                self.folded_tt_DP_batch = (self.folded_tt_DP_batch * (self.Counter - 1) + self.folded_tt_DP) / self.Counter
                self.folded_tt_FS_batch = (self.folded_tt_FS_batch * (self.Counter - 1) + self.folded_tt_FS) / self.Counter
                self.folded_tt_S_directional_batch = (self.folded_tt_S_directional_batch * (self.Counter - 1)
                                                      + self.folded_tt_S_directional) / self.Counter
                self.folded_tt_N_directional_batch = (self.folded_tt_N_directional_batch * (self.Counter - 1)
                                                      + self.folded_tt_N_directional) / self.Counter
                self.folded_tt_BP_timebins_batch = (self.folded_tt_BP_timebins_batch * (self.Counter - 1)
                                                    + self.folded_tt_BP_timebins) / self.Counter
                self.folded_tt_DP_timebins_batch = (self.folded_tt_DP_timebins_batch * (self.Counter - 1)
                                                    + self.folded_tt_DP_timebins) / self.Counter

                # get the average number of photons in detection pulse
                self.avg_num_of_photons_per_pulse_S = self.get_avg_num_of_photons_in_seq_pulses(
                    self.folded_tt_S_directional_batch, self.pulses_location_in_seq_S,[])
                self.avg_num_of_photons_per_pulse_N = self.get_avg_num_of_photons_in_seq_pulses(
                    self.folded_tt_N_directional_batch, self.pulses_location_in_seq_N,[])
                self.avg_num_of_photons_per_pulse_A = self.get_avg_num_of_photons_in_seq_pulses(
                    (np.array(self.folded_tt_S_directional_batch) + np.array(self.folded_tt_BP_timebins_batch)
                     + np.array(self.folded_tt_DP_timebins_batch)).tolist(), self.pulses_location_in_seq_A,[])
                self.avg_num_of_photons_per_pulse_BP = self.get_avg_num_of_photons_in_seq_pulses(
                    self.folded_tt_BP_timebins_batch, self.pulses_location_in_seq[-2:],[])
                self.avg_num_of_photons_per_pulse_DP = self.get_avg_num_of_photons_in_seq_pulses(
                    self.folded_tt_DP_timebins_batch, self.pulses_location_in_seq[-2:],[])
                self.avg_num_of_photons_per_pulse = self.avg_num_of_photons_per_pulse_S + \
                                                    self.avg_num_of_photons_per_pulse_N + \
                                                    self.avg_num_of_photons_per_pulse_A
                self.avg_num_of_photons_per_pulse_MZ = [[x] + [y] for x, y in
                                                             zip(self.avg_num_of_photons_per_pulse_BP,
                                                                 self.avg_num_of_photons_per_pulse_DP)]
                # self.avg_num_of_photons_per_pulse_MZ = self.avg_num_of_photons_per_pulse_BP + \
                #                                        self.avg_num_of_photons_per_pulse_DP
                # get box location on the y-axis:
                self.max_value_per_pulse_S = self.get_max_value_in_seq_pulses(self.folded_tt_S_directional_batch,
                                                                              self.pulses_location_in_seq_S)
                self.max_value_per_pulse_N = self.get_max_value_in_seq_pulses(self.folded_tt_N_directional_batch,
                                                                              self.pulses_location_in_seq_N)
                self.max_value_per_pulse_A = self.get_max_value_in_seq_pulses(
                    (np.array(self.folded_tt_S_directional_batch) + np.array(self.folded_tt_BP_timebins_batch) +
                     np.array(self.folded_tt_DP_timebins_batch)).tolist(), self.pulses_location_in_seq_A)
                self.max_value_per_pulse_BP = self.get_max_value_in_seq_pulses(
                    self.folded_tt_BP_timebins_batch, self.pulses_location_in_seq[-2:])
                self.max_value_per_pulse_DP = self.get_max_value_in_seq_pulses(
                    self.folded_tt_DP_timebins_batch, self.pulses_location_in_seq[-2:])
                self.Num_of_photons_txt_box_y_loc = self.max_value_per_pulse_S + self.max_value_per_pulse_N + \
                                                    self.max_value_per_pulse_A
                self.Num_of_photons_txt_box_y_loc_MZ = [[x] + [y] for x, y in zip(self.max_value_per_pulse_BP,
                                                                                       self.max_value_per_pulse_DP)]
                # self.Num_of_photons_txt_box_y_loc_MZ = self.max_value_per_pulse_BP + self.max_value_per_pulse_DP

                # fold for different detectors: # TODO: delete after everything works
                for i in range(len(Num_Of_dets)):
                    for x in [elem for elem in self.tt_measure[i]]:
                        self.single_det_folded[i][x % self.QRAM_sequence_len] += 1
                        self.single_det_folded_accumulated[i][x % self.QRAM_sequence_len] += 1

                FLR_measurement = FLR_measurement[-(N - 1):] + [self.FLR_res.tolist()]
                Exp_timestr_batch = Exp_timestr_batch[-(N - 1):] + [timest]
                self.save_tt_to_batch(Num_Of_dets, N)
                if self.Counter == N:
                    print('\033[94m' + f'finished {N} Runs, {"with" if with_atoms else "without"} atoms' + '\033[0m')  # print blue
                    self.updateValue("QRAM_Exp_switch", False)
                    self.update_parameters()
                    # Other actions can be added here
                    break
                if self.Counter < N:
                    self.Counter += 1


        ############################################## END WHILE LOOP #################################################

        # For debuging:
        # self.most_common([x for vec in self.tt_S_measure_batch + self.tt_N_measure_batch for x in vec])

        ## Adding comment to measurement [prompt whether stopped or finished regularly]
        if self.Counter < N:
            aftComment = pymsgbox.prompt('Add comment to measurement: ', default='', timeout=int(30e3))
        else:
            aftComment = ''

        # if aftComment == 'Timeout': aftComment = None

        #### Handle file-names and directories #####
        ## Saving: np.savez(filedir, data = x) #note: @filedir is string of full directory; data is the queyword used to read @x from the file:
        ## Loading: file = np.load(f, allow_pickle = True)
        ##          x = file['data']

        #### ------ Save results ------
        #  -------   Create dir
        root_dirname = f'U:\\Lab_2023\\Experiment_results\\QRAM\\{datest}\\'
        dirname = root_dirname + f'{timest}_Photon_TimeTags\\'  # Specific experiment dir
        dirname_Det = dirname + 'AllDetectors\\'
        dirname_N = dirname + 'North(8)\\'
        dirname_S = dirname + 'South(5)\\'
        dirname_D = dirname + 'Dark(3,4)\\'
        dirname_B = dirname + 'Bright(1,2)\\'
        dirname_FS = dirname + 'FastSwitch(6,7)\\'
        dirname_balancing = dirname + 'BalancingRes\\'
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        if not os.path.exists(dirname_Det):
            os.makedirs(dirname_Det)
        if not os.path.exists(dirname_N):
            os.makedirs(dirname_N)
        if not os.path.exists(dirname_S):
            os.makedirs(dirname_S)
        if not os.path.exists(dirname_D):
            os.makedirs(dirname_D)
        if not os.path.exists(dirname_B):
            os.makedirs(dirname_B)
        if not os.path.exists(dirname_FS):
            os.makedirs(dirname_FS)
        if not os.path.exists(dirname_balancing):
            os.makedirs(dirname_balancing)

        # ----  msmnt files names  -----
        # Counter_str = (self.Counter)
        filename_Det_tt = []
        for i in Num_Of_dets:
            filename_Det_tt.append(f'Det'+str(i)+f'_timetags.npz')
        filename_S_tt = f'South_timetags.npz'
        filename_N_tt = f'North_timetags.npz'
        filename_DP_tt = f'Dark_timetags.npz'
        filename_BP_tt = f'Bright_timetags.npz'
        filename_FS_tt = f'FS_timetags.npz'
        filename_N_folded = f'North_timetags_folded_to_seq.npz'
        filename_S_folded = f'South_timetags_folded_to_seq.npz'
        filename_DP_folded = f'Dark_timetags_folded_to_seq.npz'
        filename_BP_folded = f'Bright_timetags_folded_to_seq.npz'
        filename_FS_folded = f'FS_timetags_folded_to_seq.npz'
        filename_FLR = f'Flouresence.npz'
        filename_timestamp = f'Drops_time_stamps.npz'
        filname_sequence_S = f'South_sequence_vector.npz'
        filname_sequence_N = f'North_sequence_vector.npz'
        filname_sequence_Early = f'Early_sequence_vector.npz'
        filname_sequence_Late = f'Late_sequence_vector.npz'
        filname_sequence_FS = f'FS_sequence_vector.npz'
        filename_reflection_averaged = f'reflection_from_detection_pulses_per_seq_averaged.npz'
        filename_transits_events = f'seq_transit_events_batched.npz'
        filename_transits_batched = f'all_transits_seq_indx_batch.npz'
        filename_MZ_BP_balancing_batched = f'MZ_BrightPort_counts_balancing_batch.npz'
        filename_MZ_BP_balancing_check_batched = f'MZ_BrightPort_counts_balancing_check_batch.npz'
        filename_MZ_DP_balancing_batched = f'MZ_DarkPort_counts_balancing_batch.npz'
        filename_MZ_DP_balancing_check_batched = f'MZ_DarkPort_counts_balancing_check_batch.npz'
        filename_MZ_balancing_Phase_Correction_vec_batch = f'MZ_balancing_Phase_Correction_vec_batch.npz'
        filename_MZ_balancing_Phase_Correction_min_vec_batch = f'MZ_balancing_Phase_Correction_min_vec_batch.npz'
        # filename_SPRINT_reflections_batched = f'all_SPRINT_pulses_reflections_during_transits.npz'
        # filename_SPRINT_transmissions_batched = f'all_SPRINT_pulses_transmissions_during_transits.npz'
        filename_experimentPlot = f'Experiment_plot.png'

        if len(FLR_measurement) > 0:
            np.savez(dirname + filename_FLR, FLR_measurement)
        if len(Exp_timestr_batch) > 0:
            np.savez(dirname + filename_timestamp, Exp_timestr_batch)
            np.savez(dirname + filname_sequence_S, Config.QRAM_Exp_Gaussian_samples_S)
            np.savez(dirname + filname_sequence_N, Config.QRAM_Exp_Gaussian_samples_N)
            np.savez(dirname + filname_sequence_Early, Config.QRAM_Exp_Square_samples_Early)
            np.savez(dirname + filname_sequence_Late, Config.QRAM_Exp_Square_samples_Late)
            np.savez(dirname + filname_sequence_FS, (1-np.array(Config.QRAM_Exp_Square_samples_FS)).tolist())
            plt.savefig(dirname + filename_experimentPlot, bbox_inches='tight')
        for i in range(len(Num_Of_dets)):
            if len(self.tt_measure_batch[i]) > 0:
                np.savez(dirname_Det + filename_Det_tt[i], self.tt_measure_batch[i])
        if len(self.folded_tt_N_batch) > 0:
            np.savez(dirname + filename_N_folded, self.folded_tt_N_batch)
        if len(self.folded_tt_S_batch) > 0:
            np.savez(dirname + filename_S_folded, self.folded_tt_S_batch)
        if len(self.folded_tt_DP_batch) > 0:
            np.savez(dirname + filename_DP_folded, self.folded_tt_DP_batch)
        if len(self.folded_tt_BP_batch) > 0:
            np.savez(dirname + filename_BP_folded, self.folded_tt_BP_batch)
        if len(self.folded_tt_FS_batch) > 0:
            np.savez(dirname + filename_FS_folded, self.folded_tt_FS_batch)
        if len(self.tt_S_measure_batch) > 0:
            np.savez(dirname_S + filename_S_tt, self.tt_S_measure_batch)
        if len(self.tt_N_measure_batch) > 0:
            np.savez(dirname_N + filename_N_tt, self.tt_N_measure_batch)
        if len(self.tt_DP_measure_batch) > 0:
            np.savez(dirname_D + filename_DP_tt, self.tt_DP_measure_batch)
        if len(self.tt_BP_measure_batch) > 0:
            np.savez(dirname_B + filename_BP_tt, self.tt_BP_measure_batch)
        if len(self.tt_FS_measure_batch) > 0:
            np.savez(dirname_FS + filename_FS_tt, self.tt_FS_measure_batch)
        if len(self.num_of_det_reflections_per_seq_accumulated) > 0:
            np.savez(dirname + filename_reflection_averaged, self.num_of_det_reflections_per_seq_accumulated/self.Counter)
        if len(self.seq_transit_events_batched) > 0:
            np.savez(dirname + filename_transits_events, self.seq_transit_events_batched)
        if len(self.all_transits_seq_indx_batch) > 0:
            np.savez(dirname + filename_transits_batched, self.all_transits_seq_indx_batch)
            # np.savez(dirname + filename_SPRINT_reflections_batched, self.reflection_SPRINT_data_per_transit_batch)
            # np.savez(dirname + filename_SPRINT_transmissions_batched, self.transmission_SPRINT_data_per_transit_batch)
        if len(self.MZ_BP_counts_balancing_check_batch) > 0:
            np.savez(dirname_balancing + filename_MZ_BP_balancing_batched, self.MZ_BP_counts_balancing_batch)
            np.savez(dirname_balancing + filename_MZ_BP_balancing_check_batched, self.MZ_BP_counts_balancing_check_batch)
            np.savez(dirname_balancing + filename_MZ_DP_balancing_batched, self.MZ_DP_counts_balancing_batch)
            np.savez(dirname_balancing + filename_MZ_DP_balancing_check_batched, self.MZ_DP_counts_balancing_check_batch)
            np.savez(dirname_balancing + filename_MZ_balancing_Phase_Correction_vec_batch, self.Phase_Correction_vec_batch)
            np.savez(dirname_balancing + filename_MZ_balancing_Phase_Correction_min_vec_batch, self.Phase_Correction_min_vec_batch)

        # if len(all_transits_batch) > 0:
        #     np.savez(dirname_S + filename_S_transits, all_transits_batch)

        ### Edit comments file ####
        cmntDir = os.path.join(root_dirname, 'daily_experiment_comments.csv')
        cmnt_header = 'Date,Time,IgnoreValid,Atoms,Cycles,Comment'
        if not os.path.exists(cmntDir):
            # Write header line
            try:
                with open(cmntDir, "a") as commentsFile:
                    commentsFile.write(cmnt_header + '\n')
            except:
                print('Could not save comments, error writing to comments-file.')

        if preComment is not None: cmnt = preComment + '; '
        if aftComment is not None: cmnt = preComment + aftComment
        if preComment is None and aftComment is None: cmnt = 'No comment.'
        if 'ignore' in cmnt:
            experiment_success = 'ignore'
        else:
            experiment_success = 'valid'
        full_line = f'{datest},{timest},{experiment_success},{with_atoms},{self.Counter},{cmnt}'
        try:
            with open(cmntDir, "a") as commentsFile:
                commentsFile.write(full_line + '\n')
        except:
            print('Could not save comments, error writing to comments-file.')

        experiment_cmnt = 'transit condition: ' + str(transit_condition) + '; ' +\
                          'reflection threshold: ' + str(reflection_threshold) + '@' +\
                          str(int(reflection_threshold_time/1e6)) + 'ms'

        comments = {'comments': experiment_cmnt}
        try:
            with open(f'{dirname}experiment_comments.txt', 'w') as file:
                json.dump(comments, file, indent=4)
        except Exception:
            pass

        self.saveConfigTable(path=dirname)
        try:
            with open(f'{dirname}max_probe_counts.txt', 'w') as file: json.dump(max_probe_counts, file, indent=4)
        except Exception as e:
            print(e)
        for qrdCtrl in self.QuadRFControllers:
            qrdCtrl.saveLinesAsCSV(f'{dirname}QuadRF_table.csv')
        ## ------------------ end of saving section -------

    def Start_QRAM_Exp_with_tt(self, N=100, qram_sequence_len=int(len(Config.QRAM_Exp_Gaussian_samples_S)),
                               transit_condition=[2, 2, 2], preComment=None, lock_err_threshold=0.01, filter_delay=[0,0,0],
                               reflection_threshold=100, reflection_threshold_time=1e6,
                               photons_per_det_pulse_threshold=12, FLR_threshold=0.11, MZ_inidelity_threshold=0.05,
                               Exp_flag=True, with_atoms=True):
        # Max_probe_counts = self.Get_Max_Probe_counts(3)  # return the average maximum probe counts of 3 cycles.
        Max_probe_counts = None  # return the average maximum probe counts of 3 cycles.
        self.QRAM_Exp_switch(True)
        self.MOT_switch(with_atoms)
        self.update_parameters()
        self.Save_SNSPDs_QRAM_Measurement_with_tt(N, qram_sequence_len, preComment, lock_err_threshold,
                                                  transit_condition, Max_probe_counts, filter_delay,
                                                  reflection_threshold, reflection_threshold_time,
                                                  photons_per_det_pulse_threshold, FLR_threshold, exp_flag=Exp_flag,
                                                  with_atoms=with_atoms,MZ_inidelity_threshold=MZ_inidelity_threshold)
    ## MW spectroscopy variable update functions: ##

    def MW_spec_frequency(self, freq):
        self.update_io_parameter(51, int(freq))  # In unit of [Hz]

    def MW_spec_MW_pulse_duration(self, p_length):
        self.update_io_parameter(52, p_length * int(1e3 / 4))  # In units of [us]

    def MW_spec_OD_pulse_duration(self, p_length):
        self.update_io_parameter(53, p_length * int(1e3 / 4))  # In units of [us]

    def MW_spec_Repetition_times(self, reps):
        self.update_io_parameter(54, reps)

    def MW_spec_Delta_freq(self, D_f):
        self.update_io_parameter(55, int(D_f))  # In units of [Hz]

    def MW_spec_n_MOT_cycles(self, N_snaps):
        return self.Film_graber(int(N_snaps))

    def MW_spec_detuning(self, det):
        self.MW_spec_frequency(det + self.MW_start_frequency)

    def MW_spec_scan(self, min_f, max_f, delta_f):
        N_snaps = int((max_f - min_f) / delta_f)  # Number of MOT cycles required
        self.MW_spec_n_MOT_cycles(N_snaps)
        self.MW_spec_Delta_freq(int(delta_f))
        self.MW_spec_detuning(int(min_f))
        return N_snaps

    def on_key_press(self, key):
        if key == keyboard.Key.esc:
            self.keyPress = 'ESC'

    def most_common(self, lst):
        # get an iterable of (item, iterable) pairs
        SL = sorted((x, i) for i, x in enumerate(lst))
        # print 'SL:', SL
        groups = itertools.groupby(SL, key=operator.itemgetter(0))

        # auxiliary function to get "quality" for an item
        def _auxfun(g):
            item, iterable = g
            count = 0
            min_index = len(lst)
            for _, where in iterable:
                count += 1
                min_index = min(min_index, where)
            if count > 20:
                print('item %r, count %r, minind %r' % (item, count, min_index))
            return count, -min_index

        # pick the highest-count/earliest item
        return max(groups, key=_auxfun)[0]

    def run_daily_experiment(self, day_experiment, transit_condition, preComment, lock_err_threshold, filter_delay,
                             reflection_threshold, reflection_threshold_time, FLR_threshold, MZ_inidelity_threshold,
                             Exp_flag=True):
        with_atoms_bool = True
        for i in range(len(day_experiment)):
            if with_atoms_bool:
                Comment = preComment + ' with atoms'
            else:
                Comment = preComment + ' without atoms'
            self.Start_QRAM_Exp_with_tt(N=day_experiment[i], transit_condition=transit_condition, preComment=Comment,
                                        lock_err_threshold=lock_err_threshold, filter_delay=filter_delay,
                                        reflection_threshold=reflection_threshold,
                                        reflection_threshold_time=reflection_threshold_time,
                                        FLR_threshold=FLR_threshold,
                                        MZ_inidelity_threshold=MZ_inidelity_threshold,
                                        Exp_flag=Exp_flag,
                                        with_atoms=with_atoms_bool)  # except KeyboardInterrupt:
            with_atoms_bool = not with_atoms_bool
            if self.Stop_run_daily_experiment:
                break

    # @Values_Factor hold factoring for values - each key has an array:
    # So, as Prep_duration is written in ms, we have to factor it into an int value, but in 4ns, for the OPX.
    # These are only factored in function "updateValue" in OPX control
    Values_Factor = {
        'Transit_Exp_switch': [bool, None, Transit_Exp_switch],
        'Spectrum_Exp_switch': [bool, None, Spectrum_Exp_switch],
        "QRAM_Exp_switch": [bool, None, QRAM_Exp_switch],
        'MOT_duration': [int, 1e6 / Config.MOT_pulse_len],
        'AntiHelmholtz_delay': [int, 1e6 / 4],  # [msec]
        'Post_MOT_delay': [int, 1e6 / 4],  # [msec]
        'N_Snaps': [int, None, update_N_Snaps],
        'Buffer_Cycles': [int, 1],
        ## PGC parameters ##
        "PGC_duration": [int, 1e6 / 4],  # [msec]
        "PGC_prep_duration": [int, None, update_PGC_prep_time],

        ## Fountain parameters ##
        "Pre_PGC_Fountain_duration": [int, 1e6 / 4],
        "Fountain_duration": [int, 1e6 / 4],
        "Fountain_prep_duration": [int, None, update_fountain_prep_time],
        "Fountain_final_Delta_freq": [float, None, update_fountain_Delta_freq],
        "Fountain_aom_chirp_rate": [float, None, update_fountain_Delta_freq],
        # The frequency in the frame of reference where the atoms have just enough velocity to reach the toroids [Hz].
        # Free-Fall parameters:
        # The time take the atoms to reach the Toroids and fall back = 2 * sqrt(2 * (Toroid height[m] - MOT height[m]) / g[m/sec^2]) * 1e3[msec]
        'FreeFall_duration': [int, 1e6 / 4],
        # We use 1.4 because there is a limit for the pulse duration. = 63ms
        'Coil_timing': [int, 1e6 / 4],  # Zeeman coils turn on time from start of free fall[msec]
        # Imaging parameters:
        'Snapshot_Intervals': [int, 1e6 / 4],  # [msec]
        'PrePulse_duration': [int, 1e6 / 4],  # [msec]
        'Trigger_delay': [int, 1e6 / 4],  # [msec]
        'Pulse_1_duration': [int, 1e6 / 4],  # [msec]
        'Pulse_1_decay_duration': [int, 1e6 / 4],  # [msec]
        'InterPulses_duration': [int, 1e6 / 4],
        'Pulse_2_duration': [int, 1e6 / 4],  # [msec]
        # OD parameters:
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


if __name__ == "__main__":
    # try:
    experiment = OPX(Config.config)
    # experiment.Start_QRAM_Exp_with_tt(N=1000, transit_condition=[2, 2, 2], preComment='ignore', lock_err_threshold=0.01,
    #                                   filter_delay=[0, 0, 0],
    #                                   reflection_threshold=10000, reflection_threshold_time=1e6, FLR_threshold=-0.11)   # except KeyboardInterrupt:
    experiment.run_daily_experiment([10000, 20] * 2, transit_condition=[2, 1, 2], preComment='"|1c, (0 + 1)t> experiment"', lock_err_threshold=0.005,
                                      filter_delay=[0, 9, 0],
                                      reflection_threshold=400, reflection_threshold_time=10e6, FLR_threshold=0.06, Exp_flag=False)

    #     experiment.job.halt()
    #     experiment.qmm.reset_data_processing()
    # finally:
    #     experiment.job.halt()
    #     experiment.qmm.reset_data_processing()
    # experiment.updateValue('PrePulse_duration', 5)
    # experiment.update_parameters()
    # experiment.Repeat_Measurement(10)
    # experiment.Repeat_Measurement(10)
