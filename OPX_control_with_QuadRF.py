import Config_with_SNSPDs_and_QuadRF as Config
from Config_Table import Initial_Values, Phases_Names  # , Values_Factor
from quadRFMOTController import QuadRFMOTController

from qm.QuantumMachinesManager import QuantumMachinesManager
from qm.qua import *
from qm.qua import lib
from scipy import signal
from qm import SimulationConfig
import matplotlib.pyplot as plt
import numpy as np
import time, json
import sys
import os
from pynput import keyboard
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


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


all_elements = ["Cooling_Sequence", "MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+", "AntiHelmholtz_Coils",
                "FLR_detection", "Measurement", "AOM_2-2/3'", "AOM_2-2/3'_detuned", "AOM_2-3'_for_interference"]  #, "PULSER_N", "AOM_S", "AOM_LO")


def MOT(mot_repetitions):
    """
    The MOT function is used to play the MOT. To that end, we send RF signal to AOM_0, AOM_+, AOM_- for the duration of the MOT.

    Parameters
    ----------
    :param mot_repetitions: derived from the MOT duration, and is calculated in the Experiment parameters section
    """
    FLR = declare(fixed)
    # update_frequency("AOM_2-2/3'", 180e6)
    align("Cooling_Sequence", "MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+", "AntiHelmholtz_Coils", "Zeeman_Coils",
          "AOM_2-2/3'", "AOM_2-3'_for_interference", "FLR_detection", "Measurement")  #, "PULSER_N", "AOM_S", "AOM_LO")

    ## MOT build-up ##
    n = declare(int)
    m = declare(int)
    play("Detection" * amp(FLR * 0), "FLR_detection", duration=4)  # we dont know why this works, but Yoav from QM made us write this line to solve an alignment problem we had in the next 2 for loops
    with for_(n, 1, n <= mot_repetitions, n + 1):
        play("MOT" * amp(Config.AOM_0_Attenuation), "MOT_AOM_0")
        play("MOT" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-")
        play("MOT" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+")
        play("AntiHelmholtz_MOT", "AntiHelmholtz_Coils")
    with for_(m, 1, m <= mot_repetitions, m + 1):
        measure("Detection", "FLR_detection", None, integration.full("Detection_opt", FLR, "out1"))
    #     # play("OD_FS" * amp(0.03), "AOM_2-3'_for_interference")

    align("Cooling_Sequence", "MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+", "AntiHelmholtz_Coils", "Zeeman_Coils",
          "AOM_2-2/3'", "AOM_2-3'_for_interference", "FLR_detection", "Measurement")  #, "PULSER_N", "AOM_S", "AOM_LO")

    return FLR


def MOT_no_Antihelmholtz(mot_repetitions):
    """
    The MOT function is used to play the MOT. To that end, we send RF signal to AOM_0, AOM_+, AOM_- for the duration of the MOT.

    Parameters
    ----------
    :param mot_repetitions: derived from the MOT duration, and is calculated in the Experiment parameters section
    """
    FLR_OFF = declare(fixed)

    align("Cooling_Sequence", "MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+", "AntiHelmholtz_Coils", "Zeeman_Coils",
          "AOM_2-2/3'", "AOM_2-2/3'_detuned", "FLR_detection", "Measurement", "Dig_detectors")  #, "PULSER_N", "AOM_S")

    ## MOT build-up ##
    n = declare(int)
    m = declare(int)
    with for_(n, 1, n <= mot_repetitions, n + 1):
        play("MOT" * amp(Config.AOM_0_Attenuation), "MOT_AOM_0")
        play("MOT" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-")
        play("MOT" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+")

    with for_(m, 1, m <= mot_repetitions, m + 1):
        measure("Detection", "FLR_detection", None, integration.full("Detection_opt", FLR_OFF, "out1"))
        play("OD_FS" * amp(0.05), "AOM_2-3'_for_interference")

    align("Cooling_Sequence", "MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+", "AntiHelmholtz_Coils", "Zeeman_Coils",
          "AOM_2-2/3'", "AOM_2-2/3'_detuned", "FLR_detection", "Measurement", "Dig_detectors")  #, "PULSER_N", "AOM_S")

    return FLR_OFF


def Pulse_const(total_pulse_duration, final_amp_0, final_amp_minus, final_amp_plus):
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
    play("Const" * amp(final_amp_0 * Config.AOM_0_Attenuation), "MOT_AOM_0", duration=total_pulse_duration)
    play("Const" * amp(final_amp_minus * Config.AOM_Minus_Attenuation), "MOT_AOM_-", duration=total_pulse_duration)
    play("Const" * amp(final_amp_plus * Config.AOM_Plus_Attenuation), "MOT_AOM_+", duration=total_pulse_duration)


def Pulse_with_prep(total_pulse_duration, prep_duration, initial_amp_0, initial_amp_minus, initial_amp_plus,
                    final_amp_0, final_amp_minus, final_amp_plus, zero_pulse_duration, plus_pulse_duration,
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
        play("Linear" * amp(initial_amp_0 * Config.AOM_0_Attenuation), "MOT_AOM_0", duration=zero_pulse_duration,
             truncate=prep_duration)
        play("Linear" * amp(initial_amp_minus * Config.AOM_Minus_Attenuation), "MOT_AOM_-",
             duration=minus_pulse_duration, truncate=prep_duration)
        play("Linear" * amp(initial_amp_plus * Config.AOM_Plus_Attenuation), "MOT_AOM_+", duration=plus_pulse_duration,
             truncate=prep_duration)

    ## Playing the pulses to the AOMs for the constant part. (Qua) ##
    align("MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+")
    # align("MOT_AOM_-", "MOT_AOM_+")
    with if_(total_pulse_duration > prep_duration):
        play("Const" * amp(final_amp_0 * Config.AOM_0_Attenuation), "MOT_AOM_0",
             duration=(total_pulse_duration - prep_duration))
        play("Const" * amp(final_amp_minus * Config.AOM_Minus_Attenuation), "MOT_AOM_-",
             duration=(total_pulse_duration - prep_duration))
        play("Const" * amp(final_amp_plus * Config.AOM_Plus_Attenuation), "MOT_AOM_+",
             duration=(total_pulse_duration - prep_duration))


def Pulse_with_prep_with_chirp(total_pulse_duration, prep_duration, initial_amp_0, initial_amp_minus, initial_amp_plus,
                               final_amp_0, final_amp_minus, final_amp_plus, zero_pulse_duration, plus_pulse_duration,
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
    # align("MOT_AOM_-", "MOT_AOM_+")
    with if_(prep_duration > 0):
        with if_(zero_pulse_duration == prep_duration):
            play("Const" * amp(initial_amp_0 * Config.AOM_0_Attenuation), "MOT_AOM_0", duration=prep_duration)
            play("Const" * amp(initial_amp_minus * Config.AOM_Minus_Attenuation), "MOT_AOM_-", duration=prep_duration,
                 chirp=(aom_chirp_rate, "mHz/nsec"))
            play("Const" * amp(initial_amp_plus * Config.AOM_Plus_Attenuation), "MOT_AOM_+", duration=prep_duration,
                 chirp=(-aom_chirp_rate, "mHz/nsec"))
        with else_():
            play("Linear" * amp(initial_amp_0 * Config.AOM_0_Attenuation), "MOT_AOM_0", duration=zero_pulse_duration,
                 truncate=prep_duration)
            play("Linear" * amp(initial_amp_minus * Config.AOM_Minus_Attenuation), "MOT_AOM_-",
                 duration=minus_pulse_duration, truncate=prep_duration,
                 chirp=(aom_chirp_rate, "mHz/nsec"))
            play("Linear" * amp(initial_amp_plus * Config.AOM_Plus_Attenuation), "MOT_AOM_+",
                 duration=plus_pulse_duration, truncate=prep_duration,
                 chirp=(-aom_chirp_rate, "mHz/nsec"))

    ## Playing the pulses to the AOMs for the constant part. (Qua) ##
    align("MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+")
    # align("MOT_AOM_-", "MOT_AOM_+")
    update_frequency("MOT_AOM_-", Config.IF_AOM_MOT + delta_f)
    update_frequency("MOT_AOM_+", Config.IF_AOM_MOT - delta_f)
    with if_(total_pulse_duration > prep_duration):
        play("Const" * amp(final_amp_0 * Config.AOM_0_Attenuation), "MOT_AOM_0",
             duration=(total_pulse_duration - prep_duration))
        play("Const" * amp(final_amp_minus * Config.AOM_Minus_Attenuation), "MOT_AOM_-",
             duration=(total_pulse_duration - prep_duration))
        play("Const" * amp(final_amp_plus * Config.AOM_Plus_Attenuation), "MOT_AOM_+",
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

    update_frequency("MOT_AOM_0", Config.IF_AOM_MOT)
    update_frequency("MOT_AOM_-", Config.IF_AOM_MOT)
    update_frequency("MOT_AOM_+", Config.IF_AOM_MOT)

    ## Aligning all the different elements used during the freefall time of the experiment ##
    align("Cooling_Sequence", "MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+", "Zeeman_Coils", "AOM_2-2/3'",
          "AOM_2-2/3'_detuned", "Measurement")  #, "PULSER_N", "AOM_S")

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

#
# def OD_Measure_SNSPDs(M_delay, OD_delay, rep, OD_rep_pulse1, OD_rep_pulse2, OD_sleep, m_time, m_window,
#                       ON_counts_st1, ON_counts_st2, ON_counts_st3, ON_counts_st4,
#                       ON_counts_st5, ON_counts_st6, ON_counts_st7, ON_counts_st8,
#                       OFF_counts_st1, OFF_counts_st2, OFF_counts_st3, OFF_counts_st4, OFF_counts_st5, OFF_counts_st6,
#                       OFF_counts_st7, OFF_counts_st8, OD_st,
#                       antihelmholtz_ON):
#     """
#     Measuring the OD with and without atoms.
#
#     Parameters
#     ----------
#     :param OD_delay: The time delay from end of PGC until the first OD pulse (OD + trigger).
#     :param rep: The number of measuring window repetitions, derived from OD measuring duration (M_time/M_window).
#     :param m_time: The duration of the entire measurement time.
#     :param m_window: The duration of each measuring window - fpr each window there is 28 nsec "deadtime".
#     :param counts_st1: The stream array for the number of photons counted at each measuring window for detector 1.
#     :param counts_st2: The stream array for the number of photons counted at each measuring window for detector 2.
#     :param counts_st3: The stream array for the number of photons counted at each measuring window for detector 3.
#     :param counts_st4: The stream array for the number of photons counted at each measuring window for detector 4.
#     :param counts_st5: The stream array for the number of photons counted at each measuring window for detector 5.
#     :param counts_st6: The stream array for the number of photons counted at each measuring window for detector 6.
#     :param counts_st7: The stream array for the number of photons counted at each measuring window for detector 7.
#     :param counts_st8: The stream array for the number of photons counted at each measuring window for detector 8.
#     """
#     counts1 = declare(int)
#     counts2 = declare(int)
#     counts3 = declare(int)
#     counts4 = declare(int)
#     counts5 = declare(int)
#     counts6 = declare(int)
#     counts7 = declare(int)
#     counts8 = declare(int)
#     OD = declare(fixed)
#     n = declare(int)
#     m = declare(int)
#     l = declare(int)
#
#     align("AOM_2-2/3'", "Dig_detectors")
#     wait(OD_delay + 4, "AOM_2-2/3'")
#     with for_(m, 0, m < OD_rep_pulse1, m + 1):
#         measure("OD", "AOM_2-2/3'", None, integration.full("Detection_opt", OD, "out1"))
#         save(OD, OD_st)
#     with if_(OD_rep_pulse2 > 0):
#         wait(OD_sleep, "AOM_2-2/3'")
#         with for_(l, 0, l < OD_rep_pulse2, l + 1):
#             measure("OD", "AOM_2-2/3'", None, integration.full("Detection_opt", OD, "out1"))
#             save(OD, OD_st)
#     wait(M_delay + 4, "Dig_detectors")
#     with for_(n, 0, n < rep, n + 1):
#         measure("readout", "Dig_detectors", None,
#                 counting.digital(counts1, m_window, element_outputs="out1"),
#                 counting.digital(counts2, m_window, element_outputs="out2"),
#                 counting.digital(counts3, m_window, element_outputs="out3"),
#                 counting.digital(counts4, m_window, element_outputs="out4"),
#                 counting.digital(counts5, m_window, element_outputs="out5"),
#                 counting.digital(counts6, m_window, element_outputs="out6"),
#                 counting.digital(counts7, m_window, element_outputs="out7"),
#                 counting.digital(counts8, m_window, element_outputs="out8"),
#                 )
#
#         ## Save Data: ##
#
#         ## Number of Photons (NOP) Count stream for each detector: ##
#         with if_(antihelmholtz_ON):
#             save(counts1, ON_counts_st1)
#             save(counts2, ON_counts_st2)
#             save(counts3, ON_counts_st3)
#             save(counts4, ON_counts_st4)
#             save(counts5, ON_counts_st5)
#             save(counts6, ON_counts_st6)
#             save(counts7, ON_counts_st7)
#             save(counts8, ON_counts_st8)
#         with else_():
#             save(counts1, OFF_counts_st1)
#             save(counts2, OFF_counts_st2)
#             save(counts3, OFF_counts_st3)
#             save(counts4, OFF_counts_st4)
#             save(counts5, OFF_counts_st5)
#             save(counts6, OFF_counts_st6)
#             save(counts7, OFF_counts_st7)
#             save(counts8, OFF_counts_st8)
#
#
# def Transits_Exp(M_delay, rep, m_time, m_window, antihelmholtz_ON,
#                  ON_counts_st1, ON_counts_st2, ON_counts_st3, ON_counts_st4,
#                  ON_counts_st5, ON_counts_st6, ON_counts_st7, ON_counts_st8,
#                  OFF_counts_st1, OFF_counts_st2, OFF_counts_st3, OFF_counts_st4,
#                  OFF_counts_st5, OFF_counts_st6, OFF_counts_st7, OFF_counts_st8):
#     """
#     Transit measurement with and without atoms.
#
#     Parameters
#     ----------
#     :param M_delay: The time delay from end of PGC until the first 2-3' pulse (OD + trigger).
#     :param rep: The number of measuring window repetitions, derived from OD measuring duration (M_time/M_window).
#     :param m_time: The duration of the entire measurement time.
#     :param m_window: The duration of each measuring window - fpr each window there is 28 nsec "deadtime".
#     :param counts_st1: The stream array for the number of photons counted at each measuring window for detector 1.
#     :param counts_st2: The stream array for the number of photons counted at each measuring window for detector 2.
#     :param counts_st3: The stream array for the number of photons counted at each measuring window for detector 3.
#     :param counts_st4: The stream array for the number of photons counted at each measuring window for detector 4.
#     :param counts_st5: The stream array for the number of photons counted at each measuring window for detector 5.
#     :param counts_st6: The stream array for the number of photons counted at each measuring window for detector 6.
#     :param counts_st7: The stream array for the number of photons counted at each measuring window for detector 7.
#     :param counts_st8: The stream array for the number of photons counted at each measuring window for detector 8.
#     """
#     counts1 = declare(int)
#     counts2 = declare(int)
#     counts3 = declare(int)
#     counts4 = declare(int)
#     counts5 = declare(int)
#     counts6 = declare(int)
#     counts7 = declare(int)
#     counts8 = declare(int)
#
#     n = declare(int)
#
#     align("AOM_2-2/3'", "Dig_detectors")
#     wait(M_delay + 4, "Dig_detectors", "AOM_2-2/3'")
#     play("OD", "AOM_2-2/3'", duration=m_time)
#     with for_(n, 0, n < rep, n + 1):
#         measure("readout", "Dig_detectors", None,
#                 counting.digital(counts1, m_window, element_outputs="out1"),
#                 counting.digital(counts2, m_window, element_outputs="out2"),
#                 counting.digital(counts3, m_window, element_outputs="out3"),
#                 counting.digital(counts4, m_window, element_outputs="out4"),
#                 counting.digital(counts5, m_window, element_outputs="out5"),
#                 counting.digital(counts6, m_window, element_outputs="out6"),
#                 counting.digital(counts7, m_window, element_outputs="out7"),
#                 counting.digital(counts8, m_window, element_outputs="out8"),
#                 )
#
#         ## Save Data: ##
#
#         ## Number of Photons (NOP) Count stream for each detector: ##
#         with if_(antihelmholtz_ON):
#             save(counts1, ON_counts_st1)
#             save(counts2, ON_counts_st2)
#             save(counts3, ON_counts_st3)
#             save(counts4, ON_counts_st4)
#             save(counts5, ON_counts_st5)
#             save(counts6, ON_counts_st6)
#             save(counts7, ON_counts_st7)
#             save(counts8, ON_counts_st8)
#         with else_():
#             save(counts1, OFF_counts_st1)
#             save(counts2, OFF_counts_st2)
#             save(counts3, OFF_counts_st3)
#             save(counts4, OFF_counts_st4)
#             save(counts5, OFF_counts_st5)
#             save(counts6, OFF_counts_st6)
#             save(counts7, OFF_counts_st7)
#             save(counts8, OFF_counts_st8)
#
#
# def Transits_Exp_TT(M_delay, rep, m_time, m_window, shutter_open_time,
#                     ON_counts_st1, ON_counts_st2, ON_counts_st3, ON_counts_st4,
#                     ON_counts_st5, ON_counts_st6, ON_counts_st7, ON_counts_st8,
#                     tt_st_N, tt_st_S, rep_st):
#     """
#      Transit measurement with and without atoms.
#
#      Parameters
#      ----------
#      :param M_delay: The time delay from end of PGC until the first 2-3' pulse (OD + trigger).
#      :param rep: The number of measuring window repetitions, derived from OD measuring duration (M_time/M_window).
#      :param m_time: The duration of the entire measurement time.
#      :param m_window: The duration of each measuring window - fpr each window there is 28 nsec "deadtime".
#      :param counts_st1: The stream array for the number of photons counted at each measuring window for detector 1.
#      :param counts_st2: The stream array for the number of photons counted at each measuring window for detector 2.
#      :param counts_st3: The stream array for the number of photons counted at each measuring window for detector 3.
#      :param counts_st4: The stream array for the number of photons counted at each measuring window for detector 4.
#      :param counts_st5: The stream array for the number of photons counted at each measuring window for detector 5.
#      :param counts_st6: The stream array for the number of photons counted at each measuring window for detector 6.
#      :param counts_st7: The stream array for the number of photons counted at each measuring window for detector 7.
#      :param counts_st8: The stream array for the number of photons counted at each measuring window for detector 8.
#      """
#     counts1 = declare(int)
#     counts2 = declare(int)
#     counts3 = declare(int)
#     counts4 = declare(int)
#     counts5 = declare(int)
#     counts6 = declare(int)
#     counts7 = declare(int)
#     counts8 = declare(int)
#
#     tt_vec1 = declare(int, size=1250)
#     tt_vec2 = declare(int, size=1250)
#     tt_vec3 = declare(int, size=1250)
#     tt_vec4 = declare(int, size=1250)
#     tt_vec5 = declare(int, size=1250)
#     tt_vec6 = declare(int, size=1250)
#     tt_vec7 = declare(int, size=1250)
#     tt_vec8 = declare(int, size=1250)
#
#     n = declare(int)
#     m = declare(int)
#     i = declare(int)
#     j = declare(int)
#     m1 = declare(int)
#     m2 = declare(int)
#     m3 = declare(int)
#     m4 = declare(int)
#     m5 = declare(int)
#     m6 = declare(int)
#     m7 = declare(int)
#     m8 = declare(int)
#     zero = declare(int, value=0)
#
#     align("AOM_2-2/3'", "Dig_detectors")
#     wait(M_delay - shutter_open_time, "Dig_detectors", "AOM_2-2/3'")
#     play("OD", "AOM_2-2/3'", duration=shutter_open_time)
#     align("AOM_2-2/3'", "Dig_detectors")
#     play("OD", "AOM_2-2/3'", duration=m_time)
#     with for_(n, 0, n < (m_time * 4), n + m_window):
#         measure("readout", "Dig_detectors", None,
#                 time_tagging.digital(tt_vec1, m_window, element_output="out1", targetLen=counts1),
#                 time_tagging.digital(tt_vec2, m_window, element_output="out2", targetLen=counts2),
#                 time_tagging.digital(tt_vec3, m_window, element_output="out3", targetLen=counts3),
#                 time_tagging.digital(tt_vec4, m_window, element_output="out4", targetLen=counts4),
#                 time_tagging.digital(tt_vec5, m_window, element_output="out5", targetLen=counts5),
#                 time_tagging.digital(tt_vec6, m_window, element_output="out6", targetLen=counts6),
#                 time_tagging.digital(tt_vec7, m_window, element_output="out7", targetLen=counts7),
#                 time_tagging.digital(tt_vec8, m_window, element_output="out8", targetLen=counts8),
#                 )
#
#         ## Save Data: ##
#
#         ## Number of Photons (NOP) Count stream for each detector: ##
#         save(counts1, ON_counts_st1)
#         save(counts2, ON_counts_st2)
#         save(counts3, ON_counts_st3)
#         save(counts4, ON_counts_st4)
#         save(counts5, ON_counts_st5)
#         save(counts6, ON_counts_st6)
#         save(counts7, ON_counts_st7)
#         save(counts8, ON_counts_st8)
#
#         with for_(m, 0, m1 < counts1, m1 + 1):
#             save(tt_vec1[m1], tt_st_N)
#         with for_(m2, 0, m2 < counts2, m2 + 1):
#             save(tt_vec2[m2], tt_st_N)
#         with for_(m3, 0, m3 < counts3, m3 + 1):
#             save(tt_vec3[m3], tt_st_N)
#         with for_(m4, 0, m4 < counts4, m4 + 1):
#             save(tt_vec4[m4], tt_st_N)
#         with for_(m5, 0, m5 < counts5, m5 + 1):
#             save(tt_vec5[m5], tt_st_S)
#         with for_(m6, 0, m6 < counts6, m6 + 1):
#             save(tt_vec6[m6], tt_st_S)
#         with for_(m7, 0, m7 < counts7, m7 + 1):
#             save(tt_vec7[m7], tt_st_S)
#         with for_(m8, 0, m8 < counts8, m8 + 1):
#             save(tt_vec8[m8], tt_st_S)
#         with for_(i, 0, i < (5000 - counts1 - counts2 - counts3 - counts4), i + 1):
#             save(zero, tt_st_N)
#         with for_(j, 0, j < (5000 - counts5 - counts6 - counts7 - counts8), j + 1):
#             save(zero, tt_st_S)
#         with for_(m, 0, m < 5000, m + 1):
#             save(n, rep_st)
#
# def Probe_counts_Measure_SNSPDs(M_delay, rep, m_time, m_window, shutter_open_time,
#                                 ON_counts_st1, ON_counts_st2, ON_counts_st3, ON_counts_st4,
#                                 ON_counts_st5, ON_counts_st6, ON_counts_st7, ON_counts_st8,):
#     """
#     Measuring the OD with and without atoms.
#
#     Parameters
#     ----------
#     :param M_delay: The time delay from end of PGC until the first Probe pulse.
#     :param rep: The number of measuring window repetitions, derived from OD measuring duration (M_time/M_window).
#     :param m_time: The duration of the entire measurement time.
#     :param m_window: The duration of each measuring window - fpr each window there is 28 nsec "deadtime".
#     :param counts_st1: The stream array for the number of photons counted at each measuring window for detector 1.
#     :param counts_st2: The stream array for the number of photons counted at each measuring window for detector 2.
#     :param counts_st3: The stream array for the number of photons counted at each measuring window for detector 3.
#     :param counts_st4: The stream array for the number of photons counted at each measuring window for detector 4.
#     :param counts_st5: The stream array for the number of photons counted at each measuring window for detector 5.
#     :param counts_st6: The stream array for the number of photons counted at each measuring window for detector 6.
#     :param counts_st7: The stream array for the number of photons counted at each measuring window for detector 7.
#     :param counts_st8: The stream array for the number of photons counted at each measuring window for detector 8.
#     """
#     counts1 = declare(int)
#     counts2 = declare(int)
#     counts3 = declare(int)
#     counts4 = declare(int)
#     counts5 = declare(int)
#     counts6 = declare(int)
#     counts7 = declare(int)
#     counts8 = declare(int)
#
#     n = declare(int)
#
#     align("AOM_2-2/3'", "Dig_detectors")
#     wait(M_delay - shutter_open_time, "Dig_detectors", "AOM_2-2/3'")
#     play("OD", "AOM_2-2/3'", duration=shutter_open_time)
#     align("AOM_2-2/3'", "Dig_detectors")
#     play("OD", "AOM_2-2/3'", duration=m_time)
#     with for_(n, 0, n < (m_time * 4), n + m_window):
#         measure("readout", "Dig_detectors", None,
#                 counting.digital(counts1, m_window, element_outputs="out1"),
#                 counting.digital(counts2, m_window, element_outputs="out2"),
#                 counting.digital(counts3, m_window, element_outputs="out3"),
#                 counting.digital(counts4, m_window, element_outputs="out4"),
#                 counting.digital(counts5, m_window, element_outputs="out5"),
#                 counting.digital(counts6, m_window, element_outputs="out6"),
#                 counting.digital(counts7, m_window, element_outputs="out7"),
#                 counting.digital(counts8, m_window, element_outputs="out8"),
#                 )
#
#         ## Save Data: ##
#
#         ## Number of Photons (NOP) Count stream for each detector: ##
#
#         save(counts1, ON_counts_st1)
#         save(counts2, ON_counts_st2)
#         save(counts3, ON_counts_st3)
#         save(counts4, ON_counts_st4)
#         save(counts5, ON_counts_st5)
#         save(counts6, ON_counts_st6)
#         save(counts7, ON_counts_st7)
#         save(counts8, ON_counts_st8)
#
#
#
# def Spectrum_Experiment(M_delay, rep, m_time, m_window, shutter_open_time,
#                         ON_counts_st1, ON_counts_st2, ON_counts_st3, ON_counts_st4,
#                         ON_counts_st5, ON_counts_st6, ON_counts_st7, ON_counts_st8,
#                         tt_st_N, tt_st_S, rep_st, p_duration):
#     """
#        Transit measurement with and without atoms.
#
#        Parameters
#        ----------
#        :param M_delay: The time delay from end of PGC until the first 2-3' pulse (OD + trigger).
#        :param rep: The number of measuring window repetitions, derived from OD measuring duration (M_time/M_window).
#        :param m_time: The duration of the entire measurement time.
#        :param m_window: The duration of each measuring window - fpr each window there is 28 nsec "deadtime".
#        :param counts_st1: The stream array for the number of photons counted at each measuring window for detector 1.
#        :param counts_st2: The stream array for the number of photons counted at each measuring window for detector 2.
#        :param counts_st3: The stream array for the number of photons counted at each measuring window for detector 3.
#        :param counts_st4: The stream array for the number of photons counted at each measuring window for detector 4.
#        :param counts_st5: The stream array for the number of photons counted at each measuring window for detector 5.
#        :param counts_st6: The stream array for the number of photons counted at each measuring window for detector 6.
#        :param counts_st7: The stream array for the number of photons counted at each measuring window for detector 7.
#        :param counts_st8: The stream array for the number of photons counted at each measuring window for detector 8.
#        """
#
#     counts5 = declare(int)
#     counts6 = declare(int)
#     counts7 = declare(int)
#     counts8 = declare(int)
#
#     tt_vec5 = declare(int, size=700)
#     tt_vec6 = declare(int, size=700)
#     tt_vec7 = declare(int, size=700)
#     tt_vec8 = declare(int, size=700)
#
#     n = declare(int)
#     m = declare(int)
#     j = declare(int)
#     m5 = declare(int)
#     m6 = declare(int)
#     m7 = declare(int)
#     m8 = declare(int)
#     t_1 = declare(int)
#     t_2 = declare(int)
#
#     zero = declare(int, value=0)
#
#     align("AOM_2-2/3'", "AOM_2-2/3'_detuned", "Dig_detectors_spectrum")
#     wait(M_delay - shutter_open_time, "Dig_detectors_spectrum", "AOM_2-2/3'", "AOM_2-2/3'_detuned")
#     # play("OD", "AOM_2-2/3'", duration=shutter_open_time)
#     align("AOM_2-2/3'", "Dig_detectors_spectrum", "AOM_2-2/3'_detuned")
#     with for_(t_1, 0, t_1 < (m_time * 4), t_1 + 2 * 60):
#         play("OD_Gaussian_Pulse", "AOM_2-2/3'")
#         wait(60, "AOM_2-2/3'", "AOM_2-2/3'_detuned")
#         play("OD_Gaussian_Pulse", "AOM_2-2/3'_detuned")
#
#     with for_(n, 0, n < (m_time * 4), n + m_window):
#         measure("readout", "Dig_detectors_spectrum", None,
#                 time_tagging.digital(tt_vec5, m_window, element_output="out5", targetLen=counts5),
#                 time_tagging.digital(tt_vec6, m_window, element_output="out6", targetLen=counts6),
#                 time_tagging.digital(tt_vec7, m_window, element_output="out7", targetLen=counts7),
#                 time_tagging.digital(tt_vec8, m_window, element_output="out8", targetLen=counts8))
#
#     ## Save Data: ##
#
#         ## Number of Photons (NOP) Count stream for each detector: ##
#         save(counts5, ON_counts_st5)
#         save(counts6, ON_counts_st6)
#         save(counts7, ON_counts_st7)
#         save(counts8, ON_counts_st8)
#
#         with for_(m5, 0, m5 < counts5, m5 + 1):
#             save(tt_vec5[m5], tt_st_S)
#         with for_(m6, 0, m6 < counts6, m6 + 1):
#             save(tt_vec6[m6], tt_st_S)
#         with for_(m7, 0, m7 < counts7, m7 + 1):
#             save(tt_vec7[m7], tt_st_S)
#         with for_(m8, 0, m8 < counts8, m8 + 1):
#             save(tt_vec8[m8], tt_st_S)
#         with for_(j, 0, j < (5000 - counts5 - counts6 - counts7 - counts8), j + 1):
#             save(zero, tt_st_S)
#         with for_(m, 0, m < 5000, m + 1):
#             save(n, rep_st)


# def Main_exp(calibration_time, pulse_length, rep):
#     """
#
#     """
#     NREP = 20
#     n = declare(int)
#     frac = declare(fixed, value=0.5 / NREP)
#
#     align("PULSER_N", "AOM_S", "AOM_LO", "Pulse_EOM")
#     play("Detection_pulses", "Pulse_EOM")
#
#     # First pulse, phase = 0
#     reset_frame('AOM_LO')
#     align('AOM_LO', 'AOM_S', 'AOM_N')
#     wait(66, "AOM_LO")
#     play('Homodyne_Pulse', 'AOM_LO')
#     play('Homodyne_Pulse', 'AOM_S')
#     play('Homodyne_Pulse', 'AOM_N')
#
#     # Second pulse, phase = 2*pi/4
#     frame_rotation_2pi(0.25, 'AOM_LO')
#     play('Homodyne_Pulse', 'AOM_LO')
#     play('Homodyne_Pulse', 'AOM_S')
#     play('Homodyne_Pulse', 'AOM_N')

    # we shouldn't use "for_" loops: https://qm-docs.qualang.io/guides/timing_in_qua#the-implicit-align
    # we should use "for" loops, like in the last example here: https://qm-docs.qualang.io/guides/features.html#switch-case
    # All other pulses, phase = pi/n * k


def opx_control(obj, qm):
    with program() as opx_control_prog:
        ## declaring program variables: ##
        i = declare(int)
        k = declare(int)
        j = declare(int)
        filler = declare(int, value=1)
        Trigger_Phase = declare(int, value=Phases_Names.index(obj.Exp_Values['Triggering_Phase']))
        Imaging_Phase = declare(int, value=Phases_Names.index(obj.Exp_Values['Imaging_Phase']))

        # Boolean variables:
        MOT_ON = declare(bool, value=True)
        AntiHelmholtz_delay_ON = declare(bool, value=False)
        AntiHelmholtz_ON = declare(bool, value=True)
        Linear_PGC_ON = declare(bool, value=True)
        Transits_Exp_ON = declare(bool, value=False)
        Probe_max_counts_Exp_ON = declare(bool, value=False)
        Spectrum_Exp_ON = declare(bool, value=False)
        MW_Spectroscopy_ON = declare(bool, value=False)
        Main_ON = declare(bool, value=False)

        # MOT variables
        MOT_Repetitions = declare(int, value=obj.Exp_Values['MOT_rep'])
        MOT_OFF_Repetitions = declare(int, value=int(obj.Exp_Values['MOT_rep']))
        post_MOT_delay = declare(int, value=int(obj.Exp_Values['Post_MOT_delay'] * 1e6 / 4))

        # AntiHelmholtz delay after MOT:
        antihelmholtz_delay = declare(int, value=int(obj.Exp_Values['AntiHelmholtz_delay'] * 1e6 / 4))

        # PGC variables:
        pgc_duration = declare(int, value=int(obj.pgc_duration))
        pgc_prep_time = declare(int, value=int(obj.pgc_prep_duration))
        pgc_pulse_duration_0 = declare(int, value=int(
            obj.pgc_pulse_duration_0))  # The relative duration to reach the desired amplitude
        pgc_pulse_duration_minus = declare(int, value=int(
            obj.pgc_pulse_duration_minus))  # The relative duration to reach the desired amplitude
        pgc_pulse_duration_plus = declare(int, value=int(
            obj.pgc_pulse_duration_plus))  # The relative duration to reach the desired amplitude
        pgc_initial_amp_0 = declare(fixed, value=obj.pgc_initial_amp_0)
        pgc_initial_amp_minus = declare(fixed, value=obj.pgc_initial_amp_minus)
        pgc_initial_amp_plus = declare(fixed, value=obj.pgc_initial_amp_plus)
        pgc_final_amp_0 = declare(fixed, value=obj.pgc_final_amp_0)
        pgc_final_amp_minus = declare(fixed, value=obj.pgc_final_amp_minus)
        pgc_final_amp_plus = declare(fixed, value=obj.pgc_final_amp_plus)
        # pgc_aom_chirp_rate = declare(int, value=int(
        #     obj.pgc_aom_chirp_rate))  # [mHz/nsec], If needed pgc preparation duration must be constant!!!

        # Fountain variables:
        fountain_duration = declare(int, value=int(obj.fountain_duration))
        # TODO: pre_PGC_fountain_duration = declare(int, value=int(obj.pre_PGC_fountain_duration))
        fountain_prep_time = declare(int, value=int(obj.fountain_prep_duration))  # Can't be used with Chirp!!!
        fountain_pulse_duration_0 = declare(int, value=int(obj.fountain_pulse_duration_0))  # The relative duration to reach the desired amplitude
        fountain_pulse_duration_minus = declare(int, value=int(obj.fountain_pulse_duration_minus))  # The relative duration to reach the desired amplitude
        fountain_pulse_duration_plus = declare(int, value=int(obj.fountain_pulse_duration_plus))  # The relative duration to reach the desired amplitude
        fountain_initial_amp_0 = declare(fixed, value=obj.fountain_initial_amp_0)
        fountain_initial_amp_minus = declare(fixed, value=obj.fountain_initial_amp_minus)
        fountain_initial_amp_plus = declare(fixed, value=obj.fountain_initial_amp_plus)
        fountain_final_amp_0 = declare(fixed, value=obj.fountain_final_amp_0)
        fountain_final_amp_minus = declare(fixed, value=obj.fountain_final_amp_minus)
        fountain_final_amp_plus = declare(fixed, value=obj.fountain_final_amp_plus)
        fountain_aom_chirp_rate = declare(int, value=int(obj.fountain_aom_chirp_rate))  # mHz/nsec
        fountain_delta_f = declare(int, value=int(obj.Exp_Values['Fountain_final_Delta_freq']))

        # Free fall duration:
        FreeFall_duration = declare(int, value=int(obj.Exp_Values['FreeFall_duration'] * 1e6 / 4))
        coils_timing = declare(int, value=int(obj.Exp_Values['Coil_timing'] * 1e6 / 4))

        # Imaging variables:
        N_Snaps = declare(int, value=obj.Exp_Values['N_Snaps'])
        Buffer_Cycles = declare(int, value=obj.Exp_Values['Buffer_Cycles'])

        # General measurements variables:
        Trigger_delay = declare(int, value=int(obj.Exp_Values['Trigger_delay'] * 1e6 / 4))
        PrePulse_duration = declare(int, value=int(obj.Exp_Values['PrePulse_duration'] * 1e6 / 4))

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
        pulse_1_initial_amp_0 = declare(fixed, value=obj.Exp_Values['Pulse_1_initial_amp_0'])
        pulse_1_initial_amp_minus = declare(fixed, value=obj.Exp_Values['Pulse_1_initial_amp_minus'])
        pulse_1_initial_amp_plus = declare(fixed, value=obj.Exp_Values['Pulse_1_initial_amp_plus'])
        pulse_1_final_amp_0 = declare(fixed, value=obj.Exp_Values['Pulse_1_final_amp_0'])
        pulse_1_final_amp_minus = declare(fixed, value=obj.Exp_Values['Pulse_1_final_amp_minus'])
        pulse_1_final_amp_plus = declare(fixed, value=obj.Exp_Values['Pulse_1_final_amp_plus'])

        # Depump measurement variables:
        Depump_pulse_duration = declare(int, value=int(obj.Depump_pulse_duration * 1e6 / 4))
        Depump_pulses_spacing = declare(int, value=int(obj.Depump_pulses_spacing * 1e6 / 4))
        Depump_start = declare(int, value=int(obj.Depump_Start * 1e6 / 4))

        # OD measurement variables:
        ## Free Space:
        OD_FS_pulse_duration = declare(int, value=int(obj.OD_FS_pulse_duration * 1e6 / 4))
        OD_FS_pulses_spacing = declare(int, value=int(obj.OD_FS_pulses_spacing * 1e6 / 4))
        OD_FS_start = declare(int, value=int(obj.OD_FS_Start * 1e6 / 4))
        OD_FS_sleep = declare(int, value=int(obj.OD_FS_sleep * 1e6 / 4))
        ## In fiber:
        M_Delay = declare(int, value=int(obj.M_delay * 1e6 / 4))
        M_duration = declare(int, value=int(obj.M_time / 4)) # From [nsec] to [4 nsec]
        OD_Delay = declare(int, value=int(obj.OD_delay * 1e6 / 4))
        OD_Rep_pulse1 = declare(int, value=int(obj.OD_duration_pulse1 * 1e6 / Config.OD_pulse_len))
        OD_Rep_pulse2 = declare(int, value=int(obj.OD_duration_pulse2 * 1e6 / Config.OD_pulse_len))
        OD_sleep = declare(int, value=int(obj.OD_sleep * 1e6 / 4))
        shutter_open_time = declare(int, value=int(obj.Shutter_open_time * 1e6 / 4))
        Rep = declare(int, value=obj.rep)

        # Spectrum:
        p_duration = declare(int, value=int(200/4))

        # MW spectroscopy variables:
        MW_freq = declare(int, value=obj.MW_start_frequency)
        pulse_length_MW = declare(int, value=(obj.Pulse_Length_MW * int(1e3 / 4)))
        pulse_length_OD = declare(int, value=(obj.Pulse_Length_OD * int(1e3 / 4)))
        Repetitions = declare(int, value=1)
        Delta_f = declare(int, value=0)


        # Stream processing:
        ON_counts_st1 = declare_stream()
        ON_counts_st2 = declare_stream()
        ON_counts_st3 = declare_stream()
        ON_counts_st4 = declare_stream()
        ON_counts_st5 = declare_stream()
        ON_counts_st6 = declare_stream()
        ON_counts_st7 = declare_stream()
        ON_counts_st8 = declare_stream()
        tt_st_N = declare_stream()
        tt_st_S = declare_stream()
        rep_st = declare_stream()
        AntiHelmholtz_ON_st = declare_stream()
        FLR_st = declare_stream()

        assign(IO1, 0)
        assign(IO2, 0)

        wait(100)

        with infinite_loop_():
            assign(i, IO1)
            with for_(k, 1, k <= N_Snaps + Buffer_Cycles, k + 1):
                with for_(j, 1, j <= filler, j + 1):
                    ##########################
                    ## Cooling Sequence ##
                    ##########################

                    with if_(MOT_ON & AntiHelmholtz_ON):
                        FLR = MOT(MOT_Repetitions)
                        save(FLR, FLR_st)
                        align(*all_elements)
                    with if_(antihelmholtz_delay > 0):
                        play("AntiHelmholtz_MOT", "AntiHelmholtz_Coils", duration=antihelmholtz_delay)
                    with if_(Trigger_Phase == 1):  # Trigger on PGC
                        ## Trigger QuadRF Sequence #####################
                        play("C_Seq", "Cooling_Sequence", duration=2500)
                        ################################################
                    # align(*all_elements, "AOM_2-2/3'")
                    with if_(post_MOT_delay > 0):
                        wait(post_MOT_delay, "Cooling_Sequence")
                        align(*all_elements, "Zeeman_Coils", "AOM_2-2/3'")
                    with if_(fountain_duration > 0):
                        wait(fountain_duration, "Cooling_Sequence")
                        Pulse_with_prep_with_chirp(fountain_duration, obj.fountain_prep_duration,
                                                   fountain_initial_amp_0,
                                                   fountain_initial_amp_minus, fountain_initial_amp_plus,
                                                   fountain_final_amp_0,
                                                   fountain_final_amp_minus, fountain_final_amp_plus,
                                                   fountain_pulse_duration_0,
                                                   fountain_pulse_duration_minus, fountain_pulse_duration_plus,
                                                   fountain_aom_chirp_rate, fountain_delta_f)
                        align(*all_elements, "Zeeman_Coils", "AOM_2-2/3'")
                    with if_(pgc_duration > 0):
                        wait(pgc_duration, "Cooling_Sequence")
                        with if_(pgc_initial_amp_minus == pgc_final_amp_minus):
                            Pulse_const(pgc_duration, pgc_final_amp_0, pgc_final_amp_minus, pgc_final_amp_plus)
                        with else_():
                            Pulse_with_prep(pgc_duration, pgc_prep_time, pgc_initial_amp_0, pgc_initial_amp_minus,
                                            pgc_initial_amp_plus, pgc_final_amp_0, pgc_final_amp_minus, pgc_final_amp_plus,
                                            pgc_pulse_duration_0, pgc_pulse_duration_minus, pgc_pulse_duration_plus)
                        align(*all_elements, "Zeeman_Coils", "AOM_2-2/3'")
                    with if_((k >= 1) & (Imaging_Phase == 4)):  # 4 means imaging phase on pulse_1
                        with if_(k <= Buffer_Cycles):
                            FreeFall(FreeFall_duration, coils_timing)
                        with else_():
                            FreeFall(FreeFall_duration + PrePulse_duration * (k - 1 - Buffer_Cycles), coils_timing)

                    ##########################
                    ## Measurement Sequence ##
                    ##########################

                    ## Measurement start time:
                    with if_((k > Buffer_Cycles) & (Imaging_Phase == 4)):  # 4 means imaging phase on pulse_1
                        wait(PrePulse_duration * (k - Buffer_Cycles), "Cooling_Sequence", "Measurement")
                    with else_():
                        wait(PrePulse_duration, "Cooling_Sequence", "Measurement")
                    align(*all_elements)

                    with if_(Trigger_Phase == 4):  # when trigger on pulse 1
                        ## Trigger QuadRF Sequence #####################
                        play("C_Seq", "Cooling_Sequence", duration=250)
                        ################################################
                        with if_((Imaging_Phase == 4) & (Pulse_1_duration > 0)):  # 4 means imaging phase on pulse_1
                            ## For Depump measurement:
                            with if_(Depump_pulse_duration > 0):
                                wait(Depump_start + 4, "AOM_2-2/3'")
                                Depump_Measure(Depump_pulse_duration, Depump_pulses_spacing)
                                align(*all_elements, "AOM_2-2/3'")
                            with else_():
                                ## For free-space OD measurement:
                                with if_(OD_FS_pulse_duration > 0):
                                    wait(OD_FS_start + 4, "AOM_2-2/3'")
                                    OD_Measure(OD_FS_pulse_duration, OD_FS_pulses_spacing, OD_FS_sleep)
                                    align(*all_elements, "AOM_2-2/3'")
                                with else_():
                                    align(*all_elements, "AOM_2-2/3'")
                                    Pulse_with_prep(Pulse_1_duration, Pulse_1_decay_time, pulse_1_initial_amp_0,
                                                    pulse_1_initial_amp_minus, pulse_1_initial_amp_plus,
                                                    pulse_1_final_amp_0, pulse_1_final_amp_minus,
                                                    pulse_1_final_amp_plus, pulse_1_duration_0,
                                                    pulse_1_duration_minus, pulse_1_duration_plus)
                                    # play("OD", "AOM_2-2/3'", duration=Pulse_1_duration)
                                    Measure(Pulse_1_duration)  # This triggers camera (Control 7)
                                    align(*all_elements, "AOM_2-2/3'")
                        align(*all_elements, "AOM_2-2/3'")

            with if_(~(AntiHelmholtz_ON | Transits_Exp_ON)):
                assign(AntiHelmholtz_ON, True)
                # assign(IO2, False)

            assign(N_Snaps, 1)
            assign(Buffer_Cycles, 0)
            assign(i, IO1)

            ## PARAMETERS UPDATE ##
            with if_(i > 0):
                pause()
            with while_(i > 0):
                ## Boolean variables control: ##
                with if_(i == 1):
                    assign(MOT_ON, IO2)
                with if_(i == 2):
                    assign(Linear_PGC_ON, IO2)
                with if_(i == 3):
                    assign(Transits_Exp_ON, IO2)
                with if_(i == 4):
                    assign(Spectrum_Exp_ON, IO2)
                with if_(i == 5):
                    assign(AntiHelmholtz_delay_ON, IO2)
                with if_(i == 6):
                    assign(Probe_max_counts_Exp_ON, IO2)
                with if_(i == 7):
                    assign(filler, IO2)
                ## AntiHelmholtz control ##
                with if_(i == 10):
                    assign(antihelmholtz_delay, IO2)
                ## MOT variables control ##
                with if_(i == 11):  # Live control over the
                    assign(MOT_Repetitions, IO2)
                with if_(i == 12):  # Live control over the
                    assign(post_MOT_delay, IO2)

                ## PGC variables control ##
                with if_(i == 20):  # Live control over the PGC duration
                    assign(pgc_duration, IO2)
                with if_(i == 21):  # Live control over the preparation time of the PGC
                    assign(pgc_prep_time, IO2)
                with if_(i == 22):  # Live control over the initial amplitude of the PGC AOM 0
                    assign(pgc_initial_amp_0, IO2)
                with if_(i == 23):  # Live control over the initial amplitude of the PGC AOM -
                    assign(pgc_initial_amp_minus, IO2)
                with if_(i == 24):  # Live control over the initial amplitude of the PGC AOM +
                    assign(pgc_initial_amp_plus, IO2)
                with if_(i == 25):  # Live control over the final amplitude of the PGC AOM 0
                    assign(pgc_final_amp_0, IO2)
                with if_(i == 26):  # Live control over the final amplitude of the PGC AOM -
                    assign(pgc_final_amp_minus, IO2)
                with if_(i == 27):  # Live control over the final amplitude of the PGC AOM +
                    assign(pgc_final_amp_plus, IO2)
                with if_(i == 28):  # Live control over the final amplitude of the PGC AOM 0
                    assign(pgc_pulse_duration_0, IO2)
                with if_(i == 29):  # Live control over the final amplitude of the PGC AOM -
                    assign(pgc_pulse_duration_minus, IO2)
                with if_(i == 30):  # Live control over the final amplitude of the PGC AOM +
                    assign(pgc_pulse_duration_plus, IO2)

                ## Fountain variables control: ##
                with if_(i == 31):  # Live control over the fountain duration
                    assign(fountain_duration, IO2)
                with if_(i == 32):  # Live control over the preparation time of the Fountain
                    assign(fountain_prep_time, IO2)
                # TODO: note this is the right place for changing fountain initial amps. However, index order is wrong.
                # with if_(i == 73):  # Live control over the fountain duration
                #     assign(pre_PGC_fountain_duration, IO2)
                # however this should work.
                with if_(i == 70):  # Live control over the initial amplitude of the fountain AOM 0
                    assign(fountain_initial_amp_0, IO2)
                ## Fountain variables control: ##
                with if_(i == 71):  # Live control over the initial amplitude of the fountain AOM -
                    assign(fountain_initial_amp_minus, IO2)
                with if_(i == 72):  # Live control over the initial amplitude of the fountain AOM +
                    assign(fountain_initial_amp_plus, IO2)

                with if_(i == 33):  # Live control over the final amplitude of the fountain AOM 0
                    assign(fountain_final_amp_0, IO2)
                with if_(i == 34):  # Live control over the final amplitude of the fountain AOM -
                    assign(fountain_final_amp_minus, IO2)
                with if_(i == 35):  # Live control over the final amplitude of the fountain AOM +
                    assign(fountain_final_amp_plus, IO2)
                with if_(i == 36):  # Live control over the final amplitude of the fountain AOM 0
                    assign(fountain_pulse_duration_0, IO2)
                with if_(i == 37):  # Live control over the final amplitude of the fountain AOM -
                    assign(fountain_pulse_duration_minus, IO2)
                with if_(i == 38):  # Live control over the final amplitude of the fountain AOM +
                    assign(fountain_pulse_duration_plus, IO2)
                with if_(i == 39):  # Live control over the fountain final frequency of the MOT + and - AOMs
                    assign(fountain_aom_chirp_rate, IO2)

                ## Measurement variables control: ##
                with if_(i == 41):
                    assign(Trigger_delay, IO2)
                with if_(i == 42):
                    assign(PrePulse_duration, IO2)
                with if_(i == 43):
                    assign(Pulse_1_duration, IO2)
                with if_(i == 44):
                    assign(Pulse_1_decay_time, IO2)
                with if_(i == 45):  # The number of frames(snaps) to take for a video with growing snap time intervals
                    assign(N_Snaps, IO2)
                with if_(i == 46):
                    assign(Buffer_Cycles, IO2)

                ## OD and N_atoms measuring variable control ##
                with if_(i == 51):  # Live control of the Free Space OD measurement start time
                    assign(OD_FS_start, IO2)
                with if_(i == 52):  # Live control of the Free Space OD measurement pulses duration
                    assign(OD_FS_pulse_duration, IO2)
                with if_(i == 53):  # Live control of the Free Space OD measurement wait duration between the 2 pulses
                    assign(OD_FS_pulses_spacing, IO2)
                with if_(i == 54):  # Live control of the SNSPDs OD measurement start time
                    assign(OD_Delay, IO2)
                with if_(i == 55):  # Live control of the SNSPDs OD measurement duration
                    assign(Rep, IO2)
                with if_(i == 56):  # Live control of the Depump measurement start time
                    assign(Depump_start, IO2)
                with if_(i == 57):  # Live control of the Depump measurement pulses duration
                    assign(Depump_pulse_duration, IO2)
                with if_(i == 58):  # Live control of the Depump measurement wait duration between the 2 pulses
                    assign(Depump_pulses_spacing, IO2)
                with if_(i == 59):  # Live control of the delay due to shutter opening time.
                    assign(shutter_open_time, IO2)
                ## MW spectroscopy control ##
                with if_(i == 61):  # Live control of the MW spectroscopy MW frequency
                    assign(MW_freq, IO2)
                with if_(i == 62):  # Live control of the MW spectroscopy MW pulse length
                    assign(pulse_length_MW, IO2)
                with if_(i == 63):  # Live control of the MW spectroscopy OD pulse length
                    assign(pulse_length_OD, IO2)
                with if_(i == 64):  # Live control of the MW spectroscopy MW pulse repetitions
                    assign(Repetitions, IO2)
                with if_(i == 65):  # Live control of the MW spectroscopy MW pulse frequency change
                    assign(Delta_f, IO2)

                pause()
                assign(i, IO1)

        with stream_processing():
            # (ON_counts_st1 + ON_counts_st2 + ON_counts_st3 + ON_counts_st4).buffer(obj.rep).save('North_Probe')
            # (ON_counts_st5 + ON_counts_st6 + ON_counts_st7 + ON_counts_st8).buffer(obj.rep).save('South_Probe')
            # (tt_st_N + rep_st).buffer(5000 * obj.rep).save('North_Probe_TT')
            # (tt_st_S + rep_st).buffer(5000 * obj.rep).save('South_Probe_TT')
            FLR_st.save('FLR_measure')
            # AntiHelmholtz_ON_st.save("antihelmholtz_on")

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

        self.Exp_Values = Initial_Values  # Initialize experiment values to be as in Config_Table.py
        QuadRFMOTController(initialValues=self.Exp_Values, updateChannels=(1, 4), topticaLockWhenUpdating = False, debugging=True,
                            continuous=False)  # updates values on QuadRF (uploads table)
        # QuadRFMOTController(initialValues={'Operation_Mode': 'Continuous',  'CH3_freq': '133MHz', 'CH3_amp': '31dbm'},
        #                     updateChannels=[3], debugging=False, continuous=False)  # updates values on QuadRF (uploads table)
        QuadRFMOTController(initialValues={'Operation_Mode': 'Continuous', 'CH2_freq': '113MHz', 'CH2_amp': '15.25dbm'},
                            updateChannels=[2], debugging=True, continuous=False)  # updates values on QuadRF (uploads table)

        self.Update_QuadRF_channels = set({})  # Only update these channels on QuadRF when UpdateParameters method is called [note: this is a python set]

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
            self.pgc_pulse_duration_0 = int((self.Exp_Values['PGC_initial_amp_0'] * self.Exp_Values['PGC_prep_duration'] / (self.Exp_Values['PGC_initial_amp_0'] - self.Exp_Values['PGC_final_amp_0'])) * 1e6 / 4)  # The relative duration to reach the desired amplitude
            self.pgc_pulse_duration_minus = int((self.Exp_Values['PGC_initial_amp_minus'] * self.Exp_Values['PGC_prep_duration'] / (self.Exp_Values['PGC_initial_amp_minus'] - self.Exp_Values['PGC_final_amp_minus'])) * 1e6 / 4)  # The relative duration to reach the desired amplitude
            self.pgc_pulse_duration_plus = int((self.Exp_Values['PGC_initial_amp_plus'] * self.Exp_Values['PGC_prep_duration'] / (self.Exp_Values['PGC_initial_amp_plus'] - self.Exp_Values['PGC_final_amp_plus'])) * 1e6 / 4)  # The relative duration to reach the desired amplitude
            if self.pgc_pulse_duration_0 > int(60 * 1e6 / 4):  # longer then 60ms
                printError('The values for PGC_initial_amp_0 and PGC_final_amp_0 are too close or PGC_prep_duration is too long, might cause an ERROR!!!')
            if self.pgc_pulse_duration_minus > int(60 * 1e6 / 4):  # longer then 60ms
                printError('The values for PGC_initial_amp_minus and PGC_final_amp_minus are too close or PGC_prep_duration is too long, might cause an ERROR!!!')
            if self.pgc_pulse_duration_plus > int(60 * 1e6 / 4):  # longer then 60ms
                printError('The values for PGC_initial_amp_plus and PGC_final_amp_plus are too close or PGC_prep_duration is too long, might cause an ERROR!!!')
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
            self.fountain_pulse_duration_minus = int(self.Exp_Values['Fountain_initial_amp_minus'] * self.Exp_Values['Fountain_prep_duration'] * 1e6 / 4 / (
                        self.Exp_Values['Fountain_initial_amp_minus'] - self.Exp_Values['Fountain_final_amp_minus']))  # The relative duration to reach the desired amplitude
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
        self.fountain_aom_chirp_rate = int(self.Exp_Values['Fountain_final_Delta_freq'] * 1e3 / (self.Exp_Values['Fountain_prep_duration'] * 1e6))  # mHz/nsec

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
        self.M_delay = self.Exp_Values['M_delay']  # [msec]
        self.OD_delay = self.Exp_Values['OD_delay']  # [msec]
        self.M_window = self.Exp_Values['M_window']  # [nsec]
        self.OD_duration_pulse1 = self.Exp_Values['OD_duration_pulse1']  # [msec]
        self.OD_sleep = self.Exp_Values['OD_sleep']  # [msec]
        self.OD_duration_pulse2 = self.Exp_Values['OD_duration_pulse2']  # [msec]
        self.M_time = int(self.Exp_Values['M_time'] * 1e6) # [nsec]
        self.Shutter_open_time = 3  # [msec]
        # self.rep = int(self.M_time / (self.M_window + 28 + 170))
        self.rep = int(self.M_time / self.M_window)

        # MW spectroscopy parameters:
        self.MW_start_frequency = int(100e6)  # [Hz]
        self.Pulse_Length_MW = 400  # [usec]
        self.Pulse_Length_OD = 20  # [usec]

        # Main Experiment:
        self.TOP2_pulse_len = int(Config.Probe_pulse_len / 4)  # [nsec]
        self.Calibration_time = 10  # [msec]

        self.qmm = QuantumMachinesManager(host='132.77.54.230', port='80')
        self.qmm.clear_all_job_results()
        self.qm = self.qmm.open_qm(config)
        self.job = opx_control(self, self.qm)

        self.io1_list = []
        self.io2_list = []

        self.keyPress = None  # initing key listener

    def __del__(self):
        # self.job.halt()
        self.qm.close()
        self.qmm.close()

    def updateValue(self, key, value, update_parameters = False):
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

    def saveConfigTable(self, path = '.\\Config_Table Archive\\'):
        timeStamp = time.strftime("%d-%m-%Y %H_%M_%S", time.localtime())
        helmCoilsState = 'Working on it'#self.getHelmholtzCoilsState()
        expValuesToSave = dict(self.Exp_Values)
        if 'opticalPowerCalibrationFunction' in expValuesToSave:
            expValuesToSave['opticalPowerCalibrationFunction'] = 'Experiment values contained power calibration function. It is, however, impossible to save a function as JSON. Thus, this line is added as a warning'
        saveData = {'Time Stamp': timeStamp, 'Exp_Values': expValuesToSave, 'Helmholtz Coils': helmCoilsState}
        with open(f'{path}{timeStamp}.json', 'w') as file:
            json.dump(saveData, file, indent=4)

    def getHelmholtzCoilsState(self):
        printGreen('Getting Helmholtz coils state...')
        res = {}
        try:
            HCoilsController = HMP4040Visa()
            for ch in [1,2,3]:
                HCoilsController.setOutput(ch = ch)
                state = HCoilsController.getOutputState()
                current = HCoilsController.getCurrent()
                voltage = HCoilsController.getVoltage()
                res['Channel %s' %str(ch)] = {'Current':current, 'Voltage':voltage, 'state': state}
            res['Output state'] = HCoilsController.getGeneralOutputState()
            HCoilsController = None # closing controller
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

    def MOT_switch(self, Bool):
        self.update_io_parameter(1, Bool)

    def Linear_PGC_switch(self, Bool):
        self.update_io_parameter(2, Bool)

    def Transit_Exp_switch(self, Bool):
        self.Transit_switch = Bool
        self.update_io_parameter(3, Bool)

    def Spectrum_Exp_switch(self, Bool):
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
            self.update_io_parameter(46, 3) # Update 3 buffer cycles
        self.update_io_parameter(45, N_Snaps) # Update N_snaps

    def update_PGC_prep_time(self, prep_time):  # In units of [msec]
        self.pgc_prep_duration = int(prep_time * 1e6 / 4)
        self.update_io_parameter(21, int(prep_time * 1e6 / 4))

        # self.update_PGC_AOM_0_final_amplitude(self.pgc_final_amp_0)
        # self.update_PGC_AOM_minus_final_amplitude(self.pgc_final_amp_minus)
        # self.update_PGC_AOM_plus_final_amplitude(self.pgc_final_amp_plus)

    def update_PGC_AOMs_initial_amplitude(self, pgc_initial_amp):
        self.update_PGC_AOM_0_initial_amplitude(pgc_initial_amp)
        self.update_PGC_AOM_minus_initial_amplitude(pgc_initial_amp)
        self.update_PGC_AOM_plus_initial_amplitude(pgc_initial_amp)

    def update_PGC_AOM_0_initial_amplitude(self, pgc_initial_amp_0):
        self.pgc_initial_amp_0 = pgc_initial_amp_0
        self.update_io_parameter(22, float(pgc_initial_amp_0))

    def update_PGC_AOM_plus_initial_amplitude(self, pgc_initial_amp_minus):
        self.pgc_initial_amp_minus = pgc_initial_amp_minus
        self.update_io_parameter(23, float(pgc_initial_amp_minus))

    def update_PGC_AOM_minus_initial_amplitude(self, pgc_initial_amp_plus):
        self.pgc_initial_amp_plus = pgc_initial_amp_plus
        self.update_io_parameter(24, float(pgc_initial_amp_plus))

    def update_PGC_AOMs_final_amplitude(self, pgc_final_amp):
        self.update_PGC_AOM_0_final_amplitude(pgc_final_amp)
        self.update_PGC_AOM_minus_final_amplitude(pgc_final_amp)
        self.update_PGC_AOM_plus_final_amplitude(pgc_final_amp)

    def update_PGC_AOM_0_final_amplitude(self, pgc_final_amp_0):
        self.pgc_final_amp_0 = pgc_final_amp_0
        self.update_io_parameter(25, float(pgc_final_amp_0))
        # change pgc pulse duration
        # The relative duration to reach the desired amplitude, for all AOMs (0, +, -) this keeps the relative detuning constant
        self.update_io_parameter(28, int(
            self.pgc_prep_duration * self.pgc_initial_amp_0 / (self.pgc_initial_amp_0 - pgc_final_amp_0)))

    def update_PGC_AOM_minus_final_amplitude(self, pgc_final_amp_minus):
        self.pgc_final_amp_minus = pgc_final_amp_minus
        self.update_io_parameter(26, float(pgc_final_amp_minus))
        # change pgc_pulse_duration_minus
        # The relative duration to reach the desired amplitude, for all AOMs (0, +, -) this keeps the relative detuning constant
        self.update_io_parameter(29, int(
            self.pgc_prep_duration * self.pgc_initial_amp_minus / (self.pgc_initial_amp_minus - pgc_final_amp_minus)))

    def update_PGC_AOM_plus_final_amplitude(self, pgc_final_amp_plus):
        self.pgc_final_amp_plus = pgc_final_amp_plus
        self.update_io_parameter(27, float(pgc_final_amp_plus))
        # change pgc_pulse_duration_plus
        # The relative duration to reach the desired amplitude, for all AOMs (0, +, -) this keeps the relative detuning constant
        self.update_io_parameter(30, int(
            self.pgc_prep_duration * self.pgc_initial_amp_plus / (self.pgc_initial_amp_plus - pgc_final_amp_plus)))

    ## Fountain variable update functions: ##

    def update_fountain_prep_time(self, prep_time):  # In units of [msec]
        self.fountain_prep_duration = int(prep_time * 1e6 / 4)
        self.update_io_parameter(32, int(prep_time * 1e6 / 4))
        # self.update_fountain_AOM_0_final_amplitude(self.fountain_final_amp_0)
        # self.update_fountain_AOM_minus_final_amplitude(self.fountain_final_amp_minus)
        # self.update_fountain_AOM_plus_final_amplitude(self.fountain_final_amp_plus)

    def update_fountain_AOMs_final_amplitude(self, fountain_final_amp):
        self.update_fountain_AOM_0_final_amplitude(fountain_final_amp)
        self.update_fountain_AOM_minus_final_amplitude(fountain_final_amp)
        self.update_fountain_AOM_plus_final_amplitude(fountain_final_amp)

    def update_fountain_AOM_0_final_amplitude(self, fountain_final_amp_0):
        self.fountain_final_amp_0 = fountain_final_amp_0
        self.update_io_parameter(33, float(fountain_final_amp_0))
        # The relative duration to reach the desired amplitude, for all AOMs (0, +, -) this keeps the relative detuning constant
        self.update_io_parameter(36, int(self.fountain_prep_duration * self.fountain_initial_amp_0 / (
                self.fountain_initial_amp_0 - fountain_final_amp_0)))

    def update_fountain_AOM_minus_final_amplitude(self, fountain_final_amp_minus):
        self.fountain_final_amp_minus = fountain_final_amp_minus
        self.update_io_parameter(34, float(fountain_final_amp_minus))
        # The relative duration to reach the desired amplitude, for all AOMs (0, +, -) this keeps the relative detuning constant
        self.update_io_parameter(37, int(self.fountain_prep_duration * self.fountain_initial_amp_minus / (
                self.fountain_initial_amp_minus - fountain_final_amp_minus)))

    def update_fountain_AOM_plus_final_amplitude(self, fountain_final_amp_plus):
        self.fountain_final_amp_plus = fountain_final_amp_plus
        self.update_io_parameter(35, float(fountain_final_amp_plus))
        # The relative duration to reach the desired amplitude, for all AOMs (0, +, -) this keeps the relative detuning constant
        self.update_io_parameter(38, int(self.fountain_prep_duration * self.fountain_initial_amp_plus / (
                self.fountain_initial_amp_plus - fountain_final_amp_plus)))

    def update_fountain_AOM_0_initial_amplitude(self, fountain_initial_amp_0):
        self.fountain_initial_amp_0 = fountain_initial_amp_0
        self.update_io_parameter(70, float(fountain_initial_amp_0))
        # The relative duration to reach the desired amplitude, for all AOMs (0, +, -) this keeps the relative detuning constant
        self.update_io_parameter(36, int(self.fountain_prep_duration * self.fountain_initial_amp_0 / (
                fountain_initial_amp_0 - self.fountain_final_amp_0)))

    def update_fountain_AOM_minus_initial_amplitude(self, fountain_initial_amp_minus):
        self.fountain_initial_amp_minus = fountain_initial_amp_minus
        self.update_io_parameter(71, float(fountain_initial_amp_minus))
        # The relative duration to reach the desired amplitude, for all AOMs (0, +, -) this keeps the relative detuning constant
        self.update_io_parameter(37, int(self.fountain_prep_duration * self.fountain_initial_amp_minus / (
                fountain_initial_amp_minus - self.fountain_final_amp_minus)))

    def update_fountain_AOM_plus_initial_amplitude(self, fountain_initial_amp_plus):
        self.fountain_initial_amp_plus = fountain_initial_amp_plus
        self.update_io_parameter(72, float(fountain_initial_amp_plus))
        # The relative duration to reach the desired amplitude, for all AOMs (0, +, -) this keeps the relative detuning constant
        self.update_io_parameter(38, int(self.fountain_prep_duration * self.fountain_initial_amp_plus / (
                fountain_initial_amp_plus - self.fountain_final_amp_plus)))

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

    def Save_SNSPDs_Measurement(self, N, bin_size, preComment, counts_threshold):
        # if preComment is True:
        if not preComment:
            preComment = pymsgbox.prompt('Add comment to measurement: ', default='', timeout=int(30e3))
        aftComment = None

        ### fetching data from server
        ### saving to file
        ###
        # avg_OD_N_handle = self.job.result_handles.get("North_OD_avg")
        # avg_OD_S_handle = self.job.result_handles.get("South_OD_avg")
        OD_N_handle = self.job.result_handles.get("North_Probe")
        OD_S_handle = self.job.result_handles.get("South_Probe")
        FLR_handle = self.job.result_handles.get("FLR_measure")

        # avg_OD_N_no_MOT_handle = self.job.result_handles.get("North_OD_avg_no_MOT")
        # avg_OD_S_no_MOT_handle = self.job.result_handles.get("South_OD_avg_no_MOT")
        OD_N_no_MOT_handle = self.job.result_handles.get("North_OD_no_MOT")
        OD_S_no_MOT_handle = self.job.result_handles.get("South_OD_no_MOT")
        FLR_no_MOT_handle = self.job.result_handles.get("FLR_measure_no_MOT")

        antihelmholtz_handle = self.job.result_handles.get("antihelmholtz_on")

        # avg_OD_N_handle.wait_for_values(1)
        # avg_OD_S_handle.wait_for_values(1)
        OD_N_handle.wait_for_values(1)
        OD_S_handle.wait_for_values(1)
        FLR_handle.wait_for_values(1)

        antihelmholtz_handle.wait_for_values(1)

        t = np.linspace(0, self.M_time, self.rep)
        t_shutter_open = 3  # [msec]
        shutter_open_indx = find_nearest(t, value=t_shutter_open)

        pass_threshold = False

        FLR_measurement = []
        FLR_no_MOT_measurement = []
        N_OD_measurement = []
        N_OD_no_MOT_measurement = []
        S_OD_measurement = []
        S_OD_no_MOT_measurement = []

        # avg_OD_N_res = avg_OD_N_handle.fetch_all()
        # avg_OD_S_res = avg_OD_S_handle.fetch_all()
        OD_N_res = OD_N_handle.fetch_all()
        OD_S_res = OD_S_handle.fetch_all()
        FLR_res = -FLR_handle.fetch_all()
        FLR_measurement.append(FLR_res.tolist())

        # avg_OD_N_no_MOT_res = avg_OD_N_res  # Wrong data. using as place-holder
        # avg_OD_S_no_MOT_res = avg_OD_S_res  # Wrong data. using as place-holder
        OD_N_no_MOT_res = OD_N_res * 0  # Wrong data. using as place-holder
        OD_S_no_MOT_res = OD_S_res * 0  # Wrong data. using as place-holder
        FLR_no_MOT_res = []  # Wrong data. using as place-holder

        ## Prepare plots
        fig = plt.figure()
        ax1 = plt.subplot2grid((5, 2), (0, 0), colspan=1, rowspan=2)
        ax2 = plt.subplot2grid((5, 2), (0, 1), colspan=1, rowspan=2)
        ax3 = plt.subplot2grid((5, 2), (2, 0), colspan=1, rowspan=2)
        ax4 = plt.subplot2grid((5, 2), (2, 1), colspan=1, rowspan=2)
        ax5 = plt.subplot2grid((5, 2), (4, 0), colspan=2, rowspan=1)

        ## record time
        timest = time.strftime("%Y%m%d-%H%M%S")
        datest = time.strftime("%Y%m%d")

        ## Listen for keyboard
        listener = keyboard.Listener(on_press=self.on_key_press)
        listener.start()  # start to listen on a separate thread
        self.keyPress = None
        print('\033[94m' + 'Press ESC to stop measurement.' + '\033[0m')  # print blue
        reps = 1  # a counter, number of repeats actually made.

        # Place holders for results
        OD_N_res_batch = [OD_N_res.tolist()]
        OD_N_no_MOT_res_batch = []
        OD_S_res_batch = [OD_S_res.tolist()]
        OD_S_no_MOT_res_batch = []

        # Photon counts per sec/msec converters:

        Per_sec = 1e9 / Config.readout_pulse_len
        Per_usec = 1e3 / Config.readout_pulse_len

        # while OD_N_handle.is_processing():  # and self.qm.get_io2_value()['boolean_value']:
        while True:
            if self.keyPress == 'ESC':
                print('\033[94m' + 'ESC pressed. Stopping measurement.' + '\033[0m')  # print blue
                self.updateValue("Transit_Exp_switch", False)
                self.update_parameters()
                # Other actions can be added here
                break
            if reps < N:
                reps += 1

            ## For online plot of the moving average of the OD measurement: ##

            antihelmholtz_res = antihelmholtz_handle.fetch_all()

            if antihelmholtz_res:
                # avg_OD_N_res = avg_OD_N_handle.fetch_all()
                # avg_OD_S_res = avg_OD_S_handle.fetch_all()
                OD_N_res = OD_N_handle.fetch_all()
                OD_S_res = OD_S_handle.fetch_all()
                FLR_res = -FLR_handle.fetch_all()
                # FLR_measurement.append(FLR_res.tolist())
                FLR_measurement = FLR_measurement[-(N - 1):] + [FLR_res.tolist()]
                if np.mean(OD_S_res[shutter_open_indx:(
                        shutter_open_indx + (shutter_open_indx // t_shutter_open))]) * Per_usec < counts_threshold:
                    pass_threshold = True
                    N_OD_measurement.append(OD_N_res.tolist())
                    S_OD_measurement.append(OD_S_res.tolist())
            else:
                # avg_OD_N_no_MOT_res = avg_OD_N_no_MOT_handle.fetch_all()
                # avg_OD_S_no_MOT_res = avg_OD_S_no_MOT_handle.fetch_all()
                OD_N_no_MOT_res = OD_N_no_MOT_handle.fetch_all()
                OD_S_no_MOT_res = OD_S_no_MOT_handle.fetch_all()
                FLR_no_MOT_res = -FLR_no_MOT_handle.fetch_all()
                # FLR_no_MOT_measurement.append(FLR_no_MOT_res.tolist())
                FLR_no_MOT_measurement = FLR_no_MOT_measurement[-(N - 1):] + [FLR_no_MOT_res.tolist()]
                if np.mean(OD_S_res[shutter_open_indx:(
                        shutter_open_indx + (shutter_open_indx // t_shutter_open))]) * Per_usec < counts_threshold:
                    pass_threshold = True
                    N_OD_no_MOT_measurement.append(OD_N_no_MOT_res.tolist())
                    S_OD_no_MOT_measurement.append(OD_S_no_MOT_res.tolist())

            # measurements come in a whole vector at at time, and we want to save the last 100 measurements
            if pass_threshold:
                OD_N_res_batch = OD_N_res_batch[-(N - 1):] + [OD_N_res.tolist()]
                OD_N_no_MOT_res_batch = OD_N_no_MOT_res_batch[-(N - 1):] + [OD_N_no_MOT_res.tolist()]
                OD_S_res_batch = OD_S_res_batch[-(N - 1):] + [OD_S_res.tolist()]
                OD_S_no_MOT_res_batch = OD_S_no_MOT_res_batch[-(N - 1):] + [OD_S_no_MOT_res.tolist()]

            # if reps % 1 == 0:
            ax1.clear()
            ax2.clear()
            ax3.clear()
            ax4.clear()
            ax5.clear()

            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            textstr_N = r'$Probe_N = %.3f$' % (np.mean(OD_N_res[shutter_open_indx:]) * Per_usec,) + '[MPhotons/sec]'
            textstr_N_avg = r'$\overline{Probe}_N = %.3f$' % (
                np.mean([x[shutter_open_indx:] for x in OD_N_res_batch]) * Per_usec,) + '[MPhotons/sec]'
            textstr_S = r'$Probe_S = %.3f$' % (np.mean(OD_S_res[shutter_open_indx:]) * Per_usec,) + '[MPhotons/sec]'
            textstr_S_avg = r'$\overline{Probe}_S = %.3f$' % (
                np.mean([x[shutter_open_indx:] for x in OD_S_res_batch]) * Per_usec,) + '[MPhotons/sec]'
            textstr_FLR = r'$\overline{FLR}_{MAX} = %.1f$' % (np.mean(FLR_measurement) * 1e5,) + r'$\times 10^{-5}$'

            ax1.plot(t, [x * Per_usec for x in OD_N_res_batch[-1]], label='MOT ON - Inst', color='b')
            # ax1.plot(t, OD_N_res * Per_usec, label='MOT ON - Inst', color='b')
            ax1.plot(t, OD_N_no_MOT_res * Per_usec, label='MOT OFF - Inst', color='r')
            ax1.set_title('North', fontweight="bold")
            ax1.set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
            ax1.text(0.05, 0.95, textstr_N, transform=ax1.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
            ax1.legend(loc='upper right')

            ax2.plot(t, [x * Per_usec for x in OD_S_res_batch[-1]], label='MOT ON - Inst', color='b')
            # ax2.plot(t, OD_S_res * Per_usec, label='MOT ON - Inst', color='b')
            ax2.plot(t, OD_S_no_MOT_res * Per_usec, label='MOT OFF - Inst', color='r')
            ax2.set_title('South', fontweight="bold")
            ax2.set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
            ax2.text(0.05, 0.95, textstr_S, transform=ax2.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
            ax2.legend(loc='upper right')

            ax3.plot(t[bin_size // 2:bin_size * (len(t) // bin_size):bin_size],
                     np.mean(np.reshape(np.mean(OD_N_res_batch, 0)[:bin_size * (len(OD_N_res) // bin_size)]
                                        * Per_usec, (-1, bin_size)), 1), label='MOT ON - Avg', marker='*',
                     color='k')
            # ax3.plot(t[bin_size // 2:bin_size * (len(t) // bin_size):bin_size],
            #          np.mean(np.reshape(np.mean(OD_N_no_MOT_res_batch, 0)[:bin_size * (len(OD_N_res) // bin_size)]
            #                             * Per_usec, (-1, bin_size)), 1), label='MOT OFF - Avg', marker='*',
            #          color='#2ca02c')
            # ax3.set_title('North Average')
            ax3.set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
            ax3.text(0.05, 0.95, textstr_N_avg, transform=ax3.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
            ax3.legend(loc='upper right')

            ax4.plot(t[bin_size // 2:bin_size * (len(t) // bin_size):bin_size],
                     np.mean(np.reshape(np.mean(OD_S_res_batch, 0)[:bin_size * (len(OD_S_res) // bin_size)]
                                        * Per_usec, (-1, bin_size)), 1), label='MOT ON - Avg', marker='*',
                     color='k')
            # ax4.plot(t[bin_size // 2:bin_size * (len(t) // bin_size):bin_size],
            #          np.mean(np.reshape(np.mean(OD_S_no_MOT_res_batch, 0)[:bin_size * (len(OD_S_res) // bin_size)]
            #                             * Per_usec, (-1, bin_size)), 1), label='MOT OFF - Avg', marker='*',
            #          color='#2ca02c')
            # ax4.set_title('South Average')
            ax4.set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
            ax4.text(0.05, 0.95, textstr_S_avg, transform=ax4.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
            ax4.legend(loc='upper right')

            ax5.plot(FLR_measurement, label='MOT ON')
            ax5.plot(FLR_no_MOT_measurement, label='MOT OFF')
            ax5.set_title('Flouresence')
            ax5.set(xlabel='Repetitions', ylabel='Flouresence')
            ax5.text(0.05, 0.95, textstr_FLR, transform=ax5.transAxes, fontsize=12, verticalalignment='top',
                     bbox=props)
            ax5.legend(loc='upper right')

            ax1.set_ylim(0, 8)
            ax2.set_ylim(0, 8)
            # ax3.set_xlim(self.OD_delay - self.OD_duration_pulse1, self.OD_delay + 2 * self.OD_duration_pulse1)

            pass_threshold = False
            plt.show()
            plt.pause(0.9)

            ###################################################################

        ## Adding comment to measurement [prompt whether stopped or finished regularly]
        aftComment = pymsgbox.prompt('Add comment to measurement: ', default='', timeout=int(30e3))
        if aftComment == 'Timeout': aftComment = None

        #### Handle file-names and directories #####

        filename_N_avg = timest + '_' + 'Avg_PhotonCounts_North_rep=' + str(reps)
        filename_N = timest + '_' + 'PhotonCounts_North_rep=' + str(reps)
        filename_S_avg = timest + '_' + 'Avg_PhotonCounts_South_rep=' + str(reps)
        filename_S = timest + '_' + 'PhotonCounts_South_rep=' + str(reps)
        filename_FLR = timest + '_' + 'Flouresence_rep=' + str(reps)
        filename_FLR_no_MOT = timest + '_' + 'Flouresence_no_MOT_rep=' + str(reps)
        new_dir = 'PhotonCounts_' + datest
        Main_dir = 'U:\\Lab_2021-2022\\Experiment_results\\Transits\\'
        dirname = Main_dir + new_dir
        dirname_N = dirname + '\\' + 'North'
        dirname_N_no_MOT = dirname + '\\' + 'North_no_MOT'
        dirname_S = dirname + '\\' + 'South'
        dirname_S_no_MOT = dirname + '\\' + 'South_no_MOT'
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        if not os.path.exists(dirname_N):
            os.makedirs(dirname_N)
        if not os.path.exists(dirname_S):
            os.makedirs(dirname_S)
        if not os.path.exists(dirname_N_no_MOT):
            os.makedirs(dirname_N_no_MOT)
        if not os.path.exists(dirname_S_no_MOT):
            os.makedirs(dirname_S_no_MOT)

        if len(OD_N_res_batch) > 0:
            np.savetxt(dirname_N + '\\' + filename_N_avg + '.csv', np.vstack(OD_N_res_batch),
                       delimiter=',',
                       fmt='%s')
            np.savetxt(dirname_N + '\\' + filename_N + '.csv', np.transpose(N_OD_measurement), delimiter=',',
                       fmt='%s')
        if len(OD_S_res_batch) > 0:
            np.savetxt(dirname_S + '\\' + filename_S_avg + '.csv', np.vstack(OD_S_res_batch),
                       delimiter=',', fmt='%s')
            np.savetxt(dirname_S + '\\' + filename_S + '.csv', np.transpose(S_OD_measurement), delimiter=',',
                       fmt='%s')

        if FLR_measurement.size > 0:
            np.savetxt(dirname + '\\' + filename_FLR + '.csv', np.vstack(FLR_measurement), delimiter=',', fmt='%s')
        if FLR_no_MOT_measurement.size > 0:
            np.savetxt(dirname + '\\' + filename_FLR_no_MOT + '.csv', np.vstack(FLR_no_MOT_measurement),
                       delimiter=',', fmt='%s')

        # ## Edit comments file
        cmntDir = dirname + '\\experiment_comments.txt'
        cmnt = timest + ' - '
        if preComment is not None: cmnt = cmnt + preComment + '; '
        if aftComment is not None: cmnt = cmnt + aftComment
        if preComment is None and aftComment is None: cmnt = cmnt + 'No comment. '
        try:
            with open(cmntDir, "a") as commentsFile:
                commentsFile.write(cmnt + '\n')
        except:
            print('Could not save comments, error writing to comments-file.')
            print(cmnt)

    def Get_Max_Probe_counts(self, repetitions):

        self.Max_Probe_counts_switch(True)
        self.update_parameters()

        Probe_N_handle = self.job.result_handles.get("North_Probe")
        Probe_S_handle = self.job.result_handles.get("South_Probe")

        Probe_N_handle.wait_for_values(1)
        Probe_S_handle.wait_for_values(1)

        Probe_N_res = Probe_N_handle.fetch_all()
        Probe_S_res = Probe_S_handle.fetch_all()

        Probe_counts_North = [Probe_N_res.tolist()]
        Probe_counts_South = [Probe_S_res.tolist()]
        print(Probe_counts_South)

        for n in range(repetitions - 1):
            while Probe_S_res.tolist() == Probe_counts_South[-1]:
                 Probe_N_res = Probe_N_handle.fetch_all()
                 Probe_S_res = Probe_S_handle.fetch_all()
            Probe_counts_North.append(Probe_N_res.tolist())
            Probe_counts_South.append(Probe_S_res.tolist())

        print(r'$Max Probe_N = %.1f$' % ((np.average(Probe_counts_North) * 1000) / self.M_time,) + '[MPhotons/sec]')
        print(r'$Max Probe_S = %.1f$' % ((np.average(Probe_counts_South) * 1000) / self.M_time,) + '[MPhotons/sec]')

        return [(np.average(Probe_counts_North) * 1000) / self.M_time, (np.average(Probe_counts_South) * 1000) / self.M_time]

    def Save_SNSPDs_Transit_Measurement_with_tt(self, N, histogram_bin_size, Transit_profile_bin_size, preComment, counts_threshold, intensity_threshold, max_probe_counts):
        """
        Function for analyzing and saving the time tags data measured from the SNSPDs using the OPX. In this specific
         program we are looking for transits of atoms next to the toroid and record them.
        :param N: Number of maximum experiments (free throws) saved and displayed.
        :param histogram_bin_size: The bin size for the general experiment histogram which means - dividing the length
                                   of the measuring time (m_time) to bins and counting the number of photon detections
                                   at each bin.
        :param Transit_profile_bin_size: The bin size for the transit histogram which means - dividing the length of the
                                         transit time to bins and counting the number of transits for which there has
                                         been a detection of photon at each bin.
        :param preComment: The comment added at the start of the experiment, usually consisting of unique experiment
                           parameters.
        :param counts_threshold: The minimum number of MCounts / sec for which the data will be displayed and collected.
        :param intensity_threshold: The number of photons at each time bin for which suspect we as a transit.
        :param max_probe_counts: The maximum counts probe counts at each direction measured when cavity isn't locked.
        :return:
        """
        # if preComment is True:
        if not preComment:
            preComment = pymsgbox.prompt('Add comment to measurement: ', default='', timeout=int(30e3))
        aftComment = None

        # print(r'$Max Probe_N = %.1f$' % ((max_probe_counts[0] * 1000) / self.M_time,) + '[MPhotons/sec]')
        # print(r'$Max Probe_S = %.1f$' % ((max_probe_counts[1] * 1000) / self.M_time,) + '[MPhotons/sec]')

        ### fetching data from server
        ### saving to file
        ###

        histogram_bin_number = self.M_time // histogram_bin_size
        time_bins = np.linspace(0, self.M_time, histogram_bin_number)
        # time_threshold = int(histogram_bin_size / intensity_threshold)  # The minimum time between two time tags to be counted for a transit. # TODO: might need a factor of 2???
        time_threshold = int(histogram_bin_size * 0.5)  # The minimum time between two time tags to be counted for a transit. # TODO: might need a factor of 2???

        ## Listen for keyboard
        listener = keyboard.Listener(on_press=self.on_key_press)
        listener.start()  # start to listen on a separate thread
        self.keyPress = None
        print('\033[94m' + 'Press ESC to stop measurement.' + '\033[0m')  # print blue
        reps = 1  # a counter, number of repeats actually made.

        Probe_N_handle = self.job.result_handles.get("North_Probe")
        Probe_S_handle = self.job.result_handles.get("South_Probe")
        tt_N_handle = self.job.result_handles.get("North_Probe_TT")
        tt_S_handle = self.job.result_handles.get("South_Probe_TT")
        FLR_handle = self.job.result_handles.get("FLR_measure")

        Probe_N_handle.wait_for_values(1)
        Probe_S_handle.wait_for_values(1)
        tt_N_handle.wait_for_values(1)
        tt_S_handle.wait_for_values(1)
        FLR_handle.wait_for_values(1)

        Probe_N_res = Probe_N_handle.fetch_all()
        Probe_S_res = Probe_S_handle.fetch_all()
        tt_N_res = tt_N_handle.fetch_all()
        tt_S_res = tt_S_handle.fetch_all()
        FLR_res = -FLR_handle.fetch_all()

        tt_N_measure = [i for i in tt_N_res if (i % self.M_window) != 0]
        tt_S_measure = [i for i in tt_S_res if (i % self.M_window) != 0]
        tt_N_measure.sort()
        tt_S_measure.sort()

        tt_N_measure_batch = []
        tt_N_binning_batch = []
        self.tt_S_measure_batch = []
        tt_S_binning_batch = []
        all_transits_batch = []
        all_transits_aligned_first_batch = []
        transit_histogram_batch = []
        FLR_measurement = []
        Exp_timestr_batch = []

        tt_N_binning = np.zeros(histogram_bin_number)
        tt_S_binning = np.zeros(histogram_bin_number)
        tt_N_transit_events = np.zeros(histogram_bin_number)
        tt_S_transit_events = np.zeros(histogram_bin_number)

        # Place holders for results
        # while (number of photons * 10^-6 [Mcounts] / Measuring time [nsec] * 10^-9 [sec/nsec])  > counts_threshold [Mcounts /sec]:
        while ((len(tt_S_measure) * 1000) / self.M_time) > counts_threshold:
            if self.keyPress == 'ESC':
                print('\033[94m' + 'ESC pressed. Stopping measurement.' + '\033[0m')  # print blue
                self.updateValue("Transit_Exp_switch", False)
                self.update_parameters()
                # Other actions can be added here
                break

            print('Above Threshold')
            Probe_N_handle.wait_for_values(1)
            Probe_S_handle.wait_for_values(1)
            tt_N_handle.wait_for_values(1)
            tt_S_handle.wait_for_values(1)
            FLR_handle.wait_for_values(1)

            Probe_N_res = Probe_N_handle.fetch_all()
            Probe_S_res = Probe_S_handle.fetch_all()
            tt_N_res = tt_N_handle.fetch_all()
            tt_S_res = tt_S_handle.fetch_all()
            FLR_res = -FLR_handle.fetch_all()

            tt_N_measure = [i for i in tt_N_res if (i % self.M_window) != 0]
            tt_S_measure = [i for i in tt_S_res if (i % self.M_window) != 0]
            tt_N_measure.sort()
            tt_S_measure.sort()

        self.tt_S_measure = tt_S_measure
        ## record time
        timest = time.strftime("%Y%m%d-%H%M%S")
        datest = time.strftime("%Y%m%d")

        FLR_measurement = FLR_measurement[-(N - 1):] + [FLR_res.tolist()]
        Exp_timestr_batch = Exp_timestr_batch[-(N - 1):] + [timest]

        for x in tt_N_measure:
            tt_N_binning[x // histogram_bin_size] += 1
        for x in tt_S_measure:
            tt_S_binning[x // histogram_bin_size] += 1

        if len(self.tt_S_measure_batch) == N:
            tt_N_transit_events[[i for i, x in enumerate(tt_N_binning_batch[0]) if x > intensity_threshold]] -= 1
            tt_S_transit_events[[i for i, x in enumerate(tt_S_binning_batch[0]) if x > intensity_threshold]] -= 1

        tt_N_measure_batch = tt_N_measure_batch[-(N - 1):] + [tt_N_measure]
        tt_N_binning_batch = tt_N_binning_batch[-(N - 1):] + [tt_N_binning]
        self.tt_S_measure_batch = self.tt_S_measure_batch[-(N - 1):] + [tt_S_measure]
        tt_S_binning_batch = tt_S_binning_batch[-(N - 1):] + [tt_S_binning]
        Counter = 1

        tt_N_transit_events[[i for i, x in enumerate(tt_N_binning) if x > intensity_threshold]] += 1
        tt_S_transit_events[[i for i, x in enumerate(tt_S_binning) if x > intensity_threshold]] += 1

        # Find transits and build histogram:
        current_transit = []
        all_transits = []
        all_transits_aligned_first = []
        t_transit = []
        t_transit_batch = []
        transit_histogram = []
        for t in tt_S_measure:
            if not current_transit:  # if the array is empty
                current_transit.append(t)
            elif (t - current_transit[-1]) < time_threshold:
                current_transit.append(t)
            elif len(current_transit) > intensity_threshold:
                all_transits.append(current_transit)
                all_transits_aligned_first.append([x - current_transit[0] for x in current_transit])
                current_transit = [t]
            else:
                current_transit = [t]

        if all_transits:
            if len(all_transits_aligned_first_batch) == N:
                for n in range(len(all_transits_aligned_first_batch[0])):
                    for m in range(len(all_transits_aligned_first_batch[0][n]) - 1):
                        transit_histogram_batch[all_transits_aligned_first_batch[0][n][m + 1]//Transit_profile_bin_size] -= 1

            all_transits_batch = all_transits_batch[-(N - 1):] + [all_transits]
            all_transits_aligned_first_batch = all_transits_aligned_first_batch[-(N - 1):] + [
                all_transits_aligned_first]
            transit_histogram = np.zeros((max([vec for elem in all_transits_aligned_first for vec in elem])
                                          + histogram_bin_size)//Transit_profile_bin_size)
            t_transit = np.linspace(0, len(transit_histogram) * Transit_profile_bin_size, len(transit_histogram))
            transit_histogram_batch = np.zeros((max(
                [vec for elem in [vec for elem in all_transits_aligned_first_batch for vec in elem] for vec in
                 elem]) + histogram_bin_size)//Transit_profile_bin_size)
            t_transit_batch = np.linspace(0, len(transit_histogram_batch) * Transit_profile_bin_size, len(transit_histogram_batch))
            for n in range(len(all_transits_aligned_first)):
                for m in range(len(all_transits_aligned_first[n]) - 1):
                    transit_histogram[all_transits_aligned_first[n][m + 1]//Transit_profile_bin_size] += 1
                    transit_histogram_batch[all_transits_aligned_first[n][m + 1]//Transit_profile_bin_size] += 1


        ## Prepare plots
        fig = plt.figure()
        ax1 = plt.subplot2grid((6, 2), (0, 0), colspan=1, rowspan=2)
        ax2 = plt.subplot2grid((6, 2), (0, 1), colspan=1, rowspan=2)
        ax3 = plt.subplot2grid((6, 2), (2, 0), colspan=1, rowspan=2)
        ax4 = plt.subplot2grid((6, 2), (2, 1), colspan=1, rowspan=2)
        ax5 = plt.subplot2grid((6, 2), (4, 0), colspan=1, rowspan=2)
        ax6 = plt.subplot2grid((6, 2), (4, 1), colspan=1, rowspan=2)

        while True:
            if self.keyPress == 'ESC':
                print('\033[94m' + 'ESC pressed. Stopping measurement.' + '\033[0m')  # print blue
                self.updateValue("Transit_Exp_switch", False)
                self.update_parameters()
                # Other actions can be added here
                break
            if reps < N:
                reps += 1

            ########################## PLOT!!! ########################################################################

            ax1.clear()
            ax2.clear()
            ax3.clear()
            ax4.clear()
            ax5.clear()
            ax6.clear()

            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            textstr_N = r'$Probe_N = %.3f$' % ((len(tt_N_measure) * 1000) / self.M_time,) + '[MPhotons/sec]\n' \
                        + '$\overline{Probe}_N = %.3f$' % ((np.mean([len(x) for x in tt_N_measure_batch]) * 1000)
                                                           / self.M_time,) + '[MPhotons/sec]'
            textstr_S = r'$Probe_S = %.3f$' % ((len(tt_S_measure) * 1000) / self.M_time,) + '[MPhotons/sec]\n' \
                        + '$\overline{Probe}_S = %.3f$' % ((np.mean([len(x) for x in self.tt_S_measure_batch]) * 1000)
                                                           / self.M_time,) + '[MPhotons/sec]'
            textstr_FLR = r'$\overline{FLR}_{MAX} = %.1f$' % (np.mean(FLR_measurement) * 1e5,) + r'$\times 10^{-5}$'
            textstr_No_transits = 'NO TRANSITS YET!!!'

            ax1.plot(time_bins, tt_N_binning, label='Counts histogram', color='b')
            ax1.set_title('North', fontweight="bold")
            ax1.set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
            ax1.text(0.05, 0.95, textstr_N, transform=ax1.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
            ax1.legend(loc='upper right')

            ax2.plot(time_bins, tt_S_binning, label='Counts histogram', color='b')
            ax2.set_title('South', fontweight="bold")
            ax2.set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
            ax2.text(0.05, 0.95, textstr_S, transform=ax2.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
            ax2.legend(loc='upper right')

            ax3.plot(time_bins, tt_N_transit_events, label='Transit events histogram', marker='*', color='k')
            ax3.set(xlabel='Time [msec]', ylabel='Transits events [#Number]')
            ax3.legend(loc='upper right')

            ax4.plot(time_bins, tt_S_transit_events, label='Transit events histogram', marker='*', color='k')
            ax4.set(xlabel='Time [msec]', ylabel='Transits events [#Number]')
            ax4.legend(loc='upper right')

            if len(transit_histogram) > 0:
                textstr_transit_counts = r'$N_{Transits} = %s $' % (len(all_transits_aligned_first),) + r'$[Counts]$'
                textstr_avg_transit_counts = r'$\overline{N}_{Transits} = %.1f $' % (np.average([len(vec) for vec in all_transits_aligned_first]),) + r'$[Counts]$'

                ax5.plot(t_transit, transit_histogram, color='b')
                ax5.set_title('Drop transits profile', fontweight="bold")
                ax5.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
                ax5.text(0.05, 0.95, textstr_transit_counts, transform=ax5.transAxes, fontsize=12,
                         verticalalignment='top', bbox=props)

                ax6.plot(t_transit_batch, transit_histogram_batch, color='b')
                ax6.set_title('Accumulated drop transits profile', fontweight="bold")
                ax6.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
                ax6.text(0.05, 0.95, textstr_avg_transit_counts, transform=ax6.transAxes, fontsize=12,
                         verticalalignment='top', bbox=props)
            else:
                ax5.plot(t_transit, transit_histogram, color='b')
                ax5.set_title('Drop transits profile', fontweight="bold")
                ax5.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
                ax5.text(0.25, 0.5, textstr_No_transits, transform=ax5.transAxes, fontsize=24,
                         verticalalignment='center', bbox=props)

                ax6.plot(t_transit_batch, transit_histogram_batch, color='b')
                ax6.set_title('Accumulated drop transits profile', fontweight="bold")
                ax6.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
                ax6.text(0.25, 0.5, textstr_No_transits, transform=ax6.transAxes, fontsize=24,
                         verticalalignment='center', bbox=props)

            ax1.set_ylim(0, 8)
            ax2.set_ylim(0, 8)

            # plt.tight_layout()
            plt.show()
            plt.pause(1.5)

            ###########################################################################################################

            # while True:
                # record time:
            timest = time.strftime("%Y%m%d-%H%M%S")
            datest = time.strftime("%Y%m%d")

            # get measures:
            Probe_N_res = Probe_N_handle.fetch_all()
            Probe_S_res = Probe_S_handle.fetch_all()
            tt_N_res = tt_N_handle.fetch_all()
            tt_S_res = tt_S_handle.fetch_all()
            FLR_res = -FLR_handle.fetch_all()

            tt_N_measure = [i for i in tt_N_res if (i % self.M_window) != 0]
            tt_S_measure = [i for i in tt_S_res if (i % self.M_window) != 0]
            tt_N_measure.sort()
            tt_S_measure.sort()

                # if self.tt_S_measure != self.tt_S_measure_batch[-1]:
                #     break

            if ((len(tt_S_measure) * 1000) / self.M_time) < counts_threshold:

                FLR_measurement = FLR_measurement[-(N - 1):] + [FLR_res.tolist()]
                Exp_timestr_batch = Exp_timestr_batch[-(N - 1):] + [timest]
                Counter += 1
                print(timest, Counter)

                tt_N_binning = np.zeros(histogram_bin_number)
                tt_S_binning = np.zeros(histogram_bin_number)

                for x in tt_N_measure:
                    tt_N_binning[x // histogram_bin_size] += 1
                for x in tt_S_measure:
                    tt_S_binning[x // histogram_bin_size] += 1

                if len(self.tt_S_measure_batch) == N:
                    tt_N_transit_events[[i for i, x in enumerate(tt_N_binning_batch[0]) if x > intensity_threshold]] -= 1
                    tt_S_transit_events[[i for i, x in enumerate(tt_S_binning_batch[0]) if x > intensity_threshold]] -= 1

                tt_N_measure_batch = tt_N_measure_batch[-(N - 1):] + [tt_N_measure]
                tt_N_binning_batch = tt_N_binning_batch[-(N - 1):] + [tt_N_binning]
                self.tt_S_measure_batch = self.tt_S_measure_batch[-(N - 1):] + [tt_S_measure]
                tt_S_binning_batch = tt_S_binning_batch[-(N - 1):] + [tt_S_binning]

                tt_N_transit_events[[i for i, x in enumerate(tt_N_binning) if x > intensity_threshold]] += 1
                tt_S_transit_events[[i for i, x in enumerate(tt_S_binning) if x > intensity_threshold]] += 1

                # Find transits and build histogram:
                current_transit = []
                all_transits = []
                all_transits_aligned_first = []
                for t in tt_S_measure:
                    if not current_transit:  # if the array is empty
                        current_transit.append(t)
                    elif (t - current_transit[-1]) < time_threshold:
                        current_transit.append(t)
                    elif len(current_transit) > intensity_threshold:
                        all_transits.append(current_transit)
                        all_transits_aligned_first.append([x - current_transit[0] for x in current_transit])
                        current_transit = [t]
                    else:
                        current_transit = [t]

                if all_transits:
                    if len(all_transits_aligned_first_batch) == N:
                        for n in range(len(all_transits_aligned_first_batch[0])):
                            for m in range(len(all_transits_aligned_first_batch[0][n]) - 1):
                                transit_histogram_batch[
                                    all_transits_aligned_first_batch[0][n][m + 1] // Transit_profile_bin_size] -= 1

                    all_transits_batch = all_transits_batch[-(N - 1):] + [all_transits]
                    all_transits_aligned_first_batch = all_transits_aligned_first_batch[-(N - 1):] + [
                        all_transits_aligned_first]

                    transit_histogram = np.zeros((max([vec for elem in all_transits_aligned_first for vec in elem])
                                                  + histogram_bin_size) // Transit_profile_bin_size)
                    t_transit = np.linspace(0, len(transit_histogram) * Transit_profile_bin_size,
                                            len(transit_histogram))
                    transit_histogram_batch = np.zeros((max(
                        [vec for elem in [vec for elem in all_transits_aligned_first_batch for vec in elem] for vec in
                         elem]) + histogram_bin_size) // Transit_profile_bin_size)
                    # np.pad(transit_histogram_batch, (0, (max(
                    #     [vec for elem in [vec for elem in all_transits_aligned_first_batch for vec in elem] for vec in
                    #      elem]) + histogram_bin_size) // Transit_profile_bin_size - len(transit_histogram_batch)),
                    #        'constant')
                    t_transit_batch = np.linspace(0, len(transit_histogram_batch) * Transit_profile_bin_size,
                                                  len(transit_histogram_batch))
                    for n in range(len(all_transits_aligned_first)):
                        for m in range(len(all_transits_aligned_first[n]) - 1):
                            transit_histogram[all_transits_aligned_first[n][m + 1] // Transit_profile_bin_size] += 1
                            transit_histogram_batch[all_transits_aligned_first[n][m + 1] // Transit_profile_bin_size] += 1

        ############################################## END WHILE LOOP #################################################

        ## Adding comment to measurement [prompt whether stopped or finished regularly]
        aftComment = pymsgbox.prompt('Add comment to measurement: ', default='', timeout=int(30e3))
        if aftComment == 'Timeout': aftComment = None

        #### Handle file-names and directories #####
        ## Saving: np.savez(filedir, data = x) #note: @filedir is string of full directory; data is the queyword used to read @x from the file:
        ## Loading: file = np.load(f, allow_pickle = True)
        ##          x = file['data']

        #### ------ Save results ------
        #  -------   Create dir
        root_dirname = f'U:\\Lab_2021-2022\\Experiment_results\\Transits\\{datest}\\'
        dirname = root_dirname + f'{timest}_Photon_TimeTags\\'  # Specific experiment dir
        dirname_N = dirname + 'North\\'
        dirname_S = dirname + 'South\\'
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        if not os.path.exists(dirname_N):
            os.makedirs(dirname_N)
        if not os.path.exists(dirname_S):
            os.makedirs(dirname_S)

        # ----  msmnt files names  -----
        # Counter_str = (Counter)
        filename_N_tt = f'North_timetags.npz'
        filename_S_tt = f'South_timetags.npz'
        filename_S_transits = f'South_Transits.npz'
        filename_FLR = f'Flouresence.npz'
        filename_timestamp = f'Drops_time_stamps.npz'

        if len(FLR_measurement) > 0:
            np.savez(dirname + filename_FLR, FLR_measurement)
        if len(Exp_timestr_batch) > 0:
            np.savez(dirname + filename_timestamp, Exp_timestr_batch)
        if len(tt_N_measure_batch) > 0:
            np.savez(dirname_N + filename_N_tt, tt_N_measure_batch)
        if len(self.tt_S_measure_batch) > 0:
            np.savez(dirname_S + filename_S_tt, self.tt_S_measure_batch)
        if len(all_transits_batch) > 0:
            np.savez(dirname_S + filename_S_transits, all_transits_batch)


        ### Edit comments file ####
        cmntDir = root_dirname + '\\daily_experiment_comments.txt'
        cmnt = timest + ' - '
        if preComment is not None: cmnt = cmnt + preComment + '; '
        if aftComment is not None: cmnt = cmnt + aftComment
        if preComment is None and aftComment is None: cmnt = cmnt + 'No comment. '
        try:
            with open(cmntDir, "a") as commentsFile:
                commentsFile.write(cmnt + '\n')
        except:
            print('Could not save comments, error writing to comments-file.')
            print(cmnt)

        comments = {'comments': cmnt}
        try:
            with open(f'{dirname}experiment_comments.txt', 'w') as file:
                json.dump(comments, file, indent=4)
        except Exception:
            pass

        self.saveConfigTable(path = dirname)
        ## ------------------ end of saving section -------

    def Save_SNSPDs_Spectrum_Measurement_with_tt(self, N, histogram_bin_size, Transit_profile_bin_size, preComment, counts_threshold, intensity_threshold, max_probe_counts):
        """
        Function for analyzing and saving the time tags data measured from the SNSPDs using the OPX. In this specific
         program we are looking for transits of atoms next to the toroid and record them.
        :param N: Number of maximum experiments (free throws) saved and displayed.
        :param histogram_bin_size: The bin size for the general experiment histogram which means - dividing the length
                                   of the measuring time (m_time) to bins and counting the number of photon detections
                                   at each bin.
        :param Transit_profile_bin_size: The bin size for the transit histogram which means - dividing the length of the
                                         transit time to bins and counting the number of transits for which there has
                                         been a detection of photon at each bin.
        :param preComment: The comment added at the start of the experiment, usually consisting of unique experiment
                           parameters.
        :param counts_threshold: The minimum number of MCounts / sec for which the data will be displayed and collected.
        :param intensity_threshold: The number of photons at each time bin for which suspect we as a transit.
        :param max_probe_counts: The maximum counts probe counts at each direction measured when cavity isn't locked.
        :return:
        """
        # if preComment is True:
        if not preComment:
            preComment = pymsgbox.prompt('Add comment to measurement: ', default='', timeout=int(30e3))
        aftComment = None

        # print(r'$Max Probe_N = %.1f$' % ((max_probe_counts[0] * 1000) / self.M_time,) + '[MPhotons/sec]')
        # print(r'$Max Probe_S = %.1f$' % ((max_probe_counts[1] * 1000) / self.M_time,) + '[MPhotons/sec]')

        ### fetching data from server
        ### saving to file
        ###

        histogram_bin_number = self.M_time // histogram_bin_size
        time_bins = np.linspace(0, self.M_time, histogram_bin_number)
        # time_threshold = int(histogram_bin_size / intensity_threshold)  # The minimum time between two time tags to be counted for a transit. # TODO: might need a factor of 2???
        time_threshold = int(histogram_bin_size * 0.5)  # The minimum time between two time tags to be counted for a transit. # TODO: might need a factor of 2???

        ## Listen for keyboard
        listener = keyboard.Listener(on_press=self.on_key_press)
        listener.start()  # start to listen on a separate thread
        self.keyPress = None
        print('\033[94m' + 'Press ESC to stop measurement.' + '\033[0m')  # print blue
        reps = 1  # a counter, number of repeats actually made.

        Probe_N_handle = self.job.result_handles.get("North_Probe")
        Probe_S_handle = self.job.result_handles.get("South_Probe")
        tt_N_handle = self.job.result_handles.get("North_Probe_TT")
        tt_S_handle = self.job.result_handles.get("South_Probe_TT")
        FLR_handle = self.job.result_handles.get("FLR_measure")

        Probe_N_handle.wait_for_values(1)
        Probe_S_handle.wait_for_values(1)
        tt_N_handle.wait_for_values(1)
        tt_S_handle.wait_for_values(1)
        FLR_handle.wait_for_values(1)

        Probe_N_res = Probe_N_handle.fetch_all()
        Probe_S_res = Probe_S_handle.fetch_all()
        tt_N_res = tt_N_handle.fetch_all()
        tt_S_res = tt_S_handle.fetch_all()
        FLR_res = -FLR_handle.fetch_all()

        tt_N_measure = [i for i in tt_N_res if (i % self.M_window) != 0]
        tt_S_measure = [i for i in tt_S_res if (i % self.M_window) != 0]
        tt_N_measure.sort()
        tt_S_measure.sort()

        tt_N_measure_batch = []
        tt_N_binning_batch = []
        self.tt_S_measure_batch = []
        tt_S_binning_batch = []
        all_transits_batch = []
        all_transits_aligned_first_batch = []
        transit_histogram_batch = []
        FLR_measurement = []
        Exp_timestr_batch = []

        tt_N_binning = np.zeros(histogram_bin_number)
        tt_S_binning = np.zeros(histogram_bin_number)
        tt_N_transit_events = np.zeros(histogram_bin_number)
        tt_S_transit_events = np.zeros(histogram_bin_number)

        # Place holders for results
        # while (number of photons * 10^-6 [Mcounts] / Measuring time [nsec] * 10^-9 [sec/nsec])  > counts_threshold [Mcounts /sec]:
        while ((len(tt_S_measure) * 1000) / self.M_time) > counts_threshold:

            print('Above Threshold')
            Probe_N_handle.wait_for_values(1)
            Probe_S_handle.wait_for_values(1)
            tt_N_handle.wait_for_values(1)
            tt_S_handle.wait_for_values(1)
            FLR_handle.wait_for_values(1)

            Probe_N_res = Probe_N_handle.fetch_all()
            Probe_S_res = Probe_S_handle.fetch_all()
            tt_N_res = tt_N_handle.fetch_all()
            tt_S_res = tt_S_handle.fetch_all()
            FLR_res = -FLR_handle.fetch_all()

            tt_N_measure = [i for i in tt_N_res if (i % self.M_window) != 0]
            tt_S_measure = [i for i in tt_S_res if (i % self.M_window) != 0]
            tt_N_measure.sort()
            tt_S_measure.sort()

        self.tt_S_measure = tt_S_measure
        ## record time
        timest = time.strftime("%Y%m%d-%H%M%S")
        datest = time.strftime("%Y%m%d")

        FLR_measurement = FLR_measurement[-(N - 1):] + [FLR_res.tolist()]
        Exp_timestr_batch = Exp_timestr_batch[-(N - 1):] + [timest]

        for x in tt_N_measure:
            tt_N_binning[x // histogram_bin_size] += 1
        for x in tt_S_measure:
            tt_S_binning[x // histogram_bin_size] += 1

        if len(self.tt_S_measure_batch) == N:
            tt_N_transit_events[[i for i, x in enumerate(tt_N_binning_batch[0]) if x > intensity_threshold]] -= 1
            tt_S_transit_events[[i for i, x in enumerate(tt_S_binning_batch[0]) if x > intensity_threshold]] -= 1

        tt_N_measure_batch = tt_N_measure_batch[-(N - 1):] + [tt_N_measure]
        tt_N_binning_batch = tt_N_binning_batch[-(N - 1):] + [tt_N_binning]
        self.tt_S_measure_batch = self.tt_S_measure_batch[-(N - 1):] + [tt_S_measure]
        tt_S_binning_batch = tt_S_binning_batch[-(N - 1):] + [tt_S_binning]
        Counter = 1

        tt_N_transit_events[[i for i, x in enumerate(tt_N_binning) if x > intensity_threshold]] += 1
        tt_S_transit_events[[i for i, x in enumerate(tt_S_binning) if x > intensity_threshold]] += 1

        # Find transits and build histogram:
        current_transit = []
        all_transits = []
        all_transits_aligned_first = []
        t_transit = []
        t_transit_batch = []
        transit_histogram = []
        for t in tt_S_measure:
            if not current_transit:  # if the array is empty
                current_transit.append(t)
            elif (t - current_transit[-1]) < time_threshold:
                current_transit.append(t)
            elif len(current_transit) > intensity_threshold:
                all_transits.append(current_transit)
                all_transits_aligned_first.append([x - current_transit[0] for x in current_transit])
                current_transit = [t]
            else:
                current_transit = [t]

        if all_transits:
            if len(all_transits_aligned_first_batch) == N:
                for n in range(len(all_transits_aligned_first_batch[0])):
                    for m in range(len(all_transits_aligned_first_batch[0][n]) - 1):
                        transit_histogram_batch[all_transits_aligned_first_batch[0][n][m + 1]//Transit_profile_bin_size] -= 1

            all_transits_batch = all_transits_batch[-(N - 1):] + [all_transits]
            all_transits_aligned_first_batch = all_transits_aligned_first_batch[-(N - 1):] + [
                all_transits_aligned_first]
            transit_histogram = np.zeros((max([vec for elem in all_transits_aligned_first for vec in elem])
                                          + histogram_bin_size)//Transit_profile_bin_size)
            t_transit = np.linspace(0, len(transit_histogram) * Transit_profile_bin_size, len(transit_histogram))
            transit_histogram_batch = np.zeros((max(
                [vec for elem in [vec for elem in all_transits_aligned_first_batch for vec in elem] for vec in
                 elem]) + histogram_bin_size)//Transit_profile_bin_size)
            t_transit_batch = np.linspace(0, len(transit_histogram_batch) * Transit_profile_bin_size, len(transit_histogram_batch))
            for n in range(len(all_transits_aligned_first)):
                for m in range(len(all_transits_aligned_first[n]) - 1):
                    transit_histogram[all_transits_aligned_first[n][m + 1]//Transit_profile_bin_size] += 1
                    transit_histogram_batch[all_transits_aligned_first[n][m + 1]//Transit_profile_bin_size] += 1


        ## Prepare plots
        fig = plt.figure()
        ax1 = plt.subplot2grid((6, 2), (0, 0), colspan=1, rowspan=2)
        ax2 = plt.subplot2grid((6, 2), (0, 1), colspan=1, rowspan=2)
        ax3 = plt.subplot2grid((6, 2), (2, 0), colspan=1, rowspan=2)
        ax4 = plt.subplot2grid((6, 2), (2, 1), colspan=1, rowspan=2)
        ax5 = plt.subplot2grid((6, 2), (4, 0), colspan=1, rowspan=2)
        ax6 = plt.subplot2grid((6, 2), (4, 1), colspan=1, rowspan=2)

        while True:
            if self.keyPress == 'ESC':
                print('\033[94m' + 'ESC pressed. Stopping measurement.' + '\033[0m')  # print blue
                self.updateValue("Transit_Exp_switch", False)
                self.update_parameters()
                # Other actions can be added here
                break
            if reps < N:
                reps += 1

            ########################## PLOT!!! ########################################################################

            ax1.clear()
            ax2.clear()
            ax3.clear()
            ax4.clear()
            ax5.clear()
            ax6.clear()

            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            textstr_N = r'$Probe_N = %.3f$' % ((len(tt_N_measure) * 1000) / self.M_time,) + '[MPhotons/sec]\n' \
                        + '$\overline{Probe}_N = %.3f$' % ((np.mean([len(x) for x in tt_N_measure_batch]) * 1000)
                                                           / self.M_time,) + '[MPhotons/sec]'
            textstr_S = r'$Probe_S = %.3f$' % ((len(tt_S_measure) * 1000) / self.M_time,) + '[MPhotons/sec]\n' \
                        + '$\overline{Probe}_S = %.3f$' % ((np.mean([len(x) for x in self.tt_S_measure_batch]) * 1000)
                                                           / self.M_time,) + '[MPhotons/sec]'
            textstr_FLR = r'$\overline{FLR}_{MAX} = %.1f$' % (np.mean(FLR_measurement) * 1e5,) + r'$\times 10^{-5}$'
            textstr_No_transits = 'NO TRANSITS YET!!!'

            ax1.plot(time_bins, tt_N_binning, label='Counts histogram', color='b')
            ax1.set_title('North', fontweight="bold")
            ax1.set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
            ax1.text(0.05, 0.95, textstr_N, transform=ax1.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
            ax1.legend(loc='upper right')

            ax2.plot(time_bins, tt_S_binning, label='Counts histogram', color='b')
            ax2.set_title('South', fontweight="bold")
            ax2.set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
            ax2.text(0.05, 0.95, textstr_S, transform=ax2.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
            ax2.legend(loc='upper right')

            ax3.plot(time_bins, tt_N_transit_events, label='Transit events histogram', marker='*', color='k')
            ax3.set(xlabel='Time [msec]', ylabel='Transits events [#Number]')
            ax3.legend(loc='upper right')

            ax4.plot(time_bins, tt_S_transit_events, label='Transit events histogram', marker='*', color='k')
            ax4.set(xlabel='Time [msec]', ylabel='Transits events [#Number]')
            ax4.legend(loc='upper right')

            if len(transit_histogram) > 0:
                textstr_transit_counts = r'$N_{Transits} = %s $' % (len(all_transits_aligned_first),) + r'$[Counts]$'
                textstr_avg_transit_counts = r'$\overline{N}_{Transits} = %.1f $' % (np.average([len(vec) for vec in all_transits_aligned_first]),) + r'$[Counts]$'

                ax5.plot(t_transit, transit_histogram, color='b')
                ax5.set_title('Drop transits profile', fontweight="bold")
                ax5.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
                ax5.text(0.05, 0.95, textstr_transit_counts, transform=ax5.transAxes, fontsize=12,
                         verticalalignment='top', bbox=props)

                ax6.plot(t_transit_batch, transit_histogram_batch, color='b')
                ax6.set_title('Accumulated drop transits profile', fontweight="bold")
                ax6.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
                ax6.text(0.05, 0.95, textstr_avg_transit_counts, transform=ax6.transAxes, fontsize=12,
                         verticalalignment='top', bbox=props)
            else:
                ax5.plot(t_transit, transit_histogram, color='b')
                ax5.set_title('Drop transits profile', fontweight="bold")
                ax5.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
                ax5.text(0.25, 0.5, textstr_No_transits, transform=ax5.transAxes, fontsize=24,
                         verticalalignment='center', bbox=props)

                ax6.plot(t_transit_batch, transit_histogram_batch, color='b')
                ax6.set_title('Accumulated drop transits profile', fontweight="bold")
                ax6.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
                ax6.text(0.25, 0.5, textstr_No_transits, transform=ax6.transAxes, fontsize=24,
                         verticalalignment='center', bbox=props)

            ax1.set_ylim(0, 8)
            ax2.set_ylim(0, 8)

            # plt.tight_layout()
            plt.show()
            plt.pause(1)

            ###########################################################################################################

            while tt_S_measure == self.tt_S_measure_batch[-1]:
                # record time:
                timest = time.strftime("%Y%m%d-%H%M%S")
                datest = time.strftime("%Y%m%d")

                # get measures:
                Probe_N_res = Probe_N_handle.fetch_all()
                Probe_S_res = Probe_S_handle.fetch_all()
                tt_N_res = tt_N_handle.fetch_all()
                tt_S_res = tt_S_handle.fetch_all()
                FLR_res = -FLR_handle.fetch_all()

                tt_N_measure = [i for i in tt_N_res if (i % self.M_window) != 0]
                tt_S_measure = [i for i in tt_S_res if (i % self.M_window) != 0]
                tt_N_measure.sort()
                tt_S_measure.sort()

            if ((len(tt_S_measure) * 1000) / self.M_time) < counts_threshold:

                FLR_measurement = FLR_measurement[-(N - 1):] + [FLR_res.tolist()]
                Exp_timestr_batch = Exp_timestr_batch[-(N - 1):] + [timest]
                Counter += 1
                print(timest, Counter)

                tt_N_binning = np.zeros(histogram_bin_number)
                tt_S_binning = np.zeros(histogram_bin_number)

                for x in tt_N_measure:
                    tt_N_binning[x // histogram_bin_size] += 1
                for x in tt_S_measure:
                    tt_S_binning[x // histogram_bin_size] += 1

                if len(self.tt_S_measure_batch) == N:
                    tt_N_transit_events[[i for i, x in enumerate(tt_N_binning_batch[0]) if x > intensity_threshold]] -= 1
                    tt_S_transit_events[[i for i, x in enumerate(tt_S_binning_batch[0]) if x > intensity_threshold]] -= 1

                tt_N_measure_batch = tt_N_measure_batch[-(N - 1):] + [tt_N_measure]
                tt_N_binning_batch = tt_N_binning_batch[-(N - 1):] + [tt_N_binning]
                self.tt_S_measure_batch = self.tt_S_measure_batch[-(N - 1):] + [tt_S_measure]
                tt_S_binning_batch = tt_S_binning_batch[-(N - 1):] + [tt_S_binning]

                tt_N_transit_events[[i for i, x in enumerate(tt_N_binning) if x > intensity_threshold]] += 1
                tt_S_transit_events[[i for i, x in enumerate(tt_S_binning) if x > intensity_threshold]] += 1

                # Find transits and build histogram:
                current_transit = []
                all_transits = []
                all_transits_aligned_first = []
                for t in tt_S_measure:
                    if not current_transit:  # if the array is empty
                        current_transit.append(t)
                    elif (t - current_transit[-1]) < time_threshold:
                        current_transit.append(t)
                    elif len(current_transit) > intensity_threshold:
                        all_transits.append(current_transit)
                        all_transits_aligned_first.append([x - current_transit[0] for x in current_transit])
                        current_transit = [t]
                    else:
                        current_transit = [t]

                if all_transits:
                    if len(all_transits_aligned_first_batch) == N:
                        for n in range(len(all_transits_aligned_first_batch[0])):
                            for m in range(len(all_transits_aligned_first_batch[0][n]) - 1):
                                transit_histogram_batch[
                                    all_transits_aligned_first_batch[0][n][m + 1] // Transit_profile_bin_size] -= 1

                    all_transits_batch = all_transits_batch[-(N - 1):] + [all_transits]
                    all_transits_aligned_first_batch = all_transits_aligned_first_batch[-(N - 1):] + [
                        all_transits_aligned_first]

                    transit_histogram = np.zeros((max([vec for elem in all_transits_aligned_first for vec in elem])
                                                  + histogram_bin_size) // Transit_profile_bin_size)
                    t_transit = np.linspace(0, len(transit_histogram) * Transit_profile_bin_size,
                                            len(transit_histogram))
                    # transit_histogram_batch = np.zeros((max(
                    #     [vec for elem in [vec for elem in all_transits_aligned_first_batch for vec in elem] for vec in
                    #      elem]) + histogram_bin_size) // Transit_profile_bin_size)
                    np.pad(transit_histogram_batch, (0, (max(
                        [vec for elem in [vec for elem in all_transits_aligned_first_batch for vec in elem] for vec in
                         elem]) + histogram_bin_size) // Transit_profile_bin_size - len(transit_histogram_batch)),
                           'constant')
                    t_transit_batch = np.linspace(0, len(transit_histogram_batch) * Transit_profile_bin_size,
                                                  len(transit_histogram_batch))
                    for n in range(len(all_transits_aligned_first)):
                        for m in range(len(all_transits_aligned_first[n]) - 1):
                            transit_histogram[all_transits_aligned_first[n][m + 1] // Transit_profile_bin_size] += 1
                            transit_histogram_batch[all_transits_aligned_first[n][m + 1] // Transit_profile_bin_size] += 1

        ############################################## END WHILE LOOP #################################################

        ## Adding comment to measurement [prompt whether stopped or finished regularly]
        aftComment = pymsgbox.prompt('Add comment to measurement: ', default='', timeout=int(30e3))
        if aftComment == 'Timeout': aftComment = None

        #### Handle file-names and directories #####
        ## Saving: np.savez(filedir, data = x) #note: @filedir is string of full directory; data is the queyword used to read @x from the file:
        ## Loading: file = np.load(f, allow_pickle = True)
        ##          x = file['data']

        #### ------ Save results ------
        #  -------   Create dir
        root_dirname = f'U:\\Lab_2021-2022\\Experiment_results\\Transits\\{datest}\\'
        dirname = root_dirname + f'{timest}_Photon_TimeTags\\'  # Specific experiment dir
        dirname_N = dirname + 'North\\'
        dirname_S = dirname + 'South\\'
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        if not os.path.exists(dirname_N):
            os.makedirs(dirname_N)
        if not os.path.exists(dirname_S):
            os.makedirs(dirname_S)

        # ----  msmnt files names  -----
        # Counter_str = (Counter)
        filename_N_tt = f'North_timetags.npz'
        filename_S_tt = f'South_timetags.npz'
        filename_S_transits = f'South_Transits.npz'
        filename_FLR = f'Flouresence.npz'
        filename_timestamp = f'Drops_time_stamps.npz'

        if len(FLR_measurement) > 0:
            np.savez(dirname + filename_FLR, FLR_measurement)
        if len(Exp_timestr_batch) > 0:
            np.savez(dirname + filename_timestamp, Exp_timestr_batch)
        if len(tt_N_measure_batch) > 0:
            np.savez(dirname_N + filename_N_tt, tt_N_measure_batch)
        if len(self.tt_S_measure_batch) > 0:
            np.savez(dirname_S + filename_S_tt, self.tt_S_measure_batch)
        if len(all_transits_batch) > 0:
            np.savez(dirname_S + filename_S_transits, all_transits_batch)


        ### Edit comments file ####
        cmntDir = root_dirname + '\\daily_experiment_comments.txt'
        cmnt = timest + ' - '
        if preComment is not None: cmnt = cmnt + preComment + '; '
        if aftComment is not None: cmnt = cmnt + aftComment
        if preComment is None and aftComment is None: cmnt = cmnt + 'No comment. '
        try:
            with open(cmntDir, "a") as commentsFile:
                commentsFile.write(cmnt + '\n')
        except:
            print('Could not save comments, error writing to comments-file.')
            print(cmnt)

        comments = {'comments': cmnt}
        try:
            with open(f'{dirname}experiment_comments.txt', 'w') as file:
                json.dump(comments, file, indent=4)
        except Exception:
            pass

        self.saveConfigTable(path = dirname)
        ## ------------------ end of saving section -------

    def Start_Transit_Exp(self, N=100, bin_size=100, preComment=None, threshold=0.1):
        self.Transit_Exp_switch(True)
        self.update_parameters()
        self.Save_SNSPDs_Measurement(N, bin_size, preComment, threshold)

    def Start_Transit_Exp_with_tt(self, N=100, Histogram_bin_size=1000, Transit_profile_bin_size=100, preComment=None, threshold=1, intensity_threshold=5):
        # Max_probe_counts = self.Get_Max_Probe_counts(3)  # return the average maximum probe counts of 3 cycles.
        Max_probe_counts = None
        self.Transit_Exp_switch(True)
        self.update_parameters()
        self.Save_SNSPDs_Transit_Measurement_with_tt(N, Histogram_bin_size, Transit_profile_bin_size, preComment, threshold, intensity_threshold, Max_probe_counts)

    def Start_Spectrum_Exp_with_tt(self, N=100, Histogram_bin_size=1000, Transit_profile_bin_size=100, preComment=None, threshold=1, intensity_threshold=5):
        Max_probe_counts = self.Get_Max_Probe_counts(3)  # return the average maximum probe counts of 3 cycles.
        self.Spectrum_Exp_switch(True)
        self.update_parameters()
        self.Save_SNSPDs_Spectrum_Measurement_with_tt(N, Histogram_bin_size, Transit_profile_bin_size, preComment, threshold, intensity_threshold, Max_probe_counts)

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

    # @Values_Factor hold factoring for values - each key has an array:
    # So, as Prep_duration is written in ms, we have to factor it into an int value, but in 4ns, for the OPX.
    # These are only factored in function "updateValue" in OPX control
    Values_Factor = {
        'Transit_Exp_switch': [bool, None, Transit_Exp_switch],
        'MOT_duration': [int, 1e6 / Config.MOT_pulse_len],
        'AntiHelmholtz_delay': [int, 1e6 / 4],  # [msec]
        'Post_MOT_delay': [int, 1e6 / 4],  # [msec]
        'N_Snaps':  [int, None, update_N_Snaps],
        'Buffer_Cycles': [int, 1],
        ## PGC parameters ##
        "PGC_duration": [int, 1e6 / 4],  # [msec]
        "PGC_prep_duration": [int, None, update_PGC_prep_time],
        "PGC_initial_amp_0": [int, None, update_PGC_AOM_0_initial_amplitude],
        "PGC_initial_amp_minus": [int, None, update_PGC_AOM_minus_initial_amplitude],
        "PGC_initial_amp_plus": [int, None, update_PGC_AOM_plus_initial_amplitude],
        "PGC_final_amp": [float, 1],  # update_PGC_AOMs_final_amplitude],
        "PGC_final_amp_0": [int, None, update_PGC_AOM_0_final_amplitude],
        "PGC_final_amp_minus": [int, None, update_PGC_AOM_minus_initial_amplitude],
        "PGC_final_amp_plus": [int, None, update_PGC_AOM_plus_initial_amplitude],
        ## Fountain parameters ##
        "Pre_PGC_Fountain_duration": [int, 1e6 / 4],
        "Fountain_duration": [int, 1e6 / 4],
        "Fountain_initial_amp_0": [float, None, update_fountain_AOM_0_initial_amplitude],
        "Fountain_initial_amp_minus": [float, None, update_fountain_AOM_minus_initial_amplitude],
        "Fountain_initial_amp_plus": [float, None, update_fountain_AOM_plus_initial_amplitude],
        "Fountain_prep_duration": [int, None, update_fountain_prep_time],
        "Fountain_final_amp_0": [float, None, update_fountain_AOM_0_final_amplitude],
        "Fountain_final_amp_minus": [float, None, update_fountain_AOM_minus_final_amplitude],
        "Fountain_final_amp_plus": [float, None, update_fountain_AOM_plus_final_amplitude],
        "Fountain_final_Delta_freq":[float, None, update_fountain_Delta_freq],
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
        'M_time': [int, 1e6 / 4],  # [msec]'
    }


if __name__ == "__main__":
    experiment = OPX(Config.config)
    # experiment.updateValue('PrePulse_duration', 5)
    # experiment.update_parameters()
    # experiment.Repeat_Measurement(10)
    # experiment.Repeat_Measurement(10)

