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
from scipy.io import savemat
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


def moving_average(a, n=3):
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
          "AOM_2-2/3'", "AOM_2-3'_for_interference", "FLR_detection", "Measurement")  # , "Dig_detectors") # , "PULSER_N", "PULSER_S")
    # update_frequency("AOM_Spectrum", 80000000)
    ## MOT build-up ##
    n = declare(int)
    m = declare(int)
    play("Detection" * amp(FLR * 0), "FLR_detection",
         duration=4)  # we don't know why this works, but Yoav from QM made us write this line to solve an alignment problem we had in the next 2 for loops
    with for_(n, 1, n <= mot_repetitions, n + 1):
        play("MOT_with_EOM" * amp(Config.AOM_0_Attenuation), "MOT_AOM_0")
        play("MOT" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-")
        play("MOT" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+")
        # with if_(n == mot_repetitions):
        #     play("Const_open", "AOM_Spectrum", duration=50000)
        # play("Const_open", "AOM_Spectrum")
        # play("Const_open", "PULSER_N")
        # play("OD_FS" * amp(0.1), "AOM_2-2/3'")
        play("AntiHelmholtz_MOT", "AntiHelmholtz_Coils")
    with for_(m, 1, m <= (mot_repetitions - 1), m + 1):
        measure("Detection", "FLR_detection", None, integration.full("Detection_opt", FLR, "out1"))
        # play("OD_FS" * amp(0.1), "AOM_2-3'_for_interference")

    align("Cooling_Sequence", "MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+", "AntiHelmholtz_Coils", "Zeeman_Coils",
          "AOM_2-2/3'", "FLR_detection", "Measurement")  # , "Dig_detectors") #, "PULSER_N", "PULSER_S")

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
          "Measurement", "AOM_Early", "AOM_Late", "PULSER_ANCILLA", "PULSER_N", "PULSER_S")  # , "Dig_detectors")

    ## Zeeman Coils turn-on sequence ##
    wait(coils_timing, "Zeeman_Coils")
    play("ZeemanOFF", "Zeeman_Coils", duration=(freefall_duration - coils_timing))
    # play("ZeemanSplit", "Zeeman_Coils", duration=(freefall_duration - coils_timing))


def Measure(measuring_duration):
    """
    The Measure function is an all purpose tool for any measurement is needed.

    Parameters
    ----------
    :param measuring_duration: Given in [msec] and plays a trigger to the measuring component.
    :return:
    """

    play('Measure', "Measurement", duration=measuring_duration)


def Transit_Exp(m_off_time, m_time, m_window, shutter_open_time,
                ON_counts_st1, ON_counts_st2, ON_counts_st3,
                ON_counts_st4, ON_counts_st5, ON_counts_st6,
                ON_counts_st7, ON_counts_st8,
                tt_st_1, tt_st_2, tt_st_3, tt_st_4, tt_st_5, tt_st_6, tt_st_7, tt_st_8, rep_st):
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
    m = declare(int)

    # assign_variables_to_element("Dig_detectors", tt_vec1[0], counts1, m_window)
    align("AOM_Late", "PULSER_N", "PULSER_S", "Dig_detectors", "PULSER_ANCILLA", "AOM_Spectrum")

    # Analog signal is disabled (multiplied by 0), only digital is sent
    # Sent to shutters
    play("Const_open_triggered" * amp(0), "PULSER_N", duration=shutter_open_time)
    play("Const_open" * amp(0), "PULSER_S", duration=shutter_open_time)

    align("AOM_Late", "PULSER_N", "PULSER_S", "Dig_detectors", "PULSER_ANCILLA", "AOM_Spectrum")

    play("Const_open_triggered", "PULSER_N", duration=(m_time + m_off_time))
    play("Const_open_triggered", "PULSER_S", duration=(m_time + m_off_time))
    play("Const_open", "AOM_Late", duration=(m_time + m_off_time))
    # play("OD_FS", "AOM_2-2/3'", duration=(m_time + m_off_time))

    # main experiment - send on resonance during whole experiment
    play("Spectrum_pulse", "AOM_Spectrum", duration=(m_time + m_off_time))

    # with for_(t, 0, t < (m_time + m_off_time) * 4, t + int(Config.frequency_sweep_duration)):  # assaf comment debbuging
    #     play("Spectrum_pulse", "AOM_Spectrum")
        # update_frequency("AOM_Spectrum", )

    # wait(298, "Dig_detectors")
    wait(330, "Dig_detectors")
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


def Spectrum_Exp(m_off_time, m_time, m_window, shutter_open_time,
                 frequency_sweep_rep, same_frequency_rep, num_of_different_frequencies, frequency_diff, frequency_start,
                 ON_counts_st1, ON_counts_st2, ON_counts_st3,
                 ON_counts_st4, ON_counts_st5, ON_counts_st6,
                 ON_counts_st7, ON_counts_st8,
                 tt_st_1, tt_st_2, tt_st_3, tt_st_4, tt_st_5, tt_st_6, tt_st_7, tt_st_8, rep_st):
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
    sweep_num = declare(int)
    freq_rep_num = declare(int)
    freq = declare(int)
    m = declare(int)

    update_frequency("PULSER_ANCILLA", Config.IF_AOM_Spectrum)

    # assign_variables_to_element("Dig_detectors", tt_vec1[0], counts1, m_window)
    align("AOM_Late", "PULSER_N", "PULSER_S", "Dig_detectors", "PULSER_ANCILLA", "AOM_Spectrum")

    play("Const_open_triggered" * amp(0), "PULSER_N", duration=shutter_open_time)
    play("Const_open" * amp(0), "PULSER_S", duration=shutter_open_time)

    align("AOM_Late", "PULSER_N", "PULSER_S", "Dig_detectors", "PULSER_ANCILLA", "AOM_Spectrum")

    play("Const_open_triggered", "PULSER_N", duration=(m_time + m_off_time))
    play("Const_open_triggered", "PULSER_S", duration=(m_time + m_off_time))
    play("Const_open", "AOM_Late", duration=(m_time + m_off_time))
    # play("OD_FS", "AOM_2-2/3'", duration=(m_time + m_off_time))
    # play("Spectrum_pulse", "AOM_Spectrum", duration=(m_time + m_off_time))
    # play("Spectrum_pulse", "PULSER_ANCILLA", duration=(m_time + m_off_time))

    with for_(sweep_num, 0, sweep_num < frequency_sweep_rep, sweep_num + 1):
        with for_(freq, frequency_start, freq < (frequency_start + num_of_different_frequencies * frequency_diff), freq + frequency_diff):
            update_frequency("AOM_Spectrum", freq)
            with for_(freq_rep_num, 0, freq_rep_num < same_frequency_rep, freq_rep_num + 1):
                play("Spectrum_pulse", "PULSER_ANCILLA")
                # wait(int(Config.frequency_sweep_duration / 4), "AOM_Spectrum")
                align("PULSER_ANCILLA", "AOM_Spectrum")
                # wait(int(Config.frequency_sweep_duration / 4), "PULSER_ANCILLA")
                play("Spectrum_pulse", "AOM_Spectrum")


    wait(339, "Dig_detectors")
    # wait(330, "Dig_detectors")
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
        fountain_pulse_duration_0 = declare(int, value=int(
            obj.fountain_pulse_duration_0))  # The relative duration to reach the desired amplitude
        fountain_pulse_duration_minus = declare(int, value=int(
            obj.fountain_pulse_duration_minus))  # The relative duration to reach the desired amplitude
        fountain_pulse_duration_plus = declare(int, value=int(
            obj.fountain_pulse_duration_plus))  # The relative duration to reach the desired amplitude
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

        # push beam params
        PushBeam_duration = declare(int, value=int(20000/4))
        OD_freq = declare(int, value=int(Config.IF_AOM_OD))
        PushBeam_Amp = declare(fixed, value=1)

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

        # OD measurement variables:
        ## In fiber:
        M_off_time = declare(int, value=int(obj.M_off_time / 4))
        shutter_open_time = declare(int, value=int(obj.Shutter_open_time * 1e6 / 4))
        OD_freq = declare(int, value=int(Config.IF_AOM_Spectrum))

        ## Spectrum experiment:
        frequency_sweep_rep = declare(int, value=obj.frequency_sweep_rep)
        same_frequency_rep = declare(int, value=obj.same_frequency_rep)
        num_of_different_frequencies = declare(int, value=int(obj.num_of_different_frequencies))
        frequency_start = declare(int, value=obj.frequency_start)
        frequency_diff = declare(int, value=obj.frequency_diff)

        # Stream processing:
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
                assign(x, (
                            30678780 - 3106 + 656000 * 2 + 4) // 4)  # TODO -  added 38688900 to fix new delay due to wait(1000) in saving sprint data with vector size 10000, should be fixed as well
            with else_():
                assign(x, - 3106)
            FreeFall(FreeFall_duration - x, coils_timing)

            ##########################
            ## Measurement Sequence ##
            ##########################

            with if_(Trigger_Phase == 3):  # when trigger on PrePulse
                ## Trigger QuadRF Sequence #####################
                play("C_Seq", "Cooling_Sequence", duration=2500)
                ################################################

            align("Cooling_Sequence", "AOM_Early", "AOM_Late")
            with if_(QRAM_Exp_ON):
                wait(PrePulse_duration - shutter_open_time, "Cooling_Sequence")
            with else_():
                wait(PrePulse_duration, "Cooling_Sequence")
            # push beam
            play("OD_FS" * amp(PushBeam_Amp), "AOM_2-2/3'", duration=PushBeam_duration)

            align(*all_elements, "AOM_2-2/3'", "AOM_Early", "AOM_Late", "PULSER_ANCILLA", "PULSER_N", "PULSER_S",
                  "Dig_detectors", "AOM_Spectrum")

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
                align("Dig_detectors", "AOM_Early", "AOM_Late", "PULSER_ANCILLA", "PULSER_N", "PULSER_S", "AOM_Spectrum")
                # Transit_Exp(M_off_time, Pulse_1_duration, obj.M_window, shutter_open_time,
                #             ON_counts_st1, ON_counts_st2, ON_counts_st3,
                #             ON_counts_st4, ON_counts_st5, ON_counts_st6,
                #             ON_counts_st7, ON_counts_st8,
                #             tt_st_1, tt_st_2, tt_st_3, tt_st_4, tt_st_5, tt_st_6, tt_st_7, tt_st_8, rep_st)
                Spectrum_Exp(M_off_time, Pulse_1_duration, obj.M_window, shutter_open_time,
                             frequency_sweep_rep, same_frequency_rep, num_of_different_frequencies, frequency_diff,
                             frequency_start,
                             ON_counts_st1, ON_counts_st2, ON_counts_st3,
                             ON_counts_st4, ON_counts_st5, ON_counts_st6,
                             ON_counts_st7, ON_counts_st8,
                             tt_st_1, tt_st_2, tt_st_3, tt_st_4, tt_st_5, tt_st_6, tt_st_7, tt_st_8, rep_st)
                align("Dig_detectors", "AOM_Early", "AOM_Late", "PULSER_ANCILLA", "PULSER_N", "PULSER_S", "AOM_Spectrum")

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
                with if_(i == 5):
                    # update_frequency("AOM_Spectrum", Config.IF_AOM_Spectrum + int(30e6))
                    assign(OD_freq, IO2)
                    update_frequency("AOM_Spectrum", OD_freq)
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

                with if_(i == 47):
                    assign(PushBeam_duration, IO2)
                with if_(i == 48):
                    assign(PushBeam_Amp, IO2)
                with if_(i == 49):
                    assign(OD_freq, IO2)

                pause()
                assign(i, IO1)

        with stream_processing():
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
        qrfContr = QuadRFMOTController(initialValues=self.Exp_Values, updateChannels=(1, 2, 4),
                                       topticaLockWhenUpdating=False,
                                       debugging=False, continuous=False)
        self.QuadRFControllers.append(qrfContr)  # updates values on QuadRF (uploads table)
        self.QuadRFControllers.append(QuadRFMOTController(MOGdevice=qrfContr.dev,
                                                          initialValues={'Operation_Mode': 'Continuous',
                                                                         'CH3_freq': '90MHz', 'CH3_amp': '31dbm'},
                                                          updateChannels=[3], debugging=False,
                                                          continuous=False))  # updates values on QuadRF (uploads table)
        # self.QuadRFControllers.append(QuadRFFrequencyScannerController(MOGdevice = qrfContr.dev, channel=2, debugging=False))  # updates values on QuadRF (uploads table)

        self.Update_QuadRF_channels = set(
            {})  # Only update these channels on QuadRF when UpdateParameters method is called [note: this is a python set]
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
        self.Balancing_check_window = 2  # [msec]
        self.rep_MZ_scan = int(0.5 * (self.prepulse_duration - self.Shutter_open_time - self.Balancing_check_window)
                               * 1e6 / len(Config.QRAM_MZ_balance_pulse_Late))
        self.rep_MZ_fast_scan = int(
            0.3 * (self.prepulse_duration - self.Shutter_open_time - self.Balancing_check_window)
            * 1e6 / len(Config.QRAM_MZ_balance_pulse_Late))
        self.rep_MZ_slow_scan = int(
            0.4 * (self.prepulse_duration - self.Shutter_open_time - self.Balancing_check_window)
            * 1e6 / len(Config.QRAM_MZ_balance_pulse_Late))
        self.points_for_sum = 5
        self.points_for_sum_fast = 3
        self.points_for_sum_slow = 4
        if self.rep_MZ_scan == 0:
            self.phase_rep_MZ = 1
        else:
            self.phase_rep_MZ = int(self.rep_MZ_scan / self.points_for_sum)
        if self.rep_MZ_slow_scan == 0:
            self.phase_rep_MZ_slow_scan = 1
        else:
            self.phase_rep_MZ_slow_scan = int(self.rep_MZ_slow_scan / self.points_for_sum_slow) - 1
        if self.rep_MZ_fast_scan == 0:
            self.phase_rep_MZ_fast_scan = 1
        else:
            self.phase_rep_MZ_fast_scan = int(self.rep_MZ_fast_scan / self.points_for_sum_fast)
        self.total_phase_rep_MZ = int(2 * self.phase_rep_MZ)
        self.total_phase_rep_MZ_scan = 2 * self.phase_rep_MZ_fast_scan + self.phase_rep_MZ_slow_scan
        self.rep_MZ_check = int(self.Balancing_check_window * 1e6 / len(Config.QRAM_MZ_balance_pulse_North))

        # MW spectroscopy parameters:
        self.MW_start_frequency = int(100e6)  # [Hz]
        self.Pulse_Length_MW = 400  # [usec]
        self.Pulse_Length_OD = 20  # [usec]

        # Main Experiment:
        self.TOP2_pulse_len = int(Config.Probe_pulse_len / 4)  # [nsec]
        self.Calibration_time = 10  # [msec]

        # Spectrum Experiment:
        self.frequency_sweep_rep = 5
        # self.same_frequency_rep = 100
        # self.spectrum_bandwidth = int(30e6)
        # self.spectrum_bandwidth = int(39e6)
        self.spectrum_bandwidth = int(48e6) # experiment spe13/11/23
        self.frequency_start = int(80e6 - self.spectrum_bandwidth/2)
        # self.num_of_different_frequencies = self.Exp_Values['Pulse_1_duration'] * 1e6 / \
        #                                   (2 * Config.frequency_sweep_duration) / \
        #                                   (self.frequency_sweep_rep * self.same_frequency_rep)
        # self.frequency_diff = int(self.spectrum_bandwidth / self.num_of_different_frequencies)
        self.frequency_diff = int(1.5e6) # difference between frequencies
        self.num_of_different_frequencies = int(self.spectrum_bandwidth // self.frequency_diff) + 1
        self.same_frequency_rep = int(self.Exp_Values['Pulse_1_duration'] * 1e6 /
                                      (2 * Config.frequency_sweep_duration) /
                                      (self.num_of_different_frequencies * self.frequency_sweep_rep)) # total time[10ms]/((1us sequence - 2 pulses time)*(num of sweeps)*(num of different freqs) )
        
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

    def Change_OD_freq(self, freq):
        self.update_io_parameter(5, int(Config.IF_AOM_Spectrum + freq))
        # self.update_io_parameter(5, int(74.65e6))


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

        # push beam params
        def update_PushBeam_duration(self, duration):
            self.update_io_parameter(47, int(duration * 1e3 / 4))  # In [us]

        def update_PushBeam_amp(self, Amp):
            self.update_io_parameter(48, float(Amp))  #

        def update_PushBeam_frequency(self, freq):
            self.update_io_parameter(49, int(Config.IF_AOM_OD + freq))  # In [us]

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

    def get_avg_num_of_photons_in_det_pulse(self, det_pulse_len, sprint_sequence_delay, num_of_det_pulses,
                                            num_of_sprint_sequences):
        self.avg_num_of_photons_in_det_pulse = np.zeros(num_of_det_pulses)
        for i in range(num_of_det_pulses):
            detection_puls_ind = \
                list(np.arange((sprint_sequence_delay + i * det_pulse_len),
                               (sprint_sequence_delay + (i + 1) * det_pulse_len)))
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
            Counts_handle.append(self.job.result_handles.get("Det" + str(i) + "_Counts"))
            tt_handle.append(self.job.result_handles.get("Det" + str(i) + "_Probe_TT"))
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

        # add the tt and counts to python vars
        for i in range(len(Num_Of_dets)):
            self.counts_res.append(Counts_handle[i].fetch_all())
            self.tt_res.append(tt_handle[i].fetch_all())

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
                                      ((elm + Config.detector_delays[
                                          i]) <= self.M_window))]  # Due to an unresolved bug in the OPD there are "ghost" readings of timetags equal to the maximum time of measuring window.
            self.tt_measure[i].sort()
        # unify detectors and windows within detectors and create vector of tt's for each direction (bright port, dark port, north, south and from FS) and sort them
        self.tt_BP_measure = sorted(sum(self.tt_measure[:2], []))  # unify detectors 1-3 and windows within detectors
        self.tt_DP_measure = sorted(sum(self.tt_measure[2:4], []))  # unify detectors 1-3 and windows within detectors
        self.tt_N_measure = sorted(sum(self.tt_measure[7:], []))  # unify detectors 1-3 and windows within detectors
        self.tt_S_measure = sorted(sum(self.tt_measure[4:5], []))  # unify detectors 6-8 and windows within detectors
        self.tt_FS_measure = sorted(sum(self.tt_measure[5:7], []))  # unify detectors 6-8 and windows within detectors
        self.tt_N_directional_measure = sorted(self.tt_N_measure + self.tt_BP_measure + self.tt_DP_measure)
        # self.tt_S_directional_measure = sorted(self.tt_N_measure + self.tt_BP_measure + self.tt_DP_measure)
        self.tt_S_directional_measure = sorted(self.tt_S_measure + self.tt_FS_measure)
        # self.tt_N_directional_measure = sorted(self.tt_S_measure + self.tt_FS_measure)
        # No need to fix as we're now not switching frequencies. This fix is relevant only for Spectrum experiment
        self.fix_gaps_spectrum_exp_tts()

    def fix_gaps_spectrum_exp_tts(self):
        '''
        fixes delays of 320 ns every 100 us in spectrum experiment
        :return:

        '''
        self.real_M_window = ((2 * Config.frequency_sweep_duration) * (self.num_of_different_frequencies
                                                                       * self.same_frequency_rep
                                                                       * self.frequency_sweep_rep))
        self.total_real_time_of_freq_sweep = int(self.real_M_window // self.frequency_sweep_rep) + \
                                             self.num_of_different_frequencies * 320 + 124
        self.total_real_time_of_same_freq_rep = (Config.frequency_sweep_duration * 2 * self.same_frequency_rep) + 320
        self.tt_S_no_gaps = [x-(x//self.total_real_time_of_freq_sweep)*124 for x in self.tt_S_directional_measure]
        self.tt_S_no_gaps = [x-(x//self.total_real_time_of_same_freq_rep)*320 for x in self.tt_S_no_gaps]

        # self.tt_S_no_gaps = [x-(x//((Config.frequency_sweep_duration * 2 * self.same_frequency_rep) + 320))*320 for x in self.tt_S_directional_measure]
        # self.tt_S_no_gaps = [x-(x//(int(self.real_M_window // self.frequency_sweep_rep) + 124))*124 for x in self.tt_S_no_gaps]
        # self.tt_S_no_gaps = [int(x) for x in self.tt_S_no_gaps if x < self.real_M_window]

        # self.tt_S_no_gaps = []
        # for rep in range(self.frequency_sweep_rep):
        #     self.start_tt = rep * self.total_real_time_of_freq_sweep
        #     self.end_tt = (rep + 1) * self.total_real_time_of_freq_sweep
        #     self.tt_S_no_gaps.append([x-(x//((Config.frequency_sweep_duration * 2 * self.same_frequency_rep) + 320))*320
        #                               for x in self.tt_S_directional_measure[]])
                   
    def latched_detectors(self):
        latched_detectors = []
        for indx, det_tt_vec in enumerate(self.tt_measure[4:7]):  # for different detectors
            if not det_tt_vec:
                latched_detectors.append(indx)
        return latched_detectors

    def get_pulses_bins(self, det_pulse_len, sprint_pulse_len, num_of_det_pulses, num_of_sprint_pulses,
                        sprint_sequence_delay, num_of_sprint_sequences, num_init_zeros, num_fin_zeros,
                        num_between_zeros):
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
        pulses_bins = [np.maximum(sprint_sequence_delay, 0)]
        for i in range(num_of_sprint_sequences):
            pulses_bins += (pulses_bins[-1] + np.arange(det_pulse_len, (num_of_det_pulses + 1) * det_pulse_len,
                                                        det_pulse_len)).tolist()
            pulses_bins += (pulses_bins[-1] + np.arange(sprint_pulse_len, (num_of_sprint_pulses + 1) * sprint_pulse_len,
                                                        sprint_pulse_len)).tolist()
            pulses_bins[-1] += num_fin_zeros + num_init_zeros - num_between_zeros
        return pulses_bins

    def divide_to_reflection_trans(self, det_pulse_len, sprint_pulse_len, sprint_sequence_delay_S,
                                   sprint_sequence_delay_N,
                                   num_of_det_pulses,
                                   num_of_sprint_pulses, num_of_sprint_sequences):
        '''
        dividing south and north tt's vectros into reflection and transmission vectors, by:
        1. converting them into counts vector in same time spacing as the initial pulse.
        2. taking the relevant pulses from the sequence to each vector (reflection/transmission)
        :return:
        '''
        # create histogram with self.M_window bins
        self.pulses_bins_S = self.get_pulses_bins(det_pulse_len, sprint_pulse_len, num_of_det_pulses,
                                                  num_of_sprint_pulses
                                                  , sprint_sequence_delay_S, num_of_sprint_sequences,
                                                  Config.num_init_zeros_S
                                                  , Config.num_fin_zeros_S, Config.num_between_zeros)
        self.pulses_bins_N = self.get_pulses_bins(det_pulse_len, sprint_pulse_len, num_of_det_pulses,
                                                  num_of_sprint_pulses
                                                  , sprint_sequence_delay_N, num_of_sprint_sequences,
                                                  Config.num_init_zeros_N
                                                  , Config.num_fin_zeros_N, Config.num_between_zeros)
        self.tt_histogram_N, _ = np.histogram(self.tt_N_measure, self.pulses_bins_N)
        self.tt_histogram_S, _ = np.histogram(self.tt_S_measure, self.pulses_bins_S)

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

        self.num_of_SPRINT_reflections_per_seq_S = np.zeros(
            [self.number_of_QRAM_sequences, self.number_of_SPRINT_pulses_per_seq])
        self.num_of_SPRINT_reflections_per_seq_N = np.zeros(
            [self.number_of_QRAM_sequences, self.number_of_SPRINT_pulses_per_seq])
        self.num_of_SPRINT_transmissions_per_seq_S = np.zeros(
            [self.number_of_QRAM_sequences, self.number_of_SPRINT_pulses_per_seq])
        self.num_of_SPRINT_transmissions_per_seq_N = np.zeros(
            [self.number_of_QRAM_sequences, self.number_of_SPRINT_pulses_per_seq])

        # tt_small_perturb = []
        for element in self.tt_N_measure + self.tt_BP_measure + self.tt_DP_measure:
            if (element > int(0.6e6)) and (element < int(9.6e6)):
                tt_inseq = element % self.QRAM_sequence_len
                if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                    self.num_of_det_reflections_per_seq_S[(element - 1) // self.QRAM_sequence_len] += self.filter_S[
                        tt_inseq]
                    self.num_of_det_transmissions_per_seq_N[(element - 1) // self.QRAM_sequence_len] += self.filter_N[
                        tt_inseq]
            # else:  # The part of the SPRINT pulses in the sequence
            #     SPRINT_pulse_num = (tt_inseq - self.end_of_det_pulse_in_seq) // (sprint_pulse_len + Config.num_between_zeros)
            #     if SPRINT_pulse_num < self.number_of_SPRINT_pulses_per_seq:
            #         self.num_of_SPRINT_reflections_per_seq_S[(element-1) // self.QRAM_sequence_len][SPRINT_pulse_num] += 1
            #         self.num_of_SPRINT_transmissions_per_seq_N[(element-1) // self.QRAM_sequence_len][SPRINT_pulse_num] += 1
            #         tt_small_perturb+=[element]

        for element in self.tt_S_measure + self.tt_FS_measure:
            if (element > int(0.6e6)) and (element < int(9.6e6)):
                tt_inseq = element % self.QRAM_sequence_len
                if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                    self.num_of_det_reflections_per_seq_N[(element - 1) // self.QRAM_sequence_len] += self.filter_N[
                        tt_inseq]
                    self.num_of_det_transmissions_per_seq_S[(element - 1) // self.QRAM_sequence_len] += self.filter_S[
                        tt_inseq]
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
        self.num_of_BP_counts_per_n_sequences = np.zeros(self.number_of_QRAM_sequences // num_of_seq_per_count)
        self.num_of_DP_counts_per_n_sequences = np.zeros(self.number_of_QRAM_sequences // num_of_seq_per_count)
        self.num_of_S_counts_per_n_sequences = np.zeros(self.number_of_QRAM_sequences // num_of_seq_per_count)

        for element in self.tt_BP_measure:
            tt_inseq = element % self.QRAM_sequence_len
            if tt_inseq > self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                self.num_of_BP_counts_per_n_sequences[
                    (element - 1) // (self.QRAM_sequence_len * num_of_seq_per_count)] += \
                    np.ceil(Config.QRAM_Exp_Square_samples_Late[tt_inseq])

        for element in self.tt_DP_measure:
            tt_inseq = element % self.QRAM_sequence_len
            if tt_inseq > self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                self.num_of_DP_counts_per_n_sequences[
                    (element - 1) // (self.QRAM_sequence_len * num_of_seq_per_count)] += \
                    np.ceil(Config.QRAM_Exp_Square_samples_Late[tt_inseq])

        for element in self.tt_FS_measure + self.tt_S_measure:
            tt_inseq = element % self.QRAM_sequence_len
            if tt_inseq > self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                self.num_of_S_counts_per_n_sequences[
                    (element - 1) // (self.QRAM_sequence_len * num_of_seq_per_count)] += \
                    np.ceil(Config.QRAM_Exp_Square_samples_Late[tt_inseq])

    def plot_folded_tt_histogram(self):

        plt.figure()
        plt.plot(sum(np.reshape(self.tt_histogram_N, [int(len(self.tt_histogram_N) / 9), 9])), label='tt_hist_N')
        plt.plot(sum(np.reshape(self.tt_histogram_S, [int(len(self.tt_histogram_S) / 9), 9])), label='tt_hist_S')
        plt.legend()

    def find_transit_events_transitexp(self, N, transit_time_threshold,transit_counts_threshold):
        # Find transits and build histogram:
        current_transit = []
        self.all_transits = []

        for t in self.tt_S_directional_measure:
            if not current_transit:  # if the array is empty
                current_transit.append(t)
            elif (t - current_transit[-1]) < transit_time_threshold:
                current_transit.append(t)
            elif len(current_transit) > transit_counts_threshold:
                self.all_transits.append(current_transit)
                current_transit = [t]
            else:
                current_transit = [t]

            self.tt_S_transit_events[[i for i in [vec for elem in self.all_transits for vec in elem]]] += 1
            self.tt_S_transit_events_accumulated = self.tt_S_transit_events_accumulated + self.tt_S_transit_events

            if self.all_transits:
                self.all_transits_batch = self.all_transits_batch[-(N - 1):] + [self.all_transits]

    def find_transit_events_spectrum(self, N, transit_time_threshold, transit_counts_threshold):
        current_transit = []
        self.all_transits = []
        all_transits_aligned_first = []
        t_transit = []
        transit_histogram = []
        for index, value in enumerate(self.tt_S_binning_resonance):
            if not current_transit and value:  # if the array is empty
                current_transit.append(index)
            elif value:
                if ((index - current_transit[0]) * Config.frequency_sweep_duration * 2) < transit_time_threshold:
                    current_transit.append(index)
                elif sum([self.tt_S_binning_resonance[i] for i in current_transit]) > transit_counts_threshold:
                    self.all_transits.append(current_transit)
                    all_transits_aligned_first.append([x - current_transit[0] for x in current_transit])
                    # tt_S_transit_events[tuple(current_transit)] += 1
                    self.tt_per_frequency = np.zeros(self.spectrum_bin_number)
                    #building cavity atom spectrum
                    for i in current_transit[:-1]:
                        # building cavit-atom spectrum with the counts of the detuned time bins between the start and
                        # finish of the transit:
                        self.tt_per_frequency[i % (self.num_of_different_frequencies * self.same_frequency_rep)
                                              // self.same_frequency_rep] += self.tt_S_binning_detuned[i]
                    self.Cavity_atom_spectrum += self.tt_per_frequency
                    self.Transits_per_freuency += (self.tt_per_frequency != 0).astype(int)
                    current_transit = [index]
                else:
                    # Finding if there any index that was saved to current transit and is close enough to the new index
                    t = [i for i, elem in enumerate(current_transit) if ((index - elem) * Config.frequency_sweep_duration*2) < transit_time_threshold]
                    if t:
                        current_transit = current_transit[t[0]:] + [index]
                    else:
                        current_transit = [index]

        #building cavity spectrum
        for index, value in enumerate(self.tt_S_binning_detuned):
            if index not in [vec for elem in self.all_transits for vec in elem]:
                self.Cavity_spectrum[index % (self.num_of_different_frequencies * self.same_frequency_rep)
                                          // self.same_frequency_rep] += value

        self.Cavity_spectrum_normalized = (self.Cavity_spectrum / self.power_per_freq_weight)
        # self.Cavity_spectrum_normalized = (self.Cavity_spectrum / (self.power_per_freq_weight * self.Transits_per_freuency))
        self.Cavity_spectrum_normalized = self.Cavity_spectrum_normalized / max(self.Cavity_spectrum_normalized)
        self.Cavity_atom_spectrum_normalized = (self.Cavity_atom_spectrum / self.power_per_freq_weight)
        # self.Cavity_atom_spectrum_normalized = (self.Cavity_atom_spectrum / (self.power_per_freq_weight * self.Transits_per_freuency))
        self.Cavity_atom_spectrum_normalized = self.Cavity_atom_spectrum_normalized / max(self.Cavity_atom_spectrum_normalized)

        self.tt_S_transit_events[[i for i in [vec for elem in self.all_transits for vec in elem]]] += 1
        self.tt_S_transit_events_accumulated = self.tt_S_transit_events_accumulated + self.tt_S_transit_events

        if self.all_transits:
            self.all_transits_batch = self.all_transits_batch[-(N - 1):] + [self.all_transits]
            # transit_histogram = np.zeros((max([vec for elem in all_transits_aligned_first for vec in elem])
            #                               + self.histogram_bin_size) // self.Transit_profile_bin_size)
            # t_transit = np.linspace(0, len(transit_histogram) * self.Transit_profile_bin_size, len(transit_histogram))

    def find_transits_events_spectrum_exp(self, tt_resonance_binning, N, cond=[2,1,2], minimum_number_of_seq_detected=2):
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
        tt_resonance_binning = np.array(tt_resonance_binning)

        for i in range(len(tt_resonance_binning[:]) - len(cond) + 1):
            cond_check = (tt_resonance_binning[i:(i + len(cond))] >= cond).astype(int)
            if sum(cond_check) >= minimum_number_of_seq_detected:
                current_transit = np.unique(
                    current_transit + [*range(i + np.where(cond_check != 0)[0][0], (i + len(cond)))]).tolist()
            elif len(current_transit) > 1:
                current_transit = current_transit[:np.where(tt_resonance_binning[current_transit] >= min(cond))[0][-1] + 1]
                if self.all_transits_seq_indx:
                    if bool(set(current_transit) & set(self.all_transits_seq_indx[-1])):
                        current_transit = self.all_transits_seq_indx[-1] + current_transit[1:]
                        self.all_transits_seq_indx = self.all_transits_seq_indx[:-1]
                self.all_transits_seq_indx.append(current_transit)
                current_transit = []
        if len(current_transit) > 1:
            current_transit = current_transit[
                              :np.where(tt_resonance_binning[current_transit] >= min(cond))[0][-1] + 1]
            if self.all_transits_seq_indx:
                if bool(set(current_transit) & set(self.all_transits_seq_indx[-1])):
                    current_transit = self.all_transits_seq_indx[-1] + current_transit[1:]
                    self.all_transits_seq_indx = self.all_transits_seq_indx[:-1]
            self.all_transits_seq_indx.append(current_transit)

        # building cavity atom spectrum
        for current_transit in self.all_transits_seq_indx:
            self.tt_per_frequency = np.zeros(self.spectrum_bin_number)
            for i in current_transit[:-1]:
                # building cavit-atom spectrum with the counts of the detuned time bins between the start and
                # finish of the transit:
                self.tt_per_frequency[i % (self.num_of_different_frequencies * self.same_frequency_rep)
                                      // self.same_frequency_rep] += self.tt_S_binning_detuned[i]
            self.Cavity_atom_spectrum += self.tt_per_frequency
            self.Transits_per_freuency += (self.tt_per_frequency != 0).astype(int)

        # building cavity spectrum
        for index, value in enumerate(self.tt_S_binning_detuned):
            if index not in [vec for elem in self.all_transits_seq_indx for vec in elem]:
                self.Cavity_spectrum[index % (self.num_of_different_frequencies * self.same_frequency_rep)
                                     // self.same_frequency_rep] += value

        self.Cavity_spectrum_normalized = (self.Cavity_spectrum / self.power_per_freq_weight)
        self.Cavity_spectrum_normalized = self.Cavity_spectrum_normalized / max(self.Cavity_spectrum_normalized)
        self.Cavity_atom_spectrum_normalized = (self.Cavity_atom_spectrum / self.power_per_freq_weight)
        # self.Cavity_atom_spectrum_normalized = (self.Cavity_atom_spectrum / (self.power_per_freq_weight * self.Transits_per_freuency))
        self.Cavity_atom_spectrum_normalized = self.Cavity_atom_spectrum_normalized / max(self.Cavity_atom_spectrum_normalized)

        self.tt_S_transit_events[[i for i in [vec for elem in self.all_transits_seq_indx for vec in elem]]] += 1
        self.tt_S_transit_events_accumulated = self.tt_S_transit_events_accumulated + self.tt_S_transit_events

        if self.all_transits_seq_indx:
            self.all_transits_indx_batch = self.all_transits_indx_batch[-(N - 1):] + [self.all_transits_seq_indx]

    def get_pulses_location_in_seq(self, delay, seq=Config.QRAM_Exp_Gaussian_samples_S,
                                   smearing=int(Config.num_between_zeros / 2)):
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
                if seq_indx[i] - seq_indx[i - 1] > 1:
                    pulses_loc.append((start_indx - int(smearing), seq_indx[i - 1] + int(smearing)))
                    for j in range(start_indx - int(smearing), seq_indx[i - 1] + int(smearing)):
                        seq_filter_with_smearing[j] = 1
                    start_indx = seq_indx[i]
            pulses_loc.append((start_indx - int(smearing), seq_indx[-1] + int(smearing)))
            for j in range(start_indx - int(smearing), seq_indx[-1] + int(smearing)):
                seq_filter_with_smearing[j] = 1
        return pulses_loc, seq_filter_with_smearing

    def get_avg_num_of_photons_in_seq_pulses(self, seq, pulse_loc, tt_measure):
        avg_num_of_photons_in_seq_pulses = []
        try:
            real_number_of_seq = math.ceil(max(tt_measure) / len(Config.QRAM_Exp_Gaussian_samples_S))
            # print('Real number of seq = %d' %real_number_of_seq)
        except:
            real_number_of_seq = self.number_of_QRAM_sequences
            # print('Max number of seq')
        for t in pulse_loc:
            avg_num_of_photons_in_seq_pulses.append((sum(seq[t[0]:t[1]]) + seq[t[1]]) / (
                        real_number_of_seq * 0.167))  # Sagnac configuiration efficiency 16.7%
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
            box_loc = (np.sum(pulse_loc, axis=1) / 2).astype(int)
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
            + np.array(self.folded_tt_BP[:self.end_of_det_pulse_in_seq]) \
            + np.array(self.folded_tt_DP[:self.end_of_det_pulse_in_seq])

        self.folded_tt_BP_timebins[self.end_of_det_pulse_in_seq:] = self.folded_tt_BP[self.end_of_det_pulse_in_seq:]
        self.folded_tt_DP_timebins[self.end_of_det_pulse_in_seq:] = self.folded_tt_DP[self.end_of_det_pulse_in_seq:]
        if self.pulses_location_in_seq_A or (Config.sprint_pulse_amp_N[0] > 0):
            self.folded_tt_N_directional[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]] = \
                (np.array(
                    self.folded_tt_N_directional[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]])
                 + (Config.sprint_pulse_amp_Early[1]
                    * np.array(self.folded_tt_BP[
                               self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]]))
                 + (Config.sprint_pulse_amp_Early[1]
                    * np.array(self.folded_tt_DP[
                               self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]]))
                 )
            self.folded_tt_BP_timebins[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]] = \
                (np.array(
                    self.folded_tt_BP_timebins[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]])
                 - (Config.sprint_pulse_amp_Early[1]
                    * np.array(self.folded_tt_BP[
                               self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]]))
                 )
            self.folded_tt_DP_timebins[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]] = \
                (np.array(
                    self.folded_tt_DP_timebins[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]])
                 - (Config.sprint_pulse_amp_Early[1]
                    * np.array(self.folded_tt_DP[
                               self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]]))
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

    def plot_transit_figures(self, fig, ax, Num_Of_dets):
        ########################## PLOT!!! ########################################################################

        ax1 = ax[0]
        ax2 = ax[1]
        ax3 = ax[2]
        ax4 = ax[3]
        ax5 = ax[4]
        ax6 = ax[5]

        ax1.clear()
        ax2.clear()
        ax3.clear()
        ax4.clear()
        ax5.clear()
        ax6.clear()

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        textstr_N = r'$Probe_N = %.3f$' % ((len(self.tt_N_measure) * 1000) / self.M_time,) + '[MPhotons/sec]\n' \
                    + '$\overline{Probe}_N = %.3f$' % ((np.mean([len(x) for x in self.tt_N_measure_batch]) * 1000)
                                                       / self.M_time,) + '[MPhotons/sec]'
        textstr_S = r'$Probe_S = %.3f$' % ((len(self.tt_S_measure) * 1000) / self.M_time,) + '[MPhotons/sec]\n' \
                    + '$\overline{Probe}_S = %.3f$' % ((np.mean([len(x) for x in self.tt_S_measure_batch]) * 1000)
                                                       / self.M_time,) + '[MPhotons/sec]'
        textstr_FLR = r'$\overline{FLR}_{MAX} = %.1f$' % (np.mean(self.FLR_measurement) * 1e5,) + r'$\times 10^{-5}$'
        textstr_No_transits = 'NO TRANSITS YET!!!'

        ax1.plot(self.time_bins, self.tt_N_binning, label='Counts histogram', color='b')
        ax1.set_title('North', fontweight="bold")
        ax1.set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
        ax1.text(0.05, 0.95, textstr_N, transform=ax1.transAxes, fontsize=12,
                 verticalalignment='top', bbox=props)
        ax1.legend(loc='upper right')

        ax2.plot(self.time_bins, self.tt_S_binning, label='Counts histogram', color='b')
        ax2.set_title('South', fontweight="bold")
        ax2.set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
        ax2.text(0.05, 0.95, textstr_S, transform=ax2.transAxes, fontsize=12,
                 verticalalignment='top', bbox=props)
        ax2.legend(loc='upper right')

        ax3.plot(self.time_bins, self.tt_N_transit_events, label='Transit events histogram', marker='*', color='k')
        ax3.set(xlabel='Time [msec]', ylabel='Transits events [#Number]')
        ax3.legend(loc='upper right')

        ax4.plot(self.time_bins, self.tt_S_transit_events, label='Transit events histogram', marker='*', color='k')
        ax4.set(xlabel='Time [msec]', ylabel='Transits events [#Number]')
        ax4.legend(loc='upper right')

        if len(self.all_transits_batch) > 0:
            if self.all_transits:
                textstr_transit_counts = r'$N_{Transits} = %s $' % (len(self.all_transits),) + r'$[Counts]$'
            textstr_transit_event_counter = r'$N_{Transits Total} = %s $' % (
            len([vec for elem in self.all_transits_batch for vec in elem]),) + r'$[Counts]$'

        #     ax5.plot(self.t_transit, self.transit_histogram, label='Transit profile', color='b')
        #     ax5.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
        #     ax5.text(0.05, 0.95, textstr_transit_counts, transform=ax5.transAxes, fontsize=12,
        #              verticalalignment='top', bbox=props)
        #
        #     ax6.plot(self.time_bins, self.tt_S_transit_events, label='Transit events histogram', marker='*', color='b')
        #     ax6.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
        #     ax6.text(0.05, 0.95, textstr_transit_event_counter, transform=ax6.transAxes, fontsize=12,
        #              verticalalignment='top', bbox=props)
        # else:
        #     ax5.plot(self.t_transit, self.transit_histogram, label='Transit profile', color='b')
        #     ax5.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
        #     ax5.text(0.25, 0.5, textstr_No_transits, transform=ax5.transAxes, fontsize=24,
        #              verticalalignment='center', bbox=props)
        #
        #     ax6.plot(self.time_bins, self.tt_S_transit_events, label='Transit events histogram', marker='*', color='b')
        #     ax6.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
        #     ax6.text(0.25, 0.5, textstr_No_transits, transform=ax6.transAxes, fontsize=24,
        #              verticalalignment='center', bbox=props)

        # if len(transit_histogram) > 0:
        #     textstr_transit_counts = r'$N_{Transits} = %s $' % (len(all_transits_aligned_first),) + r'$[Counts]$'
        #     textstr_avg_transit_counts = r'$\overline{N}_{Transits} = %.1f $' % (
        #     np.average([len(vec) for vec in all_transits_aligned_first]),) + r'$[Counts]$'
        #
        #     ax5.plot(t_transit, transit_histogram, color='b')
        #     ax5.set_title('Drop transits profile', fontweight="bold")
        #     ax5.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
        #     ax5.text(0.05, 0.95, textstr_transit_counts, transform=ax5.transAxes, fontsize=12,
        #              verticalalignment='top', bbox=props)
        #
        #     ax6.plot(t_transit_batch, transit_histogram_batch, color='b')
        #     ax6.set_title('Accumulated drop transits profile', fontweight="bold")
        #     ax6.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
        #     ax6.text(0.05, 0.95, textstr_avg_transit_counts, transform=ax6.transAxes, fontsize=12,
        #              verticalalignment='top', bbox=props)
        # else:
        #     ax5.plot(t_transit, transit_histogram, color='b')
        #     ax5.set_title('Drop transits profile', fontweight="bold")
        #     ax5.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
        #     ax5.text(0.25, 0.5, textstr_No_transits, transform=ax5.transAxes, fontsize=24,
        #              verticalalignment='center', bbox=props)
        #
        #     ax6.plot(t_transit_batch, transit_histogram_batch, color='b')
        #     ax6.set_title('Accumulated drop transits profile', fontweight="bold")
        #     ax6.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
        #     ax6.text(0.25, 0.5, textstr_No_transits, transform=ax6.transAxes, fontsize=24,
        #              verticalalignment='center', bbox=props)

        ax1.set_ylim(0, 8)
        ax2.set_ylim(0, 8)

        # plt.tight_layout()
        plt.show()
        plt.pause(1)

    def plot_spectrum_figures(self, fig, ax, Num_Of_dets):

        ax[0].clear()
        ax[1].clear()
        ax[2].clear()
        ax[3].clear()
        ax[4].clear()
        ax[5].clear()

        # Detectors display:
        Detectors = []
        for i, det in enumerate(Num_Of_dets):
            Detectors.append(["Det %d" % det, -0.15, 1.9 - i * 0.4, 1])

        # Threshold Box:
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        if self.acquisition_flag:
            flag_color = 'green'
        else:
            flag_color = 'red'
            playsound('C:/Windows/Media/cricket-1.wav')
        props_thresholds = dict(boxstyle='round', edgecolor=flag_color, linewidth=2, facecolor=flag_color, alpha=0.5)

        pause_str = ' , PAUSED!' if self.pause_flag else ''
        textstr_thresholds = '# %d - ' % self.Counter + 'Transmission: %d, ' % self.sum_for_threshold + \
                             'Efficiency: %.2f, ' % self.lockingEfficiency + \
                             'Flr: %.2f, ' % (1000 * np.average(self.FLR_res.tolist())) + \
                             'Lock Error: %.3f' % self.lock_err + pause_str

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        textstr_detuned = r'$S_{detuned} = %.3f$' % (
        (sum(self.tt_S_binning_detuned) * 1000) / (self.M_time / 2),) + '[MPhotons/sec]\n' \
                          + '$\overline{S}_{detuned} = %.3f$' % (
                          (np.mean([sum(x) for x in self.tt_S_binning_detuned_batch]) * 1000)
                          / (self.M_time / 2),) + '[MPhotons/sec]'
        textstr_resonance = r'$S_{res} = %.3f$' % (
        (sum(self.tt_S_binning_resonance) * 1000) / (self.M_time / 2),) + '[MPhotons/sec]\n' \
                            + '$\overline{S}_{res} = %.3f$' % (
                            (np.mean([sum(x) for x in self.tt_S_binning_resonance_batch]) * 1000)
                            / (self.M_time / 2),) + '[MPhotons/sec]'
        textstr_FLR = r'$\overline{FLR}_{MAX} = %.1f$' % (np.mean(self.FLR_measurement) * 1e5,) + r'$\times 10^{-5}$'
        textstr_No_transits = 'NO TRANSITS YET!!!'

        ax[0].plot(self.time_bins[::2], self.S_bins_res_acc, label='Counts histogram', color='b')
        # ax[1].plot([sum(time_bins[i:i+4]) for i in range(0, len(time_bins), 4)],
        #          [sum(tt_S_binning_detuned[i:i+4]) for i in range(0, len(tt_S_binning_detuned), 4)],
        #          label='Counts histogram', color='b')
        ax[0].set_title('On resonance histogram accumulated', fontweight="bold")
        ax[0].set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
        ax[0].text(0.05, 0.95, textstr_resonance, transform=ax[0].transAxes, fontsize=12,
                 verticalalignment='top', bbox=props)
        ax[0].text(0.05, 1.4, textstr_thresholds, transform=ax[0].transAxes, fontsize=28,
                   verticalalignment='top', bbox=props_thresholds)
        ax[0].legend(loc='upper right')

        ax[1].plot(self.time_bins[::2], self.S_bins_detuned_acc, label='Counts histogram', color='b')
        # ax[1].plot([sum(time_bins[i:i+4]) for i in range(0, len(time_bins), 4)],
        #          [sum(tt_S_binning_detuned[i:i+4]) for i in range(0, len(tt_S_binning_detuned), 4)],
        #          label='Counts histogram', color='b')
        ax[1].set_title('Detuned histogram accumulated', fontweight="bold")
        ax[1].set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
        ax[1].text(0.05, 0.95, textstr_detuned, transform=ax[1].transAxes, fontsize=12,
                 verticalalignment='top', bbox=props)
        ax[1].legend(loc='upper right')

        ax[2].plot(self.time_bins[::2], self.tt_S_binning_resonance, label='Counts histogram', color='b')
        # ax[2].plot(self.time_bins[::2], self.tt_S_binning_detuned, label='Counts histogram detuned', color='g')
        # ax[2].plot([sum(time_bins[i:i+4]) for i in range(0, len(time_bins), 4)],
        #          [sum(tt_S_binning_resonance[i:i+4]) for i in range(0, len(tt_S_binning_resonance), 4)],
        #          label='Counts histogram', color='b')
        ax[2].set_title('On resonant counts (live)', fontweight="bold")
        ax[2].set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
        for indx, det_circle in enumerate(Detectors):
            if indx in self.latched_detectors():
                det_color = 'red'
            else:
                det_color = 'green'
            ax[2].text(det_circle[1], det_circle[2], det_circle[0], ha="center", va="center", transform=ax[2].transAxes,
                     bbox=dict(boxstyle=f"circle,pad={det_circle[3]}", edgecolor=det_color, linewidth=2, facecolor=det_color, alpha=0.5))
        ax[2].legend(loc='upper right')

        # ax[3].plot(self.freq_bins/1e6, self.Cavity_atom_spectrum, label='Cavity-atom Spectrum', color='k')
        # ax[3].plot(self.freq_bins/1e6, self.Cavity_spectrum, label='Cavity Spectrum', color='b')
        ax[3].plot(self.freq_bins/1e6, self.Cavity_atom_spectrum_normalized, label='Cavity-atom Spectrum', color='k')
        ax[3].plot(self.freq_bins/1e6, self.Cavity_spectrum_normalized, label='Cavity Spectrum', color='b')
        ax[3].set(xlabel='Frequency [MHz]', ylabel='Transmission Normalized')
        ax[3].legend(loc='upper right')

        ax[4].plot(np.linspace(0, self.histogram_bin_size-1, self.histogram_bin_size), self.folded_tt_S_acc, label='pulses folded', color='k')
        ax[4].plot(np.linspace(0, self.histogram_bin_size-1, self.histogram_bin_size), self.folded_tt_S_acc_2, label='pulses folded_2', color='b')
        ax[4].plot(np.linspace(0, self.histogram_bin_size-1, self.histogram_bin_size), self.folded_tt_S_acc_3, label='pulses folded_2', color='g')
        ax[4].set(xlabel='Time [nsec]', ylabel='Counts [#Number]')
        ax[4].legend(loc='upper right')
        # ax[4].plot(self.freq_bins, self.Cavity_spectrum, label='Cavity Spectrum', color='k')
        # ax[4].set(xlabel='Time [msec]', ylabel='Counts [#Number]')
        # ax[4].legend(loc='upper right')

        # if len(self.all_transits_batch) > 0:
        #     # if self.all_transits:
        #     textstr_transit_counts = r'$N_{Transits} = %s $' % (len(self.all_transits),) + r'$[Counts]$'
        #     textstr_transit_event_counter = r'$N_{Transits Total} = %s $' % (
        #                                     len([vec for elem in self.all_transits_batch for vec in elem]),) \
        #                                     + r'$[Counts]$' + '\n' + textstr_transit_counts

        if len(self.all_transits_indx_batch) > 0:
            # if self.all_transits:
            textstr_transit_counts = r'$N_{Transits} = %s $' % (len(self.all_transits_seq_indx),) + r'$[Counts]$'
            textstr_transit_event_counter = r'$N_{Transits Total} = %s $' % (
                                            len([vec for elem in self.all_transits_indx_batch for vec in elem]),) \
                                            + r'$[Counts]$' + '\n' + textstr_transit_counts

            # ax[5].plot(self.t_transit, self.transit_histogram, label='Transit profile', color='b')
            # ax[5].set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
            # ax[5].text(0.05, 0.95, textstr_transit_counts, transform=ax[5].transAxes, fontsize=12,
            #          verticalalignment='top', bbox=props)

            ax[5].plot(self.time_bins[::2], self.tt_S_transit_events_accumulated, label='Transit events histogram', marker='*', color='b')
            ax[5].set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
            ax[5].text(0.05, 0.95, textstr_transit_event_counter, transform=ax[5].transAxes, fontsize=12,
                       verticalalignment='top', bbox=props)
        else:
            # ax[5].plot(self.t_transit, self.transit_histogram, label='Transit profile', color='b')
            # ax[5].set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
            # ax[5].text(0.25, 0.5, textstr_No_transits, transform=ax[5].transAxes, fontsize=24,
            #          verticalalignment='center', bbox=props)

            ax[5].plot(self.time_bins[::2], self.tt_S_transit_events_accumulated, label='Transit events histogram', marker='*', color='b')
            ax[5].set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
            ax[5].text(0.25, 0.5, textstr_No_transits, transform=ax[5].transAxes, fontsize=24,
                       verticalalignment='center', bbox=props)

        # ax[1].set_ylim(0, 8)
        # ax[2].set_ylim(0, 8)

        # plt.tight_layout()
        plt.show()
        plt.pause(0.5)
        ###########################################################################################################

    def init_params_for_save_sprint(self, Num_Of_dets):
        # define empty variables
        self.tt_measure = []
        self.tt_measure_batch = [[]] * len(Num_Of_dets)
        self.tt_N_measure_batch = []
        self.tt_N_binning_batch = []
        self.tt_N_resonance_batch = []
        self.tt_N_detuned_batch = []
        self.tt_S_measure_batch = []
        self.tt_S_binning_batch = []
        self.tt_S_binning_detuned_batch = []
        self.tt_S_binning_resonance_batch = []
        self.tt_N_binning_detuned_batch = []
        self.tt_N_binning_resonance_batch = []
        self.tt_S_resonance_batch = []
        self.tt_S_detuned_batch = []
        self.tt_BP_measure_batch = []
        self.tt_DP_measure_batch = []
        self.tt_FS_measure_batch = []
        self.all_transits_batch = []
        self.all_transits_indx_batch = []
        self.transit_histogram_batch = []
        self.all_transits_aligned_first_batch = []
        self.transit_histogram_batch = []
        self.FLR_measurement = []
        self.Exp_timestr_batch = []

    def AOM_power_per_freq_calibration(self, dirname):
        if dirname is None:
            return np.full((1, self.spectrum_bin_number), 1)[0]
        if not os.path.exists(dirname):
            return np.full((1, self.spectrum_bin_number), 1)[0]
        self.calibration_spectrum = np.load(dirname)
        return self.calibration_spectrum['arr_0']/max(self.calibration_spectrum['arr_0'])

    # def Save_SNSPDs_Transit_Measurement_with_tt(self, N, histogram_bin_size, transit_profile_bin_size, preComment,
    #                                             total_counts_threshold, transit_counts_threshold, FLR_threshold,
    #                                             lock_err_threshold, exp_flag, with_atoms):
    #     """
    #     Function for analyzing,saving and displaying data from sprint experiment.
    #     :param N: Number of maximum experiments (free throws) saved and displayed.
    #              program we are looking for transits of atoms next to the toroid and record them.
    #     :param qram_sequence_len: the number of sprint sequences (detection and sprint pulses combination)
    #     :param preComment: The comment added at the start of the experiment, usually consisting of unique experiment
    #                        parameters.
    #     :param transit_condition:
    #     :param lock_err_threshold: The maximal error in locking resonator to rb line (in ms, depends on the scan width/vpp)
    #     :param max_probe_counts: The maximum counts probe counts at each direction measured when cavity isn't locked.
    #     :param filter_delay: Delay of the filter window compared to original location (can be taken by comparing the
    #                          time of the 1st detection pulse pick location to the original sequence)
    #     :param reflection_threshold:
    #     :param reflection_threshold_time:
    #     :return:
    #     """
    #     # if not preComment:
    #     #     preComment = pymsgbox.prompt('Add comment to measurement: ', default='', timeout=int(30e3))
    #
    #     # set constant parameters for the function
    #     Num_Of_dets = [1, 2, 3, 4, 5, 6, 7, 8]
    #     self.histogram_bin_size = histogram_bin_size
    #     histogram_bin_number = self.M_time // histogram_bin_size
    #     self.time_bins = np.linspace(0, self.M_time, histogram_bin_number)
    #     self.Transit_profile_bin_size = transit_profile_bin_size
    #     time_threshold = int(
    #         histogram_bin_size * 0.8)  # The minimum time between two time tags to be counted for a transit. # TODO: might need a factor of 2???
    #
    #     ## Listen for keyboard
    #     self.listener = keyboard.Listener(on_press=self.on_press, on_release=self.on_release)
    #     self.listener.start()
    #     self.keyPress = None
    #     print('\033[94m' + 'Press ESC to stop measurement.' + '\033[0m')  # print blue
    #
    #     self.init_params_for_save_sprint(Num_Of_dets)
    #
    #     ####     get tt and counts from OPX to python   #####
    #     Counts_handle, tt_handle, FLR_handle = self.get_handles_from_OPX_Server(Num_Of_dets)
    #
    #     ## take data only if
    #     start = True
    #     # take threshold from npz ( error from resonator lock PID)
    #     self.lock_err = 0.001  # ZIV change back to None
    #     if exp_flag:
    #         while self.lock_err == None:
    #             try:
    #                 self.lock_err = np.abs(np.load(
    #                     'U:\Lab_2023\Experiment_results\QRAM\Locking_PID_Error\locking_err.npy',
    #                     allow_pickle=True))  # the error of locking the resontor to Rb line
    #             except:
    #                 print('error in loading file')
    #             time.sleep(0.01)  # in seconds
    #     else:
    #         self.lock_err = lock_err_threshold / 2
    #
    #     self.sum_for_threshold = total_counts_threshold
    #     cycle = 0
    #     # Place holders for results # TODO: ask dor - is it works as we expect?
    #     while ((self.lock_err > lock_err_threshold) or
    #            (self.sum_for_threshold >= total_counts_threshold) and exp_flag) or start:
    #         if self.keyPress == 'ESC':
    #             print('\033[94m' + 'ESC pressed. Stopping measurement.' + '\033[0m')  # print blue
    #             self.updateValue("Sprint_Exp_switch", False)
    #             self.MOT_switch(True)
    #             self.update_parameters()
    #             self.Stop_run_daily_experiment = True
    #             # Other actions can be added here
    #             break
    #         if start:
    #             if cycle > 2:
    #                 start = False
    #             else:
    #                 cycle += 1
    #         else:
    #             print('Above Threshold')
    #         self.get_tt_from_handles(Num_Of_dets, Counts_handle, tt_handle, FLR_handle)
    #         # sum_for_threshold [Mcounts /sec] = (number of photons * 10^-6 [Mcounts] / Measuring time [nsec] * 10^-9 [sec/nsec]):
    #         self.sum_for_threshold = (len(self.tt_S_directional_measure) * 1000) / self.M_time
    #         print(self.lock_err, self.lock_err > lock_err_threshold, self.sum_for_threshold)
    #         if exp_flag:
    #             try:
    #                 self.lock_err = np.abs(np.load(
    #                     'U:\Lab_2023\Experiment_results\QRAM\Locking_PID_Error\locking_err.npy'))  # the error of locking the resontor to Rb line
    #             except:
    #                 pass
    #                 print('error in loading file')
    #     ####    end get tt and counts from OPX to python   #####
    #
    #     self.tt_N_binning = np.zeros(histogram_bin_number)
    #     self.tt_S_binning = np.zeros(histogram_bin_number)
    #     self.tt_N_transit_events = np.zeros(histogram_bin_number)
    #     self.tt_S_transit_events = np.zeros(histogram_bin_number)
    #     self.tt_S_transit_events_accumulated = np.zeros(histogram_bin_number)
    #     self.all_transits_accumulated = np.zeros(histogram_bin_number)
    #
    #     # fold reflections and transmission
    #     for x in self.tt_N_directional_measure:
    #         self.tt_N_binning[(x - 1) // histogram_bin_size] += 1
    #     for x in self.tt_S_directional_measure:
    #         self.tt_S_binning[(x - 1) // histogram_bin_size] += 1
    #     self.tt_S_binning_avg = self.tt_S_binning
    #
    #     if len(self.tt_S_measure_batch) == N:
    #         self.tt_N_transit_events[
    #             [i for i, x in enumerate(self.tt_N_binning_batch[0]) if x > transit_counts_threshold]] -= 1
    #         self.tt_S_transit_events[
    #             [i for i, x in enumerate(self.tt_S_binning_batch[0]) if x > transit_counts_threshold]] -= 1
    #
    #     self.tt_N_measure_batch = self.tt_N_measure_batch[-(N - 1):] + [self.tt_N_directional_measure]
    #     self.tt_N_binning_batch = self.tt_N_binning_batch[-(N - 1):] + [self.tt_N_binning]
    #     self.tt_S_measure_batch = self.tt_S_measure_batch[-(N - 1):] + [self.tt_S_directional_measure]
    #     self.tt_S_binning_batch = self.tt_S_binning_batch[-(N - 1):] + [self.tt_S_binning]
    #
    #     ## record time
    #     timest = time.strftime("%Y%m%d-%H%M%S")
    #     datest = time.strftime("%Y%m%d")
    #     FLR_measurement = []
    #     Exp_timestr_batch = []
    #     lock_err_batch = []
    #
    #     FLR_measurement = FLR_measurement[-(N - 1):] + [self.FLR_res.tolist()]
    #     Exp_timestr_batch = Exp_timestr_batch[-(N - 1):] + [timest]
    #     lock_err_batch = lock_err_batch[-(N - 1):] + [self.lock_err]
    #
    #     self.tt_N_transit_events[[i for i, x in enumerate(self.tt_N_binning) if x > transit_counts_threshold]] += 1
    #     self.tt_S_transit_events[[i for i, x in enumerate(self.tt_S_binning) if x > transit_counts_threshold]] += 1
    #     self.tt_S_transit_events_accumulated = self.tt_S_transit_events
    #
    #     self.find_transit_events_transitexp(N, transit_time_threshold=time_threshold,
    #                              transit_counts_threshold=transit_counts_threshold)
    #
    #     self.Counter = 1  # Total number of successful cycles
    #     self.repitions = 1  # Total number of cycles
    #     self.acquisition_flag = True
    #     self.threshold_flag = True
    #     self.pause_flag = False
    #     #
    #     # create figures template
    #     fig = plt.figure()
    #     ax1 = plt.subplot2grid((6, 2), (0, 0), colspan=1, rowspan=2)
    #     ax2 = plt.subplot2grid((6, 2), (0, 1), colspan=1, rowspan=2)
    #     ax3 = plt.subplot2grid((6, 2), (2, 0), colspan=1, rowspan=2)
    #     ax4 = plt.subplot2grid((6, 2), (2, 1), colspan=1, rowspan=2)
    #     ax5 = plt.subplot2grid((6, 2), (4, 0), colspan=1, rowspan=2)
    #     ax6 = plt.subplot2grid((6, 2), (4, 1), colspan=1, rowspan=2)
    #
    #     figManager = plt.get_current_fig_manager()
    #     figManager.window.showMaximized()
    #
    #     while True:
    #         if self.keyPress == 'ESC':
    #             print('\033[94m' + 'ESC pressed. Stopping measurement.' + '\033[0m')  # print blue
    #             self.updateValue("QRAM_Exp_switch", False)
    #             self.MOT_switch(True)
    #             self.Stop_run_daily_experiment = True
    #             self.update_parameters()
    #             # Other actions can be added here
    #             break
    #         if self.keyPress == 'SPACE' and not self.pause_flag:
    #             print('\033[94m' + 'SPACE pressed. Pausing measurement.' + '\033[0m')  # print blue
    #             self.pause_flag = True
    #             self.keyPress = None
    #             # Other actions can be added here
    #         if self.keyPress == 'SPACE' and self.pause_flag:
    #             print('\033[94m' + 'SPACE pressed. Continuing measurement.' + '\033[0m')  # print blue
    #             self.pause_flag = False
    #             self.keyPress = None
    #             # Other actions can be added here
    #         self.lockingEfficiency = self.Counter / self.repitions
    #         print(timest, self.Counter, 'Eff: %.2f' % self.lockingEfficiency,
    #               'Flr: %.2f' % (1000 * np.average(self.FLR_res.tolist())))
    #         self.repitions += 1
    #         ######################################## PLOT!!! ###########################################################
    #         ax = [ax1, ax2, ax3, ax4, ax5, ax6]
    #         self.plot_transit_figures(fig, ax, Num_Of_dets)
    #         ############################################################################################################
    #         if self.Counter == N:
    #             print(
    #                 '\033[94m' + f'finished {N} Runs, {"with" if with_atoms else "without"} atoms' + '\033[0m')  # print blue
    #             # Other actions can be added here
    #             break
    #         count = 1
    #         while True:
    #             # record time:
    #             timest = time.strftime("%H%M%S")
    #             datest = time.strftime("%Y%m%d")
    #             self.get_tt_from_handles(Num_Of_dets, Counts_handle, tt_handle, FLR_handle)
    #
    #             # Check if new tt's arrived:
    #             lenS = min(len(self.tt_S_directional_measure), len(self.tt_S_measure_batch[-1]))
    #
    #             # Check if the number of same values in the new and last vector are less than 1/2 of the total number of values.
    #             is_new_tts_S = sum(np.array(self.tt_S_directional_measure[:lenS]) == np.array(
    #                 self.tt_S_measure_batch[-1][:lenS])) < lenS / 2
    #
    #             print(count)
    #             count += 1
    #             time.sleep(0.01)
    #             if is_new_tts_S:
    #                 break
    #             if self.keyPress == 'ESC':
    #                 print('\033[94m' + 'ESC pressed. Stopping measurement.' + '\033[0m')  # print blue
    #                 self.updateValue("QRAM_Exp_switch", False)
    #                 self.MOT_switch(True)
    #                 self.update_parameters()
    #                 # Other actions can be added here
    #                 break
    #         # assaf - if x=self.M_window the index is out of range so i added 1
    #         try:
    #             self.lock_err = np.abs(np.load(
    #                 'U:\Lab_2023\Experiment_results\QRAM\Locking_PID_Error\locking_err.npy'))  # the error of locking the resontor to Rb line
    #         except:
    #             pass
    #             print('error in loading file')
    #
    #         self.sum_for_threshold = (len(self.tt_S_directional_measure) * 1000) / self.M_time
    #
    #         if exp_flag:
    #             if (self.lock_err > lock_err_threshold) or \
    #                     ((1000 * np.average(self.FLR_res.tolist()) < FLR_threshold) and with_atoms) or \
    #                     self.latched_detectors():
    #                 self.acquisition_flag = False
    #             else:
    #                 self.acquisition_flag = True
    #
    #         self.threshold_flag = (self.sum_for_threshold < total_counts_threshold)
    #
    #         if (self.threshold_flag or not exp_flag) and self.acquisition_flag and not self.pause_flag:
    #
    #             if self.Counter < N:
    #                 self.Counter += 1
    #
    #             self.tt_N_binning = np.zeros(histogram_bin_number)
    #             self.tt_S_binning = np.zeros(histogram_bin_number)
    #             self.tt_S_transit_events = np.zeros(histogram_bin_number)
    #
    #             self.tt_S_binning_avg = self.tt_S_binning_avg * (1 - 1 / self.Counter)
    #
    #             for x in self.tt_N_directional_measure:
    #                 self.tt_N_binning[(x - 1) // histogram_bin_size] += 1
    #             for x in self.tt_S_directional_measure:
    #                 self.tt_S_binning[(x - 1) // histogram_bin_size] += 1
    #                 self.tt_S_binning_avg[(x - 1) // histogram_bin_size] += 1 / self.Counter
    #
    #             if len(self.tt_S_measure_batch) == N:
    #                 self.tt_N_transit_events[
    #                     [i for i, x in enumerate(self.tt_N_binning_batch[0]) if x > transit_counts_threshold]] -= 1
    #                 self.tt_S_transit_events[
    #                     [i for i, x in enumerate(self.tt_S_binning_batch[0]) if x > transit_counts_threshold]] -= 1
    #
    #             self.tt_N_measure_batch = self.tt_N_measure_batch[-(N - 1):] + [self.tt_N_directional_measure]
    #             self.tt_N_binning_batch = self.tt_N_binning_batch[-(N - 1):] + [self.tt_N_binning]
    #             self.tt_S_measure_batch = self.tt_S_measure_batch[-(N - 1):] + [self.tt_S_directional_measure]
    #             self.tt_S_binning_batch = self.tt_S_binning_batch[-(N - 1):] + [self.tt_S_binning]
    #
    #             self.tt_N_transit_events[
    #                 [i for i, x in enumerate(self.tt_N_binning) if x > transit_counts_threshold]] += 1
    #             self.tt_S_transit_events[
    #                 [i for i, x in enumerate(self.tt_S_binning) if x > transit_counts_threshold]] += 1
    #             self.tt_S_transit_events_accumulated = self.tt_S_transit_events_accumulated + self.tt_S_transit_events
    #
    #             self.find_transit_events(N, minimum_transit_time=time_threshold,
    #                                      total_transit_reflections=transit_counts_threshold)
    #
    #             FLR_measurement = FLR_measurement[-(N - 1):] + [self.FLR_res.tolist()]
    #             Exp_timestr_batch = Exp_timestr_batch[-(N - 1):] + [timest]
    #             lock_err_batch = lock_err_batch[-(N - 1):] + [self.lock_err]
    #             self.save_tt_to_batch(Num_Of_dets, N)
    #
    #     ############################################## END WHILE LOOP #################################################
    #
    #     # For debuging:
    #     # self.most_common([x for vec in self.tt_S_measure_batch + self.tt_N_measure_batch for x in vec])
    #
    #     ## Adding comment to measurement [prompt whether stopped or finished regularly]
    #     if exp_flag:
    #         if self.Counter < N:
    #             aftComment = pymsgbox.prompt('Add comment to measurement: ', default='', timeout=int(30e3))
    #         else:
    #             aftComment = ''
    #     else:
    #         aftComment = 'ignore'
    #
    #     # if aftComment == 'Timeout': aftComment = None
    #
    #     #### Handle file-names and directories #####
    #     ## Saving: np.savez(filedir, data = x) #note: @filedir is string of full directory; data is the queyword used to read @x from the file:
    #     ## Loading: file = np.load(f, allow_pickle = True)
    #     ##          x = file['data']
    #
    #     #### ------ Save results ------
    #     #  -------   Create dir
    #     root_dirname = f'U:\\Lab_2023\\Experiment_results\\QRAM\\{datest}\\'
    #     dirname = root_dirname + f'{timest}_Photon_TimeTags\\'  # Specific experiment dir
    #     dirname_Det = dirname + 'AllDetectors\\'
    #     dirname_N = dirname + 'North(8)\\'
    #     dirname_S = dirname + 'South(5)\\'
    #     dirname_D = dirname + 'Dark(3,4)\\'
    #     dirname_B = dirname + 'Bright(1,2)\\'
    #     dirname_FS = dirname + 'FastSwitch(6,7)\\'
    #     dirname_expfigures = root_dirname + 'Experiment_figures\\'
    #     if not os.path.exists(dirname):
    #         os.makedirs(dirname)
    #     if not os.path.exists(dirname_Det):
    #         os.makedirs(dirname_Det)
    #     if not os.path.exists(dirname_N):
    #         os.makedirs(dirname_N)
    #     if not os.path.exists(dirname_S):
    #         os.makedirs(dirname_S)
    #     if not os.path.exists(dirname_D):
    #         os.makedirs(dirname_D)
    #     if not os.path.exists(dirname_B):
    #         os.makedirs(dirname_B)
    #     if not os.path.exists(dirname_FS):
    #         os.makedirs(dirname_FS)
    #     if not os.path.exists(dirname_expfigures):
    #         os.makedirs(dirname_expfigures)
    #
    #     # ----  msmnt files names  -----
    #     # Counter_str = (self.Counter)
    #     filename_Det_tt = []
    #     for i in Num_Of_dets:
    #         filename_Det_tt.append(f'Det' + str(i) + f'_timetags.npz')
    #     filename_S_tt = f'South_timetags.npz'
    #     filename_N_tt = f'North_timetags.npz'
    #     filename_DP_tt = f'Dark_timetags.npz'
    #     filename_BP_tt = f'Bright_timetags.npz'
    #     filename_FS_tt = f'FS_timetags.npz'
    #     filename_FLR = f'Flouresence.npz'
    #     filename_timestamp = f'Drops_time_stamps.npz'
    #     filename_lockerr = f'Cavity_Lock_Error.npz'
    #     filname_sequence_S = f'South_sequence_vector.npz'
    #     filname_sequence_N = f'North_sequence_vector.npz'
    #     filname_sequence_Early = f'Early_sequence_vector.npz'
    #     filname_sequence_Late = f'Late_sequence_vector.npz'
    #     filname_sequence_FS = f'FS_sequence_vector.npz'
    #     filename_transits_events = f'seq_transit_events_batched.npz'
    #     filename_transits_batched = f'all_transits_seq_indx_batch.npz'
    #     filename_experimentPlot = f'Experiment_plot.png'
    #     filename_experimentfigure = f'{timest}_Experiment_Figure.png'
    #
    #     if len(FLR_measurement) > 0:
    #         np.savez(dirname + filename_FLR, FLR_measurement)
    #     if len(Exp_timestr_batch) > 0:
    #         np.savez(dirname + filename_timestamp, Exp_timestr_batch)
    #         np.savez(dirname + filename_lockerr, lock_err_batch)
    #         np.savez(dirname + filname_sequence_S, Config.QRAM_Exp_Gaussian_samples_S)
    #         np.savez(dirname + filname_sequence_N, Config.QRAM_Exp_Gaussian_samples_N)
    #         np.savez(dirname + filname_sequence_Early, Config.QRAM_Exp_Square_samples_Early)
    #         np.savez(dirname + filname_sequence_Late, Config.QRAM_Exp_Square_samples_Late)
    #         np.savez(dirname + filname_sequence_FS, (1 - np.array(Config.QRAM_Exp_Square_samples_FS)).tolist())
    #         plt.savefig(dirname + filename_experimentPlot, bbox_inches='tight')
    #         plt.savefig(dirname_expfigures + filename_experimentfigure, bbox_inches='tight')
    #     for i in range(len(Num_Of_dets)):
    #         if len(self.tt_measure_batch[i]) > 0:
    #             np.savez(dirname_Det + filename_Det_tt[i], self.tt_measure_batch[i])
    #     if len(self.tt_S_measure_batch) > 0:
    #         np.savez(dirname_S + filename_S_tt, self.tt_S_measure_batch)
    #     if len(self.tt_N_measure_batch) > 0:
    #         np.savez(dirname_N + filename_N_tt, self.tt_N_measure_batch)
    #     if len(self.tt_DP_measure_batch) > 0:
    #         np.savez(dirname_D + filename_DP_tt, self.tt_DP_measure_batch)
    #     if len(self.tt_BP_measure_batch) > 0:
    #         np.savez(dirname_B + filename_BP_tt, self.tt_BP_measure_batch)
    #     if len(self.tt_FS_measure_batch) > 0:
    #         np.savez(dirname_FS + filename_FS_tt, self.tt_FS_measure_batch)
    #     if len(self.all_transits_batch) > 0:
    #         np.savez(dirname + filename_transits_events, self.all_transits_batch)
    #     # if len(all_transits_batch) > 0:
    #     #     np.savez(dirname_S + filename_S_transits, all_transits_batch)
    #
    #     ### Edit comments file ####
    #     cmntDir = os.path.join(root_dirname, 'daily_experiment_comments.csv')
    #     cmnt_header = 'Date,Time,IgnoreValid,Atoms,Cycles,Comment'
    #     if not os.path.exists(cmntDir):
    #         # Write header line
    #         try:
    #             with open(cmntDir, "a") as commentsFile:
    #                 commentsFile.write(cmnt_header + '\n')
    #         except:
    #             print('Could not save comments, error writing to comments-file.')
    #
    #     if preComment is not None: cmnt = preComment + '; '
    #     if aftComment is not None: cmnt = preComment + aftComment
    #     if preComment is None and aftComment is None: cmnt = 'No comment.'
    #     if 'ignore' in cmnt:
    #         experiment_success = 'ignore'
    #     else:
    #         experiment_success = 'valid'
    #     full_line = f'{datest},{timest},{experiment_success},{with_atoms},{self.Counter},{cmnt}'
    #     try:
    #         with open(cmntDir, "a") as commentsFile:
    #             commentsFile.write(full_line + '\n')
    #     except:
    #         print('Could not save comments, error writing to comments-file.')
    #
    #     experiment_cmnt = 'transit condition: minimum time between reflection = ' + str(time_threshold) + \
    #                       'with at least ' + str(transit_counts_threshold) + ' reflections ; ' + \
    #                       'reflection threshold: ' + str(total_counts_threshold) + 'MCounts / sec'
    #
    #     comments = {'comments': experiment_cmnt}
    #     try:
    #         with open(f'{dirname}experiment_comments.txt', 'w') as file:
    #             json.dump(comments, file, indent=4)
    #     except Exception:
    #         pass
    #
    #     self.saveConfigTable(path=dirname)
    #     for qrdCtrl in self.QuadRFControllers:
    #         qrdCtrl.saveLinesAsCSV(f'{dirname}QuadRF_table.csv')
    #
    #     self.updateValue("QRAM_Exp_switch", False)
    #     self.update_parameters()
    #
    #     ## ------------------ end of saving section -------

    def Save_SNSPDs_Spectrum_Measurement_with_tt(self, N, transit_profile_bin_size, preComment, transit_cond,
                                                 total_counts_threshold, transit_counts_threshold, FLR_threshold,
                                                 lock_err_threshold, exp_flag, with_atoms, calibration_dirname):
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
        # if not preComment:
        #     preComment = pymsgbox.prompt('Add comment to measurement: ', default='', timeout=int(30e3))

        # set constant parameters for the function
        Num_Of_dets = [1, 2, 3, 4, 5, 6, 7, 8]
        self.histogram_bin_size = Config.frequency_sweep_duration * 2
        self.sequence_duration = Config.frequency_sweep_duration * 2
        self.histogram_bin_number = self.M_time // self.histogram_bin_size # num of bins in cycle of frequency sweep
        self.time_bins = np.linspace(0, self.M_time, self.histogram_bin_number*2)
        self.spectrum_bin_number = self.num_of_different_frequencies
        self.freq_bins = np.linspace((self.frequency_start - Config.IF_AOM_Spectrum) * 2,
                                     (self.frequency_start + self.spectrum_bandwidth - Config.IF_AOM_Spectrum) * 2,
                                     self.spectrum_bin_number)
        self.Transit_profile_bin_size = transit_profile_bin_size

        time_threshold = int(self.histogram_bin_size * 1.6)  # The minimum time between two time tags to be counted for a transit. # TODO: might need a factor of 2???

        ## Listen for keyboard
        self.listener = keyboard.Listener(on_press=self.on_press, on_release=self.on_release)
        self.listener.start()
        self.keyPress = None
        print('\033[94m' + 'Press ESC to stop measurement.' + '\033[0m')  # print blue

        self.init_params_for_save_sprint(Num_Of_dets)

        ####     get tt and counts from OPX to python   #####
        Counts_handle, tt_handle, FLR_handle = self.get_handles_from_OPX_Server(Num_Of_dets)

        ## take data only if
        start = True
        # take threshold from npz ( error from resonator lock PID)
        self.lock_err = 0.001  # ZIV change back to None
        if exp_flag:
            while self.lock_err == None:
                try:
                    self.lock_err = np.abs(np.load(
                        'U:\Lab_2023\Experiment_results\QRAM\Locking_PID_Error\locking_err.npy',
                        allow_pickle=True))  # the error of locking the resontor to Rb line
                except:
                    print('error in loading file')
                time.sleep(0.01)  # in seconds
        else:
            self.lock_err = lock_err_threshold / 2
            self.sum_for_threshold = total_counts_threshold

        self.sum_for_threshold = total_counts_threshold # is the resonance criticaly coupled enough
        cycle = 0
        # Place holders for results # TODO: ask dor - is it works as we expect?
        while ((self.lock_err > lock_err_threshold) or
               (self.sum_for_threshold >= total_counts_threshold) and exp_flag) or start:
            if self.keyPress == 'ESC':
                print('\033[94m' + 'ESC pressed. Stopping measurement.' + '\033[0m')  # print blue
                self.updateValue("Spectrum_Exp_switch", False)
                self.MOT_switch(True)
                self.update_parameters()
                self.Stop_run_daily_experiment = True
                # Other actions can be added here
                break
            if start:
                if cycle > 2:
                    start = False
                else:
                    cycle += 1
            else:
                print('Above Threshold')
                print(self.sum_for_threshold)
            self.get_tt_from_handles(Num_Of_dets, Counts_handle, tt_handle, FLR_handle)
            self.tt_S_binning = np.zeros(self.histogram_bin_number * 2)

            for x in self.tt_S_no_gaps:
                self.tt_S_binning[(x - 1) // Config.frequency_sweep_duration] += 1
            self.tt_S_binning_resonance = [self.tt_S_binning[x] for x in range(len(self.tt_S_binning)) if
                                           not x % 2]  # even



            print(self.lock_err, self.lock_err > lock_err_threshold, self.sum_for_threshold)
            if exp_flag:
                try:
                    self.lock_err = np.abs(np.load(
                        'U:\Lab_2023\Experiment_results\QRAM\Locking_PID_Error\locking_err.npy'))  # the error of locking the resontor to Rb line
                    # sum_for_threshold [Mcounts /sec] = (number of photons * 10^-6 [Mcounts] / Measuring time [nsec] * 10^-9 [sec/nsec]):
                    self.sum_for_threshold = (np.sum(self.tt_S_binning_resonance) * 1000) / self.M_time
                except:
                    print('error in loading file')
                    pass
        ####    end get tt and counts from OPX to python   #####

        ## record time
        timest = time.strftime("%Y%m%d-%H%M%S")
        datest = time.strftime("%Y%m%d")
        FLR_measurement = []
        Exp_timestr_batch = []
        lock_err_batch = []

        # initilaization
        self.tt_N_binning = np.zeros(self.histogram_bin_number * 2)
        self.tt_N_transit_events = np.zeros(self.histogram_bin_number)
        self.tt_S_transit_events = np.zeros(self.histogram_bin_number)
        self.tt_S_transit_events_accumulated = np.zeros(self.histogram_bin_number)
        self.all_transits_accumulated = np.zeros(self.histogram_bin_number)
        self.folded_tt_S_acc = np.zeros(self.histogram_bin_size, dtype=int)
        self.folded_tt_S_acc_2 = np.zeros(self.histogram_bin_size, dtype=int)
        self.folded_tt_S_acc_3 = np.zeros(self.histogram_bin_size, dtype=int)


        self.Cavity_atom_spectrum = np.zeros(self.spectrum_bin_number)
        self.Transits_per_freuency = np.zeros(self.spectrum_bin_number)
        self.Cavity_atom_spectrum_normalized = np.zeros(self.spectrum_bin_number)
        self.Cavity_spectrum = np.zeros(self.spectrum_bin_number)
        self.Cavity_spectrum_normalized = np.zeros(self.spectrum_bin_number)

        # calibration_dirname = 'U:\\Lab_2023\\Experiment_results\\QRAM\\20231109\\160748_Photon_TimeTags\\Cavity_spectrum.npz'
        # calibration_dirname = r'U:\Lab_2023\Experiment_results\QRAM\20231114\104323_Photon_TimeTags\Cavity_spectrum.npz'  # Widening spectrum
        calibration_dirname = r'U:\Lab_2023\Experiment_results\QRAM\20231116\133237_Photon_TimeTags\Cavity_spectrum.npz'  # Widening spectrum
        # calibration_dirname = None  # Set to None when you are running calibration process
        self.power_per_freq_weight = self.AOM_power_per_freq_calibration(calibration_dirname)

        # fold reflections and transmission
        for x in self.tt_N_directional_measure:
            self.tt_N_binning[(x - 1) // Config.frequency_sweep_duration] += 1 #TODO: x-1?? in spectrum its x
        self.tt_N_binning_avg = self.tt_N_binning

        for x in self.tt_S_no_gaps:
            if x < 2e6:
                self.folded_tt_S_acc[(x - 1) % self.histogram_bin_size] += 1
            if 6e6 < x < 8e6:
                self.folded_tt_S_acc_2[(x - 1) % self.histogram_bin_size] += 1
            if 8e6 < x < 10e6:
                self.folded_tt_S_acc_3[(x - 1) % self.histogram_bin_size] += 1
        self.tt_S_binning_avg = self.tt_S_binning

        # split the binning vector to odd and even - on and off resonance pulses
        self.tt_N_binning_detuned = [self.tt_N_binning[x] for x in range(len(self.tt_N_binning)) if x % 2]  # odd
        self.tt_N_binning_resonance = [self.tt_N_binning[x] for x in range(len(self.tt_N_binning)) if not x % 2]  # even
        self.tt_S_binning_detuned = [self.tt_S_binning[x] for x in range(len(self.tt_S_binning)) if x % 2]  # odd

        # put into batches
        self.tt_N_measure_batch = self.tt_N_measure_batch[-(N - 1):] + [self.tt_N_directional_measure]
        self.tt_N_binning_batch = self.tt_N_binning_batch[-(N - 1):] + [self.tt_N_binning]
        self.tt_S_measure_batch = self.tt_S_measure_batch[-(N - 1):] + [self.tt_S_no_gaps]
        self.tt_S_binning_batch = self.tt_S_binning_batch[-(N - 1):] + [self.tt_S_binning]
        self.tt_S_binning_resonance_batch = self.tt_S_binning_resonance_batch[-(N - 1):] + [self.tt_S_binning_resonance]
        self.tt_S_binning_detuned_batch = self.tt_S_binning_detuned_batch[-(N - 1):] + [self.tt_S_binning_detuned]
        self.S_bins_res_acc = np.sum(np.array(self.tt_S_binning_resonance_batch),
                                     0)  # tt_S_binning_resonance accumulated (sum over the batch)
        self.S_bins_detuned_acc = np.sum(np.array(self.tt_S_binning_detuned_batch),
                                     0)  # tt_S_binning_resonance accumulated (sum over the batch)

        FLR_measurement = FLR_measurement[-(N - 1):] + [self.FLR_res.tolist()]
        Exp_timestr_batch = Exp_timestr_batch[-(N - 1):] + [timest]
        lock_err_batch = lock_err_batch[-(N - 1):] + [self.lock_err]
        self.save_tt_to_batch(Num_Of_dets, N)

        # self.find_transit_events_spectrum(N, transit_time_threshold=time_threshold,
        #                                   transit_counts_threshold=transit_counts_threshold)
        self.find_transits_events_spectrum_exp(self.tt_S_binning_resonance, N, transit_cond)

        self.Counter = 1  # Total number of successful cycles
        self.repitions = 1  # Total number of cycles
        self.acquisition_flag = True
        self.threshold_flag = True
        self.pause_flag = False
        #
        # create figures template
        fig = plt.figure()
        ax1 = plt.subplot2grid((6, 2), (0, 0), colspan=1, rowspan=2)
        ax2 = plt.subplot2grid((6, 2), (0, 1), colspan=1, rowspan=2)
        ax3 = plt.subplot2grid((6, 2), (2, 0), colspan=1, rowspan=2)
        ax4 = plt.subplot2grid((6, 2), (2, 1), colspan=1, rowspan=2)
        ax5 = plt.subplot2grid((6, 2), (4, 0), colspan=1, rowspan=2)
        ax6 = plt.subplot2grid((6, 2), (4, 1), colspan=1, rowspan=2)

        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()

        while True:
            if self.keyPress == 'ESC':
                print('\033[94m' + 'ESC pressed. Stopping measurement.' + '\033[0m')  # print blue
                self.updateValue("QRAM_Exp_switch", False)
                self.MOT_switch(True)
                self.Stop_run_daily_experiment = True
                self.update_parameters()
                # Other actions can be added here
                break
            if self.keyPress == 'SPACE' and not self.pause_flag:
                print('\033[94m' + 'SPACE pressed. Pausing measurement.' + '\033[0m')  # print blue
                self.pause_flag = True
                self.keyPress = None
                # Other actions can be added here
            if self.keyPress == 'SPACE' and self.pause_flag:
                print('\033[94m' + 'SPACE pressed. Continuing measurement.' + '\033[0m')  # print blue
                self.pause_flag = False
                self.keyPress = None
                # Other actions can be added here
            self.lockingEfficiency = self.Counter / self.repitions
            print(timest, self.Counter, 'Eff: %.2f' % self.lockingEfficiency,
                  'Flr: %.2f' % (1000 * np.average(self.FLR_res.tolist())))
            self.repitions += 1
            ######################################## PLOT!!! ###########################################################
            ax = [ax1, ax2, ax3, ax4, ax5, ax6]
            self.plot_spectrum_figures(fig, ax, Num_Of_dets)
            ############################################################################################################
            if self.Counter == N:
                print(
                    '\033[94m' + f'finished {N} Runs, {"with" if with_atoms else "without"} atoms' + '\033[0m')  # print blue
                # Other actions can be added here
                break
            count = 1
            while True:
                # record time:
                timest = time.strftime("%H%M%S")
                datest = time.strftime("%Y%m%d")
                self.get_tt_from_handles(Num_Of_dets, Counts_handle, tt_handle, FLR_handle)

                # Check if new tt's arrived:
                lenS = min(len(self.tt_S_no_gaps), len(self.tt_S_measure_batch[-1]))

                # Check if the number of same values in the new and last vector are less than 1/2 of the total number of values.
                is_new_tts_S = sum(np.array(self.tt_S_no_gaps[:lenS]) == np.array(
                    self.tt_S_measure_batch[-1][:lenS])) < lenS / 2

                # make sure that measurement did not fill 95% or more the buffer
                self.tt_S_binning = np.zeros(self.histogram_bin_number * 2)
                for x in self.tt_S_no_gaps:
                    self.tt_S_binning[(x - 1) // Config.frequency_sweep_duration] += 1
                S_det_isnt_full = sum(np.array(self.tt_S_binning[int(0.95*len(self.tt_S_binning)):])) != 0

                print(count)
                count += 1
                time.sleep(0.01)
                if is_new_tts_S and S_det_isnt_full:
                    break
                if self.keyPress == 'ESC':
                    print('\033[94m' + 'ESC pressed. Stopping measurement.' + '\033[0m')  # print blue
                    self.updateValue("QRAM_Exp_switch", False)
                    self.MOT_switch(True)
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

            # fold reflections and transmission
            self.tt_N_binning = np.zeros(self.histogram_bin_number*2)

            for x in self.tt_N_directional_measure:
                self.tt_N_binning[(x - 1) // Config.frequency_sweep_duration] += 1  # TODO: x-1?? in spectrum its x
            self.tt_N_binning_avg = self.tt_N_binning

            self.tt_S_binning_avg = self.tt_S_binning
            # for x in self.tt_S_directional_measure:
            for x in self.tt_S_no_gaps:
                # if x < (2e6 + 6400):
                self.folded_tt_S_acc[(x-1) % self.histogram_bin_size] += 1
                # if (2e6 + 6400 + 120) < x < (4e6 + 6400*2 + 120):
                if 6e6 < x < 8e6:
                    self.folded_tt_S_acc_2[(x - 1) % self.histogram_bin_size] += 1
                # if (8e6 + 6400*4 + 120*4) < x < (10e6):
                if 8e6 < x < 10e6:
                    self.folded_tt_S_acc_3[(x - 1) % self.histogram_bin_size] += 1

            # split the binning vector to odd and even - on and off resonance pulses
            self.tt_N_binning_detuned = [self.tt_N_binning[x] for x in range(len(self.tt_N_binning)) if x % 2]  # odd
            self.tt_N_binning_resonance = [self.tt_N_binning[x] for x in range(len(self.tt_N_binning)) if
                                           not x % 2]  # even
            self.tt_S_binning_detuned = [self.tt_S_binning[x] for x in range(len(self.tt_S_binning)) if x % 2]  # odd
            self.tt_S_binning_resonance = [self.tt_S_binning[x] for x in range(len(self.tt_S_binning)) if
                                           not x % 2]  # even

            self.sum_for_threshold = (np.sum(self.tt_S_binning_resonance) * 1000) / (self.M_time / 2)
            print(f'num of photons in [us] {self.sum_for_threshold}')

            if exp_flag:
                if (self.lock_err > lock_err_threshold) or \
                        ((1000 * np.average(self.FLR_res.tolist()) < FLR_threshold) and with_atoms) or \
                        self.latched_detectors():
                    self.acquisition_flag = False
                else:
                    self.acquisition_flag = True

            self.threshold_flag = (self.sum_for_threshold < total_counts_threshold)

            if (self.threshold_flag or not exp_flag) and self.acquisition_flag and not self.pause_flag:

                if self.Counter < N:
                    self.Counter += 1

                self.tt_N_binning = np.zeros(self.histogram_bin_number*2)
                self.tt_S_binning = np.zeros(self.histogram_bin_number*2)
                self.tt_S_transit_events = np.zeros(self.histogram_bin_number)

                self.tt_S_binning_avg = self.tt_S_binning_avg * (1 - 1 / self.Counter)

                for x in self.tt_N_directional_measure:
                    self.tt_N_binning[(x - 1) // self.histogram_bin_size] += 1
                for x in self.tt_S_no_gaps:
                    self.tt_S_binning[(x - 1) // self.histogram_bin_size] += 1
                    self.tt_S_binning_avg[(x - 1) // self.histogram_bin_size] += 1 / self.Counter

                self.tt_N_measure_batch = self.tt_N_measure_batch[-(N - 1):] + [self.tt_N_directional_measure]
                self.tt_N_binning_batch = self.tt_N_binning_batch[-(N - 1):] + [self.tt_N_binning]
                self.tt_N_binning_resonance_batch = self.tt_N_binning_resonance_batch[-(N - 1):] + [self.tt_N_binning_resonance]
                self.tt_N_binning_detuned_batch = self.tt_N_binning_detuned_batch[-(N - 1):] + [self.tt_N_binning_detuned]
                self.tt_S_measure_batch = self.tt_S_measure_batch[-(N - 1):] + [self.tt_S_no_gaps]
                self.tt_S_binning_batch = self.tt_S_binning_batch[-(N - 1):] + [self.tt_S_binning]
                self.tt_S_binning_resonance_batch = self.tt_S_binning_resonance_batch[-(N - 1):] + [self.tt_S_binning_resonance]
                self.tt_S_binning_detuned_batch = self.tt_S_binning_detuned_batch[-(N - 1):] + [self.tt_S_binning_detuned]
                self.S_bins_res_acc = np.sum(np.array(self.tt_S_binning_resonance_batch), 0) # tt_S_binning_resonance accumulated (sum over the batch)
                self.S_bins_detuned_acc = np.sum(np.array(self.tt_S_binning_detuned_batch),
                                                 0)  # tt_S_binning_resonance accumulated (sum over the batch)

                # self.find_transit_events_spectrum(N, transit_time_threshold=time_threshold,
                #                                   transit_counts_threshold=transit_counts_threshold)
                self.find_transits_events_spectrum_exp(self.tt_S_binning_resonance, N, transit_cond)

                FLR_measurement = FLR_measurement[-(N - 1):] + [self.FLR_res.tolist()]
                Exp_timestr_batch = Exp_timestr_batch[-(N - 1):] + [timest]
                lock_err_batch = lock_err_batch[-(N - 1):] + [self.lock_err]
                self.save_tt_to_batch(Num_Of_dets, N)

        ############################################## END WHILE LOOP #################################################

        # For debuging:
        # self.most_common([x for vec in self.tt_S_measure_batch + self.tt_N_measure_batch for x in vec])

        ## Adding comment to measurement [prompt whether stopped or finished regularly]
        if exp_flag:
            if self.Counter < N:
                aftComment = pymsgbox.prompt('Add comment to measurement: ', default='', timeout=int(30e3))
            else:
                aftComment = ''
        else:
            aftComment = 'ignore'

        # if aftComment == 'Timeout': aftComment = None

        #### Handle file-names and directories #####
        ## Saving: np.savez(filedir, data = x) #note: @filedir is string of full directory; data is the queyword used to read @x from the file:
        ## Loading: file = np.load(f, allow_pickle = True)
        ##          x = file['data']

        #### ------ Save results ------
        #  -------   Create dir
        root_dirname = f'U:\\Lab_2023\\Experiment_results\\Spectrum\\{datest}\\'
        dirname = root_dirname + f'{timest}_Photon_TimeTags\\'  # Specific experiment dir
        dirname_Det = dirname + 'AllDetectors\\'
        dirname_N = dirname + 'North(8)\\'
        dirname_S = dirname + 'South(5)\\'
        dirname_D = dirname + 'Dark(3,4)\\'
        dirname_B = dirname + 'Bright(1,2)\\'
        dirname_FS = dirname + 'FastSwitch(6,7)\\'
        dirname_expfigures = root_dirname + 'Experiment_figures\\'
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
        if not os.path.exists(dirname_expfigures):
            os.makedirs(dirname_expfigures)

        # ----  msmnt files names  -----
        # Counter_str = (self.Counter)
        data_to_save = ["frequency_sweep_rep", "spectrum_bandwidth", "frequency_start", "frequency_diff", 
                        "num_of_different_frequencies", "same_frequency_rep", "freq_bins", "sequence_duration"]
        filename_Det_tt = []
        for i in Num_Of_dets:
            filename_Det_tt.append(f'Det' + str(i) + f'_timetags.npz')
        filename_S_tt = f'South_timetags.npz'
        filename_N_tt = f'North_timetags.npz'
        filename_DP_tt = f'Dark_timetags.npz'
        filename_BP_tt = f'Bright_timetags.npz'
        filename_FS_tt = f'FS_timetags.npz'
        filename_FLR = f'Flouresence.npz'
        filename_timestamp = f'Drops_time_stamps.npz'
        filename_lockerr = f'Cavity_Lock_Error.npz'
        filname_sequence_S = f'South_sequence_vector.npz'
        filname_sequence_N = f'North_sequence_vector.npz'
        filname_sequence_Early = f'Early_sequence_vector.npz'
        filname_sequence_Late = f'Late_sequence_vector.npz'
        filname_sequence_FS = f'FS_sequence_vector.npz'
        filename_transits_events = f'seq_transit_events_batched.npz'
        filename_transits_batched = f'all_transits_seq_indx_batch.npz'
        filename_Cavity_spectrum = f'Cavity_spectrum.npz'
        filename_Cavity_atom_spectrum = f'Cavity_atom_spectrum.npz'
        filename_freq_vector = f'frequency_vector.npz'
        filename_experimentPlot = f'Experiment_plot.png'
        filename_experimentfigure = f'{timest}_Experiment_Figure.png'
        filaname_spectrum_parameters = f'Spectrum_parameters.mat'

        if len(FLR_measurement) > 0:
            np.savez(dirname + filename_FLR, FLR_measurement)
        if len(Exp_timestr_batch) > 0:
            np.savez(dirname + filename_timestamp, Exp_timestr_batch)
            np.savez(dirname + filename_lockerr, lock_err_batch)
            np.savez(dirname + filname_sequence_S, Config.QRAM_Exp_Gaussian_samples_S)
            np.savez(dirname + filname_sequence_N, Config.QRAM_Exp_Gaussian_samples_N)
            np.savez(dirname + filname_sequence_Early, Config.QRAM_Exp_Square_samples_Early)
            np.savez(dirname + filname_sequence_Late, Config.QRAM_Exp_Square_samples_Late)
            np.savez(dirname + filname_sequence_FS, (1 - np.array(Config.QRAM_Exp_Square_samples_FS)).tolist())
            plt.savefig(dirname + filename_experimentPlot, bbox_inches='tight')
            plt.savefig(dirname_expfigures + filename_experimentfigure, bbox_inches='tight')
            self.save_matlab_file(data_to_save, "spectrum_parameters", dirname + filaname_spectrum_parameters)
        for i in range(len(Num_Of_dets)):
            if len(self.tt_measure_batch[i]) > 0:
                np.savez(dirname_Det + filename_Det_tt[i], self.tt_measure_batch[i])
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
        if len(self.all_transits_indx_batch) > 0:
            np.savez(dirname + filename_transits_events, self.all_transits_indx_batch)
        if len(self.Cavity_spectrum) > 0:
            np.savez(dirname + filename_Cavity_spectrum, self.Cavity_spectrum)
            np.savez(dirname + filename_freq_vector, self.freq_bins)
        if len(self.Cavity_atom_spectrum) > 0:
            np.savez(dirname + filename_Cavity_atom_spectrum, self.Cavity_atom_spectrum)
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
        # if aftComment is not None: cmnt = preComment + aftComment
        if aftComment is not None: cmnt = aftComment
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

        # experiment_cmnt = f' transit condition: \n minimum time between reflection = ' + str(time_threshold) + \
        #                   f'\n with at least ' + str(transit_counts_threshold) + ' reflections ' + \
        #                   f'\n reflection threshold: ' + str(total_counts_threshold) + ' [MCounts / sec]'

        experiment_cmnt = f'Transit condition: ' + str(transit_cond) + \
                          f'\nWhich means - for ' + str(len(transit_cond)) + \
                          ' consecutive on-resonance pulses at least 2 of the conditions must apply.' + \
                          '\nEach element represents the minimum number of photon required per pulse to be regarded as transit.'


        # comments = {'comments': experiment_cmnt}
        try:
            with open(f'{dirname}experiment_comments.txt', 'w') as file:
                # json.dump(comments, file, indent=4)
                file.write(experiment_cmnt)
        except Exception:
            pass

        self.saveConfigTable(path=dirname)
        for qrdCtrl in self.QuadRFControllers:
            qrdCtrl.saveLinesAsCSV(f'{dirname}QuadRF_table.csv')

        self.updateValue("QRAM_Exp_switch", False)
        self.update_parameters()

        ## ------------------ end of saving section -------

    def Start_Transit_Exp_with_tt(self, N=500, Histogram_bin_size=1000, Transit_profile_bin_size=100, preComment=None,
                                  total_counts_threshold=1, transit_counts_threshold=5, FLR_threshold=0.11,
                                  lock_err_threshold=0.003, Exp_flag=True, with_atoms=True):
        self.QRAM_Exp_switch(True)
        self.MOT_switch(with_atoms)
        self.update_parameters()
        self.Save_SNSPDs_Transit_Measurement_with_tt(N, Histogram_bin_size, Transit_profile_bin_size, preComment,
                                                     total_counts_threshold, transit_counts_threshold, FLR_threshold,
                                                     lock_err_threshold, exp_flag=Exp_flag,
                                                     with_atoms=with_atoms)

    def Start_Spectrum_Exp_with_tt(self, N=1000, Transit_profile_bin_size=100, preComment=None, transit_cond=[2, 1, 2],
                                  total_counts_threshold=0.1, transit_counts_threshold=3, FLR_threshold=0.03,
                                  lock_err_threshold=0.002, Exp_flag=True, with_atoms=True, Calibration_dirname=None):
        self.Spectrum_Exp_switch(True)
        self.MOT_switch(with_atoms)
        self.update_parameters()
        self.Save_SNSPDs_Spectrum_Measurement_with_tt(N, Transit_profile_bin_size, preComment, transit_cond,
                                                     total_counts_threshold, transit_counts_threshold, FLR_threshold,
                                                     lock_err_threshold, exp_flag=Exp_flag,
                                                     with_atoms=with_atoms, calibration_dirname=Calibration_dirname)

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

    def on_press(self, key):
        # print(f'{key} pressed')
        if key == keyboard.Key.shift or key == keyboard.Key.esc:
            self.alt_modifier = True

    def on_release(self, key):
        # print(f'{key} released')
        if key == keyboard.Key.esc and self.alt_modifier:
            # print('Alt-ESC released')
            self.keyPress = 'ESC'

        if key == keyboard.Key.space and self.alt_modifier:
            # print('Alt-Space released')
            self.keyPress = 'SPACE'

        if key == keyboard.Key.shift:
            self.alt_modifier = False

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
    
    def save_matlab_file(self, data_to_save, name, path):
        # Ensure there is the ".mat" extension
        if not path.endswith('.mat'):
            path = path + '.mat'
        # Build dictionary with values
        dict = {
                name: {}
                }
        for key in data_to_save:
            dict[name][key] = self.__dict__[key]
        # Save the dictionary to file
        savemat(path, mdict=dict)
        pass
        

    def run_daily_experiment(self, day_experiment, Histogram_bin_size, Transit_profile_bin_size, preComment,
                             total_counts_threshold, transit_counts_threshold, FLR_threshold, lock_err_threshold,
                             Exp_flag=True):
        with_atoms_bool = False
        for i in range(len(day_experiment)):
            if with_atoms_bool:
                Comment = preComment + ' with atoms'
            else:
                Comment = preComment + ' without atoms'
            self.Start_Transit_Exp_with_tt(N=day_experiment[i], Histogram_bin_size=Histogram_bin_size,
                                           Transit_profile_bin_size=Transit_profile_bin_size,
                                           preComment=Comment, total_counts_threshold=total_counts_threshold,
                                           transit_counts_threshold=transit_counts_threshold,
                                           FLR_threshold=FLR_threshold, lock_err_threshold=lock_err_threshold,
                                           Exp_flag=Exp_flag, with_atoms=with_atoms_bool)
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
    # experiment.Start_Spectrum_Exp_with_tt(preComment="Spectrum experiment; bandwidth 48 MHz, 1.5MHz jumps")
    # experiment.Start_Spectrum_Exp_with_tt(N=200, preComment="Spectrum experiment; On resonance max counts",
    #                                       total_counts_threshold=10, with_atoms=False)
    # experiment.Start_Spectrum_Exp_with_tt(Exp_flag=False)



    # experiment.Start_Spectrum_Exp_with_tt(total_counts_threshold=0.3, Exp_flag=False)

    # experiment.run_daily_experiment([50, 500], Histogram_bin_size=1000, Transit_profile_bin_size=10,
    #                                 preComment='Transit experinment with 2uW of 731.30nm light coupled',
    #                                 total_counts_threshold=1, transit_counts_threshold=5, FLR_threshold=0.08,
    #                                 lock_err_threshold=0.004, Exp_flag=False)
    #     experiment.job.halt()
    #     experiment.qmm.reset_data_processing()
    # finally:
    #     experiment.job.halt()
    #     experiment.qmm.reset_data_processing()
    # experiment.updateValue('PrePulse_duration', 5)
    # experiment.update_parameters()
    # experiment.Repeat_Measurement(10)
    # experiment.Repeat_Measurement(10)
