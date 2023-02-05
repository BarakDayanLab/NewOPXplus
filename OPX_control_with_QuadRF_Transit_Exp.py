# from Config import config
import Config_with_SNSPDs_and_QuadRF as Config
from Config_Table import Initial_Values, Phases_Names  # , Values_Factor
from quadRFMOTController import QuadRFMOTController
from quadRFFrequencyScannerController import QuadRFFrequencyScannerController

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
          "AOM_2-2/3'", "AOM_2-3'_for_interference", "AOM_2-2'", "FLR_detection", "Measurement") # , "Dig_detectors_spectrum", "Dig_detectors") # , "PULSER_N", "PULSER_S")

    ## MOT build-up ##
    n = declare(int)
    m = declare(int)
    play("Detection" * amp(FLR * 0), "FLR_detection", duration=4)  # we dont know why this works, but Yoav from QM made us write this line to solve an alignment problem we had in the next 2 for loops
    with for_(n, 1, n <= mot_repetitions, n + 1):
        play("MOT" * amp(Config.AOM_0_Attenuation), "MOT_AOM_0")
        play("MOT" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-")
        play("MOT" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+")
        # play("OD_FS" * amp(0.5), "AOM_2-2/3'")
        play("AntiHelmholtz_MOT", "AntiHelmholtz_Coils")
    with for_(m, 1, m <= (mot_repetitions - 1), m + 1):
        measure("Detection", "FLR_detection", None, integration.full("Detection_opt", FLR, "out1"))
        # play("OD_FS" * amp(0.1), "AOM_2-3'_for_interference")

    align("Cooling_Sequence", "MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+", "AntiHelmholtz_Coils", "Zeeman_Coils",
          "AOM_2-2/3'", "AOM_2-2'", "FLR_detection", "Measurement") # , "Dig_detectors", "Dig_detectors_spectrum") #, "PULSER_N", "PULSER_S")

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

    update_frequency("MOT_AOM_0", Config.IF_AOM_MOT)
    update_frequency("MOT_AOM_-", Config.IF_AOM_MOT)
    update_frequency("MOT_AOM_+", Config.IF_AOM_MOT)

    ## Aligning all the different elements used during the freefall time of the experiment ##
    align("Cooling_Sequence", "MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+", "Zeeman_Coils", "AOM_2-2/3'", "AOM_2-2'",
          "Measurement", "PULSER_N", "PULSER_S") # , "Dig_detectors", "Dig_detectors_spectrum")

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


def Probe_counts_Measure_SNSPDs(m_off_time, m_time, m_window, shutter_open_time,
                                ON_counts_st1, ON_counts_st2, ON_counts_st3, ON_counts_st4,
                                ON_counts_st5, ON_counts_st6, ON_counts_st7, ON_counts_st8 ):
    """
    Measuring the OD with and without atoms.

    Parameters
    ----------
    :param M_delay: The time delay from end of PGC until the first Probe pulse.
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
    counts1 = declare(int)
    counts2 = declare(int)
    counts3 = declare(int)
    counts4 = declare(int)
    counts5 = declare(int)
    counts6 = declare(int)
    counts7 = declare(int)
    counts8 = declare(int)

    n = declare(int)

    align("AOM_2-2/3'", "Dig_detectors")
    # wait(M_delay - shutter_open_time, "Dig_detectors", "AOM_2-2/3'")
    play("OD", "AOM_2-2/3'", duration=shutter_open_time)
    align("AOM_2-2/3'", "Dig_detectors")
    play("OD", "AOM_2-2/3'", duration=m_time)
    with for_(n, 0, n < (m_time * 4), n + m_window):
        measure("readout", "Dig_detectors", None,
                counting.digital(counts1, m_window, element_outputs="out1"),
                counting.digital(counts2, m_window, element_outputs="out2"),
                counting.digital(counts3, m_window, element_outputs="out3"),
                counting.digital(counts4, m_window, element_outputs="out4"),
                counting.digital(counts5, m_window, element_outputs="out5"),
                counting.digital(counts6, m_window, element_outputs="out6"),
                counting.digital(counts7, m_window, element_outputs="out7"),
                counting.digital(counts8, m_window, element_outputs="out8"),
                )
        # wait(m_off_time, "Dig_detectors")

        ## Save Data: ##

        ## Number of Photons (NOP) Count stream for each detector: ##

        save(counts1, ON_counts_st1)
        save(counts2, ON_counts_st2)
        save(counts3, ON_counts_st3)
        save(counts4, ON_counts_st4)
        save(counts5, ON_counts_st5)
        save(counts6, ON_counts_st6)
        save(counts7, ON_counts_st7)
        save(counts8, ON_counts_st8)


def Sprint_Exp(m_off_time, m_time, m_window, shutter_open_time,
               ON_counts_st1, ON_counts_st2, ON_counts_st3, ON_counts_st4,
               ON_counts_st5, ON_counts_st6, ON_counts_st7, ON_counts_st8, ON_counts_st9, ON_counts_st10,
               tt_st_1, tt_st_2, tt_st_3, tt_st_4, tt_st_5, tt_st_6, tt_st_7, tt_st_8,tt_st_9, tt_st_10, rep_st):
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
    counts9 = declare(int)
    counts10 = declare(int)

    tt_vec1 = declare(int, size=vec_size)
    tt_vec2 = declare(int, size=vec_size)
    tt_vec3 = declare(int, size=vec_size)
    tt_vec4 = declare(int, size=vec_size)
    tt_vec5 = declare(int, size=vec_size)
    tt_vec6 = declare(int, size=vec_size)
    tt_vec7 = declare(int, size=vec_size)
    tt_vec8 = declare(int, size=vec_size)
    tt_vec9 = declare(int, size=vec_size)
    tt_vec10 = declare(int, size=vec_size)

    n = declare(int)
    t = declare(int)
    m = declare(int)

    align("PULSER_N", "PULSER_S", "Dig_detectors")
    play("Const_open", "PULSER_S", duration=shutter_open_time)
    # wait(shutter_open_time, "PULSER_S")
    align("PULSER_N", "PULSER_S", "Dig_detectors")

    # play("OD", "AOM_2-2/3'", duration=m_time)  # CRUS with pulses from AWG and constantly ON-resonance

    with for_(t, 0, t < (m_time + m_off_time) * 4, t + int(len(Config.Sprint_Exp_Gaussian_samples_S))): #assaf comment debbuging
        play("Sprint_experiment_pulses_S", "PULSER_S")
        play("Sprint_experiment_pulses_N", "PULSER_N")

    with for_(n, 0, n < m_time * 4, n + m_window):
        measure("readout_SPRINT", "Dig_detectors", None,
                time_tagging.digital(tt_vec1, m_window, element_output="out1", targetLen=counts1),
                time_tagging.digital(tt_vec2, m_window, element_output="out2", targetLen=counts2),
                time_tagging.digital(tt_vec3, m_window, element_output="out3", targetLen=counts3),
                time_tagging.digital(tt_vec4, m_window, element_output="out4", targetLen=counts4),
                time_tagging.digital(tt_vec5, m_window, element_output="out5", targetLen=counts5),
                time_tagging.digital(tt_vec6, m_window, element_output="out6", targetLen=counts6),
                time_tagging.digital(tt_vec7, m_window, element_output="out7", targetLen=counts7),
                time_tagging.digital(tt_vec8, m_window, element_output="out8", targetLen=counts8),
                time_tagging.digital(tt_vec9, m_window, element_output="out8", targetLen=counts9),
                time_tagging.digital(tt_vec10, m_window, element_output="out8", targetLen=counts10),
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
        save(counts9, ON_counts_st9)
        save(counts10, ON_counts_st10)

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
            save(tt_vec9[m], tt_st_9)
            save(tt_vec10[m], tt_st_10)
            save(n, rep_st)


def opx_control(obj, qm):
    with program() as opx_control_prog:
        ## declaring program variables: ##
        i = declare(int)
        Trigger_Phase = declare(int, value=Phases_Names.index(obj.Exp_Values['Triggering_Phase']))
        Imaging_Phase = declare(int, value=Phases_Names.index(obj.Exp_Values['Imaging_Phase']))

        # Boolean variables:
        AntiHelmholtz_ON = declare(bool, value=True)
        SPRINT_Exp_ON = declare(bool, value=False)
        Probe_max_counts_Exp_ON = declare(bool, value=False)

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
        ON_counts_st1 = declare_stream()
        ON_counts_st2 = declare_stream()
        ON_counts_st3 = declare_stream()
        ON_counts_st4 = declare_stream()
        ON_counts_st5 = declare_stream()
        ON_counts_st6 = declare_stream()
        ON_counts_st7 = declare_stream()
        ON_counts_st8 = declare_stream()
        ON_counts_st9 = declare_stream()
        ON_counts_st10 = declare_stream()
        tt_st_1 = declare_stream()
        tt_st_2 = declare_stream()
        tt_st_3 = declare_stream()
        tt_st_4 = declare_stream()
        tt_st_5 = declare_stream()
        tt_st_6 = declare_stream()
        tt_st_7 = declare_stream()
        tt_st_8 = declare_stream()
        tt_st_9 = declare_stream()
        tt_st_10 = declare_stream()
        rep_st = declare_stream()
        AntiHelmholtz_ON_st = declare_stream()
        FLR_st = declare_stream()
        x = declare(int)
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
            with if_(SPRINT_Exp_ON):
                assign(x, (23000000+656000 * 2) // 4) # TODO -  added 23000000 to fix new delay due to wait(1000) in saving sprint data, should be fixed as well
            with else_():
                assign(x, 0)
            FreeFall(FreeFall_duration - x, coils_timing)
            ##########################
            ## Measurement Sequence ##
            ##########################

            with if_(Trigger_Phase == 3):  # when trigger on PrePulse
                ## Trigger QuadRF Sequence #####################
                play("C_Seq", "Cooling_Sequence", duration=2500)
                ################################################

            with if_(SPRINT_Exp_ON):
                play("Depump", "AOM_2-2'", duration=(PrePulse_duration - shutter_open_time))
            with else_():
                wait(PrePulse_duration, "Cooling_Sequence")
            align(*all_elements, "AOM_2-2/3'", "AOM_2-2'", "PULSER_N", "PULSER_S" , "Dig_detectors")

            with if_(Trigger_Phase == 4):  # when trigger on pulse 1
                ## Trigger QuadRF Sequence #####################
                play("C_Seq", "Cooling_Sequence", duration=2500)
                ################################################
            with if_((Imaging_Phase == 4) & (Pulse_1_duration > 0)):  # 4 means imaging phase on pulse_1
                with if_(SPRINT_Exp_ON):
                    align("Dig_detectors", "PULSER_N", "PULSER_S")
                    Sprint_Exp(M_off_time, Pulse_1_duration, obj.M_window, shutter_open_time,
                               ON_counts_st1, ON_counts_st2, ON_counts_st3, ON_counts_st4,
                               ON_counts_st5, ON_counts_st6, ON_counts_st7,ON_counts_st8, ON_counts_st9, ON_counts_st10,
                               tt_st_1, tt_st_2, tt_st_3, tt_st_4, tt_st_5, tt_st_6, tt_st_7, tt_st_8, tt_st_9, tt_st_10, rep_st)
                    save(AntiHelmholtz_ON, AntiHelmholtz_ON_st)
                    with if_(AntiHelmholtz_ON):
                        save(FLR, FLR_st)
                    align("Dig_detectors", "PULSER_N", "PULSER_S")
                with else_():
                    with if_(Probe_max_counts_Exp_ON):
                        align("Dig_detectors", "AOM_2-2/3'")
                        Probe_counts_Measure_SNSPDs(M_off_time, Pulse_1_duration, obj.M_window,
                                                    shutter_open_time, ON_counts_st1, ON_counts_st2,
                                                    ON_counts_st3, ON_counts_st4, ON_counts_st5,
                                                    ON_counts_st6, ON_counts_st7, ON_counts_st8)
                        align("Dig_detectors", "AOM_2-2/3'")
                    ## For taking an image:
                    with else_():
                        align(*all_elements)
                        Pulse_with_prep(Pulse_1_duration, Pulse_1_decay_time, pulse_1_duration_0,
                                        pulse_1_duration_minus, pulse_1_duration_plus)
                        Measure(Pulse_1_duration)  # This triggers camera (Control 7)
                        align(*all_elements)

            assign(N_Snaps, 1)
            assign(Buffer_Cycles, 0)
            assign(i, IO1)

            ## PARAMETERS UPDATE ##
            with if_(i > 0):
                pause()
            with while_(i > 0):
                ## Boolean variables control: ##
                with if_(i == 4):
                    assign(SPRINT_Exp_ON, IO2)
                with if_(i == 6):
                    assign(Probe_max_counts_Exp_ON, IO2)
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

                ## Fountain variables control: ##
                with if_(i == 31):  # Live control over the fountain duration
                    assign(fountain_duration, IO2)
                with if_(i == 36):  # Live control over the final amplitude of the fountain AOM 0
                    assign(fountain_pulse_duration_0, IO2)
                with if_(i == 37):  # Live control over the final amplitude of the fountain AOM -
                    assign(fountain_pulse_duration_minus, IO2)
                with if_(i == 38):  # Live control over the final amplitude of the fountain AOM +
                    assign(fountain_pulse_duration_plus, IO2)
                with if_(i == 39):  # Live control over the fountain final frequency of the MOT + and - AOMs
                    assign(fountain_aom_chirp_rate, IO2)

                ## Measurement variables control: ##
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
                with if_(i == 56):  # Live control of the Depump measurement start time
                    assign(Depump_start, IO2)
                with if_(i == 57):  # Live control of the Depump measurement pulses duration
                    assign(Depump_pulse_duration, IO2)
                with if_(i == 58):  # Live control of the Depump measurement wait duration between the 2 pulses
                    assign(Depump_pulses_spacing, IO2)
                with if_(i == 59):  # Live control of the delay due to shutter opening time.
                    assign(shutter_open_time, IO2)


                pause()
                assign(i, IO1)

        with stream_processing():
            # (ON_counts_st1 + ON_counts_st2 + ON_counts_st3 + ON_counts_st4).buffer(obj.rep).save('North_Probe')
            # (ON_counts_st1 + ON_counts_st2 ).buffer(obj.rep).save('North_Probe')
            # (ON_counts_st5 + ON_counts_st6 + ON_counts_st7 + ON_counts_st8).buffer(obj.rep).save('South_Probe')
            # (ON_counts_st5 + ON_counts_st6+ ON_counts_st7 + ON_counts_st8).buffer(obj.rep).save('South_Probe')
            # ON_counts_st1.buffer(obj.rep).save('Det1_Counts')
            # ON_counts_st2.buffer(obj.rep).save('Det2_Counts')
            ON_counts_st3.buffer(obj.rep).save('Det6_Counts')
            ON_counts_st4.buffer(obj.rep).save('Det5_Counts')
            ON_counts_st5.buffer(obj.rep).save('Det1_Counts') # assaf - renamed to det3 to ease further calculation
            ON_counts_st6.buffer(obj.rep).save('Det2_Counts') # assaf - renamed to det4 to ease further calculation
            ON_counts_st7.buffer(obj.rep).save('Det3_Counts')
            ON_counts_st8.buffer(obj.rep).save('Det4_Counts')
            # (tt_st_1 + rep_st).buffer(obj.vec_size * obj.rep).save('Det1_Probe_TT')
            # (tt_st_2 + rep_st).buffer(obj.vec_size * obj.rep).save('Det2_Probe_TT')
            (tt_st_3 + rep_st).buffer(obj.vec_size * obj.rep).save('Det6_Probe_TT')
            (tt_st_4 + rep_st).buffer(obj.vec_size * obj.rep).save('Det5_Probe_TT')
            (tt_st_5 + rep_st).buffer(obj.vec_size * obj.rep).save('Det1_Probe_TT') # assaf - renamed to det3 to ease further calculation
            (tt_st_6 + rep_st).buffer(obj.vec_size * obj.rep).save('Det2_Probe_TT') # assaf - renamed to det4 to ease further calculation
            (tt_st_7 + rep_st).buffer(obj.vec_size * obj.rep).save('Det3_Probe_TT')
            (tt_st_8 + rep_st).buffer(obj.vec_size * obj.rep).save('Det4_Probe_TT')
            FLR_st.save('FLR_measure')
            AntiHelmholtz_ON_st.save("antihelmholtz_on")

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
        qrfContr = QuadRFMOTController(initialValues=self.Exp_Values, updateChannels=(1, 4), topticaLockWhenUpdating=False,
                                        debugging=False, continuous=False)
        self.QuadRFControllers.append(qrfContr)  # updates values on QuadRF (uploads table)
        self.QuadRFControllers.append(QuadRFMOTController(MOGdevice=qrfContr.dev, initialValues={'Operation_Mode': 'Continuous', 'CH3_freq': '100MHz', 'CH3_amp': '31dbm'},
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

        # MW spectroscopy parameters:
        self.MW_start_frequency = int(100e6)  # [Hz]
        self.Pulse_Length_MW = 400  # [usec]
        self.Pulse_Length_OD = 20  # [usec]

        # Main Experiment:
        self.TOP2_pulse_len = int(Config.Probe_pulse_len / 4)  # [nsec]
        self.Calibration_time = 10  # [msec]

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

    def MOT_switch(self, Bool):
        self.update_io_parameter(1, Bool)

    def Linear_PGC_switch(self, Bool):
        self.update_io_parameter(2, Bool)

    def Transit_Exp_switch(self, Bool):
        self.Transit_switch = Bool
        self.update_io_parameter(3, Bool)

    def CRUS_Exp_switch(self, Bool):
        self.update_io_parameter(4, Bool)

    def Spectrum_Exp_switch(self, Bool):
        self.update_io_parameter(4, Bool)

    def SPRINT_Exp_switch(self, Bool):
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

    # def Save_SNSPDs_Transit_Measurement_with_tt(self, N, histogram_bin_size, Transit_profile_bin_size, preComment,
    #                                             total_counts_threshold,  transit_counts_threshold, max_probe_counts):
    #     """
    #     Function for analyzing and saving the time tags data measured from the SNSPDs using the OPX. In this specific
    #      program we are looking for transits of atoms next to the toroid and record them.
    #     :param N: Number of maximum experiments (free throws) saved and displayed.
    #     :param histogram_bin_size: The bin size for the general experiment histogram which means - dividing the length
    #                                of the measuring time (m_time) to bins and counting the number of photon detections
    #                                at each bin.
    #     :param Transit_profile_bin_size: The bin size for the transit histogram which means - dividing the length of the
    #                                      transit time to bins and counting the number of transits for which there has
    #                                      been a detection of photon at each bin.
    #     :param preComment: The comment added at the start of the experiment, usually consisting of unique experiment
    #                        parameters.
    #     :param total_counts_threshold: The minimum number of MCounts / sec for which the data will be displayed and collected.
    #     :param intensity_threshold: The number of photons at each time bin for which suspect we as a transit.
    #     :param max_probe_counts: The maximum counts probe counts at each direction measured when cavity isn't locked.
    #     :return:
    #     """
    #     # if preComment is True:
    #     if not preComment:
    #         preComment = pymsgbox.prompt('Add comment to measurement: ', default='', timeout=int(30e3))
    #     aftComment = None
    #
    #     ### fetching data from server
    #     ### saving to file
    #     ###
    #
    #     histogram_bin_number = self.M_time // histogram_bin_size
    #     time_bins = np.linspace(0, self.M_time, histogram_bin_number)
    #     # time_threshold = int(histogram_bin_size / intensity_threshold)  # The minimum time between two time tags to be counted for a transit. # TODO: might need a factor of 2???
    #     time_threshold = int(
    #         histogram_bin_size * 0.8)  # The minimum time between two time tags to be counted for a transit. # TODO: might need a factor of 2???
    #
    #     ## Listen for keyboard
    #     listener = keyboard.Listener(on_press=self.on_key_press)
    #     listener.start()  # start to listen on a separate thread
    #     self.keyPress = None
    #     print('\033[94m' + 'Press ESC to stop measurement.' + '\033[0m')  # print blue
    #     reps = 1  # a counter, number of repeats actually made.
    #
    #     Probe_N_handle = self.job.result_handles.get("North_Probe")
    #     Probe_S_handle = self.job.result_handles.get("South_Probe")
    #     tt_N_handle = self.job.result_handles.get("North_Probe_TT")
    #     tt_S_handle = self.job.result_handles.get("South_Probe_TT")
    #     FLR_handle = self.job.result_handles.get("FLR_measure")
    #
    #     Probe_N_handle.wait_for_values(1)
    #     Probe_S_handle.wait_for_values(1)
    #     tt_N_handle.wait_for_values(1)
    #     tt_S_handle.wait_for_values(1)
    #     FLR_handle.wait_for_values(1)
    #
    #     Probe_N_res = Probe_N_handle.fetch_all()
    #     Probe_S_res = Probe_S_handle.fetch_all()
    #     tt_N_res = tt_N_handle.fetch_all()
    #     tt_S_res = tt_S_handle.fetch_all()
    #     FLR_res = -FLR_handle.fetch_all()
    #
    #     tt_N_measure = [i for i in tt_N_res if (i % self.M_window) != 0]
    #     tt_S_measure = [i for i in tt_S_res if (i % self.M_window) != 0]
    #     tt_N_measure.sort()
    #     tt_S_measure.sort()
    #
    #     tt_N_measure_batch = []
    #     tt_N_binning_batch = []
    #     self.tt_S_measure_batch = []
    #     tt_S_binning_batch = []
    #     all_transits_batch = []
    #     all_transits_aligned_first_batch = []
    #     transit_histogram_batch = []
    #     FLR_measurement = []
    #     Exp_timestr_batch = []
    #
    #     tt_N_binning = np.zeros(histogram_bin_number)
    #     tt_S_binning = np.zeros(histogram_bin_number)
    #     tt_N_transit_events = np.zeros(histogram_bin_number)
    #     tt_S_transit_events = np.zeros(histogram_bin_number)
    #
    #     start = True
    #
    #     # Place holders for results
    #     # while (number of photons * 10^-6 [Mcounts] / Measuring time [nsec] * 10^-9 [sec/nsec])  > total_counts_threshold [Mcounts /sec]:
    #     while ((len(tt_S_measure) * 1000) / self.M_time) > total_counts_threshold or start:
    #         if self.keyPress == 'ESC':
    #             print('\033[94m' + 'ESC pressed. Stopping measurement.' + '\033[0m')  # print blue
    #             self.updateValue("Transit_Exp_switch", False)
    #             self.update_parameters()
    #             # Other actions can be added here
    #             break
    #         if start:
    #             start = False
    #         else:
    #             print('Above Threshold')
    #
    #         Probe_N_res = Probe_N_handle.fetch_all()
    #         Probe_S_res = Probe_S_handle.fetch_all()
    #         tt_N_res = tt_N_handle.fetch_all()
    #         tt_S_res = tt_S_handle.fetch_all()
    #         FLR_res = -FLR_handle.fetch_all()
    #
    #         tt_N_measure = [i for i in tt_N_res if (i % self.M_window) != 0]
    #         tt_S_measure = [i for i in tt_S_res if (i % self.M_window) != 0]
    #         tt_N_measure.sort()
    #         tt_S_measure.sort()
    #
    #     self.tt_S_measure = tt_S_measure
    #     ## record time
    #     timest = time.strftime("%Y%m%d-%H%M%S")
    #     datest = time.strftime("%Y%m%d")
    #
    #     FLR_measurement = FLR_measurement[-(N - 1):] + [FLR_res.tolist()]
    #     Exp_timestr_batch = Exp_timestr_batch[-(N - 1):] + [timest]
    #
    #     for x in tt_N_measure:
    #         tt_N_binning[x // histogram_bin_size] += 1
    #     for x in tt_S_measure:
    #         tt_S_binning[x // histogram_bin_size] += 1
    #
    #     if len(self.tt_S_measure_batch) == N:
    #         tt_N_transit_events[[i for i, x in enumerate(tt_N_binning_batch[0]) if x >  transit_counts_threshold]] -= 1
    #         tt_S_transit_events[[i for i, x in enumerate(tt_S_binning_batch[0]) if x >  transit_counts_threshold]] -= 1
    #
    #     tt_N_measure_batch = tt_N_measure_batch[-(N - 1):] + [tt_N_measure]
    #     tt_N_binning_batch = tt_N_binning_batch[-(N - 1):] + [tt_N_binning]
    #     self.tt_S_measure_batch = self.tt_S_measure_batch[-(N - 1):] + [tt_S_measure]
    #     tt_S_binning_batch = tt_S_binning_batch[-(N - 1):] + [tt_S_binning]
    #     Counter = 1
    #
    #     tt_N_transit_events[[i for i, x in enumerate(tt_N_binning) if x > transit_counts_threshold]] += 1
    #     tt_S_transit_events[[i for i, x in enumerate(tt_S_binning) if x > transit_counts_threshold]] += 1
    #
    #     # Find transits and build histogram:
    #     current_transit = []
    #     all_transits = []
    #     all_transits_aligned_first = []
    #     t_transit = []
    #     t_transit_batch = []
    #     transit_histogram = []
    #     for t in tt_S_measure:
    #         if not current_transit:  # if the array is empty
    #             current_transit.append(t)
    #         elif (t - current_transit[-1]) < time_threshold:
    #             current_transit.append(t)
    #         elif len(current_transit) > transit_counts_threshold:
    #             all_transits.append(current_transit)
    #             all_transits_aligned_first.append([x - current_transit[0] for x in current_transit])
    #             current_transit = [t]
    #         else:
    #             current_transit = [t]
    #
    #     if all_transits:
    #         if len(all_transits_aligned_first_batch) == N:
    #             for n in range(len(all_transits_aligned_first_batch[0])):
    #                 for m in range(len(all_transits_aligned_first_batch[0][n]) - 1):
    #                     transit_histogram_batch[
    #                         all_transits_aligned_first_batch[0][n][m + 1] // Transit_profile_bin_size] -= 1
    #
    #         all_transits_batch = all_transits_batch[-(N - 1):] + [all_transits]
    #         all_transits_aligned_first_batch = all_transits_aligned_first_batch[-(N - 1):] + [
    #             all_transits_aligned_first]
    #         transit_histogram = np.zeros((max([vec for elem in all_transits_aligned_first for vec in elem])
    #                                       + histogram_bin_size) // Transit_profile_bin_size)
    #         t_transit = np.linspace(0, len(transit_histogram) * Transit_profile_bin_size, len(transit_histogram))
    #         transit_histogram_batch = np.zeros((max(
    #             [vec for elem in [vec for elem in all_transits_aligned_first_batch for vec in elem] for vec in
    #              elem]) + histogram_bin_size) // Transit_profile_bin_size)
    #         t_transit_batch = np.linspace(0, len(transit_histogram_batch) * Transit_profile_bin_size,
    #                                       len(transit_histogram_batch))
    #         for n in range(len(all_transits_aligned_first)):
    #             for m in range(len(all_transits_aligned_first[n]) - 1):
    #                 transit_histogram[all_transits_aligned_first[n][m + 1] // Transit_profile_bin_size] += 1
    #                 transit_histogram_batch[all_transits_aligned_first[n][m + 1] // Transit_profile_bin_size] += 1
    #
    #     ## Prepare plots
    #     fig = plt.figure()
    #     ax1 = plt.subplot2grid((6, 2), (0, 0), colspan=1, rowspan=2)
    #     ax2 = plt.subplot2grid((6, 2), (0, 1), colspan=1, rowspan=2)
    #     ax3 = plt.subplot2grid((6, 2), (2, 0), colspan=1, rowspan=2)
    #     ax4 = plt.subplot2grid((6, 2), (2, 1), colspan=1, rowspan=2)
    #     ax5 = plt.subplot2grid((6, 2), (4, 0), colspan=1, rowspan=2)
    #     ax6 = plt.subplot2grid((6, 2), (4, 1), colspan=1, rowspan=2)
    #
    #     while True:
    #         if self.keyPress == 'ESC':
    #             print('\033[94m' + 'ESC pressed. Stopping measurement.' + '\033[0m')  # print blue
    #             self.updateValue("Transit_Exp_switch", False)
    #             self.update_parameters()
    #             # Other actions can be added here
    #             break
    #         if reps < N:
    #             reps += 1
    #
    #         ########################## PLOT!!! ########################################################################
    #
    #         ax1.clear()
    #         ax2.clear()
    #         ax3.clear()
    #         ax4.clear()
    #         ax5.clear()
    #         ax6.clear()
    #
    #         props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #         textstr_N = r'$Probe_N = %.3f$' % ((len(tt_N_measure) * 1000) / self.M_time,) + '[MPhotons/sec]\n' \
    #                     + '$\overline{Probe}_N = %.3f$' % ((np.mean([len(x) for x in tt_N_measure_batch]) * 1000)
    #                                                        / self.M_time,) + '[MPhotons/sec]'
    #         textstr_S = r'$Probe_S = %.3f$' % ((len(tt_S_measure) * 1000) / self.M_time,) + '[MPhotons/sec]\n' \
    #                     + '$\overline{Probe}_S = %.3f$' % ((np.mean([len(x) for x in self.tt_S_measure_batch]) * 1000)
    #                                                        / self.M_time,) + '[MPhotons/sec]'
    #         textstr_FLR = r'$\overline{FLR}_{MAX} = %.1f$' % (np.mean(FLR_measurement) * 1e5,) + r'$\times 10^{-5}$'
    #         textstr_No_transits = 'NO TRANSITS YET!!!'
    #
    #         ax1.plot(time_bins, tt_N_binning, label='Counts histogram', color='b')
    #         ax1.set_title('North', fontweight="bold")
    #         ax1.set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
    #         ax1.text(0.05, 0.95, textstr_N, transform=ax1.transAxes, fontsize=12,
    #                  verticalalignment='top', bbox=props)
    #         ax1.legend(loc='upper right')
    #
    #         ax2.plot(time_bins, tt_S_binning, label='Counts histogram', color='b')
    #         ax2.set_title('South', fontweight="bold")
    #         ax2.set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
    #         ax2.text(0.05, 0.95, textstr_S, transform=ax2.transAxes, fontsize=12,
    #                  verticalalignment='top', bbox=props)
    #         ax2.legend(loc='upper right')
    #
    #         ax3.plot(time_bins, tt_N_transit_events, label='Transit events histogram', marker='*', color='k')
    #         ax3.set(xlabel='Time [msec]', ylabel='Transits events [#Number]')
    #         ax3.legend(loc='upper right')
    #
    #         ax4.plot(time_bins, tt_S_transit_events, label='Transit events histogram', marker='*', color='k')
    #         ax4.set(xlabel='Time [msec]', ylabel='Transits events [#Number]')
    #         ax4.legend(loc='upper right')
    #
    #         if len(transit_histogram) > 0:
    #             textstr_transit_counts = r'$N_{Transits} = %s $' % (len(all_transits_aligned_first),) + r'$[Counts]$'
    #             textstr_avg_transit_counts = r'$\overline{N}_{Transits} = %.1f $' % (
    #             np.average([len(vec) for vec in all_transits_aligned_first]),) + r'$[Counts]$'
    #
    #             ax5.plot(t_transit, transit_histogram, color='b')
    #             ax5.set_title('Drop transits profile', fontweight="bold")
    #             ax5.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
    #             ax5.text(0.05, 0.95, textstr_transit_counts, transform=ax5.transAxes, fontsize=12,
    #                      verticalalignment='top', bbox=props)
    #
    #             ax6.plot(t_transit_batch, transit_histogram_batch, color='b')
    #             ax6.set_title('Accumulated drop transits profile', fontweight="bold")
    #             ax6.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
    #             ax6.text(0.05, 0.95, textstr_avg_transit_counts, transform=ax6.transAxes, fontsize=12,
    #                      verticalalignment='top', bbox=props)
    #         else:
    #             ax5.plot(t_transit, transit_histogram, color='b')
    #             ax5.set_title('Drop transits profile', fontweight="bold")
    #             ax5.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
    #             ax5.text(0.25, 0.5, textstr_No_transits, transform=ax5.transAxes, fontsize=24,
    #                      verticalalignment='center', bbox=props)
    #
    #             ax6.plot(t_transit_batch, transit_histogram_batch, color='b')
    #             ax6.set_title('Accumulated drop transits profile', fontweight="bold")
    #             ax6.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
    #             ax6.text(0.25, 0.5, textstr_No_transits, transform=ax6.transAxes, fontsize=24,
    #                      verticalalignment='center', bbox=props)
    #
    #         ax1.set_ylim(0, 8)
    #         ax2.set_ylim(0, 8)
    #
    #         # plt.tight_layout()
    #         plt.show()
    #         plt.pause(1)
    #
    #         ###########################################################################################################
    #
    #         while tt_S_measure == self.tt_S_measure_batch[-1]:
    #             # record time:
    #             timest = time.strftime("%Y%m%d-%H%M%S")
    #             datest = time.strftime("%Y%m%d")
    #
    #             # get measures:
    #             Probe_N_res = Probe_N_handle.fetch_all()
    #             Probe_S_res = Probe_S_handle.fetch_all()
    #             tt_N_res = tt_N_handle.fetch_all()
    #             tt_S_res = tt_S_handle.fetch_all()
    #             FLR_res = -FLR_handle.fetch_all()
    #
    #             tt_N_measure = [i for i in tt_N_res if (i % self.M_window) != 0]
    #             tt_S_measure = [i for i in tt_S_res if (i % self.M_window) != 0]
    #             tt_N_measure.sort()
    #             tt_S_measure.sort()
    #
    #         if ((len(tt_S_measure) * 1000) / self.M_time) < total_counts_threshold:
    #
    #             FLR_measurement = FLR_measurement[-(N - 1):] + [FLR_res.tolist()]
    #             Exp_timestr_batch = Exp_timestr_batch[-(N - 1):] + [timest]
    #             Counter += 1
    #             print(timest, Counter)
    #
    #             tt_N_binning = np.zeros(histogram_bin_number)
    #             tt_S_binning = np.zeros(histogram_bin_number)
    #
    #             for x in tt_N_measure:
    #                 tt_N_binning[x // histogram_bin_size] += 1
    #             for x in tt_S_measure:
    #                 tt_S_binning[x // histogram_bin_size] += 1
    #
    #             if len(self.tt_S_measure_batch) == N:
    #                 tt_N_transit_events[
    #                     [i for i, x in enumerate(tt_N_binning_batch[0]) if x >  transit_counts_threshold]] -= 1
    #                 tt_S_transit_events[
    #                     [i for i, x in enumerate(tt_S_binning_batch[0]) if x >  transit_counts_threshold]] -= 1
    #
    #             tt_N_measure_batch = tt_N_measure_batch[-(N - 1):] + [tt_N_measure]
    #             tt_N_binning_batch = tt_N_binning_batch[-(N - 1):] + [tt_N_binning]
    #             self.tt_S_measure_batch = self.tt_S_measure_batch[-(N - 1):] + [tt_S_measure]
    #             tt_S_binning_batch = tt_S_binning_batch[-(N - 1):] + [tt_S_binning]
    #
    #             tt_N_transit_events[[i for i, x in enumerate(tt_N_binning) if x >  transit_counts_threshold]] += 1
    #             tt_S_transit_events[[i for i, x in enumerate(tt_S_binning) if x >  transit_counts_threshold]] += 1
    #
    #             # Find transits and build histogram:
    #             current_transit = []
    #             all_transits = []
    #             all_transits_aligned_first = []
    #             for t in tt_S_measure:
    #                 if not current_transit:  # if the array is empty
    #                     current_transit.append(t)
    #                 elif (t - current_transit[-1]) < time_threshold:
    #                     current_transit.append(t)
    #                 elif len(current_transit) >  transit_counts_threshold:
    #                     all_transits.append(current_transit)
    #                     all_transits_aligned_first.append([x - current_transit[0] for x in current_transit])
    #                     current_transit = [t]
    #                 else:
    #                     current_transit = [t]
    #
    #             if all_transits:
    #                 if len(all_transits_aligned_first_batch) == N:
    #                     for n in range(len(all_transits_aligned_first_batch[0])):
    #                         for m in range(len(all_transits_aligned_first_batch[0][n]) - 1):
    #                             transit_histogram_batch[
    #                                 all_transits_aligned_first_batch[0][n][m + 1] // Transit_profile_bin_size] -= 1
    #
    #                 all_transits_batch = all_transits_batch[-(N - 1):] + [all_transits]
    #                 all_transits_aligned_first_batch = all_transits_aligned_first_batch[-(N - 1):] + [
    #                     all_transits_aligned_first]
    #
    #                 transit_histogram = np.zeros((max([vec for elem in all_transits_aligned_first for vec in elem])
    #                                               + histogram_bin_size) // Transit_profile_bin_size)
    #                 t_transit = np.linspace(0, len(transit_histogram) * Transit_profile_bin_size,
    #                                         len(transit_histogram))
    #                 # transit_histogram_batch = np.zeros((max(
    #                 #     [vec for elem in [vec for elem in all_transits_aligned_first_batch for vec in elem] for vec in
    #                 #      elem]) + histogram_bin_size) // Transit_profile_bin_size)
    #                 transit_histogram_batch = np.pad(transit_histogram_batch, (0, (max(
    #                     [vec for elem in [vec for elem in all_transits_aligned_first_batch for vec in elem] for vec in
    #                      elem]) + histogram_bin_size) // Transit_profile_bin_size - len(transit_histogram_batch)),
    #                        'constant')
    #                 t_transit_batch = np.linspace(0, len(transit_histogram_batch) * Transit_profile_bin_size,
    #                                               len(transit_histogram_batch))
    #                 for n in range(len(all_transits_aligned_first)):
    #                     for m in range(len(all_transits_aligned_first[n]) - 1):
    #                         transit_histogram[all_transits_aligned_first[n][m + 1] // Transit_profile_bin_size] += 1
    #                         transit_histogram_batch[
    #                             all_transits_aligned_first[n][m + 1] // Transit_profile_bin_size] += 1
    #
    #     ############################################## END WHILE LOOP #################################################
    #
    #     ## Adding comment to measurement [prompt whether stopped or finished regularly]
    #     aftComment = pymsgbox.prompt('Add comment to measurement: ', default='', timeout=int(30e3))
    #     if aftComment == 'Timeout': aftComment = None
    #
    #     #### Handle file-names and directories #####
    #     ## Saving: np.savez(filedir, data = x) #note: @filedir is string of full directory; data is the queyword used to read @x from the file:
    #     ## Loading: file = np.load(f, allow_pickle = True)
    #     ##          x = file['data']
    #
    #     #### ------ Save results ------
    #     #  -------   Create dir
    #     root_dirname = f'U:\\Lab_2021-2022\\Experiment_results\\Transits\\{datest}\\'
    #     dirname = root_dirname + f'{timest}_Photon_TimeTags\\'  # Specific experiment dir
    #     dirname_N = dirname + 'North\\'
    #     dirname_S = dirname + 'South\\'
    #     if not os.path.exists(dirname):
    #         os.makedirs(dirname)
    #     if not os.path.exists(dirname_N):
    #         os.makedirs(dirname_N)
    #     if not os.path.exists(dirname_S):
    #         os.makedirs(dirname_S)
    #
    #     # ----  msmnt files names  -----
    #     # Counter_str = (Counter)
    #     filename_N_tt = f'North_timetags.npz'
    #     filename_S_tt = f'South_timetags.npz'
    #     filename_S_transits = f'South_Transits.npz'
    #     filename_FLR = f'Flouresence.npz'
    #     filename_timestamp = f'Drops_time_stamps.npz'
    #
    #     if len(FLR_measurement) > 0:
    #         np.savez(dirname + filename_FLR, FLR_measurement)
    #     if len(Exp_timestr_batch) > 0:
    #         np.savez(dirname + filename_timestamp, Exp_timestr_batch)
    #     if len(tt_N_measure_batch) > 0:
    #         np.savez(dirname_N + filename_N_tt, tt_N_measure_batch)
    #     if len(self.tt_S_measure_batch) > 0:
    #         np.savez(dirname_S + filename_S_tt, self.tt_S_measure_batch)
    #     if len(all_transits_batch) > 0:
    #         np.savez(dirname_S + filename_S_transits, all_transits_batch)
    #
    #     ### Edit comments file ####
    #     Counter_str = str(Counter)
    #     cmntDir = root_dirname + '\\daily_experiment_comments.txt'
    #     cmnt = timest + ' - '
    #     if preComment is not None: cmnt = cmnt + preComment + '; '
    #     if aftComment is not None: cmnt = cmnt + aftComment
    #     if preComment is None and aftComment is None: cmnt = cmnt + 'No comment. '
    #     try:
    #         with open(cmntDir, "a") as commentsFile:
    #             commentsFile.write(cmnt + '\n')
    #     except:
    #         print('Could not save comments, error writing to comments-file.')
    #         print(cmnt)
    #
    #     comments = {'comments': cmnt}
    #     try:
    #         with open(f'{dirname}experiment_comments.txt', 'w') as file:
    #             json.dump(comments, file, indent=4)
    #     except Exception:
    #         pass
    #
    #     self.saveConfigTable(path=dirname)
    #     ## ------------------ end of saving section -------
    def Save_SNSPDs_CRUS_Measurement_with_tt(self, N, histogram_bin_size, Transit_profile_bin_size, preComment,
                                             total_counts_threshold, transit_counts_threshold, transit_time_threshold,
                                             bandwidth, freq_step, max_probe_counts, Mock=False):
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
        :param total_counts_threshold: The minimum number of MCounts / sec for which the data will be displayed and collected.
        :param transit_counts_threshold: The number of photons at each time bin for which suspect we as a transit.
        :param transit_time_threshold: The maximum time between the start and finish of the transit in [nsec]
        :param max_probe_counts: The maximum counts probe counts at each direction measured when cavity isn't locked.
        :return:
        """
        # if preComment is True:
        if not preComment:
            preComment = pymsgbox.prompt('Add comment to measurement: ', default='', timeout=int(30e3))
        aftComment = None

        ### fetching data from server
        ### saving to file
        ###

        histogram_bin_number = self.M_window // (histogram_bin_size)
        time_bins = np.linspace(0, self.M_window, histogram_bin_number)
        spectrum_bin_number = bandwidth // freq_step + 1  # TODO: ask natan how many freq steps he uses
        freq_bins = np.linspace(-int(bandwidth / 2), int(bandwidth / 2), spectrum_bin_number)
        CRUS_pulse_time = np.arange(512)

        ## Listen for keyboard
        listener = keyboard.Listener(on_press=self.on_key_press)
        listener.start()  # start to listen on a separate thread
        self.keyPress = None
        print('\033[94m' + 'Press ESC to stop measurement.' + '\033[0m')  # print blue
        reps = 1  # a counter, number of repeats actually made.

        # Probe_S_handle = self.job.result_handles.get("South_Probe")
        # tt_S_handle = self.job.result_handles.get("South_Probe_TT")
        Counts_5_handle = self.job.result_handles.get("Det1_Counts")
        Counts_6_handle = self.job.result_handles.get("Det2_Counts")
        Counts_7_handle = self.job.result_handles.get("Det3_Counts")
        Counts_8_handle = self.job.result_handles.get("Det4_Counts")
        tt_5_handle = self.job.result_handles.get("Det1_Probe_TT")
        tt_6_handle = self.job.result_handles.get("Det2_Probe_TT")
        tt_7_handle = self.job.result_handles.get("Det3_Probe_TT")
        tt_8_handle = self.job.result_handles.get("Det4_Probe_TT")
        FLR_handle = self.job.result_handles.get("FLR_measure")

        # define empty variables
        self.tt_5_measure_batch = []
        self.tt_6_measure_batch = []
        self.tt_7_measure_batch = []
        self.tt_8_measure_batch = []
        self.tt_S_measure_batch = []
        tt_S_binning_batch = []
        tt_S_binning_resonance_batch = []
        tt_S_binning_detuned_batch = []
        all_transits_batch = []
        FLR_measurement = []
        Exp_timestr_batch = []

        tt_S_binning_resonance = np.zeros(histogram_bin_number)
        tt_S_binning_detuned = np.zeros(histogram_bin_number)
        tt_S_transit_events = np.zeros(histogram_bin_number)
        self.tt_S_CRUS_events_batch = np.zeros(histogram_bin_size)
        Cavity_atom_spectrum = np.zeros(spectrum_bin_number)
        Cavity_spectrum = np.zeros(spectrum_bin_number)

        start = True

        # Place holders for results
        # while (number of photons * 10^-6 [Mcounts] / Measuring time [nsec] * 10^-9 [sec/nsec])  > counts_threshold [Mcounts /sec]:
        while (((sum(tt_S_binning_resonance) * 1000) / (self.M_window / 2)) > total_counts_threshold) or start:
            if self.keyPress == 'ESC':
                print('\033[94m' + 'ESC pressed. Stopping measurement.' + '\033[0m')  # print blue
                self.updateValue("CRUS_Exp_switch", False)
                self.update_parameters()
                # Other actions can be added here
                break
            if start:
                start = False
            else:
                print('Above Threshold')

            tt_5_handle.wait_for_values(1)
            tt_6_handle.wait_for_values(1)
            tt_7_handle.wait_for_values(1)
            tt_8_handle.wait_for_values(1)
            FLR_handle.wait_for_values(1)

            counts_res5 = Counts_5_handle.fetch_all()
            counts_res6 = Counts_6_handle.fetch_all()
            counts_res7 = Counts_7_handle.fetch_all()
            counts_res8 = Counts_8_handle.fetch_all()
            tt_5_res = tt_5_handle.fetch_all()
            tt_6_res = tt_6_handle.fetch_all()
            tt_7_res = tt_7_handle.fetch_all()
            tt_8_res = tt_8_handle.fetch_all()
            FLR_res = -FLR_handle.fetch_all()

            self.tt_5_measure = [tt_5_res[(index * Config.vec_size): (index * Config.vec_size + counts)].tolist() for
                                 index, counts in
                                 enumerate(counts_res5)]
            self.tt_6_measure = [tt_6_res[(index * Config.vec_size): (index * Config.vec_size + counts)].tolist() for
                                 index, counts in
                                 enumerate(counts_res6)]
            self.tt_7_measure = [tt_7_res[(index * Config.vec_size): (index * Config.vec_size + counts)].tolist() for
                                 index, counts in
                                 enumerate(counts_res7)]
            self.tt_8_measure = [tt_8_res[(index * Config.vec_size): (index * Config.vec_size + counts)].tolist() for
                                 index, counts in
                                 enumerate(counts_res8)]
            self.tt_S_measure = self.tt_5_measure[0] + self.tt_6_measure[0] + self.tt_7_measure[0] + self.tt_8_measure[
                0]
            self.tt_5_measure.sort()
            self.tt_6_measure.sort()
            self.tt_7_measure.sort()
            self.tt_8_measure.sort()
            self.tt_S_measure.sort()

            self.tt_S_binning = np.zeros(histogram_bin_number * 2)
            self.tt_S_CRUS_events = np.zeros(histogram_bin_size)

            for x in [elem for elem in self.tt_S_measure if elem < self.M_window]:
                self.tt_S_binning[x // int(histogram_bin_size / 2)] += 1
                self.tt_S_CRUS_events[x % histogram_bin_size] += 1
                self.tt_S_CRUS_events_batch[x % histogram_bin_size] += 1

            # split the binning vector to odd and even - on and off resonance pulses
            tt_S_binning_resonance = [self.tt_S_binning[x] for x in range(len(self.tt_S_binning)) if x % 2]  # odd
            tt_S_binning_detuned = [self.tt_S_binning[x] for x in range(len(self.tt_S_binning)) if not x % 2]  # even

        ## record time
        timest = time.strftime("%Y%m%d-%H%M%S")
        datest = time.strftime("%Y%m%d")

        FLR_measurement = FLR_measurement[-(N - 1):] + [FLR_res.tolist()]
        Exp_timestr_batch = Exp_timestr_batch[-(N - 1):] + [timest]

        self.tt_5_measure_batch = self.tt_5_measure_batch[-(N - 1):] + [self.tt_5_measure]
        self.tt_6_measure_batch = self.tt_6_measure_batch[-(N - 1):] + [self.tt_6_measure]
        self.tt_7_measure_batch = self.tt_7_measure_batch[-(N - 1):] + [self.tt_7_measure]
        self.tt_8_measure_batch = self.tt_8_measure_batch[-(N - 1):] + [self.tt_8_measure]
        self.tt_S_measure_batch = self.tt_S_measure_batch[-(N - 1):] + [self.tt_S_measure]
        tt_S_binning_batch = tt_S_binning_batch[-(N - 1):] + [self.tt_S_binning]
        tt_S_binning_resonance_batch = tt_S_binning_resonance_batch[-(N - 1):] + [tt_S_binning_resonance]
        tt_S_binning_detuned_batch = tt_S_binning_detuned_batch[-(N - 1):] + [tt_S_binning_detuned]
        Counter = 1

        # Find transits and build histogram:
        current_transit = []
        all_transits = []
        all_transits_aligned_first = []
        t_transit = []
        t_transit_batch = []
        transit_histogram = []
        for index, value in enumerate(tt_S_binning_resonance):
            if not current_transit and value:  # if the array is empty
                current_transit.append(index)
            if value:
                if ((index - current_transit[0]) * histogram_bin_size) < transit_time_threshold:
                    current_transit.append(index)
                elif len(current_transit) > transit_counts_threshold:
                    all_transits.append(current_transit)
                    all_transits_aligned_first.append([x - current_transit[0] for x in current_transit])
                    # tt_S_transit_events[tuple(current_transit)] += 1
                    current_transit = [index]
                else:
                    # Finding if there any index that was saved to current transit and is close enough to the new index
                    t = [i for i, elem in enumerate(current_transit) if
                         ((index - elem) * histogram_bin_size) < transit_time_threshold]
                    if t:
                        current_transit = current_transit[t[0]:] + [index]
                    else:
                        current_transit = [index]

        tt_S_transit_events[[i for i in [vec for elem in all_transits for vec in elem]]] += 1
        #
        if all_transits:
            all_transits_batch = all_transits_batch[-(N - 1):] + [all_transits]

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
                self.updateValue("CRUS_Exp_switch", False)
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
            textstr_detuned = r'$Probe_N = %.3f$' % (
            (sum(tt_S_binning_detuned) * 1000) / (self.M_window / 2),) + '[MPhotons/sec]\n' \
                              + '$\overline{Probe}_S = %.3f$' % (
                              (np.mean([sum(x) for x in tt_S_binning_detuned_batch]) * 1000)
                              / (self.M_window / 2),) + '[MPhotons/sec]'
            textstr_resonance = r'$Probe_S = %.3f$' % (
            (sum(tt_S_binning_resonance) * 1000) / (self.M_window / 2),) + '[MPhotons/sec]\n' \
                                + '$\overline{Probe}_S = %.3f$' % (
                                (np.mean([sum(x) for x in tt_S_binning_resonance_batch]) * 1000)
                                / (self.M_window / 2),) + '[MPhotons/sec]'
            textstr_FLR = r'$\overline{FLR}_{MAX} = %.1f$' % (np.mean(FLR_measurement) * 1e5,) + r'$\times 10^{-5}$'
            textstr_No_transits = 'NO TRANSITS YET!!!'

            ax1.plot(time_bins, tt_S_binning_detuned_batch[-1], label='Counts histogram', color='b')

            ax1.set_title('Detuned counts', fontweight="bold")
            ax1.set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
            ax1.text(0.05, 0.95, textstr_detuned, transform=ax1.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
            ax1.legend(loc='upper right')
            print('ok')
            ax2.plot(time_bins, tt_S_binning_resonance_batch[-1], label='Counts histogram', color='b')

            ax2.set_title('On resonant counts', fontweight="bold")
            ax2.set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
            ax2.text(0.05, 0.95, textstr_resonance, transform=ax2.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
            ax2.legend(loc='upper right')

            ax3.plot(CRUS_pulse_time, self.tt_S_CRUS_events, label='Folded photon detection events per cycle',
                     color='k')
            ax3.set(xlabel='Time [nsec]', ylabel='Counts [#Number]')
            ax3.legend(loc='upper right')

            ax4.plot(CRUS_pulse_time, self.tt_S_CRUS_events_batch, label='Folded photon detection events per run',
                     color='k')
            ax4.set(xlabel='Time [msec]', ylabel='Counts [#Number]')
            ax4.legend(loc='upper right')

            if len(all_transits_batch) > 0:
                if all_transits:
                    textstr_transit_counts = r'$N_{Transits} = %s $' % (len(all_transits),) + r'$[Counts]$'
                textstr_transit_event_counter = r'$N_{Transits Total} = %s $' % (
                len([vec for elem in all_transits_batch for vec in elem]),) + r'$[Counts]$'

                ax5.plot(t_transit, transit_histogram, label='Transit profile', color='b')
                ax5.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
                ax5.text(0.05, 0.95, textstr_transit_counts, transform=ax5.transAxes, fontsize=12,
                         verticalalignment='top', bbox=props)

                ax6.plot(time_bins, tt_S_transit_events, label='Transit events histogram', marker='*', color='b')
                ax6.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
                ax6.text(0.05, 0.95, textstr_transit_event_counter, transform=ax6.transAxes, fontsize=12,
                         verticalalignment='top', bbox=props)
            else:
                ax5.plot(t_transit, transit_histogram, label='Transit profile', color='b')
                ax5.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
                ax5.text(0.25, 0.5, textstr_No_transits, transform=ax5.transAxes, fontsize=24,
                         verticalalignment='center', bbox=props)

                ax6.plot(time_bins, tt_S_transit_events, label='Transit events histogram', marker='*', color='b')
                ax6.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
                ax6.text(0.25, 0.5, textstr_No_transits, transform=ax6.transAxes, fontsize=24,
                         verticalalignment='center', bbox=props)

            ax1.set_ylim(0, 8)
            ax2.set_ylim(0, 8)

            # plt.tight_layout()
            plt.show()
            plt.pause(0.5)

            ###########################################################################################################

            while True:
                # record time:
                timest = time.strftime("%Y%m%d-%H%M%S")
                datest = time.strftime("%Y%m%d")

                tt_5_handle.wait_for_values(1)
                tt_6_handle.wait_for_values(1)
                tt_7_handle.wait_for_values(1)
                tt_8_handle.wait_for_values(1)
                FLR_handle.wait_for_values(1)

                counts_res5 = Counts_5_handle.fetch_all()
                counts_res6 = Counts_6_handle.fetch_all()
                counts_res7 = Counts_7_handle.fetch_all()
                counts_res8 = Counts_8_handle.fetch_all()
                tt_5_res = tt_5_handle.fetch_all()
                tt_6_res = tt_6_handle.fetch_all()
                tt_7_res = tt_7_handle.fetch_all()
                tt_8_res = tt_8_handle.fetch_all()
                FLR_res = -FLR_handle.fetch_all()

                self.tt_5_measure = [tt_5_res[(index * Config.vec_size): (index * Config.vec_size + counts)].tolist()
                                     for index, counts in enumerate(counts_res5)]
                self.tt_6_measure = [tt_6_res[(index * Config.vec_size): (index * Config.vec_size + counts)].tolist()
                                     for index, counts in enumerate(counts_res6)]
                self.tt_7_measure = [tt_7_res[(index * Config.vec_size): (index * Config.vec_size + counts)].tolist()
                                     for index, counts in enumerate(counts_res7)]
                self.tt_8_measure = [tt_8_res[(index * Config.vec_size): (index * Config.vec_size + counts)].tolist()
                                     for index, counts in enumerate(counts_res8)]
                self.tt_S_measure = self.tt_5_measure[0] + self.tt_6_measure[0] + self.tt_7_measure[0] + \
                                    self.tt_8_measure[0]
                self.tt_5_measure.sort()
                self.tt_6_measure.sort()
                self.tt_7_measure.sort()
                self.tt_8_measure.sort()
                self.tt_S_measure.sort()

                if self.tt_S_measure != self.tt_S_measure_batch[-1]:
                    break

            self.tt_S_binning = np.zeros(histogram_bin_number * 2)

            for x in [elem for elem in self.tt_S_measure if elem < self.M_window]:
                self.tt_S_binning[x // int(histogram_bin_size / 2)] += 1

            # split the binning vector to odd and even - on and off resonance pulses
            tt_S_binning_resonance = [self.tt_S_binning[x] for x in range(len(self.tt_S_binning)) if x % 2]  # odd
            tt_S_binning_detuned = [self.tt_S_binning[x] for x in range(len(self.tt_S_binning)) if not x % 2]  # even

            if ((sum(tt_S_binning_resonance) * 1000) / (self.M_window / 2)) < total_counts_threshold:

                self.tt_S_CRUS_events = np.zeros(histogram_bin_size)

                for x in [elem for elem in self.tt_S_measure if elem <= self.M_window]:
                    self.tt_S_CRUS_events[x % histogram_bin_size] += 1
                    self.tt_S_CRUS_events_batch[x % histogram_bin_size] += 1

                FLR_measurement = FLR_measurement[-(N - 1):] + [FLR_res.tolist()]
                Exp_timestr_batch = Exp_timestr_batch[-(N - 1):] + [timest]
                if Counter < N:
                    Counter += 1
                print(timest, Counter)

                self.tt_5_measure_batch = self.tt_5_measure_batch[-(N - 1):] + [self.tt_5_measure]
                self.tt_6_measure_batch = self.tt_6_measure_batch[-(N - 1):] + [self.tt_6_measure]
                self.tt_7_measure_batch = self.tt_7_measure_batch[-(N - 1):] + [self.tt_7_measure]
                self.tt_8_measure_batch = self.tt_8_measure_batch[-(N - 1):] + [self.tt_8_measure]
                self.tt_S_measure_batch = self.tt_S_measure_batch[-(N - 1):] + [self.tt_S_measure]
                tt_S_binning_batch = tt_S_binning_batch[-(N - 1):] + [self.tt_S_binning]
                tt_S_binning_resonance_batch = tt_S_binning_resonance_batch[-(N - 1):] + [tt_S_binning_resonance]
                tt_S_binning_detuned_batch = tt_S_binning_detuned_batch[-(N - 1):] + [tt_S_binning_detuned]

                # Find transits and build histogram:
                current_transit = []
                all_transits = []
                all_transits_aligned_first = []
                for index, value in enumerate(tt_S_binning_resonance):
                    if not current_transit and value:  # if the array is empty
                        current_transit.append(index)
                    if value:
                        if ((index - current_transit[0]) * 480) < transit_time_threshold:
                            current_transit.append(index)
                        elif len(current_transit) > transit_counts_threshold:
                            all_transits.append(current_transit)
                            all_transits_aligned_first.append([x - current_transit[0] for x in current_transit])
                            current_transit = [index]
                        else:
                            # Finding if there any index that was saved to current transit and is close enough to the new index
                            t = [i for i, elem in enumerate(current_transit) if
                                 ((index - elem) * 480) < transit_time_threshold]
                            if t:
                                current_transit = current_transit[t[0]:] + [index]
                            else:
                                current_transit = [index]

                tt_S_transit_events[[i for i in [vec for elem in all_transits for vec in elem]]] += 1

                if all_transits:
                    all_transits_batch = all_transits_batch[-(N - 1):] + [all_transits]
            else:
                tt_S_binning_resonance = tt_S_binning_resonance_batch[-1]  # odd
                tt_S_binning_detuned = tt_S_binning_detuned_batch[-1]  # even

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
        filename_Det5_tt = f'Det5_timetags.npz'
        filename_Det6_tt = f'Det6_timetags.npz'
        filename_Det7_tt = f'Det7_timetags.npz'
        filename_Det8_tt = f'Det8_timetags.npz'
        filename_S_tt = f'South_timetags.npz'
        filename_S_transits = f'South_Transits.npz'
        filename_S_folded = f'South_timetags_folded_512ns.npz'
        filename_FLR = f'Flouresence.npz'
        filename_timestamp = f'Drops_time_stamps.npz'

        if len(FLR_measurement) > 0:
            np.savez(dirname + filename_FLR, FLR_measurement)
        if len(Exp_timestr_batch) > 0:
            np.savez(dirname + filename_timestamp, Exp_timestr_batch)
        if len(self.tt_S_CRUS_events_batch) > 0:
            np.savez(dirname + filename_S_folded, self.tt_S_CRUS_events_batch)
        if len(self.tt_5_measure_batch) > 0:
            np.savez(dirname_S + filename_Det5_tt, self.tt_5_measure_batch)
        if len(self.tt_6_measure_batch) > 0:
            np.savez(dirname_S + filename_Det6_tt, self.tt_6_measure_batch)
        if len(self.tt_7_measure_batch) > 0:
            np.savez(dirname_S + filename_Det7_tt, self.tt_7_measure_batch)
        if len(self.tt_8_measure_batch) > 0:
            np.savez(dirname_S + filename_Det8_tt, self.tt_8_measure_batch)
        if len(self.tt_S_measure_batch) > 0:
            np.savez(dirname_S + filename_S_tt, self.tt_S_measure_batch)
        if len(all_transits_batch) > 0:
            np.savez(dirname_S + filename_S_transits, all_transits_batch)

        ### Edit comments file ####
        cmntDir = root_dirname + '\\daily_experiment_comments.txt'
        cmnt = timest + ' - ' + 'max probe counts:'  # +max_probe_counts+'-'
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

        self.saveConfigTable(path=dirname)
        try:
            with open(f'{dirname}max_probe_counts.txt', 'w') as file:
                json.dump(max_probe_counts, file, indent=4)
        except Exception as e:
            print(e)
        for qrdCtrl in self.QuadRFControllers:
            qrdCtrl.saveLinesAsCSV(f'{dirname}QuadRF_table.csv')
        ## ------------------ end of saving section -------

    def Save_SNSPDs_Sprint_Measurement_with_tt(self, N, histogram_bin_size, Transit_profile_bin_size, preComment,
                                                 lock_err_threshold, transit_counts_threshold, transit_time_threshold,
                                                 bandwidth, freq_step, max_probe_counts, Mock = False):
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
        :param lock_err_threshold: The maximalerror in locking resonator to rb line.
        :param transit_counts_threshold: The number of photons at each time bin for which suspect we as a transit.
        :param transit_time_threshold: The maximum time between the start and finish of the transit in [nsec]
        :param max_probe_counts: The maximum counts probe counts at each direction measured when cavity isn't locked.
        :return:
        """
        # if preComment is True:
        if not preComment:
            preComment = pymsgbox.prompt('Add comment to measurement: ', default='', timeout=int(30e3))
        aftComment = None

        ### fetching data from server
        ### saving to file
        ###
        Num_Of_dets = 6
        # detector_delay = [5,0,0,15] # For detectors 1-4 "N"
        detector_delay = [0,0,0,0] # For detectors 5-8 "S"


        histogram_bin_number = self.M_window // (histogram_bin_size)
        time_bins = np.linspace(0, self.M_window, histogram_bin_number)
        spectrum_bin_number = bandwidth // freq_step + 1  # TODO: ask natan how many freq steps he uses
        freq_bins = np.linspace(-int(bandwidth / 2), int(bandwidth / 2), spectrum_bin_number)
        SPRINT_pulse_time = np.arange(512)

        ## Listen for keyboard
        listener = keyboard.Listener(on_press=self.on_key_press)
        listener.start()  # start to listen on a separate thread
        self.keyPress = None
        print('\033[94m' + 'Press ESC to stop measurement.' + '\033[0m')  # print blue
        reps = 1  # a counter, number of repeats actually made.

        ####     get tt and counts from OPX to python   #####

        # Probe_S_handle = self.job.result_handles.get("South_Probe")
        # tt_S_handle = self.job.result_handles.get("South_Probe_TT")
        Counts_handle = []
        tt_handle = []
        for i in range(Num_Of_dets):
            Counts_handle.append(self.job.result_handles.get("Det"+str(i+1)+"_Counts"))
            tt_handle.append(self.job.result_handles.get("Det"+str(i+1)+"_Probe_TT"))

        FLR_handle = self.job.result_handles.get("FLR_measure")

        # define empty variables
        tt_S_binning_batch = []

        all_transits_batch = []
        FLR_measurement = []
        Exp_timestr_batch = []

        tt_S_SPRINT_events = np.zeros(histogram_bin_number)
        self.tt_S_SPRINT_events_batch = np.zeros(histogram_bin_size)

        start = True

        # take threshold from npz ( error from resonator lock PID)
        lock_err = np.load(
            'U:\Lab_2021-2022\Experiment_results\Sprint\Locking_PID_Error\locking_err.npy')  # the error of locking the resontor to Rb line
        # lock_err = lock_err_threshold/2

        # Place holders for results
        # while (number of photons * 10^-6 [Mcounts] / Measuring time [nsec] * 10^-9 [sec/nsec])  > counts_threshold [Mcounts /sec]:
        while lock_err > lock_err_threshold or start:
            if self.keyPress == 'ESC':
                print('\033[94m' + 'ESC pressed. Stopping measurement.' + '\033[0m')  # print blue
                self.updateValue("Sprint_Exp_switch", False)
                self.update_parameters()
                # Other actions can be added here
                break

            if start:
                start = False
            else:
                print('Above Threshold')

            # handles wait for values
            for i in range(Num_Of_dets):
                tt_handle[i].wait_for_values(1)
            FLR_handle.wait_for_values(1)

            # add the tt and counts to python vars
            counts_res = []
            tt_res = []
            for i in range(Num_Of_dets):
                counts_res.append(Counts_handle[i].fetch_all())
                tt_res.append(tt_handle[i].fetch_all())
            FLR_res = -FLR_handle.fetch_all()

            self.tt_measure = []
            self.tt_S_measure = []
            for i in range(Num_Of_dets): # for different detectors
                self.tt_measure.append([tt_res[i][(index * Config.vec_size): (index * Config.vec_size + counts)].tolist() for index, counts in
                                        enumerate(counts_res[i])])
                # for window in range(len(self.tt_measure[i])):
                #     self.tt_measure[i][window] = [elem + detector_delay[i] for elem in self.tt_measure[i][window]]
                self.tt_S_measure += self.tt_measure[i][-1]
                self.tt_measure[i].sort()
            self.tt_S_measure.sort()
            ####    end get tt and counts from OPX to python   #####

            self.tt_S_binning = np.zeros(histogram_bin_number +1)
            self.tt_S_SPRINT_events = np.zeros(histogram_bin_size)
            self.tt_S_SPRINT_events_batch = np.zeros(histogram_bin_size)
            self.tt_Single_det_SPRINT_events = np.zeros((Num_Of_dets,histogram_bin_size))
            self.tt_Single_det_SPRINT_events_batch = np.zeros((Num_Of_dets,histogram_bin_size))

            # fold South:
            # for x in [elem for elem in self.tt_S_measure if elem < self.M_window]: - for debugging assaf
            for x in [elem for elem in self.tt_S_measure]:
                self.tt_S_binning[x // int(histogram_bin_size)] += 1
                self.tt_S_SPRINT_events[x % histogram_bin_size] += 1
                self.tt_S_SPRINT_events_batch[x % histogram_bin_size] += 1

            # fold for different detectors:
            for i in range(Num_Of_dets):
                # for x in [elem for elem in self.tt_S_measure if elem < self.M_window]: - for debugging assaf
                for x in [elem for elem in self.tt_measure[i][-1]]:
                    self.tt_Single_det_SPRINT_events[i][x % histogram_bin_size] += 1
                    self.tt_Single_det_SPRINT_events_batch[i][x % histogram_bin_size] += 1

        ## record time
        timest = time.strftime("%Y%m%d-%H%M%S")
        datest = time.strftime("%Y%m%d")

        FLR_measurement = FLR_measurement[-(N - 1):] + [FLR_res.tolist()]
        Exp_timestr_batch = Exp_timestr_batch[-(N - 1):] + [timest]

        self.tt_measure_batch = []
        self.tt_S_measure_batch = []
        for i in range(Num_Of_dets):
            self.tt_measure_batch.append([self.tt_measure[i]])
        self.tt_S_measure_batch = self.tt_S_measure_batch[-(N - 1):] + [self.tt_S_measure]
        tt_S_binning_batch = tt_S_binning_batch[-(N - 1):] + [self.tt_S_binning]
        Counter = 1
        #
        # ### Find transits and build histogram:  ###
        # current_transit = []
        # all_transits = []
        # all_transits_aligned_first = []
        # t_transit = []
        # t_transit_batch = []
        # transit_histogram = []
        # for index, value in enumerate(tt_S_binning_resonance):
        #     if not current_transit and value:  # if the array is empty
        #         current_transit.append(index)
        #     if value:
        #         if ((index - current_transit[0]) * histogram_bin_size) < transit_time_threshold:
        #             current_transit.append(index)
        #         elif len(current_transit) > transit_counts_threshold:
        #             all_transits.append(current_transit)
        #             all_transits_aligned_first.append([x - current_transit[0] for x in current_transit])
        #             # tt_S_transit_events[tuple(current_transit)] += 1
        #             current_transit = [index]
        #         else:
        #             # Finding if there any index that was saved to current transit and is close enough to the new index
        #             t = [i for i, elem in enumerate(current_transit) if ((index - elem) * histogram_bin_size) < transit_time_threshold]
        #             if t:
        #                 current_transit = current_transit[t[0]:] + [index]
        #             else:
        #                 current_transit = [index]
        #
        #
        # tt_S_transit_events[[i for i in [vec for elem in all_transits for vec in elem]]] += 1
        # #
        # if all_transits:
        #     all_transits_batch = all_transits_batch[-(N - 1):] + [all_transits]
        #
        #

        #
        # ## Prepare plots
        fig = plt.figure()
        # fig2 = plt.figure()
        #
        # ax1 =
        # ax1 = plt.subplot2grid((6, 2), (0, 0), colspan=1, rowspan=2)
        # ax2 = plt.subplot2grid((6, 2), (0, 1), colspan=1, rowspan=2)
        # ax3 = plt.subplot2grid((6, 2), (2, 0), colspan=1, rowspan=2)
        # ax4 = plt.subplot2grid((6, 2), (2, 1), colspan=1, rowspan=2)
        # # ax5 = plt.subplot2grid((6, 2), (4, 0), colspan=1, rowspan=2)
        # ax6 = plt.subplot2grid((6, 2), (4, 1), colspan=1, rowspan=2)
        #
        while True:
            if self.keyPress == 'ESC':
                print('\033[94m' + 'ESC pressed. Stopping measurement.' + '\033[0m')  # print blue
                self.updateValue("CRUS_Exp_switch", False)
                self.update_parameters()
                # Other actions can be added here
                break
            if reps < N:
                reps += 1
            #
            #     ########################## PLOT!!! ########################################################################
            #
            fig.clear()
            #     ax2.clear()
            #     ax3.clear()
            #     ax4.clear()
            #     ax5.clear()
            #     ax6.clear()
            #
            #     props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            #     textstr_detuned = r'$Probe_N = %.3f$' % ((sum(tt_S_binning_detuned) * 1000) / (self.M_window / 2),) + '[MPhotons/sec]\n' \
            #                + '$\overline{Probe}_S = %.3f$' % ((np.mean([sum(x) for x in tt_S_binning_detuned_batch]) * 1000)
            #                                                   / (self.M_window / 2),) + '[MPhotons/sec]'
            #     textstr_resonance = r'$Probe_S = %.3f$' % ((sum(tt_S_binning_resonance) * 1000) / (self.M_window / 2),) + '[MPhotons/sec]\n' \
            #                 + '$\overline{Probe}_S = %.3f$' % ((np.mean([sum(x) for x in tt_S_binning_resonance_batch]) * 1000)
            #                                                    / (self.M_window / 2),) + '[MPhotons/sec]'
            #     textstr_FLR = r'$\overline{FLR}_{MAX} = %.1f$' % (np.mean(FLR_measurement) * 1e5,) + r'$\times 10^{-5}$'
            #     textstr_No_transits = 'NO TRANSITS YET!!!'

            # plt.plot(self.tt_S_SPRINT_events_batch, label='Counts histogram', color='b')
            for i in range(4):
                plt.plot(self.tt_Single_det_SPRINT_events_batch[i], label='detector'+str(i+1+2*(i//2)))
            plt.legend(loc='upper right')
            plt.show()
            plt.pause(0.5)

            # fig2.clear()
            # plt.figure()
            # plt.plot(tt_S_binning_batch, label='unfolded data')
            # plt.show()
            # plt.pause(0.5)
            # ax1.plot(time_bins, tt_S_binning_batch[-1], label='Counts histogram', color='b')
            #
            #     ax1.set_title('Detuned counts', fontweight="bold")
            #     ax1.set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
            #     ax1.text(0.05, 0.95, textstr_detuned, transform=ax1.transAxes, fontsize=12,
            #              verticalalignment='top', bbox=props)
            #     ax1.legend(loc='upper right')
            #     print('ok')
            # ax2.plot(time_bins, tt_S_binning_batch[-1], label='Counts histogram', color='b')
            #
            #     ax2.set_title('On resonant counts', fontweight="bold")
            #     ax2.set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
            #     ax2.text(0.05, 0.95, textstr_resonance, transform=ax2.transAxes, fontsize=12,
            #              verticalalignment='top', bbox=props)
            #     ax2.legend(loc='upper right')
            #
            # ax3.plot(SPRINT_pulse_time, self.tt_S_SPRINT_events, label='Folded photon detection events per cycle', color='k')
            #     ax3.set(xlabel='Time [nsec]', ylabel='Counts [#Number]')
            #     ax3.legend(loc='upper right')
            #
            # ax4.plot(SPRINT_pulse_time, self.tt_S_SPRINT_events, label='Folded photon detection events per run', color='k')
            #     ax4.set(xlabel='Time [msec]', ylabel='Counts [#Number]')
            #     ax4.legend(loc='upper right')
            #
            #     if len(all_transits_batch) > 0:
            #         if all_transits:
            #             textstr_transit_counts = r'$N_{Transits} = %s $' % (len(all_transits),) + r'$[Counts]$'
            #         textstr_transit_event_counter = r'$N_{Transits Total} = %s $' % (len([vec for elem in all_transits_batch for vec in elem]),) + r'$[Counts]$'
            #
            #         ax5.plot(t_transit, transit_histogram, label='Transit profile', color='b')
            #         ax5.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
            #         ax5.text(0.05, 0.95, textstr_transit_counts, transform=ax5.transAxes, fontsize=12,
            #                  verticalalignment='top', bbox=props)
            #
            #         ax6.plot(time_bins, tt_S_transit_events, label='Transit events histogram', marker='*', color='b')
            #         ax6.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
            #         ax6.text(0.05, 0.95, textstr_transit_event_counter, transform=ax6.transAxes, fontsize=12,
            #                  verticalalignment='top', bbox=props)
            #     else:
            #         ax5.plot(t_transit, transit_histogram, label='Transit profile', color='b')
            #         ax5.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
            #         ax5.text(0.25, 0.5, textstr_No_transits, transform=ax5.transAxes, fontsize=24,
            #                  verticalalignment='center', bbox=props)
            #
            #         ax6.plot(time_bins, tt_S_transit_events, label='Transit events histogram', marker='*', color='b')
            #         ax6.set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
            #         ax6.text(0.25, 0.5, textstr_No_transits, transform=ax6.transAxes, fontsize=24,
            #                  verticalalignment='center', bbox=props)
            #
            #     ax1.set_ylim(0, 8)
            #     ax2.set_ylim(0, 8)
            #
            #     # plt.tight_layout()

            # #
            #     ###########################################################################################################

            while True:
                # record time:
                timest = time.strftime("%Y%m%d-%H%M%S")
                datest = time.strftime("%Y%m%d")

                # handles wait for values
                for i in range(Num_Of_dets):
                    tt_handle[i].wait_for_values(1)
                FLR_handle.wait_for_values(1)

                tt_res = []
                counts_res = []
                # add the tt and counts to python vars
                for i in range(Num_Of_dets):
                    counts_res.append(Counts_handle[i].fetch_all())
                    tt_res.append(tt_handle[i].fetch_all())
                FLR_res = -FLR_handle.fetch_all()

                self.tt_S_measure = []
                self.tt_measure = []
                for i in range(Num_Of_dets):  # for different detectors
                    self.tt_measure.append(
                        [tt_res[i][(index * Config.vec_size):
                                   (index * Config.vec_size + counts)].tolist() for index, counts in
                         enumerate(counts_res[i])])
                    #add delay to a detctors tt's
                    # for window in range(len(self.tt_measure[i])):
                    #     self.tt_measure[i][window] = [elem + detector_delay[i] for elem in self.tt_measure[i][window]]
                    self.tt_S_measure += self.tt_measure[i][-1]
                    self.tt_measure[i].sort()
                self.tt_S_measure.sort()

                if self.tt_S_measure != self.tt_S_measure_batch[-1]:
                    break
            # assaf - if x=self.M_window the index is out of range so i added 1
            self.tt_S_binning = np.zeros(histogram_bin_number+1) #  self.tt_S_binning = np.zeros(histogram_bin_number * 2)

            for x in [elem for elem in self.tt_S_measure if elem < self.M_window]:
                self.tt_S_binning[x // int(histogram_bin_size)] += 1

            if lock_err < lock_err_threshold:

                self.tt_S_SPRINT_events = np.zeros(histogram_bin_size)
                # for x in [elem for elem in self.tt_S_measure if elem <= self.M_window]:
                for x in [elem for elem in self.tt_S_measure]:
                    self.tt_S_SPRINT_events[x % histogram_bin_size] += 1
                    self.tt_S_SPRINT_events_batch[x % histogram_bin_size] += 1

                self.tt_Single_det_SPRINT_events = np.zeros((Num_Of_dets, histogram_bin_size))
                # fold for different detectors:
                for i in range(Num_Of_dets):
                    # for x in [elem for elem in self.tt_S_measure if elem < self.M_window]:
                    for x in [elem for elem in self.tt_measure[i][-1]]:
                        self.tt_Single_det_SPRINT_events[i][x % histogram_bin_size] += 1
                        self.tt_Single_det_SPRINT_events_batch[i][x % histogram_bin_size] += 1

                FLR_measurement = FLR_measurement[-(N - 1):] + [FLR_res.tolist()]
                Exp_timestr_batch = Exp_timestr_batch[-(N - 1):] + [timest]
                if Counter < N:
                    Counter += 1
                print(timest, Counter)

                for i in range(Num_Of_dets):
                    self.tt_measure_batch[i] = self.tt_measure_batch[i][-(N - 1):] + [self.tt_measure[i]]
                self.tt_S_measure_batch = self.tt_S_measure_batch[-(N - 1):] + [self.tt_S_measure]
                tt_S_binning_batch = tt_S_binning_batch[-(N - 1):] + [self.tt_S_binning]

            #
            #     # Find transits and build histogram:
            #     current_transit = []
            #     all_transits = []
            #     all_transits_aligned_first = []
            #     for index, value in enumerate(tt_S_binning_resonance):
            #         if not current_transit and value:  # if the array is empty
            #             current_transit.append(index)
            #         if value:
            #             if ((index - current_transit[0]) * 480) < transit_time_threshold:
            #                 current_transit.append(index)
            #             elif len(current_transit) > transit_counts_threshold:
            #                 all_transits.append(current_transit)
            #                 all_transits_aligned_first.append([x - current_transit[0] for x in current_transit])
            #                 current_transit = [index]
            #             else:
            #                 # Finding if there any index that was saved to current transit and is close enough to the new index
            #                 t = [i for i, elem in enumerate(current_transit) if
            #                      ((index - elem) * 480) < transit_time_threshold]
            #                 if t:
            #                     current_transit = current_transit[t[0]:] + [index]
            #                 else:
            #                     current_transit = [index]
            #
            #     tt_S_transit_events[[i for i in [vec for elem in all_transits for vec in elem]]] += 1
            #
            #     if all_transits:
            #         all_transits_batch = all_transits_batch[-(N - 1):] + [all_transits]
            # else:
            #     tt_S_binning_resonance = tt_S_binning_resonance_batch[-1]  # odd
            #     tt_S_binning_detuned = tt_S_binning_detuned_batch[-1]  # even

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
        root_dirname = f'U:\\Lab_2021-2022\\Experiment_results\\Sprint\\{datest}\\'
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
        filename_Det_tt = []
        for i in range(Num_Of_dets):
            filename_Det_tt.append(f'Det'+str(i+1)+f'_timetags.npz')
        filename_S_tt = f'South_timetags.npz'
        filename_S_transits = f'South_Transits.npz'
        filename_S_folded = f'South_timetags_folded_512ns.npz'
        filename_FLR = f'Flouresence.npz'
        filename_timestamp = f'Drops_time_stamps.npz'

        if len(FLR_measurement) > 0:
            np.savez(dirname + filename_FLR, FLR_measurement)
        if len(Exp_timestr_batch) > 0:
            np.savez(dirname + filename_timestamp, Exp_timestr_batch)
        if len(self.tt_S_SPRINT_events_batch) > 0:
            np.savez(dirname + filename_S_folded, self.tt_S_SPRINT_events_batch)
        for i in range(Num_Of_dets):
            if len(self.tt_measure_batch[i]) > 0:
                np.savez(dirname_S + filename_Det_tt[i], self.tt_measure_batch[i])
        if len(self.tt_S_measure_batch) > 0:
            np.savez(dirname_S + filename_S_tt, self.tt_S_measure_batch)
        # if len(all_transits_batch) > 0:
        #     np.savez(dirname_S + filename_S_transits, all_transits_batch)

        ### Edit comments file ####
        cmntDir = root_dirname + '\\daily_experiment_comments.txt'
        cmnt = timest + ' - '+'max probe counts:'#+max_probe_counts+'-'
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

        self.saveConfigTable(path=dirname)
        try:
            with open(f'{dirname}max_probe_counts.txt', 'w') as file: json.dump(max_probe_counts, file, indent=4)
        except Exception as e:
            print(e)
        for qrdCtrl in self.QuadRFControllers:
            qrdCtrl.saveLinesAsCSV(f'{dirname}QuadRF_table.csv')
        ## ------------------ end of saving section -------

    def Start_Transit_Exp(self, N=100, bin_size=100, preComment=None, threshold=0.1):
        self.Transit_Exp_switch(True)
        self.update_parameters()
        self.Save_SNSPDs_Measurement(N, bin_size, preComment, threshold)

    def Start_Transit_Exp_with_tt(self, N=100, Histogram_bin_size=1000, Transit_profile_bin_size=100, preComment=None,
                                  total_counts_threshold=1, transit_counts_threshold=5):
        # Max_probe_counts = self.Get_Max_Probe_counts(3)  # return the average maximum probe counts of 3 cycles.
        Max_probe_counts = None  # return the average maximum probe counts of 3 cycles.
        self.Transit_Exp_switch(True)
        self.update_parameters()
        self.Save_SNSPDs_Transit_Measurement_with_tt(N, Histogram_bin_size, Transit_profile_bin_size, preComment,
                                                     total_counts_threshold, transit_counts_threshold, Max_probe_counts)

    def Start_Sprint_Exp_with_tt(self, N=100, Histogram_bin_size=720,#int(len(Config.CRUS_pulser_samples)),
                                   Transit_profile_bin_size=100, preComment=None, lock_err_threshold=1,
                                   transit_counts_threshold=5, transit_time_threshold=6000, bandwidth=80,
                                   freq_step=4):
        # Max_probe_counts = self.Get_Max_Probe_counts(3)  # return the average maximum probe counts of 3 cycles.
        Max_probe_counts = None  # return the average maximum probe counts of 3 cycles.
        self.SPRINT_Exp_switch(True)
        self.update_parameters()
        self.Save_SNSPDs_Sprint_Measurement_with_tt(N, Histogram_bin_size, Transit_profile_bin_size, preComment,
                                                      lock_err_threshold, transit_counts_threshold,
                                                      transit_time_threshold, bandwidth, freq_step, Max_probe_counts)

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
        'Spectrum_Exp_switch': [bool, None, Spectrum_Exp_switch],
        "CRUS_Exp_switch": [bool, None, CRUS_Exp_switch],
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
        # experiment.Start_Sprint_Exp_with_tt(N=500, preComment='test')
    # except KeyboardInterrupt:
    #     experiment.job.halt()
    #     experiment.qmm.reset_data_processing()
    # finally:
    #     experiment.job.halt()
    #     experiment.qmm.reset_data_processing()
    # experiment.updateValue('PrePulse_duration', 5)
    # experiment.update_parameters()
    # experiment.Repeat_Measurement(10)
    # experiment.Repeat_Measurement(10)
