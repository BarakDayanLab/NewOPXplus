# from Config import config
import Config_with_SNSPDs_and_QuadRF_Sprint as Config
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
import math
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
          "AOM_2-2/3'", "AOM_2-3'_for_interference", "AOM_2-2'", "FLR_detection", "Measurement") # , "Dig_detectors") # , "PULSER_N", "PULSER_S")

    ## MOT build-up ##
    n = declare(int)
    m = declare(int)
    play("Detection" * amp(FLR * 0), "FLR_detection", duration=4)  # we dont know why this works, but Yoav from QM made us write this line to solve an alignment problem we had in the next 2 for loops
    with for_(n, 1, n <= mot_repetitions, n + 1):
        play("MOT" * amp(Config.AOM_0_Attenuation), "MOT_AOM_0")
        play("MOT" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-")
        play("MOT" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+")
        # play("Const_open", "PULSER_N")
        # play("OD_FS" * amp(0.5), "AOM_2-2/3'")
        play("AntiHelmholtz_MOT", "AntiHelmholtz_Coils")
    with for_(m, 1, m <= (mot_repetitions - 1), m + 1):
        measure("Detection", "FLR_detection", None, integration.full("Detection_opt", FLR, "out1"))
        # play("OD_FS" * amp(0.1), "AOM_2-3'_for_interference")

    align("Cooling_Sequence", "MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+", "AntiHelmholtz_Coils", "Zeeman_Coils",
          "AOM_2-2/3'", "AOM_2-2'", "FLR_detection", "Measurement") # , "Dig_detectors") #, "PULSER_N", "PULSER_S")

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
          "Measurement", "PULSER_N", "PULSER_S", "PULSER_LO") # , "Dig_detectors")

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


def Sprint_Exp(m_off_time, m_time, m_window, shutter_open_time,
               ON_counts_st1, ON_counts_st2, ON_counts_st3,
               ON_counts_st6, ON_counts_st7, ON_counts_st8,
               tt_st_1, tt_st_2, tt_st_3, tt_st_6, tt_st_7, tt_st_8, rep_st):
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
    counts6 = declare(int)
    counts7 = declare(int)
    counts8 = declare(int)

    tt_vec1 = declare(int, size=vec_size)
    tt_vec2 = declare(int, size=vec_size)
    tt_vec3 = declare(int, size=vec_size)
    tt_vec6 = declare(int, size=vec_size)
    tt_vec7 = declare(int, size=vec_size)
    tt_vec8 = declare(int, size=vec_size)

    n = declare(int)
    t = declare(int)
    m = declare(int)

    align("PULSER_N", "PULSER_S", "Dig_detectors", "AOM_2-2/3'", "PULSER_LO")
    play("OD", "AOM_2-2/3'", duration=shutter_open_time)
    play("Const_open", "PULSER_LO", duration=shutter_open_time)
    align("PULSER_N", "PULSER_S", "Dig_detectors", "AOM_2-2/3'", "PULSER_LO")
    play("OD", "AOM_2-2/3'", duration=(m_time + m_off_time))
    play("Const_open", "PULSER_LO", duration=(m_time + m_off_time))
    with for_(t, 0, t < (m_time + m_off_time) * 4, t + int(len(Config.Sprint_Exp_Gaussian_samples_S))): #assaf comment debbuging
        play("Sprint_experiment_pulses_S", "PULSER_S")
        play("Sprint_experiment_pulses_N", "PULSER_N")

    # wait(137, "Dig_detectors")
    # wait(117, "Dig_detectors") # For 20ns between pulses in sequence
    wait(297, "Dig_detectors") # For 20ns between pulses in sequence of only detections
    with for_(n, 0, n < m_time * 4, n + m_window):
        measure("readout_SPRINT", "Dig_detectors", None,
        # measure("readout_QRAM", "Dig_detectors", None,
                time_tagging.digital(tt_vec1, m_window, element_output="out1", targetLen=counts1),
                time_tagging.digital(tt_vec2, m_window, element_output="out2", targetLen=counts2),
                time_tagging.digital(tt_vec3, m_window, element_output="out3", targetLen=counts3),
                time_tagging.digital(tt_vec6, m_window, element_output="out6", targetLen=counts6),
                time_tagging.digital(tt_vec7, m_window, element_output="out7", targetLen=counts7),
                time_tagging.digital(tt_vec8, m_window, element_output="out8", targetLen=counts8),
                )

        ## Save Data: ##

        # Number of Photons (NOP) Count stream for each detector: ##
        save(counts1, ON_counts_st1)
        save(counts2, ON_counts_st2)
        save(counts3, ON_counts_st3)
        save(counts6, ON_counts_st6)
        save(counts7, ON_counts_st7)
        save(counts8, ON_counts_st8)


        with for_(m, 0, m < vec_size, m + 1):
            wait(1000)
            save(tt_vec1[m], tt_st_1)
            save(tt_vec2[m], tt_st_2)
            save(tt_vec3[m], tt_st_3)
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
        SPRINT_Exp_ON = declare(bool, value=False)

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
        ON_counts_st6 = declare_stream()
        ON_counts_st7 = declare_stream()
        ON_counts_st8 = declare_stream()
        tt_st_1 = declare_stream()
        tt_st_2 = declare_stream()
        tt_st_3 = declare_stream()
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
            with if_(SPRINT_Exp_ON):
                assign(x, (38688900 - 3106 + 656000 * 2) // 4) # TODO -  added 38688900 to fix new delay due to wait(1000) in saving sprint data with vector size 10000, should be fixed as well
            with else_():
                assign(x,  - 3106)
            FreeFall(FreeFall_duration - x, coils_timing)

            ##########################
            ## Measurement Sequence ##
            ##########################

            with if_(Trigger_Phase == 3):  # when trigger on PrePulse
                ## Trigger QuadRF Sequence #####################
                play("C_Seq", "Cooling_Sequence", duration=2500)
                ################################################

            wait(2500, "AOM_2-2'")
            with if_(SPRINT_Exp_ON):
                play("Depump", "AOM_2-2'", duration=(PrePulse_duration - shutter_open_time))
            with else_():
                play("Depump", "AOM_2-2'", duration=PrePulse_duration)
                # wait(PrePulse_duration, "Cooling_Sequence")
            align(*all_elements, "AOM_2-2/3'", "AOM_2-2'", "PULSER_N", "PULSER_S", "Dig_detectors")

            with if_(Trigger_Phase == 4):  # when trigger on pulse 1
                ## Trigger QuadRF Sequence #####################
                play("C_Seq", "Cooling_Sequence", duration=2500)
                ################################################
            ## For taking an image:
            with if_((Pulse_1_duration > 0) & ~SPRINT_Exp_ON):
                align(*all_elements)
                Pulse_with_prep(Pulse_1_duration, Pulse_1_decay_time, pulse_1_duration_0,
                                pulse_1_duration_minus, pulse_1_duration_plus)
                Measure(Pulse_1_duration)  # This triggers camera (Control 7)
                align(*all_elements)

            with if_((Pulse_1_duration > 0) & SPRINT_Exp_ON):
                align("Dig_detectors", "PULSER_N", "PULSER_S", "PULSER_LO")
                Sprint_Exp(M_off_time, Pulse_1_duration, obj.M_window, shutter_open_time,
                           ON_counts_st1, ON_counts_st2, ON_counts_st3,
                           ON_counts_st6, ON_counts_st7, ON_counts_st8,
                           tt_st_1, tt_st_2, tt_st_3, tt_st_6, tt_st_7, tt_st_8, rep_st)
                save(AntiHelmholtz_ON, AntiHelmholtz_ON_st)
                with if_(AntiHelmholtz_ON):
                    save(FLR, FLR_st)
                align("Dig_detectors", "PULSER_N", "PULSER_S", "PULSER_LO")

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
            ON_counts_st1.buffer(obj.rep).save('Det1_Counts')
            ON_counts_st2.buffer(obj.rep).save('Det2_Counts')
            ON_counts_st3.buffer(obj.rep).save('Det3_Counts')
            ON_counts_st6.buffer(obj.rep).save('Det6_Counts')
            ON_counts_st7.buffer(obj.rep).save('Det7_Counts')
            ON_counts_st8.buffer(obj.rep).save('Det8_Counts')
            (tt_st_1 + rep_st).buffer(obj.vec_size * obj.rep).save('Det1_Probe_TT')
            (tt_st_2 + rep_st).buffer(obj.vec_size * obj.rep).save('Det2_Probe_TT')
            (tt_st_3 + rep_st).buffer(obj.vec_size * obj.rep).save('Det3_Probe_TT')
            (tt_st_6 + rep_st).buffer(obj.vec_size * obj.rep).save('Det6_Probe_TT')
            (tt_st_7 + rep_st).buffer(obj.vec_size * obj.rep).save('Det7_Probe_TT')
            (tt_st_8 + rep_st).buffer(obj.vec_size * obj.rep).save('Det8_Probe_TT')
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
        # self.QuadRFControllers.append(QuadRFMOTController(MOGdevice=qrfContr.dev, initialValues={'Operation_Mode': 'Continuous', 'CH3_freq': '90MHz', 'CH3_amp': '31dbm'},
        #                                                   updateChannels=[3], debugging=False, continuous=False))  # updates values on QuadRF (uploads table)
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

    def get_avg_num_of_photons_in_det_pulse(self, det_pulse_len, sprint_sequence_delay, num_of_det_pulses, num_of_sprint_sequences):
        self.avg_num_of_photons_in_det_pulse = np.zeros(num_of_det_pulses)
        for i in range(num_of_det_pulses):
            detection_puls_ind =\
            list(np.arange((sprint_sequence_delay + i * det_pulse_len),(sprint_sequence_delay + (i + 1) * det_pulse_len)))
            self.avg_num_of_photons_in_det_pulse[i] = \
                np.sum(self.tt_S_SPRINT_events[detection_puls_ind]) / num_of_sprint_sequences

    def get_handles_from_OPX_Server(self,Num_Of_dets):
        '''
         gets handles of timetags and counts from OPX
        :return:
        '''
        Counts_handle = []
        tt_handle = []

        for i in Num_Of_dets:
            Counts_handle.append(self.job.result_handles.get("Det"+str(i)+"_Counts"))
            tt_handle.append(self.job.result_handles.get("Det"+str(i)+"_Probe_TT"))
        FLR_handle = self.job.result_handles.get("FLR_measure")
        return Counts_handle,tt_handle,FLR_handle

    def get_tt_from_handles(self, Num_Of_dets, Counts_handle, tt_handle, FLR_handle):
        counts_res = []
        tt_res = []

        # handles wait for values
        for i in range(len(Num_Of_dets)):
            tt_handle[i].wait_for_values(1)
        FLR_handle.wait_for_values(1)

        # add the tt and counts to python vars
        for i in range(len(Num_Of_dets)):
            counts_res.append(Counts_handle[i].fetch_all())
            tt_res.append(tt_handle[i].fetch_all())
        self.FLR_res = -FLR_handle.fetch_all()

        # clear tt's vecs from last data
        self.tt_measure = []
        self.tt_N_measure = []
        self.tt_S_measure = []

        # remove zero-padding from tt-res and append into tt_measure
        for i in range(len(Num_Of_dets)): # for different detectors
            self.tt_measure.append([tt_res[i][(index * Config.vec_size): (index * Config.vec_size + counts)].tolist() for index, counts in
                                    enumerate(counts_res[i])])
            self.tt_measure[i].sort()
        # create vector of tt's for each direction (north and south) and sort them
        self.tt_N_measure = sorted(sum(sum(self.tt_measure[:3],[]),[])) # unify detectors 1-3 and windows within detectors
        self.tt_S_measure = sorted(sum(sum(self.tt_measure[3:],[]),[])) # unify detectors 6-8 and windows within detectors


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
        self.num_of_det_reflections_per_seq_S = np.zeros(self.number_of_sprint_sequences)
        self.num_of_det_reflections_per_seq_N = np.zeros(self.number_of_sprint_sequences)
        self.num_of_det_transmissions_per_seq_S = np.zeros(self.number_of_sprint_sequences)
        self.num_of_det_transmissions_per_seq_N = np.zeros(self.number_of_sprint_sequences)

        self.num_of_SPRINT_reflections_per_seq_S = np.zeros([self.number_of_sprint_sequences, self.number_of_SPRINT_pulses_per_seq])
        self.num_of_SPRINT_reflections_per_seq_N = np.zeros([self.number_of_sprint_sequences, self.number_of_SPRINT_pulses_per_seq])
        self.num_of_SPRINT_transmissions_per_seq_S = np.zeros([self.number_of_sprint_sequences, self.number_of_SPRINT_pulses_per_seq])
        self.num_of_SPRINT_transmissions_per_seq_N = np.zeros([self.number_of_sprint_sequences, self.number_of_SPRINT_pulses_per_seq])

        end_of_det_pulse_in_seq = max(self.pulses_location_in_seq_N[num_of_det_pulses - 1][1],
                                      self.pulses_location_in_seq_S[num_of_det_pulses - 1][1])

        for element in self.tt_N_measure:
            tt_inseq = element % self.sprint_sequence_len
            if tt_inseq <= end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                self.num_of_det_reflections_per_seq_S[(element-1)//self.sprint_sequence_len] += self.filter_S[tt_inseq]
                self.num_of_det_transmissions_per_seq_N[(element-1)//self.sprint_sequence_len] += self.filter_N[tt_inseq]
            else:  # The part of the SPRINT pulses in the sequence
                SPRINT_pulse_num = (tt_inseq - end_of_det_pulse_in_seq) // (sprint_pulse_len + Config.num_between_zeros)
                if SPRINT_pulse_num < self.number_of_SPRINT_pulses_per_seq:
                    self.num_of_SPRINT_reflections_per_seq_S[(element-1) // self.sprint_sequence_len][SPRINT_pulse_num] += self.filter_S[tt_inseq]
                    self.num_of_SPRINT_transmissions_per_seq_N[(element-1) // self.sprint_sequence_len][SPRINT_pulse_num] += self.filter_N[tt_inseq]

        for element in self.tt_S_measure:
            tt_inseq = element % self.sprint_sequence_len
            if tt_inseq <= end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                self.num_of_det_reflections_per_seq_N[(element-1)//self.sprint_sequence_len] += self.filter_N[tt_inseq]
                self.num_of_det_transmissions_per_seq_S[(element-1)//self.sprint_sequence_len] += self.filter_S[tt_inseq]
            else:  # The part of the SPRINT pulses in the sequence
                SPRINT_pulse_num = (tt_inseq - end_of_det_pulse_in_seq) // (sprint_pulse_len + Config.num_between_zeros)
                if SPRINT_pulse_num < self.number_of_SPRINT_pulses_per_seq:
                    self.num_of_SPRINT_reflections_per_seq_N[(element-1) // self.sprint_sequence_len][SPRINT_pulse_num] += self.filter_N[tt_inseq]
                    self.num_of_SPRINT_transmissions_per_seq_S[(element-1) // self.sprint_sequence_len][SPRINT_pulse_num] += self.filter_S[tt_inseq]

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
        detection_pulse_range = np.zeros([self.number_of_sprint_sequences, num_of_det_pulses], dtype=int)
        sprint_pulse_range = np.zeros([self.number_of_sprint_sequences, num_of_sprint_pulses], dtype=int)

        #define ranges for tt_histogram (a vector that counts number of photons in each pulse)
        for i in range(self.number_of_sprint_sequences):
            detection_pulse_range[i] =\
                list(range(i*num_of_pulses, i*num_of_pulses+num_of_det_pulses))
            # sprint pulse range starts from last detection pulse andfo num_of_sprint_pulses
            sprint_pulse_range[i] = \
                list(range(int(detection_pulse_range[i][-1])+1, int(detection_pulse_range[i][-1])+1+num_of_sprint_pulses))
        # find transits and sprint events
        for j in range(self.number_of_sprint_sequences - 2):
            if \
                sum(self.tt_histogram_reflection[detection_pulse_range[j]]) >= detection_condition[0] and \
                        sum(self.tt_histogram_reflection[detection_pulse_range[j + 1]]) >= detection_condition[1]:
                    num_of_detected_atom += 1
                    transit_sequences_init_tt.append(j * self.sprint_sequence_len)
                    if self.tt_histogram_reflection[sprint_pulse_range[j][0]]:
                        num_of_sprints += 1
                        sprints_events.append(self.tt_histogram_reflection[sprint_pulse_range[j]])
        sprints_data = [sprints_events, num_of_sprints]
        atom_detect_data = [transit_sequences_init_tt, num_of_detected_atom]
        return atom_detect_data, sprints_data

    def find_transits_and_sprint_events_changed(self, cond=[2,2,2], minimum_number_of_seq_detected=2):
        '''
        find transits of atoms by searching events with detection_condition[0] photons before
        the sprint sequence and detection_condition[1] after.
        :param minimum_number_of_seq_detected:
        :param cond:
        :param detection_condition:
        :param num_of_det_pulses:
        :param num_of_sprint_pulses:
        :return:
        '''

        current_transit = []
        self.all_transits_seq_indx = []
        self.reflection_SPRINT_data_per_transit = []
        self.transmission_SPRINT_data_per_transit = []

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
                        self.reflection_SPRINT_data_per_transit = self.reflection_SPRINT_data_per_transit[:-1]
                        self.transmission_SPRINT_data_per_transit = self.transmission_SPRINT_data_per_transit[:-1]
                self.all_transits_seq_indx.append(current_transit)
                self.reflection_SPRINT_data_per_transit.append([self.num_of_SPRINT_reflections_per_seq[elem].tolist()
                                                                for elem in current_transit])
                self.transmission_SPRINT_data_per_transit.append([self.num_of_SPRINT_transmissions_per_seq[elem].tolist()
                                                                  for elem in current_transit])
                current_transit = []
        if len(current_transit) > 1:
            current_transit = current_transit[
                              :np.where(self.num_of_det_reflections_per_seq[current_transit] >= min(cond))[0][-1] + 1]
            if self.all_transits_seq_indx:
                if bool(set(current_transit) & set(self.all_transits_seq_indx[-1])):
                    current_transit = self.all_transits_seq_indx[-1] + current_transit[1:]
                    self.all_transits_seq_indx = self.all_transits_seq_indx[:-1]
                    self.reflection_SPRINT_data_per_transit = self.reflection_SPRINT_data_per_transit[:-1]
                    self.transmission_SPRINT_data_per_transit = self.transmission_SPRINT_data_per_transit[:-1]
            self.all_transits_seq_indx.append(current_transit)
            self.reflection_SPRINT_data_per_transit.append([self.num_of_SPRINT_reflections_per_seq[elem].tolist()
                                                            for elem in current_transit])
            self.transmission_SPRINT_data_per_transit.append([self.num_of_SPRINT_transmissions_per_seq[elem].tolist()
                                                              for elem in current_transit])

    def get_pulses_location_in_seq(self, delay, seq=Config.Sprint_Exp_Gaussian_samples_S, smearing = int(Config.num_between_zeros/2)):
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
        # seq_filter = np.append(np.zeros(delay), seq_filter[delay:])
        seq_filter = np.roll(seq_filter, delay)
        seq_indx = np.where(seq_filter > 0)[0]
        seq_filter_with_smearing = seq_filter
        pulses_loc = []
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


    def get_avg_num_of_photons_in_seq_pulses(self, seq, pulse_loc):
        avg_num_of_photons_in_seq_pulses = []
        try:
            real_number_of_seq = math.ceil(max(self.tt_S_measure)/len(Config.Sprint_Exp_Gaussian_samples_S))
            # print('Real number of seq = %d' %real_number_of_seq)
        except:
            real_number_of_seq = self.number_of_sprint_sequences
            # print('Max number of seq')
        for t in pulse_loc:
            # avg_num_of_photons_in_seq_pulses.append((sum(seq[t[0]:t[1]]) + seq[t[1]]) / (real_number_of_seq * 0.147))
            avg_num_of_photons_in_seq_pulses.append((sum(seq[t[0]:t[1]]) + seq[t[1]]) / (real_number_of_seq * 0.13))  # Sagnac configuiration efficiency 0.
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
        return (np.sum(pulse_loc, axis=1)/2).astype(int)


    # def get_avg_num_of_photons_in_seq_pulse(self, delay, seq=Config.Sprint_Exp_Gaussian_samples_S):
    #     seq_filter = (np.array(seq) > 0).astype(int)
    #     seq_filter = np.append(np.zeros(delay), seq_filter[delay:])
    #     seq_indx = np.where(seq_filter > 0)[0]
    #     self.avg_num_of_photons_in_seq_pulse = []
    #     self.pulse_loc = []
    #     start_indx = seq_indx[0]
    #     number_of_photons_in_current_pulse = self.tt_N_det_SPRINT_events_batch[seq_indx[0]] + self.tt_S_det_SPRINT_events_batch[seq_indx[0]]
    #     for i in range(1, len(seq_indx)):
    #         if seq_indx[i] - seq_indx[i-1] > 1:
    #             self.avg_num_of_photons_in_seq_pulse.append(number_of_photons_in_current_pulse)
    #             self.pulse_loc.append((start_indx, seq_indx[i-1]))
    #             start_indx = seq_indx[i]
    #             number_of_photons_in_current_pulse = 0
    #         number_of_photons_in_current_pulse += self.tt_N_det_SPRINT_events_batch[seq_indx[i]] + self.tt_S_det_SPRINT_events_batch[seq_indx[i]]
    #     self.avg_num_of_photons_in_seq_pulse.append(number_of_photons_in_current_pulse)
    #     self.pulse_loc.append((start_indx, seq_indx[-1]))


    def fold_tt_histogram(self, exp_sequence_len, delay_S,delay_N,delay_det_N=11,delay_det_S=6,delay_rf_N=15,delay_rf_S=13):
        # fold north and south
        self.folded_tt_S = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_N= np.zeros(exp_sequence_len, dtype=int)

        for x in [elem for elem in self.tt_S_measure]:
            self.folded_tt_S[x % exp_sequence_len] += 1
        for x in [elem for elem in self.tt_N_measure]:
            self.folded_tt_N[x % exp_sequence_len] += 1



        # a vector of ones at pulses inexes and zeros else, used to take only pulses location
        S_pulses_location = np.asarray(Config.Sprint_Exp_Gaussian_samples_S).astype(bool).astype(int)
        N_pulses_location = np.asarray(Config.Sprint_Exp_Gaussian_samples_N).astype(bool).astype(int)
        self.S_pulses_loc_delayed = np.roll(S_pulses_location, delay_rf_S)
        self.N_pulses_loc_delayed = np.roll(N_pulses_location, delay_rf_N)


        self.folded_tt_S_delayed = np.roll(self.folded_tt_S, Config.num_init_zeros_S - delay_S+delay_det_S)
        self.folded_tt_N_delayed = np.roll(self.folded_tt_N, Config.num_init_zeros_N - delay_N+delay_det_N)

        folded_transmission = self.folded_tt_S_delayed*self.S_pulses_loc_delayed + self.folded_tt_N_delayed *self.N_pulses_loc_delayed
        folded_reflection = self.N_pulses_loc_delayed*self.folded_tt_S_delayed + self.S_pulses_loc_delayed*self.folded_tt_N_delayed
        return folded_transmission, folded_reflection

    def plt_pulses_roll(self, delay_det_N, delay_det_S, delay_rf_N, delay_rf_S):
        plt.figure()
        plt.plot(np.roll(self.folded_tt_S_delayed, delay_det_S), label='folded_S')
        plt.plot(np.roll(self.folded_tt_N_delayed, delay_det_N), label='folded_N')
        plt.plot(np.roll(self.S_pulses_loc_delayed, delay_rf_S), label='sent_S')
        plt.plot(np.roll(self.N_pulses_loc_delayed, delay_rf_N), label='sent_N')
        plt.legend()
        plt.title([delay_det_N, delay_det_S, delay_rf_N, delay_rf_S])

    def save_tt_to_batch(self,Num_Of_dets,N):
        for i in range(len(Num_Of_dets)):
            self.tt_measure_batch[i] = self.tt_measure_batch[i][-(N - 1):] + [self.tt_measure[i]]
        self.tt_S_measure_batch = self.tt_S_measure_batch[-(N - 1):] + [self.tt_S_measure]
        self.tt_N_measure_batch = self.tt_N_measure_batch[-(N - 1):] + [self.tt_N_measure]

        # self.transit_sequences_batch += self.transit_sequences


    def plot_sprint_figures(self,ax,Num_Of_dets):

        ax[0].clear()
        ax[1].clear()
        ax[2].clear()
        ax[3].clear()
        ax[4].clear()
        ax[5].clear()
        #
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        textstr_total_reflections = 'Total reflections per cycle "N"=%d \n' % (sum(self.num_of_det_reflections_per_seq_N),)\
                                    + 'Total reflections per cycle "S"=%d' % (sum(self.num_of_det_reflections_per_seq_S),)
        textstr_avg_reflections = r'Average reflections per cycle = %.2f' % (sum(self.num_of_det_reflections_per_seq_accumulated/self.Counter),)
        #     textstr_detuned = r'$Probe_N = %.3f$' % ((sum(tt_S_binning_detuned) * 1000) / (self.M_window / 2),) + '[MPhotons/sec]\n' \
        #                + '$\overline{Probe}_S = %.3f$' % ((np.mean([sum(x) for x in tt_S_binning_detuned_batch]) * 1000)
        #                                                   / (self.M_window / 2),) + '[MPhotons/sec]'
        #     textstr_resonance = r'$Probe_S = %.3f$' % ((sum(tt_S_binning_resonance) * 1000) / (self.M_window / 2),) + '[MPhotons/sec]\n' \
        #                 + '$\overline{Probe}_S = %.3f$' % ((np.mean([sum(x) for x in tt_S_binning_resonance_batch]) * 1000)
        #                                                    / (self.M_window / 2),) + '[MPhotons/sec]'
        #     textstr_FLR = r'$\overline{FLR}_{MAX} = %.1f$' % (np.mean(FLR_measurement) * 1e5,) + r'$\times 10^{-5}$'
        # textstr_No_transits = 'NO TRANSITS YET!!!'

        ax[0].plot(self.folded_tt_N_batch, label='"N" detectors')
        ax[0].plot(self.folded_tt_S_batch, label='"S" detectors')
        ax[0].plot((self.filter_N) * max(self.folded_tt_N_batch + self.folded_tt_S_batch), '--k', label='Filter "N"')
        for i in range(len(self.Num_of_photons_txt_box_y_loc)):
            ax[0].text(self.Num_of_photons_txt_box_x_loc.tolist()[i], self.Num_of_photons_txt_box_y_loc[i],
                     '%.2f' % self.avg_num_of_photons_per_pulse[i],
                     horizontalalignment='center', fontsize=12, fontweight='bold', family=['Comic Sans MS'])
        ax[0].set_title('binned timetags from all detectors folded (Averaged)', fontweight="bold")
        ax[0].legend(loc='upper right')

        ax[1].plot(self.folded_tt_N, label='"N" detectors')
        ax[1].plot(self.folded_tt_S, label='"S" detectors')
        ax[1].plot((self.filter_S) * max(self.folded_tt_N + self.folded_tt_S), '--k', label='Filter "S"')
        for i in range(len(self.Num_of_photons_txt_box_y_loc_live)):
            ax[1].text(self.Num_of_photons_txt_box_x_loc.tolist()[i], self.Num_of_photons_txt_box_y_loc_live[i],
                       '%.2f' % self.avg_num_of_photons_per_pulse_live[i],
                       horizontalalignment='center', fontsize=12, fontweight='bold', family=['Comic Sans MS'])
        ax[1].set_title('binned timetags from all detectors folded (Live)', fontweight="bold")
        ax[1].legend(loc='upper right')

        # ax[2].plot(self.folded_transmission, label='folded transmission')
        # ax[2].plot(self.folded_reflection, label='folded reflection')
        # ax[2].set_title('folded reflection and transmission(Live)', fontweight="bold")
        # ax[2].legend(loc='upper right')

        ax[2].plot(self.num_of_det_reflections_per_seq, label='Num of reflections per sequence (Live)')
        ax[2].set_title('Num of reflections per sequence (live)', fontweight="bold")
        ax[2].text(0.1, 0.9*max(self.num_of_det_reflections_per_seq), textstr_total_reflections, fontsize=14,
                   verticalalignment='top', bbox=props)

        ax[3].plot(self.num_of_det_reflections_per_seq_accumulated/self.Counter, label='Num of reflections per sequence (Live)')
        ax[3].set_title('Num of reflections per sequence (Average)', fontweight="bold")
        ax[3].text(0.1, 0.9*max(self.num_of_det_reflections_per_seq_accumulated/self.Counter), textstr_avg_reflections, fontsize=14,
                   verticalalignment='top', bbox=props)


        # ax[4].plot(np.histogram(self.transit_sequences,bins=self.number_of_sprint_sequences,
        #                         range=(0, self.M_window))[0], label='transits(Live)')
        # ax[4].set_title('transit histogram(live)', fontweight="bold")
        # ax[4].legend(loc='upper right')
        # ax[4].text(0.5, 0.5, 'number of transits='+str(self.atom_detect_data[1]))
        #
        # ax[5].plot(np.histogram(self.transit_sequences_batch,bins=self.number_of_sprint_sequences
        #                         ,range=(0,self.M_window))[0], label='transits(Live)')
        # ax[5].set_title('transit histogram(accumulated)', fontweight="bold")
        # ax[5].legend(loc='upper right')

        # if len(self.all_transits_seq_indx_batch) > 0:
        if self.number_of_transits_live:
            textstr_transit_counts = r'$N_{Transits} = %s $' % (self.number_of_transits_live,) + r'$[Counts]$'
        else:
            textstr_transit_counts = r'$N_{Transits} = %s $' % (0,) + r'$[Counts]$'
        textstr_transit_event_counter = r'$N_{Transits Total} = %s $' % (self.number_of_transits_total,) + r'$[Counts]$'

        ax[4].plot(range(self.number_of_sprint_sequences), self.seq_transit_events_live, label='Current transit events')
        ax[4].set(xlabel='Sequence [#]', ylabel='Counts [Photons]')
        ax[4].set_title('Transits per sequence (live)', fontweight="bold")
        ax[4].text(0.05, 0.95, textstr_transit_counts, transform=ax[4].transAxes, fontsize=14,
                   verticalalignment='top', bbox=props)

        ax[5].plot(range(self.number_of_sprint_sequences), self.seq_transit_events_batched, label='Transit events histogram')
        ax[5].set(xlabel='Sequence [#]', ylabel='Counts [Photons]')
        ax[5].set_title('Transits per sequence (accumulated)', fontweight="bold")
        ax[5].text(0.05, 0.95, textstr_transit_event_counter, transform=ax[5].transAxes, fontsize=14,
                   verticalalignment='top', bbox=props)
        # else:
        #     ax[4].plot(range(self.number_of_sprint_sequences), self.seq_transit_events_live, label='Current transit events')
        #     ax[4].set(xlabel='Sequence [#]', ylabel='Counts [Photons]')
        #     ax[4].text(0.25, 0.5, textstr_No_transits, transform=ax[4].transAxes, fontsize=24,
        #                verticalalignment='center', bbox=props)
        #
        #     ax[5].plot(range(self.number_of_sprint_sequences), self.seq_transit_events_batched, label='Transit events histogram')
        #     ax[5].set(xlabel='Sequence [#]', ylabel='Counts [Photons]')
        #     ax[5].text(0.25, 0.5, textstr_No_transits, transform=ax[5].transAxes, fontsize=24,
        #                verticalalignment='center', bbox=props)

        # plt.show()
        plt.pause(1.0)
        ###########################################################################################################

    def init_params_for_save_sprint(self,sprint_sequence_len,Num_Of_dets):
        # define empty variables
        self.sprint_sequence_len=sprint_sequence_len
        # self.number_of_sprint_sequences = self.M_window // (self.sprint_sequence_len)
        self.number_of_sprint_sequences = math.ceil(experiment.M_window / experiment.sprint_sequence_len)
        self.number_of_SPRINT_pulses_per_seq = len(Config.sprint_pulse_amp_S)
        time_bins = np.linspace(0, self.M_window, self.number_of_sprint_sequences)

        tt_S_SPRINT_events = np.zeros(self.number_of_sprint_sequences)
        self.tt_S_SPRINT_events_batch = np.zeros(self.sprint_sequence_len)

        #
        FLR_measurement = []
        Exp_timestr_batch = []
        self.tt_measure = []
        self.tt_S_measure = []
        self.tt_measure_batch = [[], [], [], [], [], []]
        self.tt_S_measure_batch = []
        self.tt_N_measure_batch = []
        self.transit_sequences_batch = []
        self.all_transits_seq_indx_batch = []

        self.folded_transmission = np.zeros(len(Config.Sprint_Exp_Gaussian_samples_S))
        self.folded_reflection = np.zeros(len(Config.Sprint_Exp_Gaussian_samples_S))

        self.tt_S_binning = np.zeros(self.number_of_sprint_sequences + 1)
        self.seq_transit_events_live = np.zeros(self.number_of_sprint_sequences)
        self.seq_transit_events_batched = np.zeros(self.number_of_sprint_sequences)
        self.tt_S_SPRINT_events = np.zeros(self.sprint_sequence_len)
        self.tt_S_SPRINT_events_batch = np.zeros(self.sprint_sequence_len)
        self.Single_det_foldeded = np.zeros((len(Num_Of_dets), self.sprint_sequence_len))
        self.single_det_folded_accumulated = np.zeros((len(Num_Of_dets), self.sprint_sequence_len))
        self.num_of_det_reflections_per_seq_accumulated = np.zeros(self.number_of_sprint_sequences)

    def Save_SNSPDs_Sprint_Measurement_with_tt(self, N, sprint_sequence_len, preComment, lock_err_threshold, transit_condition,
                                               max_probe_counts, filter_delay, reflection_threshold, reflection_threshold_time):
        """
        Function for analyzing,saving and displaying data from sprint experiment.
        :param N: Number of maximum experiments (free throws) saved and displayed.
                 program we are looking for transits of atoms next to the toroid and record them.
        :param sprint_sequence_len: the number of sprint sequences (detection and sprint pulses combination)
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
        Num_Of_dets = [1, 2, 3, 6, 7, 8]
        delay_in_detection_N = 30 # choose the correct delay in samples to the first detection pulse # TODO: 40?
        delay_in_detection_S = 20 # choose the correct delay in samples to the first detection pulse # TODO: 40?
        det_pulse_len = Config.det_pulse_len+Config.num_between_zeros
        sprint_pulse_len = Config.sprint_pulse_len+Config.num_between_zeros

        # initialize parameters - set zeros vectors and empty lists
        self.init_params_for_save_sprint(sprint_sequence_len, Num_Of_dets)

        # get pulses location south and north
        self.pulses_location_in_seq_S, self.filter_S = self.get_pulses_location_in_seq(filter_delay[0],
                                                                                       Config.Sprint_Exp_Gaussian_samples_S,
                                                                                       smearing=5)  # smearing=int(Config.num_between_zeros/2))
        self.pulses_location_in_seq_N, self.filter_N = self.get_pulses_location_in_seq(filter_delay[1],
                                                                                       Config.Sprint_Exp_Gaussian_samples_N,
                                                                                       smearing=5)  # smearing=int(Config.num_between_zeros/2))
        self.Num_of_photons_txt_box_x_loc = np.concatenate((self.num_of_photons_txt_box_loc(self.pulses_location_in_seq_S),
                                                           self.num_of_photons_txt_box_loc(self.pulses_location_in_seq_N)))

        ## Listen for keyboard
        listener = keyboard.Listener(on_press=self.on_key_press)
        listener.start()  # start to listen on a separate thread
        self.keyPress = None
        print('\033[94m' + 'Press ESC to stop measurement.' + '\033[0m')  # print blue

        reps = 1  # a counter, number of repeats actually made.

        ####     get tt and counts from OPX to python   #####
        Counts_handle, tt_handle, FLR_handle = self.get_handles_from_OPX_Server(Num_Of_dets)

        ## take data only if
        start = True
        # take threshold from npz ( error from resonator lock PID)
        # lock_err = np.abs(np.load(
        #     'U:\Lab_2021-2022\Experiment_results\Sprint\Locking_PID_Error\locking_err.npy', allow_pickle=True)) # the error of locking the resontor to Rb line
        lock_err = lock_err_threshold/2
        sum_for_threshold = 0

        # Place holders for results # TODO: ask dor - is it works as we expect?
        # while (number of photons * 10^-6 [Mcounts] / Measuring time [nsec] * 10^-9 [sec/nsec])  > counts_threshold [Mcounts /sec]:
        while ((lock_err > lock_err_threshold) and (sum_for_threshold < reflection_threshold)) or start:
            # lock_err = np.abs(np.load(
            #     'U:\Lab_2021-2022\Experiment_results\Sprint\Locking_PID_Error\locking_err.npy'))  # the error of locking the resontor to Rb line
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
            self.get_tt_from_handles(Num_Of_dets, Counts_handle, tt_handle, FLR_handle)
            self.divide_tt_to_reflection_trans(sprint_pulse_len, len([x for x in Config.det_pulse_amp_S if x > 0]))
            self.num_of_det_reflections_per_seq = self.num_of_det_reflections_per_seq_S \
                                                  + self.num_of_det_reflections_per_seq_N
            self.num_of_SPRINT_reflections_per_seq = self.num_of_SPRINT_reflections_per_seq_S \
                                                     + self.num_of_SPRINT_reflections_per_seq_N
            self.num_of_SPRINT_transmissions_per_seq = self.num_of_SPRINT_transmissions_per_seq_S \
                                                       + self.num_of_SPRINT_transmissions_per_seq_N
            sum_for_threshold = sum(self.num_of_det_reflections_per_seq[-int(reflection_threshold_time//len(Config.Sprint_Exp_Gaussian_samples_S)):])  # summing over the reflection from detection pulses of each sequence corresponding the the reflection_threshold_time
        ####    end get tt and counts from OPX to python   #####

        self.num_of_det_reflections_per_seq_accumulated += self.num_of_det_reflections_per_seq_S \
                                                           + self.num_of_det_reflections_per_seq_N

        # divide south and north into reflection and transmission
        self.tt_histogram_transmission, self.tt_histogram_reflection = \
        self.divide_to_reflection_trans(det_pulse_len=det_pulse_len, sprint_pulse_len=sprint_pulse_len,
                                        sprint_sequence_delay_S=delay_in_detection_S, sprint_sequence_delay_N=delay_in_detection_N,
                                        num_of_det_pulses=len(Config.det_pulse_amp_S),
                                        num_of_sprint_pulses=len(Config.sprint_pulse_amp_S),
                                        num_of_sprint_sequences=self.number_of_sprint_sequences)

        # TODO: needed?
        self.tt_S_SPRINT_events = np.zeros(sprint_sequence_len)
        self.tt_S_SPRINT_events_batch = np.zeros(sprint_sequence_len)
        self.tt_Single_det_SPRINT_events = np.zeros((len(Num_Of_dets), sprint_sequence_len))
        self.tt_Single_det_SPRINT_events_batch = np.zeros((len(Num_Of_dets), sprint_sequence_len))
        self.tt_N_det_SPRINT_events_batch = np.zeros(sprint_sequence_len)
        self.tt_S_det_SPRINT_events_batch = np.zeros(sprint_sequence_len)

        # fold reflections and transmission
        self.folded_transmission, self.folded_reflection = \
            self.fold_tt_histogram(exp_sequence_len=self.sprint_sequence_len, delay_S=delay_in_detection_S, delay_N=delay_in_detection_N)
        self.folded_transmission_accumulated,self.folded_reflection_accumulated = \
            self.folded_transmission, self.folded_reflection

        # fold data from different detectors: # TODO: is needed? @Dor
        for i in range(len(Num_Of_dets)):
            # for x in [elem for elem in self.tt_S_measure if elem < self.M_window]: - for debugging assaf
            for x in [elem for elem in self.tt_measure[i][-1]]:
                self.Single_det_foldeded[i][x % self.sprint_sequence_len] += 1
                self.single_det_folded_accumulated[i][x % self.sprint_sequence_len] += 1

        # Batch folded tt "N" and "S"
        self.folded_tt_N_batch = self.folded_tt_N
        self.folded_tt_S_batch = self.folded_tt_S

        # get the average number of photons in detection pulse
        self.avg_num_of_photons_per_pulse_S_live = self.get_avg_num_of_photons_in_seq_pulses(self.folded_tt_S, self.pulses_location_in_seq_S)
        self.avg_num_of_photons_per_pulse_N_live = self.get_avg_num_of_photons_in_seq_pulses(self.folded_tt_N, self.pulses_location_in_seq_N)
        self.avg_num_of_photons_per_pulse_live = self.avg_num_of_photons_per_pulse_S_live + self.avg_num_of_photons_per_pulse_N_live
        self.avg_num_of_photons_per_pulse = self.avg_num_of_photons_per_pulse_live

        # get box location on the y-axis:
        self.max_value_per_pulse_S_live = self.get_max_value_in_seq_pulses(self.folded_tt_S, self.pulses_location_in_seq_S)
        self.max_value_per_pulse_N_live = self.get_max_value_in_seq_pulses(self.folded_tt_N, self.pulses_location_in_seq_N)
        self.Num_of_photons_txt_box_y_loc_live = self.max_value_per_pulse_S_live + self.max_value_per_pulse_N_live
        self.Num_of_photons_txt_box_y_loc = self.Num_of_photons_txt_box_y_loc_live

        # self.get_avg_num_of_photons_in_det_pulse(det_pulse_len=(Config.det_pulse_len+Config.num_between_zeros),
        #                                          sprint_sequence_delay=delay_in_detection_S, num_of_det_pulses=len(Config.det_pulse_amp_S),
        #                                          num_of_sprint_sequences=self.number_of_sprint_sequences)
        # print('average number of photons in detection pulses is:', self.avg_num_of_photons_in_det_pulse)

        ## record time
        timest = time.strftime("%Y%m%d-%H%M%S")
        datest = time.strftime("%Y%m%d")
        FLR_measurement = []
        Exp_timestr_batch = []

        FLR_measurement = FLR_measurement[-(N - 1):] + [self.FLR_res.tolist()]
        Exp_timestr_batch = Exp_timestr_batch[-(N - 1):] + [timest]

        # ### Find transits and build histogram:  ###
        # [self.atom_detect_data, self.sprints_data] = \
        #     self.find_transits_and_sprint_events(detection_condition=transit_condition,
        #                                          num_of_det_pulses=len(Config.det_pulse_amp_S),
        #                                          num_of_sprint_pulses=len(Config.sprint_pulse_amp_S))
        # self.transit_sequences = self.atom_detect_data[0]

        ### Dor version ###
        self.find_transits_and_sprint_events_changed(cond=transit_condition, minimum_number_of_seq_detected=2)
        self.seq_transit_events_live[[elem for vec in self.all_transits_seq_indx for elem in vec]] += 1
        self.seq_transit_events_batched[[elem for vec in self.all_transits_seq_indx for elem in vec]] += 1
        self.all_transits_seq_indx_batch = self.all_transits_seq_indx_batch[-(N-1):] + [self.all_transits_seq_indx]
        self.number_of_transits_live = len(self.all_transits_seq_indx)
        self.number_of_transits_total = len([vec for lst in self.all_transits_seq_indx_batch for vec in lst])
        ###################

        self.save_tt_to_batch(Num_Of_dets, N)
        self.Counter = 1
        #
        # create figures template
        ax1 = plt.subplot2grid((3, 2), (0, 0), colspan=1, rowspan=1)
        ax2 = plt.subplot2grid((3, 2), (0, 1), colspan=1, rowspan=1)
        ax3 = plt.subplot2grid((3, 2), (1, 0), colspan=1, rowspan=1)
        ax4 = plt.subplot2grid((3, 2), (1, 1), colspan=1, rowspan=1)
        ax5 = plt.subplot2grid((3, 2), (2, 0), colspan=1, rowspan=1)
        ax6 = plt.subplot2grid((3, 2), (2, 1), colspan=1, rowspan=1)
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
            ######################################## PLOT!!! ###########################################################
            ax = [ax1, ax2, ax3, ax4, ax5, ax6]
            self.plot_sprint_figures(ax, Num_Of_dets)
            ############################################################################################################
            while True:
                # record time:
                timest = time.strftime("%Y%m%d-%H%M%S") # TODO: is it needed? already writen above..
                datest = time.strftime("%Y%m%d")

                self.get_tt_from_handles(Num_Of_dets, Counts_handle, tt_handle, FLR_handle)

                if self.tt_S_measure != self.tt_S_measure_batch[-1]:
                    break
            # assaf - if x=self.M_window the index is out of range so i added 1
            # lock_err = np.abs(np.load(
            #     'U:\Lab_2021-2022\Experiment_results\Sprint\Locking_PID_Error\locking_err.npy'))  # the error of locking the resontor to Rb line
            self.divide_tt_to_reflection_trans(sprint_pulse_len, len([x for x in Config.det_pulse_amp_S if x > 0]))
            self.num_of_det_reflections_per_seq = self.num_of_det_reflections_per_seq_S \
                                                  + self.num_of_det_reflections_per_seq_N
            self.num_of_SPRINT_reflections_per_seq = self.num_of_SPRINT_reflections_per_seq_N \
                                                     + self.num_of_SPRINT_reflections_per_seq_S
            self.num_of_SPRINT_transmissions_per_seq = self.num_of_SPRINT_transmissions_per_seq_N \
                                                       + self.num_of_SPRINT_transmissions_per_seq_S
            sum_for_threshold = sum(self.num_of_det_reflections_per_seq[-int(reflection_threshold_time//len(Config.Sprint_Exp_Gaussian_samples_S)):])  # summing over the reflection from detection pulses of each sequence corresponding the the reflection_threshold_time

            if (lock_err < lock_err_threshold) and (sum_for_threshold < reflection_threshold):
                print('Sum of reflections: %d' %sum_for_threshold)
                self.num_of_det_reflections_per_seq_accumulated += self.num_of_det_reflections_per_seq_S \
                                                                   + self.num_of_det_reflections_per_seq_N

                self.Single_det_foldeded = np.zeros((len(Num_Of_dets), self.sprint_sequence_len))
                self.seq_transit_events_live = np.zeros(self.number_of_sprint_sequences)

                # divide south and north into reflection and transmission
                self.tt_histogram_transmission, self.tt_histogram_reflection = \
                    self.divide_to_reflection_trans(det_pulse_len=det_pulse_len, sprint_pulse_len=sprint_pulse_len,
                                                    sprint_sequence_delay_S=delay_in_detection_S,
                                                    sprint_sequence_delay_N=delay_in_detection_N,
                                                    num_of_det_pulses=len(Config.det_pulse_amp_S),
                                                    num_of_sprint_pulses=len(Config.sprint_pulse_amp_S),
                                                    num_of_sprint_sequences=self.number_of_sprint_sequences)

                # fold reflections and transmission
                self.folded_transmission, self.folded_reflection = \
                    self.fold_tt_histogram(exp_sequence_len=self.sprint_sequence_len, delay_S=delay_in_detection_S,
                                           delay_N=delay_in_detection_N)
                self.folded_transmission_accumulated += self.folded_transmission
                self.folded_reflection_accumulated += self.folded_reflection

                ### Find transits and build histogram:  ###
                # [self.atom_detect_data, self.sprints_data] = \
                #     self.find_transits_and_sprint_events(detection_condition=transit_condition,
                #                                          num_of_det_pulses=len(Config.det_pulse_amp_S),
                #                                          num_of_sprint_pulses=len(Config.sprint_pulse_amp_S))
                # self.transit_sequences = self.atom_detect_data[0]

                ### Dor version ###
                self.find_transits_and_sprint_events_changed(cond=transit_condition, minimum_number_of_seq_detected=2)
                self.seq_transit_events_live[[vec for elem in self.all_transits_seq_indx for vec in elem]] += 1
                self.seq_transit_events_batched[[vec for elem in self.all_transits_seq_indx for vec in elem]] += 1
                self.all_transits_seq_indx_batch = self.all_transits_seq_indx_batch[-(N - 1):] + [
                    self.all_transits_seq_indx]
                self.number_of_transits_live = len(self.all_transits_seq_indx)
                self.number_of_transits_total = len([vec for lst in self.all_transits_seq_indx_batch for vec in lst])
                ###################

                # Batch folded tt "N" and "S"
                self.folded_tt_N_batch = (self.folded_tt_N_batch * (self.Counter - 1) + self.folded_tt_N) / self.Counter
                self.folded_tt_S_batch = (self.folded_tt_S_batch * (self.Counter - 1) + self.folded_tt_S) / self.Counter

                # get the average number of photons in detection pulse
                self.avg_num_of_photons_per_pulse_S_live = self.get_avg_num_of_photons_in_seq_pulses(self.folded_tt_S, self.pulses_location_in_seq_S)
                self.avg_num_of_photons_per_pulse_N_live = self.get_avg_num_of_photons_in_seq_pulses(self.folded_tt_N, self.pulses_location_in_seq_N)
                self.avg_num_of_photons_per_pulse_live = self.avg_num_of_photons_per_pulse_S_live + self.avg_num_of_photons_per_pulse_N_live
                self.avg_num_of_photons_per_pulse_S = self.get_avg_num_of_photons_in_seq_pulses(self.folded_tt_S_batch, self.pulses_location_in_seq_S)
                self.avg_num_of_photons_per_pulse_N = self.get_avg_num_of_photons_in_seq_pulses(self.folded_tt_N_batch, self.pulses_location_in_seq_N)
                self.avg_num_of_photons_per_pulse = self.avg_num_of_photons_per_pulse_S + self.avg_num_of_photons_per_pulse_N

                # get box location on the y-axis:
                self.max_value_per_pulse_S_live = self.get_max_value_in_seq_pulses(self.folded_tt_S, self.pulses_location_in_seq_S)
                self.max_value_per_pulse_N_live = self.get_max_value_in_seq_pulses(self.folded_tt_N, self.pulses_location_in_seq_N)
                self.Num_of_photons_txt_box_y_loc_live = self.max_value_per_pulse_S_live + self.max_value_per_pulse_N_live
                self.max_value_per_pulse_S = self.get_max_value_in_seq_pulses(self.folded_tt_S_batch, self.pulses_location_in_seq_S)
                self.max_value_per_pulse_N = self.get_max_value_in_seq_pulses(self.folded_tt_N_batch, self.pulses_location_in_seq_N)
                self.Num_of_photons_txt_box_y_loc = self.max_value_per_pulse_S + self.max_value_per_pulse_N

                # fold for different detectors: # TODO: delete after everything works
                for i in range(len(Num_Of_dets)):
                    # for x in [elem for elem in self.tt_S_measure if elem < self.M_window]:
                    for x in [elem for elem in self.tt_measure[i][-1]]:
                        self.Single_det_foldeded[i][x % self.sprint_sequence_len] += 1
                        self.single_det_folded_accumulated[i][x % self.sprint_sequence_len] += 1

                FLR_measurement = FLR_measurement[-(N - 1):] + [self.FLR_res.tolist()]
                Exp_timestr_batch = Exp_timestr_batch[-(N - 1):] + [timest]
                if self.Counter < N:
                    self.Counter += 1
                print(timest, self.Counter)

                self.save_tt_to_batch(Num_Of_dets, N)
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
        # Counter_str = (self.Counter)
        filename_N_tt = f'North_timetags.npz'
        filename_Det_tt = []
        for i in Num_Of_dets:
            filename_Det_tt.append(f'Det'+str(i)+f'_timetags.npz')
        filename_S_tt = f'South_timetags.npz'
        filename_N_tt = f'North_timetags.npz'
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
        for i in range(len(Num_Of_dets)):
            if len(self.tt_measure_batch[i]) > 0:
                np.savez(dirname_S + filename_Det_tt[i], self.tt_measure_batch[i])
        if len(self.tt_S_measure_batch) > 0:
            np.savez(dirname_S + filename_S_tt, self.tt_S_measure_batch)
        if len(self.tt_N_measure_batch) > 0:
            np.savez(dirname_N + filename_N_tt, self.tt_N_measure_batch)
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

    def Start_Sprint_Exp_with_tt(self, N=100, sprint_sequence_len=int(len(Config.Sprint_Exp_Gaussian_samples_S)),
                                 transit_condition=[2,2], preComment=None, lock_err_threshold=0.01, filter_delay=[-7,2],
                                 reflection_threshold=100, reflection_threshold_time=1e6):
        # Max_probe_counts = self.Get_Max_Probe_counts(3)  # return the average maximum probe counts of 3 cycles.
        Max_probe_counts = None  # return the average maximum probe counts of 3 cycles.
        self.SPRINT_Exp_switch(True)
        self.update_parameters()
        self.Save_SNSPDs_Sprint_Measurement_with_tt(N, sprint_sequence_len, preComment,lock_err_threshold,
                                                    transit_condition, Max_probe_counts, filter_delay,
                                                    reflection_threshold, reflection_threshold_time)

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
        # experiment.Start_Sprint_Exp_with_tt(N=1000, transit_condition=[2,2,2], preComment='test', filter_delay=[0,0], reflection_threshold=35)
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


