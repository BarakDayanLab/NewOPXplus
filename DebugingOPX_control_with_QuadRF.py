# from Config import config
import Config_with_SNSPDs_and_QuadRF as Config
from Config_Table import Initial_Values, Phases_Names  # , Values_Factor
# from quadRFMOTController import QuadRFMOTController

from qm.QuantumMachinesManager import QuantumMachinesManager
from qm.qua import *
from qm import generate_qua_script
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
                "FLR_detection", "Measurement", "AOM_2-2/3'", "AOM_2-2/3'_detuned", "AOM_2-3'_for_interference"]  #, "PULSER_N", "PULSER_S", "AOM_LO")


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
          "AOM_2-2/3'", "AOM_2-3'_for_interference", "FLR_detection", "Measurement")  #, "PULSER_N", "PULSER_S", "AOM_LO")

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
          "AOM_2-2/3'", "AOM_2-3'_for_interference", "FLR_detection", "Measurement")  #, "PULSER_N", "PULSER_S", "AOM_LO")

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
          "AOM_2-2/3'", "AOM_2-2/3'_detuned", "FLR_detection", "Measurement", "Dig_detectors")  #, "PULSER_N", "PULSER_S")

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
          "AOM_2-2/3'", "AOM_2-2/3'_detuned", "FLR_detection", "Measurement", "Dig_detectors")  #, "PULSER_N", "PULSER_S")

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
          "AOM_2-2/3'_detuned", "Measurement")  #, "PULSER_N", "PULSER_S")

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


        # MOT variables
        MOT_Repetitions = declare(int, value=obj.Exp_Values['MOT_rep'])
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
        OD_Delay = declare(int, value=int(obj.OD_delay * 1e6 / 4))
        shutter_open_time = declare(int, value=int(obj.Shutter_open_time * 1e6 / 4))
        Rep = declare(int, value=obj.rep)

        # MW spectroscopy variables:
        MW_freq = declare(int, value=obj.MW_start_frequency)
        pulse_length_MW = declare(int, value=(obj.Pulse_Length_MW * int(1e3 / 4)))
        pulse_length_OD = declare(int, value=(obj.Pulse_Length_OD * int(1e3 / 4)))
        Repetitions = declare(int, value=1)
        Delta_f = declare(int, value=0)

        # Stream processin
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
                                    Measure(Pulse_1_duration)  # This triggers camera (Control 7)
                                    align(*all_elements, "AOM_2-2/3'")
                        align(*all_elements, "AOM_2-2/3'")

            with if_(~(AntiHelmholtz_ON | Transits_Exp_ON)):
                assign(AntiHelmholtz_ON, True)

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
            FLR_st.save('FLR_measure')

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
    sourceFile = open('debug.py', 'w')
    print(generate_qua_script(experiment.opx_control_prog, Config.config), file=sourceFile)
    sourceFile.close()


