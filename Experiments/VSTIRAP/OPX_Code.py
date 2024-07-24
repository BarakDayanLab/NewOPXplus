from Experiments.VSTIRAP import VSTIRAP_Config_Experiment as Config
from Experiments.Enums.Phases import Phases
from Experiments.Enums.IOParameters import IOParameters as IOP
from qm.qua import *
from Utilities.OPX_Utils import OPX_Utils

# TODO: Move this elsewhere...
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
        # play("OD_FS"*amp(0.2), "AOM_2-2/3'") #for when we want the 2-2/3' to play doring the MOT
        # play("Const_open" * amp(Config.AOM_Late_Attenuation_From_Const), "AOM_Late")
        # play("Const_open","PULSER_N")
        # play("Const_open","PULSER_S")

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

def Pulse_with_prep(total_pulse_duration, prep_duration, zero_pulse_duration, plus_pulse_duration, minus_pulse_duration):
    """
    This pulse id divided to two parts:
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
          "PULSER_VSTIRAP_1_1", "PULSER_VSTIRAP_1_0",
          "Measurement", "AOM_Early", "AOM_Late", "PULSER_N", "PULSER_S") # , "Dig_detectors")

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

# From one direction
def MZ_balancing(m_time, m_window, shutter_open_time, phase_rep, points_for_sum,
                 counts_st_B, counts_st_D, phase_correction_st, phase_for_min_st, tt_st_B, tt_st_D):

    counts_B = declare(int)
    # tt_vec_B2 = declare(int, size=Config.vec_size)
    counts_B_sum = declare(int, value=0)
    counts_D = declare(int)
    # tt_vec_D2 = declare(int, size=Config.vec_size)
    counts_D_sum = declare(int, value=0)
    tot_counts_MZ_S = declare(int, value=0)

    OPX_Utils.assign_variables_to_element("AOM_Early", counts_B_sum)
    OPX_Utils.assign_variables_to_element("AOM_Late", counts_D_sum)
    OPX_Utils.assign_variables_to_element("PULSER_N", tot_counts_MZ_S)

    t = declare(int)
    i = declare(fixed)
    j = declare(int)
    k = declare(int)

    phase_scan = declare(fixed, value=[1, 0.25])
    # phase_scan = declare(fixed, value=[1, 0.1])
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
    assign(tot_counts_MZ_S, 0)
    assign(phase_correction_min, 0.0)
    assign(phase_correction_start, 0.0)

    # align("AOM_Early", "AOM_Late", "PULSER_N", "PULSER_S", "PULSER_ANCILLA")

    align("AOM_Early", "AOM_Late", "PULSER_N", "PULSER_S")

    with for_each_(i, phase_scan):
        assign(phase_correction_per_scan, phase_correction * i)
        assign(phase_correction_start, phase_correction_min - i * 0.5)
        assign(max_counts_B, 0)
        reset_frame("AOM_Early", "AOM_Late") #, "PULSER_S", "PULSER_N")
        frame_rotation_2pi(phase_correction_start, "AOM_Early")
        with for_(j, 0, j < phase_rep, j + 1):
            assign(counts_B_sum, 0)
            assign(counts_D_sum, 0)
            with for_(k, 1, k <= points_for_sum, k + 1):
                play("MZ_balancing_pulses_N", "PULSER_N")
                play("MZ_balancing_pulses_N", "PULSER_S")
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
            assign(tot_counts_MZ_S, tot_counts_MZ_S + counts_D_sum + counts_B_sum)

            with if_((counts_B_sum - counts_D_sum) > max_counts_B):
                assign(max_counts_B, (counts_B_sum - counts_D_sum))
                assign(phase_correction_min, phase_correction_value)

            save(phase_correction_min, phase_for_min_st)
            frame_rotation_2pi(phase_correction_per_scan, "AOM_Early")

    return phase_correction_min, tot_counts_MZ_S


def MZ_balancing_2(m_time, m_window, shutter_open_time, phase_rep_fast, phase_rep_slow, points_for_sum_fast,
                   points_for_sum_slow, counts_st_B, counts_st_D, phase_correction_st, phase_for_min_st,
                   phase_diff_last, max_tot_counts_MZ):

    counts_B = declare(int)
    counts_B_sum = declare(int, value=0)
    counts_D = declare(int)
    counts_D_sum = declare(int, value=0)
    tot_counts_MZ_S = declare(int, value=0)
    tot_counts_MZ_N = declare(int, value=0)

    OPX_Utils.assign_variables_to_element("AOM_Early", counts_B_sum)
    OPX_Utils.assign_variables_to_element("AOM_Late", counts_D_sum)
    OPX_Utils.assign_variables_to_element("PULSER_S", tot_counts_MZ_S)
    OPX_Utils.assign_variables_to_element("PULSER_N", tot_counts_MZ_N)
    tot_phase = declare(fixed)
    phase_correction = declare(fixed)
    ratio = declare(fixed)
    phase_rep = declare(int)
    points_for_sum = declare(int)
    i = declare(int)
    j = declare(int)

    phase_scan = declare(fixed, value=[1, 0.25])
    # phase_scan = declare(fixed, value=[1.0, 1.0])
    points_for_sum_vec = declare(int, value=[points_for_sum_fast, points_for_sum_slow])
    phase_correction_fast = declare(fixed, value=(1/phase_rep_fast))
    phase_correction_slow = declare(fixed, value=(1/phase_rep_slow))
    phase_correction_vec = declare(fixed, value=[(1/phase_rep_fast), (1/phase_rep_slow)])
    phase_rep_vec = declare(int, value=[phase_rep_fast, phase_rep_slow])
    phase_correction_per_scan = declare(fixed)
    phase_correction_value = declare(fixed)
    phase_correction_min = declare(fixed, value=0.0)
    phase_correction_min_south = declare(fixed, value=0.0)
    phase_correction_min_north = declare(fixed, value=0.0)
    phase_correction_diff = declare(fixed, value=0.8)
    phase_correction_diff_default = declare(fixed, value=0.8)
    phase_correction_diff_std = declare(fixed, value=0.06)
    phase_correction_start = declare(fixed, value=0.0)
    min_counts_D = declare(int, value=0)
    max_counts_B = declare(int, value=0)

    assign(max_counts_B, 0)
    assign(min_counts_D, 0)
    assign(phase_correction_min, 0.0)
    assign(phase_correction_start, 0.0)

    align("AOM_Early", "AOM_Late", "PULSER_N", "PULSER_S")
    # play("Const_open_triggered" * amp(0), "PULSER_N", duration=shutter_open_time)
    # play("Const_open" * amp(0), "PULSER_S", duration=shutter_open_time)

    align("AOM_Early", "AOM_Late", "PULSER_N", "PULSER_S")
    assign(phase_correction_per_scan, phase_correction_fast)
    # assign(phase_correction_start, -0.5)
    assign(max_counts_B, 0)
    assign(tot_counts_MZ_S, 0)
    reset_frame("AOM_Early", "AOM_Late", "PULSER_S", "PULSER_N")
    # frame_rotation_2pi(phase_correction_start, "AOM_Early")
    with for_(j, 0, j < phase_rep_fast, j + 1):
        assign(counts_B_sum, 0)
        assign(counts_D_sum, 0)
        with for_(i, 1, i <= points_for_sum_fast, i + 1):
            play("MZ_balancing_pulses_N", "PULSER_N")
            play("MZ_balancing_pulses_N", "PULSER_S")
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
        assign(tot_counts_MZ_S, tot_counts_MZ_S + counts_D_sum + counts_B_sum)

        with if_((counts_B_sum - counts_D_sum) > max_counts_B):
            assign(max_counts_B, (counts_B_sum - counts_D_sum))
            assign(phase_correction_min, phase_correction_value)

        save(phase_correction_min, phase_for_min_st)
        frame_rotation_2pi(phase_correction_per_scan, "AOM_Early")

    assign(phase_correction_min_south, phase_correction_min)
    assign(phase_correction_min, 0)

    align("AOM_Early", "AOM_Late", "PULSER_N", "PULSER_S")
    assign(phase_correction_per_scan, phase_correction_fast * 1)
    # assign(phase_correction_start, phase_correction_min - 1 * 0.5)
    assign(max_counts_B, 0)
    assign(tot_counts_MZ_N, 0)
    reset_frame("AOM_Early", "AOM_Late", "PULSER_S", "PULSER_N")
    # frame_rotation_2pi(phase_correction_start, "AOM_Early")
    with for_(j, 0, j < phase_rep_fast, j + 1):
        assign(counts_B_sum, 0)
        assign(counts_D_sum, 0)
        with for_(i, 1, i <= points_for_sum_fast, i + 1):
            play("MZ_balancing_pulses_N", "PULSER_N")
            play("MZ_balancing_pulses_N", "PULSER_S")
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
        assign(tot_counts_MZ_N, tot_counts_MZ_N + counts_D_sum + counts_B_sum)

        with if_((counts_B_sum - counts_D_sum) > max_counts_B):
            assign(max_counts_B, (counts_B_sum - counts_D_sum))
            assign(phase_correction_min, phase_correction_value)

        # assign(tot_counts_MZ_N, max_counts_B)
        save(phase_correction_min, phase_for_min_st)
        frame_rotation_2pi(phase_correction_per_scan, "AOM_Early")

    align("AOM_Early", "AOM_Late", "PULSER_N", "PULSER_S")
    assign(phase_correction_per_scan, phase_correction_slow * 0.25)
    assign(phase_correction_start, phase_correction_min - 0.25 * 0.5)
    assign(max_counts_B, 0)
    reset_frame("AOM_Early", "AOM_Late", "PULSER_S", "PULSER_N")
    frame_rotation_2pi(phase_correction_start, "AOM_Early")
    with for_(j, 0, j < phase_rep_slow, j + 1):
        assign(counts_B_sum, 0)
        assign(counts_D_sum, 0)
        with for_(i, 1, i <= points_for_sum_slow, i + 1):
            play("MZ_balancing_pulses_N", "PULSER_N")
            play("MZ_balancing_pulses_N", "PULSER_S")
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

        with if_((counts_B_sum - counts_D_sum) > max_counts_B):
            assign(max_counts_B, (counts_B_sum - counts_D_sum))
            assign(phase_correction_min, phase_correction_value)

        save(phase_correction_min, phase_for_min_st)
        frame_rotation_2pi(phase_correction_per_scan, "AOM_Early")

    assign(phase_correction_min_north, phase_correction_min)
    assign(phase_correction_diff, Util.cond(tot_counts_MZ_S > 0.20 * tot_counts_MZ_N,
                                            phase_correction_min_south - phase_correction_min_north,
                                            phase_correction_diff))

    assign(ratio, tot_counts_MZ_S/tot_counts_MZ_N)
    return phase_correction_min, phase_correction_diff, ratio

def MZ_balancing_check(m_time, m_window, rep,
                       counts_st_B_balanced, counts_st_D_balanced):

    counts_B = declare(int)
    counts_D = declare(int)

    t = declare(int)

    align("AOM_Early", "AOM_Late", "PULSER_N", "PULSER_S")
    # play("Const_open_triggered" * amp(0), "PULSER_ANCILLA", duration=m_time)
    with for_(t, 0, t < rep, t + 1):
        play("MZ_balancing_pulses_N", "PULSER_N")
        play("MZ_balancing_pulses_N", "PULSER_S")
        # play("MZ_balancing_pulses_S", "PULSER_N")
        # play("MZ_balancing_pulses_S", "PULSER_S")
        measure("MZ_balancing_pulses", "AOM_Early", None,
                counting.digital(counts_B, m_window, element_outputs="OutBright1"))
        measure("MZ_balancing_pulses", "AOM_Late", None,
                counting.digital(counts_D, m_window, element_outputs="OutDark2"))
        save(counts_B, counts_st_B_balanced)
        save(counts_D, counts_st_D_balanced)

def VSTIRAP_Exp(m_off_time, m_time, m_window, shutter_open_time,
             ON_counts_st1, ON_counts_st2, ON_counts_st3,
             ON_counts_st4, ON_counts_st5, ON_counts_st6,
             ON_counts_st7, ON_counts_st8,
             tt_st_1, tt_st_2, tt_st_3, tt_st_4, tt_st_5, tt_st_6, tt_st_7, tt_st_8, rep_st,
             balancing_check_window, length,
             rep_MZ_check, counts_st_B_balanced, counts_st_D_balanced, phase_correction):
    # TODO: the below comment looks irrelevant
    """
     Generates train of 8 pulses (to be configured in config - whether North or South)
     for SPRINT experiments

       Parameters
       ----------
       :param M_delay: The time delay from end of PGC until the first 2-3' pulse (OD + trigger).
       :param rep: The number of measuring window repetitions, derived from OD measuring duration (M_time/M_window).
       :param m_time: The duration of the entire measurement time.
       :param m_window: The duration of each measuring window - for each window there is 28 nsec "deadtime".
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

    OPX_Utils.assign_variables_to_element("Dig_detectors", counts1)

    n = declare(int, value=0)
    t = declare(int)
    m = declare(int)

    # OPX_Utils.assign_variables_to_element("Dig_detectors", tt_vec1[0], counts1, m_window)
    # align("AOM_Early", "AOM_Late", "PULSER_N", "PULSER_S", "Dig_detectors")
    # play("Const_open_triggered" * amp(0), "PULSER_N", duration=shutter_open_time)
    # play("Const_open" * amp(0), "PULSER_S", duration=shutter_open_time)
    # play("Const_open_triggered" * amp(0), "AOM_Early", duration=shutter_open_time)
    # # play("OD_FS" * amp(0), "AOM_2-2/3'")
    # # play("Const_open_triggered" * amp(0), "PULSER_ANCILLA", duration=shutter_open_time)

    align("AOM_Early", "AOM_Late", "PULSER_N", "PULSER_S", "Dig_detectors", "AOM_2-2/3'")
    #Sends trigger to AOM_Late
    play("Const_open", "AOM_Late", duration=(m_time + m_off_time))
    # # Sends trigger to FS
    play("Const_open_triggered", "PULSER_S", duration=(m_time + m_off_time))
    # #Sends trigger to shutters
    play("Const_open_triggered"*amp(0), "PULSER_N", duration=(m_time + m_off_time))
    play("OD_FS"*amp(0.2), "AOM_2-2/3'", duration=(m_time + m_off_time))

    # with for_(t, 0, t < (m_time + m_off_time) * 4, t + int(len(Config.VSTIRAP_Gaussian_pulse_samples))):  ## for VSTIRAP experiment ##
    # with for_(t, 0, t < (m_time + m_off_time) * 4, t + int(Config.vstirap_transits_pulse_len)): ## for transit experiment ##
    #     ## for transit experiment ##
    #

        # play("VSTIRAP_experiment_pulse", "PULSER_VSTIRAP_1_1")

        # wait(int(Config.frequency_sweep_duration / 4), "AOM_Spectrum")
        # align("PULSER_VSTIRAP_1_1", "PULSER_VSTIRAP_1_0")
        # wait(int(Config.frequency_sweep_duration / 4), "PULSER_ANCILLA")
        # play("VSTIRAP_experiment_pulse", "PULSER_VSTIRAP_1_0")

    # wait(298, "Dig_detectors")
    # wait(300-12, "Dig_detectors")
    wait(288, "Dig_detectors")
    # with for_(n, 0, n < m_time * 4, n + m_window):
    measure("readout_PNSA", "Dig_detectors", None,
            time_tagging.digital(tt_vec1, m_window, element_output="out1", targetLen=counts1),
            time_tagging.digital(tt_vec2, m_window, element_output="out2", targetLen=counts2),
            time_tagging.digital(tt_vec3, m_window, element_output="out3", targetLen=counts3),
            time_tagging.digital(tt_vec4, m_window, element_output="out4", targetLen=counts4),
            time_tagging.digital(tt_vec5, m_window, element_output="out5", targetLen=counts5),
            time_tagging.digital(tt_vec6, m_window, element_output="out6", targetLen=counts6),
            time_tagging.digital(tt_vec7, m_window, element_output="out7", targetLen=counts7),
            time_tagging.digital(tt_vec8, m_window, element_output="out8", targetLen=counts8),
            )

    reset_frame("AOM_Early", "AOM_Late", "PULSER_S", "PULSER_N")
    frame_rotation_2pi(phase_correction, "AOM_Early")
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
        Trigger_Phase = declare(int, value=obj.Exp_Values['Triggering_Phase'])
        Imaging_Phase = declare(int, value=obj.Exp_Values['Imaging_Phase'])
        x = declare(int)

        # Boolean variables:
        AntiHelmholtz_ON = declare(bool, value=True)
        Exp_ON = declare(bool, value=False)

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

        # MZ balancing vaiables:
        phase_correction = declare(fixed)
        phase_correction_diff = declare(fixed)
        max_counts_MZ = declare(fixed, value=0)

        # Stream processing:
        counts_st_B = declare_stream()
        counts_st_B_balanced = declare_stream()
        tt_st_B = declare_stream()
        counts_st_D = declare_stream()
        counts_st_D_balanced = declare_stream()
        tt_st_D = declare_stream()
        phase_correction_st = declare_stream()
        phase_for_min_st = declare_stream()
        phase_correction_diff_st = declare_stream()
        max_counts_MZ_st = declare_stream()
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
            # TODO: What is this calc?
            with if_(Exp_ON):
                assign(x, (30678780 - 3106 + 656000 * 2 + 4) // 4)  # TODO - added 38688900 to fix new delay due to wait(1000) in saving sprint data with vector size 10000, should be fixed as well
            with else_():
                assign(x, - 3106)
            # FreeFall(FreeFall_duration, coils_timing)
            FreeFall(FreeFall_duration - x, coils_timing)
            # play("Const_open_triggered" * amp(0), "PULSER_ANCILLA", duration=(PrePulse_duration - shutter_open_time - Balancing_check_window))

            ##########################
            ## Measurement Sequence ##
            ##########################

            # TODO: if we use Enums, it will be like this:
            #with if_(Trigger_Phase == Phases.FREE_FALL.value)
            # TODO: but I see that in Phases_Names, #3 is "Free_Fall" - is it the same like "PrePulse"?
            with if_(Trigger_Phase == Phases.FREE_FALL):  # when trigger on PrePulse
                ## Trigger QuadRF Sequence #####################
                play("C_Seq", "Cooling_Sequence", duration=2500)
                ################################################

            align("Cooling_Sequence", "AOM_Early", "AOM_Late")
            # wait(2500, "AOM_2-2'")
            with if_(Exp_ON):
                wait(PrePulse_duration - shutter_open_time, "Cooling_Sequence")
                # play("Depump", "AOM_2-2'", duration=(PrePulse_duration - shutter_open_time))
                phase_correction, max_counts_MZ = MZ_balancing((PrePulse_duration - shutter_open_time - Balancing_check_window), len(Config.PNSA_MZ_balance_pulse_Late),
                                                               shutter_open_time, obj.phase_rep_MZ, obj.points_for_sum,
                                                               counts_st_B, counts_st_D, phase_correction_st, phase_for_min_st, tt_st_B, tt_st_D)
                # phase_correction, phase_correction_diff, max_counts_MZ = \
                #     MZ_balancing_2((PrePulse_duration - shutter_open_time - Balancing_check_window),
                #                    len(Config.PNSA_MZ_balance_pulse_Late), shutter_open_time,
                #                    obj.phase_rep_MZ_fast_scan, obj.phase_rep_MZ_slow_scan, obj.points_for_sum_fast,
                #                    obj.points_for_sum_slow, counts_st_B, counts_st_D, phase_correction_st,
                #                    phase_for_min_st, phase_correction_diff, max_counts_MZ)
                save(phase_correction_diff, phase_correction_diff_st)
                save(max_counts_MZ, max_counts_MZ_st)
                reset_frame("AOM_Early", "AOM_Late", "PULSER_S", "PULSER_N")
                frame_rotation_2pi(phase_correction, "AOM_Early")
                MZ_balancing_check(Balancing_check_window, len(Config.PNSA_MZ_balance_pulse_Late),
                                   obj.rep_MZ_check, counts_st_B_balanced, counts_st_D_balanced)
            with else_():
                # play("Depump", "AOM_2-2'", duration=PrePulse_duration)
                wait(PrePulse_duration, "Cooling_Sequence")
            align(*all_elements, "AOM_2-2/3'", "AOM_Early", "AOM_Late", "PULSER_N", "PULSER_S", "Dig_detectors")

            with if_(Trigger_Phase == Phases.PULSE_1):  # when trigger on pulse 1
                ## Trigger QuadRF Sequence #####################
                play("C_Seq", "Cooling_Sequence", duration=2500)
                ################################################
            ## For taking an image:
            with if_((Pulse_1_duration > 0) & ~Exp_ON):
                align(*all_elements)
                Pulse_with_prep(Pulse_1_duration, Pulse_1_decay_time, pulse_1_duration_0,
                                pulse_1_duration_minus, pulse_1_duration_plus)
                Measure(Pulse_1_duration)  # This triggers camera (Control 7)
                align(*all_elements)

            with if_((Pulse_1_duration > 0) & Exp_ON):
                align("Dig_detectors", "AOM_Early", "AOM_Late", "PULSER_N", "PULSER_S","AOM_2-2/3'")
                reset_frame("AOM_Early", "AOM_Late", "PULSER_S", "PULSER_N")
                # frame_rotation_2pi(phase_correction + phase_correction_diff, "AOM_Early")
                frame_rotation_2pi(phase_correction, "AOM_Early")

                # Run the experiment part!
                VSTIRAP_Exp(m_off_time=M_off_time, m_time=Pulse_1_duration, m_window=obj.M_window, shutter_open_time=shutter_open_time,
                         ON_counts_st1=ON_counts_st1, ON_counts_st2=ON_counts_st2, ON_counts_st3=ON_counts_st3,
                         ON_counts_st4=ON_counts_st4, ON_counts_st5=ON_counts_st5, ON_counts_st6=ON_counts_st6,
                         ON_counts_st7=ON_counts_st7, ON_counts_st8=ON_counts_st8,
                         tt_st_1=tt_st_1, tt_st_2=tt_st_2, tt_st_3=tt_st_3, tt_st_4=tt_st_4,
                         tt_st_5=tt_st_5, tt_st_6=tt_st_6, tt_st_7=tt_st_7, tt_st_8=tt_st_8,
                         rep_st=rep_st,
                         balancing_check_window=Balancing_check_window,
                         length=len(Config.PNSA_MZ_balance_pulse_Late),
                         rep_MZ_check=obj.rep_MZ_check,
                         counts_st_B_balanced=counts_st_B_balanced,
                         counts_st_D_balanced=counts_st_D_balanced,
                         phase_correction=phase_correction)
                align("Dig_detectors", "AOM_Early", "AOM_Late", "PULSER_N", "PULSER_S","AOM_2-2/3'")

            save(AntiHelmholtz_ON, AntiHelmholtz_ON_st)
            with if_(AntiHelmholtz_ON):
                save(FLR, FLR_st)

            assign(N_Snaps, 1)
            assign(Buffer_Cycles, 0)
            assign(i, IO1)

            ## PARAMETERS UPDATE ##
            with if_(i > 0):
                pause()

            # We will break the loop only when receiving a value of 0 - this marks the end of parameters to be updated
            with while_(i > 0):
                ## Boolean variables control: ##
                with if_(i == IOP.MOT_SWITCH_ON.value):
                    update_frequency("MOT_AOM_0", Config.IF_AOM_MOT)
                with if_(i == IOP.MOT_SWITCH_OFF.value):
                    update_frequency("MOT_AOM_0", Config.IF_AOM_MOT_OFF)
                with if_(i == IOP.EXPERIMENT_SWITCH.value):
                    assign(Exp_ON, IO2)
                # ## AntiHelmholtz control ##
                with if_(i == IOP.ANTIHELMHOLTZ_DELAY.value):
                    assign(antihelmholtz_delay, IO2)
                # ## MOT variables control ##
                # with if_(i == 11):  # Live control over the
                #     assign(MOT_Repetitions, IO2)
                with if_(i == IOP.POST_MOT_DELAY.value):  # Live control over the
                    assign(post_MOT_delay, IO2)

                # ## PGC variables control ##
                with if_(i == IOP.PGC_DURATION.value):  # Live control over the PGC duration
                    assign(pgc_duration, IO2)

                # ## Fountain variables control: ##
                with if_(i == IOP.FOUNTAIN_DURATION.value):  # Live control over the fountain duration
                    assign(fountain_duration, IO2)
                # with if_(i == 36):  # Live control over the final amplitude of the fountain AOM 0
                #     assign(fountain_pulse_duration_0, IO2)
                # with if_(i == 37):  # Live control over the final amplitude of the fountain AOM -
                #     assign(fountain_pulse_duration_minus, IO2)
                # with if_(i == 38):  # Live control over the final amplitude of the fountain AOM +
                #     assign(fountain_pulse_duration_plus, IO2)
                with if_(i == IOP.FOUNTAIN_FINAL_DELTA_FREQ.value):  # Live control over the fountain final frequency of the MOT + and - AOMs
                    assign(fountain_aom_chirp_rate, IO2)
                #
                # ## Measurement variables control: ##
                with if_(i == IOP.PREPULSE_DURATION.value):
                    assign(PrePulse_duration, IO2)
                with if_(i == IOP.PULSE_1_DURATION.value):
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

        # TODO: we want to go over steams and use them to "save" the data to the streams
        with stream_processing():
            # counts_st_B.buffer(obj.total_phase_rep_MZ).save('Bright_Port_Counts')
            # counts_st_B_balanced.save('Bright_Port_Counts_check')
            # counts_st_B.buffer(obj.total_phase_rep_MZ_scan).zip(counts_st_B_balanced.buffer(obj.rep_MZ_check*2)).save('Bright_Port_Counts')
            counts_st_B.buffer(obj.total_phase_rep_MZ).zip(counts_st_B_balanced.buffer(obj.rep_MZ_check*2)).save('Bright_Port_Counts')
            # counts_st_B.buffer(obj.total_phase_rep_MZ_scan).save('Bright_Port_Counts')
            # counts_st_D.buffer(obj.total_phase_rep_MZ).save('Dark_Port_Counts')
            # counts_st_D_balanced.save('Dark_Port_Counts_check')
            counts_st_D.buffer(obj.total_phase_rep_MZ).zip(counts_st_D_balanced.buffer(obj.rep_MZ_check*2)).save('Dark_Port_Counts')
            # counts_st_D.buffer(obj.total_phase_rep_MZ_scan).zip(counts_st_D_balanced.buffer(obj.rep_MZ_check*2)).save('Dark_Port_Counts')

            # ZIP together 'phase_correction_st' and 'phase_for_min_st'
            phase_correction_st.buffer(obj.total_phase_rep_MZ).zip(phase_for_min_st.buffer(obj.total_phase_rep_MZ)).save('Phase_Correction_array')
            # phase_correction_st.buffer(obj.total_phase_rep_MZ_scan).zip(phase_for_min_st.buffer(obj.total_phase_rep_MZ_scan)).save('Phase_Correction_array')

            # phase_correction_diff_st.save('Phase_diff')
            max_counts_MZ_st.save('Max_counts')
            # ON_counts_st1.buffer(obj.rep).save('Detector_1_Counts')
            # ON_counts_st2.buffer(obj.rep).save('Detector_2_Counts')
            # ON_counts_st3.buffer(obj.rep).save('Detector_3_Counts')
            # ON_counts_st4.buffer(obj.rep).save('Detector_4_Counts')
            # ON_counts_st5.buffer(obj.rep).save('Detector_5_Counts')
            # ON_counts_st6.buffer(obj.rep).save('Detector_6_Counts')
            # ON_counts_st7.buffer(obj.rep).save('Detector_7_Counts')
            # ON_counts_st8.buffer(obj.rep).save('Detector_8_Counts')
            (tt_st_1 + rep_st).buffer(obj.vec_size * obj.rep + 1).save('Detector_1_Timetags')
            (tt_st_2 + rep_st).buffer(obj.vec_size * obj.rep + 1).save('Detector_2_Timetags')
            (tt_st_3 + rep_st).buffer(obj.vec_size * obj.rep + 1).save('Detector_3_Timetags')
            (tt_st_4 + rep_st).buffer(obj.vec_size * obj.rep + 1).save('Detector_4_Timetags')
            (tt_st_5 + rep_st).buffer(obj.vec_size * obj.rep + 1).save('Detector_5_Timetags')
            (tt_st_6 + rep_st).buffer(obj.vec_size * obj.rep + 1).save('Detector_6_Timetags')
            (tt_st_7 + rep_st).buffer(obj.vec_size * obj.rep + 1).save('Detector_7_Timetags')
            (tt_st_8 + rep_st).buffer(obj.vec_size * obj.rep + 1).save('Detector_8_Timetags')
            FLR_st.save('FLR_measure')

            # TODO: this is not in the streams definitions. So where does it come from?
            AntiHelmholtz_ON_st.save("antihelmholtz_on")

    # sourceFile = open('serialized_code.py', 'w')
    # print(generate_qua_script(opx_control_prog, Config.config), file=sourceFile)
    # sourceFile.close()

    job = qm.execute(opx_control_prog, flags=['auto-element-thread'])

    return job
