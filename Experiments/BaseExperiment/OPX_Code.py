import Config_Experiment as Config
from Experiments.BaseExperiment.Config_Table import Phases_Names, Values_Factor
from qm.qua import *

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
        # play("OD_FS" * amp(0.5), "AOM_2-3'_for_interference")

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
          "Measurement") # , "Dig_detectors", "Dig_detectors_spectrum")

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


def Transits_Exp_TT(m_off_time, m_time, m_window, shutter_open_time,
                    ON_counts_st1, ON_counts_st2, ON_counts_st3,
                    ON_counts_st6, ON_counts_st7, ON_counts_st8,
                    tt_st_1, tt_st_2, tt_st_3, tt_st_6, tt_st_7, tt_st_8, rep_st):
    """
     Transit measurement with and without atoms.

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

    align("AOM_2-2/3'", "Dig_detectors")
    play("OD", "AOM_2-2/3'", duration=shutter_open_time)
    align("AOM_2-2/3'", "Dig_detectors")
    play("OD", "AOM_2-2/3'", duration=(m_time + m_off_time))
    with for_(n, 0, n < (m_time * 4), n + m_window):
        measure("readout_CRUS", "Dig_detectors", None,
                time_tagging.digital(tt_vec1, m_window, element_output="out1", targetLen=counts1),
                time_tagging.digital(tt_vec2, m_window, element_output="out2", targetLen=counts2),
                time_tagging.digital(tt_vec3, m_window, element_output="out3", targetLen=counts3),
                time_tagging.digital(tt_vec6, m_window, element_output="out6", targetLen=counts6),
                time_tagging.digital(tt_vec7, m_window, element_output="out7", targetLen=counts7),
                time_tagging.digital(tt_vec8, m_window, element_output="out8", targetLen=counts8),
                )

        ## Save Data: ##

        ## Number of Photons (NOP) Count stream for each detector: ##
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


def Probe_counts_Measure_SNSPDs(m_off_time, m_time, m_window, shutter_open_time,
                                ON_counts_st1, ON_counts_st2, ON_counts_st3, ON_counts_st4,
                                ON_counts_st5, ON_counts_st6, ON_counts_st7, ON_counts_st8, ):
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


def opx_control(obj, qm):
    with program() as opx_control_prog:
        ## declaring program variables: ##
        i = declare(int)
        filler = declare(int, value=1)
        Trigger_Phase = declare(int, value=Phases_Names.index(obj.Exp_Values['Triggering_Phase']))
        Imaging_Phase = declare(int, value=Phases_Names.index(obj.Exp_Values['Imaging_Phase']))

        # Boolean variables:
        MOT_ON = declare(bool, value=True)
        AntiHelmholtz_delay_ON = declare(bool, value=False)
        AntiHelmholtz_ON = declare(bool, value=True)
        Linear_PGC_ON = declare(bool, value=True)
        Transits_Exp_ON = declare(bool, value=False)
        CRUS_Exp_ON = declare(bool, value=False)
        Probe_max_counts_Exp_ON = declare(bool, value=False)

        # MOT variables
        MOT_Repetitions = declare(int, value=obj.Exp_Values['MOT_rep'])
        post_MOT_delay = declare(int, value=int(obj.Exp_Values['Post_MOT_delay'] * 1e6 / 4))

        # AntiHelmholtz delay after MOT:
        antihelmholtz_delay = declare(int, value=int(obj.Exp_Values['AntiHelmholtz_delay'] * 1e6 / 4))

        # PGC variables:
        pgc_duration = declare(int, value=int(obj.pgc_duration))
        pgc_prep_time = declare(int, value=int(obj.pgc_prep_duration))
        pgc_pulse_duration_0 = declare(int, value=int(obj.pgc_pulse_duration_0))  # The relative duration to reach the desired amplitude
        pgc_pulse_duration_minus = declare(int, value=int(obj.pgc_pulse_duration_minus))  # The relative duration to reach the desired amplitude
        pgc_pulse_duration_plus = declare(int, value=int(obj.pgc_pulse_duration_plus))  # The relative duration to reach the desired amplitude

        # Fountain variables:
        fountain_duration = declare(int, value=int(obj.fountain_duration))
        # TODO: pre_PGC_fountain_duration = declare(int, value=int(obj.pre_PGC_fountain_duration))
        fountain_prep_time = declare(int, value=int(obj.fountain_prep_duration))  # Can't be used with Chirp!!!
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
        M_duration = declare(int, value=int(obj.M_time / 4))  # From [nsec] to [4 nsec]
        M_off_time = declare(int, value=int(obj.M_off_time / 4))
        OD_Delay = declare(int, value=int(obj.OD_delay * 1e6 / 4))
        shutter_open_time = declare(int, value=int(obj.Shutter_open_time * 1e6 / 4))
        Rep = declare(int, value=obj.rep)

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
            with if_(Transits_Exp_ON):
                assign(x, 766000 // 4)
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

            with if_(Transits_Exp_ON):
                wait(PrePulse_duration - shutter_open_time, "Cooling_Sequence")
            with else_():
                wait(PrePulse_duration, "Cooling_Sequence")
            align(*all_elements, "AOM_2-2/3'", "AOM_2-2'", "Dig_detectors") # , "Dig_detectors"

            with if_(Trigger_Phase == 4):  # when trigger on pulse 1
                ## Trigger QuadRF Sequence #####################
                play("C_Seq", "Cooling_Sequence", duration=2500)
                ################################################
            ## For taking an image:
            with if_((Pulse_1_duration > 0) & ~Transits_Exp_ON):
                align(*all_elements)
                Pulse_with_prep(Pulse_1_duration, Pulse_1_decay_time, pulse_1_duration_0,
                                pulse_1_duration_minus, pulse_1_duration_plus)
                Measure(Pulse_1_duration)  # This triggers camera (Control 7)
                align(*all_elements)

            with if_((Pulse_1_duration > 0) & Transits_Exp_ON):
                align(*all_elements, "Dig_detectors", "AOM_2-2/3'")
                Transits_Exp_TT(M_off_time, Pulse_1_duration, obj.M_window, shutter_open_time,
                                ON_counts_st1, ON_counts_st2, ON_counts_st3,
                                ON_counts_st6, ON_counts_st7, ON_counts_st8,
                                tt_st_1, tt_st_2, tt_st_3, tt_st_6, tt_st_7, tt_st_8, rep_st)
                save(AntiHelmholtz_ON, AntiHelmholtz_ON_st)
                with if_(AntiHelmholtz_ON):
                    save(FLR, FLR_st)
                align(*all_elements, "Dig_detectors", "AOM_2-2/3'")

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
                    assign(CRUS_Exp_ON, IO2)
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
