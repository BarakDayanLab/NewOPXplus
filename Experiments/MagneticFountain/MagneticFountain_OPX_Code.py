import Experiments.BaseExperiment.Config_Experiment as Config
from Experiments.Enums.IOParameters import IOParameters as IOP
from Experiments.Enums.Phases import Phases
from Utilities.OPX_Utils import OPX_Utils
from qm.qua import *

all_elements = ["Cooling_Sequence", "MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+", "AntiHelmholtz_Coils", "Measurement", "Magnetic_Fountain"]


def MOT(mot_repetitions, OD_Attenuation):
    """
    The MOT function is used to play the MOT. To that end, we send RF signal to AOM_0, AOM_+, AOM_- for the duration of the MOT.

    Parameters
    ----------
    :param mot_repetitions: derived from the MOT duration, and is calculated in the Experiment parameters section
    """
    FLR = declare(fixed)
    # OPX_Utils.assign_variables_to_element("FLR_detection", FLR)

    align("Cooling_Sequence", "MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+", "AntiHelmholtz_Coils", "Zeeman_Coils",
          "AOM_2-2/3'", "AOM_2-3'_for_interference", "AOM_2-2'", "FLR_detection", "Measurement",
          "Magnetic_Fountain") # , "Dig_detectors_spectrum", "Dig_detectors") # , "PULSER_N", "PULSER_S")

    ## MOT build-up ##
    n = declare(int)
    m = declare(int)
    play("Detection" * amp(FLR * 0), "FLR_detection", duration=4)  # we don't know why this works, but Yoav from QM made us write this line to solve an alignment problem we had in the next 2 for loops
    with for_(n, 1, n <= mot_repetitions, n + 1):
        play("MOT" * amp(Config.AOM_0_Attenuation), "MOT_AOM_0")
        play("MOT" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-")
        play("MOT" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+")

        # OD beam
        # play("OD_FS" * amp(OD_Attenuation), "AOM_2-2/3'")
        # In the case of configuring Funnel, we use the OD laser to see if we hit the atom cloud. Therefore, we reduce a bit the amp to the MOT, increasing what goes to the OD
        # (this will "hurt" the MOT a bit)
        # play("OD_FS" * amp(0.3), "AOM_2-2/3'")
        play("AntiHelmholtz_MOT", "AntiHelmholtz_Coils")
    with for_(m, 1, m <= (mot_repetitions - 1), m + 1):
        measure("Detection", "FLR_detection", None, integration.full("Detection_opt", FLR, "out1"))

        # play("OD_FS" * amp(OD_Attenuation), "AOM_2-2/3'")

        # play("OD_FS" * amp(0.5), "AOM_2-3'_for_interference")

    align("Cooling_Sequence", "MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+", "AntiHelmholtz_Coils", "Zeeman_Coils",
          "AOM_2-2/3'", "AOM_2-2'", "FLR_detection", "Measurement") # , "Dig_detectors", "Dig_detectors_spectrum") #, "PULSER_N", "PULSER_S")

    return FLR


def PGC(pgc_duration, pgc_prep_duration, pgc_beams_0_off_duration, fountain_aom_chirp_rate, magnetic_fountain_duration,AOM_0_att_3rd_stage):
    """
    The PGC process is controlled both by OPX and QuadRF.
    - The OPX code controls the AOMs (either by constant pulse, or turning on/off the 0 beams
    - The Quad code controls the gradient and plateau for the Eilam Pulser shifting the frequency
    """

    # There are 3 options to run PGC:

    # --- Option 1 ---- => Plateau only
    # Pulse_const(pgc_duration)

    # --- Option 2 ---- => Gradient Descent and then Plateau
    # Pulse_with_prep_without_0_at_pgc(pgc_duration, pgc_prep_duration, pgc_prep_duration, pgc_prep_duration,
    #                                  pgc_prep_duration, fountain_aom_chirp_rate)

    # --- Option 3 ---- => Gradient Descent and then Plateau with a period at the end turning off 0-beams
    Three_stage_pgc(pgc_duration, pgc_prep_duration, pgc_beams_0_off_duration, magnetic_fountain_duration,AOM_0_att_3rd_stage)


def Pulse_const(total_pulse_duration):
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

    ## Playing the pulses to the AOMs for the constant part. (Qua) ##
    align("MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+")
    play("Const" * amp(Config.AOM_0_Attenuation), "MOT_AOM_0", duration=total_pulse_duration)
    play("Const" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-", duration=total_pulse_duration)
    play("Const" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+", duration=total_pulse_duration)


def Pulse_with_prep(total_pulse_duration, prep_duration, zero_pulse_duration, plus_pulse_duration,
                    minus_pulse_duration):
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
        play("Const" * amp(Config.AOM_0_Attenuation_pulse_1), "MOT_AOM_0",
             duration=(total_pulse_duration - prep_duration))
        play("Const" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-",
             duration=(total_pulse_duration - prep_duration))
        play("Const" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+",
             duration=(total_pulse_duration - prep_duration))


def Magnetic_Fountain(fountain_duration, prep_duration, zero_pulse_duration, plus_pulse_duration,
                               minus_pulse_duration, aom_chirp_rate, magnetic_fountain_duration):

    align("MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+", "Magnetic_Fountain")
    play("Const_open" * amp(1.0), "Magnetic_Fountain", duration=magnetic_fountain_duration)
    pass

def Pulse_with_prep_with_chirp(fountain_duration, prep_duration, zero_pulse_duration, plus_pulse_duration,
                               minus_pulse_duration, aom_chirp_rate, delta_f):
    """
    This pulse id divided to two parts:
        1) Preparation sequence where different parameters changes such as amplitudes and frequencies of the 0 - + AOMs (scan each separately)
        2) Part where we keep the different values constant for the 0 - + AOMs

    Parameters
    ----------
    :param fountain_duration: duration of the whole pulse.
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
    with if_(fountain_duration > prep_duration):
        play("Const" * amp(Config.AOM_0_Attenuation), "MOT_AOM_0",
             duration=(fountain_duration - prep_duration))
        play("Const" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-",
             duration=(fountain_duration - prep_duration))
        play("Const" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+",
             duration=(fountain_duration - prep_duration))

def Pulse_with_prep_without_0_at_pgc(pulse_duration, prep_duration, zero_pulse_duration, plus_pulse_duration,
                               minus_pulse_duration, aom_chirp_rate):
    """
    This pulse id divided to two parts:
        1) Preparation sequence where different parameters changes such as amplitudes and frequencies of the 0 - + AOMs (scan each separately)
        2) Part where we keep the different values constant for the 0 - + AOMs

    Parameters
    ----------
    :param fountain_duration: duration of the whole pulse.
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
            play("Const" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-", duration=prep_duration)
            play("Const" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+", duration=prep_duration)
        with else_():
            play("Linear" * amp(Config.AOM_0_Attenuation), "MOT_AOM_0", duration=zero_pulse_duration,
                 truncate=prep_duration)
            play("Linear" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-", duration=minus_pulse_duration,
                 truncate=prep_duration, chirp=(aom_chirp_rate, "mHz/nsec"))
            play("Linear" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+", duration=plus_pulse_duration,
                 truncate=prep_duration, chirp=(-aom_chirp_rate, "mHz/nsec"))


    ## Playing the pulses to the AOMs for the constant part. (Qua) ##
    align("MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+")
    update_frequency("MOT_AOM_-", Config.IF_AOM_MOT)
    update_frequency("MOT_AOM_+", Config.IF_AOM_MOT)
    with if_(pulse_duration > prep_duration):
        play("Const" * amp(0.15), "MOT_AOM_0",
             duration=(pulse_duration - prep_duration))
        play("Const" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-",
             duration=(pulse_duration - prep_duration))
        play("Const" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+",
             duration=(pulse_duration - prep_duration))


def Three_stage_pgc(pgc_duration, pgc_prep_duration, pgc_beams_off_duration, magnetic_fountain_duration,AOM_0_att_3rd_stage):

    """
    This pulse id divided to two parts:
        1) Preparation sequence where different parameters changes such as amplitudes and frequencies of the 0 - + AOMs (scan each separately)
        2) Part where we keep the different values constant for the 0 - + AOMs

    Parameters
    ----------
    :param pgc_duration: duration of the whole pulse.
    :param pgc_prep_duration: duration of the preparation part of the pulse
    """

    # pgc_duration => The total time of the PGC process
    # - pgc_prep_duration => the gradient descent part
    # - pgc_beams_off_duration => the plateau part - 0 beams are turned OFF
    # - pgc_duration - pgc_prep_duration - pgc_beams_off_duration -> the plateau part - 0 beams are turned ON

    ## Playing the pulses to the AOMs for the preparation sequence. (Qua) ##
    align("MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+")
    with if_(pgc_prep_duration > 0):
        play("Const" * amp(Config.AOM_0_Attenuation), "MOT_AOM_0", duration=pgc_prep_duration)
        play("Const" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-", duration=pgc_prep_duration)
        play("Const" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+", duration=pgc_prep_duration)

    align("MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+")
    with if_(pgc_duration > pgc_prep_duration+pgc_beams_off_duration):
        play("Const" * amp(Config.AOM_0_Attenuation), "MOT_AOM_0", duration=(pgc_duration - pgc_prep_duration - pgc_beams_off_duration))
        play("Const" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-", duration=(pgc_duration - pgc_prep_duration - pgc_beams_off_duration))
        play("Const" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+", duration=(pgc_duration - pgc_prep_duration - pgc_beams_off_duration))

    align("MOT_AOM_0", "MOT_AOM_-", "MOT_AOM_+", "Magnetic_Fountain")

    wait(us(50), "Magnetic_Fountain")
    play("Const_open" * amp(1.0), "Magnetic_Fountain", duration=magnetic_fountain_duration)

    with if_(pgc_beams_off_duration > 0):
        play("Const" * amp(Config.AOM_0_Attenuation), "MOT_AOM_0", duration=pgc_beams_off_duration)
        play("Const" * amp(Config.AOM_Minus_Attenuation), "MOT_AOM_-", duration=pgc_beams_off_duration)
        play("Const" * amp(Config.AOM_Plus_Attenuation), "MOT_AOM_+", duration=pgc_beams_off_duration)


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
          "Measurement", "Magnetic_Fountain") # , "Dig_detectors", "Dig_detectors_spectrum")

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


def opx_control_dbg(obj, qm):

    with program() as opx_control_prog:
        i = declare(int)
        assign(i, 5)

    job = qm.execute(opx_control_prog, flags=['auto-element-thread'])

    return job

def ms(val):
    return ms(val*1000)
def us(val):
    return ns(val*1000)
def ns(val):
    return val >> 2  # divide by 4

def opx_control(obj, qm):
    with program() as opx_control_prog:

        # ----------------------------------
        # Declaring program variables
        # ----------------------------------

        # TODO: when mechanism is complete, uncomment the below and test @@@ - 1
        # Populate parameters vector with all experiment values
        #p = OPX_Utils.assign_experiment_variables(obj)

        i = declare(int)
        filler = declare(int, value=1)
        Trigger_Phase = declare(int, value=int(obj.Exp_Values['Triggering_Phase']))
        Imaging_Phase = declare(int, value=int(obj.Exp_Values['Imaging_Phase']))

        # Boolean variables:
        Magnetic_Fountain_ON = declare(bool, value=True)

        # MOT variables
        MOT_Repetitions = declare(int, value=obj.Exp_Values['MOT_rep'])

        # TODO: why do we have here * 1e6 / 4 - if the Python code factors it before sending?
        post_MOT_delay = declare(int, value=int(obj.Exp_Values['Post_MOT_delay'] * 1e6 / 4))

        # AntiHelmholtz delay after MOT:
        antihelmholtz_delay = declare(int, value=int(obj.Exp_Values['AntiHelmholtz_delay'] * 1e6 / 4))

        # PGC variables:
        pgc_duration = declare(int, value=int(obj.pgc_duration))
        pgc_prep_duration = declare(int, value=int(obj.pgc_prep_duration))
        pgc_beams_0_off_duration = declare(int, value=int(obj.pgc_beams_0_off_duration))

        pgc_pulse_duration_0 = declare(int, value=int(obj.pgc_pulse_duration_0))  # The relative duration to reach the desired amplitude
        pgc_pulse_duration_minus = declare(int, value=int(obj.pgc_pulse_duration_minus))  # The relative duration to reach the desired amplitude
        pgc_pulse_duration_plus = declare(int, value=int(obj.pgc_pulse_duration_plus))  # The relative duration to reach the desired amplitude

        # Fountain variables:
        fountain_duration = declare(int, value=int(obj.fountain_duration))
        magnetic_fountain_duration = declare(int, value=int(obj.magnetic_fountain_duration))
        AOM_0_att_3rd_stage = declare(fixed,value=obj.Exp_Values['AOM_0_att_3rd_stage'])

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
        OD_attenuation = declare(fixed, value=obj.Exp_Values['OD_continuous_attenuation'])
        shutter_open_time = declare(int, value=int(obj.Shutter_open_time * 1e6 / 4))
        Rep = declare(int, value=obj.rep)

        # MW spectroscopy variables:
        MW_freq = declare(int, value=obj.MW_start_frequency)
        pulse_length_MW = declare(int, value=(obj.Pulse_Length_MW * int(1e3 / 4)))
        pulse_length_OD = declare(int, value=(obj.Pulse_Length_OD * int(1e3 / 4)))
        Repetitions = declare(int, value=1)
        Delta_f = declare(int, value=0)

        # TODO: can remove the following - not in use by the generic MOT sequence
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
            FLR = MOT(MOT_Repetitions, OD_attenuation)
            # play("AntiHelmholtz_MOT", "AntiHelmholtz_Coils", duration=antihelmholtz_delay)

            # Delay before fountain:

            wait(post_MOT_delay, "Cooling_Sequence")
            #wait(params[IOP.POST_MOT_DELAY], "Cooling_Sequence")  # @@@ - 2  TODO: Generalize this - for all parameters :-)

            align(*all_elements)

            # Optical Fountain sequence:
            with if_(fountain_duration > 0):
                wait(fountain_duration, "Cooling_Sequence")
                Pulse_with_prep_with_chirp(fountain_duration, obj.fountain_prep_duration,
                                           fountain_pulse_duration_0, fountain_pulse_duration_minus,
                                           fountain_pulse_duration_plus, fountain_aom_chirp_rate,
                                           fountain_delta_f)

            align(*all_elements)

            # wait(pgc_duration, "Cooling_Sequence")
            PGC(pgc_duration, pgc_prep_duration, pgc_beams_0_off_duration, fountain_aom_chirp_rate,
                magnetic_fountain_duration,AOM_0_att_3rd_stage)

            align(*all_elements)

            # TODO: make the opx_quadrf_misalignment_delay a parameter
            # FreeFall sequence:
            misalignment_delay = 0 - obj.Exp_Values['OPX_Quad_Misalignment_Delay']
            assign(x, ns(misalignment_delay))
            #assign(x, (0 - obj.Exp_Values['OPX_Quad_Misalignment_Delay']) // 4)
            FreeFall(FreeFall_duration - x, coils_timing)

            ##########################
            ## Measurement Sequence ##
            ##########################

            with if_(Trigger_Phase == Phases.FREE_FALL):  # when trigger on PrePulse
                ## Trigger QuadRF Sequence #####################
                play("C_Seq", "Cooling_Sequence", duration=us(10))
                ################################################

            wait(PrePulse_duration, "Cooling_Sequence")
            align(*all_elements, "AOM_2-2/3'", "AOM_2-2'", "Dig_detectors")  # , "Dig_detectors"

            with if_(Trigger_Phase == Phases.PULSE_1):  # when trigger on pulse 1
                ## Trigger QuadRF Sequence #####################
                play("C_Seq", "Cooling_Sequence", duration=us(10))
                ################################################

            # For taking an image:
            with if_(Pulse_1_duration > 0):
                align(*all_elements, "AOM_2-2/3'")
                Pulse_with_prep(Pulse_1_duration, Pulse_1_decay_time, pulse_1_duration_0,
                                pulse_1_duration_minus, pulse_1_duration_plus)
                Measure(Pulse_1_duration)  # This triggers camera (Control 7)
                align(*all_elements, "AOM_2-2/3'")

            assign(N_Snaps, 1)
            assign(Buffer_Cycles, 0)
            assign(i, IO1)

            ## PARAMETERS UPDATE ##
            with if_(i > 0):
                pause()
            with while_(i > 0):
                ## Boolean variables control: ##
                with if_(i == IOP.MOT_SWITCH_ON.value):
                    update_frequency("MOT_AOM_0", Config.IF_AOM_MOT)
                with if_(i == IOP.MOT_SWITCH_OFF.value):
                    update_frequency("MOT_AOM_0", Config.IF_AOM_MOT_OFF)
                with if_(i == IOP.MAGNETIC_FOUNTAIN_ON.value):
                    assign(Magnetic_Fountain_ON, IO2)
                ## AntiHelmholtz control ##
                with if_(i == IOP.ANTIHELMHOLTZ_DELAY.value):
                    assign(antihelmholtz_delay, IO2)
                ## MOT variables control ##
                with if_(i == IOP.MOT_DURATION.value):  # Live control over the
                    assign(MOT_Repetitions, IO2)
                with if_(i == IOP.POST_MOT_DELAY.value):  # Live control over the post MOT delay
                    assign(post_MOT_delay, IO2)
                    #assign(params[IOP.POST_MOT_DELAY], IO2)  # @@@ - 4  TODO: Generalize this! assign(p[i], IO2)

                ## PGC variables control ##
                with if_(i == IOP.PGC_DURATION.value):  # Live control over the PGC duration
                    assign(pgc_duration, IO2)
                with if_(i == IOP.PGC_PREP_DURATION.value):  # Live control over the preparation time of the PGC
                    assign(pgc_prep_duration, IO2)
                with if_(i == IOP.PGC_PULSE_DURATION_0.value):  # Live control over the final amplitude of the PGC AOM 0
                    assign(pgc_pulse_duration_0, IO2)
                with if_(i == IOP.PGC_PULSE_DURATION_MINUS.value):  # Live control over the final amplitude of the PGC AOM -
                    assign(pgc_pulse_duration_minus, IO2)
                with if_(i == IOP.PGC_PULSE_DURATION_PLUS.value):  # Live control over the final amplitude of the PGC AOM +
                    assign(pgc_pulse_duration_plus, IO2)
                with if_(i == IOP.PGC_BEAMS_0_OFF_DURATION.value):
                    assign(pgc_beams_0_off_duration, IO2)

                ## Fountain variables control: ##
                with if_(i == IOP.FOUNTAIN_DURATION.value):  # Live control over the fountain duration
                    assign(fountain_duration, IO2)
                with if_(i == IOP.FOUNTAIN_PREP_TIME.value):  # Live control over the preparation time of the Fountain
                    assign(fountain_prep_time, IO2)
                with if_(i == IOP.FOUNTAIN_PULSE_DURATION_0.value):  # Live control over the final amplitude of the fountain AOM 0
                    assign(fountain_pulse_duration_0, IO2)
                with if_(i == IOP.FOUNTAIN_PULSE_DURATION_MINUS.value):  # Live control over the final amplitude of the fountain AOM -
                    assign(fountain_pulse_duration_minus, IO2)
                with if_(i == IOP.FOUNTAIN_PULSE_DURATION_PLUS.value):  # Live control over the final amplitude of the fountain AOM +
                    assign(fountain_pulse_duration_plus, IO2)
                with if_(i == IOP.FOUNTAIN_FINAL_DELTA_FREQ.value):  # Live control over the fountain final frequency of the MOT + and - AOMs
                    assign(fountain_aom_chirp_rate, IO2)

                ## Measurement variables control: ##
                with if_(i == IOP.TRIGGER_DELAY.value):
                    assign(Trigger_delay, IO2)
                with if_(i == IOP.PREPULSE_DURATION.value):
                    assign(PrePulse_duration, IO2)
                with if_(i == IOP.MAGNETIC_FOUNTAIN_DURATION.value):
                    assign(magnetic_fountain_duration, IO2)
                with if_(i == IOP.AOM_0_ATT_3RD_STAGE.value):
                    assign(AOM_0_att_3rd_stage, IO2)

                with if_(i == IOP.PULSE_1_DURATION.value):
                    assign(Pulse_1_duration, IO2)
                with if_(i == IOP.PULSE_1_DECAY_DURATION.value):
                    assign(Pulse_1_decay_time, IO2)
                with if_(i == IOP.N_SNAPS.value):  # The number of frames(snaps) to take for a video with growing snap time intervals
                    assign(N_Snaps, IO2)
                with if_(i == IOP.BUFFER_CYCLES.value):
                    assign(Buffer_Cycles, IO2)

                ## OD and N_atoms measuring variable control ##
                with if_(i == IOP.OD_FS_START.value):  # Live control of the Free Space OD measurement start time
                    assign(OD_FS_start, IO2)
                with if_(i == IOP.OD_FS_PULSE_DURATION.value):  # Live control of the Free Space OD measurement pulses duration
                    assign(OD_FS_pulse_duration, IO2)
                with if_(i == IOP.OD_FS_PULSES_SPACING.value):  # Live control of the Free Space OD measurement wait duration between the 2 pulses
                    assign(OD_FS_pulses_spacing, IO2)
                with if_(i == IOP.OD_DELAY.value):  # Live control of the SNSPDs OD measurement start time
                    assign(OD_Delay, IO2)
                with if_(i == IOP.OD_SNSPDS_DURATION.value):  # Live control of the SNSPDs OD measurement duration
                    assign(Rep, IO2)
                with if_(i == IOP.OD_CONTINUOUS_ATTENUATION.value):  # Turn on continuous OD and attenuate factor
                    assign(OD_attenuation, IO2)
                with if_(i == IOP.DEPUMP_START.value):  # Live control of the Depump measurement start time
                    assign(Depump_start, IO2)
                with if_(i == IOP.DEPUMP_PULSE_DURATION.value):  # Live control of the Depump measurement pulses duration
                    assign(Depump_pulse_duration, IO2)
                with if_(i == IOP.DEPUMP_PULSES_SPACING.value):  # Live control of the Depump measurement wait duration between the 2 pulses
                    assign(Depump_pulses_spacing, IO2)
                with if_(i == IOP.SHUTTER_OPENING_TIME.value):  # Live control of the delay due to shutter opening time.
                    assign(shutter_open_time, IO2)

                ## MW spectroscopy control ##
                # with if_(i == 61):  # Live control of the MW spectroscopy MW frequency
                #     assign(MW_freq, IO2)
                # with if_(i == 62):  # Live control of the MW spectroscopy MW pulse length
                #     assign(pulse_length_MW, IO2)
                # with if_(i == 63):  # Live control of the MW spectroscopy OD pulse length
                #     assign(pulse_length_OD, IO2)
                # with if_(i == 64):  # Live control of the MW spectroscopy MW pulse repetitions
                #     assign(Repetitions, IO2)
                # with if_(i == 65):  # Live control of the MW spectroscopy MW pulse frequency change
                #     assign(Delta_f, IO2)

                pause()
                assign(i, IO1)

        # TODO: remove the following - not in use by the generic MOT sequence
        # with stream_processing():
        #     ON_counts_st1.buffer(obj.rep).save('Det1_Counts')
        #     ON_counts_st2.buffer(obj.rep).save('Det2_Counts')
        #     ON_counts_st3.buffer(obj.rep).save('Det3_Counts')
        #     ON_counts_st6.buffer(obj.rep).save('Det6_Counts')
        #     ON_counts_st7.buffer(obj.rep).save('Det7_Counts')
        #     ON_counts_st8.buffer(obj.rep).save('Det8_Counts')
        #     (tt_st_1 + rep_st).buffer(obj.vec_size * obj.rep).save('Det1_Probe_TT')
        #     (tt_st_2 + rep_st).buffer(obj.vec_size * obj.rep).save('Det2_Probe_TT')
        #     (tt_st_3 + rep_st).buffer(obj.vec_size * obj.rep).save('Det3_Probe_TT')
        #     (tt_st_6 + rep_st).buffer(obj.vec_size * obj.rep).save('Det6_Probe_TT')
        #     (tt_st_7 + rep_st).buffer(obj.vec_size * obj.rep).save('Det7_Probe_TT')
        #     (tt_st_8 + rep_st).buffer(obj.vec_size * obj.rep).save('Det8_Probe_TT')
        #     FLR_st.save('FLR_measure')
        #     AntiHelmholtz_ON_st.save("antihelmholtz_on")

    job = qm.execute(opx_control_prog, flags=['auto-element-thread'])

    return job
