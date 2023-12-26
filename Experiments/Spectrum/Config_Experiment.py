import numpy as np
from scipy import signal
import json
import math
from Utilities.Utils import Utils

Spectrum_Exp_Gaussian_samples = [0] * 10 + (signal.gaussian(220, std=30) * 0.4).tolist() + [0] * 10
readout_pulse_spectrum_len = int(len(Spectrum_Exp_Gaussian_samples) * 1000 * 21 * 4)

# TODO - need to repalce all 'QRAM' to 'Spectrum', or a generic name (also rename parameters)
def QRAM_Exp_Gaussian_samples(sprint_pulse_len=110, det_pulse_len=30, det_pulses_amp=[0.4]*6, sprint_pulses_amp=[0.4]*4,
                              num_between_zeros=10, num_init_zeros=12, num_mid_zeros=12, num_fin_zeros=0):
    qram_exp_gaussian_samples = [0] * (num_init_zeros - 2 + 42 + 16)  # -2 for echos from south
    for n in det_pulses_amp:
        qram_exp_gaussian_samples += (signal.gaussian(det_pulse_len, std=(det_pulse_len * 0.5 / 2.355)) * n).tolist() + [0] * (num_between_zeros - 3) # -3 for echos from south
    # qram_exp_gaussian_samples += [0] * (num_mid_zeros - 10) # due to unresolved reflections +40 for S echos
    qram_exp_gaussian_samples = qram_exp_gaussian_samples[:-(num_between_zeros - 3)] + [0] * (num_between_zeros - 3 + num_mid_zeros - 16) # due to unresolved reflections +40 for S echos
    # for m in sprint_pulses_amp:
    #     qram_exp_gaussian_samples += (signal.gaussian((sprint_pulse_len-4), std=((sprint_pulse_len-4) / 2.355)) * m).tolist() + [0] * (num_between_zeros + 4)
    # # qram_exp_gaussian_samples += [0] * num_fin_zeros
    # qram_exp_gaussian_samples = qram_exp_gaussian_samples[:-(num_between_zeros + 4)] + [0] * (num_between_zeros + 4 + num_fin_zeros - 16) # -14 due to S echos
    return qram_exp_gaussian_samples


def QRAM_Exp_Square_samples(amp=0.45, sprint_pulse_len=110, det_pulse_len=30, det_pulses_amp=[0.4]*6, sprint_pulses_amp=[0.4]*4,
                            num_between_vals=10, num_init_val=12, num_mid_val=12, num_fin_val=0):

    qram_exp_gaussian_samples = [det_pulses_amp[0] * amp] * (num_init_zeros - 2 + 42 + 16)

    for n in range(len(det_pulses_amp)-1):
        qram_exp_gaussian_samples += [det_pulses_amp[n] * amp] * det_pulse_len + \
                                     [det_pulses_amp[n] * det_pulses_amp[n+1] * amp] * (num_between_vals - 3)
    qram_exp_gaussian_samples += [det_pulses_amp[-1] * amp] * det_pulse_len + \
                                 [det_pulses_amp[-1] * sprint_pulses_amp[0] * amp] * (num_between_vals - 16 - 3)

    qram_exp_gaussian_samples += [sprint_pulses_amp[0] * amp] * (num_mid_val)  # due to unresolved reflections

    # for m in range(len(sprint_pulses_amp)-1):
    #     qram_exp_gaussian_samples += [sprint_pulses_amp[m] * amp] * (sprint_pulse_len - 4) + \
    #                                  [sprint_pulses_amp[m] * sprint_pulses_amp[m + 1] * amp] * (num_between_vals + 4)
    # qram_exp_gaussian_samples += [sprint_pulses_amp[-1] * amp] * (sprint_pulse_len - 4) + \
    #                              [sprint_pulses_amp[-1] * det_pulses_amp[0] * amp] * (num_between_vals + 4)
    #
    # qram_exp_gaussian_samples += [det_pulses_amp[0] * amp] * num_fin_val

    return qram_exp_gaussian_samples

def QRAM_Exp_samples(delta=240, pulse_len=10000000):
    QRAM_exp_samples = []
    for n in range(int(pulse_len//(3*delta))):
        QRAM_exp_samples += [0.4]*delta + [0]*2*delta
    return QRAM_exp_samples

def get_pulses_location_in_seq(delay, seq, smearing = 0):
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
    seq_indx = np.where(seq_filter > 0)[0] + 1
    pulses_loc = []
    if seq_indx.any():
        start_indx = seq_indx[0]
        for i in range(1, len(seq_indx)):
            if (seq_indx[i] - seq_indx[i-1]) > 1:
                if (start_indx > 1) and pulses_loc:
                    pulses_loc.append((0, start_indx - 1 - pulses_loc[-1][-1] - int(smearing)))
                pulses_loc.append((1, seq_indx[i-1] - start_indx + 1 + int(smearing)))
                start_indx = seq_indx[i]
        if (start_indx > 1) and pulses_loc:
            pulses_loc.append((0, start_indx - 1 - pulses_loc[-1][-1] - int(smearing)))
        pulses_loc.append((1, seq_indx[-1] - start_indx + 1 + int(smearing)))
    pulses_loc.append((0, 0))
    return pulses_loc

controller = 'con1'

# Parameters:
# delays [ns]:
detector_delays = [26, 30, 33, 27, 8, 8, 8, 0]  # For detectors [1,2,3,9,15,6,7,8]
AOM_Late_delay = 670
AOM_Early_delay = 525
AOM_S_to_N_delay = 40
# time tags vector size
# parameters of sizes
vec_size = 8000
num_of_detectors = 2
opx_max_per_window = vec_size*num_of_detectors/2

# Pulse_durations
readout_pulse_len = 10000000
# MOT_pulse_len = 2e6
MOT_pulse_len = 2e6
Short_pulse_len = 40
PGC_pulse_len = 40
Fountain_pulse_len = 40
FreeFall_pulse_len = 40
Probe_pulse_len = 5e6
Depump_pulse_len = 40
Measuring_pulse_len = 40
OD_pulse_len = 2e6
Repump_pulse_len = 10e6
north_const_pulse_len = 500
south_const_pulse_len = 500
analyzer_const_pulse_len = 500
Detection_pulse_len = 1000
Flash_pulse_len = 10000
spectrum_pulse = 240
frequency_sweep_duration = 500

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
IF_AOM_N = 129e6
# IF_AOM_S = 89.2368e6
# IF_AOM_S = 129.2368e6
IF_AOM_S = 129e6
# IF_AOM_LO = 89.2368e6
IF_AOM_ANCILLA = 129.2368e6 + 110e6/2
# IF_AOM_LO = 129.2368e6
IF_AOMs_MZ = 110e6
IF_AOM_SigmaPlus = 114.58e6
IF_AOM_SigmaMinus = 114.58e6
IF_AOM_Pi = 75.34e6
IF_CRUS_pulser = 125e6
IF_AOM_Spectrum = 133.325e6/2

IF_Divert = 20e6
# IF_AOM_Analyzer = np.abs(IF_AOM_N - IF_AOM_S) * 2

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

num_of_det_pulses_S = 4
# num_of_det_pulses_S = 2
num_of_sprint_pulses_S = 1
# num_of_sprint_pulses_S = 0  # For only det pulses sequence
num_of_det_pulses_N = 4
# num_of_det_pulses_N = 2
num_of_sprint_pulses_N = 3
# num_of_sprint_pulses_N = 0  # For only det pulses sequence

# parameters for window len
# efficiency = 0.5 # the efficiency of the system
# efficiency = 0.29 # the efficiency of the system
efficiency = 0.168  # the efficiency of the system with sagnac
num_of_photons_per_sequence_S = num_of_photons_det_pulses * num_of_det_pulses_S + num_of_photons_sprint_pulses * num_of_sprint_pulses_S
num_of_photons_per_sequence_N = num_of_photons_det_pulses * num_of_det_pulses_N + num_of_photons_sprint_pulses * num_of_sprint_pulses_N


det_pulse_len = 40
sprint_pulse_len = 105
num_between_zeros = 20

# For general sequence pulses shape
num_init_zeros = 10  # For only det pulses sequence
num_mid_zeros = 10
num_fin_zeros = 0  # For only det pulses sequence
det_pulse_amp_General = [1, 1, 1, 1, 1, 1, 1, 1]
sprint_pulse_amp_General = [1, 1, 1, 1]


num_init_zeros_S = 10  # For only det pulses sequence
num_mid_zeros_S = 10
num_fin_zeros_S = 0  # For only det pulses sequence
det_pulse_amp_S = [0, 0, 0, 0, 0, 0, 0, 0]
# sprint_pulse_amp_S = [0, 0.45, 0, 0.45]
# For pulse sync
# det_pulse_amp_S = [0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45]
# sprint_pulse_amp_S = [0, 0, 0, 0]
# sprint_pulse_amp_S = [1, 1, 1, 1]
# For Bell |(0 + 1)c, 1t>
# det_pulse_amp_S = [0.45, 0, 0.45, 0, 0.45, 0, 0.45, 0]
# sprint_pulse_amp_S = [0.075, 0, 0, 0]
# |0c, (0 + 1)t>
# det_pulse_amp_S = [0, 0, 0.45, 0, 0.45, 0, 0.45, 0]
# det_pulse_amp_S = [0, 0.45, 0, 0.45, 0, 0.45, 0, 0]
sprint_pulse_amp_S = [0, 0, 0, 0]
# |1c, (0 + 1)t>
# det_pulse_amp_S = [0, 0.45, 0, 0.45, 0, 0.45, 0, 0.45]
# sprint_pulse_amp_S = [0, 0.07, 0, 0.07]
# sprint_pulse_amp_S = [0, 0.15, 0, 0.15]
# det_pulse_amp_S = [0.155, 0, 0.155, 0, 0.155, 0, 0.155, 0]
# sprint_pulse_amp_S = [0, 0.078, 0, 0.078]
# sprint_pulse_amp_S = [0.095, 0, 0, 0]
# sprint_pulse_amp_S = [0.078, 0, 0.078, 0]

# Sprint_Exp_Gaussian_samples_S = Sprint_Exp_Gaussian_samples(sprint_pulse_len=sprint_pulse_len,
#                                                             det_pulse_len=det_pulse_len,
#                                                             det_pulses_amp=det_pulse_amp_S,
#                                                             sprint_pulses_amp=sprint_pulse_amp_S,
#                                                             num_between_zeros=num_between_zeros,
#                                                             num_init_zeros=num_init_zeros_S,
#                                                             num_fin_zeros=num_fin_zeros_S)

num_init_zeros_N = 10  # For only det pulses sequence
num_mid_zeros_N = 10
num_fin_zeros_N = 0  # For only det pulses sequence
# det_pulse_amp_N = [0.45, 0, 0.45, 0, 0.45, 0, 0.45, 0]
# det_pulse_amp_N = [0, 0.45, 0, 0.45, 0, 0.45, 0, 0.45]
# sprint_pulse_amp_N = [0.45, 0, 0.45, 0]
# sprint_pulse_amp_N = [0, 0.105, 0, 0.105]
# sprint_pulse_amp_N = [0, 0.25, 0, 0.25]
# For pulse sync
det_pulse_amp_N = [0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45]
# det_pulse_amp_N = [0, 0, 0, 0, 0, 0, 0, 0]
# det_pulse_amp_N = [1, 1, 1, 1, 1, 1, 1, 1]
# sprint_pulse_amp_N = [0.45, 0.45, 0.45, 0]
# sprint_pulse_amp_N = [0.45, 0, 0.45, 0]
# sprint_pulse_amp_N = [1, 1, 1, 1]
# For Bell |(0 + 1)c, 1t>
# det_pulse_amp_N = [0, 0.45, 0, 0.45, 0, 0.45, 0, 0.45]
# sprint_pulse_amp_N = [0, 0, 0, 0.085]
# |0c, (0 + 1)t>
# det_pulse_amp_N = [0, 0.45, 0, 0.45, 0, 0.45, 0, 0.45]
# det_pulse_amp_N = [0.45, 0, 0.45, 0, 0.45, 0, 0.45, 0]
# det_pulse_amp_N = [0, 0.155, 0, 0.155, 0, 0.155, 0, 0.155]
# sprint_pulse_amp_N = [0, 0.078, 0, 0.078]
# # |1c, (0 + 1)t>
# det_pulse_amp_N = [0.45, 0, 0.45, 0, 0.45, 0, 0.45, 0]
# sprint_pulse_amp_N = [0.095, 0, 0, 0]
sprint_pulse_amp_N = [0, 0, 0, 0]

# Sprint_Exp_Gaussian_samples_N = Sprint_Exp_Gaussian_samples(sprint_pulse_len=sprint_pulse_len,
#                                                             det_pulse_len=det_pulse_len,
#                                                             det_pulses_amp=det_pulse_amp_N,
#                                                             sprint_pulses_amp=sprint_pulse_amp_N,
#                                                             num_between_zeros=num_between_zeros,
#                                                             num_init_zeros=num_init_zeros_N,
#                                                             # num_mid_zeros=num_mid_zeros_N,
#                                                             num_fin_zeros=num_fin_zeros_N)

num_init_zeros_Ancilla = 10  # For only det pulses sequence
num_mid_zeros_Ancilla = 10
num_fin_zeros_Ancilla = 0  # For only det pulses sequence

# For pulse sync
# det_pulse_amp_Ancilla = [0, 0, 0, 0, 0, 0, 0, 0]
# sprint_pulse_amp_Ancilla = [0, 0, 0, 0]
# For Bell |(0 + 1)c, 1t>
# det_pulse_amp_Ancilla = [0, 0, 0, 0, 0, 0, 0, 0]
# sprint_pulse_amp_Ancilla = [0, 0, 0.085, 0]
# |0c, (0 + 1)t>
# det_pulse_amp_Ancilla = [0, 0, 0, 0, 0, 0, 0, 0]
# sprint_pulse_amp_Ancilla = [0, 0, 0.085, 0]
# # |1c, (0 + 1)t>
det_pulse_amp_Ancilla = [0, 0, 0, 0, 0, 0, 0, 0]
sprint_pulse_amp_Ancilla = [0, 0, 0, 0]

# Qram

QRAM_Exp_Gaussian_samples_General = QRAM_Exp_Gaussian_samples(sprint_pulse_len=sprint_pulse_len,
                                                              det_pulse_len=det_pulse_len,
                                                              det_pulses_amp=det_pulse_amp_General,
                                                              sprint_pulses_amp=sprint_pulse_amp_General,
                                                              num_between_zeros=num_between_zeros,
                                                              num_init_zeros=num_init_zeros,
                                                              num_mid_zeros=num_mid_zeros,
                                                              num_fin_zeros=num_fin_zeros)

QRAM_Exp_Gaussian_samples_Ancilla = QRAM_Exp_Gaussian_samples(sprint_pulse_len=sprint_pulse_len,
                                                              det_pulse_len=det_pulse_len,
                                                              det_pulses_amp=det_pulse_amp_Ancilla,
                                                              sprint_pulses_amp=sprint_pulse_amp_Ancilla,
                                                              num_between_zeros=num_between_zeros,
                                                              num_init_zeros=num_init_zeros_Ancilla,
                                                              num_mid_zeros=num_mid_zeros_Ancilla,
                                                              num_fin_zeros=num_fin_zeros_Ancilla)

QRAM_Exp_Gaussian_samples_S = QRAM_Exp_Gaussian_samples(sprint_pulse_len=sprint_pulse_len,
                                                        det_pulse_len=det_pulse_len,
                                                        det_pulses_amp=det_pulse_amp_S,
                                                        sprint_pulses_amp=sprint_pulse_amp_S,
                                                        num_between_zeros=num_between_zeros,
                                                        num_init_zeros=num_init_zeros_S,
                                                        num_mid_zeros=num_mid_zeros_S,
                                                        num_fin_zeros=num_fin_zeros_S)

QRAM_Exp_Gaussian_samples_N = QRAM_Exp_Gaussian_samples(sprint_pulse_len=sprint_pulse_len,
                                                        det_pulse_len=det_pulse_len,
                                                        det_pulses_amp=det_pulse_amp_N,
                                                        sprint_pulses_amp=sprint_pulse_amp_N,
                                                        num_between_zeros=num_between_zeros,
                                                        num_init_zeros=num_init_zeros_N,
                                                        num_mid_zeros=num_mid_zeros_N,
                                                        num_fin_zeros=num_fin_zeros_N)

num_between_vals = 20
num_init_val_Early = 10  # For only det pulses sequence
num_mid_val_Early = 10
num_fin_val_Early = 0  # For only det pulses sequence
Pulses_Amp = 0.45
Pulses_Amp_Early = 0.495
# For pulse sync
# det_pulse_amp_Early = [0, 0, 0, 0, 0, 0, 1, 1]
# sprint_pulse_amp_Early = [1, 0, 0, 0]
# For Bell |(0 + 1)c, 1t>s
# det_pulse_amp_Early = [0, 0, 0, 0, 0, 0, 1, 1]
# sprint_pulse_amp_Early = [1, 0, 0, 0]
# |0c, (0 + 1)t>
# det_pulse_amp_Early = [0, 0, 0, 0, 0, 0, 0, 0]
# sprint_pulse_amp_Early = [1, 1, 0, 0]
# |1c, (0 + 1)t>
det_pulse_amp_Early = [0, 0, 0, 0, 0, 0, 0, 0]
# sprint_pulse_amp_Early = [1, 1, 0, 0]
sprint_pulse_amp_Early = [0, 0, 0, 0]


QRAM_Exp_Square_samples_Early = QRAM_Exp_Square_samples(amp=Pulses_Amp_Early,
                                                        sprint_pulse_len=sprint_pulse_len,
                                                        det_pulse_len=det_pulse_len,
                                                        det_pulses_amp=det_pulse_amp_Early,
                                                        sprint_pulses_amp=sprint_pulse_amp_Early,
                                                        num_between_vals=num_between_vals,
                                                        num_init_val=num_init_val_Early,
                                                        num_mid_val=num_mid_val_Early,
                                                        num_fin_val=num_fin_val_Early)
QRAM_Exp_Square_samples_Early_delayed = np.roll(QRAM_Exp_Square_samples_Early, AOM_Early_delay)

num_init_val_Late = 10  # For only det pulses sequence
num_mid_val_Late = 10
num_fin_val_Late = 0  # For only det pulses sequence
Pulses_Amp_Late = 0.35  # For balancing with AOM Early (since added the switch)
# Pulses_Amp_Late = 0.45  # For balancing with AOM Early (since added the switch)
# For pulse sync
det_pulse_amp_Late = [1, 1, 1, 1, 1, 1, 1, 1]
# sprint_pulse_amp_Late = [0, 0, 0, 0]
sprint_pulse_amp_Late = [1, 1, 1, 1]
# sprint_pulse_amp_Late = [0, 1, 1, 0]
# For Bell |(0 + 1)c, 1t>
# det_pulse_amp_Late = [0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0, 0]
# sprint_pulse_amp_Late = [0, 0, 0, 0]
# |0c, (0 + 1)t>
# det_pulse_amp_Late = [0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8]
# sprint_pulse_amp_Late = [0, 0, 1, 1]
# # |1c, (0 + 1)t>
# det_pulse_amp_Late = [1, 1, 1, 1, 1, 1, 1, 1]
# sprint_pulse_amp_Late = [0, 0, 1, 1]
# sprint_pulse_amp_Late = [1, 1, 1, 1]

QRAM_Exp_Square_samples_Late = QRAM_Exp_Square_samples(amp=Pulses_Amp_Late,
                                                       sprint_pulse_len=sprint_pulse_len,
                                                       det_pulse_len=det_pulse_len,
                                                       det_pulses_amp=det_pulse_amp_Late,
                                                       sprint_pulses_amp=sprint_pulse_amp_Late,
                                                       num_between_vals=num_between_vals,
                                                       num_init_val=num_init_val_Late,
                                                       num_mid_val=num_mid_val_Late,
                                                       num_fin_val=num_fin_val_Late)
QRAM_Exp_Square_samples_Late_delayed = np.roll(QRAM_Exp_Square_samples_Late, AOM_Late_delay)

num_init_val_FS = 10  # For only det pulses sequence
num_mid_val_FS = 10
num_fin_val_FS = 0  # For only det pulses sequence
# For pulse sync
# det_pulse_FS = [1, 1, 1, 1, 1, 1, 1, 1]
# sprint_pulse_FS = [1, 1, 1, 1]
# For Bell |(0 + 1)c, 1t>
# det_pulse_amp_FS = [1, 1, 1, 1, 1, 1, 1, 1]
# sprint_pulse_amp_FS = [1, 0, 0, 0]
# |0c, (0 + 1)t>
# det_pulse_amp_FS = [1, 1, 1, 1, 1, 1, 1, 1]
# sprint_pulse_FS = [1, 1, 1, 1]
# # |1c, (0 + 1)t>
det_pulse_amp_FS = [1, 1, 1, 1, 1, 1, 1, 1]
sprint_pulse_FS = [1, 1, 1, 1]

QRAM_Exp_Square_samples_FS = QRAM_Exp_Square_samples(amp=Pulses_Amp,
                                                     sprint_pulse_len=sprint_pulse_len,
                                                     det_pulse_len=det_pulse_len,
                                                     det_pulses_amp=det_pulse_amp_FS,
                                                     sprint_pulses_amp=sprint_pulse_FS,
                                                     num_between_vals=num_between_vals,
                                                     num_init_val=num_init_val_FS,
                                                     num_mid_val=num_mid_val_FS,
                                                     num_fin_val=num_fin_val_FS)

QRAM_Exp_digital_samples_FS = get_pulses_location_in_seq(delay=0,
                                                         seq=QRAM_Exp_Square_samples_FS,
                                                         smearing=0)


MZ_delay = int(len(QRAM_Exp_Gaussian_samples_N) / 4)
# Pulses_Amp_balance = 0.18
Pulses_Amp_balance = 0.3
AOM_risetime = 120
AOM_risetime_pulsers = 130
MZ_balancing_seq_rep = 40
# QRAM_MZ_balance_pulse_North = ([Pulses_Amp*0.5] * (MZ_delay - AOM_risetime_pulsers) + [0] * MZ_delay + [0] * AOM_risetime_pulsers) * MZ_balancing_seq_rep
QRAM_MZ_balance_pulse_North = ([0] * (MZ_delay - AOM_risetime_pulsers) + [Pulses_Amp_balance] * (AOM_risetime_pulsers - 20) + [0] * 20
                               + [0] * (MZ_delay - AOM_risetime_pulsers) + [Pulses_Amp_balance] * (AOM_risetime_pulsers - 20) + [0] * 20) * \
                              MZ_balancing_seq_rep
QRAM_MZ_balance_pulse_South = ([0] * (MZ_delay - AOM_risetime_pulsers) + [0] * (AOM_risetime_pulsers - 20) + [0] * 20
                               + [0] * (MZ_delay - AOM_risetime_pulsers) + [0] * (AOM_risetime_pulsers - 20) + [0] * 20) * \
                              MZ_balancing_seq_rep
QRAM_MZ_balance_pulse_Early = ([Pulses_Amp_Early] * MZ_delay + [0] * MZ_delay) * MZ_balancing_seq_rep
QRAM_MZ_balance_pulse_Early_delayed = np.roll(QRAM_MZ_balance_pulse_Early, AOM_Early_delay-10)
# QRAM_MZ_balance_pulse_Early_delayed = np.roll(QRAM_MZ_balance_pulse_Early, AOM_Early_delay-10+AOM_S_to_N_delay)
QRAM_MZ_balance_pulse_Late = ([0] * MZ_delay + [Pulses_Amp_Late] * (MZ_delay)) * MZ_balancing_seq_rep
QRAM_MZ_balance_pulse_Late_delayed = np.roll(QRAM_MZ_balance_pulse_Late, AOM_Late_delay-10)
# QRAM_MZ_balance_pulse_Late_delayed = np.roll(QRAM_MZ_balance_pulse_Late, AOM_Late_delay-10+AOM_S_to_N_delay)

QRAM_Exp_Gaussian_samples_N = ([0] * (MZ_delay - AOM_risetime_pulsers) + [0] * (AOM_risetime_pulsers - 20) + [0] * 20
                               + [0] * (MZ_delay - AOM_risetime_pulsers) + [0] * (AOM_risetime_pulsers - 20) + [0] * 20) * 2
# QRAM_Exp_Gaussian_samples_N = ([0] * (MZ_delay - AOM_risetime_pulsers) + [Pulses_Amp_balance] * (AOM_risetime_pulsers - 20) + [0] * 20
#                                + [0] * (MZ_delay - AOM_risetime_pulsers) + [Pulses_Amp_balance] * (AOM_risetime_pulsers - 20) + [0] * 20) * 2
QRAM_Exp_Gaussian_samples_S = ([Pulses_Amp_balance] * (MZ_delay - AOM_risetime_pulsers) + [Pulses_Amp_balance] * (AOM_risetime_pulsers - 20) + [Pulses_Amp_balance] * 20
                               + [Pulses_Amp_balance] * (MZ_delay - AOM_risetime_pulsers) + [Pulses_Amp_balance] * (AOM_risetime_pulsers - 20) + [Pulses_Amp_balance] * 20) * 2
# QRAM_Exp_Gaussian_samples_S = ([0] * (MZ_delay - AOM_risetime_pulsers) + [Pulses_Amp_balance] * (AOM_risetime_pulsers - 20) + [0] * 20
#                                + [0] * (MZ_delay - AOM_risetime_pulsers) + [Pulses_Amp_balance] * (AOM_risetime_pulsers - 20) + [0] * 20) * 2
QRAM_Exp_Square_samples_Late = ([Pulses_Amp_Late] * MZ_delay + [Pulses_Amp_Late] * MZ_delay) * 2
# QRAM_Exp_Square_samples_Early = ([Pulses_Amp_Early] * MZ_delay + [0] * MZ_delay) * 2
QRAM_Exp_Square_samples_Early = ([0] * MZ_delay + [0] * MZ_delay) * 2
QRAM_Exp_Square_samples_Late_delayed = np.roll(QRAM_Exp_Square_samples_Late, AOM_Late_delay-10-AOM_S_to_N_delay)
QRAM_Exp_Square_samples_Early_delayed = np.roll(QRAM_Exp_Square_samples_Early, AOM_Early_delay-10-AOM_S_to_N_delay)


# readout_pulse_sprint_len_N = math.ceil(((opx_max_per_window/1.5)/(efficiency*1e6*num_of_photons_per_sequence_N))*len(Sprint_Exp_Gaussian_samples_N))*1e6# [ns] length of the measurment window for North, the 4's are for division in 4
readout_pulse_sprint_len_N = 10*1e6# [ns] length of the measurment window for North, the 4's are for division in 4
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
Depump_pulse_samples = np.convolve(square_samples, Utils.gaussian(x, 600, 20))[900 // 4 * 4:1500 // 4 * 4]
OD_pulse_samples = np.convolve(square_samples, Utils.gaussian(x, 600, 20))[900 // 4 * 4:1500 // 4 * 4]

square_samples2 = ([0] * 400 + [0.8] * 400 + [0] * 400)
xx = np.linspace(0, len(square_samples2), len(square_samples2))
Gaussian_pulse_samples = (signal.gaussian(500, std=(300 / 2.355)) * 0.2).tolist()
Gaussian_pulse_samples2 = (signal.gaussian(500, std=(300 / 2.355)) * 0.2).tolist()
# Gaussian_pulse_samples = Utils.gauss_adaptive(0.45, 500)

## Attenuators (global) for AOMs 0, + & - of MOT sequence
# Factor 0.0-1.0
AOM_0_Attenuation = 1 # 0.85 seems to be about right, for 0.8 the slopes are exactly the same
AOM_Plus_Attenuation = 0.42
# AOM_Plus_Attenuation = 0.38
AOM_Minus_Attenuation = 1 # 0.85 seems to be about right, for 0.8 the slopes are exactly the same

## For Homodyne ##
LO_pulse_samples = ([0.4] * 20) + [0] * 100


#------------------------------------------------------
# NOTES - OPX Limitations:
# - In 'elements', OPX allows only N intermediate channels
# - Elements: avoid having multiple elements playing on the same port at the same time - may cause problems to the OPX
# - Pulses: must be a multiplication of 4
# - Waveforms: sample (amplitude) must not exceed 0.5 v (OPX will not alert, but perform modulu)
#------------------------------------------------------
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
                5: {'offset': +0.0},  # AOM 2-2/3' (Depump/OD)
                6: {'offset': +0.0},  # AOM early
                7: {'offset': +0.0},  # AOM late
                8: {'offset': +0.0},  # AOM Ancilla
                9: {'offset': +0.0},  # AOM N
                10: {'offset': +0.0}, # AOM S
            },

            'digital_outputs': {
                1: {},  # Trigger TOP1 - QuadRF
                2: {},  # OD/Depump Switch
                3: {},  # Switch ON/OFF anti-Helmholtz coils (for MOT)
                4: {},  #
                5: {},  # Switch Shutters to SNSPDs
                6: {},  #
                7: {},  # Measurement Trigger
                8: {},  #
                9: {},  # Switch ON/OFF Helmholtz coils (for Super-SPRINT)
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
            'digitalInputs': {
                "AWG_Switch": {
                    "port": (controller, 10),
                    "delay": 210,
                    "buffer": 0,
                },
            },
            'operations': {
                'MOT': "MOT_lock",
                'MOT_with_EOM': "MOT_with_EOM_pulses",
                'Linear': "Linear_pulse",
                'Const': "Const_pulse",
            },
            'intermediate_frequency': IF_AOM_MOT,
        },

        "MOT_AOM_-": {
            'singleInput': {
                "port": (controller, 3)
            },
            'digitalInputs': {
                "RedPitaya": {
                    "port": (controller, 4),
                    "delay": 0,
                    "buffer": 0,
                },
            },
            'operations': {
                # 'MOT': "MOT_lock",
                'MOT': "MOT_with_Trigger",
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
                # "Shutter_Switch": {
                #     "port": (controller, 5),
                #     "delay": 0,
                #     "buffer": 0,
                # },
            },
            'operations': {
                'OD': "OD_pulse",
                'OD_FS': "OD_FS_pulse",
                'Depump': "Depump_pulse",
            },
            'intermediate_frequency': IF_AOM_OD,
        },

        # "AOM_2-2'": {
        #     'singleInput': {
        #         "port": (controller, 5)
        #     },
        #     'operations': {
        #         'Depump': "Depump_pulse",
        #     },
        #     'intermediate_frequency': IF_AOM_Depump,
        # },

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

        "Dig_detectors": {
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
                "out1": (controller, 1),
                "out2": (controller, 2),
                "out3": (controller, 3),
                "out4": (controller, 4),
                "out5": (controller, 5),
                "out6": (controller, 6),
                "out7": (controller, 7),
                "out8": (controller, 8),
                "out9": (controller, 9),
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
            # Fake Input:
            'digitalInputs': {
                "Helmholtz": {
                    "port": (controller, 4),
                    "delay": 0,
                    "buffer": 0,
                },
            },
            'operations': {
                'ZeemanSplit': "Zeeman_on",
                'ZeemanOFF': "Zeeman_off",
            },
        },

        "AOM_Early": {
            "singleInput": {
                "port": (controller, 6),
            },
            'digitalInputs': {
                "APD_Switch": {
                    "port": (controller, 6),
                    "delay": 0,
                    "buffer": 0,
                },
                # "AWG_Switch": {
                #     "port": (controller, 10),
                #     # "delay": 176, # OPX control EOM
                #     "delay": 160, # OPX control EOM double pass
                #     # "delay": 400, # AWG control EOM
                #     "buffer": 0,
                # },
            },
            "digitalOutputs": {
                "OutBright1": (controller, 1),
                "OutBright2": (controller, 2),
            },
            'outputs': {
                'out1': (controller, 1)
            },
            'operations': {
                'Const_open': "MOT_lock",
                'MZ_balancing_pulses': "MZ_balance_pulses_Early",
                'QRAM_experiment_pulses_Early': "Square_pulse_seq_MZ_Early",
            },
            'time_of_flight': 36,
            'smearing': 0,
            'intermediate_frequency': IF_AOMs_MZ,
        },

        "AOM_Late": {
            "singleInput": {
                "port": (controller, 7),
            },
            'digitalInputs': {
                "AWG_Switch": {
                    "port": (controller, 10),
                    # "delay": 176, # OPX control EOM
                    "delay": 160,  # OPX control EOM double pass
                    # "delay": 400, # AWG control EOM
                    "buffer": 0,
                },
            },
            "digitalOutputs": {
                "OutDark1": (controller, 3),
                "OutDark2": (controller, 4),
            },
            'outputs': {
                'out1': (controller, 1)
            },
            'operations': {
                'Const_open': "MOT_lock",
                'MZ_balancing_pulses': "MZ_balance_pulses_Late",
                'QRAM_experiment_pulses_Late': "Square_pulse_seq_MZ_Late",
            },
            'time_of_flight': 36,
            'smearing': 0,
            'intermediate_frequency': IF_AOMs_MZ,
        },

        "PULSER_ANCILLA": {
            "singleInput": {
                "port": (controller, 8),
            },
            'digitalInputs': {
                # "Shutter_Switch": {
                #     "port": (controller, 5),
                #     "delay": 0,
                #     "buffer": 0,
                # },
                "SouthtoNorth_Shutter": {
                    "port": (controller, 9),
                    "delay": 0,
                    "buffer": 0,
                },
            },
            'operations': {
                'Const_open': "MOT_lock",
                'Const_open_triggered': "MOT_lock_ON",
                'Detection_pulses': "Square_detection_pulses",
                'Homodyne_Pulse': "Homodyne_Pulse",
                # 'Balancing_support': "Balancing_support_Ancilla",
                'QRAM_experiment_pulses_Ancilla': "QRAM_seq_pulse_Ancilla",
                'Spectrum_pulse': "Frequency_Sweep"
            },
            'intermediate_frequency': IF_AOM_ANCILLA,
        },

        "AOM_Spectrum": {
            'singleInput': {
                "port": (controller, 8)
            },
            'operations': {
                'Spectrum_pulse': "Frequency_Sweep2",
            },
            'intermediate_frequency': IF_AOM_Spectrum,
        },

        "PULSER_N": {
            "singleInput": {
                "port": (controller, 9
                         ),
            },
            'digitalInputs': {
                "Shutter_Switch": {
                    "port": (controller, 5),
                    "delay": 0,
                    "buffer": 0,
                },
                # "FS_South": {
                #     "port": (controller, 8),
                #     "delay": 500,
                #     "buffer": 0,
                # },
            },
            'operations': {
                'Const_open': "MOT_lock",
                'Const_open_triggered': "MOT_lock_ON",
                'Detection_pulses': "Square_detection_pulses",
                'Homodyne_Pulse': "Homodyne_Pulse",
                'MZ_balancing_pulses': "MZ_balance_pulses_N",
                'MZ_balancing_pulses_N': "MZ_balance_pulses_N",
                'MZ_balancing_pulses_S': "MZ_balance_pulses_S",
                'QRAM_experiment_pulses_N': "QRAM_seq_pulse_N",
            },
            'intermediate_frequency': IF_AOM_N,
        },

        "PULSER_S": {
            "singleInput": {
                "port": (controller, 10),
            },
            'digitalInputs': {
                # "SouthtoNorth_Shutter": {
                #     "port": (controller, 9),
                #     "delay": 0,
                #     "buffer": 0,
                # },
                "FS_South": {
                    "port": (controller, 8),
                    "delay": 500,
                    "buffer": 0,
                },
            },
            'operations': {
                'Const_open': "MOT_lock",
                'Const_open_triggered': "MOT_lock_ON",
                'Detection_pulses': "Square_detection_pulses",
                'Homodyne_Pulse': "Homodyne_Pulse",
                'MZ_balancing_pulses_N': "MZ_balance_pulses_S",
                'MZ_balancing_pulses_S': "MZ_balance_pulses_N",
                # 'MZ_balancing_pulses_S_check': "MZ_balance_pulses_N_check",
                'QRAM_experiment_pulses_S': "QRAM_seq_pulse_S",
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
        },

        "MOT_with_EOM_pulses": {
            'operation': 'control',
            'length': MOT_pulse_len,
            'waveforms': {
                'single': 'const_wf'
            },
            'digital_marker': 'Trig_EOM_MOT',
            # 'digital_marker': 'ON',
        },

        "MOT_with_Trigger": {
            'operation': 'control',
            'length': MOT_pulse_len,
            'waveforms': {
                'single': 'const_wf'
            },
            'digital_marker': 'Trig_AWG_MOT',
            # 'digital_marker': 'ON',
        },

        "MOT_lock_ON": {
            'operation': 'control',
            'length': MOT_pulse_len,
            'waveforms': {
                'single': 'const_wf'
            },
            'digital_marker': 'ON',
            # 'digital_marker': 'ON_pulse_switch'
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

        # "Gaussian_Sprint_pulse_S": {
        #     'operation': 'control',
        #     'length': len(Sprint_Exp_Gaussian_samples_S),
        #     'waveforms': {
        #         'single': 'Sprint_Gaussian_wf_S'
        #     },
        #     'digital_marker': 'ON'
        # },

        # "Gaussian_Sprint_pulse_N": {
        #     'operation': 'control',
        #     'length': len(Sprint_Exp_Gaussian_samples_N),
        #     'waveforms': {
        #         'single': 'Sprint_Gaussian_wf_N'
        #     }
        # },

        "Square_pulse_seq_MZ_Early": {
            'operation': 'control',
            'length': len(QRAM_Exp_Square_samples_Early),
            'waveforms': {
                'single': 'QRAM_Square_wf_Early'
            },
            # 'digital_marker': 'Trig_EOM'
            'digital_marker': 'ON'
        },

        "Square_pulse_seq_MZ_Late": {
            'operation': 'control',
            'length': len(QRAM_Exp_Square_samples_Late),
            'waveforms': {
                'single': 'QRAM_Square_wf_Late'
            },
            'digital_marker': 'Trig_EOM'
        },

        "QRAM_seq_pulse_Ancilla": {
            'operation': 'control',
            'length': len(QRAM_Exp_Gaussian_samples_Ancilla),
            'waveforms': {
                'single': 'QRAM_Gaussian_wf_Ancilla'
            },
            # 'digital_marker': 'QRAM_TOP2_pulsar'
            'digital_marker': 'ON'
        },

        "QRAM_seq_pulse_N": {
            'operation': 'control',
            'length': len(QRAM_Exp_Gaussian_samples_N),
            'waveforms': {
                'single': 'QRAM_Gaussian_wf_N'
            },
            'digital_marker': 'ON'
            # 'digital_marker': 'FS_North'
        },

        "QRAM_seq_pulse_S": {
            'operation': 'control',
            'length': len(QRAM_Exp_Gaussian_samples_S),
            'waveforms': {
                'single': 'QRAM_Gaussian_wf_S'
            },
            'digital_marker': 'FS_North'
            # 'digital_marker': 'ON'
        },

        "Frequency_Sweep": {
            'operation': 'control',
            'length': len(Gaussian_pulse_samples),
            'waveforms': {
                # 'single': 'const_wf'
                'single': 'wf_gaus'
            },
        },

        "Frequency_Sweep2": {
            'operation': 'control',
            'length': len(Gaussian_pulse_samples),
            'waveforms': {
                # 'single': 'const_wf'
                'single': 'wf_gaus2'
            },
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
            'length': int(max(readout_pulse_sprint_len_N, readout_pulse_sprint_len_S)),
            'operation': 'control',
            'waveforms': {
                'single': 'zero_wf'
            },
            # 'digital_marker': 'ON'
        },

        "MZ_balance_pulses_N": {
            'length': len(QRAM_MZ_balance_pulse_North),
            'operation': 'control',
            'waveforms': {
                'single': 'north_south_wf'
            },
            # 'digital_marker': 'Trig_EOM_MZ'
            'digital_marker': 'ON'
        },

        "MZ_balance_pulses_N_check": {
            'length': len(QRAM_MZ_balance_pulse_North),
            'operation': 'control',
            'waveforms': {
                'single': 'north_south_wf'
            },
            # 'digital_marker': 'Trig_EOM_MZ'
            # 'digital_marker': 'ON'
        },

        "MZ_balance_pulses_S": {
            'length': len(QRAM_MZ_balance_pulse_South),
            'operation': 'control',
            'waveforms': {
                'single': 'south_north_wf'
            },
            # 'digital_marker': 'Trig_EOM_MZ'
            'digital_marker': 'ON'
        },

        "MZ_balance_pulses_Early": {
            'length': len(QRAM_MZ_balance_pulse_Early),
            'operation': 'control',
            'waveforms': {
                'single': 'early_late_wf'
            },
            # 'digital_marker': 'Trig_EOM_MZ'
            'digital_marker': 'ON'
        },

        "MZ_balance_pulses_Late": {
            'length': len(QRAM_MZ_balance_pulse_Late),
            'operation': 'control',
            'waveforms': {
                'single': 'late_early_wf'
            },
            'digital_marker': 'Trig_EOM_MZ'
            # 'digital_marker': 'ON'
        },

        # TODO: Remove? Belongs to MW Spectroscopy
        # 'pulse1_in': {
        #     'operation': 'control',
        #     'length': pulse_time,
        #     'waveforms': {
        #         'I': 'wf1',
        #         'Q': 'wf2'
        #     },
        # },
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
        # 'Sprint_Gaussian_wf_S': {
        #     'type': 'arbitrary',
        #     'samples': Sprint_Exp_Gaussian_samples_S
        # },
        # 'Sprint_Gaussian_wf_N': {
        #     'type': 'arbitrary',
        #     'samples': Sprint_Exp_Gaussian_samples_N
        # },
        'QRAM_Gaussian_wf_Ancilla': {
            'type': 'arbitrary',
            'samples': QRAM_Exp_Gaussian_samples_Ancilla
        },
        'QRAM_Gaussian_wf_S': {
            'type': 'arbitrary',
            'samples': QRAM_Exp_Gaussian_samples_S
        },
        'QRAM_Gaussian_wf_N': {
            'type': 'arbitrary',
            'samples': QRAM_Exp_Gaussian_samples_N
        },
        'QRAM_Square_wf_Early': {
            'type': 'arbitrary',
            'samples': QRAM_Exp_Square_samples_Early_delayed
        },
        'QRAM_Square_wf_Late': {
            'type': 'arbitrary',
            'samples': QRAM_Exp_Square_samples_Late_delayed
        },
        'qram_wf': {
            'type': 'arbitrary',
            'samples': QRAM_Exp_TOP2_samples
        },
        'north_south_wf': {
            'type': 'arbitrary',
            'samples': QRAM_MZ_balance_pulse_North
        },
        'south_north_wf': {
            'type': 'arbitrary',
            'samples': QRAM_MZ_balance_pulse_South
        },
        'early_late_wf': {
            'type': 'arbitrary',
            'samples': QRAM_MZ_balance_pulse_Early_delayed
        },
        'late_early_wf': {
            'type': 'arbitrary',
            'samples': QRAM_MZ_balance_pulse_Late_delayed
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
        # TODO: Remove? Belongs to MW Spectroscopy
        # 'wf1': {
        #     'type': 'constant',
        #     'sample': wf_samples[0]
        # },
        # TODO: Remove? Belongs to MW Spectroscopy
        # 'wf2': {
        #     'type': 'constant',
        #     'sample': wf_samples[1]
        # },
        'wf_gaus': {
            'type': 'arbitrary',
            'samples': Gaussian_pulse_samples
        },
        'wf_gaus2': {
            'type': 'arbitrary',
            'samples': Gaussian_pulse_samples2
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
            "samples": [(1, 250), (0, 0)]
        },
        "Trig_EOM": {
            "samples": [(0, 250), (1, 250), (0, 250), (1, 0)]
            # "samples": [(0, 270), (1, 230), (0, 270), (1, 0)]
        },
        "Trig_EOM_MZ": {
            "samples": [(0, 250), (1, 250)] * MZ_balancing_seq_rep + [(0, 0)]
            # "samples": [(0, 270), (1, 230)] * MZ_balancing_seq_rep + [(0, 0)]
        },
        "Trig_EOM_MOT": {
            "samples": [(0, 250), (1, 250)] * int(MOT_pulse_len / (MZ_delay * 2)) + [(0, 0)]
            # "samples": [(0, 270), (1, 230)] * int(MOT_pulse_len / (MZ_delay * 2)) + [(0, 0)]
        },
        "Trig_AWG_MOT": {
            "samples": [(1, 20000), (0, 0)]
        },
        "FS_North": {
            "samples": QRAM_Exp_digital_samples_FS
        }
    },

    # TODO: Remove? Belongs to MW Spectroscopy
    # "mixers": {
    #     "my_mixer": [
    #         {"intermediate_frequency": im_freq, "lo_frequency": lo_freq, "correction": correction_matrix}
    #     ]
    # }
}

streams = {
    # "Detector_1_Counts": {
    #     "number": 1,
    #     "name": "Detector_1_Counts"
    # },
    # "Detector_2_Counts": {
    #     "number": 2,
    #     "name": "Detector_2_Counts"
    # },
    # "Detector_3_Counts": {
    #     "number": 3,
    #     "name": "Detector_3_Counts"
    # },
    # "Detector_4_Counts": {
    #     "number": 4,
    #     "name": "Detector_4_Counts"
    # },
    # "Detector_5_Counts": {
    #     "number": 5,
    #     "name": "Detector_5_Counts"
    # },
    # "Detector_6_Counts": {
    #     "number": 6,
    #     "name": "Detector_6_Counts"
    # },
    # "Detector_7_Counts": {
    #     "number": 7,
    #     "name": "Detector_7_Counts"
    # },
    # "Detector_8_Counts": {
    #     "number": 8,
    #     "name": "Detector_8_Counts"
    # },
    "Detector_1_Timetags": {
        "number": 1,
        "name": "Detector_1_Timetags",
        "type": "int",
        "binary": "I",  # Unsigned int
        "playback": "Det1_timetags.npz"
    },
    "Detector_2_Timetags": {
        "number": 2,
        "name": "Detector_2_Timetags",
        "type": "int",
        "binary": "I",  # Unsigned int
        "playback": "Det2_timetags.npz"
    },
    "Detector_3_Timetags": {
        "number": 3,
        "name": "Detector_3_Timetags",
        "type": "int",
        "binary": "I",  # Unsigned int
        "playback": "Det3_timetags.npz"
    },
    "Detector_4_Timetags": {
        "number": 4,
        "name": "Detector_4_Timetags",
        "type": "int",
        "binary": "I",  # Unsigned int
        "playback": "Det4_timetags.npz"
    },
    "Detector_5_Timetags": {
        "number": 5,
        "name": "Detector_5_Timetags",
        "type": "int",
        "binary": "I",  # Unsigned int
        "playback": "Det5_timetags.npz"
    },
    "Detector_6_Timetags": {
        "number": 6,
        "name": "Detector_6_Timetags",
        "type": "int",
        "binary": "I",  # Unsigned int
        "playback": "Det6_timetags.npz"
    },
    "Detector_7_Timetags": {
        "number": 7,
        "name": "Detector_7_Timetags",
        "type": "int",
        "binary": "I",  # Unsigned int
        "playback": "Det7_timetags.npz"
    },
    "Detector_8_Timetags": {
        "number": 8,
        "name": "Detector_8_Timetags",
        "type": "int",
        "binary": "I",  # Unsigned int
        "playback": "Det8_timetags.npz"
    },
    "FLR_measure": {
        "name": "FLR_measure",
        "type": "float",
        "binary": "d",  # double
        "playback": "Flouresence.npz"
    }
}
