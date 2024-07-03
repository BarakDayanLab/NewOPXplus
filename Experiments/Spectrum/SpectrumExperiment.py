from Experiments.BaseExperiment.BaseExperiment import BaseExperiment
from Experiments.BaseExperiment.BaseExperiment import TerminationReason
from Experiments.Spectrum import Spectrum_Config_Experiment as Config

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import time
import pymsgbox
from Utilities.BDSound import SOUNDS
from Utilities.Utils import Utils


class SpectrumExperiment(BaseExperiment):

    def __init__(self, playback=False, save_raw_data=False):
        # Invoking BaseClass constructor. It will initiate OPX, QuadRF, BDLogger, Camera, BDResults, KeyEvents etc.
        super().__init__(playback, save_raw_data)
        pass

    def __del__(self):
        print('**** SpectrumExperiment Destructor ****')
        super(type(self), self).__del__()
        pass

    # Override base-class method with the variables/values this experiment wants to use
    def initialize_experiment_variables(self):

        # -----------------------------------------------------------
        # Handle Free-Fall Variables
        # -----------------------------------------------------------

        self.FreeFall_duration = int(self.Exp_Values['FreeFall_duration'] * 1e6 / 4)
        self.Coils_timing = int(self.Exp_Values['Coil_timing'] * 1e6 / 4)
        if (self.Exp_Values['FreeFall_duration'] - self.Exp_Values['Coil_timing']) > 60:
            raise ValueError("FreeFall_duration - Coils_timing can't be larger than 60 ms")

        # -----------------------------------------------------------
        # Handle PGC Variables
        # -----------------------------------------------------------

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
                self.logger.error(
                    'The values for PGC_initial_amp_0 and PGC_final_amp_0 are too close or PGC_prep_duration is too long, might cause an ERROR!!!')
            if self.pgc_pulse_duration_minus > int(60 * 1e6 / 4):  # longer then 60ms
                self.logger.error(
                    'The values for PGC_initial_amp_minus and PGC_final_amp_minus are too close or PGC_prep_duration is too long, might cause an ERROR!!!')
            if self.pgc_pulse_duration_plus > int(60 * 1e6 / 4):  # longer then 60ms
                self.logger.error(
                    'The values for PGC_initial_amp_plus and PGC_final_amp_plus are too close or PGC_prep_duration is too long, might cause an ERROR!!!')
        self.pgc_initial_amp_0 = self.Exp_Values['PGC_initial_amp_0']
        self.pgc_initial_amp_minus = self.Exp_Values['PGC_initial_amp_minus']
        self.pgc_initial_amp_plus = self.Exp_Values['PGC_initial_amp_plus']
        self.pgc_final_amp_0 = self.Exp_Values['PGC_final_amp_0']
        self.pgc_final_amp_minus = self.Exp_Values['PGC_final_amp_minus']
        self.pgc_final_amp_plus = self.Exp_Values['PGC_final_amp_plus']
        # self.pgc_aom_chirp_rate = int(self.Exp_Values['PGC_final_Delta_freq'] * 1e3 / (self.Exp_Values['PGC_prep_duration'] * 1e6))  # [mHz/nsec], If needed pgc preparation duration must be constant!!!

        # -----------------------------------------------------------
        # Handle Fountain Variables
        # -----------------------------------------------------------

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
                self.logger.error(
                    'The values for Fountain_initial_amp_0 and Fountain_final_amp_0 are too close or Fountain_prep_duration is too long, might cause an ERROR!!!')
            if self.fountain_pulse_duration_minus > int(60 * 1e6 / 4):  # longer then 60ms
                self.logger.error(
                    'The values for Fountain_initial_amp_minus and Fountain_final_amp_minus are too close or Fountain_prep_duration is too long, might cause an ERROR!!!')
            if self.fountain_pulse_duration_plus > int(60 * 1e6 / 4):  # longer then 60ms
                self.logger.error(
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

        # -----------------------------------------------------------
        # MW spectroscopy parameters:
        # -----------------------------------------------------------

        self.MW_start_frequency = int(100e6)  # [Hz]
        self.Pulse_Length_MW = 400  # [usec]
        self.Pulse_Length_OD = 20  # [usec]

        # -----------------------------------------------------------
        # Main Experiment:
        # -----------------------------------------------------------

        self.TOP2_pulse_len = int(Config.Probe_pulse_len / 4)  # [nsec]
        self.Calibration_time = 10  # [msec]

        # Spectrum Experiment:
        self.frequency_sweep_rep = 5
        self.spectrum_bandwidth = int(48e6)  # experiment spe13/11/23
        self.frequency_start = int(80e6 - self.spectrum_bandwidth/2)
        # self.frequency_diff = int(self.spectrum_bandwidth / self.num_of_different_frequencies)
        self.frequency_diff = int(1.5e6)  # difference between frequencies
        self.num_of_different_frequencies = int(self.spectrum_bandwidth // self.frequency_diff) + 1
        self.same_frequency_rep = int(self.Exp_Values['Pulse_1_duration'] * 1e6 /
                                      (2 * Config.frequency_sweep_duration) /
                                      (self.num_of_different_frequencies * self.frequency_sweep_rep)) # total time[10ms]/((1us sequence - 2 pulses time)*(num of sweeps)*(num of different freqs) )
        pass

    # ## Boolean functions: ##
    def update_pgc_final_amplitude(self, final_amplitude):
        x = np.log(int(1 / (1 - final_amplitude)))
        self.logger.info('Called update_pgc_final_amplitude', final_amplitude)
        # self.update_lin_pgc_final_amplitude(final_amplitude)
        # self.update_exp_pgc_final_amplitude(final_amplitude)
        return x


    # TODO: Not in use for this experiment
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
        self.logger.info(Probe_counts_South)

        for n in range(repetitions - 1):
            while np.sum(Probe_S_res.tolist()) == Probe_counts_South[-1]:
                Probe_N_res = Probe_N_handle.fetch_all()
                Probe_S_res = Probe_S_handle.fetch_all()
            Probe_counts_North.append(np.sum(Probe_N_res.tolist()))
            Probe_counts_South.append(np.sum(Probe_S_res.tolist()))

        self.logger.info(r'$Max Probe_N = %.1f$' % ((np.average(Probe_counts_North) * 1000) / self.M_time,) + '[MPhotons/sec]')
        self.logger.info(r'$Max Probe_S = %.1f$' % ((np.average(Probe_counts_South) * 1000) / self.M_time,) + '[MPhotons/sec]')

        self.Max_Probe_counts_switch(False)
        self.update_parameters()

        return [(np.average(Probe_counts_North) * 1000) / self.M_time,
                (np.average(Probe_counts_South) * 1000) / self.M_time]

    def get_handles_from_OPX_Server(self):
        """
         Gets handles of time-tags and counts from OPX
        :return:
        """
        Counts_handle = []
        tt_handle = []

        for i in self.Num_Of_dets:
            Counts_handle.append(self.job.result_handles.get("Det" + str(i) + "_Counts"))
            tt_handle.append(self.job.result_handles.get("Det" + str(i) + "_Probe_TT"))
        FLR_handle = self.job.result_handles.get("FLR_measure")
        return Counts_handle, tt_handle, FLR_handle

    def ingest_time_tags(self):
        """
        Takes all the raw results we got from the streams and does some processing on them - preparing "measures"
        """
        # Get fluorescence data
        self.FLR_res = -float(self.streams['FLR_measure']['results'])
        self.fluorescence_average = 1000 * self.FLR_res

        # Ensure we clear time-tags vectors from data from last iteration
        self.tt_measure = []

        # Normalize the data coming from the detectors
        for detector_index in range(len(self.Num_Of_dets)):  # for different detectors
            tt_res = self.streams[f'Detector_{detector_index + 1}_Timetags']['results']
            normalized_stream = self.bdstreams.normalize_stream(tt_res, detector_index, Config.detector_delays, self.M_window)
            self.tt_measure.append(normalized_stream)

        # ------------------------------------------
        # Unify detectors and windows within detectors and create vector of tt's for each direction (bright port, dark port, north, south and from FS) and sort them
        # ------------------------------------------

        self.tt_BP_measure = []
        self.tt_DP_measure = []
        self.tt_N_measure = []
        self.tt_S_measure = []
        self.tt_FS_measure = []

        # Unify detectors 0 & 1 and windows within detectors
        self.tt_BP_measure = sorted(sum(self.tt_measure[0:2], []))

        # Unify detectors 2 & 3 and windows within detectors
        self.tt_DP_measure = sorted(sum(self.tt_measure[2:4], []))

        # Unify detectors 7 & [] and windows within detectors
        self.tt_N_measure = sorted(self.tt_measure[7])

        # Unify detectors 4 & [] and windows within detectors
        self.tt_S_measure = sorted(self.tt_measure[4])

        # Unify detectors 5 & 6 and windows within detectors
        self.tt_FS_measure = sorted(sum(self.tt_measure[5:7], []))

        self.tt_N_directional_measure = sorted(self.tt_N_measure + self.tt_BP_measure + self.tt_DP_measure)
        # self.tt_N_directional_measure = sorted(self.tt_S_measure + self.tt_FS_measure)

        # self.tt_S_directional_measure = sorted(self.tt_N_measure + self.tt_BP_measure + self.tt_DP_measure)
        self.tt_S_directional_measure = sorted(self.tt_S_measure + self.tt_FS_measure)

        # Fix gaps created by OPX (take tt_S_directional_measure and convert it to self.tt_S_no_gaps)
        self.fix_gaps_spectrum_exp_tts()

    # TODO: we need to bring this NEW method up-to-speed with the OLD one (Below). Ensure the .copy part.
    def fix_gaps_spectrum_exp_tts_NEW(self):
        self.tt_S_no_gaps = self.tt_S_directional_measure.copy()
        for i, x in enumerate(self.tt_S_no_gaps):
            y = x-(x//(1e5 + 320))*320
            z = y-(y//(2e6 + 124))*124
            self.tt_S_no_gaps[i] = int(z)
        pass

    def fix_gaps_spectrum_exp_tts(self):
        """
        Fixes delays of 320 ns every 100 us in spectrum experiment
        These gaps are caused by OPX when it is switching frequencies

        Input:   self.tt_S_directional_measure
        Output:  self.tt_S_no_gaps
        """
        self.real_M_window = ((2 * Config.frequency_sweep_duration) * (self.num_of_different_frequencies
                                                                       * self.same_frequency_rep
                                                                       * self.frequency_sweep_rep))
        self.total_real_time_of_freq_sweep = int(self.real_M_window // self.frequency_sweep_rep) + \
                                             self.num_of_different_frequencies * 320 + 124
        self.total_real_time_of_same_freq_rep = (Config.frequency_sweep_duration * 2 * self.same_frequency_rep) + 320
        self.tt_S_no_gaps = [x - (x // self.total_real_time_of_freq_sweep) * 124 for x in self.tt_S_directional_measure]
        self.tt_S_no_gaps = [x - (x // self.total_real_time_of_same_freq_rep) * 320 for x in self.tt_S_no_gaps]

    # TODO: In the case of Spectrum, we're checking only specific detectors. We need to make this configurable
    def latched_detectors(self):
        latched_detectors = []
        for indx, det_tt_vec in enumerate(self.tt_measure[4:7]):  # for different detectors
            if not det_tt_vec:
                latched_detectors.append(indx)
        return latched_detectors

    def get_pulses_bins(self, det_pulse_len, sprint_pulse_len, num_of_det_pulses, num_of_sprint_pulses,
                        sprint_sequence_delay, num_of_sprint_sequences, num_init_zeros, num_fin_zeros,
                        num_between_zeros):
        """
        generating bin vwctor with edges at the efges of pulses (such that the histogram would give the number of
        clicks within pulse).
        :param det_pulse_len:
        :param sprint_pulse_len:
        :param num_of_det_pulses:
        :param num_of_sprint_pulses:
        :param sprint_sequence_delay:
        :param num_of_sprint_sequences:
        :param num_init_zeros:
        :param num_fin_zeros:
        :param num_between_zeros:
        :return:
        """
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
        """
        dividing south and north tt's vectros into reflection and transmission vectors, by:
        1. converting them into counts vector in same time spacing as the initial pulse.
        2. taking the relevant pulses from the sequence to each vector (reflection/transmission)
        :return:
        """
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

        ### build reflection and transmission pulses indices, by setting non-zero at indices elements ###
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


    def plot_folded_tt_histogram(self):

        plt.figure()
        plt.plot(sum(np.reshape(self.tt_histogram_N, [int(len(self.tt_histogram_N) / 9), 9])), label='tt_hist_N')
        plt.plot(sum(np.reshape(self.tt_histogram_S, [int(len(self.tt_histogram_S) / 9), 9])), label='tt_hist_S')
        plt.legend()

    def find_transits_events_spectrum_exp(self, tt_resonance_binning, N, cond=None, minimum_number_of_seq_detected=2):
        """
        Find transits of atoms by searching events that satisfy the number of reflected photons per sequence required at
        each cond place with minimum number of conditions needed to be satisfied, defined by the minimum_number_of_seq_detected.
        For example:
        given cond=[2,1,2] and minimum_number_of_seq_detected=2, if in 3 consecutive sequences we get either 2-1-0 or
        0-2-1 or 2-0-2, the condition is satisfied, and it is defined as a transit.
        :param cond: The condition that need to be met for number of reflections per detection pulses in sequence.
        :param minimum_number_of_seq_detected: The number of terms needed to be satisfied per cond vector.
        :return:
        """

        if cond is None:
            cond = [2, 1, 2]

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
            self.Transits_per_frequency += (self.tt_per_frequency != 0)  # TODO: removed this: .astype(int). Test it.

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
            self.batcher.add("all_transits_index_batch", self.all_transits_seq_indx)

    def get_pulses_location_in_seq(self, delay, seq=Config.QRAM_Exp_Gaussian_samples_S,
                                   smearing=int(Config.num_between_zeros / 2)):
        '''
        A function that uses the original sequence samples that the OPX uses, in order to obtain the location of the
        pulses in the sequence and build a filter. The user may add smearing which is the value that is added before and
        after each pulse in the sequence to match the filter to the performance of the physical system (AOMs).
        :param delay: Between the actual sequence pulse location to the location of the folded data
        :param seq: The sequence of pulses from which the filter is generated.
        :param smearing: The value that is added to the filter before and after each pulse in the sequence.
        :return:
        '''
        seq_filter = (np.array(seq) > 0).astype(int)
        seq_filter = np.roll(seq_filter, delay)
        seq_index = np.where(seq_filter > 0)[0]
        seq_filter_with_smearing = seq_filter
        pulses_loc = []
        if seq_index.any():
            start_indx = seq_index[0]
            for i in range(1, len(seq_index)):
                if seq_index[i] - seq_index[i - 1] > 1:
                    pulses_loc.append((start_indx - int(smearing), seq_index[i - 1] + int(smearing)))
                    for j in range(start_indx - int(smearing), seq_index[i - 1] + int(smearing)):
                        seq_filter_with_smearing[j] = 1
                    start_indx = seq_index[i]
            pulses_loc.append((start_indx - int(smearing), seq_index[-1] + int(smearing)))
            for j in range(start_indx - int(smearing), seq_index[-1] + int(smearing)):
                seq_filter_with_smearing[j] = 1
        return pulses_loc, seq_filter_with_smearing

    def load_data_for_playback(self):
        """
        This function loads relevant data saved in previous experiments
        and re-constructs the original streams raw data.
        (the raw data will serve for the playback runs)
        """
        playback_folder = self.bd_results.get_custom_root('playback_data')

        # Iterate over all stream and load data from the relevant file
        for stream in self.streams.values():
            self.info(f"Loading stream for {stream['name']}...")
            if stream['handler'] is not None:
                stream['handler'] = 'playback: data loaded!'  # Mark the stream as loaded
                if stream['name'].startswith('Detector_'):
                    # Load the file and add first element with the size
                    data_file = os.path.join(playback_folder, 'AllDetectors', stream['playback'])
                    zipped_data = np.load(data_file, allow_pickle=True)
                    stream['all_rows'] = zipped_data['arr_0'].tolist()

                    # Add the time-tags count as the first element of each array
                    if len(stream['all_rows']) == 0:
                        stream['all_rows'] = [0]
                    else:
                        # Reduce the detector delay that was added originally before experiment saved the data
                        delay = Config.detector_delays[stream['number'] - 1]
                        for i in range(len(stream['all_rows'])):
                            if len(stream['all_rows'][i]) > 0:
                                stream['all_rows'][i] = [time_tag - delay for time_tag in stream['all_rows'][i]]
                        # Add the count before the data
                        stream['all_rows'] = [[len(arr)] + arr for arr in stream['all_rows']]

        # Load flouresence
        data_file = os.path.join(playback_folder, 'Flouresence.npz')
        zipped_data = np.load(data_file, allow_pickle=True)
        self.streams['FLR_measure']['all_rows'] = (-zipped_data['arr_0']).tolist()

        pass

    def print_experiment_information(self):
        """
        Prints during the experiment a formatted line that brings all concise information
        """
        locking_efficiency = self.counter / self.repetitions
        locking_efficiency_str = '%.2f' % locking_efficiency
        fluorescence_str = '%.2f' % self.fluorescence_average
        time_formatted = time.strftime("%Y/%m/%d, %H:%M:%S")
        status_str = f'[Warm Up: {self.warm_up_cycles}]' if self.warm_up else f'# {self.counter} ({self.repetitions})'

        self.logger.info(f'{time_formatted}: {status_str}, Eff: {locking_efficiency_str}, Flr: {fluorescence_str}, #Photons/[us]: {self.sum_for_threshold}')

    def experiment_mainloop_delay(self):
        if self.playback['active']:
            if self.playback['delay'] != -1:
                time.sleep(self.playback['delay'])
        else:
            time.sleep(0.01)  # TODO: do we need this delay?

    def prepare_figures(self):

        super().prepare_figures()

        self.subplots = []

        self.subplots.append(plt.subplot2grid((6, 2), (0, 0), colspan=1, rowspan=2))
        self.subplots.append(plt.subplot2grid((6, 2), (0, 1), colspan=1, rowspan=2))
        self.subplots.append(plt.subplot2grid((6, 2), (2, 0), colspan=1, rowspan=2))
        self.subplots.append(plt.subplot2grid((6, 2), (2, 1), colspan=1, rowspan=2))
        self.subplots.append(plt.subplot2grid((6, 2), (4, 0), colspan=1, rowspan=2))
        self.subplots.append(plt.subplot2grid((6, 2), (4, 1), colspan=1, rowspan=2))

        #self.maximize_figure()  # TODO: uncomment

        return

    def plot_figures(self):
        try:
            self._plot_figures()
        except Exception as err:
            self.warn(f'Failed on plot_figures: {err}')
        pass

    def _plot_figures(self):

        if self.playback['active']:
            if self.playback['plot'] == 'NONE':
                return
            if self.playback['plot'] == 'LAST' and self.counter < (self.N - 10):
                return

        # If it is the first time, plot.
        if not self.plot_shown:
            plt.show(block=False)
            self.plot_shown = True

        plot_switches = {
            "playsound": True,
            "header": True,
            "detectors": True,
            "graph-0": True,
            "graph-1": True,
            "graph-2": True,
            "graph-3": True,
            "graph-4": True,
            "graph-5": True
        }

        # Take time
        start_plot_time = time.time()

        # Used mainly as a shortcut for shorter lines
        ax = self.subplots

        ax[0].clear()
        ax[1].clear()
        ax[2].clear()
        ax[3].clear()
        ax[4].clear()
        ax[5].clear()

        # Play sound:
        if plot_switches['playsound']:
            if not self.acquisition_flag:
                self.bdsound.play(SOUNDS.INDICATOR_1)

        # Prepare header text
        pause_str = ' , PAUSED!' if self.pause_flag else ''
        trs_str = '%.2f' % self.sum_for_threshold
        flr_str = '%.2f' % self.fluorescence_average
        eff_str = '%.2f' % (self.counter / self.repetitions)
        lck_str = '%.3f' % self.lock_err
        status_str = f'[Warm Up: {self.warm_up_cycles}]' if self.warm_up else f'# {self.counter} ({self.repetitions})'
        playback_str = 'PLAYBACK: ' if self.playback['active'] else ''
        header_text = f'{playback_str} {status_str} - Trans: {trs_str}, Eff: {eff_str}, Flr: {flr_str}, Lock Error: {lck_str} {pause_str}'

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        textstr_detuned = r'$S_{detuned} = %.3f$' % (
        (sum(self.tt_S_binning_detuned) * 1000) / (self.M_time / 2),) + '[MPhotons/sec]\n' \
                          + '$\overline{S}_{detuned} = %.3f$' % (
                          (np.mean([sum(x) for x in self.batcher['tt_S_binning_detuned_batch']]) * 1000)
                          / (self.M_time / 2),) + '[MPhotons/sec]'
        textstr_resonance = r'$S_{res} = %.3f$' % (
        (sum(self.tt_S_binning_resonance) * 1000) / (self.M_time / 2),) + '[MPhotons/sec]\n' \
                            + '$\overline{S}_{res} = %.3f$' % (
                            (np.mean([sum(x) for x in self.batcher['tt_S_binning_resonance_batch']]) * 1000)
                            / (self.M_time / 2),) + '[MPhotons/sec]'
        textstr_no_transits = 'NO TRANSITS YET!!!'

        if plot_switches['graph-0']:
            flag_color = 'green' if self.acquisition_flag else 'red'
            props_thresholds = dict(boxstyle='round', edgecolor=flag_color, linewidth=2, facecolor=flag_color, alpha=0.5)
            ax[0].plot(self.time_bins[::2], self.S_bins_res_acc, label='Counts histogram', color='b')
            ax[0].set_title('On resonance histogram accumulated', fontweight="bold")
            ax[0].set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
            ax[0].text(0.05, 0.95, textstr_resonance, transform=ax[0].transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
            # Add figures header
            ax[0].text(0.05, 1.4, header_text, transform=ax[0].transAxes, fontsize=26,
                       verticalalignment='top', bbox=props_thresholds)
            ax[0].legend(loc='upper right')
            #self.debug("Plot: Drew graph 0")

        if plot_switches['graph-1']:
            ax[1].plot(self.time_bins[::2], self.S_bins_detuned_acc, label='Counts histogram', color='b')
            ax[1].set_title('Detuned histogram accumulated', fontweight="bold")
            ax[1].set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
            ax[1].text(0.05, 0.95, textstr_detuned, transform=ax[1].transAxes, fontsize=12, verticalalignment='top', bbox=props)
            ax[1].legend(loc='upper right')
            #self.debug("Plot: Drew graph 1")

        if plot_switches['graph-2']:
            ax[2].plot(self.time_bins[::2], self.tt_S_binning_resonance, label='Counts histogram', color='b')
            ax[2].set_title('On resonant counts (live)', fontweight="bold")
            ax[2].set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')

            # Detectors status:
            if plot_switches['detectors']:
                latched_detectors = self.latched_detectors()
                for i, det in enumerate(self.Num_Of_dets):
                    x = -0.15
                    y = 1.9 - i * 0.4
                    pad = 1
                    text = 'Det %d' % det
                    det_color = 'red' if latched_detectors else 'green'
                    ax[2].text(x, y, text, ha="center", va="center", transform=ax[2].transAxes,
                             bbox=dict(boxstyle=f"circle,pad={pad}", edgecolor=det_color, linewidth=2, facecolor=det_color, alpha=0.5))
            ax[2].legend(loc='upper right')
            #self.debug("Plot: Drew graph 2")

        if plot_switches['graph-3']:
            ax[3].plot(self.freq_bins/1e6, self.Cavity_atom_spectrum_normalized, label='Cavity-atom Spectrum', color='k')
            ax[3].plot(self.freq_bins/1e6, self.Cavity_spectrum_normalized, label='Cavity Spectrum', color='b')
            ax[3].set(xlabel='Frequency [MHz]', ylabel='Transmission Normalized')
            ax[3].legend(loc='upper right')
            #self.debug("Plot: Drew graph 3")

        if plot_switches['graph-4']:
            ax[4].plot(np.linspace(0, self.histogram_bin_size-1, self.histogram_bin_size), self.folded_tt_S_acc_1, label='pulses folded', color='k')
            ax[4].plot(np.linspace(0, self.histogram_bin_size-1, self.histogram_bin_size), self.folded_tt_S_acc_2, label='pulses folded_2', color='b')
            ax[4].plot(np.linspace(0, self.histogram_bin_size-1, self.histogram_bin_size), self.folded_tt_S_acc_3, label='pulses folded_2', color='g')
            ax[4].set(xlabel='Time [nsec]', ylabel='Counts [#Number]')
            ax[4].legend(loc='upper right')
            #self.debug("Plot: Drew graph 4")

        if plot_switches['graph-5']:
            if len(self.batcher['all_transits_index_batch']) > 0:
                # if self.all_transits:
                textstr_transit_counts = r'$N_{Transits} = %s $' % (len(self.all_transits_seq_indx),) + r'$[Counts]$'
                textstr_transit_event_counter = r'$N_{Transits Total} = %s $' % (
                                                len([vec for elem in self.batcher['all_transits_index_batch'] for vec in elem]),) \
                                                + r'$[Counts]$' + '\n' + textstr_transit_counts

                ax[5].plot(self.time_bins[::2], self.tt_S_transit_events_accumulated, label='Transit events histogram', marker='*', color='b')
                ax[5].set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
                ax[5].text(0.05, 0.95, textstr_transit_event_counter, transform=ax[5].transAxes, fontsize=12,
                           verticalalignment='top', bbox=props)
            else:
                ax[5].plot(self.time_bins[::2], self.tt_S_transit_events_accumulated, label='Transit events histogram', marker='*', color='b')
                ax[5].set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
                ax[5].text(0.25, 0.5, textstr_no_transits, transform=ax[5].transAxes, fontsize=24,
                           verticalalignment='center', bbox=props)
            #self.debug("Plot: Drew graph 5")

        # End timer
        total_prep_time = time.time() - start_plot_time

        # plt.tight_layout()
        plt.pause(0.5)

    # TODO: Document this method
    # TODO: Q: Do we need to move this to BaseExperiment
    def AOM_power_per_freq_calibration(self, calibration_file):
        if calibration_file is None:
            return np.full((1, self.spectrum_bin_number), 1)[0]
        if not os.path.exists(calibration_file):
            return np.full((1, self.spectrum_bin_number), 1)[0]
        self.calibration_spectrum = np.load(calibration_file)
        return self.calibration_spectrum['arr_0']/max(self.calibration_spectrum['arr_0'])

    def new_timetags_detected(self, prev_measures, curr_measure):
        """
        Attempt to deduce if the time-tags that are in are new ones (as opposed to leftovers from prev measure)

        returns: True/False
        """

        # If we don't have any previous measures, clearly it's new data
        if len(prev_measures) == 0:
            return True

        # Take the last sample from batch of previous measurements
        prev_measure = prev_measures[-1]

        # If last measure was empty and this one is not, clearly it's new data
        if len(prev_measure) == 0 and len(curr_measure) > 0:
            return True

        # Get the minimum value between new measure (tt_S_no_gaps) and last old measure (tt_S_measure_batch)
        min_len = min(len(prev_measure), len(curr_measure))

        # Compare values of measurements
        compare_values = np.array(prev_measure[:min_len]) == np.array(curr_measure[:min_len])

        # TODO: replace the line below with the one in comment:
        #new_timetags_found = np.sum(compare_values) < min_len / 2
        new_timetags_found = sum(compare_values) < min_len / 2
        return new_timetags_found

    def should_terminate(self):
        """
        Decides if to terminate mainloop. Overrides BaseExperiment method.
        Did we complete N successful iterations? (we check for N+1 because we started with 1)

        Note we also consult with super() - if it tells us to terminate, we "listen" :-)
        """
        flag = super().should_terminate()
        return flag or self.counter == self.N+1

    def await_for_values(self):
        """
        Awaits for values to come from OPX streams. Returns the timestamp of the incoming data.
        Returns the status of the operation (enumeration). It can be either:
        - Successful - we got data
        - User terminated - user terminated the experiment
        - Too many cycles - no data arrived after a given number of cycles
        """
        WAIT_TIME = 0.01  # [sec]
        #TOO_MANY_WAITING_CYCLES = WAIT_TIME*100*20  # 5 seconds. Make this -1 to wait indefinitely
        TOO_MANY_WAITING_CYCLES = -1

        count = 1
        while True:

            self.handle_user_events()
            if self.runs_status == TerminationReason.USER:
                break

            # Get data from OPX streams
            self.get_results_from_streams()
            if self.runs_status != TerminationReason.SUCCESS:
                break

            self.ingest_time_tags()

            is_new_tts_S = True  # We'll assume it's new data. Especially true in the first run, where there's no previous data to compare
            curr_measure = self.tt_S_no_gaps

            # Check if new time tags arrived - compare previous time-tags with current:
            is_new_tts_S = self.new_timetags_detected(self.batcher['tt_S_measure_batch'], curr_measure)

            # Bin the South time-tags we just got
            self.tt_S_binning = Utils.bin_values(values=curr_measure, bin_size=Config.frequency_sweep_duration, num_of_bins=self.histogram_bin_number * 2)

            # If there are time-tags in the last 5%, it is ok!
            S_det_is_full = Utils.tail_contains_non_zeros(values=self.tt_S_binning, percent=5)

            # We break if it is new data and not full
            if is_new_tts_S and S_det_is_full:
                break

            # Are we waiting too long?
            count += 1
            if TOO_MANY_WAITING_CYCLES != -1 and count > TOO_MANY_WAITING_CYCLES:
                self.runs_status = TerminationReason.ERROR
                break

            self.info(f'Waiting for values from stream (count: {count}). No new data coming from detectors (exp_flag = {self.exp_flag})')
            time.sleep(WAIT_TIME)

        # We return the timestamp - this is when we declare we got the measurement
        timestamp = time.strftime("%Y%m%d-%H%M%S")
        return timestamp

    def experiment_calculations(self):

        self.tt_N_binning_avg = self.tt_N_binning
        self.tt_S_binning_avg = self.tt_S_binning

        # Prepare folded arrays for display
        for x in self.tt_S_no_gaps:
            # if x < (2e6 + 6400):
            self.folded_tt_S_acc_1[(x - 1) % self.histogram_bin_size] += 1
            # if (2e6 + 6400 + 120) < x < (4e6 + 6400*2 + 120):
            if 6e6 < x < 8e6:
                self.folded_tt_S_acc_2[(x - 1) % self.histogram_bin_size] += 1
            # if (8e6 + 6400*4 + 120*4) < x < (10e6):
            if 8e6 < x < 10e6:
                self.folded_tt_S_acc_3[(x - 1) % self.histogram_bin_size] += 1

        # Split the binning vector to even and odd - Off/On resonance pulses
        self.tt_N_binning_resonance, self.tt_N_binning_detuned = Utils.split_to_even_odd_elements(self.tt_N_binning)
        self.tt_S_binning_resonance, self.tt_S_binning_detuned = Utils.split_to_even_odd_elements(self.tt_S_binning)

        # Calculate threshold for experiment
        self.sum_for_threshold = (np.sum(self.tt_S_binning_resonance) * 1000) / (self.M_time / 2)

        pass

    def run(self, sequence_definitions, run_parameters):
        """
        Function for analyzing,saving and displaying data from spectrum experiment.
        :param N: Number of maximum experiments (free throws) saved and displayed.
                 program we are looking for transits of atoms next to the toroid and record them.
        :param transit_profile_bin_size:
        :param pre_comment: The comment added at the start of the experiment, usually consisting of unique experiment
                           parameters.
        :param transit_cond:
        :param total_counts_threshold:
        :param transit_counts_threshold:
        :param FLR_threshold:
        :param lock_err_threshold: The maximal error in locking resonator to rb line (in ms, depends on the scan width/vpp)
        :param exp_flag:
        :param with_atoms:

        :return: Returns the run status (Success, User-Termination, Error, etc.)
        """

        self.run_parameters = run_parameters
        self.N = run_parameters['N']
        self.with_atoms = run_parameters['with_atoms']
        self.pre_comment = run_parameters['pre_comment']
        self.lock_err_threshold = run_parameters['lock_err_threshold']
        self.FLR_threshold = run_parameters['FLR_threshold']
        self.transit_cond = run_parameters['transit_cond']
        self.transit_profile_bin_size = run_parameters['transit_profile_bin_size']
        self.exp_flag = run_parameters['exp_flag'] if 'exp_flag' in run_parameters else False

        # Handle pre-comment - keep what we got in parameters or prompt the user for a comment
        if not self.pre_comment:
            self.pre_comment = self.prompt(title='Pre Experiment Run', msg='Add pre-run comment to measurement (click Cancel for "Ignore")')
            if self.pre_comment == None:
                self.pre_comment = 'Ignore'

        self.switch_atom_no_atoms = "atoms" if self.with_atoms else "!atoms"

        # set constant parameters for the function
        self.Num_Of_dets = [1, 2, 3, 4, 5, 6, 7, 8]
        self.histogram_bin_size = Config.frequency_sweep_duration * 2
        self.sequence_duration = Config.frequency_sweep_duration * 2
        self.histogram_bin_number = self.M_time // self.histogram_bin_size # num of bins in cycle of frequency sweep
        self.time_bins = np.linspace(start=0, stop=self.M_time, num=self.histogram_bin_number*2)
        self.spectrum_bin_number = self.num_of_different_frequencies
        self.freq_bins = np.linspace((self.frequency_start - Config.IF_AOM_Spectrum) * 2,
                                     (self.frequency_start + self.spectrum_bandwidth - Config.IF_AOM_Spectrum) * 2,
                                     self.spectrum_bin_number)
        time_threshold = int(self.histogram_bin_size * 1.6)  # The minimum time between two time tags to be counted for a transit. # TODO: might need a factor of 2???

        self.logger.blue('Press ESC to stop measurement.')

        # Initialize the batcher
        self.batcher.set_batch_size(self.N)
        self.batcher.empty_all()

        # Associate the streams filled in OPX (FPGA) code with result handles
        self.get_handles_from_OPX_server()

        # TODO: move this into the "real"/"united" initialization method - pre_run ?
        # Initialization
        self.tt_N_binning = np.zeros(self.histogram_bin_number * 2)
        self.tt_N_transit_events = np.zeros(self.histogram_bin_number)
        self.tt_S_transit_events = np.zeros(self.histogram_bin_number)
        self.tt_S_transit_events_accumulated = np.zeros(self.histogram_bin_number)
        self.all_transits_accumulated = np.zeros(self.histogram_bin_number)
        self.folded_tt_S_acc_1 = np.zeros(self.histogram_bin_size, dtype=int)
        self.folded_tt_S_acc_2 = np.zeros(self.histogram_bin_size, dtype=int)
        self.folded_tt_S_acc_3 = np.zeros(self.histogram_bin_size, dtype=int)

        self.Cavity_atom_spectrum = np.zeros(self.spectrum_bin_number)
        self.Transits_per_frequency = np.zeros(self.spectrum_bin_number)
        self.Cavity_atom_spectrum_normalized = np.zeros(self.spectrum_bin_number)
        self.Cavity_spectrum = np.zeros(self.spectrum_bin_number)
        self.Cavity_spectrum_normalized = np.zeros(self.spectrum_bin_number)

        # Get calibration weights from file
        calibration_folder = self.bd_results.get_custom_root('calibration_folder')
        calibration_file = os.path.join(calibration_folder, 'Cavity_spectrum.npz')
        self.power_per_freq_weight = self.AOM_power_per_freq_calibration(calibration_file)

        # Prepare the sub figures needed for the plotting phase
        self.prepare_figures()

        ############################################ START MAIN LOOP #################################################

        # Initialize the variables before loop starts
        self.counter = 1  # Total number of successful cycles
        self.repetitions = 1  # Total number of cycles
        self.acquisition_flag = True
        self.warm_up = True
        self.post_warm_up_completed = False  # Did we finish the post warm-up process

        self.pause_flag = False
        self.runs_status = TerminationReason.SUCCESS  # default status of the run is Success

        # ---------------------------------------------------------------------------------
        # Experiment Mainloop
        # - We iterate few warm-up rounds
        # - We iterate until we got N successful runs, or user terminated
        # ---------------------------------------------------------------------------------
        while True:

            # Handle cases where user terminates or pauses experiment
            self.handle_user_events()
            if self.runs_status == TerminationReason.USER:
                break

            # Get locking error
            self.lock_err = self._read_locking_error()

            # Experiment delay
            self.experiment_mainloop_delay()

            # Handle warm-up phase
            self.warm_up = self.handle_warm_up_phase()
            if self.warm_up:
                continue

            # Await for values from OPX
            self.sampling_timestamp = self.await_for_values()
            if self.runs_status != TerminationReason.SUCCESS:
                break

            # Bin North time-tags (South time-tags were binned in await_for_values)
            self.tt_N_binning = Utils.bin_values(values=self.tt_N_directional_measure, bin_size=Config.frequency_sweep_duration, num_of_bins=self.histogram_bin_number*2)

            # Informational printing
            self.print_experiment_information()

            # Perform all analytics and calculations needed for display
            self.experiment_calculations()

            # Plot figures
            self.plot_figures()

            # Determine if we're "acquired" - e.g., worthy to save data :-)
            self.acquisition_flag = self.is_acquired()
            if self.acquisition_flag and not self.pause_flag:

                if self.counter < self.N:
                    self.counter += 1

                self.tt_N_binning = np.zeros(self.histogram_bin_number*2)
                self.tt_S_binning = np.zeros(self.histogram_bin_number*2)
                self.tt_S_transit_events = np.zeros(self.histogram_bin_number)

                self.tt_S_binning_avg = self.tt_S_binning_avg * (1 - 1 / self.counter)

                for x in self.tt_N_directional_measure:
                    self.tt_N_binning[(x - 1) // self.histogram_bin_size] += 1
                for x in self.tt_S_no_gaps:
                    self.tt_S_binning[(x - 1) // self.histogram_bin_size] += 1
                    self.tt_S_binning_avg[(x - 1) // self.histogram_bin_size] += 1 / self.counter

                # Batch all data samples
                self.batcher.batch_all(self)

                self.S_bins_res_acc = np.sum(np.array(self.batcher['tt_S_binning_resonance_batch']),0)  # tt_S_binning_resonance accumulated (sum over the batch)
                self.S_bins_detuned_acc = np.sum(np.array(self.batcher['tt_S_binning_detuned_batch']), 0)  # tt_S_binning_resonance accumulated (sum over the batch)

                self.find_transits_events_spectrum_exp(self.tt_S_binning_resonance, self.N, self.transit_cond)

            # Did user request we change with/without atoms state?
            self.handle_user_atoms_on_off_switch()

            # Did we complete N successful iterations?
            if self.counter == self.N+1:
                break

            # We completed another repetition.
            self.repetitions += 1

        ############################################## END MAIN-LOOP #################################################

        # TODO: move all the below to post_run()

        if self.runs_status == TerminationReason.SUCCESS:
            self.logger.info(f'Finished {self.N} Runs, {"with" if self.with_atoms else "without"} atoms')
        elif self.runs_status == TerminationReason.PLAYBACK_END:
            self.logger.info(f'Playback completed. Finished {self.counter} Runs, {"with" if self.with_atoms else "without"} atoms')
        elif self.runs_status == TerminationReason.USER:
            self.logger.info(f'User Terminated the measurements after {self.counter} Runs, {"with" if self.with_atoms else "without"} atoms')
        elif self.runs_status == TerminationReason.ERROR:
            self.logger.info(f'Error Terminated the measurements after {self.counter} Runs, {"with" if self.with_atoms else "without"} atoms')

        # Adding comment to measurement [prompt whether stopped or finished regularly]
        if self.exp_flag:
            if self.counter < self.N:
                aftComment = pymsgbox.prompt('Add comment to measurement: ', default='', timeout=int(30e3))
            else:
                aftComment = ''
        else:
            aftComment = 'ignore'


        experiment_comment = f'Transit condition: {self.transit_cond}\nWhich means - for {len(self.transit_cond)} consecutive on-resonance pulses at least 2 of the conditions must apply.\n Each element represents the minimum number of photon required per pulse to be regarded as transit.'
        daily_experiment_comments = self.generate_experiment_summary_line(self.pre_comment, aftComment, self.with_atoms, self.counter)

        self.save_experiment_results(experiment_comment, daily_experiment_comments)

        # End of Experiment! Ensure we turn the experiment flag off and return the MOT
        self.updateValue("Experiment_Switch", False)
        self.MOT_switch(with_atoms=True, update_parameters=True)

        # Close streams writer
        self.bdstreams.close_file()

        return self.runs_status

    def save_experiment_results(self, experiment_comment, daily_experiment_comments):
        """
        This method is responsible for saving all the results gathered in the experiment and required by analysis
        """

        # Save Quad RF controllers commands
        # TODO: re-implement in QuadRF class - to get the data - and BDResults will save...
        if not self.playback['active']:
            dirname = self.bd_results.get_root()
            for qrdCtrl in self.QuadRFControllers:
                qrdCtrl.saveLinesAsCSV(f'{dirname}\\QuadRF_table.csv')

        # Save all other files
        results = {
            "early_sequence": Config.QRAM_Exp_Square_samples_Early,  # TODO: Why are these called "QRAM" - need to rename them
            "late_sequence": Config.QRAM_Exp_Square_samples_Late,  # TODO: Why are these called "QRAM" - need to rename them
            "north_sequence": Config.QRAM_Exp_Gaussian_samples_N,
            "south_sequence": Config.QRAM_Exp_Gaussian_samples_S,
            "fs_sequence": Config.QRAM_Exp_Square_samples_FS,
            "all_transits_index_batch": self.batcher['all_transits_index_batch'],
            "tt_measure_batch": self.batcher['tt_measure_batch'],
            "tt_N_measure_batch": self.batcher['tt_N_measure_batch'],
            "tt_S_measure_batch": self.batcher['tt_S_measure_batch'],
            "tt_FS_measure_batch": self.batcher['tt_FS_measure_batch'],
            "tt_BP_measure_batch": self.batcher['tt_BP_measure_batch'],
            "tt_DP_measure_batch": self.batcher['tt_DP_measure_batch'],
            "cavity_spectrum": self.Cavity_spectrum,
            "cavity_atom_spectrum": self.Cavity_atom_spectrum,
            "frequency_vector": self.freq_bins,
            #"quadrf_table_ch1": self.QuadRFControllers[0].get_channel_data(0),
            "FLR_measurement": self.batcher['flr_batch'],
            "lock_error": self.batcher['lock_err_batch'],
            "exp_timestr": self.batcher['exp_timestr_batch'],
            "exp_comment": experiment_comment,
            "daily_experiment_comments": daily_experiment_comments,
            #"max_probe_counts": "TBD",  # TODO: ...
            "experiment_config_values": self.Exp_Values
        }
        self.bd_results.save_results(results)


    def is_acquired(self):
        """
        Acquired means that all conditions are met for the experiment to be valid for saving data.

        So we check the following:
        - Locking of resonator is occurring
        - Fluorescence is in accordance to threshold when there are atoms (if no atoms, we don't care)
        - The detectors are functioning and not latched
        - The sum-for-threshold is within the params set for the experiment
        """

        # If we are not in experiment, we don't care about the rest - just say we're acquired.
        if not self.exp_flag:
            return True

        # Are we properly locked on resonance?
        if self.lock_err > self.lock_err_threshold:
            return False

        # Is fluorescence strong enough? (assuming we're running with atoms)
        if self.fluorescence_average < self.FLR_threshold and self.with_atoms:
            return False

        # Are any of the detectors latched?
        if self.latched_detectors():
           return False

        # Is the sum-for-threshold outside the params set for the experiment
        if self.sum_for_threshold > self.total_counts_threshold:
            return False

        # All is well!
        return True

    def handle_warm_up_phase(self):

        if self.post_warm_up_completed:
            return False

        # Get data from OPX streams
        self.get_results_from_streams()

        # Make something out of the data we received on the streams
        self.ingest_time_tags()

        self.tt_S_binning = np.zeros(self.histogram_bin_number * 2)

        for x in self.tt_S_no_gaps:
            self.tt_S_binning[(x - 1) // Config.frequency_sweep_duration] += 1
        self.tt_S_binning_resonance = [self.tt_S_binning[x] for x in range(len(self.tt_S_binning)) if
                                       not x % 2]  # even

        # TODO: we may want to conditon this with the exp_flag (like in original code)
        self.sum_for_threshold = (np.sum(self.tt_S_binning_resonance) * 1000) / self.M_time

        # Check if conditions have been met for the completion of the warm-up phase
        warm_up_phase_complete = self.is_warm_up_phase_complete()

        # This is the first time we realize we are no loger in warm-up - so run post stage
        if warm_up_phase_complete:
            self.post_warm_up()
            self.post_warm_up_completed = True  # Mark the fact we are done

        return not warm_up_phase_complete

    # TODO: make generic in superclass
    def is_warm_up_phase_complete(self):
        """
        We are in warm-up if these conditions are met:
        - We haven't yet completed WARMUP_CYCLES
        - We haven't yet acquired a lock on resonance
        - We are in experiment, and the threshold is not yet filled
        """

        # Did we finish going through all warm-up cycles? If not, we're still in warm-up -> return True
        if self.warm_up_cycles > 0:
            self.warm_up_cycles -= 1
            return False

        # We stay in warm-up while we're not locked
        if self.lock_err > self.lock_err_threshold:
            return False

        # We stay in warm-up if we're not within threshold
        if self.exp_flag and self.sum_for_threshold >= self.total_counts_threshold:
            return False

        return True

    def post_warm_up(self):

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

        # TODO: Complete the stuff in OPX_control_with_QuadRF_Transit_Exp_QRAMconfig.py
        # TODO: from lines 2715 and probably to 2753
        # TODO: until this line: "find_transits_events_spectrum_exp(...)"
        pass

    def maximize_figure(self):
        """
        Attempt to maximize figure based on the matplotlib backend engine. It is different per backed...
        This is taken from StackOverflow - https://stackoverflow.com/questions/12439588/how-to-maximize-a-plt-show-window
        """
        figure_manager = plt.get_current_fig_manager()
        backend_name = matplotlib.get_backend()
        if backend_name == 'QT4Agg':
            figure_manager.window.showMaximized()
        elif backend_name == 'TkAgg':
            figure_manager.window.state('zoomed')
        elif backend_name == 'wxAgg':
            figure_manager.frame.Maximize(True)
        else:
            self.logger.warn(f'Unknown matplotlib backend ({backend_name}). Cannot maximize figure')

    def generate_experiment_summary_line(self, pre_comment, aft_comment, with_atoms, counter):
        time_str = time.strftime("%H%M%S")
        date_str = time.strftime("%Y%m%d")
        comment = None
        if pre_comment is not None:
            comment = pre_comment + '; '
        if aft_comment is not None:
            comment = aft_comment
        if pre_comment is None and aft_comment is None:
            comment = 'No comment.'
        experiment_success = 'ignore' if 'ignore' in comment else 'valid'
        full_line = f'{date_str},{time_str},{experiment_success},{with_atoms},{counter},{comment}'
        return full_line

        # ------------------ end of saving section -------

    def pre_run(self, sequence_definitions, run_parameters):

        super().pre_run(sequence_definitions, run_parameters)

        # Append the pre-comment - the with_atoms parameter
        suffix = ' with atoms' if run_parameters['with_atoms'] else ' without atoms'
        run_parameters['pre_comment'] += suffix

        # Turn Experiment flag on
        self.Experiment_Switch(True)

        # Turn MOT on/off according to parameter
        self.MOT_switch(run_parameters['with_atoms'])

        self.update_parameters()
        pass

    def post_run(self, sequence_definitions, run_parameters):

        # The experiment is done - we turn the experiment flag off (for OPX code)
        self.updateValue("Experiment_Switch", False)

        # Ensure we turn the MOT back on
        self.MOT_switch(with_atoms=True, update_parameters=True)

        pass


if __name__ == "__main__":

    # Playback definitions
    playback_parameters = {
        "active": False,
        'playback_files_path': r"C:\temp\playback_data\PNSA\165103_Photon_TimeTags\playback",
        "old_format": False,
        "save_results": True,
        "save_results_path": 'C:\\temp\\playback_data',
        "max_iterations": 5,  # -1 if you want playback to run through entire data
        "plot": "LIVE",  # "LIVE", "LAST", "NONE"
        "delay": -1,  # -1,  # 0.5,  # In seconds. Use -1 for not playback delay
    }

    run_parameters = {
        'N': 500,
        'transit_profile_bin_size': 100,  # TODO: Q: should this be 10 or 100?
        'pre_comment': 'Spectrum Experiment. No pre-comments defined',
        'transit_cond': [2, 1, 2],
        'total_counts_threshold': 0.4,
        'transit_counts_threshold': 3,
        'FLR_threshold': 0.03,
        'lock_err_threshold': 0.002,
        'Exp_flag': False,
        'with_atoms': True
    }
    sequence_definitions = {
        'total_cycles': 1,
        'delay_between_cycles': None,  # seconds
        'sequence': [
            {
                'name': "With atoms",
                'parameters': {
                    'N': 500,
                    'with_atoms': True
                }
            },
            # {
            #     'name': 'Without Atoms',
            #     'parameters': {
            #         'N': 50,
            #         'with_atoms': False
            #     }
            # }
        ]
    }

    experiment = SpectrumExperiment(playback=False, save_raw_data=False)
    run_status = experiment.run_sequence(sequence_definitions, run_parameters)

    pass
