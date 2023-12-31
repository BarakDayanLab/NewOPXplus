from Experiments.BaseExperiment.BaseExperiment import BaseExperiment
from Experiments.BaseExperiment.BaseExperiment import TerminationReason
from Experiments.QRAM import Config_Experiment as Config

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import time
import math
import pymsgbox
import traceback

from Utilities.BDSound import SOUNDS
from Utilities.Utils import Utils


class QRAMExperiment(BaseExperiment):
    def __init__(self):
        # Invoking BaseClass constructor. It will initiate OPX, QuadRF, BDLogger, Camera, BDResults, KeyEvents etc.
        super().__init__()
        pass

    def __del__(self):
        print('**** SpectrumExperiment Destructor ****')
        super(type(self), self).__del__()
        pass

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
        # self.pgc_aom_chirp_rate = int(self.Exp_Values['PGC_final_Delta_freq'] * 1e3 / (self.Exp_Values[

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
        # TODO: remove? Not used in code (only setting it, not reading from it)
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
        pass

    # ## Boolean functions: ##
    def update_pgc_final_amplitude(self, final_amplitude):
        x = np.log(int(1 / (1 - final_amplitude)))
        self.logger.info(f'Called update_pgc_final_amplitude {final_amplitude}')
        # self.update_lin_pgc_final_amplitude(final_amplitude)
        # self.update_exp_pgc_final_amplitude(final_amplitude)
        return x

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
        self.logger.debug(Probe_counts_South)

        for n in range(repetitions - 1):
            while np.sum(Probe_S_res.tolist()) == Probe_counts_South[-1]:
                Probe_N_res = Probe_N_handle.fetch_all()
                Probe_S_res = Probe_S_handle.fetch_all()
            Probe_counts_North.append(np.sum(Probe_N_res.tolist()))
            Probe_counts_South.append(np.sum(Probe_S_res.tolist()))

        self.logger.debug(
            r'$Max Probe_N = %.1f$' % ((np.average(Probe_counts_North) * 1000) / self.M_time,) + '[MPhotons/sec]')
        self.logger.debug(
            r'$Max Probe_S = %.1f$' % ((np.average(Probe_counts_South) * 1000) / self.M_time,) + '[MPhotons/sec]')

        self.Max_Probe_counts_switch(False)
        self.update_parameters()

        return [(np.average(Probe_counts_North) * 1000) / self.M_time,
                (np.average(Probe_counts_South) * 1000) / self.M_time]

    def ingest_time_tags(self):
        """
        Takes all the raw results we got from the streams and does some processing on them - preparing "measures"
        """

        # Get fluorescence data
        self.FLR_res = -self.streams['FLR_measure']['results']
        self.fluorescence_average = 1000 * self.FLR_res

        # Clear time-tags vectors from last data
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

        # TODO: Q: the below unification means that if 2 detectors clicked, we may have double appearings
        # Unify detectors 1 & 2
        self.tt_BP_measure = sorted(sum(self.tt_measure[0:2], []))
        # Unify detectors 3 & 4
        self.tt_DP_measure = sorted(sum(self.tt_measure[2:4], []))
        # Take detector 8
        self.tt_N_measure = sorted(sum(self.tt_measure[7:], []))
        # Detector 4
        self.tt_S_measure = sorted(sum(self.tt_measure[4:5], []))
        # Unify detectors 6 & 7
        self.tt_FS_measure = sorted(sum(self.tt_measure[5:7], []))

        # plt.plot(self.tt_S_measure)

        # Phase Correction is a result of ZIP action in opx_control, thus we have "value_0" and "value_1" for the tupples
        self.Phase_Correction_vec = self.streams['Phase_Correction_array']['results']['value_0']
        self.Phase_Correction_min_vec = self.streams['Phase_Correction_array']['results']['value_1']
        self.Phase_Correction_value = self.streams['Phase_Correction_array']['results']['value_1'][-1]
        # self.Phase_diff = self.streams['Phase_diff']['results']
        self.MZ_S_tot_counts = self.streams['Max_counts']['results']

        # TODO: can we give these a more meaningful name - what is 'value_0' and 'value_1' ?
        self.MZ_BP_counts_res_value_0 = self.streams['Bright_Port_Counts']['results']['value_0']  # batch
        self.MZ_BP_counts_res_value_1 = self.streams['Bright_Port_Counts']['results']['value_1']  # check batch

        self.MZ_DP_counts_res_value_0 = self.streams['Dark_Port_Counts']['results']['value_0']  # batch
        self.MZ_DP_counts_res_value_1 = self.streams['Dark_Port_Counts']['results']['value_1']  # check batch

        pass


    # TODO: Q: this currently does nothing. What was the intention and why isn't it in use?
    def latched_detectors(self):
        latched_detectors = []
        # for indx, det_tt_vec in enumerate(self.tt_measure):  # for different detectors
        #     if not det_tt_vec:
        #         latched_detectors.append(indx)
        return latched_detectors

    def get_pulses_bins(self, sprint_pulse_len, num_of_det_pulses, num_of_sprint_pulses,
                        sprint_sequence_delay, num_of_sprint_sequences, num_init_zeros, num_fin_zeros,
                        num_between_zeros):
        '''
        Generating bin vector with edges at the edges of pulses (such that the histogram would give the number of
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
            pulses_bins += (pulses_bins[-1] + np.arange(self.det_pulse_len, (num_of_det_pulses + 1) * self.det_pulse_len,
                                                        self.det_pulse_len)).tolist()
            pulses_bins += (pulses_bins[-1] + np.arange(sprint_pulse_len, (num_of_sprint_pulses + 1) * sprint_pulse_len,
                                                        sprint_pulse_len)).tolist()
            pulses_bins[-1] += num_fin_zeros + num_init_zeros - num_between_zeros
        return pulses_bins

    def divide_to_reflection_trans(self, sprint_pulse_len, num_of_det_pulses, num_of_sprint_pulses, num_of_sprint_sequences):
        '''
        dividing south and north tt's vectors into reflection and transmission vectors, by:
        1. converting them into counts vector in same time spacing as the initial pulse.
        2. taking the relevant pulses from the sequence to each vector (reflection/transmission)
        :return:
        '''
        # create histogram with self.M_window bins
        self.pulses_bins_S = self.get_pulses_bins(sprint_pulse_len, num_of_det_pulses,
                                                  num_of_sprint_pulses, self.delay_in_detection_S,
                                                  num_of_sprint_sequences, Config.num_init_zeros_S,
                                                  Config.num_fin_zeros_S, Config.num_between_zeros)
        self.pulses_bins_N = self.get_pulses_bins(sprint_pulse_len, num_of_det_pulses,
                                                  num_of_sprint_pulses, self.delay_in_detection_N,
                                                  num_of_sprint_sequences, Config.num_init_zeros_N,
                                                  Config.num_fin_zeros_N, Config.num_between_zeros)
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

    def divide_tt_to_reflection_trans(self):
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
            # TODO: Q: I assume this is a time-window to ignore some time on start/end, how did we decide on this time?  TODO: Q: what us the meaning of this specific time-window?
            if (element > int(0.6e6)) and (element < int(9.6e6)):
                tt_inseq = element % self.QRAM_sequence_len
                if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                    seq_num = (element - 1) // self.QRAM_sequence_len
                    self.num_of_det_reflections_per_seq_S[seq_num] += self.filter_S[tt_inseq]
                    self.num_of_det_transmissions_per_seq_N[seq_num] += self.filter_N[tt_inseq]
            # else:  # The part of the SPRINT pulses in the sequence
            #     SPRINT_pulse_num = (tt_inseq - self.end_of_det_pulse_in_seq) // (self.sprint_pulse_len + Config.num_between_zeros)
            #     if SPRINT_pulse_num < self.number_of_SPRINT_pulses_per_seq:
            #         self.num_of_SPRINT_reflections_per_seq_S[(element-1) // self.QRAM_sequence_len][SPRINT_pulse_num] += 1
            #         self.num_of_SPRINT_transmissions_per_seq_N[(element-1) // self.QRAM_sequence_len][SPRINT_pulse_num] += 1
            #         tt_small_perturb+=[element]

        for element in self.tt_S_measure + self.tt_FS_measure:
            if (element > int(0.6e6)) and (element < int(9.6e6)):
                tt_inseq = element % self.QRAM_sequence_len
                if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                    seq_num = (element - 1) // self.QRAM_sequence_len
                    self.num_of_det_reflections_per_seq_N[seq_num] += self.filter_N[tt_inseq]
                    self.num_of_det_transmissions_per_seq_S[seq_num] += self.filter_S[tt_inseq]
            # else:  # The part of the SPRINT pulses in the sequence
            #     SPRINT_pulse_num = (tt_inseq - self.end_of_det_pulse_in_seq) // (self.sprint_pulse_len + Config.num_between_zeros)
            #     if SPRINT_pulse_num < self.number_of_SPRINT_pulses_per_seq:
            #         self.num_of_SPRINT_reflections_per_seq_N[(element-1) // self.QRAM_sequence_len][SPRINT_pulse_num] += 1
            #         self.num_of_SPRINT_transmissions_per_seq_S[(element-1) // self.QRAM_sequence_len][SPRINT_pulse_num] += 1

    """
        Dor: "Method was built so we can debug coherence"
        In this method we count the clicks on the Dark-Port and Bright-Port after the detection pulses

        Sequences:19231, Sequence-Len:520, End-of-detection-Pulse:510
        +-------------+--------------+--------------+--------------+--------+---+
        +      0      |      19      |      211     |      23      |    54  |   | ...  +
        +-------------+--------------+--------------+--------------+--------+---+
        <- 520 [ns] -> <- 520 [ns] -> <- 520 [ns] -> <- 520 [ns] -> <-510-><-10->     
    """

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
                try:
                    self.num_of_DP_counts_per_n_sequences[
                        (element - 1) // (self.QRAM_sequence_len * num_of_seq_per_count)] += \
                        np.ceil(Config.QRAM_Exp_Square_samples_Late[tt_inseq])
                except Exception as err:
                    print(err)

        for element in self.tt_FS_measure + self.tt_S_measure:
            tt_inseq = element % self.QRAM_sequence_len
            if tt_inseq > self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                try:
                    # lll = (element - 1) // (self.QRAM_sequence_len * num_of_seq_per_count)
                    # if (element - 1) // (self.QRAM_sequence_len * num_of_seq_per_count) > 385:
                    #     lll = 385
                    # self.num_of_S_counts_per_n_sequences[lll] += np.ceil(Config.QRAM_Exp_Square_samples_Late[tt_inseq])
                    self.num_of_S_counts_per_n_sequences[
                        (element - 1) // (self.QRAM_sequence_len * num_of_seq_per_count)] += \
                        np.ceil(Config.QRAM_Exp_Square_samples_Late[tt_inseq])
                except Exception as err:
                    print(err)

    def plot_folded_tt_histogram(self):

        plt.figure()
        plt.plot(sum(np.reshape(self.tt_histogram_N, [int(len(self.tt_histogram_N) / 9), 9])), label='tt_hist_N')
        plt.plot(sum(np.reshape(self.tt_histogram_S, [int(len(self.tt_histogram_S) / 9), 9])), label='tt_hist_S')
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
        detection_pulse_range = np.zeros([self.number_of_QRAM_sequences, num_of_det_pulses], dtype=int)
        sprint_pulse_range = np.zeros([self.number_of_QRAM_sequences, num_of_sprint_pulses], dtype=int)

        # define ranges for tt_histogram (a vector that counts number of photons in each pulse)
        for i in range(self.number_of_QRAM_sequences):
            detection_pulse_range[i] = \
                list(range(i * num_of_pulses, i * num_of_pulses + num_of_det_pulses))
            # sprint pulse range starts from last detection pulse andfo num_of_sprint_pulses
            sprint_pulse_range[i] = \
                list(range(int(detection_pulse_range[i][-1]) + 1,
                           int(detection_pulse_range[i][-1]) + 1 + num_of_sprint_pulses))
        # find transits and sprint events
        for j in range(self.number_of_QRAM_sequences - 2):
            if \
                    sum(self.tt_histogram_reflection[detection_pulse_range[j]]) >= detection_condition[0] and \
                            sum(self.tt_histogram_reflection[detection_pulse_range[j + 1]]) >= detection_condition[1]:
                num_of_detected_atom += 1
                transit_sequences_init_tt.append(j * self.QRAM_sequence_len)
                if self.tt_histogram_reflection[sprint_pulse_range[j][0]]:
                    num_of_sprints += 1
                    sprints_events.append(self.tt_histogram_reflection[sprint_pulse_range[j]])
        sprints_data = [sprints_events, num_of_sprints]
        atom_detect_data = [transit_sequences_init_tt, num_of_detected_atom]
        return atom_detect_data, sprints_data

    def find_transits_and_sprint_events_changed(self, cond=None, minimum_number_of_seq_detected=2):
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

        if cond is None:
            cond = [2, 2, 2]

        current_transit = []
        self.all_transits_seq_indx = []  # Array of the sequence indexes of all recognized transits per cycle. The length of it will be the number of all transits at the current cycle.
        self.reflection_SPRINT_data_per_transit = []  # Array of vectors with data on the number of reflections per SPRINT pulse in sequence.
        self.transmission_SPRINT_data_per_transit = []  # Array of vectors with data on the number of transmissions per SPRINT pulse in sequence.

        for i in range(len(self.num_of_det_reflections_per_seq[:]) - len(cond) + 1):
            cond_check = (self.num_of_det_reflections_per_seq[i:(i + len(cond))] >= cond).astype(int)
            if sum(cond_check) >= minimum_number_of_seq_detected:
                current_transit = np.unique(
                    current_transit + [*range(i + np.where(cond_check != 0)[0][0], (i + len(cond)))]).tolist()
            elif len(current_transit) > 1:
                current_transit = current_transit[
                                  :np.where(self.num_of_det_reflections_per_seq[current_transit] >= min(cond))[0][
                                       -1] + 1]
                if self.all_transits_seq_indx:
                    if bool(set(current_transit) & set(self.all_transits_seq_indx[-1])):
                        current_transit = self.all_transits_seq_indx[-1] + current_transit[1:]
                        self.all_transits_seq_indx = self.all_transits_seq_indx[:-1]
                        # self.reflection_SPRINT_data_per_transit = self.reflection_SPRINT_data_per_transit[:-1]
                        # self.transmission_SPRINT_data_per_transit = self.transmission_SPRINT_data_per_transit[:-1]
                self.all_transits_seq_indx.append(current_transit)
                # self.reflection_SPRINT_data_per_transit.append([self.num_of_SPRINT_reflections_per_seq[elem].tolist()
                #                                                 for elem in current_transit[:-1]])
                # self.transmission_SPRINT_data_per_transit.append([self.num_of_SPRINT_transmissions_per_seq[elem].tolist()
                #                                                   for elem in current_transit[:-1]])
                current_transit = []
        if len(current_transit) > 1:
            current_transit = current_transit[
                              :np.where(self.num_of_det_reflections_per_seq[current_transit] >= min(cond))[0][-1] + 1]
            if self.all_transits_seq_indx:
                if bool(set(current_transit) & set(self.all_transits_seq_indx[-1])):
                    current_transit = self.all_transits_seq_indx[-1] + current_transit[1:]
                    self.all_transits_seq_indx = self.all_transits_seq_indx[:-1]
                    # self.reflection_SPRINT_data_per_transit = self.reflection_SPRINT_data_per_transit[:-1]
                    # self.transmission_SPRINT_data_per_transit = self.transmission_SPRINT_data_per_transit[:-1]
            self.all_transits_seq_indx.append(current_transit)
            # self.reflection_SPRINT_data_per_transit.append([self.num_of_SPRINT_reflections_per_seq[elem].tolist()
            #                                                 for elem in current_transit[:-1]])
            # self.transmission_SPRINT_data_per_transit.append([self.num_of_SPRINT_transmissions_per_seq[elem].tolist()
            #                                                   for elem in current_transit[:-1]])

    def get_pulses_location_in_seq(self, delay, seq, smearing):
        '''
        :description
        A function that takes a pulse sequence, (a) applies a delay to it, (b) detects the pulses inside and
        (c) widens them pulses by the "smearing" factor, to build a filter (done to match the filter to the performance
        of the physical system, e.g. AOMs).

        :param delay: Between the actual sequence pulse location to the location of the folded data
        :param seq: The sequence of pulses from which the filter is generated.
        :param smearing: The value that is added to the filter before and after each pulse in the sequence.

        :return: the pulses location in the sequence (indices tuples) and the filter signal itself
        '''

        # Create an array consisting of zeros/ones that acts as filter, based on sequence values.
        seq_filter = (np.array(seq) > 0).astype(int)
        # Shift all values ("roll") to the right, based on the delay parameter
        seq_filter = np.roll(seq_filter, delay)
        # Find signal boundries indices - comparing each element with its neighbor
        indices = np.where(seq_filter[:-1] != seq_filter[1:])[0] + 1
        # Create the pulses boundry tuples - with the smearing widening
        pulses_loc = [(indices[i] - smearing, indices[i + 1] + smearing) for i in range(0, len(indices), 2)]
        # Recreate the signal by filling it with 1's in the relevant places
        seq_filter_with_smearing = np.zeros(seq_filter.shape[0])
        for (start, end) in pulses_loc: np.put_along_axis(seq_filter_with_smearing, np.arange(start, end), 1, axis=0)
        return pulses_loc, seq_filter_with_smearing

    def get_pulses_location_in_seq_DEP(self, delay, seq=Config.QRAM_Exp_Gaussian_samples_S,
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
        # Create an array consisting of zeros/ones that acts as filter, based on sequence values.
        # Ones - if value is positive, Zero - if it is non-positive
        seq_filter = (np.array(seq) > 0).astype(int)

        # Shift all values ("roll") to the right, based on the delay parameter
        seq_filter = np.roll(seq_filter, delay)

        # Find the indices where the filter is 1 (e.g value was positive)
        seq_index = np.where(seq_filter > 0)[0]
        seq_filter_with_smearing = seq_filter
        seq_filter_with_smearing = np.copy(seq_filter)
        pulses_loc = []
        if seq_index.any():
            start_index = seq_index[0]

            # Iterate over all the indices - check if each value has "gap" to the previous one
            for i in range(1, len(seq_index)):
                if seq_index[i] - seq_index[i - 1] > 1:
                    index_range = (start_index - int(smearing), seq_index[i - 1] + int(smearing))
                    pulses_loc.append(index_range)
                    for j in index_range:
                        seq_filter_with_smearing[j] = 1
                    # Advance the start index to the end of the "gap"
                    start_index = seq_index[i]
            pulses_loc.append((start_index - int(smearing), seq_index[-1] + int(smearing)))
            for j in range(start_index - int(smearing), seq_index[-1] + int(smearing)):
                seq_filter_with_smearing[j] = 1
        self._plot(seq_filter_with_smearing)
        return pulses_loc, seq_filter_with_smearing

    def get_avg_num_of_photons_in_seq_pulses(self, seq, pulse_loc, tt_measure):
        avg_num_of_photons_in_seq_pulses = []
        try:
            real_number_of_seq = math.ceil(max(tt_measure) / len(Config.QRAM_Exp_Gaussian_samples_S))
            # self.logger.debug('Real number of seq = %d' %real_number_of_seq)
        except:
            real_number_of_seq = self.number_of_QRAM_sequences
            # self.logger.debug('Max number of seq')
        for t in pulse_loc:
            avg_num_of_photons_in_seq_pulses.append((sum(seq[t[0]:t[1]]) + seq[t[1]]) / (
                        real_number_of_seq * 0.167))  # Sagnac configuration efficiency 16.7%
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
        '''
        Finds the average value of the pulse indices - so it's the middle position. Used for the UI/UX

        :param pulse_loc:
        :return: returns an array - averages of the values of the pulses indices
        '''
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

    # TODO: move this to Base or Utils (make it static)
    def merge_playback_data_files(self, file_1, file_2):
        playback_folder = self.bd_results.get_custom_root('playback_data')

        data_file_1 = os.path.join(playback_folder, file_1)
        zipped_data_1 = np.load(data_file_1, allow_pickle=True)
        data_1 = zipped_data_1['arr_0']
        data_file_2 = os.path.join(playback_folder, file_2)
        zipped_data_2 = np.load(data_file_2, allow_pickle=True)
        data_2 = zipped_data_2['arr_0']

        if len(data_1) != len(data_2):
            raise Exception(f'Found files with different length - cannot construct playback! ({file_1}')

        merged_array = []
        for i in range(0, len(zipped_data_1['arr_0'])):
            merged_array.append({
                "value_0": data_1[i],
                "value_1": data_2[i]
            })
        return merged_array

    def load_data_for_playback(self):
        """
        This function loads relevant data saved in previous experiments
        and re-constructs the original streams raw data.
        (the raw data will serve for the playback runs)
        """
        playback_folder = self.bd_results.get_custom_root('playback_data')

        # TODO: try/raise on all - to catch cases where data file is missing!

        # Iterate over all stream and load detectors data from the relevant file
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

        # Load fluorescence
        data_file = os.path.join(playback_folder, 'Flouresence.npz')
        zipped_data = np.load(data_file, allow_pickle=True)
        self.streams['FLR_measure']['all_rows'] = (-zipped_data['arr_0']).tolist()

        # Load Dark-Port/Bright-Port/Phase-Correction data from separate files and combine them
        self.streams['Bright_Port_Counts']['all_rows'] = \
            self.merge_playback_data_files(file_1='BalancingRes\\MZ_BrightPort_counts_balancing_batch.npz',
                                           file_2='BalancingRes\\MZ_BrightPort_counts_balancing_check_batch.npz')
        self.streams['Dark_Port_Counts']['all_rows'] = \
            self.merge_playback_data_files(file_1='BalancingRes\\MZ_DarkPort_counts_balancing_batch.npz',
                                           file_2='BalancingRes\\MZ_DarkPort_counts_balancing_check_batch.npz')
        self.streams['Phase_Correction_array']['all_rows'] = \
            self.merge_playback_data_files(file_1='BalancingRes\\MZ_balancing_Phase_Correction_vec_batch.npz',
                                           file_2='BalancingRes\\MZ_balancing_Phase_Correction_min_vec_batch.npz')

        # There is no file that saved this information, so we just put zeros here
        self.streams['Max_counts']['all_rows'] = np.zeros(shape=len(zipped_data['arr_0']), dtype=int)

    def print_experiment_information(self):
        """
        Prints during the experiment a formatted line that brings all concise information
        """
        locking_efficiency = self.counter / self.repetitions
        locking_efficiency_str = '%.2f' % locking_efficiency
        fluorescence_str = '%.2f' % self.fluorescence_average
        time_formatted = time.strftime("%Y/%m/%d, %H:%M:%S")
        status_str = f'[Warm Up: {self.warm_up_cycles}]' if self.warm_up else f'# {self.counter} ({self.repetitions})'
        photons_str = '' if self.warm_up else f', #Photons/[us]: {self.sum_for_threshold}'

        self.logger.info(f'{time_formatted}: {status_str}, Eff: {locking_efficiency_str}, Flr: {fluorescence_str}{photons_str}')

    def experiment_mainloop_delay(self):
        if self.playback['active']:
            if self.playback['delay'] != -1:
                time.sleep(self.playback['delay'])
        else:
            time.sleep(0.01)  # TODO: do we need this delay?

    def prepare_figures(self):
        self.fig = plt.figure()
        self.subplots = []
        self.plot_shown = False

        # Figure 0 - Top Left - Binned Time-Tags from All Detectors
        self.subplots.append(plt.subplot2grid((3, 4), (0, 0), colspan=2, rowspan=1))

        # Figure 1 - Top Right - Binned Time-Tags from All Detectors Folded - Averaged
        self.subplots.append(plt.subplot2grid((3, 4), (0, 2), colspan=2, rowspan=1))

        # Figure 2 - Middle Left - MZ outputs while locking + Detectors Circles
        self.subplots.append(plt.subplot2grid((3, 4), (1, 0), colspan=2, rowspan=1))

        # Figure 3 - Num of reflections per sequence
        self.subplots.append(plt.subplot2grid((3, 4), (1, 2), colspan=2, rowspan=1))

        # Figure 4 - Bottom Left (left) - MZ outputs around experiment
        self.subplots.append(plt.subplot2grid((3, 4), (2, 0), colspan=1, rowspan=1))

        # Figure 5 - Bottom Right - Transits per sequence
        self.subplots.append(plt.subplot2grid((3, 4), (2, 2), colspan=2, rowspan=1))

        # Figure 6 - Middle Left - PLACED on top of Graph 2 - same X, different Y -
        self.subplots.append(self.subplots[2].twinx())

        # Figure 7 - Bottom Left (right) - MZ outputs during experiment
        self.subplots.append(plt.subplot2grid((3, 4), (2, 1), colspan=1, rowspan=1))

        #self.maximize_figure()  # TODO: uncomment

        return

    def plot_figures(self):
        try:
            self._plot_figures()
        except Exception as err:
            tb = traceback.format_exc()
            self.warn(f'Failed on plot_figures: {err}')
        pass

    def _plot_figures(self):

        if self.playback['active']:
            if self.playback['plot'] == 'NONE':
                return
            if self.playback['plot'] == 'LAST' and self.counter < (self.N - 10):
                return

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
            "graph-5": True,
            "graph-6": False,
            "graph-7": True
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
        ax[6].clear()
        ax[7].clear()

        # Play sound:
        if plot_switches['playsound']:
            if not self.acquisition_flag:
                self.bdsound.play(SOUNDS.INDICATOR_1)

        # Prepare header text
        pause_str = ' , PAUSED!' if self.pause_flag else ''
        ref_str = '%.2f' % self.sum_for_threshold
        flr_str = '%.2f' % self.fluorescence_average
        eff_str = '%.2f' % (self.counter / self.repetitions)
        lck_str = '%.3f' % self.lock_err
        status_str = f'[Warm Up: {self.warm_up_cycles}]' if self.warm_up else f'# {self.counter} ({self.repetitions})'
        playback_str = 'PLAYBACK: ' if self.playback['active'] else ''
        header_text = f'{playback_str} {status_str} - Reflections: {ref_str}, Eff: {eff_str}, Flr: {flr_str}, Lock Error: {lck_str} {pause_str}'

        # Threshold Box:
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        flag_color = 'green' if self.acquisition_flag else 'red'
        props_thresholds = dict(boxstyle='round', edgecolor=flag_color, linewidth=2, facecolor=flag_color, alpha=0.5)

        avg_BP = np.average(self.num_of_BP_counts_per_n_sequences)
        avg_DP = np.average(self.num_of_DP_counts_per_n_sequences)

        # TODO: temp workaround, REMOVE
        avg_BP = 0.001 if avg_BP == 0 else avg_BP
        textstr_BP_DP = 'Average Bright counts = %.2f \n' % (avg_BP,) + \
                        'Average Dark counts = %.2f \n' % (avg_DP,) + \
                        'Infidelity = %.2f' % (avg_DP / (avg_DP + avg_BP),)

        textstr_BP_DP_BA = 'Infidelity before = %.2f \n' % (self.Infidelity_before,) + \
                           'Infidelity after= %.2f \n' % (self.Infidelity_after,) + \
                           'S total counts MZ = %.2f ' % (self.MZ_S_tot_counts,)
        # 'phase difference S-N = %.2f \n' % (self.Phase_diff,) + \
        # 'S-N total counts ratio = %.2f \n' % (self.MZ_S_tot_counts,)

        # Binned time-tags from all detectors folded (live)
        if plot_switches['graph-0']:
            ax[0].plot(self.folded_tt_N_directional, label='"N" detectors')
            ax[0].plot(self.folded_tt_S_directional, label='"S" detectors')
            ax[0].plot(self.folded_tt_BP_timebins, label='"BP" detectors')
            ax[0].plot(self.folded_tt_DP_timebins, label='"DP" detectors')
            ax[0].plot((self.filter_S) * max(self.folded_tt_N_directional + self.folded_tt_S_directional),
                       '--', color='orange', label='Filter "S"')
            for i in range(len(self.Num_of_photons_txt_box_y_loc_live)):
                ax[0].text(self.Num_of_photons_txt_box_x_loc.tolist()[i], self.Num_of_photons_txt_box_y_loc_live[i],
                           '%.2f' % self.avg_num_of_photons_per_pulse_live[i],
                           horizontalalignment='center', fontsize=12, fontweight='bold', family=['Comic Sans MS'])
            for j in range(len(self.Num_of_photons_txt_box_y_loc_live_MZ)):
                ax[0].text(self.Num_of_photons_txt_box_x_loc_for_MZ_ports.tolist()[j],
                           self.Num_of_photons_txt_box_y_loc_live_MZ[j][0],
                           '%.2f' % self.avg_num_of_photons_per_pulse_live_MZ[j][0],
                           horizontalalignment='center', fontsize=12, fontweight='bold', family=['Comic Sans MS'],
                           color='#2ca02c')
                ax[0].text(self.Num_of_photons_txt_box_x_loc_for_MZ_ports.tolist()[j],
                           self.Num_of_photons_txt_box_y_loc_live_MZ[j][1],
                           '%.2f' % self.avg_num_of_photons_per_pulse_live_MZ[j][1],
                           horizontalalignment='center', fontsize=12, fontweight='bold', family=['Comic Sans MS'],
                           color='#d62728')
            # ax[0].set_ylim(0, 0.3)
            ax[0].set_title('Binned time-tags from all detectors folded (Live)', fontweight="bold")
            ax[0].legend(loc='upper right')
            ax[0].text(0.05, 1.4, header_text, transform=ax[0].transAxes, fontsize=26, verticalalignment='top', bbox=props_thresholds)

        # Binned time-tags from all detectors folded (Averaged)
        if plot_switches['graph-1']:
            # ax[1].plot(self.folded_tt_BP_cumulative_avg, label='"BP" detectors')
            # ax[1].plot(self.folded_tt_DP_cumulative_avg, label='"DP" detectors')
            ax[1].plot(self.folded_tt_N_directional_cumulative_avg, label='"N" detectors')
            ax[1].plot(self.folded_tt_S_directional_cumulative_avg, label='"S" detectors')
            ax[1].plot(self.folded_tt_BP_timebins_cumulative_avg, label='"BP" detectors')
            ax[1].plot(self.folded_tt_DP_timebins_cumulative_avg, label='"DP" detectors')
            ax[1].plot((self.filter_N) * max(self.folded_tt_N_directional_cumulative_avg + self.folded_tt_S_directional_cumulative_avg),
                       '--b', label='Filter "N"')
            for i in range(len(self.Num_of_photons_txt_box_y_loc)):
                ax[1].text(self.Num_of_photons_txt_box_x_loc.tolist()[i], self.Num_of_photons_txt_box_y_loc[i],
                           '%.2f' % self.avg_num_of_photons_per_pulse[i],
                           horizontalalignment='center', fontsize=12, fontweight='bold', family=['Comic Sans MS'])
            for j in range(len(self.Num_of_photons_txt_box_y_loc_MZ)):
                ax[1].text(self.Num_of_photons_txt_box_x_loc_for_MZ_ports.tolist()[j],
                           self.Num_of_photons_txt_box_y_loc_MZ[j][0],
                           '%.2f' % self.avg_num_of_photons_per_pulse_MZ[j][0],
                           horizontalalignment='center', fontsize=12, fontweight='bold', family=['Comic Sans MS'],
                           color='#2ca02c')
                ax[1].text(self.Num_of_photons_txt_box_x_loc_for_MZ_ports.tolist()[j],
                           self.Num_of_photons_txt_box_y_loc_MZ[j][1],
                           '%.2f' % self.avg_num_of_photons_per_pulse_MZ[j][1],
                           horizontalalignment='center', fontsize=12, fontweight='bold', family=['Comic Sans MS'],
                           color='#d62728')
            # ax[1].set_ylim(0, 0.3)
            ax[1].set_title('binned time-tags from all detectors folded (Averaged)', fontweight="bold")
            ax[1].legend(loc='upper right')

        # MZ outputs while locking + Detectors Circles
        if plot_switches['graph-2']:
            # self.Phase_Correction_min_diff += [(1 + (self.Phase_Correction_min_vec[-1] - self.Phase_Correction_min_vec[67])) % 1]
            ax[2].plot(self.MZ_BP_counts_res_value_0, label='MZ Bright port')
            ax[2].plot(self.MZ_DP_counts_res_value_0, label='MZ Dark port')
            ax[2].plot(self.MZ_BP_counts_res_value_0 - self.MZ_DP_counts_res_value_0, label='Dif ports')

            ax[6].tick_params(axis="y", labelcolor='#8c564b')
            ax[6].plot(self.Phase_Correction_vec, label='Phase correction values', color='#8c564b')
            ax[6].plot(self.Phase_Correction_min_vec, label='Phase correction values', color='#9467bd')

            ax[2].set_ylim(0, 1.1 * np.max([self.MZ_BP_counts_res_value_0, self.MZ_DP_counts_res_value_0]))
            ax[2].set_title('MZ outputs while locking', fontweight="bold")
            ax[2].legend(loc='upper right')

            # Detectors status:
            if plot_switches['detectors']:
                for i, det in enumerate(self.Num_Of_dets):
                    x = -0.15
                    y = 1.9 - i * 0.4
                    pad = 1
                    text = 'Det %d' % det
                    det_color = 'red' if self.latched_detectors() else 'green'
                    ax[2].text(x, y, text, ha="center", va="center", transform=ax[2].transAxes,
                             bbox=dict(boxstyle=f"circle,pad={pad}", edgecolor=det_color, linewidth=2, facecolor=det_color, alpha=0.5))

        # Num of reflections per sequence
        if plot_switches['graph-3']:
            max_reflect_avg = max(self.num_of_det_reflections_per_seq_accumulated / self.counter)
            max_reflect = max(self.num_of_det_reflections_per_seq)

            textstr_total_reflections = 'Total reflections per cycle "N" = %d \n' % (
                sum(self.num_of_det_reflections_per_seq_N),) \
                                        + 'Total reflections per cycle "S" = %d \n' % (
                                            sum(self.num_of_det_reflections_per_seq_S),) \
                                        + 'Average reflections per cycle = %.2f \n' % (
                                            sum(self.num_of_det_reflections_per_seq_accumulated / self.counter),) \
                                        + 'Average transmissions per cycle = %.2f \n' % (
                                            sum(self.num_of_det_transmissions_per_seq_accumulated / self.counter),) \
                                        + 'Average reflections precentage = %.2f' % (
                                            sum(self.num_of_det_reflections_per_seq_accumulated) /
                                            (sum(self.num_of_det_reflections_per_seq_accumulated) +
                                             sum(self.num_of_det_transmissions_per_seq_accumulated)),)

            ax[3].plot(self.num_of_det_reflections_per_seq_accumulated / self.counter, label='Num of reflections per sequence')

            # TODO: original code did not have this condition. Can it be that max_relect is never zero?
            if max_reflect > 0:
                ax[3].plot(self.num_of_det_reflections_per_seq * 0.5 * max_reflect_avg / max_reflect * 0.3, label='Num of reflections per sequence (Live)')
            ax[3].set_title('Num of reflections per sequence', fontweight="bold")
            ax[3].legend(loc='upper right')
            ax[3].text(0.1, 0.9 * max_reflect_avg, textstr_total_reflections, fontsize=14, verticalalignment='top', bbox=props)

        # MZ outputs around experiment
        if plot_switches['graph-4']:
            ax[4].plot(self.MZ_BP_counts_res_value_1, label='MZ BP counts before and after')
            ax[4].plot(self.MZ_DP_counts_res_value_1, label='MZ DP port counts before and after')
            ax[4].axvline(len(self.MZ_DP_counts_res_value_1) / 2, linestyle='--', color='red')
            ax[4].set_title('MZ outputs around experiment', fontweight="bold")
            # ax[4].legend(loc='upper right')
            ax[4].text(0.05, 0.6, textstr_BP_DP_BA, transform=ax[4].transAxes, fontsize=14, verticalalignment='top', bbox=props)

        # MZ outputs during experiment
        if plot_switches['graph-7']:
            ax[7].plot(self.num_of_BP_counts_per_n_sequences, label='MZ BP counts per %d seq' % 50)
            ax[7].plot(self.num_of_DP_counts_per_n_sequences, label='MZ DP counts per %d seq' % 50)
            ax[7].plot(self.num_of_S_counts_per_n_sequences, label='"S" transmission counts per %d seq' % 50)
            ax[7].plot(
                self.num_of_S_counts_per_n_sequences + self.num_of_DP_counts_per_n_sequences + self.num_of_BP_counts_per_n_sequences,
                label='"S" total counts per %d seq' % 50)
            ax[7].set_title('MZ outputs during experiment', fontweight="bold")
            # ax[7].legend(loc='upper right')
            ax[7].text(0.05, 0.6, textstr_BP_DP, transform=ax[7].transAxes, fontsize=14,
                       verticalalignment='top', bbox=props)

        if self.number_of_transits_live:
            textstr_transit_counts = r'$N_{Transits} = %s $' % (self.number_of_transits_live,) + r'$[Counts]$'
        else:
            textstr_transit_counts = r'$N_{Transits} = %s $' % (0,) + r'$[Counts]$'
        textstr_transit_event_counter = r'$N_{Transits Total} = %s $' % (self.number_of_transits_total,) + '[Counts]\n' \
                                        + textstr_transit_counts
        # + textstr_total_SPRINT_reflections

        # Transits per sequence
        if plot_switches['graph-5']:
            ax[5].plot(range(self.number_of_QRAM_sequences), self.seq_transit_events_batched,
                       label='Transit Events Accumulated')
            ax[5].plot(range(self.number_of_QRAM_sequences), self.seq_transit_events_live, label='Transit Events Live')
            ax[5].set(xlabel='Sequence [#]', ylabel='Counts [Photons]')
            ax[5].set_title('Transits per sequence', fontweight="bold")
            ax[5].legend(loc='upper right')
            ax[5].text(0.05, 0.95, textstr_transit_event_counter, transform=ax[5].transAxes, fontsize=14,
                       verticalalignment='top', bbox=props)

        # End timer
        total_prep_time = time.time() - start_plot_time

        # plt.tight_layout()
        plt.pause(0.2)

    def init_params_for_experiment(self):
        # define empty variables
        self.number_of_QRAM_sequences = math.ceil(self.M_window / self.QRAM_sequence_len)

        # Reformatting the above variables into an object - (a) for better clarity (b) make it "experiment-agnostic"
        self.experiment = {
            "sequence_length": self.QRAM_sequence_len,
            "measurement_window": self.M_window,  # [ns]
            "number_of_sequences": math.ceil(self.M_window / self.QRAM_sequence_len)
        }
        self.number_of_SPRINT_pulses_per_seq = len(Config.sprint_pulse_amp_S)

        # TODO: Do we need these?
        self.tt_measure = []
        self.tt_S_measure = []

        self.folded_transmission = np.zeros(len(Config.QRAM_Exp_Gaussian_samples_S))
        self.folded_reflection = np.zeros(len(Config.QRAM_Exp_Gaussian_samples_S))

        self.tt_S_binning = np.zeros(self.number_of_QRAM_sequences + 1)
        self.seq_transit_events_live = np.zeros(self.number_of_QRAM_sequences)
        self.seq_transit_events_batched = np.zeros(self.number_of_QRAM_sequences)
        self.tt_S_SPRINT_events = np.zeros(self.QRAM_sequence_len)
        self.tt_S_SPRINT_events_batch = np.zeros(self.QRAM_sequence_len)
        self.num_of_det_reflections_per_seq_accumulated = np.zeros(self.number_of_QRAM_sequences)
        self.num_of_det_transmissions_per_seq_accumulated = np.zeros(self.number_of_QRAM_sequences)

        num_of_seq_per_count = 50
        self.num_of_BP_counts_per_n_sequences = np.zeros(self.number_of_QRAM_sequences//num_of_seq_per_count)
        self.num_of_DP_counts_per_n_sequences = np.zeros(self.number_of_QRAM_sequences//num_of_seq_per_count)
        self.num_of_S_counts_per_n_sequences = np.zeros(self.number_of_QRAM_sequences//num_of_seq_per_count)

        self.num_of_det_reflections_per_seq_S = np.zeros(self.number_of_QRAM_sequences)
        self.num_of_det_reflections_per_seq_N = np.zeros(self.number_of_QRAM_sequences)
        self.num_of_det_transmissions_per_seq_S = np.zeros(self.number_of_QRAM_sequences)
        self.num_of_det_transmissions_per_seq_N = np.zeros(self.number_of_QRAM_sequences)


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

        # Get the minimum value between new measure (tt_S_no_gaps) and last old measure
        min_len = min(len(prev_measure), len(curr_measure))

        # Compare values of measurements
        compare_values = np.array(prev_measure[:min_len]) == np.array(curr_measure[:min_len])

        # TODO: replace the line below with the one in comment:
        #new_timetags = np.sum(compare_values) < min_len / 2
        new_timetags = sum(compare_values) < min_len / 2
        return new_timetags

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

            if self.new_timetags_detected(self.batcher['tt_S_measure_batch'], self.tt_S_measure):
                break

            if self.new_timetags_detected(self.batcher['tt_N_measure_batch'], self.tt_N_measure):
                break

            if self.new_timetags_detected(self.batcher['tt_DP_measure_batch'], self.tt_DP_measure):
                break

            if self.new_timetags_detected(self.batcher['tt_BP_measure_batch'], self.tt_BP_measure):
                break

            if self.new_timetags_detected(self.batcher['tt_FS_measure_batch'], self.tt_FS_measure):
                break

            # Are we waiting too long?
            count += 1
            if TOO_MANY_WAITING_CYCLES != -1 and count > TOO_MANY_WAITING_CYCLES:
                self.runs_status = TerminationReason.ERROR
                break

            time.sleep(WAIT_TIME)

        # We return the timestamp - this is when we declare we got the measurement
        timestamp = time.strftime("%Y%m%d-%H%M%S")
        return timestamp

    def experiment_calculations(self):

        self.divide_tt_to_reflection_trans()

        # TODO: replace this with my code?
        self.divide_BP_and_DP_counts(50)

        self.num_of_det_reflections_per_seq = self.num_of_det_reflections_per_seq_S \
                                              + self.num_of_det_reflections_per_seq_N
        self.num_of_det_transmissions_per_seq = self.num_of_det_transmissions_per_seq_S \
                                                + self.num_of_det_transmissions_per_seq_N
        # self.num_of_SPRINT_reflections_per_seq = self.num_of_SPRINT_reflections_per_seq_N \
        #                                          + self.num_of_SPRINT_reflections_per_seq_S
        # self.num_of_SPRINT_transmissions_per_seq = self.num_of_SPRINT_transmissions_per_seq_N \
        #                                            + self.num_of_SPRINT_transmissions_per_seq_S
        self.sum_for_threshold = sum(self.num_of_det_reflections_per_seq[-int(self.reflection_threshold_time // len(
            Config.QRAM_Exp_Gaussian_samples_S)):])  # summing over the reflection from detection pulses of each sequence corresponding the the reflection_threshold_time

        # fold reflections and transmission
        self.fold_tt_histogram(exp_sequence_len=self.QRAM_sequence_len)

        # get the average number of photons in detection pulse
        self.avg_num_of_photons_per_pulse_S_live = self.get_avg_num_of_photons_in_seq_pulses(
            self.folded_tt_S_directional, self.pulses_location_in_seq_S, self.tt_FS_measure)
        self.avg_num_of_photons_per_pulse_N_live = self.get_avg_num_of_photons_in_seq_pulses(
            self.folded_tt_N_directional, self.pulses_location_in_seq_N, self.tt_BP_measure + self.tt_DP_measure)
        self.avg_num_of_photons_per_pulse_A_live = self.get_avg_num_of_photons_in_seq_pulses(
            (np.array(self.folded_tt_S_directional) + np.array(self.folded_tt_BP_timebins)
             + np.array(self.folded_tt_DP_timebins)).tolist(), self.pulses_location_in_seq_A,
            self.tt_BP_measure + self.tt_DP_measure)
        self.avg_num_of_photons_per_pulse_BP_live = self.get_avg_num_of_photons_in_seq_pulses(
            self.folded_tt_BP_timebins, self.pulses_location_in_seq[-2:], self.tt_BP_measure)
        self.avg_num_of_photons_per_pulse_DP_live = self.get_avg_num_of_photons_in_seq_pulses(
            self.folded_tt_DP_timebins, self.pulses_location_in_seq[-2:], self.tt_DP_measure)
        self.avg_num_of_photons_per_pulse_live = self.avg_num_of_photons_per_pulse_S_live + \
                                                 self.avg_num_of_photons_per_pulse_N_live + \
                                                 self.avg_num_of_photons_per_pulse_A_live
        self.avg_num_of_photons_per_pulse_live_MZ = [[x] + [y] for x, y in
                                                     zip(self.avg_num_of_photons_per_pulse_BP_live,
                                                         self.avg_num_of_photons_per_pulse_DP_live)]

        # self.avg_num_of_photons_per_pulse_live_MZ = self.avg_num_of_photons_per_pulse_BP_live + \
        #                                             self.avg_num_of_photons_per_pulse_DP_live

        # get box location on the y-axis:
        self.max_value_per_pulse_S_live = self.get_max_value_in_seq_pulses(self.folded_tt_S_directional,
                                                                           self.pulses_location_in_seq_S)
        self.max_value_per_pulse_N_live = self.get_max_value_in_seq_pulses(self.folded_tt_N_directional,
                                                                           self.pulses_location_in_seq_N)
        self.max_value_per_pulse_A_live = self.get_max_value_in_seq_pulses(
            (np.array(self.folded_tt_S_directional) + np.array(self.folded_tt_BP_timebins) +
             np.array(self.folded_tt_DP_timebins)).tolist(), self.pulses_location_in_seq_A)
        self.max_value_per_pulse_BP_live = self.get_max_value_in_seq_pulses(
            self.folded_tt_BP_timebins, self.pulses_location_in_seq[-2:])
        self.max_value_per_pulse_DP_live = self.get_max_value_in_seq_pulses(
            self.folded_tt_DP_timebins, self.pulses_location_in_seq[-2:])
        self.Num_of_photons_txt_box_y_loc_live = self.max_value_per_pulse_S_live + \
                                                 self.max_value_per_pulse_N_live + \
                                                 self.max_value_per_pulse_A_live
        self.Num_of_photons_txt_box_y_loc_live_MZ = [[x] + [y] for x, y in zip(self.max_value_per_pulse_BP_live,
                                                                               self.max_value_per_pulse_DP_live)]
        # self.Num_of_photons_txt_box_y_loc_live_MZ = self.max_value_per_pulse_BP_live + self.max_value_per_pulse_DP_live

        avg_BP_before = np.average(self.MZ_BP_counts_res_value_1[:self.rep_MZ_check])
        avg_BP_after = np.average(self.MZ_BP_counts_res_value_1[self.rep_MZ_check:])
        avg_DP_before = np.average(self.MZ_DP_counts_res_value_1[:self.rep_MZ_check])
        avg_DP_after = np.average(self.MZ_DP_counts_res_value_1[self.rep_MZ_check:])
        self.Infidelity_before = avg_DP_before / (avg_DP_before + avg_BP_before)
        self.Infidelity_after = avg_DP_after / (avg_DP_after + avg_BP_after)

        pass

    def calculate_running_averages(self):
        self.folded_tt_S_cumulative_avg = Utils.running_average(self.folded_tt_S_cumulative_avg, self.folded_tt_S, self.counter)
        self.folded_tt_N_cumulative_avg = Utils.running_average(self.folded_tt_N_cumulative_avg, self.folded_tt_N, self.counter)
        self.folded_tt_BP_cumulative_avg = Utils.running_average(self.folded_tt_BP_cumulative_avg, self.folded_tt_BP, self.counter)
        self.folded_tt_DP_cumulative_avg = Utils.running_average(self.folded_tt_DP_cumulative_avg, self.folded_tt_DP, self.counter)
        self.folded_tt_FS_cumulative_avg = Utils.running_average(self.folded_tt_FS_cumulative_avg, self.folded_tt_FS, self.counter)
        self.folded_tt_S_directional_cumulative_avg = Utils.running_average(self.folded_tt_S_directional_cumulative_avg, self.folded_tt_S_directional, self.counter)
        self.folded_tt_N_directional_cumulative_avg = Utils.running_average(self.folded_tt_N_directional_cumulative_avg, self.folded_tt_N_directional, self.counter)
        self.folded_tt_BP_timebins_cumulative_avg = Utils.running_average(self.folded_tt_BP_timebins_cumulative_avg, self.folded_tt_BP_timebins, self.counter)
        self.folded_tt_DP_timebins_cumulative_avg = Utils.running_average(self.folded_tt_DP_timebins_cumulative_avg, self.folded_tt_DP_timebins, self.counter)
        pass

    # TODO: Refactor/Rename (this method analyzes the results)
    def Save_SNSPDs_QRAM_Measurement_with_tt(self, N, qram_sequence_len, pre_comment, lock_err_threshold,
                                             transit_condition,
                                             max_probe_counts, filter_delay, reflection_threshold,
                                             reflection_threshold_time,
                                             photons_per_det_pulse_threshold, FLR_threshold, exp_flag, with_atoms,
                                             MZ_infidelity_threshold):
        """
        Function for analyzing, saving and displaying data from SPRINT experiment.
        :param N: Number of maximum experiments (free throws) saved and displayed.
                 program we are looking for transits of atoms next to the toroid and record them.
        :param qram_sequence_len: the number of sprint sequences (detection and sprint pulses combination)
        :param pre_comment: The comment added at the start of the experiment, usually consisting of unique experiment
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

        if not pre_comment:
            pre_comment = self.prompt('Add comment to measurement: ')

        # set constant parameters for the function

        # TODO: Q: Are these just detector "names"/"Symbols"? And also used to define the number of detectors (e.g. 8)?
        self.Num_Of_dets = [1, 2, 3, 4, 5, 6, 7, 8]

        self.delay_in_detection_N = 30  # choose the correct delay in samples to the first detection pulse # TODO: 40?
        self.delay_in_detection_S = 20  # choose the correct delay in samples to the first detection pulse # TODO: 40?

        self.det_pulse_len = Config.det_pulse_len + Config.num_between_zeros

        self.sprint_pulse_len = Config.sprint_pulse_len + Config.num_between_zeros  # TODO: Q: what is num_between_zeros ?

        # TODO: These come as parameters for the experiment run - we may want to have them all under "self.params[..]" so we can
        # TODO: (a) group logically (b) pass to other functions, instead of passing them one by one or relying on them being on self
        self.N = N
        self.QRAM_sequence_len = qram_sequence_len
        self.with_atoms = with_atoms
        self.switch_atom_no_atoms = "atoms" if self.with_atoms else "!atoms"
        self.lock_err_threshold = lock_err_threshold
        self.transit_condition = transit_condition
        self.reflection_threshold = reflection_threshold
        self.reflection_threshold_time = reflection_threshold_time
        self.MZ_infidelity_threshold = MZ_infidelity_threshold
        self.photons_per_det_pulse_threshold = photons_per_det_pulse_threshold
        self.FLR_threshold = FLR_threshold
        self.exp_flag = exp_flag

        self.warm_up_cycles = 3  # Should also come as params, with some default. E.g. self.warm_up_cycles = 10 if 'warm_up_cycles' not in params else params['warm_up_cycles']

        # initialize parameters - set zeros vectors and empty lists
        self.init_params_for_experiment()

        # ----------------------------------------------------------
        # Prepare pulses location of South, North and Ancilla
        # ----------------------------------------------------------

        # TODO: Q: we used to rely on this: "int(Config.num_between_zeros/2)" - why aren't we anymore?
        self.pulses_location_in_seq, self.filter_gen = self.get_pulses_location_in_seq(0,
                                                                                       Config.QRAM_Exp_Gaussian_samples_General,
                                                                                       smearing=5)  # smearing=int(Config.num_between_zeros/2))
        self.pulses_location_in_seq_S, self.filter_S = self.get_pulses_location_in_seq(filter_delay[0],
                                                                                       Config.QRAM_Exp_Gaussian_samples_S,
                                                                                       smearing=5)  # smearing=int(Config.num_between_zeros/2))
        # TODO: Q: Why are we fixing the Gaussian here?
        self.QRAM_Exp_Gaussian_samples_N = Config.QRAM_Exp_Gaussian_samples_N
        self.QRAM_Exp_Gaussian_samples_N[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]] = \
            (np.array(Config.QRAM_Exp_Gaussian_samples_N[
                      self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]]) + \
             np.array(Config.QRAM_Exp_Gaussian_samples_General[
                      self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]]) * \
             Config.sprint_pulse_amp_Early[0]).tolist()
        self.pulses_location_in_seq_N, self.filter_N = self.get_pulses_location_in_seq(filter_delay[1],
                                                                                       self.QRAM_Exp_Gaussian_samples_N,
                                                                                       smearing=5)  # smearing=int(Config.num_between_zeros/2))
        self.pulses_location_in_seq_A, self.filter_A = self.get_pulses_location_in_seq(filter_delay[2],
                                                                                       Config.QRAM_Exp_Gaussian_samples_Ancilla,
                                                                                       smearing=5)  # smearing=int(Config.num_between_zeros/2))

        # Find the center indices of the pulses and concatenate them to one list - to be used to put text boxes in figures
        # TODO: this can be moved/calculated later when we're doing the figure work...
        self.Num_of_photons_txt_box_x_loc = np.concatenate(
            (self.num_of_photons_txt_box_loc(self.pulses_location_in_seq_S),
             self.num_of_photons_txt_box_loc(self.pulses_location_in_seq_N),
             self.num_of_photons_txt_box_loc(self.pulses_location_in_seq_A)))

        self.Num_of_photons_txt_box_x_loc_for_MZ_ports = self.num_of_photons_txt_box_loc(
            self.pulses_location_in_seq[-2:])
        self.num_of_detection_pulses = len(Config.det_pulse_amp_N)

        # Take the 2nd number in the last tuple (which is the time of the last pulse)
        self.end_of_det_pulse_in_seq = self.pulses_location_in_seq[self.num_of_detection_pulses - 1][1]

        # Start experiment flag and set MOT according to flag
        self.updateValue("Experiment_Switch", True)
        self.MOT_switch(with_atoms=self.with_atoms, update_parameters=True)

        self.logger.blue('Press ESC to stop measurement.')

        # Initialize the batcher
        self.batcher.set_batch_size(N)
        self.batcher.empty_all()

        # Associate the streams filled in OPX (FPGA) code with result handles
        self.get_handles_from_OPX_server()

        ############################# WHILE 1 - START #############################

        # TODO: remove the "old" warmup loop
        WARMUP_CYCLES = -1
        cycle = 0
        self.sum_for_threshold = self.reflection_threshold
        while exp_flag and (cycle < WARMUP_CYCLES or self.sum_for_threshold > self.reflection_threshold):

            if self.keyPress == 'ESC':
                self.logger.blue('ESC pressed. Stopping measurement.')
                self.updateValue("Experiment_Switch", False)
                self.MOT_switch(True)
                self.update_parameters()
                break

            # -------------------------------------------
            # Deal with the time-tags and counts
            # -------------------------------------------

            # Filter/Manipulate the values we got
            self.get_results_from_streams()
            self.ingest_time_tags()

            # TODO: Review:
            # TODO: The below calls come instead of:          self.divide_tt_to_reflection_trans(...)

            buckets = Utils.bucket_timetags(
                timetags=self.tt_N_measure + self.tt_BP_measure + self.tt_DP_measure,
                window_size=self.experiment["sequence_length"],  # [ns]
                buckets_number=self.experiment["number_of_sequences"],
                start_time=int(0.6e6),
                end_time=int(9.6e6),
                filters=[{"name": "reflections_south", "filter": self.filter_S},
                         {"name": "transmissions_north", "filter": self.filter_N}])

            self.num_of_det_reflections_per_seq_S = buckets["reflections_south"]["counts"]
            self.num_of_det_transmissions_per_seq_N = buckets["transmissions_north"]["counts"]

            buckets = Utils.bucket_timetags(
                timetags=self.tt_S_measure + self.tt_FS_measure,
                window_size=self.experiment["sequence_length"],  # [ns]
                buckets_number=self.experiment["number_of_sequences"],
                start_time=int(0.6e6),
                end_time=int(9.6e6),
                filters=[{"name": "reflections_north", "filter": self.filter_N},
                         {"name": "transmissions_south", "filter": self.filter_S}])

            self.num_of_det_reflections_per_seq_N = buckets["reflections_north"]["counts"]
            self.num_of_det_transmissions_per_seq_S = buckets["transmissions_south"]["counts"]

            # self.divide_tt_to_reflection_trans(self.sprint_pulse_len, self.num_of_detection_pulses)

            # TODO: Review:
            # TODO: The below calls come instead of:          self.divide_BP_and_DP_counts(50)
            buckets = Utils.bucket_timetags(
                timetags=self.tt_BP_measure,
                window_size=self.experiment["sequence_length"],  # [ns]
                buckets_number=self.experiment["number_of_sequences"] // 50,
                start_time=self.end_of_det_pulse_in_seq,
                filters=[{"name": "BP_counts_per_n_sequences", "filter": np.ceil(Config.QRAM_Exp_Square_samples_Late)}]
            )
            self.num_of_BP_counts_per_n_sequences = buckets["BP_counts_per_n_sequences"]["counts"]

            buckets = Utils.bucket_timetags(
                timetags=self.tt_DP_measure,
                window_size=self.experiment["sequence_length"],  # [ns]
                buckets_number=self.experiment["number_of_sequences"] // 50,
                start_time=self.end_of_det_pulse_in_seq,
                filters=[{"name": "DP_counts_per_n_sequences", "filter": np.ceil(Config.QRAM_Exp_Square_samples_Late)}]
            )
            self.num_of_DP_counts_per_n_sequences = buckets["DP_counts_per_n_sequences"]["counts"]

            buckets = Utils.bucket_timetags(
                timetags=self.tt_FS_measure + self.tt_S_measure,
                window_size=self.experiment["sequence_length"],  # [ns]
                buckets_number=self.experiment["number_of_sequences"] // 50,
                start_time=self.end_of_det_pulse_in_seq,
                filters=[{"name": "S_counts_per_n_sequences", "filter": np.ceil(Config.QRAM_Exp_Square_samples_Late)}]
            )
            self.num_of_S_counts_per_n_sequences = buckets["S_counts_per_n_sequences"]["counts"]

            # self.divide_BP_and_DP_counts(50)

            self.num_of_det_reflections_per_seq = self.num_of_det_reflections_per_seq_S \
                                                  + self.num_of_det_reflections_per_seq_N
            self.num_of_det_transmissions_per_seq = self.num_of_det_transmissions_per_seq_S \
                                                    + self.num_of_det_transmissions_per_seq_N
            # self.num_of_SPRINT_reflections_per_seq = self.num_of_SPRINT_reflections_per_seq_S \
            #                                          + self.num_of_SPRINT_reflections_per_seq_N
            # self.num_of_SPRINT_transmissions_per_seq = self.num_of_SPRINT_transmissions_per_seq_S \
            #                                            + self.num_of_SPRINT_transmissions_per_seq_N

            # Summing over the reflection from detection pulses of each sequence corresponding to the reflection_threshold_time
            self.sum_for_threshold = sum(self.num_of_det_reflections_per_seq[
                                         -int(self.reflection_threshold_time // len(Config.QRAM_Exp_Gaussian_samples_S)):])

            # Check locking error, break if we are above threshold
            self.lock_err = (lock_err_threshold / 2) if exp_flag else self._read_locking_error()
            self.logger.debug(f'{self.lock_err}, {self.lock_err > lock_err_threshold}, {self.sum_for_threshold}')
            if self.lock_err > lock_err_threshold:
                break

            cycle += 1

        ############################# WHILE 1 - END #############################

        ####    end get tt and counts from OPX to python   #####

        self.acquisition_flag = True
        self.threshold_flag = True
        self.pause_flag = False
        self.Phase_Correction_min_diff = []

        # Prepare the sub figures needed for the plotting phase
        self.prepare_figures()

        ############################################ START WHILE LOOP #################################################

        # Initialize the variables before loop starts
        self.counter = 1  # Total number of successful cycles
        self.repetitions = 1  # Total number of cycles
        self.acquisition_flag = True
        self.warm_up = True
        # TODO: Change this to Experiment Phase - with enum
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

            # Informational printing
            self.print_experiment_information()

            # If we are still in warm-up phase, we do not continue further and go back to beginning of loop
            # self.warm_up = self.is_warm_up_phase()
            # if self.warm_up:
            #     continue

            # Perform all analytics and calculations needed for display
            self.experiment_calculations()

            # Plot figures
            self.plot_figures()

            # Determine if we're "acquired" - e.g., worthy to save data :-)
            self.acquisition_flag = self.is_acquired()
            if self.acquisition_flag and not self.pause_flag:

                if self.counter < N:
                    self.counter += 1

                self.logger.info('Sum of reflections: %d' % self.sum_for_threshold)
                self.num_of_det_reflections_per_seq_accumulated += self.num_of_det_reflections_per_seq_S \
                                                                   + self.num_of_det_reflections_per_seq_N
                self.num_of_det_transmissions_per_seq_accumulated += self.num_of_det_transmissions_per_seq_S \
                                                                     + self.num_of_det_transmissions_per_seq_N

                self.seq_transit_events_live = np.zeros(self.number_of_QRAM_sequences)

                ### Find transits and build histogram:  ###
                self.find_transits_and_sprint_events_changed(cond=self.transit_condition, minimum_number_of_seq_detected=2)
                self.seq_transit_events_live[[vec for elem in self.all_transits_seq_indx for vec in elem]] += 1
                self.seq_transit_events_batched[[vec for elem in self.all_transits_seq_indx for vec in elem]] += 1

                self.number_of_transits_live = len(self.all_transits_seq_indx)
                self.number_of_transits_total = len([vec for lst in self.batcher['all_transits_seq_indx_batch'] for vec in lst])

                self.calculate_running_averages()

                # get the average number of photons in detection pulse
                self.avg_num_of_photons_per_pulse_S = self.get_avg_num_of_photons_in_seq_pulses(
                    self.folded_tt_S_directional_cumulative_avg, self.pulses_location_in_seq_S, [])
                self.avg_num_of_photons_per_pulse_N = self.get_avg_num_of_photons_in_seq_pulses(
                    self.folded_tt_N_directional_cumulative_avg, self.pulses_location_in_seq_N, [])
                self.avg_num_of_photons_per_pulse_A = self.get_avg_num_of_photons_in_seq_pulses(
                    (np.array(self.folded_tt_S_directional_cumulative_avg) + np.array(self.folded_tt_BP_timebins_cumulative_avg)
                     + np.array(self.folded_tt_DP_timebins_cumulative_avg)).tolist(), self.pulses_location_in_seq_A, [])
                self.avg_num_of_photons_per_pulse_BP = self.get_avg_num_of_photons_in_seq_pulses(
                    self.folded_tt_BP_timebins_cumulative_avg, self.pulses_location_in_seq[-2:], [])
                self.avg_num_of_photons_per_pulse_DP = self.get_avg_num_of_photons_in_seq_pulses(
                    self.folded_tt_DP_timebins_cumulative_avg, self.pulses_location_in_seq[-2:], [])
                self.avg_num_of_photons_per_pulse = self.avg_num_of_photons_per_pulse_S + \
                                                    self.avg_num_of_photons_per_pulse_N + \
                                                    self.avg_num_of_photons_per_pulse_A
                self.avg_num_of_photons_per_pulse_MZ = [[x] + [y] for x, y in
                                                        zip(self.avg_num_of_photons_per_pulse_BP,
                                                            self.avg_num_of_photons_per_pulse_DP)]
                # self.avg_num_of_photons_per_pulse_MZ = self.avg_num_of_photons_per_pulse_BP + \
                #                                        self.avg_num_of_photons_per_pulse_DP
                # get box location on the y-axis:
                self.max_value_per_pulse_S = self.get_max_value_in_seq_pulses(self.folded_tt_S_directional_cumulative_avg,
                                                                              self.pulses_location_in_seq_S)
                self.max_value_per_pulse_N = self.get_max_value_in_seq_pulses(self.folded_tt_N_directional_cumulative_avg,
                                                                              self.pulses_location_in_seq_N)
                self.max_value_per_pulse_A = self.get_max_value_in_seq_pulses(
                    (np.array(self.folded_tt_S_directional_cumulative_avg) + np.array(self.folded_tt_BP_timebins_cumulative_avg) +
                     np.array(self.folded_tt_DP_timebins_cumulative_avg)).tolist(), self.pulses_location_in_seq_A)
                self.max_value_per_pulse_BP = self.get_max_value_in_seq_pulses(
                    self.folded_tt_BP_timebins_cumulative_avg, self.pulses_location_in_seq[-2:])
                self.max_value_per_pulse_DP = self.get_max_value_in_seq_pulses(
                    self.folded_tt_DP_timebins_cumulative_avg, self.pulses_location_in_seq[-2:])
                self.Num_of_photons_txt_box_y_loc = self.max_value_per_pulse_S + self.max_value_per_pulse_N + \
                                                    self.max_value_per_pulse_A
                self.Num_of_photons_txt_box_y_loc_MZ = [[x] + [y] for x, y in zip(self.max_value_per_pulse_BP,
                                                                                  self.max_value_per_pulse_DP)]
                # self.Num_of_photons_txt_box_y_loc_MZ = self.max_value_per_pulse_BP + self.max_value_per_pulse_DP

                # Batch all data samples
                self.batcher.batch_all(self)

            # Did user request we change with/without atoms state?
            self.handle_user_atoms_on_off_switch()

            # Did we complete N successful iterations? (we check for N+1 because we started with 1)
            if self.counter == N+1:
                break

            # We completed another repetition.
            self.repetitions += 1

        ############################################## END WHILE-LOOP #################################################

        if self.runs_status == TerminationReason.SUCCESS:
            self.logger.info(f'Finished {N} Runs, {"with" if self.with_atoms else "without"} atoms')
        elif self.runs_status == TerminationReason.USER:
            self.logger.info(f'User Terminated the measurements after {self.counter} Runs, {"with" if self.with_atoms else "without"} atoms')
        elif self.runs_status == TerminationReason.ERROR:
            self.logger.info(f'Error Terminated the measurements after {self.counter} Runs, {"with" if self.with_atoms else "without"} atoms')

        # Adding comment to measurement [prompt whether stopped or finished regularly]
        if exp_flag:
            if self.counter < N:
                aftComment = pymsgbox.prompt('Add comment to measurement: ', default='', timeout=int(30e3))
            else:
                aftComment = ''
        else:
            aftComment = 'ignore'

        experiment_comment = f'Transit condition: {self.transit_condition}\nReflection threshold {self.reflection_threshold} @ {int(self.reflection_threshold_time/1e6)} ms'
        daily_experiment_comments = self.generate_experiment_summary_line(pre_comment, aftComment, self.with_atoms, self.counter)

        # Save all results of experiment
        self.save_experiment_results(experiment_comment, daily_experiment_comments)

        # End of Experiment! Ensure we turn the experiment flag off and return the MOT
        self.updateValue("Experiment_Switch", False)
        self.MOT_switch(with_atoms=True, update_parameters=True)

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
            "early_sequence": Config.QRAM_Exp_Gaussian_samples_Early,
            "late_sequence": Config.QRAM_Exp_Gaussian_samples_Late,
            "north_sequence": Config.QRAM_Exp_Gaussian_samples_N,
            "south_sequence": Config.QRAM_Exp_Gaussian_samples_S,
            "fs_sequence": Config.QRAM_Exp_Square_samples_FS,

            "tt_measure_batch": self.batcher['tt_measure_batch'],
            "tt_N_measure_batch": self.batcher['tt_N_measure_batch'],
            "tt_S_measure_batch": self.batcher['tt_S_measure_batch'],
            "tt_FS_measure_batch": self.batcher['tt_FS_measure_batch'],
            "tt_BP_measure_batch": self.batcher['tt_BP_measure_batch'],
            "tt_DP_measure_batch": self.batcher['tt_DP_measure_batch'],

            "MZ_BP_counts_balancing_batch": self.batcher['MZ_BP_counts_balancing_batch'],
            "MZ_BP_counts_balancing_check_batch": self.batcher['MZ_BP_counts_balancing_check_batch'],
            "MZ_DP_counts_balancing_batch": self.batcher['MZ_DP_counts_balancing_batch'],
            "MZ_DP_counts_balancing_check_batch": self.batcher['MZ_DP_counts_balancing_check_batch'],
            "Phase_Correction_vec_batch": self.batcher['Phase_Correction_vec_batch'],
            "Phase_Correction_min_vec_batch": self.batcher['Phase_Correction_min_vec_batch'],

            "FLR_measurement": self.batcher['flr_batch'],
            "lock_error": self.batcher['lock_err_batch'],
            "exp_timestr": experiment_comment,

            "exp_comment": f'transit condition: {self.transit_condition}; reflection threshold: {self.reflection_threshold} @ {int(self.reflection_threshold_time / 1e6)} ms',
            "daily_experiment_comments": daily_experiment_comments,

            "experiment_config_values": self.Exp_Values
        }
        self.bd_results.save_results(results)

    def is_acquired(self):

        # If we are not in experiment, we don't care about the rest - just say we're acquired.
        if not self.exp_flag:
            return True

        # Are we properly locked on resonance?
        if self.lock_err > self.lock_err_threshold:
            return False

        # Is fluorescence strong enough? (assuming we're running with atoms)
        if self.fluorescence_average < self.FLR_threshold and self.with_atoms:
            return False

        # TODO: should be self.avg... What does the below do?
        if len(experiment.avg_num_of_photons_per_pulse_live) == 0 or np.average(experiment.avg_num_of_photons_per_pulse_live) > self.photons_per_det_pulse_threshold:
            return False

        # Are any of the detectors latched?
        if self.latched_detectors():
           return False

        threshold_flag = (self.sum_for_threshold < self.reflection_threshold) and \
                              (self.Infidelity_before <= self.MZ_infidelity_threshold) and \
                              (self.Infidelity_after <= self.MZ_infidelity_threshold)
        if threshold_flag or not self.exp_flag:
            return True

        # All is well!
        return True

    def handle_warm_up_phase(self):

        if self.post_warm_up_completed:
            return False

        # Get data from OPX streams
        self.get_results_from_streams()

        # Make something out of the data we received on the streams
        self.ingest_time_tags()

        self.divide_tt_to_reflection_trans()
        self.divide_BP_and_DP_counts(50)
        self.num_of_det_reflections_per_seq = self.num_of_det_reflections_per_seq_S + self.num_of_det_reflections_per_seq_N
        self.num_of_det_transmissions_per_seq = self.num_of_det_transmissions_per_seq_S + self.num_of_det_transmissions_per_seq_N

        # Summing over the reflection from detection pulses of each sequence corresponding the the reflection_threshold_time
        self.sum_for_threshold = sum(self.num_of_det_reflections_per_seq[-int(self.reflection_threshold_time // len(Config.QRAM_Exp_Gaussian_samples_S)):])

        # Check if conditions have been met for the completeion of the warm-up phase
        warm_up_phase_complete = self.is_warm_up_phase_complete()

        # This is the first time we realize we are no loger in warm-up - so run post stage
        if warm_up_phase_complete:
            self.post_warm_up()
            self.post_warm_up_completed = True  # Mark the fact we are done

        return not warm_up_phase_complete

    def is_warm_up_phase_complete(self):
        """
        We are in warm-up if these conditions are met:
        - We didn't complete already the warm-up post process
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
        if self.exp_flag and self.sum_for_threshold > self.reflection_threshold:
            return False

        return True

    def post_warm_up(self):
        """
        Assuming warm-up phase has completed and we have some information, we now prepare for the display and the experiment cycles
        """

        self.num_of_det_reflections_per_seq_accumulated += self.num_of_det_reflections_per_seq_S + self.num_of_det_reflections_per_seq_N
        self.num_of_det_transmissions_per_seq_accumulated += self.num_of_det_transmissions_per_seq_S + self.num_of_det_transmissions_per_seq_N

        # Divide south and north into reflection and transmission
        self.tt_histogram_transmission, self.tt_histogram_reflection = \
        self.divide_to_reflection_trans(sprint_pulse_len=self.sprint_pulse_len,
                                        num_of_det_pulses=len(Config.det_pulse_amp_S),
                                        num_of_sprint_pulses=len(Config.sprint_pulse_amp_S),
                                        num_of_sprint_sequences=self.number_of_QRAM_sequences)

        # TODO: needed?
        self.tt_S_SPRINT_events = np.zeros(self.QRAM_sequence_len)
        self.tt_S_SPRINT_events_batch = np.zeros(self.QRAM_sequence_len)
        self.tt_Single_det_SPRINT_events = np.zeros((len(self.Num_Of_dets), self.QRAM_sequence_len))
        self.tt_Single_det_SPRINT_events_batch = np.zeros((len(self.Num_Of_dets), self.QRAM_sequence_len))
        self.tt_N_det_SPRINT_events_batch = np.zeros(self.QRAM_sequence_len)
        self.tt_S_det_SPRINT_events_batch = np.zeros(self.QRAM_sequence_len)

        # Fold time-tags: S, N, BP, DP, FS, N-directional, S-directional, BP-timebins, DP-timebins
        self.fold_tt_histogram(exp_sequence_len=self.QRAM_sequence_len)

        # Batch folded tt "N", "S", BP, DP, FS
        self.folded_tt_S_cumulative_avg = self.folded_tt_S
        self.folded_tt_N_cumulative_avg = self.folded_tt_N
        self.folded_tt_BP_cumulative_avg = self.folded_tt_BP
        self.folded_tt_DP_cumulative_avg = self.folded_tt_DP
        self.folded_tt_FS_cumulative_avg = self.folded_tt_FS
        self.folded_tt_S_directional_cumulative_avg = self.folded_tt_S_directional
        self.folded_tt_N_directional_cumulative_avg = self.folded_tt_N_directional
        self.folded_tt_BP_timebins_cumulative_avg = self.folded_tt_BP_timebins
        self.folded_tt_DP_timebins_cumulative_avg = self.folded_tt_DP_timebins

        # Get the average number of photons in detection pulse
        self.avg_num_of_photons_per_pulse_S_live = self.get_avg_num_of_photons_in_seq_pulses(self.folded_tt_S_directional, self.pulses_location_in_seq_S, self.tt_FS_measure)
        self.avg_num_of_photons_per_pulse_N_live = self.get_avg_num_of_photons_in_seq_pulses(self.folded_tt_N_directional, self.pulses_location_in_seq_N, self.tt_BP_measure + self.tt_DP_measure)
        self.avg_num_of_photons_per_pulse_A_live = self.get_avg_num_of_photons_in_seq_pulses((np.array(self.folded_tt_S_directional) + np.array(self.folded_tt_BP_timebins)+ np.array(self.folded_tt_DP_timebins)).tolist(), self.pulses_location_in_seq_A,self.tt_BP_measure + self.tt_DP_measure)
        self.avg_num_of_photons_per_pulse_BP_live = self.get_avg_num_of_photons_in_seq_pulses(self.folded_tt_BP_timebins, self.pulses_location_in_seq[-2:], self.tt_BP_measure)
        self.avg_num_of_photons_per_pulse_DP_live = self.get_avg_num_of_photons_in_seq_pulses(self.folded_tt_DP_timebins, self.pulses_location_in_seq[-2:], self.tt_DP_measure)
        self.avg_num_of_photons_per_pulse_live = self.avg_num_of_photons_per_pulse_S_live + self.avg_num_of_photons_per_pulse_N_live + self.avg_num_of_photons_per_pulse_A_live
        self.avg_num_of_photons_per_pulse_live_MZ = [[x]+[y] for x, y in zip(self.avg_num_of_photons_per_pulse_BP_live, self.avg_num_of_photons_per_pulse_DP_live)]

        self.avg_num_of_photons_per_pulse = self.avg_num_of_photons_per_pulse_live
        self.avg_num_of_photons_per_pulse_MZ = self.avg_num_of_photons_per_pulse_live_MZ

        # get box location on the y-axis:
        self.max_value_per_pulse_S_live = self.get_max_value_in_seq_pulses(self.folded_tt_S_directional, self.pulses_location_in_seq_S)
        self.max_value_per_pulse_N_live = self.get_max_value_in_seq_pulses(self.folded_tt_N_directional, self.pulses_location_in_seq_N)
        self.max_value_per_pulse_A_live = self.get_max_value_in_seq_pulses((np.array(self.folded_tt_S_directional) + np.array(self.folded_tt_BP_timebins) + np.array(self.folded_tt_DP_timebins)).tolist(), self.pulses_location_in_seq_A)
        self.max_value_per_pulse_BP_live = self.get_max_value_in_seq_pulses(self.folded_tt_BP_timebins, self.pulses_location_in_seq[-2:])
        self.max_value_per_pulse_DP_live = self.get_max_value_in_seq_pulses(self.folded_tt_DP_timebins, self.pulses_location_in_seq[-2:])

        self.Num_of_photons_txt_box_y_loc_live = self.max_value_per_pulse_S_live + self.max_value_per_pulse_N_live + self.max_value_per_pulse_A_live
        self.Num_of_photons_txt_box_y_loc_live_MZ = [[x]+[y] for x, y in zip(self.max_value_per_pulse_BP_live, self.max_value_per_pulse_DP_live)]
        # self.Num_of_photons_txt_box_y_loc_live_MZ = self.max_value_per_pulse_BP_live + self.max_value_per_pulse_DP_live

        self.Num_of_photons_txt_box_y_loc = self.Num_of_photons_txt_box_y_loc_live
        self.Num_of_photons_txt_box_y_loc_MZ = self.Num_of_photons_txt_box_y_loc_live_MZ

        avg_BP_before = np.average(self.MZ_BP_counts_res_value_1[:self.rep_MZ_check])
        avg_BP_after = np.average(self.MZ_BP_counts_res_value_1[self.rep_MZ_check:])
        avg_DP_before = np.average(self.MZ_DP_counts_res_value_1[:self.rep_MZ_check])
        avg_DP_after = np.average(self.MZ_DP_counts_res_value_1[self.rep_MZ_check:])

        self.Infidelity_before = avg_DP_before / (avg_DP_before + avg_BP_before)
        self.Infidelity_after = avg_DP_after / (avg_DP_after + avg_BP_after)

        self.find_transits_and_sprint_events_changed(cond=self.transit_condition, minimum_number_of_seq_detected=2)

        self.seq_transit_events_live[[elem for vec in self.all_transits_seq_indx for elem in vec]] += 1
        self.seq_transit_events_batched[[elem for vec in self.all_transits_seq_indx for elem in vec]] += 1

        self.number_of_transits_live = len(self.all_transits_seq_indx)

        self.batcher.batch_all(self)

        self.number_of_transits_total = len([vec for lst in self.batcher['all_transits_seq_indx_batch'] for vec in lst])
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
    def generate_experiment_summary_line(self, pre_comment, aftComment, with_atoms, counter):
        time_str = time.strftime("%H%M%S")
        date_str = time.strftime("%Y%m%d")
        cmnt = None
        if pre_comment is not None:
            cmnt = pre_comment + '; '
        if aftComment is not None:
            cmnt = pre_comment + ' After comment: ' + aftComment
        if pre_comment is None and aftComment is None:
            cmnt = 'No comment.'
        experiment_success = 'ignore' if 'ignore' in cmnt else 'valid'
        full_line = f'{date_str},{time_str},{experiment_success},{with_atoms},{counter},{cmnt}'
        return full_line

    def pre_run(self, run_parameters):
        # Change the pre comment based on the with_atoms parameter
        suffix = ' with atoms' if run_parameters['with_atoms'] else ' without atoms'
        run_parameters['pre_comment'] += suffix
        pass

    def run(self, run_parameters):
        rp = run_parameters  # Set so we can use in short - "rp", instead of "run_parameters"...

        # max_probe_counts = self.Get_Max_Probe_counts(3)  # return the average maximum probe counts of 3 cycles.
        max_probe_counts = None  # return the average maximum probe counts of 3 cycles.

        # Set switches
        # TODO: Can we move here to Enums?
        self.Experiment_Switch(True)
        self.MOT_switch(rp['with_atoms'])
        self.update_parameters()

        # TODO: Q: Config.QRAM_Exp_Gaussian_samples_S is constructed in a function, using the parameter "sprint_pulse_len" - so why not use it here?
        # TODO: Q: (a) we don't want to use duplicate variables holding the same value, (b) it mentions "samples_S" - but it's the same for "N" as well...
        run_status = self.Save_SNSPDs_QRAM_Measurement_with_tt(N=rp['N'],
                                                  #qram_sequence_len=len(Config.QRAM_Exp_Gaussian_samples_S),  # TODO: Which one of these should it be?
                                                  qram_sequence_len=len(Config.QRAM_Exp_Square_samples_Late),
                                                  pre_comment=rp['pre_comment'],
                                                  lock_err_threshold=rp['lock_err_threshold'],
                                                  transit_condition=rp['transit_condition'],
                                                  max_probe_counts=max_probe_counts,
                                                  filter_delay=rp['filter_delay'],
                                                  reflection_threshold=rp['reflection_threshold'],
                                                  reflection_threshold_time=rp['reflection_threshold_time'],
                                                  photons_per_det_pulse_threshold=rp['photons_per_det_pulse_threshold'],
                                                  FLR_threshold=rp['FLR_threshold'],
                                                  exp_flag=rp['Exp_flag'],
                                                  with_atoms=rp['with_atoms'],
                                                  MZ_infidelity_threshold=rp['MZ_infidelity_threshold'])
        return run_status

    def post_run(self, run_parameters):
        pass


if __name__ == "__main__":

    matplotlib_version = matplotlib.get_backend()
    print(f"In use: {matplotlib_version}")
    matplotlib.use("Qt5Agg")
    #plt.ion()

    run_parameters = {
        'N': 370,  # 50,
        'transit_condition': [2, 1, 2],
        'pre_comment': '"N-S-N-S-N-S-N-S experiment"',
        'lock_err_threshold': 0.005,
        'filter_delay': [0, 0, 0],
        'reflection_threshold': 2550,
        'reflection_threshold_time': 9e6,
        'FLR_threshold': 0.08,
        'MZ_infidelity_threshold': 1.12,
        'photons_per_det_pulse_threshold': 12,
        'Exp_flag': True,
        'with_atoms': True
    }
    sequence_definitions = {
        'total_cycles': 2,
        'delay_between_cycles': None,  # seconds
        'sequence': [
            {
                'parameters': {
                    'N': 500,
                    'with_atoms': True
                }
            },
            {
                'parameters': {
                    'N': 50,
                    'with_atoms': False
                }
            }
        ]
    }

    experiment = QRAMExperiment()

    # TODO: REMOVE, for debug only
    sequence_definitions = None

    if sequence_definitions is None:
        experiment.run(run_parameters)
    else:
        experiment.run_sequence(sequence_definitions, run_parameters)
