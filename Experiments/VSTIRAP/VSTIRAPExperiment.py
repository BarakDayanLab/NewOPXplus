from Experiments.BaseExperiment.BaseExperiment import BaseExperiment
from Experiments.BaseExperiment.BaseExperiment import TerminationReason
from Experiments.VSTIRAP import VSTIRAP_Config_Experiment as Config

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import time
import math
import traceback
from importlib.machinery import SourceFileLoader

from Utilities.BDSound import SOUNDS
from Utilities.BDDialog import BDDialog
from Utilities.Utils import Utils


class VSTIRAPExperiment(BaseExperiment):
    def __init__(self, playback_parameters=None, save_raw_data=False):

        # ------------------------------------------------------------------------------------------------------------
        # If we are in playback mode, we load the VSTIRAP_Config_Experiment.py file that was running during the recording
        #
        # Note:
        # - We do this BEFORE calling super().__init__ as it may invoke methods in this file - and we want the
        #   playback Config to be used!
        # - We use "global" to override the Config that was used in the import statement
        # ------------------------------------------------------------------------------------------------------------
        global Config
        if playback_parameters['active']:
            the_path = os.path.join(playback_parameters['playback_files_path'], 'Source Files', 'VSTIRAP_Config_Experiment.py')
            Config = SourceFileLoader("VSTIRAP_Config_Experiment", the_path).load_module()

        # Invoking BaseClass constructor. It will initiate OPX, QuadRF, BDLogger, Camera, BDResults, KeyEvents etc.
        super().__init__(playback_parameters, save_raw_data)

        pass

    def __del__(self):
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
        self.ignored_marginals = self.Exp_Values['ignored_marginals']  # [nsec] data at the beginning and at the end of a window to "throw away" because of shutters noise
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
                               * 1e6 / len(Config.PNSA_MZ_balance_pulse_Late))
        self.rep_MZ_fast_scan = int(
            0.3 * (self.prepulse_duration - self.Shutter_open_time - self.Balancing_check_window)
            * 1e6 / len(Config.PNSA_MZ_balance_pulse_Late))
        self.rep_MZ_slow_scan = int(
            0.4 * (self.prepulse_duration - self.Shutter_open_time - self.Balancing_check_window)
            * 1e6 / len(Config.PNSA_MZ_balance_pulse_Late))
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
        self.rep_MZ_check = int(self.Balancing_check_window * 1e6 / len(Config.PNSA_MZ_balance_pulse_North))

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

    def ingest_time_tags(self):
        """
        Takes all the raw results we got from the streams and does some processing on them - preparing "measures"
        """

        # Get fluorescence data
        self.FLR_res = -self.streams['FLR_measure']['results']
        self.fluorescence_average = 100000 * self.FLR_res

        # Clear time-tags vectors from last data
        self.tt_measure = []

        # Normalize the data coming from the detectors
        for detector_index in range(len(self.Num_Of_dets)):  # for different detectors
            tt_res = self.streams[f'Detector_{detector_index + 1}_Timetags']['results']
            if isinstance(tt_res, np.ndarray):
                tt_res = tt_res.astype(np.int32) # changing type from int64 to int32
            normalized_stream = self.bdstreams.normalize_stream(tt_res, detector_index, Config.detector_delays,
                                                                self.M_window, self.ignored_marginals)
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
        # Detector 5
        self.tt_S_measure = sorted(sum(self.tt_measure[4:5], []))
        # Unify detectors 6 & 7
        self.tt_FS_measure = sorted(sum(self.tt_measure[5:7], []))
        # tt's from all North detectors
        self.tt_N_measure_Total = self.tt_N_measure + self.tt_DP_measure + self.tt_BP_measure
        # tt's from all South detectors
        self.tt_S_measure_Total = self.tt_S_measure + self.tt_FS_measure

        # --------------------------

        # Calculate the total North/South counts in cycle
        self.total_south_clicks = len(self.tt_S_measure) + len(self.tt_FS_measure)
        self.total_north_clicks = len(self.tt_BP_measure) + len(self.tt_DP_measure) + len(self.tt_N_measure)
        self.total_clicks = self.total_north_clicks + self.total_south_clicks

        self.avg_clicks += (self.total_clicks - self.avg_clicks) / (self.counter+5)
        self.avg_north_clicks += (-self.avg_north_clicks + self.total_north_clicks) / (self.counter+5)
        self.avg_south_clicks += (-self.avg_south_clicks + self.total_south_clicks) / (self.counter+5)

        self.logger.info("total clicks in cycle is %d , average %d " % (self.total_clicks,self.avg_clicks ))
        self.logger.info("total clicks in North cycle is %d, average %d" % (self.total_north_clicks,self.avg_north_clicks))
        self.logger.info("total clicks in South cycle is %d, average %d" % (self.total_south_clicks,self.avg_south_clicks))

        # --------------------------

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

    def latched_detectors(self):
        latched_detectors = []
        for indx, det_tt_vec in enumerate(self.tt_measure):  # for different detectors
            if not det_tt_vec:
                latched_detectors.append(indx)
        return latched_detectors

    def saturated_detectors(self):
        saturated_detectors = []
        for indx, det_tt_vec in enumerate(self.tt_measure):  # for different detectors
            if len(det_tt_vec) >= (Config.vec_size * 0.99):
                saturated_detectors.append(indx)
        return saturated_detectors

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
        self.num_of_det_reflections_per_seq_S = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_det_reflections_per_seq_N = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_det_transmissions_per_seq_S = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_det_transmissions_per_seq_N = np.zeros(self.number_of_PNSA_sequences)

        self.num_of_SPRINT_reflections_per_seq_S = np.zeros(
            [self.number_of_PNSA_sequences, self.number_of_SPRINT_pulses_per_seq])
        self.num_of_SPRINT_reflections_per_seq_N = np.zeros(
            [self.number_of_PNSA_sequences, self.number_of_SPRINT_pulses_per_seq])
        self.num_of_SPRINT_transmissions_per_seq_S = np.zeros(
            [self.number_of_PNSA_sequences, self.number_of_SPRINT_pulses_per_seq])
        self.num_of_SPRINT_transmissions_per_seq_N = np.zeros(
            [self.number_of_PNSA_sequences, self.number_of_SPRINT_pulses_per_seq])

        # tt_small_perturb = []
        for element in self.tt_N_measure + self.tt_BP_measure + self.tt_DP_measure:
            # TODO: Q: I assume this is a time-window to ignore some time on start/end, how did we decide on this time?  TODO: Q: what us the meaning of this specific time-window?
            if (element > int(0.6e6)) and (element < int(self.M_time - 0.4e6)):
                tt_inseq = element % self.sequence_len
                if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                    seq_num = (element - 1) // self.sequence_len
                    self.num_of_det_reflections_per_seq_S[seq_num] += self.filter_S[tt_inseq]
                    self.num_of_det_transmissions_per_seq_N[seq_num] += self.filter_N[tt_inseq]
            # else:  # The part of the SPRINT pulses in the sequence
            #     SPRINT_pulse_num = (tt_inseq - self.end_of_det_pulse_in_seq) // (self.sprint_pulse_len + Config.num_between_zeros)
            #     if SPRINT_pulse_num < self.number_of_SPRINT_pulses_per_seq:
            #         self.num_of_SPRINT_reflections_per_seq_S[(element-1) // self.sequence_len][SPRINT_pulse_num] += 1
            #         self.num_of_SPRINT_transmissions_per_seq_N[(element-1) // self.sequence_len][SPRINT_pulse_num] += 1
            #         tt_small_perturb+=[element]

        for element in self.tt_S_measure + self.tt_FS_measure:
            if (element > int(0.6e6)) and (element < int(self.M_time - 0.4e6)):
                tt_inseq = element % self.sequence_len
                if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                    seq_num = (element - 1) // self.sequence_len
                    self.num_of_det_reflections_per_seq_N[seq_num] += self.filter_N[tt_inseq]
                    self.num_of_det_transmissions_per_seq_S[seq_num] += self.filter_S[tt_inseq]
            # else:  # The part of the SPRINT pulses in the sequence
            #     SPRINT_pulse_num = (tt_inseq - self.end_of_det_pulse_in_seq) // (self.sprint_pulse_len + Config.num_between_zeros)
            #     if SPRINT_pulse_num < self.number_of_SPRINT_pulses_per_seq:
            #         self.num_of_SPRINT_reflections_per_seq_N[(element-1) // self.sequence_len][SPRINT_pulse_num] += 1
            #         self.num_of_SPRINT_transmissions_per_seq_S[(element-1) // self.sequence_len][SPRINT_pulse_num] += 1

    """
        Dor: "Method was built so we can debug coherence"
        In this method we count the clicks on the Dark-Port and Bright-Port after the detection pulses

        Sequences:19231, Sequence-Len:520, End-of-detection-Pulse:510
        +-------------+--------------+--------------+--------------+--------+---+
        +      0      |      19      |      211     |      23      |    54  |   | ...  +
        +-------------+--------------+--------------+--------------+--------+---+
        <- 520 [ns] -> <- 520 [ns] -> <- 520 [ns] -> <- 520 [ns] -> <-510-><-10->     
    """

    def divide_tt_to_reflection_trans_extended(self):
        '''
        A function designed to count the number of photons reflected or transmitted for each sequence, such that,
        for the detection pulses the number of photons will be accumulated for each sequence and for the SPRINT
        pulses there will be number of reflected or transmitted photons for each SPRINT pulse.
        :param sprint_pulse_len: the length in [ns] of the SPRINT pulses in the sequence.
        :param num_of_det_pulses: the number of detection pulses in the sequence.
        :return:
        '''

        self.num_of_det_reflections_per_seq_S = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_det_reflections_per_seq_N = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_det_transmissions_per_seq_S = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_det_transmissions_per_seq_N = np.zeros(self.number_of_PNSA_sequences)

        self.num_of_det_reflections_per_seq_S_,\
        self.num_of_det_reflections_per_seq_N_, \
        self.num_of_det_transmissions_per_seq_S_, \
        self.num_of_det_transmissions_per_seq_N_ = [
            [
                [
                    [] for _ in range(self.number_of_detection_pulses_per_seq)
                ]
                for _ in range(self.number_of_PNSA_sequences)
            ] for _ in range(4)
        ]

        self.num_of_SPRINT_reflections_per_seq_S_, \
        self.num_of_SPRINT_reflections_per_seq_N_, \
        self.num_of_SPRINT_transmissions_per_seq_S_, \
        self.num_of_SPRINT_transmissions_per_seq_N_, \
        self.num_of_BP_counts_per_seq_in_SPRINT_pulse, \
        self.num_of_DP_counts_per_seq_in_SPRINT_pulse = [
            [
                [
                    [] for _ in range(self.number_of_SPRINT_pulses_per_seq)
                ]
                for _ in range(self.number_of_PNSA_sequences)
            ] for _ in range(6)
        ]

        # N direction:
        self.N_tt = np.array(self.tt_N_measure)
        tt_inseq_ = self.N_tt % self.sequence_len
        seq_num_ = self.N_tt // self.sequence_len
        for (element, tt_inseq, seq_num) in zip(self.N_tt, tt_inseq_, seq_num_):
            # TODO: Q: I assume this is a time-window to ignore some time on start/end, how did we decide on this time?
            #  TODO: Q: what is the meaning of this specific time-window?
            if (element > int(0.6e6)) and (element < int(self.M_time - 0.4e6)):
                for indx, tup in enumerate(self.sorted_pulses):
                    # if all(x == 0 for x in Config.PNSA_Exp_Square_samples_Early[tup[0]:tup[1]]):
                    # if (Config.PNSA_Exp_Square_samples_Early[tup[0]] == 0) or \
                    #         (self.number_of_SPRINT_pulses_per_seq > 1):
                    start_pulse_time_in_seq = tup[0]
                    end_pulse_time_in_seq = tup[1]
                    # else:
                    #     start_pulse_time_in_seq = tup[0] + Config.MZ_delay
                    #     end_pulse_time_in_seq = tup[1] + Config.MZ_delay
                    if (tt_inseq >= start_pulse_time_in_seq) and (tt_inseq <= end_pulse_time_in_seq):
                        if indx < self.number_of_detection_pulses_per_seq:
                            if tup[2] == 'N':
                                self.num_of_det_transmissions_per_seq_N_[seq_num][indx].append(element)
                            if tup[2] == 'S':
                                self.num_of_det_reflections_per_seq_S_[seq_num][indx].append(element)
                        else:
                            ind = indx - self.number_of_detection_pulses_per_seq
                            if tup[2] == 'n':
                                self.num_of_SPRINT_transmissions_per_seq_N_[seq_num][ind].append(element)
                            if tup[2] == 's':
                                self.num_of_SPRINT_reflections_per_seq_S_[seq_num][ind].append(element)

                if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                    self.num_of_det_reflections_per_seq_S[seq_num] += self.filter_S[tt_inseq]
                    self.num_of_det_transmissions_per_seq_N[seq_num] += self.filter_N[tt_inseq]

        # Bright port direction:
        self.BP_tt = np.array(self.tt_BP_measure)
        tt_inseq_ = self.BP_tt % self.sequence_len
        seq_num_ = self.BP_tt // self.sequence_len
        for (element, tt_inseq, seq_num) in zip(self.BP_tt, tt_inseq_, seq_num_):
            # TODO: Q: I assume this is a time-window to ignore some time on start/end, how did we decide on this time?
            #  TODO: Q: what is the meaning of this specific time-window?
            if (element > int(0.6e6)) and (element < int(self.M_time - 0.4e6)):
                for indx, tup in enumerate(self.sorted_pulses):
                    # if all(x == 0 for x in Config.PNSA_Exp_Square_samples_Early[tup[0]:tup[1]]):
                    if (Config.PNSA_Exp_Square_samples_Early[tup[0]] == 0) or \
                            (self.number_of_SPRINT_pulses_per_seq > 1):
                        start_pulse_time_in_seq = tup[0]
                        end_pulse_time_in_seq = tup[1]
                    else:
                        start_pulse_time_in_seq = tup[0] + Config.MZ_delay
                        end_pulse_time_in_seq = tup[1] + Config.MZ_delay
                    if (tt_inseq >= start_pulse_time_in_seq) and (tt_inseq <= end_pulse_time_in_seq):
                        if indx < self.number_of_detection_pulses_per_seq:
                            if tup[2] == 'N':
                                self.num_of_det_transmissions_per_seq_N_[seq_num][indx].append(element)
                            if tup[2] == 'S':
                                self.num_of_det_reflections_per_seq_S_[seq_num][indx].append(element)
                        else:
                            ind = indx - self.number_of_detection_pulses_per_seq
                            if tup[2] == 'n':
                                self.num_of_SPRINT_transmissions_per_seq_N_[seq_num][ind].append(element)
                            if tup[2] == 's':
                                self.num_of_SPRINT_reflections_per_seq_S_[seq_num][ind].append(element)
                            self.num_of_BP_counts_per_seq_in_SPRINT_pulse[seq_num][ind].append(element)

                if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                    self.num_of_det_reflections_per_seq_S[seq_num] += self.filter_S[tt_inseq]
                    self.num_of_det_transmissions_per_seq_N[seq_num] += self.filter_N[tt_inseq]

        # Dark port direction:
        self.DP_tt = np.array(self.tt_DP_measure)
        tt_inseq_ = self.DP_tt % self.sequence_len
        seq_num_ = self.DP_tt // self.sequence_len
        for (element, tt_inseq, seq_num) in zip(self.DP_tt, tt_inseq_, seq_num_):
            # TODO: Q: I assume this is a time-window to ignore some time on start/end, how did we decide on this time?
            # TODO: Q: what is the meaning of this specific time-window?
            if (element > int(0.6e6)) and (element < int(self.M_time - 0.4e6)):
                for indx, tup in enumerate(self.sorted_pulses):
                    # if all(x == 0 for x in Config.PNSA_Exp_Square_samples_Early[tup[0]:tup[1]]):
                    if (Config.PNSA_Exp_Square_samples_Early[tup[0]] == 0) or \
                            (self.number_of_SPRINT_pulses_per_seq > 1):
                        start_pulse_time_in_seq = tup[0]
                        end_pulse_time_in_seq = tup[1]
                    else:
                        start_pulse_time_in_seq = tup[0] + Config.MZ_delay
                        end_pulse_time_in_seq = tup[1] + Config.MZ_delay
                    if (tt_inseq >= start_pulse_time_in_seq) and (tt_inseq <= end_pulse_time_in_seq):
                        if indx < self.number_of_detection_pulses_per_seq:
                            if tup[2] == 'N':
                                self.num_of_det_transmissions_per_seq_N_[seq_num][indx].append(element)
                            if tup[2] == 'S':
                                self.num_of_det_reflections_per_seq_S_[seq_num][indx].append(element)
                        else:
                            ind = indx - self.number_of_detection_pulses_per_seq
                            if tup[2] == 'n':
                                self.num_of_SPRINT_transmissions_per_seq_N_[seq_num][ind].append(element)
                            if tup[2] == 's':
                                self.num_of_SPRINT_reflections_per_seq_S_[seq_num][ind].append(element)
                            self.num_of_DP_counts_per_seq_in_SPRINT_pulse[seq_num][ind].append(element)

                if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                    self.num_of_det_reflections_per_seq_S[seq_num] += self.filter_S[tt_inseq]
                    self.num_of_det_transmissions_per_seq_N[seq_num] += self.filter_N[tt_inseq]

        self.S_tt = np.array(self.tt_S_measure + self.tt_FS_measure)
        tt_inseq_ = self.S_tt % self.sequence_len
        seq_num_ = self.S_tt // self.sequence_len
        for (element, tt_inseq, seq_num) in zip(self.S_tt, tt_inseq_, seq_num_):
            if (element > int(0.6e6)) and (element < int(self.M_time - 0.4e6)):
                for indx, tup in enumerate(self.sorted_pulses):
                    if (tt_inseq >= tup[0]) and (tt_inseq <= tup[1]):
                        if indx < self.number_of_detection_pulses_per_seq:
                            if tup[2] == 'N':
                                self.num_of_det_reflections_per_seq_N_[seq_num][indx].append(element)
                            if tup[2] == 'S':
                                self.num_of_det_transmissions_per_seq_S_[seq_num][indx].append(element)
                        else:
                            ind = indx - self.number_of_detection_pulses_per_seq
                            if tup[2] == 'n':
                                self.num_of_SPRINT_reflections_per_seq_N_[seq_num][ind].append(element)
                            if tup[2] == 's':
                                self.num_of_SPRINT_transmissions_per_seq_S_[seq_num][ind].append(element)

                if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                    seq_num = (element - 1) // self.sequence_len
                    self.num_of_det_reflections_per_seq_N[seq_num] += self.filter_N[tt_inseq]
                    self.num_of_det_transmissions_per_seq_S[seq_num] += self.filter_S[tt_inseq]

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

        # self.num_of_BP_counts_per_n_sequences = np.zeros(self.number_of_PNSA_sequences)
        # self.num_of_DP_counts_per_n_sequences = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_BP_counts_per_n_sequences = np.zeros(self.number_of_PNSA_sequences // num_of_seq_per_count + 1)
        self.num_of_DP_counts_per_n_sequences = np.zeros(self.number_of_PNSA_sequences // num_of_seq_per_count + 1)
        self.num_of_S_counts_per_n_sequences = np.zeros(self.number_of_PNSA_sequences // num_of_seq_per_count + 1)

        for element in self.tt_BP_measure:
            tt_inseq = element % self.sequence_len
            if tt_inseq > self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                # print(element, (element - 1) // (self.sequence_len * num_of_seq_per_count))
                self.num_of_BP_counts_per_n_sequences[
                    (element - 1) // (self.sequence_len * num_of_seq_per_count)] += \
                    np.ceil(Config.PNSA_Exp_Square_samples_Late[tt_inseq])

        for element in self.tt_DP_measure:
            tt_inseq = element % self.sequence_len
            if tt_inseq > self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                # print(element, (element - 1) // (self.sequence_len * num_of_seq_per_count))
                try:
                    self.num_of_DP_counts_per_n_sequences[
                        (element - 1) // (self.sequence_len * num_of_seq_per_count)] += \
                        np.ceil(Config.PNSA_Exp_Square_samples_Late[tt_inseq])
                except Exception as err:
                    print(err)

        for element in self.tt_FS_measure + self.tt_S_measure:
            tt_inseq = element % self.sequence_len
            if tt_inseq > self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                # print(element, (element - 1) // (self.sequence_len * num_of_seq_per_count))
                try:
                    # lll = (element - 1) // (self.sequence_len * num_of_seq_per_count)
                    # if (element - 1) // (self.sequence_len * num_of_seq_per_count) > 385:
                    #     lll = 385
                    # self.num_of_S_counts_per_n_sequences[lll] += np.ceil(Config.PNSA_Exp_Square_samples_Late[tt_inseq])
                    self.num_of_S_counts_per_n_sequences[
                        (element - 1) // (self.sequence_len * num_of_seq_per_count)] += \
                        np.ceil(Config.PNSA_Exp_Square_samples_Late[tt_inseq])
                except Exception as err:
                    print(err)

    def plot_folded_tt_histogram(self):

        plt.figure()
        plt.plot(sum(np.reshape(self.tt_histogram_N, [int(len(self.tt_histogram_N) / 9), 9])), label='tt_hist_N')
        plt.plot(sum(np.reshape(self.tt_histogram_S, [int(len(self.tt_histogram_S) / 9), 9])), label='tt_hist_S')
        plt.legend()

    def find_transit_events_transitexp(self, transit_condition, tt_vector,
                                       all_transits_batch, tt_transit_events_accumulated):
        """
         :param tt_vector: The tt's vector on which we want to find transits
         :param transit_condition: consists of two integers: t - the max time allowed between clicks, n - the number of consecutive clicks

         The function checks if there are n consecutive clicks, within a time-frame that does not exceed t
        """
        transit_time_threshold = transit_condition[0]
        transit_counts_threshold = transit_condition[1]

        # Find transits and build histogram:
        current_transit = []
        all_transits = []
        tt_transit_events = np.zeros(int(self.M_window/self.bin_size))

        for t in tt_vector:
            if not current_transit:  # if the array is empty
                current_transit.append(t)
            elif (t - current_transit[-1]) < transit_time_threshold:
                current_transit.append(t)
            elif len(current_transit) >= transit_counts_threshold:
                all_transits.append(current_transit)
                current_transit = [t]
            else:
                current_transit = [t]

        tt_transit_events[[int(i//self.bin_size) for i in [vec for elem in all_transits for vec in elem]]] += 1
        tt_transit_events_accumulated = tt_transit_events_accumulated + tt_transit_events

        if all_transits:  # DROR: why do we need this?
            all_transits_batch = all_transits_batch[-(self.N - 1):] + [all_transits]

        return tt_transit_events, tt_transit_events_accumulated,all_transits_batch

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
        detection_pulse_range = np.zeros([self.number_of_PNSA_sequences, num_of_det_pulses], dtype=int)
        sprint_pulse_range = np.zeros([self.number_of_PNSA_sequences, num_of_sprint_pulses], dtype=int)

        # define ranges for tt_histogram (a vector that counts number of photons in each pulse)
        for i in range(self.number_of_PNSA_sequences):
            detection_pulse_range[i] = \
                list(range(i * num_of_pulses, i * num_of_pulses + num_of_det_pulses))
            # sprint pulse range starts from last detection pulse andfo num_of_sprint_pulses
            sprint_pulse_range[i] = \
                list(range(int(detection_pulse_range[i][-1]) + 1,
                           int(detection_pulse_range[i][-1]) + 1 + num_of_sprint_pulses))
        # find transits and sprint events
        for j in range(self.number_of_PNSA_sequences - 2):
            if \
                    sum(self.tt_histogram_reflection[detection_pulse_range[j]]) >= detection_condition[0] and \
                            sum(self.tt_histogram_reflection[detection_pulse_range[j + 1]]) >= detection_condition[1]:
                num_of_detected_atom += 1
                transit_sequences_init_tt.append(j * self.sequence_len)
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
                # TODO: ask dor (08.01.24) - what happens at [0,4,0]? and why including the middle at [2,0,2]?
                # adding to current transit the indices from first element satisfing the condition to the last element checked.
                # for example:
                # if the condition is [1,1,1] and for i=7000 the reflections were [0(i=7000),1 (i=7001),1 (i=7002)]
                # than current transit would add [7001,7002]
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

    def find_transit_events_Old(self, N, transit_time_threshold, transit_counts_threshold, tt_histogram,tt_histogram_bin_size):
        # this list appends the tt's that are distanced less than transit_time_threshold.gets empty every new transit.
        current_transit = []
        # list of all transits that go through the constrain
        self.all_transits = []
        # transits time relative to the first click in transit
        all_transits_aligned_first = []
        # runnning over all tt's binned
        for index, value in enumerate(tt_histogram):
            if not current_transit and value:  # if the array is empty start adding tt's
                current_transit.append(index)
            elif value:
                if ((index - current_transit[0]) *tt_histogram_bin_size) < transit_time_threshold:
                    current_transit.append(index)
                elif sum([tt_histogram[i] for i in current_transit]) > transit_counts_threshold:
                    self.all_transits.append(current_transit)
                    all_transits_aligned_first.append([x - current_transit[0] for x in current_transit])
                else:
                    # Finding if there any index that was saved to current transit and is close enough to the new index
                    t = [i for i, elem in enumerate(current_transit) if ((index - elem) * 480) < transit_time_threshold]
                    if t:
                        current_transit = current_transit[t[0]:] + [index]
                    else:
                        current_transit = [index]

        self.tt_S_transit_events[[i for i in [vec for elem in self.all_transits for vec in elem]]] += 1

        if self.all_transits:
            self.all_transits_batch = self.all_transits_batch[-(N - 1):] + [self.all_transits]


    def find_transit_events(self, cond=None, minimum_number_of_seq_detected=2):
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

        for i in range(self.number_of_PNSA_sequences - len(cond) + 1):
            cond_check = (self.num_of_det_reflections_per_seq[i:(i + len(cond))] >= cond).astype(int)
            if sum(cond_check) >= minimum_number_of_seq_detected:
                # TODO: ask dor (08.01.24) - what happens at [0,4,0]? and why including the middle at [2,0,2]?
                # adding to current transit the indices from first element satisfing the condition to the last element checked.
                # for example:
                # if the condition is [1,1,1] and for i=7000 the reflections were [0(i=7000),1 (i=7001),1 (i=7002)]
                # than current transit would add [7001,7002]
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
                self.all_transits_seq_indx.append(current_transit)
                current_transit = []
        if len(current_transit) > 1:
            current_transit = current_transit[
                              :np.where(self.num_of_det_reflections_per_seq[current_transit] >= min(cond))[0][-1] + 1]
            if self.all_transits_seq_indx:
                if bool(set(current_transit) & set(self.all_transits_seq_indx[-1])):
                    current_transit = self.all_transits_seq_indx[-1] + current_transit[1:]
                    self.all_transits_seq_indx = self.all_transits_seq_indx[:-1]
            self.all_transits_seq_indx.append(current_transit)

    def analyze_SPRINT_data_points(self, all_transits_seq_indx, SPRINT_pulse_number=[1], background=False):
        '''
        For a given vector of sequence indexes, find the relevant data points for SPRINT pulse and analyze results.
        :param all_transits_seq_indx: The vector of indexes of sequences that need to be analyzed.
        :param SPRINT_pulse_number: The SPRINT pulse number for which we want to check the results.
        :param background: If the check is for background data and not for actual transits, the "potential data" check,
                           for which we condition the relevance of the data if we get at least 1 reflection from the
                           last detection pulse, is irrelevant.
        '''
        reflection_SPRINT_data = []  # Array of vectors with data on the number of reflections per SPRINT pulse in sequence.
        transmission_SPRINT_data = []  # Array of vectors with data on the number of transmissions per SPRINT pulse in sequence.
        seq_with_data_points = []
        BP_counts_SPRINT_data = []
        DP_counts_SPRINT_data = []
        if len(self.sorted_pulses) > 0:
            for transit in all_transits_seq_indx:
                for seq_indx in transit[:-1]:
                    # Checking potential for data point by looking for a single photon at the reflection of the last
                    # detection pulse:
                    potential_data = False if not background else True
                    if self.sorted_pulses[self.number_of_detection_pulses_per_seq-1][2] == 'N' and \
                       len(self.num_of_det_reflections_per_seq_N_[seq_indx][-1]) >= 1:
                        potential_data = True
                        # self.potential_data += 1
                    elif self.sorted_pulses[self.number_of_detection_pulses_per_seq-1][2] == 'S' and \
                            len(self.num_of_det_reflections_per_seq_S_[seq_indx][-1]) >= 1:
                        potential_data = True
                        # self.potential_data += 1
                    # Getting SPRINT data if the SPRINT pulse has one photon in the reflection or transmission
                    if potential_data and len(self.sorted_pulses) > self.number_of_detection_pulses_per_seq:
                        if self.sorted_pulses[self.number_of_detection_pulses_per_seq-1+SPRINT_pulse_number[0]][2] == 'n':
                            transmissions = 0
                            reflections = 0
                            for sprint_pulse in SPRINT_pulse_number:
                                transmissions += len(self.num_of_SPRINT_transmissions_per_seq_N_[seq_indx][sprint_pulse-1])
                                reflections += len(self.num_of_SPRINT_reflections_per_seq_N_[seq_indx][sprint_pulse-1])
                            if (transmissions + reflections) == 1:
                                seq_with_data_points.append(seq_indx)
                                reflection_SPRINT_data.append(reflections)
                                transmission_SPRINT_data.append(transmissions)
                            if (transmissions + reflections) > 0:
                                BP_counts_SPRINT_data.append(len(self.num_of_BP_counts_per_seq_in_SPRINT_pulse[seq_indx][SPRINT_pulse_number[-1]-1]))
                                DP_counts_SPRINT_data.append(len(self.num_of_DP_counts_per_seq_in_SPRINT_pulse[seq_indx][SPRINT_pulse_number[-1]-1]))
                        elif self.sorted_pulses[self.number_of_detection_pulses_per_seq-1+SPRINT_pulse_number[0]][2] == 's':
                            transmissions = 0
                            reflections = 0
                            for sprint_pulse in SPRINT_pulse_number:
                                transmissions += len(self.num_of_SPRINT_transmissions_per_seq_S_[seq_indx][sprint_pulse-1])
                                reflections += len(self.num_of_SPRINT_reflections_per_seq_S_[seq_indx][sprint_pulse-1])
                            if (transmissions + reflections) == 1:
                                seq_with_data_points.append(seq_indx)
                                reflection_SPRINT_data.append(reflections)
                                transmission_SPRINT_data.append(transmissions)
                            if (transmissions + reflections) > 0:
                                BP_counts_SPRINT_data.append(len(self.num_of_BP_counts_per_seq_in_SPRINT_pulse[seq_indx][SPRINT_pulse_number[-1]-1]))
                                DP_counts_SPRINT_data.append(len(self.num_of_DP_counts_per_seq_in_SPRINT_pulse[seq_indx][SPRINT_pulse_number[-1]-1]))
        return seq_with_data_points, reflection_SPRINT_data, transmission_SPRINT_data, BP_counts_SPRINT_data, DP_counts_SPRINT_data


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

    def get_pulses_location_in_seq_DEP(self, delay, seq=Config.PNSA_Exp_Gaussian_samples_S,
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

    def get_avg_num_of_photons_in_seq_pulses(self, seq, pulse_loc, tt_measure, efficiency):
        avg_num_of_photons_in_seq_pulses = []
        try:
            real_number_of_seq = math.ceil(max(tt_measure) / len(Config.PNSA_Exp_Gaussian_samples_S))
            # self.logger.debug('Real number of seq = %d' %real_number_of_seq)
        except:
            real_number_of_seq = self.number_of_PNSA_sequences
            # self.logger.debug('Max number of seq')
        for t in pulse_loc:
            avg_num_of_photons_in_seq_pulses.append((sum(seq[t[0]:t[1]]) + seq[t[1]]) / (
                        real_number_of_seq * efficiency))  # Sagnac configuration efficiency 16.7%
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

    def number_of_pulses_per_seq(self):
        '''

        :return:
        '''
        detection_pulses_per_seq = 0
        SPRINT_pulses_per_seq = 0
        for tup in self.sorted_pulses:
            if (tup[2] == 'N') or (tup[2] == 'S'):
                detection_pulses_per_seq += 1
            elif (tup[2] == 'n') or (tup[2] == 's'):
                SPRINT_pulses_per_seq += 1
        return detection_pulses_per_seq, SPRINT_pulses_per_seq

    def tt_histogram(self, num_of_bins):
        self.tt_histogram_N, _ = np.histogram(self.tt_N_measure_Total, bins=num_of_bins,range =(0,self.M_window))
        self.tt_histogram_S, _ = np.histogram(self.tt_S_measure_Total, bins=num_of_bins,range =(0,self.M_window))

    def fold_tt_histogram(self, exp_sequence_len):

        st = time.time()

        # Zero all vectors
        # TODO: These are used globally, so maybe it's not the best place to zero them - should be outside this function (somewhere)
        self.folded_tt_S = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_N = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_BP = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_DP = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_FS = np.zeros(exp_sequence_len, dtype=int)

        self.folded_tt_S_directional = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_N_directional = np.zeros(exp_sequence_len, dtype=int)

        self.folded_tt_BP_timebins = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_DP_timebins = np.zeros(exp_sequence_len, dtype=int)

        # Bin all the time-tags into buckets for all the time-tag sequences we care about
        for x in self.tt_S_measure:
            self.folded_tt_S[x % exp_sequence_len] += 1
        for x in self.tt_N_measure:
            self.folded_tt_N[x % exp_sequence_len] += 1
        for x in self.tt_BP_measure:
            self.folded_tt_BP[x % exp_sequence_len] += 1
        for x in self.tt_DP_measure:
            self.folded_tt_DP[x % exp_sequence_len] += 1
        for x in self.tt_FS_measure:
            self.folded_tt_FS[x % exp_sequence_len] += 1

        # for x in [elem for elem in self.tt_S_measure]:
        #     self.folded_tt_S[x % exp_sequence_len] += 1
        # for x in [elem for elem in self.tt_N_measure]:
        #     self.folded_tt_N[x % exp_sequence_len] += 1
        # for x in [elem for elem in self.tt_BP_measure]:
        #     self.folded_tt_BP[x % exp_sequence_len] += 1
        # for x in [elem for elem in self.tt_DP_measure]:
        #     self.folded_tt_DP[x % exp_sequence_len] += 1
        # for x in [elem for elem in self.tt_FS_measure]:
        #     self.folded_tt_FS[x % exp_sequence_len] += 1

        # Directional is South + Fast Switch
        self.folded_tt_S_directional = (np.array(self.folded_tt_S) + np.array(self.folded_tt_FS))
        # self.folded_tt_N_directional = self.folded_tt_N
        # TODO: ask dor -  switched folded_tt_N with folded_tt_N_directional
        self.folded_tt_N_directional[:self.end_of_det_pulse_in_seq] = \
            np.array(self.folded_tt_N[:self.end_of_det_pulse_in_seq]) \
            + np.array(self.folded_tt_BP[:self.end_of_det_pulse_in_seq]) \
            + np.array(self.folded_tt_DP[:self.end_of_det_pulse_in_seq])

        self.folded_tt_BP_timebins[self.end_of_det_pulse_in_seq:] = self.folded_tt_BP[self.end_of_det_pulse_in_seq:]
        self.folded_tt_DP_timebins[self.end_of_det_pulse_in_seq:] = self.folded_tt_DP[self.end_of_det_pulse_in_seq:]
        if self.pulses_location_in_seq_A or ((Config.sprint_pulse_amp_N[0] > 0) & (len(Config.sprint_pulse_amp_N) > 1)):
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

        et = time.time()
        self.info(f'fold_tt_histogram took {et-st}')

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

    def load_data_for_playback__old_files(self):
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
        # fluorescence_str = '%.2f' % self.fluorescence_average
        fluorescence_str = '%d' % int(self.fluorescence_average)
        time_formatted = time.strftime("%Y/%m/%d, %H:%M:%S")
        status_str = f'[Warm Up: {self.warm_up_cycles}]' if self.warm_up else f'# {self.counter} ({self.repetitions})'
        photons_str = '' if self.warm_up else f', #Photons/[us]: {self.sum_for_threshold}'

        self.logger.info(f'{time_formatted}: {self.experiment_type} {status_str}, Eff: {locking_efficiency_str}, Flr: {fluorescence_str}{photons_str}')

    def experiment_mainloop_delay(self):
        super().experiment_mainloop_delay()
        pass

    def prepare_figures(self):

        super().prepare_figures()

        matplotlib.mathtext.SHRINK_FACTOR = 0.4
        matplotlib.mathtext.GROW_FACTOR = 1 / 0.4

        # Uncomment the below if you want the view to open maximized
        #self.maximize_figure()

        return

    # -------------------------------------------------------------------------
    # Experiment Plots
    # -------------------------------------------------------------------------

    def plots_handler__header(self, ax):
        """
        Plots a TEXT that is the experiment header line
        """

        pause_str = ' , PAUSED!' if self.pause_flag else ''
        ref_str = '%.2f' % self.sum_for_threshold
        eff_str = '%.1f%%' % (self.counter * 100 / self.repetitions)
        exp_str = r'$\bf{' + self.experiment_type + '}$'
        if hasattr(self, 'fluorescence_flag'):
            if self.fluorescence_flag:  # @@@
                # flr_str = '$Flr: %.2f$' % self.fluorescence_average
                flr_str = '$Flr: %d$' % int(self.fluorescence_average)
            else:
                # flr_str = r'$\bf{Flr: %.2f}$' % self.fluorescence_average
                flr_str = r'$\bf{Flr: %d}$' % int(self.fluorescence_average)
        else:
            flt_str = "Flr: N/A"

        if self.lock_err is not None:
            if self.lock_err_flag:
                lck_str = '$\Delta_{lock}: %.1f$' % self.lock_err
            else:
                lck_str = r'$\bf{\Delta_{lock}: %.1f}$' % self.lock_err
        else:
            lck_str = r'LckErr: N/A'

        if self.k_ex is not None:
            if self.k_ex_flag:
                k_ex_str = '$\kappa_{ex}: %.1f$' % self.k_ex
            else:
                k_ex_str = r'$\bf{\kappa_{ex}: %.1f}$' % self.k_ex
        else:
            k_ex_str = "K_ex: N/A"

        status_str = f'[Warm Up: {self.warm_up_cycles}]' if self.warm_up else f'# {self.counter} ({eff_str})'
        playback_str = 'P: ' if self.playback['active'] else ''

        header_text = f'{playback_str} {self.experiment_type}, {status_str} - {flr_str}, {lck_str}, {k_ex_str} {pause_str}'

        # Threshold Box:
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        flag_color = 'green' if self.acquisition_flag else 'red'
        props_thresholds = dict(boxstyle='round', edgecolor=flag_color, linewidth=2, facecolor=flag_color, alpha=0.5)

        ax.text(0, 1.4, header_text, transform=ax.transAxes, fontsize=26, verticalalignment='top', bbox=props_thresholds)

        pass

    def plots_handler__left_sidebar(self, ax):

        latched_detectors = self.latched_detectors()
        saturated_detectors = self.saturated_detectors()

        subplots_view = self.bdplots.get_subplots_view()

        x0 = -0.15
        y0 = 0.5  # or 1.0 in zoom mode
        w = 0.4  # or 0.2 in zoom mode
        if subplots_view > 0:  # Single Plot
            x0 = -0.08
            y0 = 1.0
            w = 0.14

        for i, det in enumerate(self.Num_Of_dets):
            pad = 1.2
            num_clicks = len(self.tt_measure[i])
            # num_clicks = len(self.streams[f'Detector_{det}_Timetags']['results'][0])
            text = f'{self.detectors_names[i]}-{det}\n({num_clicks - 1})'
            det_color = 'red' if i in latched_detectors else 'green'
            if i in latched_detectors:
                det_color = 'red'
            elif i in saturated_detectors:
                det_color = '#ffc710'
            props = dict(boxstyle=f"circle,pad={pad}", edgecolor=det_color, linewidth=2, facecolor=det_color, alpha=0.5)

            y = y0 - i * w
            ax.text(x0, y, text, ha="center", va="center", transform=ax.transAxes, fontsize=8, bbox=props)

        pass

    def plots_handler__right_sidebar(self, ax):

        x0 = 2.22
        y0 = 0.5  # or 1.0 in zoom mode
        w = 0.4  # or 0.2 in zoom mode

        subplots_view = self.bdplots.get_subplots_view()

        exp_flag = 'green' if self.exp_flag else 'yellow'
        connection_status = 'green' if (self.playback['active'] or self.bdsocket.is_connected('cavity_lock')) else 'red'
        with_atoms = 'green' if self.MOT_on else 'red'
        spectrum_save_flag = 'green' if ('cavity_lock' in self.comm_messages and self.comm_messages['cavity_lock']['save_succeeded'] and connection_status) else 'red'
        snspds_flag = 'orange'

        flags = [
            ('Exp Flag       ', exp_flag),
            ('Atoms           ', with_atoms),
            ('Locker Sync     ', connection_status),
            ('Spectrum Save', spectrum_save_flag),
            ('SNSPDs Temp', snspds_flag)
        ]

        for i in range(0, len(flags)):
            flag = flags[i]
            pad = 1.2
            text = flag[0]
            color = flag[1]
            y = y0 - i * w

            #props = dict(boxstyle=f"circle,pad={pad}", edgecolor=color, linewidth=2, facecolor=color, alpha=0.5)
            #ax.text(x0, y, text, ha="center", va="center", transform=ax.transAxes, fontsize=8, bbox=props)

            props = dict(boxstyle='round', edgecolor=color, linewidth=2, facecolor=color, alpha=0.5)
            ax.text(x0, y, text, transform=ax.transAxes, verticalalignment='top', fontsize=12, bbox=props, ha='center', va='center')

        pass

    def plots_handler__binned_tags_live(self, subplot_def):

        ax = subplot_def["ax"]
        # ax.plot(self.tt_histogram_N, label='"N" detectors')
        ax.plot(self.tt_histogram_S, label='"S" detectors')
        pass

    def plots_handler__binned_tags_acc(self, subplot_def):

        ax = subplot_def["ax"]
        ax.plot(np.average(self.batcher["tt_histogram_N_batch"][1:], axis=0), label='"N" detectors', alpha=0.6)
        ax.plot(np.average(self.batcher["tt_histogram_S_batch"][1:], axis=0), label='"S" detectors', alpha=0.6)
        pass


    def plots_handler__binned_tags_live_OLD(self, subplot_def):

        ax = subplot_def["ax"]

        ax.plot(self.folded_tt_N_directional, label='"N" detectors')
        ax.plot(self.folded_tt_S_directional, label='"S" detectors')
        ax.plot(self.folded_tt_BP_timebins, label='"BP" detectors')
        ax.plot(self.folded_tt_DP_timebins, label='"DP" detectors')
        ax.plot((self.filter_S) * max(self.folded_tt_N_directional + self.folded_tt_S_directional),
                   '--', color='orange', label='Filter "S"')
        for i in range(len(self.Num_of_photons_txt_box_y_loc_live)):
            ax.text(self.Num_of_photons_txt_box_x_loc.tolist()[i], self.Num_of_photons_txt_box_y_loc_live[i],
                       '%.2f' % self.avg_num_of_photons_per_pulse_live[i],
                       horizontalalignment='center', fontsize=12, fontweight='bold', family=['Comic Sans MS'])
        for j in range(len(self.Num_of_photons_txt_box_y_loc_live_MZ)):
            ax.text(self.Num_of_photons_txt_box_x_loc_for_MZ_ports.tolist()[j],
                       self.Num_of_photons_txt_box_y_loc_live_MZ[j][0],
                       '%.2f' % self.avg_num_of_photons_per_pulse_live_MZ[j][0],
                       horizontalalignment='center', fontsize=12, fontweight='bold', family=['Comic Sans MS'],
                       color='#2ca02c')
            ax.text(self.Num_of_photons_txt_box_x_loc_for_MZ_ports.tolist()[j],
                       self.Num_of_photons_txt_box_y_loc_live_MZ[j][1],
                       '%.2f' % self.avg_num_of_photons_per_pulse_live_MZ[j][1],
                       horizontalalignment='center', fontsize=12, fontweight='bold', family=['Comic Sans MS'],
                       color='#d62728')

        pass

    def plots_handler__binned_tags_averaged(self, subplot_def):

        ax = subplot_def["ax"]

        ax.plot(self.folded_tt_N_directional_cumulative_avg, label='"N" detectors')
        ax.plot(self.folded_tt_S_directional_cumulative_avg, label='"S" detectors')
        ax.plot(self.folded_tt_BP_timebins_cumulative_avg, label='"BP" detectors')
        ax.plot(self.folded_tt_DP_timebins_cumulative_avg, label='"DP" detectors')
        ax.plot((self.filter_N) * max(
            self.folded_tt_N_directional_cumulative_avg + self.folded_tt_S_directional_cumulative_avg),
                   '--b', label='Filter "N"')
        for i in range(len(self.Num_of_photons_txt_box_y_loc)):
            ax.text(self.Num_of_photons_txt_box_x_loc.tolist()[i], self.Num_of_photons_txt_box_y_loc[i],
                       '%.2f' % self.avg_num_of_photons_per_pulse[i],
                       horizontalalignment='center', fontsize=12, fontweight='bold', family=['Comic Sans MS'])
        for j in range(len(self.Num_of_photons_txt_box_y_loc_MZ)):
            ax.text(self.Num_of_photons_txt_box_x_loc_for_MZ_ports.tolist()[j],
                       self.Num_of_photons_txt_box_y_loc_MZ[j][0],
                       '%.2f' % self.avg_num_of_photons_per_pulse_MZ[j][0],
                       horizontalalignment='center', fontsize=12, fontweight='bold', family=['Comic Sans MS'],
                       color='#2ca02c')
            ax.text(self.Num_of_photons_txt_box_x_loc_for_MZ_ports.tolist()[j],
                       self.Num_of_photons_txt_box_y_loc_MZ[j][1],
                       '%.2f' % self.avg_num_of_photons_per_pulse_MZ[j][1],
                       horizontalalignment='center', fontsize=12, fontweight='bold', family=['Comic Sans MS'],
                       color='#d62728')

        # -----------------------------------------
        # SPRINT results box ("Soccer Results")
        # -----------------------------------------

        SPRINT_reflections_without_transits = '%d' % sum(self.batcher['num_of_total_SPRINT_reflections_batch'])
        SPRINT_transmissions_without_transits = '%d' % sum(self.batcher['num_of_total_SPRINT_transmissions_batch'])
        if (sum(self.batcher['num_of_total_SPRINT_reflections_batch']) + sum(self.batcher['num_of_total_SPRINT_transmissions_batch'])) > 0:
            SPRINT_reflections_percentage_without_transits = (
                    '%.1f' % ((sum(self.batcher['num_of_total_SPRINT_reflections_batch']) * 100) /
                              (sum(self.batcher['num_of_total_SPRINT_reflections_batch']) + sum(self.batcher['num_of_total_SPRINT_transmissions_batch']))))
            SPRINT_transmissions_percentage_without_transits = (
                    '%.1f' % ((sum(self.batcher['num_of_total_SPRINT_transmissions_batch']) * 100) /
                              (sum(self.batcher['num_of_total_SPRINT_reflections_batch']) + sum(self.batcher['num_of_total_SPRINT_transmissions_batch']))))
        else:
            SPRINT_reflections_percentage_without_transits = '%.1f' % 0
            SPRINT_transmissions_percentage_without_transits = '%.1f' % 0

        SPRINT_reflections_with_transits = '%d' % sum(sum(self.batcher['reflection_SPRINT_data_batch'], []))
        # SPRINT_reflections = f'${SPRINT_reflections_with_transits}_{{({SPRINT_reflections_without_transits})}}$'
        SPRINT_reflections = f'${SPRINT_reflections_with_transits}_{{({SPRINT_reflections_percentage_without_transits}\%)}}$'
        SPRINT_reflections_text = '$R_{SPRINT}$'
        SPRINT_transmissions_with_transits = '%d' % sum(sum(self.batcher['transmission_SPRINT_data_batch'], []))
        # SPRINT_transmissions = f'${SPRINT_transmissions_with_transits}_{{({SPRINT_transmissions_without_transits})}}$'
        SPRINT_transmissions = f'${SPRINT_transmissions_with_transits}_{{({SPRINT_transmissions_percentage_without_transits}\%)}}$'
        SPRINT_transmissions_text = '$T_{SPRINT}$'
        SPRINT_Score = f'{SPRINT_reflections} - {SPRINT_transmissions}'
        table_vals = [SPRINT_reflections_text, SPRINT_Score, SPRINT_transmissions_text]
        SPRINT_text = f'{SPRINT_reflections_text} {SPRINT_reflections} - {SPRINT_transmissions} {SPRINT_transmissions_text}'
        props_SPRINT = dict(boxstyle='round', edgecolor='gray', linewidth=2, facecolor='gray', alpha=0.5)

        ax.text(1, 1.4, SPRINT_text, transform=ax.transAxes, fontsize=26, verticalalignment='top',
                   horizontalalignment='right', bbox=props_SPRINT)
        pass

    def plots_handler__mz_outputs(self, subplot_def):

        ax = subplot_def["ax"]

        ax.plot(self.MZ_BP_counts_res_value_0, label='MZ Bright port')
        ax.plot(self.MZ_DP_counts_res_value_0, label='MZ Dark port')
        ax.plot(self.MZ_BP_counts_res_value_0 - self.MZ_DP_counts_res_value_0, label='Dif ports')

        ax.set_ylim(0, 1.1 * np.max([self.MZ_BP_counts_res_value_0, self.MZ_DP_counts_res_value_0]))

        pass

    def plots_handler__phase_correction(self, subplot_def):

        ax = subplot_def["ax"]

        ax.tick_params(axis="y", labelcolor='#8c564b')
        ax.plot(self.Phase_Correction_vec, label='Phase correction values', color='#8c564b')
        ax.plot(self.Phase_Correction_min_vec, label='Phase correction values', color='#9467bd')

        pass

    def plots_handler__num_reflections_per_sequence(self, subplot_def):

        ax = subplot_def["ax"]

        max_reflect_avg = max(self.num_of_det_reflections_per_seq_accumulated / self.counter)
        max_reflect = max(self.num_of_det_reflections_per_seq)

        textstr_total_reflections_N = 'Total reflections per cycle "N" = %d' % (
            sum(self.num_of_det_reflections_per_seq_N),)

        textstr_total_reflections_S = 'Total reflections per cycle "S" = %d' % (
            sum(self.num_of_det_reflections_per_seq_S),)

        textstr_total_reflections_avg = 'Average reflections per cycle = %.2f' % (
            sum(self.num_of_det_reflections_per_seq_accumulated / self.counter),)

        textstr_total_transmission_avg = 'Average transmissions per cycle = %.2f' % (
            sum(self.num_of_det_transmissions_per_seq_accumulated / self.counter),)

        total_photons_per_seq_accumulated = sum(self.num_of_det_reflections_per_seq_accumulated) + \
                                            (sum(self.num_of_det_transmissions_per_seq_accumulated) / self.Cavity_transmission)
        if total_photons_per_seq_accumulated > 0:
            textstr_total_reflections_percentage = 'Average reflections percentage = %.1f%%' % (
                (sum(self.num_of_det_reflections_per_seq_accumulated) * 100) / total_photons_per_seq_accumulated,)
        else:
            textstr_total_reflections_percentage = 'Average reflections percentage = %.1f%%' % 0

        textstr_total_reflections = f'{textstr_total_reflections_N} \n{textstr_total_reflections_S} \n' \
                                    f'{textstr_total_reflections_avg} \n{textstr_total_transmission_avg} \n' \
                                    f'{textstr_total_reflections_percentage}' \

        ax.plot(self.num_of_det_reflections_per_seq_accumulated / self.counter,
                   label='Num of reflections per sequence')

        # TODO: original code did not have this condition. Can it be that max_relect is never zero?
        if max_reflect > 0:
            ax.plot(self.num_of_det_reflections_per_seq * 0.5 * max_reflect_avg / max_reflect * 0.3,
                       label='Num of reflections per sequence (Live)')

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.02, 0.95, textstr_total_reflections, transform=ax.transAxes, fontsize=14,
                   verticalalignment='top', bbox=props)

        pass

    def plots_handler__mz_outputs_around_experiment(self, subplot_def):

        ax = subplot_def["ax"]

        ax.plot(self.MZ_BP_counts_res_value_1, label='MZ BP counts before and after')
        ax.plot(self.MZ_DP_counts_res_value_1, label='MZ DP port counts before and after')
        ax.axvline(len(self.MZ_DP_counts_res_value_1) / 2, linestyle='--', color='red')

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        flag_color = 'green' if self.acquisition_flag else 'red'
        props_thresholds = dict(boxstyle='round', edgecolor=flag_color, linewidth=2, facecolor=flag_color, alpha=0.5)

        if self.MZ_infidelity_flag:
            props_MZ = props
        else:
            props_MZ = props_thresholds

        textstr_BP_DP_BA = 'Infidelity before = %.2f \n' % (self.Infidelity_before,) + \
                           'Infidelity after= %.2f \n' % (self.Infidelity_after,) + \
                           'S total counts MZ = %.2f ' % (self.MZ_S_tot_counts,)

        ax.text(0.05, 0.6, textstr_BP_DP_BA, transform=ax.transAxes, fontsize=14, verticalalignment='top',
                   bbox=props_MZ)

        pass

    def plots_handler__mz_outputs_during_experiment(self, subplot_def):

        ax = subplot_def["ax"]

        ax.plot(self.num_of_BP_counts_per_n_sequences, label='MZ BP counts per %d seq' % 50)
        ax.plot(self.num_of_DP_counts_per_n_sequences, label='MZ DP counts per %d seq' % 50)
        ax.plot(self.num_of_S_counts_per_n_sequences, label='"S" transmission counts per %d seq' % 50)
        ax.plot(self.num_of_S_counts_per_n_sequences + self.num_of_DP_counts_per_n_sequences + self.num_of_BP_counts_per_n_sequences,
            label='"S" total counts per %d seq' % 50)

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        props_SPRINT_coherence = dict(boxstyle='round', edgecolor='gray', linewidth=2, facecolor='gray', alpha=0.5)

        # Coherence results box
        if (sum(self.batcher['num_of_total_SPRINT_BP_counts_batch']) + sum(self.batcher['num_of_total_SPRINT_DP_counts_batch'])) > 0:
            SPRINT_BP_percentage_without_transits = (
                    '%.1f' % ((sum(self.batcher['num_of_total_SPRINT_BP_counts_batch']) * 100) /
                              (sum(self.batcher['num_of_total_SPRINT_BP_counts_batch']) + sum(self.batcher['num_of_total_SPRINT_DP_counts_batch']))))
            SPRINT_DP_percentage_without_transits = (
                    '%.1f' % ((sum(self.batcher['num_of_total_SPRINT_DP_counts_batch']) * 100) /
                              (sum(self.batcher['num_of_total_SPRINT_BP_counts_batch']) + sum(self.batcher['num_of_total_SPRINT_DP_counts_batch']))))
        else:
            SPRINT_BP_percentage_without_transits = '%.1f' % 0
            SPRINT_DP_percentage_without_transits = '%.1f' % 0

        SPRINT_BP_counts_with_transits = '%d' % sum(sum(self.batcher['BP_counts_SPRINT_data_batch'], []))
        SPRINT_BP_counts = f'${SPRINT_BP_counts_with_transits}_{{({SPRINT_BP_percentage_without_transits}\%)}}$'
        SPRINT_BP_counts_text = '$BP_{SPRINT}$'
        SPRINT_DP_counts_with_transits = '%d' % sum(sum(self.batcher['DP_counts_SPRINT_data_batch'], []))
        SPRINT_DP_counts = f'${SPRINT_DP_counts_with_transits}_{{({SPRINT_DP_percentage_without_transits}\%)}}$'
        SPRINT_DP_counts_text = '$DP_{SPRINT}$'
        SPRINT_Coherence_Score = f'{SPRINT_BP_counts} - {SPRINT_DP_counts}'
        table_vals = [SPRINT_BP_counts_text, SPRINT_Coherence_Score, SPRINT_DP_counts_text]

        SPRINT_Coherence_text = f'{SPRINT_BP_counts_text} {SPRINT_BP_counts} - {SPRINT_DP_counts} {SPRINT_DP_counts_text}'

        avg_BP = np.average(self.num_of_BP_counts_per_n_sequences)
        avg_DP = np.average(self.num_of_DP_counts_per_n_sequences)

        # TODO: temp workaround, REMOVE
        avg_BP = 0.001 if avg_BP == 0 else avg_BP
        textstr_BP_DP = 'Average Bright counts = %.2f \n' % (avg_BP,) + \
                        'Average Dark counts = %.2f \n' % (avg_DP,) + \
                        'Infidelity = %.2f' % (avg_DP / (avg_DP + avg_BP),)

        ax.text(0.05, 0.6, textstr_BP_DP, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax.text(0.05, 0.9, SPRINT_Coherence_text, transform=ax.transAxes, fontsize=18, verticalalignment='top',
                   bbox=props_SPRINT_coherence)
        pass

    def plots_handler__transits_per_sequence(self, subplot_def):

        ax = subplot_def["ax"]

        ax.plot(range(len(self.tt_transit_events_accumulated_N)), self.tt_transit_events_accumulated_N, label='Transit North Events Accumulated')
        ax.plot(range(len(self.tt_transit_events_accumulated_S)), self.tt_transit_events_accumulated_S, label='Transit South Events Accumulated')
        ax.plot(range(len(self.tt_transit_events_N)), self.tt_transit_events_N, label='Transit North Events Live')
        ax.plot(range(len(self.tt_transit_events_S)), self.tt_transit_events_S, label='Transit South Events Live')
        ax.set(xlabel='Sequence [#]', ylabel='Counts [Photons]')

        if self.number_of_transits_live:
            textstr_transit_counts_N = r'$N_{Transits_North} = %s $' % (sum(self.tt_transit_events_N),) + r'$[Counts]$'
            textstr_transit_counts_S = r'$N_{Transits_South} = %s $' % (sum(self.tt_transit_events_S),) + r'$[Counts]$'
        else:
            textstr_transit_counts_N = r'$N_{Transits} = %s $' % (0,) + r'$[Counts]$'
            textstr_transit_counts_S = r'$N_{Transits} = %s $' % (0,) + r'$[Counts]$'

        textstr_transit_event_counter_N = r'$N_{Transits Total North} = %s $' % (sum(self.tt_transit_events_accumulated_N),) + '[Counts]\n' \
                                        + textstr_transit_counts_N
        textstr_transit_event_counter_S = r'$N_{Transits Total South} = %s $' % (sum(self.tt_transit_events_accumulated_S),) + '[Counts]\n'\
                                          + textstr_transit_counts_S

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.02, 0.92, f'{textstr_transit_event_counter_S} \n {textstr_transit_event_counter_N}', transform=ax.transAxes, fontsize=14,
                   verticalalignment='top', bbox=props)

        pass

    def plots_handler__scattering(self, subplot_def):
        """
        Plot the scatterings to the fiber over time - north and south
        """

        ax = subplot_def["ax"]
        ax.plot(self.batcher['total_north_clicks_batch'], label='Total North Clicks')
        ax.plot(self.batcher['total_south_clicks_batch'], label='Total South Clicks')

        avg_north = int(np.average(self.batcher['total_north_clicks_batch']))
        avg_south = int(np.average(self.batcher['total_south_clicks_batch']))

        ax.axhline(y=avg_north, color='b', linestyle='--', linewidth=1, label=f'N_avg:{avg_north}')
        ax.axhline(y=avg_south, color='y', linestyle='--', linewidth=1, label=f'S_avg: {avg_south}')

        pass

    def plots_handler__pulse_shapes(self, subplot_def):

        ax = subplot_def["ax"]
        ax.plot(Config.PNSA_Exp_Gaussian_samples_N, label='North Gaussian')
        ax.plot(Config.PNSA_Exp_Gaussian_samples_S, label='South Gaussian')
        ax.plot(Config.PNSA_Exp_Square_samples_Early, label='Early Square')
        ax.plot(Config.PNSA_Exp_Square_samples_Late, label='Late Square')

    def plot_figures(self):

        super().plot_figures()
        pass

    def init_params_for_experiment(self):
        # define empty variables
        self.number_of_PNSA_sequences = math.ceil(self.M_window / self.sequence_len)

        # Reformatting the above variables into an object - (a) for better clarity (b) make it "experiment-agnostic"
        self.experiment = {
            "sequence_length": self.sequence_len,
            "measurement_window": self.M_window,  # [ns]
            "number_of_sequences": math.ceil(self.M_window / self.sequence_len)
        }
        self.number_of_detection_pulses_per_seq = sum((np.array(Config.det_pulse_amp_S) != 0) |
                                                      (np.array(Config.det_pulse_amp_N) != 0))
        self.number_of_SPRINT_pulses_per_seq = sum((np.array(Config.sprint_pulse_amp_S) != 0) +
                                                   (np.array(Config.sprint_pulse_amp_N) != 0))

        # TODO: Do we need these?
        self.tt_measure = []
        self.tt_S_measure = []

        self.folded_transmission = np.zeros(len(Config.PNSA_Exp_Gaussian_samples_S))
        self.folded_reflection = np.zeros(len(Config.PNSA_Exp_Gaussian_samples_S))

        self.tt_S_binning = np.zeros(self.number_of_PNSA_sequences + 1)
        self.seq_transit_events_live = np.zeros(self.number_of_PNSA_sequences)
        self.seq_transit_events_batched = np.zeros(self.number_of_PNSA_sequences)
        self.tt_S_SPRINT_events = np.zeros(self.sequence_len)
        self.tt_S_SPRINT_events_batch = np.zeros(self.sequence_len)
        self.num_of_det_reflections_per_seq_accumulated = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_det_transmissions_per_seq_accumulated = np.zeros(self.number_of_PNSA_sequences)

        num_of_seq_per_count = 50
        self.num_of_BP_counts_per_n_sequences = np.zeros(self.number_of_PNSA_sequences//num_of_seq_per_count)
        self.num_of_DP_counts_per_n_sequences = np.zeros(self.number_of_PNSA_sequences//num_of_seq_per_count)
        self.num_of_S_counts_per_n_sequences = np.zeros(self.number_of_PNSA_sequences//num_of_seq_per_count)

        self.num_of_det_reflections_per_seq_S = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_det_reflections_per_seq_N = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_det_transmissions_per_seq_S = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_det_transmissions_per_seq_N = np.zeros(self.number_of_PNSA_sequences)

        self.seq_with_data_points = []
        self.reflection_SPRINT_data = []  # Array of vectors with data on the number of reflections per SPRINT pulse in sequence.
        self.transmission_SPRINT_data = []  # Array of vectors with data on the number of transmissions per SPRINT pulse in sequence.
        self.BP_counts_SPRINT_data = []
        self.DP_counts_SPRINT_data = []

        self.reflection_SPRINT_data_without_transits = []  # Array of vectors with data on the number of reflections per SPRINT pulse in sequence.
        self.transmission_SPRINT_data_without_transits = []  # Array of vectors with data on the number of transmissions per SPRINT pulse in sequence.
        self.BP_counts_SPRINT_data_without_transits = []
        self.DP_counts_SPRINT_data_without_transits = []

        # Average of clicks per cycle
        self.avg_clicks = 0
        self.avg_north_clicks = 0
        self.avg_south_clicks = 0


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

        # If last measure was empty and this one is not ,clearly it's new data
        if len(prev_measure) == 0 and len(curr_measure) > 0:
            return True

        # Get the minimum value between new measure (tt_S_no_gaps) and last old measure
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

            # Check if run should terminate (maybe ESC was pressed)
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

            self.info(f'Repetition #{self.repetitions}. Waiting for values from stream (wait_cycle: {count}). No new data coming from detectors (exp_flag = {self.exp_flag})')

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

        self.divide_tt_to_reflection_trans_extended()

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
            Config.PNSA_Exp_Gaussian_samples_S)):])  # summing over the reflection from detection pulses of each sequence corresponding the the reflection_threshold_time

        # fold reflections and transmission
        self.fold_tt_histogram(exp_sequence_len=self.sequence_len)
        self.tt_histogram(num_of_bins=int(self.M_window/self.bin_size))

        # get the average number of photons in detection pulse
        self.avg_num_of_photons_per_pulse_S_live = self.get_avg_num_of_photons_in_seq_pulses(
            self.folded_tt_S_directional, self.pulses_location_in_seq_S, self.tt_FS_measure,
            self.Eff_from_taper_S)
        self.avg_num_of_photons_per_pulse_N_live = self.get_avg_num_of_photons_in_seq_pulses(
            self.folded_tt_N_directional, self.pulses_location_in_seq_N,
            self.tt_BP_measure + self.tt_DP_measure, self.Eff_from_taper_N)
        self.avg_num_of_photons_per_pulse_A_live = self.get_avg_num_of_photons_in_seq_pulses(
            (np.array(self.folded_tt_S_directional) + np.array(self.folded_tt_BP_timebins)
             + np.array(self.folded_tt_DP_timebins)).tolist(), self.pulses_location_in_seq_A,
            self.tt_BP_measure + self.tt_DP_measure, self.Eff_from_taper_N)
        self.avg_num_of_photons_per_pulse_BP_live = self.get_avg_num_of_photons_in_seq_pulses(
            self.folded_tt_BP_timebins, self.pulses_location_in_seq[-4:], self.tt_BP_measure, self.Eff_from_taper_S)
        self.avg_num_of_photons_per_pulse_DP_live = self.get_avg_num_of_photons_in_seq_pulses(
            self.folded_tt_DP_timebins, self.pulses_location_in_seq[-4:], self.tt_DP_measure, self.Eff_from_taper_S)
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
            self.folded_tt_BP_timebins, self.pulses_location_in_seq[-4:])
        self.max_value_per_pulse_DP_live = self.get_max_value_in_seq_pulses(
            self.folded_tt_DP_timebins, self.pulses_location_in_seq[-4:])
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

        # for vstirap


        pass

    def construct_pulses_sequence_only_from_config(self):
        """
        Construct the pulses string representing what happened in the experiment.
        It relies on these configurations

            det_pulse_amp_N = [0.28, 0, 0.28, 0, 0.28, 0]  ==> "N_N_N_"
            det_pulse_amp_S = [0, 0.095, 0, 0.095, 0, 0.06]  ==> "_S_S_S"

            sprint_pulse_amp_N = [0.06, 0.06]
            sprint_pulse_amp_S = [0, 0]

            det_pulse_amp_Early = [0, 0, 0, 0, 0, 0]
            sprint_pulse_amp_Early = [0, 0, 0, 0]

            det_pulse_amp_Late = [1, 1, 1, 1, 1, 1]
            sprint_pulse_amp_Late = [1, 1, 1, 1]

        """

        str = ''
        for i in range(0, len(Config.det_pulse_amp_N)):
            if Config.det_pulse_amp_N[i] > 0:
                str += 'N'
            elif Config.det_pulse_amp_S[i] > 0:
                str += 'S'
            else:
                str += '_'
        for i in range(0, len(Config.sprint_pulse_amp_N)):
            if Config.sprint_pulse_amp_N[i] > 0:
                str += 'n'
            elif Config.sprint_pulse_amp_S[i] > 0:
                str += 's'
            else:
                str += '_'
        return str

    def construct_pulses_sequence(self, pulse_seq_N, pulse_seq_S, pulse_len):
        """
        Given the pulses time-bins from North and South, and the Pulse Length:
        - Constructs an array of tuples of the form (Start-tim, End-Time, Direction)
        - Constructs a string representing the pulse sequence (e.g. 'N-S-N-S-N-S-s-s')
        """

        N_detection_pulses_locations = [tup + ('N',) for tup in pulse_seq_N if (tup[1] - tup[0]) < pulse_len]
        n_sprint_pulses_locations = [tup + ('n',) for tup in pulse_seq_N if (tup[1] - tup[0]) >= pulse_len]
        S_detection_pulses_locations = [tup + ('S',) for tup in pulse_seq_S if (tup[1] - tup[0]) < pulse_len]
        s_sprint_pulses_locations = [tup + ('s',) for tup in pulse_seq_S if (tup[1] - tup[0]) >= pulse_len]
        all_pulses = N_detection_pulses_locations + n_sprint_pulses_locations + S_detection_pulses_locations + s_sprint_pulses_locations
        all_pulses = sorted(all_pulses, key=lambda tup: tup[1])

        str = '-'.join(tup[2] for tup in all_pulses)

        # Iterate over all pulses - if one of them is sprint pulse and the Early AOM is open,
        # this means we're "redirecting" the photon to do an "X" measurement.
        measurement_direction = 'Z'
        for tup in all_pulses:
            # if all(x == 0 for x in Config.PNSA_Exp_Square_samples_Early[tup[0]:tup[1]]):
            if tup[2] in 'ns' and Config.PNSA_Exp_Square_samples_Early[tup[0]] > 0:
                measurement_direction = 'X'
                break

        return all_pulses, str, measurement_direction

    def experiment_analysis_subfolder(self):
        """
        This overrides the BaseExperiment method and returns a string which encapsulates
        the Detection pulses, the SPRINT pulse and the Measurement direction


        Examples: "NSNSNSn-Z" or "SNSNSNSNs-X"
        """
        experiment_pulses_str = self.construct_pulses_sequence_only_from_config()
        subfolder = f'{experiment_pulses_str}-{self.measurement_direction}'
        return subfolder

    def run(self, sequence_definitions, run_parameters):
        """
        Function for analyzing, saving and displaying data from SPRINT experiment.
        :param N: Number of maximum experiments (free throws) saved and displayed.
                 program we are looking for transits of atoms next to the toroid and record them.
        :param sequence_len: the number of sprint sequences (detection and sprint pulses combination)
        :param pre_comment: The comment added at the start of the experiment, usually consisting of unique experiment
                           parameters.
        :param transit_condition:
        :param lock_err_threshold: The maximal error in locking resonator to rb line (in ms, depends on the scan width/vpp)
        :param filter_delay: Delay of the filter window compared to original location (can be taken by comparing the
                             time of the 1st detection pulse pick location to the original sequence)
        :param reflection_threshold:
        :param reflection_threshold_time:
        :return:
        """

        self.run_parameters = run_parameters
        self.N = run_parameters['N']
        self.sequence_len = run_parameters['sequence_len'] if 'sequence_len' in run_parameters else len(Config.PNSA_Exp_Gaussian_samples_N)
        self.pre_comment = run_parameters['pre_comment']
        self.lock_err_threshold = run_parameters['lock_err_threshold']
        self.interference_error_threshold = run_parameters['interference_error_threshold']
        self.desired_k_ex = run_parameters['desired_k_ex']
        self.k_ex_err = run_parameters['k_ex_err']
        self.transit_condition = run_parameters['transit_condition']
        self.filter_delay = run_parameters['filter_delay'] if 'filter_delay' in run_parameters else [0, 0, 0]
        self.reflection_threshold = run_parameters['reflection_threshold']
        self.reflection_threshold_time = run_parameters['reflection_threshold_time']
        self.photons_per_det_pulse_threshold = run_parameters['photons_per_det_pulse_threshold']
        self.FLR_threshold = run_parameters['FLR_threshold']
        self.exp_flag = run_parameters['exp_flag'] if 'exp_flag' in run_parameters else False
        self.with_atoms = run_parameters['with_atoms']
        self.MZ_infidelity_threshold = run_parameters['MZ_infidelity_threshold']

        # transits parameters:
        self.bin_size = 1e3 # the size in time of bins for tt's histogram
        self.all_transits_batch_N = []  # DROR: why do we need this?
        self.all_transits_batch_S = []  # DROR: why do we need this?
        self.tt_transit_events_accumulated_N = np.zeros(int(self.M_window/self.bin_size))
        self.tt_transit_events_accumulated_S = np.zeros(int(self.M_window/self.bin_size))

        # Handle pre-comment - keep what we got in parameters or prompt the user for a comment
        if not self.pre_comment:
            self.pre_comment = self.prompt(title='Pre Experiment Run', msg='Add pre-run comment to measurement (click Cancel for "Ignore")')
            if self.pre_comment == None:
                self.pre_comment = 'Ignore'

        # TODO: Q: Are these just detector "names"/"Symbols"? And also used to define the number of detectors (e.g. 8)?
        self.Num_Of_dets = [1, 2, 3, 4, 5, 6, 7, 8]
        self.detectors_names = ['BP', 'BP', 'DP', 'DP', 'S', 'FS', 'FS', 'N']

        self.delay_in_detection_N = 30  # choose the correct delay in samples to the first detection pulse # TODO: 40?
        self.delay_in_detection_S = 20  # choose the correct delay in samples to the first detection pulse # TODO: 40?

        self.det_pulse_len = Config.det_pulse_len + Config.num_between_zeros

        self.sprint_pulse_len = Config.sprint_pulse_len + Config.num_between_zeros  # TODO: Q: what is num_between_zeros ?

        # Flag used to hold the status of with/without atoms - operated by user
        # NOTE:
        #   - "self.with_atoms" is a parameter passed for the experiment
        #   - "self.MOT_on" is the current status - it can be changed (on/off) by the user during the experiment - pressing "a"
        #     this will be reset on every run (even if user turned it off in prev run)
        self.switch_atom_no_atoms = "atoms" if self.with_atoms else "!atoms"
        self.MOT_on = self.with_atoms

        self.warm_up_cycles = 5  # Should also come as params, with some default. E.g. self.warm_up_cycles = 10 if 'warm_up_cycles' not in params else params['warm_up_cycles']

        # initialize parameters - set zeros vectors and empty lists
        self.init_params_for_experiment()

        # ----------------------------------------------------------
        # Prepare pulses location of South, North and Ancilla
        # ----------------------------------------------------------

        # TODO: Q: we used to rely on this: "int(Config.num_between_zeros/2)" - why aren't we anymore?
        self.pulses_location_in_seq, self.filter_gen = self.get_pulses_location_in_seq(0,
                                                                                       Config.PNSA_Exp_Gaussian_samples_General,
                                                                                       smearing=0)  # smearing=int(Config.num_between_zeros/2))
        self.pulses_location_in_seq_S, self.filter_S = self.get_pulses_location_in_seq(self.filter_delay[0],
                                                                                       Config.PNSA_Exp_Gaussian_samples_S,
                                                                                       smearing=0)  # smearing=int(Config.num_between_zeros/2))
        # TODO: Q: Why are we fixing the Gaussian here?
        self.PNSA_Exp_Gaussian_samples_N = Config.PNSA_Exp_Gaussian_samples_N
        # TODO: what the hell did I use it for?!
        # self.PNSA_Exp_Gaussian_samples_N[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]] = \
        #     (np.array(Config.PNSA_Exp_Gaussian_samples_N[
        #               self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]]) + \
        #      np.array(Config.PNSA_Exp_Gaussian_samples_General[
        #               self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]]) * \
        #      Config.sprint_pulse_amp_Early[0]).tolist()
        self.pulses_location_in_seq_N, self.filter_N = self.get_pulses_location_in_seq(self.filter_delay[1],
                                                                                       self.PNSA_Exp_Gaussian_samples_N,
                                                                                       smearing=0)  # smearing=int(Config.num_between_zeros/2))
        self.pulses_location_in_seq_A, self.filter_A = self.get_pulses_location_in_seq(self.filter_delay[2],
                                                                                       Config.PNSA_Exp_Gaussian_samples_Ancilla,
                                                                                       smearing=0)  # smearing=int(Config.num_between_zeros/2))

        # Construct the pulses sequence and representing string for it (e.g. 'N-S-N-S-N-S-s-s')
        self.sorted_pulses, self.experiment_type, self.measurement_direction = self.construct_pulses_sequence(self.pulses_location_in_seq_N, self.pulses_location_in_seq_S, Config.sprint_pulse_len)

        # Find the center indices of the pulses and concatenate them to one list - to be used to put text boxes in figures
        # TODO: this can be moved/calculated later when we're doing the figure work...
        self.Num_of_photons_txt_box_x_loc = np.concatenate(
            (self.num_of_photons_txt_box_loc(self.pulses_location_in_seq_S),
             self.num_of_photons_txt_box_loc(self.pulses_location_in_seq_N),
             self.num_of_photons_txt_box_loc(self.pulses_location_in_seq_A)))

        self.Num_of_photons_txt_box_x_loc_for_MZ_ports = self.num_of_photons_txt_box_loc(
            self.pulses_location_in_seq[-4:])
        self.num_of_detection_pulses = len(Config.det_pulse_amp_N)

        self.number_of_detection_pulses_per_seq, self.number_of_SPRINT_pulses_per_seq = self.number_of_pulses_per_seq()

        # Take the 2nd number in the last tuple (which is the time of the last pulse)
        # TODO: uncomment this for PNSA
        self.end_of_det_pulse_in_seq = self.pulses_location_in_seq[self.num_of_detection_pulses - 1][1]
        # TODO: give dor box babeten - and then go back to the commented line
        # self.end_of_det_pulse_in_seq = self.pulses_location_in_seq[self.num_of_detection_pulses - 1][1]+6 # 6 only relevant for sprint - not PNSA!

        # Start experiment flag and set MOT according to flag
        self.updateValue("Experiment_Switch", True)
        self.MOT_switch(with_atoms=self.with_atoms, update_parameters=True)

        self.logger.blue('Press ESC to stop measurement.')

        # Initialize the batcher
        self.batcher.set_batch_size(self.N)
        self.batcher.empty_all()

        # Associate the streams filled in OPX (FPGA) code with result handles
        # TODO: (a) why are we doing this every "run" - can we do it per sequence?
        # TODO: (b) replace the self.streams with self.bdstreams.stream_data('name')
        #self.get_handles_from_OPX_server()

        # Initialize thresholds and flags before we start
        self.sum_for_threshold = self.reflection_threshold
        self.acquisition_flag = True
        self.threshold_flag = True
        self.pause_flag = False
        self.Phase_Correction_min_diff = []

        # Prepare the sub figures needed for the plotting phase
        self.prepare_figures()

        # Initialize the variables before loop starts
        self.counter = 1  # Total number of successful cycles
        self.repetitions = 1  # Total number of cycles
        self.potential_data = 0
        self.acquisition_flag = True
        self.warm_up = True
        self.post_warm_up_completed = False  # Did we finish the post warm-up process

        self.pause_flag = False
        self.runs_status = TerminationReason.SUCCESS  # default status of the run is Success

        ############################################ START WHILE LOOP #################################################

        # ---------------------------------------------------------------------------------
        # Experiment Mainloop
        # - We iterate few warm-up rounds
        # - We iterate until we got N successful runs, or user terminated
        # ---------------------------------------------------------------------------------
        while True:

            # Check if user maybe terminated the experiment
            if self.runs_status == TerminationReason.USER:
                break

            # Set lock error and kappa_ex as coming over communication line
            self.lock_err = 99.9
            self.k_i = 12  # [MHz
            self.k_ex = 99.9  # Dummy default value for a case we don't have it
            self.interference_error = 99.9
            if 'results' in self.bdstreams.streams['Lock_Error']:
                self.lock_err = self.bdstreams.streams['Lock_Error']['results']
                self.k_ex = self.bdstreams.streams['Kappa_Ex']['results']
                self.interference_error = self.bdstreams.streams['Interference']['results']

            # Define efficiencies:
            self.Cavity_transmission = Utils.cavity_transmission(0, self.k_ex, k_i=self.k_i, h=1)
            self.Eff_from_taper_N = Config.Eff_from_taper_N * \
                                    (self.Cavity_transmission / 0.5)
            self.Eff_from_taper_S = Config.Eff_from_taper_S * \
                                    (self.Cavity_transmission / 0.5)

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

            # Perform all analytics and calculations needed for display
            self.experiment_calculations()

            # Determine if we're "acquired" - e.g., worthy to save data :-)
            self.acquisition_flag = self.is_acquired()
            if self.acquisition_flag and not self.pause_flag:

                if self.counter < self.N+1:
                    self.counter += 1

                self.logger.info('Sum of reflections: %d' % self.sum_for_threshold)
                self.num_of_det_reflections_per_seq_accumulated += self.num_of_det_reflections_per_seq_S \
                                                                   + self.num_of_det_reflections_per_seq_N
                self.num_of_det_transmissions_per_seq_accumulated += self.num_of_det_transmissions_per_seq_S \
                                                                     + self.num_of_det_transmissions_per_seq_N

                self.seq_transit_events_live = np.zeros(self.number_of_PNSA_sequences)

                ### Find transits and build histogram:  ###
                self.tt_transit_events_N, self.tt_transit_events_accumulated_N, self.all_transits_batch_N = \
                    self.find_transit_events_transitexp(self.transit_condition, self.tt_N_measure_Total, self.all_transits_batch_N,
                                                        self.tt_transit_events_accumulated_N)

                self.tt_transit_events_S, self.tt_transit_events_accumulated_S, self.all_transits_batch_S = \
                    self.find_transit_events_transitexp(self.transit_condition, self.tt_S_measure_Total, self.all_transits_batch_S,
                                                        self.tt_transit_events_accumulated_S)


                # self.find_transit_events(cond=self.transit_condition, minimum_number_of_seq_detected=2)
                # self.seq_transit_events_live[[vec for elem in self.all_transits_seq_indx for vec in elem]] += 1
                # self.seq_transit_events_batched[[vec for elem in self.all_transits_seq_indx for vec in elem]] += 1
                self.Tot_transit_events = self.tt_transit_events_S +self.tt_transit_events_N
                self.number_of_transits_live = sum(self.Tot_transit_events)
                self.number_of_transits_total = sum([vec for lst in self.batcher['Tot_transit_events_batch'] for vec in lst])

                # print(self.potential_data)

                # get the average number of photons in detection pulse
                self.avg_num_of_photons_per_pulse_S = self.get_avg_num_of_photons_in_seq_pulses(
                    self.folded_tt_S_directional_cumulative_avg, self.pulses_location_in_seq_S, [],
                    self.Eff_from_taper_S)
                self.avg_num_of_photons_per_pulse_N = self.get_avg_num_of_photons_in_seq_pulses(
                    self.folded_tt_N_directional_cumulative_avg, self.pulses_location_in_seq_N,
                    [], self.Eff_from_taper_N)
                self.avg_num_of_photons_per_pulse_A = self.get_avg_num_of_photons_in_seq_pulses(
                    (np.array(self.folded_tt_S_directional_cumulative_avg) + np.array(self.folded_tt_BP_timebins_cumulative_avg)
                     + np.array(self.folded_tt_DP_timebins_cumulative_avg)).tolist(), self.pulses_location_in_seq_A,
                    [], self.Eff_from_taper_N)
                self.avg_num_of_photons_per_pulse_BP = self.get_avg_num_of_photons_in_seq_pulses(
                    self.folded_tt_BP_timebins_cumulative_avg, self.pulses_location_in_seq[-4:], [],
                    self.Eff_from_taper_S)
                self.avg_num_of_photons_per_pulse_DP = self.get_avg_num_of_photons_in_seq_pulses(
                    self.folded_tt_DP_timebins_cumulative_avg, self.pulses_location_in_seq[-4:], [],
                    self.Eff_from_taper_S)
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
                    self.folded_tt_BP_timebins_cumulative_avg, self.pulses_location_in_seq[-4:])
                self.max_value_per_pulse_DP = self.get_max_value_in_seq_pulses(
                    self.folded_tt_DP_timebins_cumulative_avg, self.pulses_location_in_seq[-4:])
                self.Num_of_photons_txt_box_y_loc = self.max_value_per_pulse_S + self.max_value_per_pulse_N + \
                                                    self.max_value_per_pulse_A
                self.Num_of_photons_txt_box_y_loc_MZ = [[x] + [y] for x, y in zip(self.max_value_per_pulse_BP,
                                                                                  self.max_value_per_pulse_DP)]
                # self.Num_of_photons_txt_box_y_loc_MZ = self.max_value_per_pulse_BP + self.max_value_per_pulse_DP

                # Batch all data samples
                self.batcher.batch_all(self)

            # Plot figures
            self.plot_figures()

            # Did user request we change with/without atoms state?
            self.handle_user_atoms_on_off_switch()

            # Did we complete N successful iterations? (we check for N+1 because we started with 1)
            if self.should_terminate():
                break

            # We completed another repetition.
            self.repetitions += 1

        ############################################## END WHILE-LOOP #################################################

        if self.runs_status == TerminationReason.SUCCESS:
            self.logger.info(f'Finished {self.N} Runs, {"with" if self.with_atoms else "without"} atoms')
        elif self.runs_status == TerminationReason.USER:
            self.logger.info(f'User Terminated the measurements after {self.counter} Runs, {"with" if self.with_atoms else "without"} atoms')
        elif self.runs_status == TerminationReason.ERROR:
            self.logger.info(f'Error Terminated the measurements after {self.counter} Runs, {"with" if self.with_atoms else "without"} atoms')
        elif self.runs_status == TerminationReason.PLAYBACK_END:
            self.logger.info(f'Run Terminated - completed all recorded playback data after {self.counter} Runs')

        # Adding comment to measurement [prompt whether stopped or finished regularly]
        aft_comment = 'Ignore'
        if self.exp_flag:

            # Format the sequence-step message
            num_sequence_steps = len(sequence_definitions['sequence'])
            sequence_step_str = f'Sequence step #{sequence_definitions["sequence_step"] + 1}/{num_sequence_steps} - {sequence_definitions["sequence_name"]}'
            self.save_results_for_analysis = None  # By default, do not save for analysis

            dialog_str = ''

            # Check if this is (A) The last sequence in the last cycle or (B) User Terminated
            if sequence_definitions['last_iteration_and_last_sequence'] or self.runs_status == TerminationReason.USER:
                bdd = BDDialog()
                text_res, button = bdd.prompt(self.settings['dialogs']['save_dialog'])
                if button['button_name'] == 'Ignore':
                    dialog_str = 'Ignore. ' + text_res
                else:
                    dialog_str = text_res
                    self.save_results_for_analysis = button['folder_name']

            aft_comment = f'{sequence_step_str}. {dialog_str}'

        experiment_comment = f'Transit condition: {self.transit_condition}\nReflection threshold {self.reflection_threshold} @ {int(self.reflection_threshold_time/1e6)} ms'
        daily_experiment_comments = self.generate_experiment_summary_line(self.pre_comment, aft_comment, self.with_atoms, self.counter)

        # Save all results of experiment
        # TODO: should move invocation of this to after end of Run by manager, not from here
        self.save_experiment_results(experiment_comment, daily_experiment_comments)

        return self.runs_status

    def save_experiment_results(self, experiment_comment, daily_experiment_comments):
        """
        Responsible for saving all the results gathered in the experiment and required by analysis
        Return a flag telling if save should occur at all.
        TODO: add 'save_experiment_results' to BaseExperiment, Make 'run' method invoke it at the end
        """

        # If we're in playback mode and we don't want to keep results, exit function
        if self.playback['active'] and not self.playback['save_results']:
            return

        # Save Quad RF controllers commands
        # TODO: re-implement in QuadRF class - to get the data - and BDResults will save...
        if not self.playback['active']:
            quadrf_files_folder = os.path.join(self.bd_results.get_sequence_folder(sequence_definitions), 'meta_data')
            for qrdCtrl in self.QuadRFControllers:
                qrdCtrl.saveLinesAsCSV(path=quadrf_files_folder, file_name='QuadRF_table.csv')

        # Save all other files
        results = {
            "early_sequence": Config.PNSA_Exp_Square_samples_Early,  #Config.PNSA_Exp_Gaussian_samples_Early,
            "late_sequence": Config.PNSA_Exp_Square_samples_Late,  #Config.PNSA_Exp_Gaussian_samples_Late,
            "north_sequence": Config.PNSA_Exp_Gaussian_samples_N,
            "south_sequence": Config.PNSA_Exp_Gaussian_samples_S,
            "fs_sequence": (1-np.array(Config.PNSA_Exp_Square_samples_FS)).tolist(),  #Config.PNSA_Exp_Square_samples_FS,
            "pulses_location": self.sorted_pulses,

            "tt_measure_batch": self.batcher['tt_measure_batch'],
            "tt_N_measure_batch": self.batcher['tt_N_measure_batch'],
            "tt_S_measure_batch": self.batcher['tt_S_measure_batch'],
            "tt_FS_measure_batch": self.batcher['tt_FS_measure_batch'],
            "tt_BP_measure_batch": self.batcher['tt_BP_measure_batch'],
            "tt_DP_measure_batch": self.batcher['tt_DP_measure_batch'],

            "folded_tt_S_cumulative_avg": self.folded_tt_S_cumulative_avg,
            "folded_tt_N_cumulative_avg": self.folded_tt_N_cumulative_avg,
            "folded_tt_BP_cumulative_avg": self.folded_tt_BP_cumulative_avg,
            "folded_tt_DP_cumulative_avg": self.folded_tt_DP_cumulative_avg,
            "folded_tt_FS_cumulative_avg": self.folded_tt_FS_cumulative_avg,
            "folded_tt_S_directional_cumulative_avg": self.folded_tt_S_directional_cumulative_avg,
            "folded_tt_N_directional_cumulative_avg": self.folded_tt_N_directional_cumulative_avg,
            "folded_tt_BP_timebins_cumulative_avg": self.folded_tt_BP_timebins_cumulative_avg,
            "folded_tt_DP_timebins_cumulative_avg": self.folded_tt_DP_timebins_cumulative_avg,

            "MZ_BP_counts_balancing_batch": self.batcher['MZ_BP_counts_balancing_batch'],
            "MZ_BP_counts_balancing_check_batch": self.batcher['MZ_BP_counts_balancing_check_batch'],
            "MZ_DP_counts_balancing_batch": self.batcher['MZ_DP_counts_balancing_batch'],
            "MZ_DP_counts_balancing_check_batch": self.batcher['MZ_DP_counts_balancing_check_batch'],
            "Phase_Correction_vec_batch": self.batcher['Phase_Correction_vec_batch'],
            "Phase_Correction_min_vec_batch": self.batcher['Phase_Correction_min_vec_batch'],
            "Phase_Correction_value": self.batcher['Phase_Correction_value'],
            "MZ_S_tot_counts": self.batcher['MZ_S_tot_counts'],

            "Index_of_Sequences_with_data_points": self.batcher['seq_with_data_points_batch'],
            "Reflections_per_data_point": self.batcher['reflection_SPRINT_data_batch'],
            "Transmissions_per_data_point": self.batcher['transmission_SPRINT_data_batch'],
            "Bright_port_counts_per_data_point": self.batcher['BP_counts_SPRINT_data_batch'],
            "Dark_port_counts_per_data_point": self.batcher['DP_counts_SPRINT_data_batch'],

            "FLR_measurement": self.batcher['flr_batch'],
            "lock_error": self.batcher['lock_err_batch'],
            "k_ex": self.batcher['k_ex_batch'],
            "interference_error": self.batcher['interference_error_batch'],
            "exp_timestr": experiment_comment,

            "exp_comment": f'transit condition: {self.transit_condition}; reflection threshold: {self.reflection_threshold} @ {int(self.reflection_threshold_time / 1e6)} ms',
            "daily_experiment_comments": daily_experiment_comments,

            "experiment_config_values": self.Exp_Values,

            "run_parameters": self.run_parameters
        }

        # Save the results
        if self.playback['active']:
            resolved_path = self.playback['save_results_path']
        else:
            resolved_path = self.bd_results.get_sequence_folder(sequence_definitions)
        self.bd_results.save_results(results, resolved_path)

    def is_acquired(self):

        self.lock_err_flag = True
        self.k_ex_flag = True
        self.fluorescence_flag = True
        self.MZ_infidelity_flag = True

        playback = self.playback['active']

        # If we are not in experiment, we don't care about the rest - just say we're acquired.
        if not self.exp_flag:
            return True

        # Are we properly locked on resonance?
        if self.lock_err is None and not playback:
            self.lock_err_flag = False
            return False

        # if abs(self.interference_error) > self.interference_error_threshold and not playback:
        #     print(self.interference_error)
        #     return False

        if abs(self.lock_err) > self.lock_err_threshold and not playback:
            self.lock_err_flag = False
            return False

        # if self.k_ex is None and not playback:
        #     self.k_ex_flag = False
        #     return False

        # Are we on the right width ?
        # if not np.abs(self.k_ex-self.desired_k_ex) < self.k_ex_err and not playback:
        #     self.k_ex_flag = False
        #     return False

        # Is fluorescence strong enough? (assuming we're running with atoms)
        if self.fluorescence_average < self.FLR_threshold and self.with_atoms and not playback:
            self.fluorescence_flag = False
            return False

        # TODO: should be self.avg... What does the below do?
        # if len(experiment.avg_num_of_photons_per_pulse_live) == 0 or np.average(experiment.avg_num_of_photons_per_pulse_live) > self.photons_per_det_pulse_threshold:
        #     return False

        # Are any of the detectors latched?
        if self.latched_detectors() and not playback:
            return False

        if self.saturated_detectors() and not playback:
            return False

        # if (self.Infidelity_before > self.MZ_infidelity_threshold) or \
        #    (self.Infidelity_after > self.MZ_infidelity_threshold):
        #     self.MZ_infidelity_flag = False
        #     return False

        # threshold_flag = (self.sum_for_threshold < self.reflection_threshold) and \
        #                       (self.Infidelity_before <= self.MZ_infidelity_threshold) and \
        #                       (self.Infidelity_after <= self.MZ_infidelity_threshold)
        # if threshold_flag or not self.exp_flag:
        #     return True

        # All is well!
        return True

    def handle_warm_up_phase(self):

        if self.post_warm_up_completed:
            return False

        self.logger.info(f'Running Warmup phase #{self.warm_up_cycles} (exp_flag: {self.exp_flag})')

        # Get data from OPX streams
        self.get_results_from_streams()
        if self.runs_status != TerminationReason.SUCCESS:
            self.post_warm_up_completed = True
            return False

        # Make something out of the data we received on the streams
        self.ingest_time_tags()

        self.divide_tt_to_reflection_trans_extended()
        self.divide_BP_and_DP_counts(50)
        self.num_of_det_reflections_per_seq = self.num_of_det_reflections_per_seq_S + self.num_of_det_reflections_per_seq_N
        self.num_of_det_transmissions_per_seq = self.num_of_det_transmissions_per_seq_S + self.num_of_det_transmissions_per_seq_N

        # Summing over the reflection from detection pulses of each sequence corresponding the the reflection_threshold_time
        self.sum_for_threshold = sum(self.num_of_det_reflections_per_seq[-int(self.reflection_threshold_time // len(Config.PNSA_Exp_Gaussian_samples_S)):])

        # Check if conditions have been met for the completion of the warm-up phase
        warm_up_phase_complete = self.is_warm_up_phase_complete()

        # This is the first time we realize we are no loger in warm-up - so run post stage
        if warm_up_phase_complete:
            self.post_warm_up()
            self.post_warm_up_completed = True  # Mark the fact we are done
            self.logger.info(f'Warmup phase completed.')

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
        if self.exp_flag and self.lock_err > self.lock_err_threshold and not self.playback['active']:
            return False

        # We stay in warm-up if we're not within threshold
        # if self.exp_flag and self.sum_for_threshold > self.reflection_threshold:
        #     return False

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
                                        num_of_sprint_sequences=self.number_of_PNSA_sequences)

        # TODO: needed?
        self.tt_S_SPRINT_events = np.zeros(self.sequence_len)
        self.tt_S_SPRINT_events_batch = np.zeros(self.sequence_len)
        self.tt_Single_det_SPRINT_events = np.zeros((len(self.Num_Of_dets), self.sequence_len))
        self.tt_Single_det_SPRINT_events_batch = np.zeros((len(self.Num_Of_dets), self.sequence_len))
        self.tt_N_det_SPRINT_events_batch = np.zeros(self.sequence_len)
        self.tt_S_det_SPRINT_events_batch = np.zeros(self.sequence_len)

        # Fold time-tags: S, N, BP, DP, FS, N-directional, S-directional, BP-timebins, DP-timebins
        self.fold_tt_histogram(exp_sequence_len=self.sequence_len)

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
        self.avg_num_of_photons_per_pulse_S_live = self.get_avg_num_of_photons_in_seq_pulses(self.folded_tt_S_directional, self.pulses_location_in_seq_S, self.tt_FS_measure, self.Eff_from_taper_S)
        self.avg_num_of_photons_per_pulse_N_live = self.get_avg_num_of_photons_in_seq_pulses(self.folded_tt_N_directional, self.pulses_location_in_seq_N, self.tt_BP_measure + self.tt_DP_measure, self.Eff_from_taper_N)
        self.avg_num_of_photons_per_pulse_A_live = self.get_avg_num_of_photons_in_seq_pulses((np.array(self.folded_tt_S_directional) + np.array(self.folded_tt_BP_timebins)+ np.array(self.folded_tt_DP_timebins)).tolist(), self.pulses_location_in_seq_A, self.tt_BP_measure + self.tt_DP_measure, self.Eff_from_taper_N)
        self.avg_num_of_photons_per_pulse_BP_live = self.get_avg_num_of_photons_in_seq_pulses(self.folded_tt_BP_timebins, self.pulses_location_in_seq[-4:], self.tt_BP_measure, self.Eff_from_taper_S)
        self.avg_num_of_photons_per_pulse_DP_live = self.get_avg_num_of_photons_in_seq_pulses(self.folded_tt_DP_timebins, self.pulses_location_in_seq[-4:], self.tt_DP_measure, self.Eff_from_taper_S)
        self.avg_num_of_photons_per_pulse_live = self.avg_num_of_photons_per_pulse_S_live + self.avg_num_of_photons_per_pulse_N_live + self.avg_num_of_photons_per_pulse_A_live
        self.avg_num_of_photons_per_pulse_live_MZ = [[x]+[y] for x, y in zip(self.avg_num_of_photons_per_pulse_BP_live, self.avg_num_of_photons_per_pulse_DP_live)]

        self.avg_num_of_photons_per_pulse = self.avg_num_of_photons_per_pulse_live
        self.avg_num_of_photons_per_pulse_MZ = self.avg_num_of_photons_per_pulse_live_MZ

        # get box location on the y-axis:
        self.max_value_per_pulse_S_live = self.get_max_value_in_seq_pulses(self.folded_tt_S_directional, self.pulses_location_in_seq_S)
        self.max_value_per_pulse_N_live = self.get_max_value_in_seq_pulses(self.folded_tt_N_directional, self.pulses_location_in_seq_N)
        self.max_value_per_pulse_A_live = self.get_max_value_in_seq_pulses((np.array(self.folded_tt_S_directional) + np.array(self.folded_tt_BP_timebins) + np.array(self.folded_tt_DP_timebins)).tolist(), self.pulses_location_in_seq_A)
        self.max_value_per_pulse_BP_live = self.get_max_value_in_seq_pulses(self.folded_tt_BP_timebins, self.pulses_location_in_seq[-4:])
        self.max_value_per_pulse_DP_live = self.get_max_value_in_seq_pulses(self.folded_tt_DP_timebins, self.pulses_location_in_seq[-4:])

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

        self.tt_transit_events_N, self.tt_transit_events_accumulated_N, self.all_transits_batch_N = \
            self.find_transit_events_transitexp(self.transit_condition, self.tt_N_measure_Total, self.all_transits_batch_N,
                                                self.tt_transit_events_accumulated_N)

        self.tt_transit_events_S, self.tt_transit_events_accumulated_S, self.all_transits_batch_S = \
            self.find_transit_events_transitexp(self.transit_condition, self.tt_S_measure_Total, self.all_transits_batch_S,
                                                self.tt_transit_events_accumulated_S)

        self.Tot_transit_events = self.tt_transit_events_S + self.tt_transit_events_N
        self.number_of_transits_live = sum(self.Tot_transit_events)

        self.batcher.batch_all(self)

        self.number_of_transits_total = sum([vec for lst in self.batcher['Tot_transit_events_batch'] for vec in lst])
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
        cmnt = self.experiment_type
        if pre_comment is not None:
            cmnt = cmnt + '; ' + pre_comment + '; '
        if aft_comment is not None:
            cmnt = cmnt + '; ' + pre_comment + ' After comment: ' + aft_comment
        experiment_success = 'ignore' if 'ignore' in cmnt or not self.exp_flag else 'valid'
        full_line = f'{date_str},{time_str},{experiment_success},{with_atoms},{counter},{cmnt}'
        return full_line

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
        #'playback_files_path': r'C:\temp\refactor_debug\Experiment_results\PNSA\20240225\173049_Photon_TimeTags\Iter_1_Seq_2__With Atoms\playback',
        #'playback_files_path': r'C:\temp\playback_data\PNSA\20240312\121917_Photon_TimeTags\Iter_1_Seq_2__With Atoms\playback',

        # 14-Aug-2024
        'playback_files_path': r'C:\temp\playback_data\STIRAP\20240808\143855_Photon_TimeTags\NSNSNSnn-X\Iter_1_Seq_2__With Atoms\playback',
        # --> 'playback_files_path': r'C:\temp\playback_data\STIRAP\20240808\143855_Photon_TimeTags\NSNSNSnn-X\Iter_1_Seq_1__Without Atoms\playback',

        #'playback_files_path': r'C:\temp\playback_data\STIRAP\131814_Photon_TimeTags\Iter_1_Seq_2__With Atoms\playback',
        #'playback_files_path': r'F:\temp\Weizmann\playback data\STIRAP\100843_Photon_TimeTags\Iter_1_Seq_1__Without Atoms\playback',
        "load_config_from_playback": True,  # True,
        "old_format": False,
        "save_results": False,
        "save_results_path": 'C:\\temp\\playback_data',
        "max_iterations": 400,  # -1 if you want playback to run through entire data
        "max_files_to_load": 400,  # -1 if you want to load all available playback files in folder
        "plot": "LIVE",  # "LIVE", "LAST", "NONE"
        "delay": -1,  # -1,  # 0.5,  # In seconds. Use -1 for not playback delay
    }

    run_parameters = {
        'N': 500,  # 50,
        'transit_condition': [1000, 3],
        'pre_comment': '',  # Put 'None' or '' if you don't want pre-comment prompt
        'lock_err_threshold': 5,  # [Mhz]
        'interference_error_threshold': 2,  # [MHz]
        'desired_k_ex': 37, # [Mhz]
        'k_ex_err': 3,  # [Mhz]
        'filter_delay': [0, 0, 0],
        'reflection_threshold': 2550,
        'reflection_threshold_time': 9e6,
        'FLR_threshold': -0.01,
        'MZ_infidelity_threshold': 0.8,
        'photons_per_det_pulse_threshold': 12,
        'exp_flag': True,
        'with_atoms': True
    }
    # do sequence of runs('total cycles') while changing parameters after defined number of runs ('N')
    # The sequence_definitions params would replace parameters from run_parameters('N','with_atoms')
    sequence_definitions = {
        'total_iterations': 1,
        'delay_between_iterations': None,  # seconds
        'sequence': [
            {
                'name': 'Without Atoms',
                'parameters': {
                    'N': 40,
                    'with_atoms': False
                }
            },
            {
                'name': 'With Atoms',
                'parameters': {
                    'N': 80,  # 500
                    'with_atoms': True
                }
            },
        ]
    }

    experiment = VSTIRAPExperiment(playback_parameters=playback_parameters, save_raw_data=True)
    run_status = experiment.run_sequence(sequence_definitions, run_parameters)

    pass
