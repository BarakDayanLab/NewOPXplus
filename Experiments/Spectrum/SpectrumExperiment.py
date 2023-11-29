from Experiments.BaseExperiment.BaseExperiment import BaseExperiment
from Experiments.Spectrum import Config_Experiment as Config
from Experiments.BaseExperiment.QuadRFMOTController import QuadRFMOTController

import matplotlib.pyplot as plt
import numpy as np
import time
import math
import pymsgbox
from UtilityResources.HMP4040Control import HMP4040Visa


class SpectrumExperiment(BaseExperiment):
    def __init__(self, config=Config.config):

        # Call parent class - BaseExperiment
        # Initiates OPX, Quad, Logger, Camera, etc.
        super().__init__()

        ##########################
        # EXPERIMENT PARAMETERS: #
        ##########################

        # -----------------------------------------------------------
        # Handle QuadRF
        # -----------------------------------------------------------

        self.QuadRFControllers = []
        # Note: So as not to connect again and again to QuadRF each time we update table, we now save the MOGDevic (actual QuadRF device) connected,
        # we hold this connection until update is finished, the we close the connection.
        # we do still hold the QuadRFController objects, for access to the table (read only!) when the experiment is running.
        qrfContr = QuadRFMOTController(initialValues=self.Exp_Values, updateChannels=(1, 2, 4),
                                       topticaLockWhenUpdating=False,
                                       debugging=False, continuous=False)
        self.QuadRFControllers.append(qrfContr)  # updates values on QuadRF (uploads table)
        self.QuadRFControllers.append(QuadRFMOTController(MOGdevice=qrfContr.dev,
                                                          initialValues={'Operation_Mode': 'Continuous',
                                                                         'CH3_freq': '90MHz', 'CH3_amp': '31dbm'},
                                                          updateChannels=[3], debugging=False,
                                                          continuous=False))  # updates values on QuadRF (uploads table)
        # self.QuadRFControllers.append(QuadRFFrequencyScannerController(MOGdevice = qrfContr.dev, channel=2, debugging=False))  # updates values on QuadRF (uploads table)

        self.Update_QuadRF_channels = set(
            {})  # Only update these channels on QuadRF when UpdateParameters method is called [note: this is a python set]
        qrfContr.disconnectQuadRF()
        # ---------- Finish handle QuadRF ------------

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
        #                                                                                    'PGC_prep_duration'] * 1e6))  # [mHz/nsec], If needed pgc preparation duration must be constant!!!

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

        # Spectrum Experiment:
        self.frequency_sweep_rep = 5
        self.same_frequency_rep = 100
        self.frequency_start = int(65e6)
        self.spectrum_bandwidth = int(60e6)
        self.num_of_different_frequncies = self.Exp_Values['Pulse_1_duration'] * 1e6 / \
                                           (2 * Config.frequency_sweep_duration) / \
                                           (self.frequency_sweep_rep * self.same_frequency_rep)
        self.frequency_diff = int(self.spectrum_bandwidth / self.num_of_different_frequncies)

        # run daily experiment
        # self.Stop_run_daily_experiment = False

        # Initialize OPX
        self.initialize_OPX()
        pass

    # This is overriding the method in base class
    def should_continue(self):
        cont = super().should_continue()
        return cont

    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

    def getHelmholtzCoilsState(self):
        self.logger.info('Getting Helmholtz coils state...')
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
            self.logger.error(res)
        return res

    # ## Boolean functions: ##
    def update_pgc_final_amplitude(self, final_amplitude):
        x = np.log(int(1 / (1 - final_amplitude)))
        self.logger.info('Called update_pgc_final_amplitude', final_amplitude)
        # self.update_lin_pgc_final_amplitude(final_amplitude)
        # self.update_exp_pgc_final_amplitude(final_amplitude)
        return x

    def MOT_switch(self, with_atoms):
        if with_atoms:
            self.update_io_parameter(1, 0)
        else:
            self.update_io_parameter(2, 0)

    def Linear_PGC_switch(self, Bool):
        self.update_io_parameter(2, Bool)

    def Transit_Exp_switch(self, Bool):
        self.Transit_switch = Bool
        self.update_io_parameter(3, Bool)

    def CRUS_Exp_switch(self, Bool):
        self.update_io_parameter(4, Bool)

    def Spectrum_Exp_switch(self, Bool):
        self.update_io_parameter(4, Bool)

    def QRAM_Exp_switch(self, Bool):
        self.update_io_parameter(4, Bool)

    def AntiHelmholtz_Delay_switch(self, Bool):
        self.update_io_parameter(5, Bool)

    def Max_Probe_counts_switch(self, Bool):
        self.update_io_parameter(6, Bool)

    # ## Antihelmholtz delay variable update functions: ##

    def update_AntiHelmholtz_delay(self, new_delay):  # In units of [msec]
        self.update_io_parameter(10, int(new_delay * 1e6 / 4))

    ## update N_snaps: ##

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # TODO: Are the 4 functions below called by user? Not only by code?
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
        self.logger.info(self.fountain_aom_chirp_rate)
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

    def get_avg_num_of_photons_in_det_pulse(self, det_pulse_len, sprint_sequence_delay, num_of_det_pulses,
                                            num_of_sprint_sequences):
        self.avg_num_of_photons_in_det_pulse = np.zeros(num_of_det_pulses)
        for i in range(num_of_det_pulses):
            detection_puls_ind = \
                list(np.arange((sprint_sequence_delay + i * det_pulse_len),
                               (sprint_sequence_delay + (i + 1) * det_pulse_len)))
            self.avg_num_of_photons_in_det_pulse[i] = \
                np.sum(self.tt_S_SPRINT_events[detection_puls_ind]) / num_of_sprint_sequences

    def get_handles_from_OPX_Server(self, Num_Of_dets):
        '''
         gets handles of timetags and counts from OPX
        :return:
        '''
        Counts_handle = []
        tt_handle = []

        for i in Num_Of_dets:
            Counts_handle.append(self.job.result_handles.get("Det" + str(i) + "_Counts"))
            tt_handle.append(self.job.result_handles.get("Det" + str(i) + "_Probe_TT"))
        FLR_handle = self.job.result_handles.get("FLR_measure")
        return Counts_handle, tt_handle, FLR_handle

    def get_tt_from_handles(self, Num_Of_dets, Counts_handle, tt_handle, FLR_handle):
        self.counts_res = []
        self.tt_res = []

        # handles wait for values
        for i in range(len(Num_Of_dets)):
            tt_handle[i].wait_for_values(1)
            Counts_handle[i].wait_for_values(1)
        FLR_handle.wait_for_values(1)

        # add the tt and counts to python vars
        for i in range(len(Num_Of_dets)):
            self.counts_res.append(Counts_handle[i].fetch_all())
            self.tt_res.append(tt_handle[i].fetch_all())

        self.FLR_res = -FLR_handle.fetch_all()

        # clear tt's vecs from last data
        self.tt_measure = []
        self.tt_BP_measure = []
        self.tt_DP_measure = []
        self.tt_N_measure = []
        self.tt_S_measure = []
        self.tt_FS_measure = []

        # remove zero-padding from tt-res and append into tt_measure
        for i in range(len(Num_Of_dets)):  # for different detectors
            # self.tt_measure.append([self.tt_res[i][(index * Config.vec_size): (index * Config.vec_size + counts)].tolist() for index, counts in
            #                         enumerate(self.counts_res[i])])
            self.tt_measure.append(self.tt_res[i][1:(self.tt_res[i][0])].tolist())
            self.tt_measure[i] = [elm + Config.detector_delays[i] for elm in self.tt_measure[i]
                                  if ((elm % self.M_window != 0) & (elm != 9999984) &
                                      ((elm + Config.detector_delays[
                                          i]) <= self.M_window))]  # Due to an unresolved bug in the OPD there are "ghost" readings of timetags equal to the maximum time of measuring window.
            self.tt_measure[i].sort()
        # unify detectors and windows within detectors and create vector of tt's for each direction (bright port, dark port, north, south and from FS) and sort them
        self.tt_BP_measure = sorted(sum(self.tt_measure[:2], []))  # unify detectors 1-3 and windows within detectors
        self.tt_DP_measure = sorted(sum(self.tt_measure[2:4], []))  # unify detectors 1-3 and windows within detectors
        self.tt_N_measure = sorted(sum(self.tt_measure[7:], []))  # unify detectors 1-3 and windows within detectors
        self.tt_S_measure = sorted(sum(self.tt_measure[4:5], []))  # unify detectors 6-8 and windows within detectors
        self.tt_FS_measure = sorted(sum(self.tt_measure[5:7], []))  # unify detectors 6-8 and windows within detectors
        self.tt_N_directional_measure = sorted(self.tt_N_measure + self.tt_BP_measure + self.tt_DP_measure)
        self.tt_S_directional_measure = sorted(self.tt_S_measure + self.tt_FS_measure)
        self.fix_gaps_spectrum_exp_tts()

    def fix_gaps_spectrum_exp_tts(self):
        '''
        fixes delays of 300 ns every 100 us in spectrum experiment
        :return:

        '''
        self.tt_S_no_gaps = [x-(x//(1e5 + 320))*320 for x in self.tt_S_directional_measure]
        self.tt_S_no_gaps = [x-(x//(2e6 + 124))*124 for x in self.tt_S_no_gaps]
        self.tt_S_no_gaps = [int(x) for x in self.tt_S_no_gaps]
    def latched_detectors(self):
        latched_detectors = []
        for indx, det_tt_vec in enumerate(self.tt_measure):  # for different detectors
            if not det_tt_vec:
                latched_detectors.append(indx)
        return latched_detectors

    def get_pulses_bins(self, det_pulse_len, sprint_pulse_len, num_of_det_pulses, num_of_sprint_pulses,
                        sprint_sequence_delay, num_of_sprint_sequences, num_init_zeros, num_fin_zeros,
                        num_between_zeros):
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
        '''
        dividing south and north tt's vectros into reflection and transmission vectors, by:
        1. converting them into counts vector in same time spacing as the initial pulse.
        2. taking the relevant pulses from the sequence to each vector (reflection/transmission)
        :return:
        '''
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
            if (element > int(0.6e6)) and (element < int(9.6e6)):
                tt_inseq = element % self.QRAM_sequence_len
                if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                    self.num_of_det_reflections_per_seq_S[(element - 1) // self.QRAM_sequence_len] += self.filter_S[
                        tt_inseq]
                    self.num_of_det_transmissions_per_seq_N[(element - 1) // self.QRAM_sequence_len] += self.filter_N[
                        tt_inseq]
            # else:  # The part of the SPRINT pulses in the sequence
            #     SPRINT_pulse_num = (tt_inseq - self.end_of_det_pulse_in_seq) // (sprint_pulse_len + Config.num_between_zeros)
            #     if SPRINT_pulse_num < self.number_of_SPRINT_pulses_per_seq:
            #         self.num_of_SPRINT_reflections_per_seq_S[(element-1) // self.QRAM_sequence_len][SPRINT_pulse_num] += 1
            #         self.num_of_SPRINT_transmissions_per_seq_N[(element-1) // self.QRAM_sequence_len][SPRINT_pulse_num] += 1
            #         tt_small_perturb+=[element]

        for element in self.tt_S_measure + self.tt_FS_measure:
            if (element > int(0.6e6)) and (element < int(9.6e6)):
                tt_inseq = element % self.QRAM_sequence_len
                if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                    self.num_of_det_reflections_per_seq_N[(element - 1) // self.QRAM_sequence_len] += self.filter_N[
                        tt_inseq]
                    self.num_of_det_transmissions_per_seq_S[(element - 1) // self.QRAM_sequence_len] += self.filter_S[
                        tt_inseq]
            # else:  # The part of the SPRINT pulses in the sequence
            #     SPRINT_pulse_num = (tt_inseq - self.end_of_det_pulse_in_seq) // (sprint_pulse_len + Config.num_between_zeros)
            #     if SPRINT_pulse_num < self.number_of_SPRINT_pulses_per_seq:
            #         self.num_of_SPRINT_reflections_per_seq_N[(element-1) // self.QRAM_sequence_len][SPRINT_pulse_num] += 1
            #         self.num_of_SPRINT_transmissions_per_seq_S[(element-1) // self.QRAM_sequence_len][SPRINT_pulse_num] += 1

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
                self.num_of_DP_counts_per_n_sequences[
                    (element - 1) // (self.QRAM_sequence_len * num_of_seq_per_count)] += \
                    np.ceil(Config.QRAM_Exp_Square_samples_Late[tt_inseq])

        for element in self.tt_FS_measure + self.tt_S_measure:
            tt_inseq = element % self.QRAM_sequence_len
            if tt_inseq > self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                self.num_of_S_counts_per_n_sequences[
                    (element - 1) // (self.QRAM_sequence_len * num_of_seq_per_count)] += \
                    np.ceil(Config.QRAM_Exp_Square_samples_Late[tt_inseq])

    def plot_folded_tt_histogram(self):

        plt.figure()
        plt.plot(sum(np.reshape(self.tt_histogram_N, [int(len(self.tt_histogram_N) / 9), 9])), label='tt_hist_N')
        plt.plot(sum(np.reshape(self.tt_histogram_S, [int(len(self.tt_histogram_S) / 9), 9])), label='tt_hist_S')
        plt.legend()

    def find_transit_events(self, N, transit_time_threshold, transit_counts_threshold):
        current_transit = []
        self.all_transits = []
        all_transits_aligned_first = []
        t_transit = []
        transit_histogram = []
        for index, value in enumerate(self.tt_S_binning_resonance):
            if not current_transit and value:  # if the array is empty
                current_transit.append(index)
            elif value:
                if ((index - current_transit[0]) * Config.frequency_sweep_duration*2) < transit_time_threshold:
                    current_transit.append(index)
                elif sum([self.tt_S_binning_resonance[i] for i in current_transit]) > transit_counts_threshold:
                    self.all_transits.append(current_transit)
                    all_transits_aligned_first.append([x - current_transit[0] for x in current_transit])
                    # tt_S_transit_events[tuple(current_transit)] += 1
                    for i in range(current_transit[0], current_transit[-1]):
                        # building cavit-atom spectrum with the counts of the detuned time bins between the start and
                        # finish of the transit:
                        self.Cavity_atom_spectrum[((((i + 1) * Config.frequency_sweep_duration * 2) %
                                               experiment.M_window) // (
                                                          Config.frequency_sweep_duration* 1000)) %
                                             self.spectrum_bin_number] += self.tt_S_binning_detuned[i + 1]
                    current_transit = [index]
                else:
                    # Finding if there any index that was saved to current transit and is close enough to the new index
                    t = [i for i, elem in enumerate(current_transit) if ((index - elem) * 480) < transit_time_threshold]
                    if t:
                        current_transit = current_transit[t[0]:] + [index]
                    else:
                        current_transit = [index]

        for index, value in enumerate(self.tt_S_binning_detuned):
            if index not in [vec for elem in self.all_transits for vec in elem]:
                self.Cavity_spectrum[(((index * Config.frequency_sweep_duration * 2) %
                                  experiment.M_window) // (Config.frequency_sweep_duration * 1000)) %
                                self.spectrum_bin_number] += value

        self.tt_S_transit_events[[i for i in [vec for elem in self.all_transits for vec in elem]]] += 1

        if self.all_transits:
            self.all_transits_batch = self.all_transits_batch[-(N - 1):] + [self.all_transits]
            # transit_histogram = np.zeros((max([vec for elem in all_transits_aligned_first for vec in elem])
            #                               + self.histogram_bin_size) // self.Transit_profile_bin_size)
            # t_transit = np.linspace(0, len(transit_histogram) * self.Transit_profile_bin_size, len(transit_histogram))


    def get_pulses_location_in_seq(self, delay, seq=Config.QRAM_Exp_Gaussian_samples_S,
                                   smearing=int(Config.num_between_zeros / 2)):
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
        seq_indx = np.where(seq_filter > 0)[0]
        seq_filter_with_smearing = seq_filter
        pulses_loc = []
        if seq_indx.any():
            start_indx = seq_indx[0]
            for i in range(1, len(seq_indx)):
                if seq_indx[i] - seq_indx[i - 1] > 1:
                    pulses_loc.append((start_indx - int(smearing), seq_indx[i - 1] + int(smearing)))
                    for j in range(start_indx - int(smearing), seq_indx[i - 1] + int(smearing)):
                        seq_filter_with_smearing[j] = 1
                    start_indx = seq_indx[i]
            pulses_loc.append((start_indx - int(smearing), seq_indx[-1] + int(smearing)))
            for j in range(start_indx - int(smearing), seq_indx[-1] + int(smearing)):
                seq_filter_with_smearing[j] = 1
        return pulses_loc, seq_filter_with_smearing

    def get_avg_num_of_photons_in_seq_pulses(self, seq, pulse_loc, tt_measure):
        avg_num_of_photons_in_seq_pulses = []
        try:
            real_number_of_seq = math.ceil(max(tt_measure) / len(Config.QRAM_Exp_Gaussian_samples_S))
            # self.logger.info('Real number of seq = %d' %real_number_of_seq)
        except:
            real_number_of_seq = self.number_of_QRAM_sequences
            # self.logger.info('Max number of seq')
        for t in pulse_loc:
            avg_num_of_photons_in_seq_pulses.append((sum(seq[t[0]:t[1]]) + seq[t[1]]) / (
                        real_number_of_seq * 0.167))  # Sagnac configuiration efficiency 16.7%
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

    def save_tt_to_batch(self, Num_Of_dets, N):
        for i in range(len(Num_Of_dets)):
            self.tt_measure_batch[i] = self.tt_measure_batch[i][-(N - 1):] + [self.tt_measure[i]]
        self.tt_S_measure_batch = self.tt_S_measure_batch[-(N - 1):] + [self.tt_S_measure]
        self.tt_N_measure_batch = self.tt_N_measure_batch[-(N - 1):] + [self.tt_N_measure]
        self.tt_BP_measure_batch = self.tt_BP_measure_batch[-(N - 1):] + [self.tt_BP_measure]
        self.tt_DP_measure_batch = self.tt_DP_measure_batch[-(N - 1):] + [self.tt_DP_measure]
        self.tt_FS_measure_batch = self.tt_FS_measure_batch[-(N - 1):] + [self.tt_FS_measure]

    def plot_sprint_figures(self, fig, ax, Num_Of_dets):

        ax[0].clear()
        ax[1].clear()
        ax[2].clear()
        ax[3].clear()
        ax[4].clear()
        ax[5].clear()

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        textstr_detuned = r'$Probe_N = %.3f$' % (
        (sum(self.tt_S_binning_detuned) * 1000) / (self.M_time / 2),) + '[MPhotons/sec]\n' \
                          + '$\overline{Probe}_S = %.3f$' % (
                          (np.mean([sum(x) for x in self.tt_S_binning_detuned_batch]) * 1000)
                          / (self.M_time / 2),) + '[MPhotons/sec]'
        textstr_resonance = r'$Probe_S = %.3f$' % (
        (sum(self.tt_S_binning_resonance) * 1000) / (self.M_time / 2),) + '[MPhotons/sec]\n' \
                            + '$\overline{Probe}_S = %.3f$' % (
                            (np.mean([sum(x) for x in self.tt_S_binning_resonance_batch]) * 1000)
                            / (self.M_time / 2),) + '[MPhotons/sec]'
        textstr_FLR = r'$\overline{FLR}_{MAX} = %.1f$' % (np.mean(self.FLR_measurement) * 1e5,) + r'$\times 10^{-5}$'
        textstr_No_transits = 'NO TRANSITS YET!!!'

        ax[0].plot(self.time_bins[::2], self.S_bins_res_acc, label='Counts histogram', color='b')
        # ax[1].plot([sum(time_bins[i:i+4]) for i in range(0, len(time_bins), 4)],
        #          [sum(tt_S_binning_detuned[i:i+4]) for i in range(0, len(tt_S_binning_detuned), 4)],
        #          label='Counts histogram', color='b')
        ax[0].set_title('On resonance histogram accumulated', fontweight="bold")
        ax[0].set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
        ax[0].text(0.05, 0.95, textstr_detuned, transform=ax[1].transAxes, fontsize=12,
                 verticalalignment='top', bbox=props)
        ax[1].legend(loc='upper right')

        ax[1].plot(self.time_bins[::2], self.S_bins_detuned_acc, label='Counts histogram', color='b')
        # ax[1].plot([sum(time_bins[i:i+4]) for i in range(0, len(time_bins), 4)],
        #          [sum(tt_S_binning_detuned[i:i+4]) for i in range(0, len(tt_S_binning_detuned), 4)],
        #          label='Counts histogram', color='b')
        ax[1].set_title('Detuned histogram accumulated', fontweight="bold")
        ax[1].set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
        ax[1].text(0.05, 0.95, textstr_detuned, transform=ax[1].transAxes, fontsize=12,
                 verticalalignment='top', bbox=props)
        ax[1].legend(loc='upper right')

        ax[2].plot(self.time_bins[::2], self.tt_S_binning_resonance, label='Counts histogram', color='b')
        # ax[2].plot([sum(time_bins[i:i+4]) for i in range(0, len(time_bins), 4)],
        #          [sum(tt_S_binning_resonance[i:i+4]) for i in range(0, len(tt_S_binning_resonance), 4)],
        #          label='Counts histogram', color='b')
        ax[2].set_title('On resonant counts (live)', fontweight="bold")
        ax[2].set(xlabel='Time [msec]', ylabel='Counts [Photons/usec]')
        ax[2].text(0.05, 0.95, textstr_resonance, transform=ax[2].transAxes, fontsize=12,
                 verticalalignment='top', bbox=props)
        ax[2].legend(loc='upper right')

        ax[3].plot(self.freq_bins, self.Cavity_atom_spectrum, label='Cavity-Atom Spectrum', color='k')
        ax[3].set(xlabel='Frequency [MHz]', ylabel='Counts [#Number]')
        ax[3].legend(loc='upper right')

        ax[4].plot(np.linspace(0,self.histogram_bin_size-1,self.histogram_bin_size), self.folded_tt_S_acc, label='pulses folded', color='k')
        ax[4].plot(np.linspace(0,self.histogram_bin_size-1,self.histogram_bin_size), self.folded_tt_S_acc_2, label='pulses folded_2', color='b')
        ax[4].plot(np.linspace(0,self.histogram_bin_size-1,self.histogram_bin_size), self.folded_tt_S_acc_3, label='pulses folded_2', color='g')
        ax[4].set(xlabel='Time [nsec]', ylabel='Counts [#Number]')
        ax[4].legend(loc='upper right')
        # ax[4].plot(self.freq_bins, self.Cavity_spectrum, label='Cavity Spectrum', color='k')
        # ax[4].set(xlabel='Time [msec]', ylabel='Counts [#Number]')
        # ax[4].legend(loc='upper right')

        if len(self.all_transits_batch) > 0:
            if self.all_transits:
                textstr_transit_counts = r'$N_{Transits} = %s $' % (len(self.all_transits),) + r'$[Counts]$'
            textstr_transit_event_counter = r'$N_{Transits Total} = %s $' % (
            len([vec for elem in self.all_transits_batch for vec in elem]),) + r'$[Counts]$'

            # ax[5].plot(self.t_transit, self.transit_histogram, label='Transit profile', color='b')
            # ax[5].set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
            # ax[5].text(0.05, 0.95, textstr_transit_counts, transform=ax[5].transAxes, fontsize=12,
            #          verticalalignment='top', bbox=props)

            ax[5].plot(self.time_bins, self.tt_S_transit_events, label='Transit events histogram', marker='*', color='b')
            ax[5].set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
            ax[5].text(0.05, 0.95, textstr_transit_event_counter, transform=ax[5].transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
        else:
            # ax[5].plot(self.t_transit, self.transit_histogram, label='Transit profile', color='b')
            # ax[5].set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
            # ax[5].text(0.25, 0.5, textstr_No_transits, transform=ax[5].transAxes, fontsize=24,
            #          verticalalignment='center', bbox=props)

            ax[5].plot(self.time_bins, self.tt_S_transit_events, label='Transit events histogram', marker='*', color='b')
            ax[5].set(xlabel='Time [nsec]', ylabel='Counts [Photons]')
            ax[5].text(0.25, 0.5, textstr_No_transits, transform=ax[5].transAxes, fontsize=24,
                     verticalalignment='center', bbox=props)

        # ax[1].set_ylim(0, 8)
        # ax[2].set_ylim(0, 8)

        # plt.tight_layout()
        plt.show()
        plt.pause(0.5)
        ###########################################################################################################

    def init_params_for_save_sprint(self, Num_Of_dets):
        # define empty variables
        self.tt_measure = []
        self.tt_measure_batch = [[]] * len(Num_Of_dets)
        self.tt_N_measure_batch = []
        self.tt_N_binning_batch = []
        self.tt_N_resonance_batch = []
        self.tt_N_detuned_batch = []
        self.tt_S_measure_batch = []
        self.tt_S_binning_batch = []
        self.tt_S_binning_detuned_batch = []
        self.tt_S_binning_resonance_batch = []
        self.tt_N_binning_detuned_batch = []
        self.tt_N_binning_resonance_batch = []
        self.tt_S_resonance_batch = []
        self.tt_S_detuned_batch = []
        self.tt_BP_measure_batch = []
        self.tt_DP_measure_batch = []
        self.tt_FS_measure_batch = []
        self.all_transits_batch = []
        self.transit_histogram_batch = []
        self.all_transits_aligned_first_batch = []
        self.transit_histogram_batch = []
        self.FLR_measurement = []
        self.Exp_timestr_batch = []

    def Save_SNSPDs_Spectrum_Measurement_with_tt(self, N, transit_profile_bin_size, pre_comment,
                                                 total_counts_threshold, transit_counts_threshold, FLR_threshold,
                                                 lock_err_threshold, exp_flag, with_atoms):
        """
        Function for analyzing,saving and displaying data from sprint experiment.
        :param N: Number of maximum experiments (free throws) saved and displayed.
                 program we are looking for transits of atoms next to the toroid and record them.
        :param qram_sequence_len: the number of sprint sequences (detection and sprint pulses combination)
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
        # if not preComment:
        #     preComment = pymsgbox.prompt('Add comment to measurement: ', default='', timeout=int(30e3))

        # set constant parameters for the function
        Num_Of_dets = [1, 2, 3, 4, 5, 6, 7, 8]
        self.histogram_bin_size = Config.frequency_sweep_duration * 2
        self.histogram_bin_number = self.M_time // self.histogram_bin_size # num of bins in cycle of frequency sweep
        self.time_bins = np.linspace(0, self.M_time, self.histogram_bin_number*2)
        self.spectrum_bin_number = self.spectrum_bandwidth // self.frequency_diff  # TODO: ask natan how many freq steps he uses
        self.freq_bins = np.linspace(self.frequency_start, self.frequency_start + self.spectrum_bandwidth,
                                self.spectrum_bin_number)
        self.Transit_profile_bin_size = transit_profile_bin_size
        time_threshold = int(
            self.histogram_bin_size * 0.8)  # The minimum time between two time tags to be counted for a transit. # TODO: might need a factor of 2???

        self.logger.blue('Press ESC to stop measurement.')

        self.init_params_for_save_sprint(Num_Of_dets)

        ####     get tt and counts from OPX to python   #####
        Counts_handle, tt_handle, FLR_handle = self.get_handles_from_OPX_Server(Num_Of_dets)

        self.get_handles_from_OPX_server()

        ############################# WHILE 1 - START #############################

        WARMUP_CYCLES = 3
        cycle = 0
        self.sum_for_threshold = transit_counts_threshold
        while exp_flag and (cycle < WARMUP_CYCLES or self.sum_for_threshold > transit_counts_threshold):

            if not self.should_continue():
                self.logger.blue('ESC pressed. Stopping measurement.')
                self.updateValue("Spectrum_Exp_switch", False)
                self.MOT_switch(True)
                self.update_parameters()
                break

            # -------------------------------------------
            # Deal with the time-tags and counts
            # -------------------------------------------

            # Filter/Manipulate the values we got
            self.ingest_time_tags(Num_Of_dets)

            # Check locking error, break if we are above threshold
            self.lock_err = (lock_err_threshold / 2) if exp_flag else self._read_locking_error()
            self.logger.debug(f'{self.lock_err}, {self.lock_err > lock_err_threshold}, {self.sum_for_threshold}')
            if self.lock_err > lock_err_threshold:
                break

            cycle += 1

        ############################# WHILE 1 - END #############################

        ####    end get tt and counts from OPX to python   #####

        ## record time
        timest = time.strftime("%Y%m%d-%H%M%S")
        datest = time.strftime("%Y%m%d")
        FLR_measurement = []
        Exp_timestr_batch = []
        lock_err_batch = []

        # initilaization
        self.tt_N_binning = np.zeros(self.histogram_bin_number * 2)
        self.tt_S_binning = np.zeros(self.histogram_bin_number * 2)
        self.tt_N_transit_events = np.zeros(self.histogram_bin_number*2)
        self.tt_S_transit_events = np.zeros(self.histogram_bin_number*2)
        self.tt_S_transit_events_accumulated = np.zeros(self.histogram_bin_number*2)
        self.all_transits_accumulated = np.zeros(self.histogram_bin_number*2)
        self.folded_tt_S_acc = np.zeros(self.histogram_bin_size, dtype=int)
        self.folded_tt_S_acc_2 = np.zeros(self.histogram_bin_size, dtype=int)
        self.folded_tt_S_acc_3 = np.zeros(self.histogram_bin_size, dtype=int)


        self.Cavity_atom_spectrum = np.zeros(self.spectrum_bin_number)
        self.Cavity_spectrum = np.zeros(self.spectrum_bin_number)

        # fold reflections and transmission
        for x in self.tt_N_directional_measure:
            self.tt_N_binning[(x - 1) // Config.frequency_sweep_duration] += 1 #TODO: x-1?? in spectrum its x
        self.tt_N_binning_avg = self.tt_N_binning
        for x in self.tt_S_no_gaps:
            self.tt_S_binning[(x - 1) // Config.frequency_sweep_duration] += 1
        for x in self.tt_S_no_gaps:
            if x < 2e6:
                self.folded_tt_S_acc[(x - 1) % self.histogram_bin_size] += 1
            if 2e6 < x < 4e6:
                self.folded_tt_S_acc_2[(x - 1) % self.histogram_bin_size] += 1
            if 8e6 < x < 10e6:
                self.folded_tt_S_acc_3[(x - 1) % self.histogram_bin_size] += 1
        self.tt_S_binning_avg = self.tt_S_binning

        # split the binning vector to odd and even - on and off resonance pulses
        self.tt_N_binning_resonance = [self.tt_N_binning[x] for x in range(len(self.tt_N_binning)) if x % 2]  # odd
        self.tt_N_binning_detuned = [self.tt_N_binning[x] for x in range(len(self.tt_N_binning)) if not x % 2]  # even
        self.tt_S_binning_resonance = [self.tt_S_binning[x] for x in range(len(self.tt_S_binning)) if x % 2]  # odd
        self.tt_S_binning_detuned = [self.tt_S_binning[x] for x in range(len(self.tt_S_binning)) if not x % 2]  # even

        if len(self.tt_S_measure_batch) == N:
            self.tt_N_transit_events[
                [i for i, x in enumerate(self.tt_N_binning_batch[0]) if x > transit_counts_threshold]] -= 1
            self.tt_S_transit_events[
                [i for i, x in enumerate(self.tt_S_binning_batch[0]) if x > transit_counts_threshold]] -= 1

        # put into batches
        self.tt_N_measure_batch = self.tt_N_measure_batch[-(N - 1):] + [self.tt_N_directional_measure]
        self.tt_N_binning_batch = self.tt_N_binning_batch[-(N - 1):] + [self.tt_N_binning]
        self.tt_S_measure_batch = self.tt_S_measure_batch[-(N - 1):] + [self.tt_S_no_gaps]
        self.tt_S_binning_batch = self.tt_S_binning_batch[-(N - 1):] + [self.tt_S_binning]
        self.tt_S_binning_resonance_batch = self.tt_S_binning_resonance_batch[-(N - 1):] + [self.tt_S_binning_resonance]
        self.tt_S_binning_detuned_batch = self.tt_S_binning_detuned_batch[-(N - 1):] + [self.tt_S_binning_detuned]
        self.S_bins_res_acc = np.sum(np.array(self.tt_S_binning_resonance_batch),
                                     0)  # tt_S_binning_resonance accumulated (sum over the batch)
        self.S_bins_detuned_acc = np.sum(np.array(self.tt_S_binning_detuned_batch),
                                     0)  # tt_S_binning_resonance accumulated (sum over the batch)







        FLR_measurement = FLR_measurement[-(N - 1):] + [self.FLR_res.tolist()]
        Exp_timestr_batch = Exp_timestr_batch[-(N - 1):] + [timest]
        lock_err_batch = lock_err_batch[-(N - 1):] + [self.lock_err]

        self.tt_N_transit_events[[i for i, x in enumerate(self.tt_N_binning) if x > transit_counts_threshold]] += 1
        self.tt_S_transit_events[[i for i, x in enumerate(self.tt_S_binning) if x > transit_counts_threshold]] += 1
        self.tt_S_transit_events_accumulated = self.tt_S_transit_events

        self.find_transit_events(N, transit_time_threshold=time_threshold,
                                 transit_counts_threshold=transit_counts_threshold)

        self.Counter = 1  # Total number of successful cycles
        self.repitions = 1  # Total number of cycles
        self.acquisition_flag = True
        self.threshold_flag = True
        self.pause_flag = False
        #
        # create figures template
        fig = plt.figure()
        ax1 = plt.subplot2grid((6, 2), (0, 0), colspan=1, rowspan=2)
        ax2 = plt.subplot2grid((6, 2), (0, 1), colspan=1, rowspan=2)
        ax3 = plt.subplot2grid((6, 2), (2, 0), colspan=1, rowspan=2)
        ax4 = plt.subplot2grid((6, 2), (2, 1), colspan=1, rowspan=2)
        ax5 = plt.subplot2grid((6, 2), (4, 0), colspan=1, rowspan=2)
        ax6 = plt.subplot2grid((6, 2), (4, 1), colspan=1, rowspan=2)

        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()

        ############################################ START WHILE LOOP #################################################

        while True:
            if self.keyPress == 'ESC':
                self.logger.blue('ESC pressed. Stopping measurement.')
                self.updateValue("QRAM_Exp_switch", False)
                self.MOT_switch(True)
                self.Stop_run_daily_experiment = True
                self.update_parameters()
                # Other actions can be added here
                break
            if self.keyPress == 'SPACE' and not self.pause_flag:
                self.logger.blue('SPACE pressed. Pausing measurement.')
                self.pause_flag = True
                self.keyPress = None
                # Other actions can be added here
            if self.keyPress == 'SPACE' and self.pause_flag:
                self.logger.blue('SPACE pressed. Continuing measurement.')
                self.pause_flag = False
                self.keyPress = None
                # Other actions can be added here
            self.lockingEfficiency = self.Counter / self.repitions
            self.logger.info(timest, self.Counter, 'Eff: %.2f' % self.lockingEfficiency,
                  'Flr: %.2f' % (1000 * np.average(self.FLR_res.tolist())))
            self.repitions += 1
            ######################################## PLOT!!! ###########################################################
            ax = [ax1, ax2, ax3, ax4, ax5, ax6]
            self.plot_sprint_figures(fig, ax, Num_Of_dets)
            ############################################################################################################
            if self.Counter == N:
                self.logger.info(f'finished {N} Runs, {"with" if with_atoms else "without"} atoms')
                # Other actions can be added here
                break
            count = 1
            while True:
                # record time:
                timest = time.strftime("%H%M%S")
                datest = time.strftime("%Y%m%d")

                self.get_values_from_streams()

                self.ingest_time_tags(Num_Of_dets)

                # Check if new tt's arrived:
                lenS = min(len(self.tt_S_directional_measure), len(self.tt_S_measure_batch[-1]))

                # Check if the number of same values in the new and last vector are less than 1/2 of the total number of values.
                is_new_tts_S = sum(np.array(self.tt_S_directional_measure[:lenS]) == np.array(
                    self.tt_S_measure_batch[-1][:lenS])) < lenS / 2

                self.logger.debug(count)
                count += 1
                time.sleep(0.01)
                if is_new_tts_S:
                    break
                if self.keyPress == 'ESC':
                    self.logger.blue('ESC pressed. Stopping measurement.')
                    self.updateValue("QRAM_Exp_switch", False)
                    self.MOT_switch(True)
                    self.update_parameters()
                    # Other actions can be added here
                    break
            # assaf - if x=self.M_window the index is out of range so i added 1
            self.lock_err = self._read_locking_error()

            self.sum_for_threshold = (len(self.tt_S_directional_measure) * 1000) / self.M_time

            # fold reflections and transmission
            self.tt_N_binning = np.zeros(self.histogram_bin_number*2)
            self.tt_S_binning = np.zeros(self.histogram_bin_number*2)

            for x in self.tt_N_directional_measure:
                self.tt_N_binning[(x - 1) // Config.frequency_sweep_duration] += 1  # TODO: x-1?? in spectrum its x
            self.tt_N_binning_avg = self.tt_N_binning
            for x in self.tt_S_no_gaps:
                self.tt_S_binning[(x - 1) // Config.frequency_sweep_duration] += 1
            self.tt_S_binning_avg = self.tt_S_binning
            for x in self.tt_S_no_gaps:
                # if x < (2e6 + 6400):
                self.folded_tt_S_acc[(x-1) % self.histogram_bin_size] += 1
                if (2e6 + 6400 + 120) < x < (4e6 + 6400*2 + 120):
                    self.folded_tt_S_acc_2[(x - 1) % self.histogram_bin_size] += 1
                if (8e6 + 6400*4 + 120*4) < x < (10e6):
                    self.folded_tt_S_acc_3[(x - 1) % self.histogram_bin_size] += 1

            # split the binning vector to odd and even - on and off resonance pulses
            self.tt_N_binning_resonance = [self.tt_N_binning[x] for x in range(len(self.tt_N_binning)) if x % 2]  # odd
            self.tt_N_binning_detuned = [self.tt_N_binning[x] for x in range(len(self.tt_N_binning)) if
                                         not x % 2]  # even
            self.tt_S_binning_resonance = [self.tt_S_binning[x] for x in range(len(self.tt_S_binning)) if x % 2]  # odd
            self.tt_S_binning_detuned = [self.tt_S_binning[x] for x in range(len(self.tt_S_binning)) if
                                         not x % 2]  # even

            if exp_flag:
                if (self.lock_err > lock_err_threshold) or \
                        ((1000 * np.average(self.FLR_res.tolist()) < FLR_threshold) and with_atoms) or \
                        self.latched_detectors():
                    self.acquisition_flag = False
                else:
                    self.acquisition_flag = True

            self.threshold_flag = (self.sum_for_threshold < total_counts_threshold)

            if (self.threshold_flag or not exp_flag) and self.acquisition_flag and not self.pause_flag:

                if self.Counter < N:
                    self.Counter += 1

                self.tt_N_binning = np.zeros(self.histogram_bin_number*2)
                self.tt_S_binning = np.zeros(self.histogram_bin_number*2)
                self.tt_S_transit_events = np.zeros(self.histogram_bin_number*2)

                self.tt_S_binning_avg = self.tt_S_binning_avg * (1 - 1 / self.Counter)

                for x in self.tt_N_directional_measure:
                    self.tt_N_binning[(x - 1) // self.histogram_bin_size] += 1
                for x in self.tt_S_no_gaps:
                    self.tt_S_binning[(x - 1) // self.histogram_bin_size] += 1
                    self.tt_S_binning_avg[(x - 1) // self.histogram_bin_size] += 1 / self.Counter

                self.tt_N_measure_batch = self.tt_N_measure_batch[-(N - 1):] + [self.tt_N_directional_measure]
                self.tt_N_binning_batch = self.tt_N_binning_batch[-(N - 1):] + [self.tt_N_binning]
                self.tt_N_binning_resonance_batch = self.tt_N_binning_resonance_batch[-(N - 1):] + [self.tt_N_binning_resonance]
                self.tt_N_binning_detuned_batch = self.tt_N_binning_detuned_batch[-(N - 1):] + [self.tt_N_binning_detuned]
                self.tt_S_measure_batch = self.tt_S_measure_batch[-(N - 1):] + [self.tt_S_no_gaps]
                self.tt_S_binning_batch = self.tt_S_binning_batch[-(N - 1):] + [self.tt_S_binning]
                self.tt_S_binning_resonance_batch = self.tt_S_binning_resonance_batch[-(N - 1):] + [self.tt_S_binning_resonance]
                self.tt_S_binning_detuned_batch = self.tt_S_binning_detuned_batch[-(N - 1):] + [self.tt_S_binning_detuned]
                self.S_bins_res_acc = np.sum(np.array(self.tt_S_binning_resonance_batch),0) # tt_S_binning_resonance accumulated (sum over the batch)
                self.S_bins_detuned_acc = np.sum(np.array(self.tt_S_binning_detuned_batch),
                                                 0)  # tt_S_binning_resonance accumulated (sum over the batch)

                self.tt_N_transit_events[
                    [i for i, x in enumerate(self.tt_N_binning) if x > transit_counts_threshold]] += 1
                self.tt_S_transit_events[
                    [i for i, x in enumerate(self.tt_S_binning) if x > transit_counts_threshold]] += 1
                self.tt_S_transit_events_accumulated = self.tt_S_transit_events_accumulated + self.tt_S_transit_events

                self.find_transit_events(N, transit_time_threshold=time_threshold,
                                 transit_counts_threshold=transit_counts_threshold)

                FLR_measurement = FLR_measurement[-(N - 1):] + [self.FLR_res.tolist()]
                Exp_timestr_batch = Exp_timestr_batch[-(N - 1):] + [timest]
                lock_err_batch = lock_err_batch[-(N - 1):] + [self.lock_err]
                self.save_tt_to_batch(Num_Of_dets, N)

        ############################################## END WHILE LOOP #################################################

        # For debugging:
        # Utils.most_common([x for vec in self.tt_S_measure_batch + self.tt_N_measure_batch for x in vec])

        ## Adding comment to measurement [prompt whether stopped or finished regularly]
        if exp_flag:
            if self.Counter < N:
                aftComment = pymsgbox.prompt('Add comment to measurement: ', default='', timeout=int(30e3))
            else:
                aftComment = ''
        else:
            aftComment = 'ignore'

        # if aftComment == 'Timeout': aftComment = None

        #### Handle file-names and directories #####
        ## Saving: np.savez(filedir, data = x) #note: @filedir is string of full directory; data is the queyword used to read @x from the file:
        ## Loading: file = np.load(f, allow_pickle = True)
        ##          x = file['data']

        # TODO: re-implement in QuadRF class - to get the data - and BDResults will save...
        # Save Quad RF controllers commands
        dirname = self.bd_results.get_root()
        for qrdCtrl in self.QuadRFControllers:
            qrdCtrl.saveLinesAsCSV(f'{dirname}\\QuadRF_table.csv')

        # ------------------------------------------------------------------------
        #    Save (a) Input sequences (b) Output results (c) Meta Data
        # ------------------------------------------------------------------------

        results = {
            "early_sequence": self.early_sequence,
            "late_sequence": self.late_sequence,
            "north_sequence": self.north_sequence,
            "south_sequence": self.south_sequence,
            "fs_sequence": self.fs_sequence,

            "tt_measure_batch_1": self.tt_measure_batch[0],
            "tt_measure_batch_2": self.tt_measure_batch[1],
            "tt_measure_batch_3": self.tt_measure_batch[2],
            "tt_measure_batch_4": self.tt_measure_batch[3],
            "tt_measure_batch_5": self.tt_measure_batch[4],
            "tt_measure_batch_6": self.tt_measure_batch[5],
            "tt_measure_batch_7": self.tt_measure_batch[6],
            "tt_measure_batch_8": self.tt_measure_batch[7],

            "tt_N_measure_batch": self.tt_N_measure_batch,
            "tt_S_measure_batch": self.tt_S_measure_batch,
            "tt_FS_measure_batch": self.tt_FS_measure_batch,
            "tt_BP_measure_batch": self.tt_BP_measure_batch,
            "tt_DP_measure_batch": self.tt_DP_measure_batch,

            "FLR_measurement": FLR_measurement,
            "lock_error": lock_err_batch,

            "exp_comment": f'transit condition: minimum time between reflection = {time_threshold} with at least {transit_counts_threshold} reflections ; reflection threshold: {total_counts_threshold} MCounts / sec',

            "daily_experiment_comments": self.generate_experiment_summary_line(pre_comment, aftComment, with_atoms, self.Counter),

            "max_probe_counts": "TBD",

            "experiment_config_values": self.Exp_Values
        }
        self.bd_results.save_results(results)

        self.updateValue("Spectrum_Exp_switch", False)
        self.update_parameters()

    def generate_experiment_summary_line(self, pre_comment, aftComment, with_atoms, counter):
        time_str = time.strftime("%H%M%S")
        date_str = time.strftime("%Y%m%d")
        cmnt = None
        if pre_comment is not None: cmnt = pre_comment + '; '
        if aftComment is not None: cmnt = pre_comment + aftComment
        if pre_comment is None and aftComment is None: cmnt = 'No comment.'
        experiment_success = 'ignore' if 'ignore' in cmnt else 'valid'
        full_line = f'{date_str},{time_str},{experiment_success},{with_atoms},{counter},{cmnt}'
        return full_line

        ## ------------------ end of saving section -------

    def Start_Transit_Exp_with_tt(self, N=500, Histogram_bin_size=1000, Transit_profile_bin_size=100, preComment=None,
                                  total_counts_threshold=1, transit_counts_threshold=5, FLR_threshold=0.11,
                                  lock_err_threshold=0.003, Exp_flag=True, with_atoms=True):
        self.QRAM_Exp_switch(True)
        self.MOT_switch(with_atoms)
        self.update_parameters()
        self.Save_SNSPDs_Transit_Measurement_with_tt(N, Histogram_bin_size, Transit_profile_bin_size, preComment,
                                                     total_counts_threshold, transit_counts_threshold, FLR_threshold,
                                                     lock_err_threshold, exp_flag=Exp_flag,
                                                     with_atoms=with_atoms)

    def pre_run(self, run_parameters):
        # Change the pre comment based on the with_atoms parameter
        suffix = ' with atoms' if run_parameters['with_atoms'] else ' without atoms'
        run_parameters['pre_comment'] += suffix
        pass

    def run(self, run_parameters):
        rp = run_parameters  # Set so we can use in short - "rp", instead of "run_parameters"...

        # Set switches
        self.Spectrum_Exp_switch(True)
        self.MOT_switch(rp['with_atoms'])

        # TODO: Q: This is the first time we call this, it then calls save_config_table, but there's still no experiment, so it cannot save the config anywhere...
        # TODO: this needs fixing
        self.update_parameters()

        # TODO: Q: Config.QRAM_Exp_Gaussian_samples_S is constructed in a function, using the parameter "sprint_pulse_len" - so why not use it here?
        # TODO: Q: (a) we don't want to use duplicate variables holding the same value, (b) it mentions "samples_S" - but it's the same for "N" as well...
        self.Save_SNSPDs_Spectrum_Measurement_with_tt(N=rp['N'],
                                                      transit_profile_bin_size=rp['Transit_profile_bin_size'],
                                                      pre_comment=rp['preComment'],
                                                      total_counts_threshold=rp['total_counts_threshold'],
                                                      transit_counts_threshold=rp['transit_counts_threshold'],
                                                      FLR_threshold=rp['FLR_threshold'],
                                                      lock_err_threshold=rp['lock_err_threshold'],
                                                      exp_flag=rp['Exp_flag'],
                                                      with_atoms=rp['with_atoms'])

    def post_run(self, run_parameters):
        pass

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


if __name__ == "__main__":

    run_parameters = {
        'histogram_bin_size': 1000,
        'transit_profile_bin_size': 10,
        'pre_comment': 'Transit experiment with 2uW of 731.30nm light coupled',
        'total_counts_threshold': 1,
        'transit_counts_threshold': 5,
        'FLR_threshold': 0.08,
        'lock_err_threshold': 0.004,
        'Exp_flag': False,
        # TODO: we need to read the below from the results folder, not put it here HARD CODED
        'experiment_folder_path': r'U:\Lab_2023\Experiment_results\Spectrum'
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

    experiment = SpectrumExperiment()

    # TODO: REMOVE, for debug only
    sequence_definitions = None

    if sequence_definitions is None:
        experiment.run(run_parameters)
    else:
        experiment.run_sequence(sequence_definitions, run_parameters)
