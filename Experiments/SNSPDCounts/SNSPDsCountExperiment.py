from Experiments.BaseExperiment.BaseExperiment import BaseExperiment
from Experiments.BaseExperiment.BaseExperiment import TerminationReason
from Experiments.QRAM import QRAM_Config_Experiment as Config

import time
import traceback
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

from Experiments.QuadRF.quadRFMOTController import QuadRFMOTController

import numpy as np

"""
TODO:
1) Why is original code using is_processing while QRAM code is not?
2) Align QRAM with this paradigm
3) Save basic data as binary - so we can rerun it
4) Remove TODOs from this file
5) Alt +/-
"""


class SNSPDsCountExperiment(BaseExperiment):
    def __init__(self, playback=False, save_raw_data=False):
        # Invoking BaseClass constructor. It will initiate OPX, QuadRF, BDLogger, Camera, BDResults, KeyEvents etc.
        super().__init__(playback, save_raw_data)
        pass

    def __del__(self):
        print('**** SpectrumExperiment Destructor ****')
        super(type(self), self).__del__()
        pass

    def initialize_experiment_variables(self):
        # TODO - probably remobe
        self.south_vals = []
        self.north_vals = []
        self.SPCMs_vals = []

    def initiliaze_QuadRF(self):
        if self._quadrf_skip: return

        self.QuadRFControllers = []

        qrfContr = QuadRFMOTController(initialValues={'Operation_Mode': 'Continuous', 'CH3_freq': '90MHz', 'CH3_amp': '31dbm'},
                            updateChannels=[3], debugging=False,
                            continuous=True)  # updates values on QuadRF (uploads table) #

        self.QuadRFControllers.append(qrfContr)  # updates values on QuadRF (uploads table)

    def experiment_mainloop_delay(self):
        """
        We override the base class method - as we don't need any delay at all
        """
        pass

    def wait_for_values_from_opx_streams(self):
        """
        Given the streams config and handles, wait and fetch the values from OPX
        """
        if self.streams is None:
            return

        # Wait_for_values
        for stream in self.streams.values():
            if stream['handler'] is not None:
                stream['handler'].wait_for_values(1)

        # Save streams to file
        # if self.save_raw_data:
        #     self.bdstreams.save_streams()

        pass

    def fetch_values_from_opx_streams(self):

        # TODO: insert a different check here
        first_stream = self.streams['Detector_1_Avg_Counts']['handler']
        while not first_stream.is_processing():
            pass

        # Fetch Values
        for stream in self.streams.values():
            if stream['handler'] is not None:
                stream['results'] = stream['handler'].fetch_all()
        pass

    def handle_user_events(self):
        """
        Handle cases where user pressed ESC to terminate or ALT_SPACE to pause/continue measurements
        """
        if self.keyPress == 'ESC':
            self.logger.blue('ESC pressed. Stopping measurement.')
            self.runs_status = TerminationReason.USER
            self.keyPress = None
        elif self.keyPress == '+':
            self.logger.blue('+ pressed - zooming in.')
            self.font_size += 10
            self.keyPress = None
        elif self.keyPress == '-':
            self.logger.blue('- pressed - zooming out.')
            self.font_size -= 10
            self.keyPress = None
        elif self.keyPress == 'ALT_0':
            self.logger.blue('ALT_0 pressed - reseting zoom.')
            self.keyPress = None

    def print_experiment_information(self):
        """
        We override the base class method - as we don't need to print anything
        """
        pass

    def prepare_figures(self):
        """
        Prepares the figure before the mainloop
        TODO: this should go into the BaseExperiment - for all experiments to enjoy
        """
        # Create figure and set the time and date as its title
        self.fig = plt.figure()
        current_date_time = time.strftime("%H:%M:%S (%d/%m/%Y)")
        self.fig.canvas.manager.set_window_title(current_date_time)

        # Set the font we want to use
        self.font_size = 16
        self.font = font_manager.FontProperties(family='Comic Sans MS', weight='bold', style='normal', size=self.font_size)

        self.plot_shown = False

        pass

    def _plot_figures(self):
        # Show the plot (if this is the first time)
        if not self.plot_shown:
            plt.show(block=False)
            self.plot_shown = True

        plt.clf()
        plt.plot(self.batcher['north_avg_counts'], label='North Counts: ' + self.N_counts + ' Hz')
        plt.plot(self.batcher['south_avg_counts'], label='South Counts: ' + self.S_counts + ' Hz')
        plt.plot(self.batcher['spcms_avg_counts'], label='SPCMs Counts: ' + self.SPCMs_counts + ' Hz')
        plt.title("Detectors Counts")

        self.font = font_manager.FontProperties(family='Comic Sans MS', weight='bold', style='normal', size=self.font_size)
        plt.legend(loc='upper left', prop=self.font)

        plt.pause(0.1)
        pass

    def plot_figures(self):
        try:
            self._plot_figures()
        except Exception as err:
            tb = traceback.format_exc()
            self.warn(f'Failed on plot_figures: {err}')
        pass

    def experiment_calculations(self):

        # TODO: start from 0 - so no need to append a dummy initial value
        avg_counts = [999]
        for stream in self.streams.values():
            avg_counts.append(stream['results'])

        # Add counts of South and North detectors and SPCMs
        # TODO: why are we sometimes getting here a None value? Check and ignore
        self.south_vals.append(sum(avg_counts[6] + avg_counts[7] + avg_counts[5]))
        self.north_vals.append(sum(avg_counts[1] + avg_counts[2] + avg_counts[3] + avg_counts[4] + avg_counts[8]))
        self.SPCMs_vals.append(sum(avg_counts[9] + avg_counts[10]))

        # Temp work around - use only batcher values
        self.south_avg_counts = self.south_vals[-1]
        self.north_avg_counts = self.north_vals[-1]
        self.spcm_avg_counts = self.SPCMs_vals[-1]

        # Format the N/S/SPCM counts
        # TODO: we mulitply by 2 since in the past, we were measuring every 500 ms, now: 100ms - calc it instead of hard-coded
        self.N_counts = "{:,}".format(self.north_vals[-1] * 2)
        self.S_counts = "{:,}".format(self.south_vals[-1] * 2)
        self.SPCMs_counts = "{:,}".format(self.SPCMs_vals[-1] * 2)

        pass

    def should_terminate(self):
        """
        This method decides if we should terminate mainloop or not
        TODO: move this to QRAM experiment and to BaseExperiment
        """
        return False

    def save_experiment_results(self):
        pass

    def pre_run(self, run_parameters):

        # Associate the streams filled in OPX (FPGA) code with result handles
        self.get_handles_from_OPX_server()

        # Await for values
        self.wait_for_values_from_opx_streams()
        pass

    def run(self, run_parameters):
        """
        Main executing function - performs pre-mainloop stuff and then runs the mainloop
        TODO: should eventually go into BaseExperiment - also for QRAM and Sprint
        """
        rp = run_parameters  # Set so we can use in short - "rp", instead of "run_parameters"...
        run_status = TerminationReason.SUCCESS

        self.logger.blue('Press ESC to terminate.')

        # Initialize the batcher
        self.batcher.set_batch_size(100)
        self.batcher.empty_all()

        # Associate the streams filled in OPX (FPGA) code with result handles
        self.get_handles_from_OPX_server()

        # Prepare the figures for the experiment
        self.prepare_figures()

        # ---------------------------------------------------------------------------------
        # Experiment Mainloop
        # - We iterate until user exists
        # ---------------------------------------------------------------------------------
        while True:

            # Handle cases where user terminates or pauses experiment
            self.handle_user_events()
            if self.runs_status == TerminationReason.USER:
                break

            # Experiment delay
            self.experiment_mainloop_delay()

            # Await for values from OPX
            self.fetch_values_from_opx_streams()

            # Informational printing
            self.print_experiment_information()

            # Perform all analytics and calculations needed for display
            self.experiment_calculations()

            # Plot figures
            self.plot_figures()

            # Batch the detectors count data
            self.batcher.batch_all(self)

            # Check if we should terminate
            if self.should_terminate():
                break

        self.save_experiment_results()

        return run_status

    def post_run(self, run_parameters):
        pass

if __name__ == "__main__":

    matplotlib_version = matplotlib.get_backend()
    print(f'In use: {matplotlib_version}')
    matplotlib.use("Qt5Agg")

    print('please switch SRS South and North directional detectors shutters to manual and then press enter.\n' +
          'This is important so the continuous laser beam wont be degraded by the shutters. ')
    input()

    experiment = SNSPDsCountExperiment(playback=False, save_raw_data=False)
    experiment.run(run_parameters={})
    pass