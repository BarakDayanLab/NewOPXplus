from Experiments.BaseExperiment.BaseExperiment import BaseExperiment
from Experiments.BaseExperiment.BaseExperiment import TerminationReason
from Experiments.SNSPDCounts import SNSPDCounts_Config_Experiment as Config


import time
import traceback
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

from Experiments.QuadRF.quadRFMOTController import QuadRFMOTController


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
        pass

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
        elif self.keyPress == '+' or self.keyPress == '=':
            self.font_size += 10
            self.line_width += 1
            self.keyPress = None
        elif self.keyPress == '-':
            self.font_size -= 10
            self.line_width -= 1
            self.keyPress = None
        else:
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

        # Set line_width of plots
        self.line_width = 2

        self.plot_shown = False

        pass

    def _plot_figures(self):
        # Show the plot (if this is the first time)
        if not self.plot_shown:
            plt.show(block=False)
            self.plot_shown = True

        plt.clf()
        plt.plot(self.batcher['north_avg_counts'], linewidth=self.line_width, label='North Counts: ' + self.N_counts + ' Hz')
        plt.plot(self.batcher['south_avg_counts'], linewidth=self.line_width, label='South Counts: ' + self.S_counts + ' Hz')
        plt.plot(self.batcher['spcms_avg_counts'], linewidth=self.line_width, label='SPCMs Counts: ' + self.SPCMs_counts + ' Hz')
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

        avg_counts = []
        for stream in self.streams.values():
            if stream['results'] is None:
                avg_counts.append(0)  # Just for safety. Sometimes opx code did not yet send values, so we dont want to crash
            else:
                avg_counts.append(stream['results'])

        # Add counts of South and North detectors and SPCMs
        self.south_avg_counts = sum(avg_counts[5] + avg_counts[6] + avg_counts[4])
        self.north_avg_counts = sum(avg_counts[0] + avg_counts[1] + avg_counts[2] + avg_counts[3] + avg_counts[7])
        self.spcm_avg_counts = sum(avg_counts[8] + avg_counts[9])

        factor = int(1000 / Config.Measure_Time)

        # Format the N/S/SPCM counts. Multiply by factor to get to Hz
        self.N_counts = "{:,}".format(self.south_avg_counts * factor)
        self.S_counts = "{:,}".format(self.north_avg_counts * factor)
        self.SPCMs_counts = "{:,}".format(self.spcm_avg_counts * factor)

        pass

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

    print('\033[96m' + 'NOTE!\nSwitch SRS South and North directional detectors shutters to MANUAL and then press <Enter>.\n' +
          '(this is important as continuous laser beam tends to degrade by the shutters)' + '\033[0m')
    input()

    experiment = SNSPDsCountExperiment(playback=False, save_raw_data=False)
    experiment.run(run_parameters={})
    pass