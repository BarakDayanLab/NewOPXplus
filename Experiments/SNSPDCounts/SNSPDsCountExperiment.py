from Experiments.BaseExperiment.BaseExperiment import BaseExperiment
from Experiments.BaseExperiment.BaseExperiment import TerminationReason
from Experiments.QRAM import QRAM_Config_Experiment as Config

import time
import traceback
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

import numpy as np

"""
TODO:
1) Why is original code using "job = qm_ss.execute(dig)" while QRAM is not?
2) Why is original code using is_processing while QRAM code is not?
2) Implement await_for_values
3) Align QRAM with this paradigm
6) Save basic data as binary - so we can rerun it
7) Remove TODOs from this file
8) Cntrl +/-
9) QuadRF handling - why is it in the "with program() as dig"
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

    def experiment_mainloop_delay(self):
        """
        We override the base class method - as we don't need any delay at all
        """
        pass

    def await_for_values(self):

        while not avg_count1_handle.is_processing():
            pass

        avg_counts_res1 = avg_count1_handle.fetch_all()
        avg_counts_res2 = avg_count2_handle.fetch_all()
        avg_counts_res3 = avg_count3_handle.fetch_all()
        avg_counts_res4 = avg_count4_handle.fetch_all()
        avg_counts_res5 = avg_count5_handle.fetch_all()
        avg_counts_res6 = avg_count6_handle.fetch_all()
        avg_counts_res7 = avg_count7_handle.fetch_all()
        avg_counts_res8 = avg_count8_handle.fetch_all()
        avg_counts_res9 = avg_count9_handle.fetch_all()
        avg_counts_res10 = avg_count10_handle.fetch_all()

        pass

    def print_experiment_information(self):
        """
        We override the base class method - as we don't need to print anything
        """
        pass

    def prepare_figures(self):
        """
        Prepares the figure before the mainloop
        TODO: this should go into the mainloop - for all experiments to enjoy
        """
        # Create figure and set the time and date as its title
        self.fig = plt.figure()
        current_date_time = time.strftime("%H:%M:%S (%d/%m/%Y)")
        self.fig.canvas.manager.set_window_title(current_date_time)

        # Set the font we want to use
        self.font = font_manager.FontProperties(family='Comic Sans MS', weight='bold', style='normal', size=16)

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
        plt.plot(self.batcher['spcm_avg_counts'], label='SPCMs Counts: ' + self.SPCMs_counts + ' Hz')
        plt.title("Detectors Counts")
        plt.legend(loc='upper left', prop=self.font)

        # TODO: remove after everything works
        #plt.show()

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
        # Add counts of South and North detectors and SPCMs
        self.south_vals.append(sum(avg_counts_res6 + avg_counts_res7 + avg_counts_res5))
        self.north_vals.append(sum(avg_counts_res1 + avg_counts_res2 + avg_counts_res3 + avg_counts_res4 + avg_counts_res8))
        self.SPCMs_vals.append(sum(avg_counts_res9 + avg_counts_res10))

        # Format the N/S/SPCM counts
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
            self.sampling_timestamp = self.await_for_values()
            if self.runs_status != TerminationReason.SUCCESS:
                break

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

    run_parameters = {
        'N': 100,  # 50,
        'transit_condition': [2, 1, 2],
        'pre_comment': '',
        'lock_err_threshold': 2, # [Mhz]
        'desired_k_ex': 36,# [Mhz]
        'k_ex_err': 3, # [Mhz]
        'filter_delay': [0, 0, 0],
        'reflection_threshold': 2550,
        'reflection_threshold_time': 9e6,
        'FLR_threshold': -0.01,
        'MZ_infidelity_threshold': 1.12,
        'photons_per_det_pulse_threshold': 12,
        'Exp_flag': True,
        'with_atoms': True
    }
    # do sequence of runs('total cycles') while changing parameters after defined number of runs ('N')
    # The sequence_definitions params would replace parameters from run_parameters('N','with_atoms')
    sequence_definitions = {
        'total_cycles': 2,
        'delay_between_cycles': None,  # seconds
        'sequence': [
            {
                'parameters': {
                    'N': 50,
                    'with_atoms': False
                }
            },
            {
                'parameters': {
                    'N': 500,
                    'with_atoms': True
                }
            }
        ]
    }

    print(r'please switch SRS South and North directional detectors shutters to manual and then press enter.\n ' +
          'This is important so the continuous laser beam wont be degraded by the shutters. ')
    input()

    experiment = SNSPDsCountExperiment(playback=True, save_raw_data=False)
    experiment.run(run_parameters)
