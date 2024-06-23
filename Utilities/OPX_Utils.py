import numpy as np
from scipy import signal
import inspect
import json
from qm.qua import *
from Experiments.Enums.IOParameters import IOParameters as IOP
from matplotlib import pyplot as plt
from matplotlib.widgets import RadioButtons
import matplotlib.ticker as mticker
import matplotlib as mpl

# -------------------------------------------------------------------
# OPX Stuff - can be used from different OPX_Code in all experiments
# -------------------------------------------------------------------

def ms(val):
    return ms(val*1000)
def us(val):
    return ns(val*1000)
def ns(val):
    return val >> 2  # divide by 4 - as the OPX clock is 4ns per cycle

def Hz(val):
    return val

def KHz(val):
    return val*1e3

def MHz(val):
    return val*1e6

def GHz(val):
    return val*1e9

class OPX_Utils:

    def __init__(self):
        #raise Exception('No need to initialize OPX_Utils class. It used for static methods only')
        pass

    @staticmethod
    def assign_experiment_variables(self_obj):
        """
        Populate the parameters with default values
        """
        MAX_SIZE_OF_PARAMS_VECTOR_IN_OPX = 100

        # Declare a vector of 100 entries for all possible parameters
        params = declare(int, size=MAX_SIZE_OF_PARAMS_VECTOR_IN_OPX)

        parameters_count = 0
        self_members = inspect.getmembers(self_obj)

        # Iterate over all "self" values
        for member in self_members:
            member_name = member[0]
            if member_name.startswith('_'):
                continue

            member_value = member[1]
            member_value_type = type(member_value)
            if member_value_type != int and member_value_type != float and member_value_type != bool:
                continue

            if IOP.has(member_name):
                index = IOP.value_from_string(member_name)
                assign(params[index], member_value)
                parameters_count += 1

        # Iterate over all experiment values
        for key, value in self_obj.Exp_Values.items():
            if IOP.has(key):
                index = IOP.value_from_string(key)
                assign(params[index], value)
                parameters_count += 1

        return params

    @staticmethod
    def parameters_update(param):
        return

    @staticmethod
    def load_config(config_path):
        """
        Given a path to a json file ('.json' extenstion should be included), the configuration is loaded and returned
        """
        with open(config_path, 'r') as fp:
           config = json.load(fp)
        return config

    def plot_pulses(self):
        # TODO: implement
        pass

    # TODO: Q:
    # TODO: Need to check what this does. It seems to do nothing as it multiplies by ZERO.
    # TODO: Dor says it's a hack by Quantum Machine that is suppose to tie variables to elements
    # TODO: Need to check it.
    #
    @staticmethod
    def assign_variables_to_element(element, *variables):
        """
        Forces the given variables to be used by the given element thread. Useful as a workaround for when the compiler
        wrongly assigns variables which can cause gaps.
        To be used at the beginning of a program, will add a 16ns wait to the given element. Use an `align()` if needed.

        Example::

            >>> with program() as program_name:
            >>>     a = declare(int)
            >>>     b = declare(fixed)
            >>>     assign_variables_to_element('resonator', a, b)
            >>>     align()
            >>>     ...

        """
        _exp = variables[0]
        for variable in variables[1:]:
            _exp += variable
        wait(4 + 0 * _exp, element)

    # Function to toggle visibility of plot lines
    def toggle_visibility(self, event):
        legend_line = event.artist
        legend_line.set_visible(not legend_line.get_visible())
        self.fig.canvas.draw()

    def digital_marker(self, operation_pulse):

        waveform_name = 'Digital Marker'
        pulse_length = int(operation_pulse['length'])
        waveform_name = operation_pulse['digital_marker']
        samples = self.config['digital_waveforms'][waveform_name]['samples']

        # Create full waveform
        y_values = np.array([])
        for sample in samples:
            if sample[1] == 0:
                f = pulse_length - len(y_values)
            else:
                f = sample[1]
            y_values = np.concatenate((y_values, np.repeat(sample[0], f)))

        x_values = np.linspace(1, pulse_length, pulse_length)
        return waveform_name, pulse_length, x_values, y_values

    def waveforms(self, wf_type):
        pass

    def plot_element(self, label_selected):
        element_name = label_selected
        element = self.config['elements'][element_name]
        operations = element['operations'].items()

        # Remove all previous axis before we redraw them
        all_curr_axis = plt.gcf().get_axes()
        for axis_index in range(1, len(all_curr_axis)):
            all_curr_axis[axis_index].remove()

        self.axs = []

        i = 0
        for operation_name, operation_pulse_name in element['operations'].items():

            # Create axis and Select-Current-Axis ("sca")
            self.axs.append(plt.axes([0.1, 0.05 + (i * 0.12), 0.7, 0.08]))  # [Left, Bottom, W, H]
            plt.sca(self.axs[i])

            # Get the pulse
            operation_pulse = self.config['pulses'][operation_pulse_name]

            # Iterate to draw both digital and waveform signal
            for signal_type in ['digital_marker', 'waveforms']:

                if signal_type == 'digital_marker' and 'digital_marker' in operation_pulse:
                    waveform_name, pulse_length, x_values, y_values = self.digital_marker(operation_pulse)
                elif signal_type == 'waveforms' and 'waveforms' in operation_pulse:

                    # TODO - Need to handle cases where it is not 'single', or there are multiple waveforms
                    waveform_name = operation_pulse['waveforms']['single']
                    waveform = self.config['waveforms'][waveform_name]

                    if waveform['type'] == 'arbitrary':
                        y_values = waveform['samples']
                        pulse_length = len(y_values)
                        x_values = np.linspace(1, pulse_length, pulse_length)
                    elif waveform['type'] == 'constant':
                        y = waveform['sample']
                        y_values = [y, y, y]
                        pulse_length = operation_pulse['length']
                        x_values = [0, int(pulse_length/2), pulse_length]
                else:
                    continue

                formatter = mticker.FormatStrFormatter('%d')
                formatter2 = mpl.ticker.StrMethodFormatter('{x:,.0f}')
                self.axs[i].xaxis.set_major_formatter(formatter2)
                #self.axs[i].set_xlabel("Time")
                self.axs[i].set_ylabel("Amp")

                self.axs[i].plot(x_values, y_values, label=f"{operation_name} ({waveform_name})")

            legend = self.axs[i].legend(loc='upper right')
            i += 1

        plt.title(f'Time [ns]')

        # Move subplots to the left, to allow the checkboxes to fit on the right
        plt.subplots_adjust(left=0.05, right=0.8, hspace=0.5)

        self.fig.canvas.draw_idle()
        pass

    def plot_config(self, config):

        self.config = config

        # Create the figure and axes
        self.fig = plt.figure()

        elements_names = list(config["elements"].keys())

        # Add radio buttons
        radio_ax = plt.axes([0.83, 0.1, 0.2, 0.6], facecolor='lightgoldenrodyellow')  # [Left, Bottom, W, H]
        radio_buttons = RadioButtons(radio_ax, elements_names)
        radio_ax.set_title('Elements')

        # Function to handle radio button selection
        radio_buttons.on_clicked(self.plot_element)

        plt.show()

        pass
