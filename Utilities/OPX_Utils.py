import numpy as np
from scipy import signal
import inspect
import json
from qm.qua import *
from Experiments.Enums.IOParameters import IOParameters as IOP
from matplotlib import pyplot as plt
from matplotlib.widgets import RadioButtons


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

    def plot_element(self, label_selected):
        element_name = label_selected
        element = self.config['elements'][element_name]

        for operation_name, operation_pulse in element['operations'].items():

            operation_pulse = self.config['pulses'][operation_pulse]
            # TODO - Need to handle cases where it is not 'single', or there are multiple waveforms
            waveform_name = operation_pulse['waveforms']['single']
            waveform = self.config['waveforms'][waveform_name]

            if waveform['type'] == 'arbitrary':
                y_values = waveform['samples']
            elif waveform['type'] == 'constant':
                y_values = np.full(shape=operation_pulse['length'], fill_value=waveform['sample'])

            pulse_length = len(y_values)

            #y_values = (signal.gaussian(500, std=(300 / 2.355)) * 0.2) * 2000
            x_values = np.linspace(1, pulse_length, pulse_length)

            #self.dots.set_ydata(y_values)
            #self.dots, = plt.plot(x_values, y_values, label=f"Pulse")

            plt.plot(x_values, y_values, label=f"{operation_name}")

            operation_desc = f'{operation_name}: Waveform: {waveform_name} - {waveform["type"]}'

        self.ax.set_title(f'{element_name}')
        self.fig.canvas.draw_idle()
        pass

    def plot_element_DEP(self, label_selected):

        element_name = label_selected
        element = self.config['elements'][element_name]

        y_values = (signal.gaussian(500, std=(300 / 2.355)) * 0.2) * 2000
        x_values = np.linspace(1, 500, 500)

        self.dots.set_ydata(y_values)

        self.dots, = plt.plot(x_values, y_values, label=f"Pulse")

        self.ax.set_title(f'{element_name}')
        self.fig.canvas.draw_idle()
        pass

    def plot_config(self, config):

        self.config = config

        y_values = (signal.gaussian(500, std=(300 / 2.355)) * 0.2) * 1000
        x_values = np.linspace(1, 500, 500)

        # Create the figure and axes
        self.fig, self.ax = plt.subplots(figsize=(8, 4))

        # Plot the dots
        self.dots, = plt.plot(x_values, y_values, label=f"Pulse")

        self.ax.set_xlabel("Time")
        self.ax.set_ylabel("Amplitude")

        self.ax.set_title("Pulse")
        # ax.set_ylim(0, 600)  # Set y-axis limits dynamically
        self.ax.grid(True)
        self.ax.legend()

        elements_names = list(config["elements"].keys())

        # Add radio buttons
        radio_ax = plt.axes([0.85, 0.2, 0.2, 0.4], facecolor='lightgoldenrodyellow')  # [Left, Bottom, W, H]
        radio_buttons = RadioButtons(radio_ax, elements_names)

        # Function to handle radio button selection
        radio_buttons.on_clicked(self.plot_element)

        plt.show()

        pass