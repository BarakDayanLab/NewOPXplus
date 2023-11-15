import os.path
import sys
import time
import json
import socket
import pathlib
import matplotlib.pyplot as plt

from Utilities.Utils import Utils
from Utilities.BDLogger import BDLogger
from Experiments.BaseExperiment import Config_Table
from Experiments.BaseExperiment.Config_Table import Initial_Values, Values_Factor
from Experiments.BaseExperiment.QuadRFMOTController import QuadRFMOTController

import logging
from logging import StreamHandler, Formatter, INFO, WARN, ERROR
from qm.QuantumMachinesManager import QuantumMachinesManager
from qm.qua import *
import pymsgbox
from pynput import keyboard


class BaseExperiment:

    def __init__(self, opx_definitions):

        # Setup console logger. We do this first, so rest of code can use logging functions.
        self.logger = BDLogger()

        # This should be set later by run parameters. If not, no error file can be read
        self.lock_error_file = None

        # Set paths map
        self.set_paths_map()

        # Debugging plot
        self.dbg_plot = None

        self.Exp_Values = Initial_Values  # Initialize experiment values to be as in Config_Table.py

        self.opx_definitions = opx_definitions

        # Set mainloop flags
        self.halt_experiment = False
        self.ignore_data = False

        # Open server-socket
        # TODO: complete implementation of this
        #self._open_socket()

        # Attempt to initialize Camera functionality
        self.connect_camera()

        # TODO: What is this?
        # Setup of MW spectroscopy
        self.mws_logger = logging.getLogger("MWSpectroscopy")
        self.mws_logger.setLevel(INFO)
        handler = StreamHandler(sys.stdout)
        handler.setFormatter(Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
        self.mws_logger.addHandler(handler)

        # Setup keyboard listener
        self.listener = keyboard.Listener(on_press=self._on_press, on_release=self._on_release)
        self.listener.start()
        self.keyPress = None

        # TODO: Listen to events - messages / keyboard

        pass

    def __del__(self):
        # Close all OPX related instances
        #if hasattr(self,'job'):
        # self.job.halt()
        if hasattr(self, 'qm'):
            self.qm.close()
        if hasattr(self, 'qmm'):
            self.qmm.close()
        pass

    # Initialize the OPX
    def initialize_OPX(self):
        try:
            self.qmm = QuantumMachinesManager(host=self.opx_definitions['connection']['host'], port=self.opx_definitions['connection']['port'])
            self.qmm.clear_all_job_results()
            self.qmm.reset_data_processing()
            self.qmm.close_all_quantum_machines()
            self.qm = self.qmm.open_qm(self.opx_definitions['config'])
            self.job = self.opx_definitions['control'](self, self.qm)
        except Exception as err:
            self.logger.warn(f'Unable to connect to OPX. {err}')
            pass
        except KeyboardInterrupt:
            self.job.halt()
            self.qmm.reset_data_processing()

        self.io1_list = []
        self.io2_list = []

    def set_paths_map(self):
        self.paths_map = {
            "cwd": os.getcwd(),
            "root": str(pathlib.Path(__file__).parent.resolve())
        }
        pass

    # Determines whether mainloop should continue (e.g. ESC was pressed)
    def should_continue(self):
        return not self.halt_experiment

    # Returns the error of the locking mechanism of the resonator to Rb line
    def _read_locking_error(self):
        if self.lock_error_file == None:
            self.logger.warn('Lock error file not defined. Not reading error signal.')
            return None

        max_time = 100 * 60  # The delay is 10 milli-seconds, so this is one minute
        lock_err = None
        ticks = 0
        while lock_err == None:
            try:
                data = np.load(self.lock_error_file, allow_pickle=True)
                lock_err = np.abs(data)
                break
            except Exception as err:
                self.logger.error(f'Error in loading lock-error file. {err}')
            time.sleep(0.01)  # in seconds
            ticks = ticks + 1
            if ticks > max_time:
                break
        return lock_err

    #--------------------------
    # https://stackoverflow.com/questions/23828264/how-to-make-a-simple-multithreaded-socket-server-in-python-that-remembers-client
    #--------------------------
    def _open_socket(self):
        # Create an INET, STREAMing socket
        serversocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        # Bind the socket to a public host, and a well-known port
        serversocket.bind((socket.gethostname(), 80))
        # Become a server socket
        serversocket.listen(5)

    # -------------------------------------------------------
    # Keyboard related methods
    # -------------------------------------------------------

    def _on_press(self, key):
        self.alt_modifier = ''
        if key == keyboard.Key.ctrl_l or key == keyboard.Key.ctrl_r:
            self.alt_modifier = self.alt_modifier + 'CTRL_'
        if key == keyboard.Key.alt or key == keyboard.Key.alt_l or key == keyboard.Key.alt_r or key == keyboard.Key.alt_gr:
            self.alt_modifier = self.alt_modifier + 'ALT_'
        if key == keyboard.Key.shift:
            self.alt_modifier = self.alt_modifier + 'SHIFT_'
        if len(self.alt_modifier) > 0:
            self.alt_modifier = self.alt_modifier[0:-1]
        else:
            self.alt_modifier = None

    def _is_modifier_key(self, key):
        if key == keyboard.Key.shift or \
                key == keyboard.Key.ctrl_l or key == keyboard.Key.ctrl_r or \
                key == keyboard.Key.alt or key == keyboard.Key.alt_l or key == keyboard.Key.alt_r or key == keyboard.Key.alt_gr:
            return True
        return False

    def _on_release(self, key):
        self.halt_experiment = False
        self.ignore_data = False

        if str(key) == "'q'" or str(key) == "'Q'":
            self.halt_experiment = True

        if key == keyboard.Key.esc and self.alt_modifier == 'SHIFT':
            self.halt_experiment = True

        # TODO: Generalize it to be in the form of "ALT-SPACE" or "CTRL-SPACE" or "SHIFT-ESC"
        if key == keyboard.Key.space and self.alt_modifier == 'ALT':
            self.ignore_data = True
            # print('Alt-Space released')
            self.keyPress = 'ALT_SPACE'

        if not self._is_modifier_key(key):
            self.alt_modifier = False

    # -------------------------------------------------------
    # Logger related methods
    # -------------------------------------------------------

    def blue(self, msg):
        self.logger.blue(msg)

    def error(self, msg):
        self.logger.error(msg)

    def warn(self, msg):
        self.logger.warn(msg)

    def info(self, msg):
        self.logger.info(msg)

    def debug(self, msg):
        self.logger.debug(msg)

    def prompt(self, msg, default='', timeout=30000):
        res = pymsgbox.prompt(msg, default=default, timeout=timeout)
        return res

    # -------------------------------------------------------
    # Camera connection methods
    # -------------------------------------------------------
    def connect_camera(self):
        """
        Attempt to connect to camera. Return if connected.
        :return: True - if connected, False - failed to connect
        """
        try:
            from UtilityResources import MvCameraController
            self.camera = MvCameraController.MvCameraController()
            return True
        except Exception as e:
            self.warn(f'Could not connect to camera ({e})')
            return False

    def disconnect_camera(self):
        self.camera = None

    def is_camera_connected(self):
        return self.camera is not None

    def is_camera_disconnected(self):
        return self.camera is None

    # -------------------------------------------------------
    # Experiment Management related methods
    # -------------------------------------------------------

    def pre_run(self, run_parameters):
        raise Exception("'pre_run' method must be implemented in subclass!")
        pass

    def run(self, run_parameters):
        raise Exception("'run' method must implemented in subclass!")

    def post_run(self, run_parameters):
        raise Exception("'post_run' method must be implemented in subclass!")

    def run_sequence(self, sequence_definitions, run_parameters):

        num_cycles = sequence_definitions['total_cycles']
        for i in range(num_cycles):
            self.info(f'Starting cycle sequence # {i}')

            for run in sequence_definitions['sequence']:
                # Modify params
                merged_params = Utils.merge_jsons(run_parameters, run['parameters'])

                self.pre_run(merged_params)
                self.run(merged_params)
                self.post_run(merged_params)

            if sequence_definitions['delay_between_cycles']:
                self.info(f'Sleeping')
                time.sleep(sequence_definitions['delay_between_cycles'])

        self.info(f'Completed {num_cycles} cycles!')

    #----------------------------------------
    # OPX Stuff
    #----------------------------------------

    # TODO: Q:
    # TODO: Need to check what this does. It seems to do nothing as it multiplies by ZERO.
    # TODO: Dor says it's a hack by Quantum Machine that is suppose to tie variables to elements
    # TODO: Need to check it.
    #
    # TODO: Move this to a different file: opx_utils.... ?
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

    """
        There are two vectors we get from the OPX/OPD: Counts and Time-Tags
        For each "cycle", we get two outputs:

        +---------------+------------------+
        +  Counts       |  Time-Tags Array |
        +---------------+------------------+

        Counts:
        +------+------+------+------+------+------+
        +  12  |  9   |   1  |  23  |  54  | ...  +
        +------+------+------+------+------+------+

        Timetags:
        +------+------+------+------+------+------+
        +  23  |  79  |  91  | 1023 | 2354 | ...  +
        +------+------+------+------+------+------+

        For example: the photon at position 0 was received at time 23ns on a time-window of 10e7 nano seconds
        The size of this array is dynamic - it holds values as the number of counts received

        NOTE: we have today an issue that the stream may still hold values from the previous detections        
    """

    def get_handles_from_OPX_server(self):
        '''
        Given the streams config, connect them to the OPX handles
        '''
        if self.streams is None:
            return
        self.streams = self.opx_definitions['streams']
        for key, value in self.streams.items():
            value['handler'] = self.job.result_handles.get(key)

    def get_values_from_streams(self):
        '''
        Given the streams config and handles, wait and fetch the values from OPX
        '''
        if self.streams is None:
            return
        # TODO: Should we first "wait_for_values" on all and only then "fetch_all"?
        for stream in self.streams.values():
            stream['handler'].wait_for_values(1)

        for stream in self.streams.values():
            stream['results'] = stream['handler'].fetch_all()

        pass

    #------------------------------------------------------------------------
    # Updates a specific value - either related to (a) Operation Mode (b) QuadRF or (c) OPX
    #
    # Usage:
    #   - self.updateValue("Operation_Mode", "QRAM_Exp")
    #   - self.updateValue("QRAM_Exp_switch", False)
    #------------------------------------------------------------------------
    def updateValue(self, key, value, update_parameters=False):
        if key in self.Exp_Values:
            self.Exp_Values[key] = value
            #------------------------------------------------
            # Special cases
            # - 'Operation_Mode' - we're switching the entire operation mode and its config tree
            #------------------------------------------------
            if key == 'Operation_Mode':
                if value not in Config_Table.Operation_Modes:
                    self.logger.warn(f'Warning: {key} is not a known operation mode!')
                    # TODO: ok, so now what? We just coninue? We need to raise an exception and abort things - or handle by caller!
                    # TODO: Shouldn't this be a fatal message? (not Warn)
                # Iterates over all key-values in the specific operation mode and sets them
                for k in Config_Table.Operation_Modes[value]:
                    self.updateValue(k, Config_Table.Operation_Modes[value][k])
            elif key == 'MOT_duration':  # This is a patch. It should work, but still.
                # TODO: Not sure I understood the above comment.... ?
                MOT_rep = Config_Table.updateMOT_rep()
                self.updateValue('MOT_rep', MOT_rep)
            elif key == 'OD_FS_pulse_duration':
                self.updateValue('Pulse_1_duration',
                                 self.OD_FS_Start + self.OD_FS_pulses_spacing + 2 * value + 2 * self.OD_FS_sleep)
                self.OD_FS_pulse_duration = value
            elif key == 'OD_FS_pulses_spacing':
                self.updateValue('Pulse_1_duration',
                                 self.OD_FS_Start + value + 2 * self.OD_FS_pulse_duration + 2 * self.OD_FS_sleep)
                self.OD_FS_pulses_spacing = value
            elif key == 'OD_FS_Start':
                self.updateValue('Pulse_1_duration',
                                 value + self.OD_FS_pulses_spacing + 2 * self.OD_FS_pulse_duration + 2 * self.OD_FS_sleep)
                self.OD_FS_Start = value
            elif key == 'Depump_pulse_duration':
                self.updateValue('Pulse_1_duration', self.Depump_Start + self.Depump_pulses_spacing + 2 * value)
                self.Depump_pulse_duration = value
            elif key == 'Depump_pulses_spacing':
                self.updateValue('Pulse_1_duration', self.Depump_Start + value + 2 * self.Depump_pulse_duration)
                self.Depump_pulses_spacing = value
            elif key == 'Depump_Start':
                self.updateValue('Pulse_1_duration',
                                 value + self.Depump_pulses_spacing + 2 * self.Depump_pulse_duration)
                self.Depump_Start = value

            # QuadRF - Keep track of which channels are updated
            if key in Config_Table.Key_to_Channel:
                for ch in Config_Table.Key_to_Channel[key]:
                    self.Update_QuadRF_channels.add(ch)
            else:
                self.logger.warn(f'{key} is not in Key_to_Channel. What channel does this key belong to? [change in Config_Table.py]')
                # TODO: raise an exception here

        # See whether OPX should also be updated
        if key in Config_Table.IOParametersMapping:
            io1 = Config_Table.IOParametersMapping[key]
            io2 = value
            if key in Values_Factor:
                valueType = Values_Factor[key][0]  # Can be one of: int, float, bool
                if len(Values_Factor[key]) == 1:
                    self.update_io_parameter(io1, io2)
                elif len(Values_Factor[key]) == 2:
                    valueFactor = Values_Factor[key][1]
                    # Perform factoring and cast to relevant type (int, float)
                    io2 = valueType(value * valueFactor)
                    self.update_io_parameter(io1, io2)
                elif len(Values_Factor[key]) == 3:  # If a function is attached to value change in Config_Table
                    valueFunction = Values_Factor[key][2]  # Get function we want to use to convert
                    valueFunction(self, value)  # Call the predefined function.
                    # Note: the function should call self.update_io_parameter(io1, io2),otherwise nothing happens
                    # io2 = valueType((valueFunction(self, value)))
            else:
                self.logger.error('%s not in Values_Factor!' % key)
            # TODO: why is the below commented out?
            # self.update_io_parameter(io1, io2)

        if update_parameters:
            self.update_parameters()

    def update_io_parameter(self, io1, io2):
        self.io1_list.append(io1)
        self.io2_list.append(io2)

    # TODO: need to better understand some of the vodoo going on in this func...
    # TODO: also - why aren't we disconnecting from QuadRF ?
    # TODO: why are we sending first pair of parameters and then the rest?
    def update_parameters(self):
        if len(self.Update_QuadRF_channels) != 0:
            quadController = QuadRFMOTController(initialValues=self.Exp_Values,
                                                 updateChannels=self.Update_QuadRF_channels, armChannels=True,
                                                 debugging=False,
                                                 continuous=False)  # updates values on QuadRF (uploads table)
            # Note at this point QuadRF is not yet armed, meanning it's not waiting for trigger, until OPX is updated
            self.Update_QuadRF_channels.clear()  # After only updating changed channels, reset this variable.

        # Don't update if there is nothing to update
        if self.io1_list.__len__() == 0:
            return

        # Sets the last value of both of them to be 0, such that the last loop execution will cause the OPX to continue
        # self.update_io_parameter(6, 3)  # Update filler to do 2 more cycles to let the MOT stabilize
        self.io1_list.append(0)
        self.io2_list.append(0)

        # Sending of 1st parameter pair
        self.qm.set_io_values(self.io1_list[0], self.io2_list[0])
        # IOs are sent one after the other, meaning that even if IO1 > 0, doesn't mean that IO2 has been updated.
        # This is solved by adding an initial pause/resume:
        while not self.job.is_paused():
            pass
        self.job.resume()

        # Go over the rest of the parameters:
        for i in range(1, self.io1_list.__len__(), 1):
            while not self.job.is_paused():
                pass
            self.qm.set_io_values(self.io1_list[i], self.io2_list[i])
            self.job.resume()

        # Reset the lists
        self.io1_list = []
        self.io2_list = []

        # Save Config_Table to a file with timestamp:
        self.save_config_table()

        # quadController.plotTables()

    def save_config_table(self, default_path='.\\Config Table Archive\\'):
        time_stamp = time.strftime("%d-%m-%Y %H_%M_%S", time.localtime())
        helm_coils_state = 'Working on it'  # self.getHelmholtzCoilsState()
        experiment_values_to_save = self.Exp_Values
        if 'opticalPowerCalibrationFunction' in experiment_values_to_save:
            experiment_values_to_save[
                'opticalPowerCalibrationFunction'] = 'Experiment values contained power calibration function. It is, however, impossible to save a function as JSON. Thus, this line is added as a warning'
        save_data = {
            'Time Stamp': time_stamp,
            'Exp_Values': self.Exp_Values,
            'Helmholtz Coils': helm_coils_state
        }
        try:
            # TODO: Once we know we want to write this file by-default, we need to uncomment the below
            # if not os.path.exists(default_path):
            #     os.makedirs(default_path)

            filename = os.path.join(default_path, time_stamp + '.json')
            with open(filename, 'w') as file:
                json.dump(save_data, file, indent=4)
        except Exception as err:
            self.logger.warn(f'Unable to open file {filename} for writing. {err}')
        pass

    def _plot(self, sequence_or_sequences, clear=True):
        if len(sequence_or_sequences) == 0:
            return

        if clear:
            plt.clf()

        if self.dbg_plot==None or clear==True:
            plt.figure(1)
            self.dbg_plot = True
        # If first element is an array, we're dealing here with sequences (plural)
        if hasattr(sequence_or_sequences[0], "__len__") or hasattr(np.array(sequence_or_sequences)[0], "__len__"):
            for seq in sequence_or_sequences:
                plt.plot(seq)
            return
        plt.plot(sequence_or_sequences)
