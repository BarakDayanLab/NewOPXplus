import os.path
import sys
import time
import json
import socket
import pathlib
import matplotlib.pyplot as plt
import os
import importlib
import numpy as np

from Utilities.Utils import Utils
from Utilities.BDLogger import BDLogger
from Utilities.BDResults import BDResults
from Utilities.BDStreams import BDStreams
from Utilities.BDBatch import BDBatch
from Utilities.BDSound import BDSound

from Experiments.Enums.TerminationReason import TerminationReason
from Experiments.Enums.IOParameters import IOParameters as IOP
from Experiments.Enums.KeyToChannel import KeyToChannel as K2C

from Experiments.BaseExperiment import Config_Table
from Experiments.BaseExperiment import Config_Experiment as Config  # Attempt to load the default config (may be overriden later)
from Experiments.BaseExperiment import OPX_Code  # Attempt to load the OPX Code (maybe overriden later)
from Experiments.BaseExperiment.Config_Table import Initial_Values
from Experiments.BaseExperiment.Config_Table import Default_Values
from Experiments.Enums.ValuesTransformer import ValuesTransformer

from Experiments.QuadRF.quadRFMOTController import QuadRFMOTController

from logging import StreamHandler, Formatter, INFO
from qm.QuantumMachinesManager import QuantumMachinesManager
from qm.qua import *
import pymsgbox
from pynput import keyboard


class BaseExperiment:
    """
    This class is the base experiment for all other experiments we have:

    It contains the following utility functionalities:
        - BDLogger - Logger
        - Keyboard handlers
        - BDSound - Sounds player
        - Camera usage
        - Identifier for computer (use self.login)
        - Sequence/Automation running activity with different parameters (pre_run, run, post_run, run_sequence)
        - Plot/Debugging <TBD>
        - BDResults - Results/Paths Services <TBD>
        - BDBatch - Batching Services (for gathering experiment samples)

    On the experiment side, it does the following:
        - Initializes OPX
        - Handles streams
        - Runs the following: MOT, PGC, Fountain, FreeFall
        - Handles Lock Error coming from Cavity Lock application

    OPX Functionality:
        - Initializing
        - Updating Parameters
        - Updating IO Parameters
    """

    def __init__(self):

        # Setup console logger. We do this first, so rest of code can use logging functions.
        self.logger = BDLogger()

        # Get login - to understand on what computer we are running
        self.login = os.getlogin()

        # Initialize sound player
        self.bdsound = BDSound()

        # Set paths map
        self._set_paths_map()

        # Debugging plot
        self.dbg_plot = None

        self.transformer = ValuesTransformer()

        self._opx_skip = False
        self._quadrf_skip = False

        # Playback definitions
        self.playback = {
            "active": False,
            "data_loaded": False,
            "save_results": False,
            "plot_figures": False,
            "delay": 0.5,  # sec
            "row_count": 0,
            "streams": {}
        }
        if self.playback["active"]:
            self._opx_skip = True
            self._quadrf_skip = True


        # Initialize the BDResults helper - for saving experiment results
        self.bd_results = BDResults(json_map_path=self.paths_map['cwd'], version="0.1")

        # Initialize the BDBatch helper - to serve us when batching experiment samples
        self.batcher = BDBatch(json_map_path=self.paths_map['cwd'])

        # Initialize the BDStreams
        self.bdstreams = BDStreams(save_path=self.bd_results.get_custom_root('temp_root'))

        # Set the location of the locking error file - this is a generic file that serves as communication between
        # the locking application that runs on a different computer. The file contains the current PID error value.
        self.lock_error_file = os.path.join(self.bd_results.get_custom_root('locking_error_root'), 'locking_err.npy')

        # Load Initial Values and Default Values - merge them together (Default Values prevails!)
        # These will be the experiment values
        self.Exp_Values = Utils.merge_multiple_jsons([Initial_Values, Default_Values])

        # Dynamically import the config-experiment and config-table and merge the values
        try:
            Config = importlib.import_module("Config_Experiment")
            ConfigTable = importlib.import_module("Config_Table")
            self.Exp_Values = Utils.merge_multiple_jsons([self.Exp_Values, ConfigTable.Experiment_Values])
        except Exception as err:
            self.warn(f'Unable to import Config file ({err})')

        # Get the opx-control method from the OPX-Code file in the BaseExperiment
        opx_control = OPX_Code.opx_control

        # Attempt to dynamically import the OPX-Code file and get opx-control function from the current experiment
        try:
            OPX_Code_Module = importlib.import_module("OPX_Code")
            opx_control = OPX_Code_Module.opx_control
        except Exception as err:
            self.warn(f'Unable to import OPX_Code file ({err}). Loading local opx code from BaseExperiment')
            # try:
            #     OPX_Code = importlib.import_module("Experiments.BaseExperiment.OPX_Code", "OPX_Code")
            # except Exception as err:
            #     self.warn(f'Unable to import OPX Code file from BaseExperiment ({err}).')

        # Initialize OPX with the relevant config and opx code
        # To Check on OPX status from Browser: http://132.77.54.230/ui/inventory/devices
        self.opx_definitions = {
            'connection': {'host': '132.77.54.230', 'port': '80'},  # Connection details
            'config': Config.config,  # OPX Configuration
            'streams': Config.streams,  # OPX streams we're using
            'control': opx_control  # OPX Control code
        }

        # Set the streams in the bdstreams - so it can write them to disk
        self.bdstreams.set_streams_definitions(self.opx_definitions['streams'])

        # Set mainloop flags
        self.ignore_data = False
        self.runs_status = None  # Uses the TerminationReason enum

        # Open server-socket
        # TODO: complete implementation of this - listen on socket for outer process/machines to communicate with us
        #self._open_socket()

        # Attempt to initialize Camera functionality
        self.connect_camera()

        # TODO: What is this? Do we need this in BaseExperiment
        #self.init_spectroscopy()

        # Setup keyboard listener
        self.listener = keyboard.Listener(on_press=self._on_press, on_release=self._on_release)
        self.listener.start()
        self.keyPress = None

        # TODO: Listen to events - messages / keyboard

        # (a) Initialize QuadRF (b) Set experiment related variables (c) Initialize OPX
        self.initialize_experiment()

        pass

    def __del__(self):
        print('**** BaseExperiment Destructor ****')

        # Close all OPX related instances
        #if hasattr(self,'job'):
        # self.job.halt()
        if hasattr(self, 'qm'):
            self.qm.close()
        if hasattr(self, 'qmm'):
            self.qmm.close()
        pass

    def init_spectroscopy(self):
        # Setup of MW spectroscopy
        self.mws_logger = logging.getLogger("MWSpectroscopy")
        self.mws_logger.setLevel(INFO)
        handler = StreamHandler(sys.stdout)
        handler.setFormatter(Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
        self.mws_logger.addHandler(handler)

    # Initialize the Base-Experiment: QuadRF/MOT/PGC/FreeFall
    # (a) Initialize QuadRF (b) Set experiment related variables (c) Initialize OPX
    def initialize_experiment(self):

        # Initialize Experiment Variables
        self.initialize_experiment_variables()

        # Initialize QuadRF
        self.initiliaze_QuadRF()

        # Initialize the OPX
        self.initialize_OPX()

        pass

    # Initialize the QuadRF
    def initiliaze_QuadRF(self):
        if self._quadrf_skip: return

        self.QuadRFControllers = []

        # Note: So as not to connect again and again to QuadRF each time we update table, we now save the MOGDevice (actual QuadRF device) connected,
        # we hold this connection until update is finished, then we close the connection.
        # we do still hold the QuadRFController objects, for access to the table (read only!) when the experiment is running.
        qrfContr = QuadRFMOTController(initialValues=self.Exp_Values,
                                       updateChannels=(1, 2, 4),
                                       topticaLockWhenUpdating=False,
                                       debugging=True,
                                       continuous=False)
        self.QuadRFControllers.append(qrfContr)  # updates values on QuadRF (uploads table)

        # DEBUG DEBUG DEBUG - Test if it works - REMOVE
        for qrdCtrl in self.QuadRFControllers:
            test = qrdCtrl.get_channel_data(1)

        # TODO: why aren't the values below part of the experiment values or configuration values?
        qrfContr2 = QuadRFMOTController(MOGdevice=qrfContr.device,
                                        # TODO: remove this - we don't want 'Operation Mode' in the code anymore :-)
                                        #initialValues={'Operation_Mode': 'Continuous', 'CH3_freq': '90MHz', 'CH3_amp': '31dbm'},
                                        initialValues={'CH3_freq': '90MHz', 'CH3_amp': '31dbm'},
                                        updateChannels=[3],
                                        debugging=True,
                                        continuous=True)  # updates values on QuadRF (uploads table)
        self.QuadRFControllers.append(qrfContr2)  # updates values on QuadRF (uploads table)

        # TODO: do we need this?
        # qrfContr3 = QuadRFFrequencyScannerController(MOGdevice=qrfContr.device,
        #                                              channel=2,
        #                                              debugging=False)  # updates values on QuadRF (uploads table)
        #self.QuadRFControllers.append(qrfContr3)  # updates values on QuadRF (uploads table)

        self.Update_QuadRF_channels = set({})  # Only update these channels on QuadRF when UpdateParameters method is called [note: this is a python set]

        # Disconnect from the QuadRF
        qrfContr.disconnectQuadRF()

        pass

    # Initialize the OPX
    def initialize_OPX(self):
        if self._opx_skip: return

        # NOTE: if OPX is failing to start, maybe it requires a restart.
        # - Go to the OPX admin console in the browser
        # - Select the "Gear" icon, go to the "Operations" tab, click "Restart"
        # - Click the "Hardware" icon - wait for the status to turn from "Warning" to "Healthy"
        try:
            self.job = None
            self.qmm = QuantumMachinesManager(host=self.opx_definitions['connection']['host'], port=self.opx_definitions['connection']['port'])
            self.qmm.clear_all_job_results()
            self.qmm.reset_data_processing()
            self.qmm.close_all_quantum_machines()
            self.qm = self.qmm.open_qm(self.opx_definitions['config'])
            self.job = self.opx_definitions['control'](self, self.qm)
        except Exception as err:
            self.logger.warn(f'Unable to connect to OPX. {err}')
            raise Exception(f'Unable to connect to OPX. {err}')
        except KeyboardInterrupt as err:
            self.close_opx()
            self.logger.warn(f'Got user interrupt. Closing OPX. Maybe stopped Pycharm or Python process? {err}')

        self.io1_list = []
        self.io2_list = []

    def close_opx(self):
        if hasattr(self, 'job'):
            self.job.halt()
        if hasattr(self, 'qm'):
            self.qm.close()
        if hasattr(self, 'qmm'):
            self.qmm.close_all_quantum_machines()
            self.qmm.close()

    def initialize_experiment_variables(self):

        # ----------------------------------------
        # Free fall variables
        # ----------------------------------------

        self.FreeFall_duration = int(self.Exp_Values['FreeFall_duration'] * 1e6 / 4)
        self.Coils_timing = int(self.Exp_Values['Coil_timing'] * 1e6 / 4)
        if (self.Exp_Values['FreeFall_duration'] - self.Exp_Values['Coil_timing']) > 60:
            raise ValueError("FreeFall_duration - Coils_timing can't be larger than 60 ms")

        # ----------------------------------------
        # PGC variables:
        # ----------------------------------------

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
                self.error('The values for PGC_initial_amp_0 and PGC_final_amp_0 are too close or PGC_prep_duration is too long, might cause an ERROR!!!')
            if self.pgc_pulse_duration_minus > int(60 * 1e6 / 4):  # longer then 60ms
                self.error('The values for PGC_initial_amp_minus and PGC_final_amp_minus are too close or PGC_prep_duration is too long, might cause an ERROR!!!')
            if self.pgc_pulse_duration_plus > int(60 * 1e6 / 4):  # longer then 60ms
                self.error('The values for PGC_initial_amp_plus and PGC_final_amp_plus are too close or PGC_prep_duration is too long, might cause an ERROR!!!')
        self.pgc_initial_amp_0 = self.Exp_Values['PGC_initial_amp_0']
        self.pgc_initial_amp_minus = self.Exp_Values['PGC_initial_amp_minus']
        self.pgc_initial_amp_plus = self.Exp_Values['PGC_initial_amp_plus']
        self.pgc_final_amp_0 = self.Exp_Values['PGC_final_amp_0']
        self.pgc_final_amp_minus = self.Exp_Values['PGC_final_amp_minus']
        self.pgc_final_amp_plus = self.Exp_Values['PGC_final_amp_plus']
        # self.pgc_aom_chirp_rate = int(self.Exp_Values['PGC_final_Delta_freq'] * 1e3 / (self.Exp_Values['PGC_prep_duration'] * 1e6))  # [mHz/nsec], If needed pgc preparation duration must be constant!!!

        # ----------------------------------------
        # Fountain variables:
        # ----------------------------------------

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
                self.error('The values for Fountain_initial_amp_0 and Fountain_final_amp_0 are too close or Fountain_prep_duration is too long, might cause an ERROR!!!')
            if self.fountain_pulse_duration_minus > int(60 * 1e6 / 4):  # longer then 60ms
                self.error('The values for Fountain_initial_amp_minus and Fountain_final_amp_minus are too close or Fountain_prep_duration is too long, might cause an ERROR!!!')
            if self.fountain_pulse_duration_plus > int(60 * 1e6 / 4):  # longer then 60ms
                self.error('The values for Fountain_initial_amp_plus and Fountain_final_amp_plus are too close or Fountain_prep_duration is too long, might cause an ERROR!!!')
        self.fountain_initial_amp_0 = self.Exp_Values['Fountain_initial_amp_0']
        self.fountain_initial_amp_minus = self.Exp_Values['Fountain_initial_amp_minus']
        self.fountain_initial_amp_plus = self.Exp_Values['Fountain_initial_amp_plus']
        self.fountain_final_amp_0 = self.Exp_Values['Fountain_final_amp_0']
        self.fountain_final_amp_minus = self.Exp_Values['Fountain_final_amp_minus']
        self.fountain_final_amp_plus = self.Exp_Values['Fountain_final_amp_plus']
        self.fountain_aom_chirp_rate = int(self.Exp_Values['Fountain_final_Delta_freq'] * 1e3 / (self.Exp_Values['Fountain_prep_duration'] * 1e6))  # mHz/nsec

        # ----------------------------------------
        # OD and Depump measurement parameters:
        # ----------------------------------------

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

        # ----------------------------------------
        # MW spectroscopy parameters:
        # ----------------------------------------

        self.MW_start_frequency = int(100e6)  # [Hz]
        self.Pulse_Length_MW = 400  # [usec]
        self.Pulse_Length_OD = 20  # [usec]

        # Main Experiment:
        self.switch_atom_no_atoms = 'atoms'
        self.TOP2_pulse_len = int(Config.Probe_pulse_len / 4)  # [nsec]
        self.Calibration_time = 10  # [msec]

        pass

    def _set_paths_map(self):
        self.paths_map = {
            "name": __name__,
            "file": __file__,
            "cwd": os.getcwd(),
            "root": str(pathlib.Path(__file__).parent.resolve())
        }
        pass

    def handle_user_events(self):
        """
        Handle cases where user pressed ESC to terminate or ALT_SPACE to pause/continue measurements
        """
        if self.keyPress == 'ESC':
            self.logger.blue('ESC pressed. Stopping measurement.')
            self.updateValue("Experiment_Switch", False)
            self.MOT_switch(True)
            self.update_parameters()
            self.runs_status = TerminationReason.USER
        elif self.keyPress == 'ALT_SPACE' and not self.pause_flag:
            self.logger.blue('SPACE pressed. Pausing measurement.')
            self.pause_flag = True
            self.keyPress = None
        elif self.keyPress == 'ALT_SPACE' and self.pause_flag:
            self.logger.blue('SPACE pressed. Continuing measurement.')
            self.pause_flag = False
            self.keyPress = None
        elif self.keyPress == 'A':
            self.logger.blue('A pressed. Switching atoms/no-atoms.')
            # Turn on the "switch" indication, but putting an "_" as a prefix
            # Note: user may have pressed A fast multiple times, so no need to add multiple "_" prefixes until the switch is handled
            if self.switch_atom_no_atoms[0] != '_':
                self.switch_atom_no_atoms = '_' + self.switch_atom_no_atoms
            self.keyPress = None
        pass

    # Returns the error of the locking mechanism of the resonator to Rb line
    def _read_locking_error(self):
        if self.lock_error_file is None:
            self.logger.warn('Lock error file not defined. Not reading error signal.')
            return None

        max_time = 10  # We will wait 0.1 sec for the file
        lock_err = None
        ticks = 0
        while lock_err is None:
            try:
                data = np.load(self.lock_error_file, allow_pickle=True)
                lock_err = np.abs(data)
                break
            except FileNotFoundError as err:
                pass  # We need to wait, the file may appear... (being written by another computer on a shared file-system)
            except Exception as err:
                self.logger.error(f'Error in loading lock-error file. {err}')
                break
            time.sleep(0.01)  # in seconds
            ticks = ticks + 1
            if ticks > max_time:
                break
        if lock_err is None:
            self.logger.error(f'Unable to find/load lock-error file.')

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
        self.ignore_data = False

        if str(key) == "'q'" or str(key) == "'Q'":
            self.keyPress = 'Q'

        if str(key) == "'A'":
            self.keyPress = 'A'

        if key == keyboard.Key.esc:
            self.keyPress = 'ESC'

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
        self.camera = None
        try:
            from UtilityResources import MvCameraController
            self.camera = MvCameraController.MvCameraController()
        except Exception as e:
            self.warn(f'Could not connect to camera ({e})')
        return self.camera is not None

    def disconnect_camera(self):
        self.camera = None

    def is_camera_connected(self):
        return self.camera is not None

    def is_camera_disconnected(self):
        return self.camera is None

    # -------------------------------------------------------
    # Experiment Management related methods
    # -------------------------------------------------------

    def get_batcher_map(self):
        raise Exception("'get_batcher_map' method must be implemented in subclass!")

    def pre_run(self, run_parameters):
        raise Exception("'pre_run' method must be implemented in subclass!")

    def run(self, run_parameters):
        raise Exception("'run' method must implemented in subclass!")

    def post_run(self, run_parameters):
        raise Exception("'post_run' method must be implemented in subclass!")

    def run_sequence(self, sequence_definitions, run_parameters):

        run_status = TerminationReason.SUCCESS
        num_cycles = sequence_definitions['total_cycles']
        for i in range(num_cycles):
            self.info(f'Starting cycle sequence # {i}')

            for run in sequence_definitions['sequence']:
                # Modify params
                merged_params = Utils.merge_jsons(run_parameters, run['parameters'])

                self.pre_run(merged_params)
                run_status = self.run(merged_params)
                self.post_run(merged_params)
                if run_status != TerminationReason.SUCCESS:
                    break

            if run_status != TerminationReason.SUCCESS:
                break

            if sequence_definitions['delay_between_cycles']:
                self.info(f'Sleeping')
                time.sleep(sequence_definitions['delay_between_cycles'])

        self.info(f'Completed {num_cycles} cycles! Run status: {run_status}')
        return run_status

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
        """
        Given the streams' config, get the handles from OPX
        Usually called from the Python code before we want to get the data values in the stream
        """
        if 'streams' not in self.opx_definitions or self.opx_definitions['streams'] == {}:
            return

        self.streams = self.opx_definitions['streams']

        for key, value in self.streams.items():
            if self.playback["active"]:
                value['handler'] = "playback: data not loaded yet..."
            else:
                value['handler'] = self.job.result_handles.get(key)

    def _get_results_from_opx_streams(self):
        """
        Given the streams config and handles, wait and fetch the values from OPX
        """
        if self.streams is None:
            return

        # TODO: Should we first "wait_for_values" on all and only then "fetch_all"?
        for stream in self.streams.values():
            if stream['handler'] is not None:
                stream['handler'].wait_for_values(1)

        for stream in self.streams.values():
            if stream['handler'] is not None:
                stream['results'] = stream['handler'].fetch_all()
            else:
                stream['results'] = None

        # Save streams to file
        self.bdstreams.save_streams()

        pass

    def _get_results_from_files(self):

        # If this is the first time this is called, load the data from the files
        if self.playback['data_loaded'] == False:
            self.playback['data_loaded'] = True
            self.playback['row_count'] = 0
            self.load_data_for_playback()

        try:  # TODO: remove this try-except

            # Advance all streams to hold the next data row from all_rows
            row = self.playback['row_count']
            for stream in self.streams.values():
                if 'all_rows' in stream:
                    stream['results'] = stream['all_rows'][row]

        except Exception as err:
            print(err)

        # Advance the row for the next time
        self.playback['row_count'] += 1

        pass

    def load_data_for_playback(self):
        raise Exception('Not implemented in base class. Implement in superclass!')

    def get_results_from_streams(self, playback=False):
        """
        Get the data from OPX or if in playback mode - from files
        """
        if self.streams is None:
            return

        if self.playback['active']:
            self._get_results_from_files()
        else:
            self._get_results_from_opx_streams()
        pass

    def handle_user_atoms_on_off_switch(self):
        """
        This function checks if the user requested to switch - this is indicated by "_"
        If there needs to be a switch:
        - The "_" is removed (since "command" is handled)
        - The switch is flipped
        - MOT param is updated to OPX accordingly
        """
        # Value can be: "atom" or "!atom" or "_atom" or "_!atom"
        if self.switch_atom_no_atoms.startswith('_'):
            # If first letter after the '_' is '!', we are at no-atoms, so now we will turn MOT ON
            if self.switch_atom_no_atoms[1] == '!':
                MOT_on = True
                self.switch_atom_no_atoms = 'atoms'
            else:
                MOT_on = False
                self.switch_atom_no_atoms = '!atoms'
            self.MOT_switch(with_atoms=MOT_on, update_parameters=True)

    #------------------------------------------------------------------------
    # Updates a specific value - either related to (a) Operation Mode (b) QuadRF or (c) OPX
    #
    # Usage:
    #   - self.updateValue("Operation_Mode", "QRAM_Exp")
    #   - self.updateValue("QRAM_Exp_switch", False)
    #
    #  (*) If updating OPX values, we go through Values Factor map - to provide weights
    #------------------------------------------------------------------------
    def updateValue(self, key, value, update_parameters=False):

        if key in self.Exp_Values:
            self.Exp_Values[key] = value
            # ------------------------------------------------------------
            # Special cases - these happen ON-TOP of the regular update
            # Example: if we invoked updateValue("MOT_duration", 2):
            # - it will first invoke updateValue("MOT_rep", ...)
            # - then, it will return after it, and continue the function to also update the "MOT_Duration", the usual way
            # ------------------------------------------------------------
            if key == 'MOT_duration':
                # TODO: The below is not working - need to find updateMOT_rep()
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
            if key in K2C:
                for ch in K2C[key]:
                    self.Update_QuadRF_channels.add(ch)
            else:
                self.logger.warn(f'{key} is not in KeyToChannel. What channel does this key belong to? [change in Config_Table.py]')
                # TODO: raise an exception here

        # -------------------------------------------------------
        # Now we check if this is a param we need to update OPX
        # -------------------------------------------------------
        if IOP.has(key):
            io1 = IOP.value_from_string(key)  # Get the Enum value from the string
            if self.transformer.knows(key):
                mode = self.transformer.mode(key)
                if mode == 'SET_VALUE':
                    self.update_io_parameter(io1, value)
                elif mode == 'FACTOR_AND_CAST':
                    io2 = self.transformer.factor_and_cast(key, value)
                    self.update_io_parameter(io1, io2)
                elif mode == 'TRANSFORM':
                    #io2 = self.transformer.transform(key, value)
                    msg = f'\n\nTransformer currently does not handle mode {mode}. Call Dror.\n\n'
                    self.logger.error(msg)
                    raise Exception(msg)
                else:
                    pass
            else:
                msg = f'"{key}" is not in Values_Factor map!'
                self.logger.error(msg)
                raise Exception(msg)

        if update_parameters:
            self.update_parameters()

    def update_io_parameter(self, io1, io2):
        if self.playback["active"]:
            return

        self.io1_list.append(io1)
        self.io2_list.append(io2)

    def update_parameters(self):
        """
        Update the OPX - "push" all [key,value] pairs that are awaiting on the io1_list and io2_list
        We do so by setting one pair, then waiting for the OPX to pause its work, signaling us we can send the next pair
        """
        if self.playback["active"]:
            return

        # If needed, update the QuadRF controller
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

        # TODO: why are we sending first pair of parameters and then the rest? Can we fold everything into a single loop?
        # TODO: the code below is what I suggest :-)
        # for i in range(0, self.io1_list.__len__()):
        #     self.qm.set_io_values(self.io1_list[i], self.io2_list[i])
        #     while not self.job.is_paused():
        #         pass
        #     self.job.resume()

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

        # quadController.plotTables()

    def update_values(self, array_of_parameters, update_parameters=False):
        """
        Iterate over array of [key,value] elements and update each
        """
        for index, key_value in enumerate(array_of_parameters):
            update_parameters_now = False if index<len(array_of_parameters) else update_parameters
            self.updateValue(key_value[0], key_value[1], update_parameters_now)
        pass


    # TODO: this function is currently not in use. It was used to save exp values in 2 cases:
    # TODO: (1) At the end of an experiment (2) Whenever we updateValue.
    # TODO: Now, (1) is done automatically. (2) is less important so it was leftout
    def save_config_table(self, default_path='.\\Config Table Archive\\'):
        time_stamp = time.strftime("%d-%m-%Y %H_%M_%S", time.localtime())
        helm_coils_state = 'Working on it'  # self.getHelmholtzCoilsState()
        experiment_values_to_save = self.Exp_Values
        if 'opticalPowerCalibrationFunction' in experiment_values_to_save:
            experiment_values_to_save['opticalPowerCalibrationFunction'] = 'Experiment values contained power calibration function. It is, however, impossible to save a function as JSON. Thus, this line is added as a warning'
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

# ------------------ Utility/Fast-Access Functions

    def MOT_switch(self, with_atoms, update_parameters=False):
        if with_atoms:
            self.update_io_parameter(IOP.MOT_SWITCH_ON.value, 0)
        else:
            self.update_io_parameter(IOP.MOT_SWITCH_OFF.value, 0)
        if update_parameters:
            self.update_parameters()

    def Linear_PGC_switch(self, Bool):
        self.update_io_parameter(IOP.LINEAR_PGC_SWITCH.value, Bool)

    def Experiment_Switch(self, experiment_on_off):
        self.update_io_parameter(IOP.EXPERIMENT_SWITCH.value, experiment_on_off)

    def AntiHelmholtz_Delay_switch(self, Bool):
        self.update_io_parameter(IOP.ANTIHELMHOLTZ_DELAY_SWITCH.value, Bool)

    def Change_OD_freq(self, freq):
        self.update_io_parameter(IOP.OD_FREQUENCY, int(Config.IF_AOM_Spectrum + freq))

    def Max_Probe_counts_switch(self, Bool):
        self.update_io_parameter(IOP.MAX_PROBE_COUNTS_SWITCH.value, Bool)

    # ## Antihelmholtz delay variable update functions: ##

    def update_AntiHelmholtz_delay(self, new_delay):  # In units of [msec]
        self.update_io_parameter(IOP.ANTIHELMHOLTZ_DELAY.value, int(new_delay * 1e6 / 4))

    ## update N_snaps: ##

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # TODO: Are the 4 functions below called by user? Not only by code? Because they use self.Exp_values
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ## PGC variable update functions: ##
    def update_N_Snaps(self, N_Snaps):  # In units of [msec]
        if self.Exp_Values['Buffer_Cycles'] < 0:
            self.update_io_parameter(IOP.BUFFER_CYCLES.value, 3)  # Update 3 buffer cycles
        self.update_io_parameter(IOP.N_SNAPS.value, N_Snaps)  # Update N_snaps

    def update_PGC_prep_time(self, prep_time):  # In units of [msec]
        self.pgc_prep_duration = int(prep_time * 1e6 / 4)
        self.update_io_parameter(IOP.PGC_PREP_DURATION.value, int(prep_time * 1e6 / 4))

    ## Fountain variable update functions: ##

    def update_fountain_prep_time(self, prep_time):  # In units of [msec]
        self.fountain_prep_duration = int(prep_time * 1e6 / 4)
        self.update_io_parameter(IOP.FOUNTAIN_PREP_TIME.value, int(prep_time * 1e6 / 4))

    def update_fountain_Delta_freq(self, df):
        self.fountain_aom_chirp_rate = int(df * 1e3 / (self.fountain_prep_duration * 4))  # mHz/nsec
        self.logger.info(self.fountain_aom_chirp_rate)
        self.update_io_parameter(IOP.FOUNTAIN_FINAL_DELTA_FREQ.value, self.fountain_aom_chirp_rate)

    # TODO: the below was not yet tested
    # def update_mot_duration(self, duration):
    #     # Update both the MOT duration and the MOT repetitions
    #     rep = int(np.ceil((duration * 1e6) / Config.MOT_pulse_len))
    #     self.logger.info(f'Updating MOT duration to {duration} and MOT repetitions to {rep}')
    #     self.update_io_parameter(IOP.MOT_DURATION.value, duration)
    #     self.update_io_parameter(IOP.MOT_REPETITION .value, rep)


        pass

    ## Push beam params
    def update_PushBeam_duration(self, duration):
        self.update_io_parameter(IOP.PUSHBEAM_DURATION, int(duration * 1e3 / 4))  # In [us]

    def update_PushBeam_amp(self, Amp):
        self.update_io_parameter(IOP.PUSHBEAM_AMPLITUDE, float(Amp))  #

    def update_PushBeam_frequency(self, freq):
        self.update_io_parameter(IOP.PUSHBEAM_FREQUENCY, int(Config.IF_AOM_OD + freq))  # In [us]

    ## Measuring variable update functions: ##

    def MeasureOD(self, OD_Time):
        self.update_io_parameter(41, int(OD_Time * 1e6 / 4))  # TODO: Q: 41 in IOParameters is TRIGGER_DELAY

    def MeasureOD_SNSPDs_delay(self, OD_Delay):
        self.update_io_parameter(44, int(OD_Delay * 1e6 / 4))  # TODO: Q: 44 in IOParameters is PULSE_1_DECAY_DURATION

    def MeasureOD_SNSPDs_duration(self, duration):
        self.update_io_parameter(45, int(duration / self.M_window * 1e6 / 4))  # TODO: Q: 45 in IOParameters is N_SNAPS

    ## MW spectroscopy variable update functions: ##

    def MW_spec_frequency(self, freq):
        self.update_io_parameter(IOP.MW_SPEC_FREQ.value, int(freq))  # In unit of [Hz]

    def MW_spec_MW_pulse_duration(self, p_length):
        self.update_io_parameter(IOP.MW_SPEC_MW_PULSE.value, p_length * int(1e3 / 4))  # In units of [us]

    def MW_spec_OD_pulse_duration(self, p_length):
        self.update_io_parameter(IOP.MW_SPEC_OD_PULSE.value, p_length * int(1e3 / 4))  # In units of [us]

    def MW_spec_Repetition_times(self, reps):
        self.update_io_parameter(IOP.MW_SPEC_REPETITION_TIMES.value, reps)

    def MW_spec_Delta_freq(self, D_f):
        self.update_io_parameter(IOP.MW_SPEC_DELTA_FREQ.value, int(D_f))  # In units of [Hz]

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
