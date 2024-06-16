import os.path
import sys
import time
import json
import pathlib
import traceback
import matplotlib
import matplotlib.pyplot as plt
import os
import importlib
from importlib.machinery import SourceFileLoader

from Utilities.Utils import Utils
from Utilities.BDLogger import BDLogger
from Utilities.BDResults import BDResults
from Utilities.BDStreams import BDStreams
from Utilities.BDBatch import BDBatch
from Utilities.BDSound import BDSound
from Utilities.BDSocket import BDSocket
from Utilities.BDKeyboard import BDKeyboard
from Utilities.BDDialog import BDDialog
from Utilities.BDPlots import BDPlots

from Experiments.Enums.TerminationReason import TerminationReason
from Experiments.Enums.IOParameters import IOParameters as IOP
from Experiments.Enums.KeyToChannel import KeyToChannel as K2C
from Experiments.Enums.ExperimentMode import ExperimentMode as ExperimentMode


from Experiments.BaseExperiment import Config_Table
from Experiments.BaseExperiment import Config_Experiment as Config  # Attempt to load the default config (may be overriden later)
from Experiments.BaseExperiment import OPX_Code  # Attempt to load the OPX Code (maybe overriden later)
from Experiments.BaseExperiment.Config_Table import Initial_Values
from Experiments.BaseExperiment.Config_Table import Default_Values
from Experiments.Enums.ValuesTransformer import ValuesTransformer

from Experiments.QuadRF.quadRFMOTController import QuadRFMOTController
#from Experiments.QuadRF.QuadRFMOTController2 import QuadRFMOTController
#from Experiments.QuadRF.QuadRFMOTController2 import QuadRFMode


from logging import StreamHandler, Formatter, INFO
from qm.QuantumMachinesManager import QuantumMachinesManager
from qm.qua import *
import pymsgbox


class BaseExperiment:
    """
    This class is the base experiment for all other experiments we have:

    It contains the following utility functionalities:
        - BDLogger - Logger
        - Keyboard handlers
        - BDSound - Sounds player
        - Camera Connectivity
        - Identifier for computer (use self.login)
        - Sequence/Automation running activity with different parameters (pre_run, run, post_run, run_sequence)
        - Plot/Debugging <TBD>
        - BDResults - Results/Paths Services <TBD>
        - BDBatch - Batching Services (for gathering experiment samples)
        - BDStreams - Loading/Saving of raw stream data in real-time

    On the experiment side, it does the following:
        - Initializes OPX
        - Handles streams
        - Runs the following: MOT, PGC, Fountain, FreeFall
        - Handles Lock Error, Kappa_Ex and other params coming from Cavity Lock application and outer processes

    OPX Functionality:
        - Initializing
        - Updating Parameters
        - Updating IO Parameters

    Experiment Mode:
        - Live - experiment runs with all devices connected (OPX, Quad, Camera, etc.)
        - Offline - does not connect to lab infrastructure - Quad, OPX, Camera, etc.
        - Playback - allows to replay a set of pre-recorded playback data (does not connect to OPX/Quad)

    """

    def __init__(self, playback_parameters=None, save_raw_data=False, connect_to_camera=False, experiment_mode=ExperimentMode.LIVE):

        # Load the settings file
        self.base_settings = Utils.load_json_from_file('../BaseExperiment/settings.json')
        self.experiment_settings = Utils.load_json_from_file('./settings.json')
        self.settings = Utils.merge_jsons(self.base_settings, self.experiment_settings)

        # Set experiment-mode
        self.experiment_mode = experiment_mode

        # Setup console logger. We do this first, so rest of code can use logging functions.
        self.logger = BDLogger()

        # Generate UUID
        self.UUID = Utils.generate_UUID()

        # Get login and MAC address (used later to know if we're running in lab or locally)
        self.login = os.getlogin()
        self.mac_address = Utils.get_mac_address()

        # Initialize sound player
        self.bdsound = BDSound()

        # Set paths map
        self._set_paths_map()

        # Debugging plot
        self.dbg_plot = None

        # Should we keep the detectors raw data
        self.save_raw_data = save_raw_data

        self.transformer = ValuesTransformer()

        # Set playback definitions and initialize it
        if playback_parameters is None:
            self.playback = {'active': False}
        else:
            self.playback = playback_parameters
            self.playback['data_loaded'] = False
            self.playback['row_count'] = 0
            self.playback['streams'] = {}

        # If the 'active' param is set to True, this means we're in Playback mode
        if self.playback['active']:
            self.experiment_mode = ExperimentMode.PLAYBACK

        # Insert a prompt to check we're not running on live lab devices (OPX, QuadRF)
        non_lab_user = self.login in self.settings['permissions']['allowed_offline']
        if self.experiment_mode == ExperimentMode.LIVE and non_lab_user:
            bdd = BDDialog()
            text_res, button = bdd.prompt(self.settings['dialogs']['live_dialog'])
            if (button['button_name'] == 'Continue' and text_res.lower() != 'confirm') or (button['button_name']=='Exit'):
                sys.exit('Aborting.')
            if button['button_name'] == 'Offline':
                self.experiment_mode = ExperimentMode.OFFLINE

        # If we're in playback mode, ensure we do not connect to OPX/Quad
        if self.experiment_mode != ExperimentMode.LIVE:
            self._opx_skip = True
            self._quadrf_skip = True
        else:
            self._opx_skip = False
            self._quadrf_skip = False

        # Initialize the BDResults helper - for saving experiment results
        self.bd_results = BDResults(json_map_path=self.paths_map['cwd'], version="0.1", logger=self.logger)
        self.bd_results.create_experiment_run_folder()

        # Check network driver availability
        network_drive_letter = 'U'
        default_network_path = r'\\isi.storwis.weizmann.ac.il\Labs\baraklab'
        network_drive_available = self.bd_results.is_network_drive_available(f'{network_drive_letter}:\\Lab_2023')
        if not network_drive_available:
            self.logger.error(f'Network drive {network_drive_letter} is not available/connected. PLEASE FIX. (in Windows Explorer, click "Map network drive" and insert this path: {default_network_path})')
            sys.exit(1)

        # Tell logger we want to save the log
        self.logger.turn_save_on(log_path=self.bd_results.experiment_run_folder)

        self.logger.info(f'Starting experiment. Experiment root folder is here: {self.bd_results.experiment_run_folder}')

        # Initialize the BDBatch helper - to serve us when batching experiment samples
        self.batcher = BDBatch(json_map_path=self.paths_map['cwd'])

        # Initialize the BDStreams
        self.bdstreams = BDStreams(save_raw_data=save_raw_data, max_files_to_load=self.playback['max_files_to_load'], logger=self.logger)

        # Load Initial Values and Default Values - merge them together (Default Values prevails!)
        # These will be the experiment values
        self.Exp_Values = Utils.merge_multiple_jsons([Initial_Values, Default_Values])

        # Dynamically import the config-experiment and config-table and merge the values
        try:
            # If we're in playback mode, we load the config files from the playback folder
            if self.playback['active'] and self.playback['load_config_from_playback']:
                the_path = os.path.join(self.playback['playback_files_path'], 'Source Files', f'{self.experiment_name}_Config_Experiment.py')
                Config = SourceFileLoader(f"{self.experiment_name}_Config_Experiment", the_path).load_module()
                the_path = os.path.join(self.playback['playback_files_path'], 'Source Files', f'{self.experiment_name}_Config_Table.py')
                ConfigTable = SourceFileLoader(f"{self.experiment_name}_Config_Table", the_path).load_module()
            else:
                Config = importlib.import_module(f'{self.experiment_name}_Config_Experiment')
                ConfigTable = importlib.import_module(f'{self.experiment_name}_Config_Table')
            self.Exp_Values = Utils.merge_multiple_jsons([self.Exp_Values, ConfigTable.Experiment_Values])
        except Exception as err:
            self.info(f'Unable to import Config file ({err}). Loading *BaseExperiment* config files instead.')

        # Get the opx-control method from the OPX-Code file in the BaseExperiment
        opx_control = OPX_Code.opx_control

        # Attempt to dynamically import the OPX-Code file and get opx-control function from the current experiment
        try:
            OPX_Code_Module = importlib.import_module("OPX_Code")
            opx_control = OPX_Code_Module.opx_control
        except Exception as err:
            self.info(f'Unable to import OPX_Code file ({err}). Loading __**BaseExperiment**__ opx code instead.')
            # try:
            #     OPX_Code = importlib.import_module("Experiments.BaseExperiment.OPX_Code", "OPX_Code")
            # except Exception as err:
            #     self.warn(f'Unable to import OPX Code file from BaseExperiment ({err}).')

        # Initialize OPX with the relevant config and opx code
        # To Check on OPX status from Browser: http://132.77.54.230/ui/inventory/devices
        self.opx_definitions = {
            'connection': {'host': '132.77.54.230', 'port': '80'},  # Connection details
            'config': Config.config,  # OPX Configuration
            'streams': self.settings['streams'],
            'control': opx_control  # OPX Control code
        }

        # Set the streams in the bdstreams - so it can write them to disk
        self.bdstreams.set_streams_definitions(self.opx_definitions['streams'])
        self.streams = self.opx_definitions['streams']

        # Set mainloop flags
        self.ignore_data = False
        self.runs_status = None  # Uses the TerminationReason enum

        # Start listening on sockets (except when in playback mode)
        self.comm_messages = {}
        if self.experiment_mode == ExperimentMode.LIVE:
            self.bdsocket = BDSocket(connections_map=self.settings['connections'], writeable=self.comm_messages, server_id=self.UUID)
            self.bdsocket.run_server()

        # Attempt to initialize Camera functionality
        if self.experiment_mode == ExperimentMode.LIVE and connect_to_camera:
            self.connect_camera()
        else:
            self.info('Not connecting to camera. Not required for this experiment.')

        # Set keyboard handler
        self.bdkeyboard = BDKeyboard(self.settings['keyboard'])

        # Set plots handler
        self.bdplots = BDPlots(self.settings['figures'], plotter=self, logger=self.logger)

        # (a) Initialize QuadRF (b) Set experiment related variables (c) Initialize OPX
        self.initialize_experiment()

    def __del__(self):
        # Close all OPX related instances
        #if hasattr(self,'job'):
        # self.job.halt()
        if hasattr(self, 'qm'):
            self.qm.close()
        if hasattr(self, 'qmm'):
            self.qmm.close()
        pass

    # Initialize the Base-Experiment: QuadRF/MOT/PGC/FreeFall
    # (a) Initialize QuadRF (b) Set experiment related variables (c) Initialize OPX
    def initialize_experiment(self):

        # Initialize Experiment Variables
        self.initialize_experiment_variables()

        # Initialize QuadRF
        self.initiliaze_QuadRF()

        # Initialize the OPX
        self.initialize_OPX()

    # Initialize the QuadRF
    def initiliaze_QuadRF(self):
        if self._quadrf_skip:
            return

        self.QuadRFControllers = []

        # Note: So as not to connect again and again to QuadRF each time we update table, we now save the MOGDevice (actual QuadRF device) connected,
        # we hold this connection until update is finished, then we close the connection.
        # we do still hold the QuadRFController objects, for access to the table (read only!) when the experiment is running.

        qrfContr = QuadRFMOTController(MOGdevice=None,
                                       initialValues=self.Exp_Values,
                                       updateChannels=[1, 2, 4],
                                       # updateChannels=(1, 4),  # For constant Depump
                                       topticaLockWhenUpdating=False,
                                       debugging=False,  # True,
                                       #mode=QuadRFMode.DYNAMIC,
                                       continuous=False)
        self.QuadRFControllers.append(qrfContr)  # updates values on QuadRF (uploads table)

        # TODO: why aren't the values below part of the experiment values or configuration values?
        qrfContr2 = QuadRFMOTController(MOGdevice=qrfContr.device,
                                        initialValues={'CH3_freq': self.Exp_Values['CH3_continuous_freq'], 'CH3_amp': self.Exp_Values['CH3_continuous_amp']},
                                        updateChannels=[3],
                                        debugging=False,  # True,
                                        #mode=QuadRFMode.CONTINUOUS,
                                        continuous=True)
        self.QuadRFControllers.append(qrfContr2)  # updates values on QuadRF (uploads table)

        # TODO: do we need this?
        # qrfContr3 = QuadRFFrequencyScannerController(MOGdevice=qrfContr.device,
        #                                              channel=2,
        #                                              debugging=False)  # updates values on QuadRF (uploads table)
        # self.QuadRFControllers.append(qrfContr3)  # updates values on QuadRF (uploads table)
        ##### For continues Depump!!! (please remove channel 2 from qrfContr) #####
        # qrfContr3 = QuadRFMOTController(MOGdevice=qrfContr.device,
        #                                 initialValues={'CH2_freq': '133.325MHz', 'CH2_amp': '31dbm'},
        #                                 updateChannels=[2],
        #                                 debugging=True,
        #                                 continuous=True)  # updates values on QuadRF (uploads table)
        # self.QuadRFControllers.append(qrfContr3)  # updates values on QuadRF (uploads table)

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

        self.bdstreams.set_opx_handlers(opx_job=self.job, playback=self.playback['active'])
        pass

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
        self.experiment_name = os.path.basename(os.path.normpath(os.getcwd()))
        pass

    # ----------------------------------------------------------------------------
    # Keyboard handler Functions
    # ----------------------------------------------------------------------------

    def keyboard_handler__pause_unpause(self, key):
        str = 'Pausing' if not self.pause_flag else 'Continuing'
        self.logger.blue(f'SPACE pressed. {str} measurement.')
        self.pause_flag = not self.pause_flag
        pass

    def keyboard_handler__stop_experiment(self, key):
        self.logger.blue('ESC pressed. Stopping measurement.')
        self.updateValue("Experiment_Switch", False)
        self.MOT_switch(True)
        self.update_parameters()
        self.runs_status = TerminationReason.USER

    def keyboard_handler__turn_on_off_atoms(self, key):
        # Turn on the "switch" indication, but putting an "_" as a prefix
        # Note: user may have pressed A fast multiple times, so no need to add multiple "_" prefixes until the switch is handled
        if self.switch_atom_no_atoms[0] != '_':
            self.switch_atom_no_atoms = '_' + self.switch_atom_no_atoms
        pass

    def keyboard_handler__zoom_in_out_plots(self, key):
        self.bdplots.set_subplots_view(int(key))
        pass

    def should_terminate(self):
        """
        Decides if to terminate mainloop. By default, we return False, e.g. endless loop.
        Experiment class should override this if it has a different condition
        """

        # If we're in playback mode, and we had max iterations defined, we check it
        if self.playback['active'] and self.playback['max_iterations'] > 0:
            self.playback['max_iterations'] -= 1
            if self.playback['max_iterations'] == 0:
                return True

        return False

    def _read_k_ex(self):
        if not hasattr(self, 'k_ex_file'):
            self.k_ex_file = os.path.join(self.bd_results.get_custom_root('k_ex_root'), 'k_ex.npy')

        max_time = 10  # We will wait 0.1 sec for the file
        k_ex = None
        ticks = 0
        while k_ex is None:
            try:
                data = np.load(self.k_ex_file, allow_pickle=True)
                k_ex = np.abs(data)
                break
            except FileNotFoundError as err:
                pass  # We need to wait, the file may appear... (being written by another computer on a shared file-system)
            except Exception as err:
                self.logger.error(f'Error in loading k_ex file. {err}')
                break
            time.sleep(0.01)  # in seconds
            ticks = ticks + 1
            if ticks > max_time:
                break
        if k_ex is None:
            self.logger.error(f'Unable to find/load k_ex file.')

        return k_ex

    # Returns the error of the locking mechanism of the resonator to Rb line
    def _read_locking_error(self):
        if not hasattr(self, 'lock_error_file'):
            # Set the location of the locking error file - this is a generic file that serves as communication between
            # the locking application that runs on a different computer. The file contains the current PID error value.
            self.lock_error_file = os.path.join(self.bd_results.get_custom_root('locking_error_root'), 'locking_err.npy')

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

    def prompt(self, msg, default='', timeout=30000, title='', buttons=None):
        """
        Raise a prompt. If buttons are defined, it will raise a confirmation box with just buttons and no text.
        Returns: the text in case of prompt -or- the selected buttons in case of confirmation dialog
        """
        if buttons:
            res = pymsgbox.confirm(text=msg, title=title, buttons=buttons)
        else:
            res = pymsgbox.prompt(text=msg, title=title, default=default, timeout=timeout)
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
    # Experiment Management - Iteration/Sequence/Run - related methods
    # -------------------------------------------------------

    def get_batcher_map(self):
        raise Exception("'get_batcher_map' method must be implemented in subclass!")

    def experiment_analysis_subfolder(self):
        """
        If the experiment does not define any, we are not creating any subfolders
        """
        return ''

    def sequence_started(self, sequence_definitions):
        pass

    def sequence_ended(self, sequence_definitions):
        if not self.playback['active']:
            # Copy experiment Python source files (*.py) into playback folder
            python_source_files_path = self.paths_map['cwd']
            playback_files_path = os.path.join(self.bd_results.get_sequence_folder(sequence_definitions), 'playback', 'Source Files')
            self.bd_results.copy_files(source=python_source_files_path, destination=playback_files_path, opt_in_filter='.py', create_folder=True)
            self.info(f'Copied all Python source and config files from {python_source_files_path} into playback folder ({playback_files_path})')
        pass

    def iterations_started(self, sequence_definitions):
        # By default, the experiment results do not go into the "for_analysis" folder
        # (unless the user prompts to do so at the end of all iterations)
        self.save_results_for_analysis = False
        pass

    def iterations_ended(self, sequence_definitions):

        # Get the root experiment path where all iteration/sequence data was saved locally
        experiment_path = self.bd_results.get_experiment_run_folder()

        # If these results should be saved for analysis, copy them to analysis folder in Network drive
        if self.save_results_for_analysis and not self.playback['active']:
            if ':\\' in self.save_results_for_analysis or ':/' in self.save_results_for_analysis:
                dst = self.save_results_for_analysis
            else:
                dst = self.bd_results.get_custom_root('network_root')
                dst = dst.replace('{analysis_type}', self.save_results_for_analysis)
                subfolder = self.experiment_analysis_subfolder()
                if subfolder != '':
                    dst = os.path.join(dst, subfolder)

            self.bd_results.copy_folder(source=experiment_path, destination=dst)

        pass

    def pre_run(self, sequence_definitions, run_parameters):
        pre_comment = '' if 'pre_comment' not in run_parameters else run_parameters['pre_comment']

        # Prefix the pre_comment with the sequence number
        run_parameters['pre_comment'] = f'Iteration #{sequence_definitions["current_iteration"]+1}, Sequence [{sequence_definitions["sequence_name"]}]: {pre_comment}'


    def run(self, sequence_definitions, run_parameters):
        raise Exception("'run' method must implemented in subclass!")

    def post_run(self, sequence_definitions, run_parameters):
        raise Exception("'post_run' method must be implemented in subclass!")

    def run_sequence(self, sequence_definitions, run_parameters):

        # If there was no sequence defined, we create a fictitious sequence with a single iteration an no parameters
        if sequence_definitions is None or self.playback['active']:
            sequence_definitions = {
                'total_iterations': 1,
                'delay_between_iterations': None,  # seconds
                'sequence': [
                    {
                        'name': 'Mainloop',
                        'parameters': {}
                    }
                ]
            }

        # Notify we started iterations
        self.iterations_started(sequence_definitions)

        run_status = TerminationReason.SUCCESS
        num_iterations = sequence_definitions['total_iterations']
        for i in range(num_iterations):
            self.info(f'Starting iteration {i+1} (out of {num_iterations})')

            # Keep the current cycle, so pre_run, run, post_run can relate to it and do things once/first/last/etc.
            sequence_definitions['current_iteration'] = i
            sequence_definitions['last_iteration'] = (i == num_iterations-1)

            # Notify we started the sequence
            self.sequence_started(sequence_definitions)

            j = 0
            for sequence in sequence_definitions['sequence']:
                # Modify params
                merged_params = Utils.merge_jsons(run_parameters, sequence['parameters'])

                sequence_definitions['sequence_name'] = sequence['name']
                sequence_definitions['sequence_step'] = j

                j += 1
                sequence_definitions['last_sequence_step'] = (j == len(sequence_definitions['sequence']))

                sequence_definitions['last_iteration_and_last_sequence'] = sequence_definitions['last_iteration'] and sequence_definitions['last_sequence_step']

                # TODO: Temp change - we may want to remove this - look the usage of it and also remove
                self.sequence = sequence_definitions

                self.pre_run(sequence_definitions, merged_params)
                run_status = self.run(sequence_definitions, merged_params)
                self.post_run(sequence_definitions, merged_params)
                if run_status != TerminationReason.SUCCESS:
                    break

            # Notify we finished the sequence
            self.sequence_ended(sequence_definitions)

            if run_status != TerminationReason.SUCCESS:
                break

            if sequence_definitions['delay_between_iterations']:
                self.info(f'Sleeping')
                time.sleep(float(sequence_definitions['delay_between_iterations']))

        # Notify we finished running all iterations
        self.iterations_ended(sequence_definitions)

        self.info(f'Completed {num_iterations} iteration(s)! Run status: {run_status}')
        return run_status

    def _get_results_from_files(self):

        # If this is the first time this is called, load the data from the files
        if self.playback['data_loaded'] == False:
            self.playback['data_loaded'] = True
            self.playback['row_count'] = 0
            self.load_data_for_playback()

        # Advance all streams to hold the next data row from all_rows
        row = self.playback['row_count']
        for stream in self.streams.values():

            if 'skip' in stream and stream['skip']:
                continue

            if 'all_rows' in stream:

                # If there's no more playback data on this stream, we request to terminate experiment
                if row == len(stream['all_rows']):
                    self.runs_status = TerminationReason.PLAYBACK_END
                    self.info(f'No more data on stream {stream["name"]}. Ending playback.')
                    return

                # Ensure we eventually have nparray data
                data = stream['all_rows'][row]
                if type(data) == list:
                    data = np.array(data)

                stream['results'] = data

        # Advance the row for the next time
        self.playback['row_count'] += 1

        pass

    def load_data_for_playback(self):

        # Get the path for the playback files - either from playback parameters or from json definitions:
        playback_files_path = self.bd_results.get_custom_root('playback_data')
        if 'playback_files_path' in self.playback and self.playback['playback_files_path'] is not None:
            playback_files_path = self.playback['playback_files_path']

        # Take time
        start_load_time = time.time()

        # Choose what format we should load:
        if self.playback['old_format']:
            self.load_data_for_playback_using_old_format(playback_files_path)
        else:
            self.bdstreams.load_entire_folder(playback_files_path)

        self.info(f'Finished load playback files. Took {time.time()-start_load_time} seconds')

    def get_results_from_streams(self, playback=False):
        """
        Get the data from OPX or if in playback mode - from files
        """
        if self.streams is None:
            return

        if self.playback['active']:
            self._get_results_from_files()
            return

        # Get results from OPX streams
        self.bdstreams.get_results_from_opx_streams()

        # Get results from COMM streams (e.g. lock_error, kappa_ex, interference etc.)
        self.bdstreams.get_results_from_comm_channels(self.comm_messages)

        # Save streams to file
        if self.save_raw_data:
            playback_files_path = os.path.join(self.bd_results.get_sequence_folder(self.sequence), 'playback')
            Utils.ensure_folder_exists(playback_files_path)
            self.bdstreams.save_streams(playback_files_path)

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

        updated_at = []

        if key in self.Exp_Values:
            updated_at.append('Exp Values')

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
                updated_at.append('QuadRF')
                for ch in K2C[key]:
                    self.Update_QuadRF_channels.add(ch)
            else:
                self.logger.info(f'{key} is not in KeyToChannel. No QuadRF channel updated.')

        # -------------------------------------------------------
        # Now we check if this is a param we need to update OPX
        # -------------------------------------------------------
        if IOP.has(key):
            io1 = IOP.value_from_string(key)  # Get the Enum value from the string
            if self.transformer.knows(key):
                updated_at.append('OPX')
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

        if len(updated_at) == 0:
            logger.warn(f'Cannot update {key}. Not supported in your IOP or K2C tables.')
            return

        logger.info(f'updateValue: {key} updated. New Value: {value}. Updated: {", ".join(updated_at)}')

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

    def format_iteration_and_sequence(self):
        # Get the iteration/sequence/name
        iteration = self.sequence['current_iteration']
        sequence_step = self.sequence['sequence_step']
        name = self.sequence['sequence'][sequence_step]['name']
        str = f'Iteration {iteration+1} | Sequence {sequence_step+1} | {name}'
        return str

    def experiment_mainloop_delay(self):
        """
        Handles the delay of the experiment mainloop
        By default, we have a very short delay of few milli-seconds
        If we're in playback mode, the playback parameters prevail
        """
        if self.playback['active']:
            if self.playback['delay'] != -1:
                time.sleep(self.playback['delay'])
        else:
            time.sleep(0.01)
        pass

    # ------------------------------------------------------------
    # Plot related Functions
    # ------------------------------------------------------------

    def prepare_figures(self):

        self.bdplots.create_figures(new_figure=True)

        # Set figure title - (a) iteration and sequence (b) current date and time
        str = self.format_iteration_and_sequence()
        current_date_time = time.strftime("%H:%M:%S (%d/%m/%Y)")
        title = f'{current_date_time} | {str} | {self.UUID}'
        self.bdplots.set_figure_title(title)

        # Set keyboard listener
        self.bdkeyboard.set_listener(self.bdplots.get_figure(), self)

        self.plot_shown = False

        pass

    def plot_figures(self):

        # TODO: Move this entire section to "plot_figures" and into BaseExperiment
        if self.playback['active']:
            if self.playback['plot'] == 'NONE':
                return
            if self.playback['plot'] == 'LAST' and self.counter < (self.N - 10):
                return

        self.bdplots.plot_figures()


    def _plot(self, sequence_or_sequences, clear=True):
        """
        Debug method for fast plotting of one or more sequences of data

        - sequence_or_sequences - the data
        - clear - clear the existing figure or overlay

        Usage:
               self._plot(np.array[1, 3, 4, 9, 2], True)
        """
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
