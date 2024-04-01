from topticaLockController import TopticaLockController
from Experiments.BaseExperiment.Config_Table import Initial_Values as config
from Experiments.Enums.Phases import Phases
from Experiments.Enums.Phases import Phases_Names
from Experiments.QuadRF.QuadRFController2 import QuadRFController
from Utilities.BDLogger import BDLogger
from Utilities.Utils import Utils
from enum import Enum
from enum import auto


class Channel(Enum):
    TOPTICA_1 = 1
    DEPUMP = 2
    TOPTICA_2 = 3
    REPUMP = 4

CHANNELS = [
    {'name': 'Toptica 1', 'limit': '22dbm'},  # 158.49 mW
    {'name': 'Depump', 'limit': '20dbm'},  # 100mW
    {'name': 'Toptica 2', 'limit': '32dbm'},  # ~1.58W
    {'name': 'Repump', 'limit': '4dbm'}  # ~2.5mW
]


class QuadRFChannel:
    """
    This class represents what happens in a given Channel at a given time
    What is the Frequency, Amplitude and Phase the Quad plays.

    It can be initialized with string values that include the unit (e.g. "100kHz", "5mW") or without units
    - and then default units are used: "MHz", "dBm", "degrees")

    If '0x0' is passed - it has special meaning - in the case of Power - do not play this signal
    """
    def __init__(self, name, freq, amp, phase=None, duration=None):
        self.name = name  # The channel name
        self.freq, self.freq_unit, self.freq_enabled = self._extract_value_and_unit(str(freq), float, 'MHz')
        self.amp, self.amp_unit, self.amp_enabled = self._extract_value_and_unit(str(amp), float, 'dBm')
        self.phase, self.phase_unit, self.phase_enabled = self._extract_value_and_unit(str(phase), float, 'degrees')
        self.duration, self.duration_unit, self.duration_enabled = self._extract_value_and_unit(str(duration), float, 'ms')
        pass

    def __getitem__(self, item):
        return getattr(self, item)

    def _extract_value_and_unit(self, value_str, value_cast, default_unit):
        """
        Extract the value and the unit from a string, then cast it to type_cast
        (e.g. "10kW" -> value: 10, unit: kW)

        If value_str is '0x0' - we return 0 as value, '' as unit BUT False as "enabled" flag

        TODO: we can make this more robust - remove spaces, use Regex to recognize specific units.
        """
        if value_str is None or value_str == '' or value_str == '0x0':
            return 0, '', False

        value = ''.join([n for n in str(value_str) if (n in['-', '.'] or n.isdigit())]) if '0x0' not in value_str else '0x0'
        if len(value_str) > len(value):
            unit = value_str[len(value):]
        else:
            unit = default_unit
        return value_cast(value), unit, True


class QuadRFPhase:
    """
    This class represents what the Quad does in a specific Phase - which has a given duration

    It holds 4 Channels and the Duration of the phase
    Each Channel holds a Freq/Amp/Phase that tells the Quad what signal to play
    """
    def __init__(self, duration, channel_values):
        # Ensure duration is typed as a string, with the proper formatting
        # if type(duration) is not str:
        #     duration = '{:.2f}ms'.format(duration)  # 2 decimal points, 'ms' added

        # Initialize the 4 channels
        self.channels = []
        for i in range(0, 4):
            freq = '0x0' if channel_values is None else channel_values[i][0]
            amp = '0x0' if channel_values is None else channel_values[i][1]
            phase = '0x0' if channel_values is None or len(channel_values[i]) == 2 else channel_values[i][2]
            name = Channel(i+1).name
            channel = QuadRFChannel(name, freq, amp, phase, duration)
            self.channels.append(channel)

        # All 4 channels will play for the exact duration, so the phase duration is the same as the first channel duration
        self.duration = self.channels[0].duration

        pass

    def channel(self, channel_number):
        """
        Get a specific QuadRFChannel by its index
        """
        return self.channels[channel_number]

    def get_channel_cmd(self, channel_number):
        """
        Return the command.
        It should be in the format:

        "TABLE,ENTRY,channel,[freq,power,phase,duration,flags]"
        """

        duration_str = '0x0' if self.duration==0 else '{:.2f}ms'.format(self.duration)  # 2 decimal points, 'ms' added
        freq_str = str(self.channels[channel_number].freq) + 'Hz'
        amp_str = str(self.channels[channel_number].amp) + 'dbm'
        cmd = f'TABLE,APPEND,{channel_number}, {freq_str}, {amp_str}, 0x0, {duration_str}'
        return cmd

    def __getitem__(self, item): # This enables calling attributes as such: self['attr']  = self.attr
        return getattr(self, item)


class QuadRFMode(Enum):
    OFF = auto()
    CONTINUOUS = auto()
    DYNAMIC = auto()


class QuadRFMOTController(QuadRFController):
    def __init__(self, MOGdevice=None, devPort='169.254.231.196', initialValues=None, updateChannels=(1,2,3,4), armChannels=True, opticalPowerCalibration=None, topticaLockWhenUpdating=False, debugging=False, mode=QuadRFMode.DYNAMIC):

        # Optical Power Calibration
        if opticalPowerCalibration is None and 'opticalPowerCalibrationFunction' in initialValues:
            opticalPowerCalibration = initialValues['opticalPowerCalibrationFunction']

        super(QuadRFMOTController, self).__init__(MOGdevice, devPort, opticalPowerCalibration=opticalPowerCalibration, debugging=debugging)
        self.initialValues = config if initialValues is None else initialValues

        self.logger.info(f'Initializing QuadRF @ {devPort}')

        # ---------------------------------------------------------------
        # Define the default amplitudes we will use for each channel
        # ---------------------------------------------------------------

        # self.AmplifiedAmp = 31  # dbm
        # self.Amp_Ch1 = 15.25  # dbm (15.25 -> red AOM connected to 2W amp, +++++with 15db attenuator on input++++.)
        self.Amp_Ch1 = 16.95  # dbm (15.25 -> red AOM connected to 2W amp, +++++with 15db attenuator on input++++.)
        # self.Amp_Ch2 = 4  # dbm, unAmplified (that is, amplified by an external amp. Maybe I should have called it the other way around)
        self.Amp_Ch2 = 31  # dbm, going straight into the AOM
        self.Amp_Ch3 = 1.5  # dbm
        self.Amp_Ch4 = 4.75  # dBm
        self.zeroAmp = '0x0'

        # self.UnAmplifiedAmp = 5.75  # dbm, unAmplified (that is, amplified by an external amp. Maybe I should have called it the other way around)
        self.topticaController = None
        self.topticaLockWhenUpdating = topticaLockWhenUpdating
        self.mode = mode

        # Load the Quad tables based on the mode we are running
        if self.mode == QuadRFMode.CONTINUOUS:  # run continuous mode then quit
            self.continuousTablesForChannels(updateChannels)
        elif self.mode == QuadRFMode.OFF:  # This is yet untested!
            self.turnChannelsOff(holdLocking=True)
        elif self.mode == QuadRFMode.DYNAMIC:
            self.upload_mot_tables(values=self.initialValues, updateChannels=updateChannels, armChannels=armChannels)
        return

    def __getitem__(self, item):  # This enables calling attributes as such: self['attr']  = self.attr
        return getattr(self, item)

    def initQuadRFValues(self, values=None):
        """
        TODO: Document this + make the phases a dictionary (as opposed to putting it on self....
        """
        if values is None:
            values = self.initialValues

        self.phases_order = values['Phases_Order']
        self.trigger_on_phase = values['Triggering_Phase']

        AOMOffFreq = values['MOT_AOM_freq'] + values['AOM_Off_Detuning']

        # -----------------  MOT ------------------------
        self.MOT = QuadRFPhase(duration=values['MOT_duration'],
                               channel_values=((values['MOT_freq'], self.Amp_Ch1),
                                               (AOMOffFreq, self.zeroAmp),
                                               (values['MOT_AOM_freq'], self.Amp_Ch3),
                                               (values['AOM_Repump_freq'], self.Amp_Ch4)))

        # MOT Delay (turn everything off, except for repump)
        self.Post_MOT_delay = QuadRFPhase(duration=values['Post_MOT_delay'],
                                          channel_values=((values['MOT_freq'], self.Amp_Ch1),
                                                          (AOMOffFreq, self.zeroAmp),
                                                          (AOMOffFreq, self.zeroAmp),
                                                          (values['AOM_Repump_freq'], self.Amp_Ch4)))
        # -----------------  PGC ------------------------
        PGC_final_amp_delta = Utils.amplitudeMultiplierToDBm(values['PGC_final_amp'])

        PGC_final_delta_amp_0 = Utils.amplitudeMultiplierToDBm(values['PGC_final_amp_0'])
        PGC_final_delta_amp_plus = Utils.amplitudeMultiplierToDBm(values['PGC_final_amp_plus'])
        PGC_initial_delta_amp_0 = Utils.amplitudeMultiplierToDBm(values['PGC_initial_amp_0'])
        PGC_initial_delta_amp_plus = Utils.amplitudeMultiplierToDBm(values['PGC_initial_amp_plus'])

        PGC_prep_duration = values['PGC_prep_duration'] if values['PGC_duration'] > 0 else 0
        PGC_duration = values['PGC_duration'] - values['PGC_prep_duration'] if values['PGC_duration'] > 0 else 0

        # PGC-Prep
        # Prep phase, by default takes initial / final values from phase before PGC
        self.PGC_prep = QuadRFPhase(duration=PGC_prep_duration,
                                    channel_values=((values['PGC_initial_freq'], self.Amp_Ch1),
                                                    (AOMOffFreq, self.zeroAmp),
                                                    (values['MOT_AOM_freq'] - values['PGC_initial_Delta_freq'], self.Amp_Ch3 + PGC_initial_delta_amp_plus),
                                                    (values['Repump_PGC_freq'], self.Amp_Ch4)))

        # PGC (constant)
        self.PGC = QuadRFPhase(duration=PGC_duration,
                               channel_values=((values['PGC_final_freq'], self.Amp_Ch1 + PGC_final_amp_delta),
                                               (AOMOffFreq, self.zeroAmp),
                                               (values['MOT_AOM_freq'] - values['PGC_final_Delta_freq'], self.Amp_Ch3 + PGC_final_delta_amp_plus),
                                               (values['Repump_PGC_freq'], self.Amp_Ch4)))

        # -----------------  Fountain ------------------------
        Fountain_duration = values['Fountain_duration'] - values['Fountain_prep_duration'] if values['Fountain_duration'] > 0 else 0
        Fountain_prep_duration = values['Fountain_prep_duration'] if values['Fountain_duration'] > 0 else 0

        Fountain_final_delta_amp_0 = Utils.amplitudeMultiplierToDBm(values['Fountain_final_amp_0'])
        Fountain_final_delta_amp_plus = Utils.amplitudeMultiplierToDBm(values['Fountain_final_amp_plus'])
        Fountain_initial_delta_amp_0 = Utils.amplitudeMultiplierToDBm(values['Fountain_initial_amp_0'])
        Fountain_initial_delta_amp_plus = Utils.amplitudeMultiplierToDBm(values['Fountain_initial_amp_plus'])

        # Fountain-Prep (dynamic)
        self.Fountain_prep = QuadRFPhase(duration=Fountain_prep_duration,
                                         channel_values=((values['Fountain_initial_freq'], self.Amp_Ch1),
                                                    # (values['MOT_AOM_freq'] + values['Fountain_initial_Delta_freq'], self.Amp_Ch2 + Fountain_initial_delta_amp_0),
                                                    (AOMOffFreq, self.zeroAmp),
                                                    (values['MOT_AOM_freq'] - values['Fountain_initial_Delta_freq'], self.Amp_Ch3 + Fountain_initial_delta_amp_plus),
                                                    (values['Repump_PGC_freq'], self.Amp_Ch4))) # prep phase, takes initial / final values from PGC & Fountain phases
        # Fountain (constant)
        self.Fountain = QuadRFPhase(duration=Fountain_duration,
                                    channel_values=((values['Fountain_final_freq'], self.Amp_Ch1),
                                                    # (values['MOT_AOM_freq'] + values['Fountain_final_Delta_freq'], self.Amp_Ch2 + Fountain_final_delta_amp_0),
                                                    (AOMOffFreq, self.zeroAmp),
                                                    (values['MOT_AOM_freq'] - values['Fountain_final_Delta_freq'], self.Amp_Ch3 + Fountain_final_delta_amp_plus),
                                                    (values['Repump_PGC_freq'], self.Amp_Ch4)))
        # ----------------- Free-fall -----------------
        PrePulse_delta_amp_Repump = float(Utils.amplitudeMultiplierToDBm(values['PrePulse_Repump_amp']))
        self.Free_Fall = QuadRFPhase(duration=values['PrePulse_duration'],
                                     channel_values=((values['PrePulse_CH1_freq'], self.Amp_Ch1),
                                     # channel_values=((values['PGC_final_freq'], self.Amp_Ch1 + PGC_final_amp_delta),
                                                     (values['PrePulse_CH2_freq'], self.Amp_Ch2),
                                                     (AOMOffFreq, self.Amp_Ch3),
                                                     (values['Pulse_1_CH4_Freq'], self.Amp_Ch4 + PrePulse_delta_amp_Repump)))

        # ----------------- Pulse 1 -----------------
        Pulse_1_delta_amp_i = float(Utils.amplitudeMultiplierToDBm(values['Pulse_1_amp_i']))
        Pulse_1_delta_amp_f = float(Utils.amplitudeMultiplierToDBm(values['Pulse_1_amp_f']))
        Pulse_1_delta_amp_Repump = float(Utils.amplitudeMultiplierToDBm(values['Pulse_1_Repump_amp']))
        Pulse_1_delta_amp_Depump = float(Utils.amplitudeMultiplierToDBm(values['Pulse_1_Depump_amp']))
        self.Pulse_1 = QuadRFPhase(duration=values['Pulse_1_duration'] + values['M_off_time'] - values['Pulse_1_decay_duration'],
                                   channel_values=((values['Pulse_1_CH1_Freq_f'], self.Amp_Ch1 + Pulse_1_delta_amp_f),
                                                   (values['PrePulse_CH2_freq'], self.Amp_Ch2 + Pulse_1_delta_amp_Depump),  # (AOMOffFreq, self.zeroAmp),
                                                   (values['Pulse_1_CH_2_3_Freq'], self.Amp_Ch3 + Pulse_1_delta_amp_f),
                                                   (values['Pulse_1_CH4_Freq'], self.Amp_Ch4 + Pulse_1_delta_amp_Repump)))

        self.Pulse_1_decay = QuadRFPhase(duration=values['Pulse_1_decay_duration'],
                                         channel_values=((values['Pulse_1_CH1_Freq_i'], self.Amp_Ch1 + Pulse_1_delta_amp_i),
                                                         (values['PrePulse_CH2_freq'], self.Amp_Ch2 + Pulse_1_delta_amp_Depump),
                                                         (values['Pulse_1_CH_2_3_Freq'], self.Amp_Ch3 + Pulse_1_delta_amp_i),
                                                         (values['Pulse_1_CH4_Freq'], self.Amp_Ch4 + Pulse_1_delta_amp_Repump)))

        # ----------------- Inter-Pulses -----------------
        self.Inter_Pulses = QuadRFPhase(duration=values['InterPulses_duration'],
                                        channel_values=((values['MOT_freq'], self.Amp_Ch1),
                                                        (AOMOffFreq, self.zeroAmp),
                                                        (AOMOffFreq, self.Amp_Ch3),
                                                        (AOMOffFreq, self.Amp_Ch4)))
        # ----------------- Pulse 2 -----------------
        self.Pulse_2 = QuadRFPhase(duration=values['Pulse_2_duration'],
                                   channel_values=((values['Pulse_2_CH1_Freq'], self.Amp_Ch1),
                                                   (AOMOffFreq, self.zeroAmp),
                                                   (values['Pulse_2_CH_2_3_Freq'], self.Amp_Ch3),
                                                   (values['Pulse_2_CH4_Freq'], self.Amp_Ch4)))

        # ----------------- Post-Pulses -----------------
        PostPulsesDuration = values['FreeFall_duration'] - values['PrePulse_duration'] - values['Pulse_1_duration'] - values['M_off_time'] - values['Pulse_2_duration'] - values['InterPulses_duration']
        if PostPulsesDuration > 0:
            # self.PostPulses = QuadRFPhase(duration=PostPulsesDuration, channel_values=((values['MOT_freq'], -20), (AOMOffFreq, self.Amp_Ch2),(AOMOffFreq, self.Amp_Ch3), (AOMOffFreq, self.Amp_Ch4)))
            self.PostPulses = QuadRFPhase(duration=PostPulsesDuration, channel_values=((values['MOT_freq'], self.Amp_Ch1),
                                                                                       (AOMOffFreq, self.zeroAmp),
                                                                                       (AOMOffFreq, self.Amp_Ch3),
                                                                                       (values['AOM_Repump_freq'], self.Amp_Ch4)))

        else:
            self.PostPulses = QuadRFPhase(channel_values=((values['MOT_freq'], self.Amp_Ch1),
                                                          (AOMOffFreq, self.zeroAmp),
                                                          (values['Flash_freq'], self.Amp_Ch3),
                                                          (values['AOM_Repump_freq'], self.Amp_Ch4)))

        self.changeLastPhaseDuration()

    def changeLastPhaseDuration(self):
        """
        Making sure the last phase duration is 0.01ms - this way this phase goes on until another trigger arrives.
        The previous phase is dependent on the experiments phase order
        """
        d = '0.01ms'

        previous_phase_index = (self.phases_order.index(self.trigger_on_phase)-1) % len(self.phases_order)
        last_phase_index = self.phases_order[previous_phase_index]

        triggeringPhase = Phases_Names[self.trigger_on_phase]
        lastPhase = Phases_Names[last_phase_index]

        if self.debugging:
            self.logger.debug(f'Triggering on {triggeringPhase}, {lastPhase} will continue until trigger.')

        # e.g., if we trigger on Pulse_1, lase phase is FreeFall
        self[lastPhase].duration = d

    def upload_mot_tables(self, values=None, updateChannels=(1,2,3,4), armChannels=True):
        if values is None:
            values = self.initialValues

        # --- Handle Toptica lock if necessary -----
        if 1 in updateChannels and self.topticaLockWhenUpdating:
            self.setTopticaLock(False)

        # Deletes all existing tables on QuadRF and sets limits [currently hard coded!] on each channel
        self.prepare_channels_for_table(prepareChannels=updateChannels)

        # Initialize all Quad Phases we will be using
        self.initQuadRFValues(values)

        # Append the Quad Phases to Quad tables
        self.prepareMOTTables(updateChannels=updateChannels)  # Load table to QuadRF

        # ---- Arm channels -----
        if armChannels:
            self.armChannels(updateChannels)  # arm channels
        # --- release Toptica lock ----
        if 1 in updateChannels and self.topticaLockWhenUpdating and armChannels: self.setTopticaLock(True)

    def prepare_channels_for_table(self, prepareChannels=(1, 2, 3, 4), reArm=True):
        """
        For the requested channels: arm/disarm and set power limit
        """
        for ch in prepareChannels:
            self.prepare_channel_for_table(ch=ch, reArm=reArm, limit=CHANNELS[ch - 1]['limit'])

    def setTopticaLock(self, flag):
        try:
            if self.topticaController is None: # if controller does not exist
                self.topticaController = TopticaLockController() # create one
                self.topticaController.setLocking(flag) # set locking on/off
            else: # if controller already existed
                self.topticaController.setLocking(flag)  # set locking on/offR
                self.topticaController = None  # destroy the controller, until next change
        except:
            self.logger.warn('Could not connect Toptica controller. Try to reset DigiLock server.')
            return

    # TODO: Remove/Deprecate the below method - we don't need it anymore
    # TODO: We're ensuring the values in freq/amp are only numeric and we cast to string only when we perpare command
    # Make all parameters are in the right units (i.e., unless specified otherwise by user, to be Hz & dBm)
    def changeSequenceParamsToStrings(self):
        for att in self.__dict__.values():  # runs over all @self attributes
            if type(att) == QuadRFPhase:
                for phaseAtt in att.__dict__.values():
                    if type(phaseAtt) == QuadRFChannel:
                        if type(phaseAtt.freq) is not str: phaseAtt.freq = str(phaseAtt.freq) + 'Hz'
                        if type(phaseAtt.amp) is not str: phaseAtt.amp = str(phaseAtt.amp) + 'dbm'

    # send a constant QuadRFPhase.ch to QuadRF
    def sendCmdConstantPhase(self, ch, phase, phase_name=None):
        assert type(phase) == QuadRFPhase # make sure we get an actual well-built phase
        channel = phase.channels[ch]
        cmd = f'TABLE,APPEND,{ch}, {channel.freq}, {channel.amp}, 0x0, {phase.duration}'
        self.sendCmd(cmd)
        # chStr = 'ch_{}'.format(ch)  # e.g., 'ch_1', 'ch_2', etc.
        # self.sendCmd('TABLE,APPEND,%d, %s, %s, 0x0, %s' % (ch, phase[chStr].freq, phase[chStr].amp, phase.duration), phase_name=phase_name)

    # send a dynamic QuadRFPhase.ch to QuadRF
    def sendCmdDynamicPhase(self, ch, startPhase, endPhase, duration, linearPowerChange = True):
        assert type(startPhase) == QuadRFPhase and type(endPhase) == QuadRFPhase  # make sure we get an actual well-built phases
        chStr = 'ch_{}'.format(ch)  # e.g., 'ch_1', 'ch_2', etc.
        self.sendCmdarEnvelopePulse(ch, startPhase[chStr].freq, endPhase[chStr].freq, startPhase[chStr].amp, endPhase[chStr].amp, duration, linearPowerChange = linearPowerChange)  # PGC changing amp & freq stage; from MOT to PGC

    def appendMOTSequence(self, ch):
        #self.changeSequenceParamsToStrings()  # First, make all parameters are in the right units (i.e., unless specified otherwise by user, to be Hz & dBm)
        ch = int(ch)

        # Iterate over all phases
        for phase_index in range(0, len(self.phases_order)):
            # Get the relevant phase according to experiment specific order
            i = (phase_index + self.trigger_on_phase) % len(self.phases_order)
            phase = self.phases_order[i]
            # Send the relevant command
            if phase == Phases.MOT:
                # -- MOT--
                if self.MOT.duration != '0.00ms':
                    # 2 seconds of @f0 in constant @a0
                    self.sendCmdConstantPhase(ch, self.MOT, phase_name='MOT')
                    if self.Post_MOT_delay.duration != '0.00ms':
                        self.sendCmdConstantPhase(ch, self.Post_MOT_delay) # PGC off stage - everything is off, except for repump
            elif phase == Phases.PGC:
                # -- PGC --
                if self.PGC_prep.duration != '0.00ms':
                    # PGC preparation: 3ms of frequency change from f0 to f1, power dropping from a0 to a1
                    self.sendCmdDynamicPhase(ch, startPhase=self.PGC_prep, endPhase=self.PGC, duration=self.PGC_prep.duration, linearPowerChange= False)  # PGC changing amp & freq stage; from MOT to PGC
                    # linearPowerChange = False --> linear drop in dbm --> exp drop in power
                if self.PGC.duration != '0.00ms':
                    # PGC: Then 9 ms of constant a1 & f1
                    self.sendCmdConstantPhase(ch, self.PGC, phase_name = 'PGC')  # PGC constant stage
            elif phase == Phases.FOUNTAIN:
                # Fountain
                if self.Fountain_prep.duration != '0.00ms':
                    # Fountain preparation: 1ms of frequency change from f1 to f2, power dropping from a0 to a1
                    self.sendCmdDynamicPhase(ch, startPhase=self.Fountain_prep, endPhase=self.Fountain, duration=self.Fountain_prep.duration)  # Fountain changing amp & freq stage; PGC to Fountain
                if self.Fountain.duration != '0.00ms':
                    # Fountain: 0.2ms of constant f2 & a2
                    self.sendCmdConstantPhase(ch, self.Fountain, phase_name = 'Fountain')  # Fountain constant stage
            elif phase == Phases.FREE_FALL:
                # Free Fall
                if self.Free_Fall.duration != '0.00ms':
                    self.sendCmdConstantPhase(ch, self.Free_Fall, phase_name='Free_Fall')  # Free-Fall constant stage
            elif phase == Phases.PULSE_1:
                # Pulse 1
                if self.Pulse_1.duration != '0.00ms':
                    if self.Pulse_1_decay.duration != '0.00ms':
                        self.sendCmdDynamicPhase(ch, startPhase=self.Pulse_1_decay, endPhase=self.Pulse_1, duration=self.Pulse_1_decay.duration)  # Decay from Pulse_1_decay amp & freq to Pulse_1 amp & freq
                    self.sendCmdConstantPhase(ch,self.Pulse_1, phase_name='Pulse_1')  # Pulse 1
            elif phase == Phases.INTER_PULSES:
                if self.Inter_Pulses.duration != '0.00ms':
                    # Inter-Pulses
                    self.sendCmdConstantPhase(ch, self.Inter_Pulses)  # In-between pulses; Note amp&freq should be same as post-pulses (off)
            elif i == Phases.PULSE_2:
                if self.Pulse_2.duration != '0.00ms':
                    # Pulse 2
                    self.sendCmdConstantPhase(ch, self.Pulse_2, phase_name='Pulse_2')  # Pulse 2
            elif phase == Phases.POST_PULSE:
                # Post pulses
                self.sendCmdConstantPhase(ch, self.PostPulses, phase_name='Post_Pulse')  # Post-Pulses pulses; Note amp&freq should be the same as inter-pulses (off)
            else:
                self.logger.warn(f'Unknown phase {phase}')

        pass
    def prepareMOTTables(self, updateChannels=(1,2,3,4)):
        for channel in updateChannels:
            self.appendMOTSequence(channel)
