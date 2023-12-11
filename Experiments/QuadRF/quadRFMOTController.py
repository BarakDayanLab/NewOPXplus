from topticaLockController import TopticaLockController
from Experiments.BaseExperiment.Config_Table import Initial_Values as config
from Experiments.BaseExperiment.Phases import Phases
from Experiments.BaseExperiment.Phases import Phases_Names
from Experiments.QuadRF.QuadRFController import QuadRFController
from Utilities.BDLogger import BDLogger
from Utilities.Utils import Utils


# Helper classes: each phase includes duration, and amplitude & frequency for TOP1 & AOMs
class QuadRFChannel:
    def __init__(self, freq='0x0', amp='0x0'):
        self.freq = freq
        self.amp = amp


class QuadRFPhase:
    def __init__(self, duration='0x0', initial_values=(('0x0', '0x0'), ('0x0', '0x0'),('0x0', '0x0'), ('0x0', '0x0'))):
        if type(duration) is not str:
            duration = '{:.2f}ms'.format(duration)  # 2 decimal points, as string
        self.duration = duration

        self.ch_1 = QuadRFChannel(initial_values[0][0], initial_values[0][1])  # Toptica 1
        self.ch_2 = QuadRFChannel(initial_values[1][0], initial_values[1][1])  # AOM -/0
        self.ch_3 = QuadRFChannel(initial_values[2][0], initial_values[2][1])  # AOM +
        self.ch_4 = QuadRFChannel(initial_values[3][0], initial_values[3][1])  # Repump (/Depump)

    def __getitem__(self, item): # This enables calling attributs as such: self['attr']  = self.attr
        return getattr(self, item)


class QuadRFMOTController(QuadRFController):
    def __init__(self, MOGdevice=None, devPort='169.254.231.196', initialValues=None, updateChannels=(1,2,3,4), armChannels=True, opticalPowerCalibration=None, topticaLockWhenUpdating=False, debugging=False, continuous=False):

        self.bdlogger = BDLogger()

        if opticalPowerCalibration is None and 'opticalPowerCalibrationFunction' in initialValues:
            opticalPowerCalibration = initialValues['opticalPowerCalibrationFunction']
        super(QuadRFMOTController, self).__init__(MOGdevice, devPort, opticalPowerCalibration=opticalPowerCalibration, debugging=debugging)
        self.initialValues = config if initialValues is None else initialValues

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
        self.continuous = continuous

        # TODO: remove this once we finish debugging and setting all experiments
        if 'Operation_Mode' in self.initialValues and self.initialValues['Operation_Mode'] == 'Off':
            raise Exception("**** Note - there is code that still uses Operation_Mode['off']. Please check it! *****")

        if continuous:  # run continuous mode then quit
            self.continuousTablesForChannels(updateChannels)
            return
        #elif self.initialValues['Operation_Mode'] == 'Off':
        #    self.turnChannelsOff(holdLocking=True)
        #    return
        else:
            self.uploadMOTTables(values=self.initialValues, updateChannels=updateChannels, armChannels=armChannels)

    def __getitem__(self, item):  # This enables calling attributes as such: self['attr']  = self.attr
        return getattr(self, item)

    def initQuadRFValues(self, values=None):
        if values is None:
            values = self.initialValues

        self.phases_order = values['Phases_Order']
        self.trigger_on_phase = values['Triggering_Phase']

        AOMOffFreq = values['MOT_AOM_freq'] + values['AOM_Off_Detuning']

        # -----------------  MOT ------------------------
        self.MOT = QuadRFPhase(duration=values['MOT_duration'],
                               initial_values=((values['MOT_freq'], self.Amp_Ch1), (AOMOffFreq, self.zeroAmp), (values['MOT_AOM_freq'], self.Amp_Ch3), (values['AOM_Repump_freq'], self.Amp_Ch4)))
        # MOT Delay (turn everything off, except for repump)
        self.Post_MOT_delay = QuadRFPhase(duration=values['Post_MOT_delay'], initial_values=((values['MOT_freq'], self.Amp_Ch1), (AOMOffFreq, self.zeroAmp),
                                                                                             (AOMOffFreq, self.zeroAmp), (values['AOM_Repump_freq'], self.Amp_Ch4)))
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
                                    initial_values=((values['PGC_initial_freq'], self.Amp_Ch1), (AOMOffFreq, self.zeroAmp),
                                                    (values['MOT_AOM_freq'] - values['PGC_initial_Delta_freq'], self.Amp_Ch3 + PGC_initial_delta_amp_plus), (values['Repump_PGC_freq'], self.Amp_Ch4)))

        # PGC (constant)
        self.PGC = QuadRFPhase(duration=PGC_duration,
                               initial_values=((values['PGC_final_freq'], self.Amp_Ch1 + PGC_final_amp_delta), (AOMOffFreq, self.zeroAmp),
                                               (values['MOT_AOM_freq'] - values['PGC_final_Delta_freq'], self.Amp_Ch3 + PGC_final_delta_amp_plus), (values['Repump_PGC_freq'], self.Amp_Ch4)))
        # -----------------  Fountain ------------------------
        Fountain_duration = values['Fountain_duration'] - values['Fountain_prep_duration'] if values['Fountain_duration'] > 0 else 0
        Fountain_prep_duration = values['Fountain_prep_duration'] if values['Fountain_duration'] > 0 else 0

        Fountain_final_delta_amp_0 = Utils.amplitudeMultiplierToDBm(values['Fountain_final_amp_0'])
        Fountain_final_delta_amp_plus = Utils.amplitudeMultiplierToDBm(values['Fountain_final_amp_plus'])
        Fountain_initial_delta_amp_0 = Utils.amplitudeMultiplierToDBm(values['Fountain_initial_amp_0'])
        Fountain_initial_delta_amp_plus = Utils.amplitudeMultiplierToDBm(values['Fountain_initial_amp_plus'])

        # Fountain-Prep (dynamic)
        self.Fountain_prep = QuadRFPhase(duration=Fountain_prep_duration, initial_values=((values['Fountain_initial_freq'], self.Amp_Ch1),
                                                    # (values['MOT_AOM_freq'] + values['Fountain_initial_Delta_freq'], self.Amp_Ch2 + Fountain_initial_delta_amp_0),
                                                    (AOMOffFreq, self.zeroAmp),
                                                    (values['MOT_AOM_freq'] - values['Fountain_initial_Delta_freq'], self.Amp_Ch3 + Fountain_initial_delta_amp_plus),
                                                    (values['Repump_PGC_freq'], self.Amp_Ch4))) # prep phase, takes initial / final values from PGC & Fountain phases
        # Fountain (constant)
        self.Fountain = QuadRFPhase(duration=Fountain_duration,
                                    initial_values=((values['Fountain_final_freq'], self.Amp_Ch1),
                                                    # (values['MOT_AOM_freq'] + values['Fountain_final_Delta_freq'], self.Amp_Ch2 + Fountain_final_delta_amp_0),
                                                    (AOMOffFreq, self.zeroAmp),
                                                    (values['MOT_AOM_freq'] - values['Fountain_final_Delta_freq'], self.Amp_Ch3 + Fountain_final_delta_amp_plus),
                                                    (values['Repump_PGC_freq'], self.Amp_Ch4)))
        # ----------------- Free-fall -----------------
        PrePulse_delta_amp_Repump = float(Utils.amplitudeMultiplierToDBm(values['PrePulse_Repump_amp']))
        self.Free_Fall = QuadRFPhase(duration=values['PrePulse_duration'],
                                     initial_values=((values['PrePulse_CH1_freq'], self.Amp_Ch1),
                                                     (values['PrePulse_CH2_freq'], self.Amp_Ch2),
                                                     (AOMOffFreq, self.Amp_Ch3),
                                                     (values['Pulse_1_CH4_Freq'], self.Amp_Ch4 + PrePulse_delta_amp_Repump)))

        # ----------------- Pulse 1 -----------------
        Pulse_1_delta_amp_i = float(Utils.amplitudeMultiplierToDBm(values['Pulse_1_amp_i']))
        Pulse_1_delta_amp_f = float(Utils.amplitudeMultiplierToDBm(values['Pulse_1_amp_f']))
        Pulse_1_delta_amp_Repump = float(Utils.amplitudeMultiplierToDBm(values['Pulse_1_Repump_amp']))
        self.Pulse_1 = QuadRFPhase(duration=values['Pulse_1_duration'] + values['M_off_time'] - values['Pulse_1_decay_duration'],
                                   initial_values=((values['Pulse_1_CH1_Freq_f'], self.Amp_Ch1 + Pulse_1_delta_amp_f), (AOMOffFreq, self.zeroAmp),
                                                   (values['Pulse_1_CH_2_3_Freq'], self.Amp_Ch3 + Pulse_1_delta_amp_f), (values['Pulse_1_CH4_Freq'], self.Amp_Ch4 + Pulse_1_delta_amp_Repump)))

        self.Pulse_1_decay = QuadRFPhase(duration=values['Pulse_1_decay_duration'],
                                         initial_values=((values['Pulse_1_CH1_Freq_i'], self.Amp_Ch1 + Pulse_1_delta_amp_i), (AOMOffFreq, self.zeroAmp),
                                                         (values['Pulse_1_CH_2_3_Freq'], self.Amp_Ch3 + Pulse_1_delta_amp_i), (values['Pulse_1_CH4_Freq'], self.Amp_Ch4 + Pulse_1_delta_amp_Repump)))

        # ----------------- Inter-Pulses -----------------
        self.Inter_Pulses = QuadRFPhase(duration=values['InterPulses_duration'], initial_values=((values['MOT_freq'], self.Amp_Ch1), (AOMOffFreq, self.zeroAmp),
                                                                                                 (AOMOffFreq, self.Amp_Ch3), (AOMOffFreq, self.Amp_Ch4)))
        # ----------------- Pulse 2 -----------------
        self.Pulse_2 = QuadRFPhase(duration=values['Pulse_2_duration'], initial_values=((values['Pulse_2_CH1_Freq'], self.Amp_Ch1), (AOMOffFreq, self.zeroAmp),
                                                                                        (values['Pulse_2_CH_2_3_Freq'], self.Amp_Ch3), (values['Pulse_2_CH4_Freq'], self.Amp_Ch4)))

        # ----------------- Post-Pulses -----------------
        PostPulsesDuration = values['FreeFall_duration'] - values['PrePulse_duration'] - values['Pulse_1_duration'] - values['M_off_time'] - values['Pulse_2_duration'] - values['InterPulses_duration']
        if PostPulsesDuration > 0:
            # self.PostPulses = QuadRFPhase(duration=PostPulsesDuration, initial_values=((values['MOT_freq'], -20), (AOMOffFreq, self.Amp_Ch2),(AOMOffFreq, self.Amp_Ch3), (AOMOffFreq, self.Amp_Ch4)))
            self.PostPulses = QuadRFPhase(duration=PostPulsesDuration, initial_values=((values['MOT_freq'], self.Amp_Ch1), (AOMOffFreq, self.zeroAmp), (AOMOffFreq, self.Amp_Ch3), (values['AOM_Repump_freq'], self.Amp_Ch4)))

        else:
            self.PostPulses = QuadRFPhase(initial_values=((values['MOT_freq'], self.Amp_Ch1), (AOMOffFreq, self.zeroAmp),
                                                          (values['Flash_freq'], self.Amp_Ch3), (values['AOM_Repump_freq'], self.Amp_Ch4)))

        self.changeLastPhaseDuration()

    def changeLastPhaseDuration(self):
        """
        Making sure the last phase duration is 0.01ms - this way this phase goes on until another trigger arrives.
        The previous phase is dependant on the experiments phase order
        """
        d = '0.01ms'

        previous_phase_index = (self.phases_order.index(self.trigger_on_phase)-1) % len(self.phases_order)
        last_phase_index = self.phases_order[previous_phase_index]

        triggeringPhase = Phases_Names[self.trigger_on_phase]
        lastPhase = Phases_Names[last_phase_index]

        if self.debugging:
            self.bdlogger.debug(f'Triggering on {triggeringPhase}, {lastPhase} will continue until trigger.')

        # e.g., if we trigger on Pulse_1, lase phase is FreeFall
        self[lastPhase].duration = d

    def uploadMOTTables(self, values=None, updateChannels=(1,2,3,4), armChannels=True):
        if values is None:
            values = self.initialValues

        # --- Handle Toptica lock if necessary -----
        if 1 in updateChannels and self.topticaLockWhenUpdating:
            self.setTopticaLock(False)

        # --- Prepare and load table ------
        self.prepareChannelsForTable(prepareChannels=updateChannels)  # Deletes all existing tables on QuadRF. Also sets limits [hard coded!] on each channel
        self.initQuadRFValues(values)
        self.prepareMOTTables(updateChannels=updateChannels)  # Load table to QuadRF

        # ---- Arm channels -----
        if armChannels:
            self.armChannels(updateChannels)  # arm channels
        # --- release Toptica lock ----
        if 1 in updateChannels and self.topticaLockWhenUpdating and armChannels: self.setTopticaLock(True)

    def prepareChannelsForTable(self, prepareChannels=(1,2,3,4), reArm=True):
        if 1 in prepareChannels: self.prepareChannelForTable(ch=1, reArm=reArm, limit='22dbm')
        if 2 in prepareChannels: self.prepareChannelForTable(ch=2, reArm=reArm, limit='20dbm')
        if 3 in prepareChannels: self.prepareChannelForTable(ch=3, reArm=reArm, limit='32dbm')
        if 4 in prepareChannels: self.prepareChannelForTable(ch=4, reArm=reArm, limit='4dbm')

    def setTopticaLock(self, flag):
        try:
            if self.topticaController is None: # if controller does not exist
                self.topticaController = TopticaLockController() # create one
                self.topticaController.setLocking(flag) # set locking on/off
            else: # if controller already existed
                self.topticaController.setLocking(flag)  # set locking on/offR
                self.topticaController = None  # destroy the controller, until next change
        except:
            self.bdlogger.warn('Could not connect Toptica controller. Try to reset DigiLock server.')
            return

    # Make all parameters are in the right units (i.e., unless specified otherwise by user, to be Hz & dBm)
    def changeSequenceParamsToStrings(self):
        for att in self.__dict__.values(): # runs over all @self attributes
            if type(att) == QuadRFPhase:
                for phaseAtt in att.__dict__.values():
                    if type(phaseAtt) == QuadRFChannel:
                        if type(phaseAtt.freq) is not str: phaseAtt.freq = str(phaseAtt.freq) + 'Hz'
                        if type(phaseAtt.amp) is not str: phaseAtt.amp = str(phaseAtt.amp) + 'dbm'

    # send a constant QuadRFPhase.ch to QuadRF
    def sendCmdConstantPhase(self, ch, phase, phase_name=None):
        assert type(phase) == QuadRFPhase # make sure we get an actual well-built phase
        chStr = 'ch_{}'.format(ch)  # e.g., 'ch_1', 'ch_2', etc.
        self.sendCmd('TABLE,APPEND,%d, %s, %s, 0x0, %s' % (ch, phase[chStr].freq, phase[chStr].amp, phase.duration), phase_name=phase_name)

    # send a dynamic QuadRFPhase.ch to QuadRF
    def sendCmdDynamicPhase(self, ch, startPhase, endPhase, duration, linearPowerChange = True):
        assert type(startPhase) == QuadRFPhase and type(endPhase) == QuadRFPhase  # make sure we get an actual well-built phases
        chStr = 'ch_{}'.format(ch)  # e.g., 'ch_1', 'ch_2', etc.
        self.sendCmdarEnvelopePulse(ch, startPhase[chStr].freq, endPhase[chStr].freq, startPhase[chStr].amp, endPhase[chStr].amp, duration, linearPowerChange = linearPowerChange)  # PGC changing amp & freq stage; from MOT to PGC

    def appendMOTSequence(self, ch):
        self.changeSequenceParamsToStrings()  # First, make all parameters are in the right units (i.e., unless specified otherwise by user, to be Hz & dBm)
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
                self.bdlogger.warn(f'Unknown phase {phase}')

        pass
    def prepareMOTTables(self, updateChannels=(1,2,3,4)):
        for channel in updateChannels:
            self.appendMOTSequence(channel)
