from quadRFController import QuadRFController
from topticaLockController import TopticaLockController
from Config_Table import Initial_Values as config
import time


# Helper classes: each phase includes duration, and amplitude & frequency for TOP1 & AOMs
class QuadRFChannel:
    def __init__(self, freq = '0x0', amp = '0x0'):
        self.freq = freq
        self.amp = amp


class QuadRFPhase:
    def __init__(self, duration = '0x0', initial_values = (('0x0','0x0'),('0x0','0x0'),('0x0','0x0'),('0x0','0x0'))):
        if type(duration) is not str:
            duration = str(duration) +'ms'
        self.duration = duration

        self.ch_1 = QuadRFChannel(initial_values[0][1],initial_values[0][1])  # Toptica 1
        self.ch_2 = QuadRFChannel(initial_values[1][1],initial_values[1][1])    # AOM -/0
        self.ch_3 = QuadRFChannel(initial_values[2][1],initial_values[2][1])    # AOM +
        self.ch_4 = QuadRFChannel(initial_values[3][1],initial_values[3][1])    # Repump (/Depump)


class QuadRFMOTController(QuadRFController):
    def __init__(self, devPort = '169.254.231.196', initialValues = None):
        self.AmplifiedAmp = 31  # dbm
        self.UnAmplifiedAmp = 5.75  # dbm, unAmplified

        self.initialValues = config if initialValues is None else initialValues

        super(QuadRFMOTController, self).__init__(devPort)
        self.topticaController = None
        self.setTopticaLock(False)
        self.prepareAllChannelsForTable() # Also sets limits on each channel
        self.initQuadRFValues()

    def initQuadRFValues(self):
        # MOT
        self.MOT = QuadRFPhase(duration = config['MOT_duration'], initial_values = ((config['TOP1_MOT'],self.AmplifiedAmp),(config['AOM_MOT'], self.UnAmplifiedAmp)
                                                                ,(config['AOM_MOT'],self.AmplifiedAmp),(config['AOM_MOT'], self.UnAmplifiedAmp)))
        # PGC-Prep
        self.PGC_prep = QuadRFPhase(duration = config['PGC_prep_duration']) # prep phase, takes initial / final values from MOT & PGC phases

        # PGC (constant)
        PGC_delta_amp = self.amplitudeMultiplierToDBm(config['PGC_final_amp'])
        self.PGC = QuadRFPhase(duration=config['PGC_duration'] - config['PGC_prep_duration'],
                               initial_values=((config['TOP1_PGC'], self.AmplifiedAmp), (config['AOM_MOT'], self.UnAmplifiedAmp + PGC_delta_amp), (config['AOM_MOT'], self.AmplifiedAmp + PGC_delta_amp),
                                               (config['Repump_PGC_freq'], self.UnAmplifiedAmp + PGC_delta_amp)))
        # Fountain-Prep
        self.Fountain_prep = QuadRFPhase(duration = config['Prep_duration']) # prep phase, takes initial / final values from PGC & Fountain phases
        # Fountain (constant)
        Fountain_delta_amp_0 = self.amplitudeMultiplierToDBm(config['Fountain_amp_0'])
        Fountain_delta_amp_plus = self.amplitudeMultiplierToDBm(config['Fountain_amp_plus'])
        self.Fountain = QuadRFPhase(duration=config['Fountain_duration'] - config['Prep_duration'],
                               initial_values=((config['TOP1_PGC'], self.AmplifiedAmp),
                                               (config['AOM_MOT'] - config['Delta_freq'], self.UnAmplifiedAmp + Fountain_delta_amp_0),
                                               (config['AOM_MOT'] + config['Delta_freq'], self.AmplifiedAmp + Fountain_delta_amp_plus),
                                               (config['Repump_PGC_freq'], self.UnAmplifiedAmp + PGC_delta_amp)))
        # Free-fall
        self.FreeFall = QuadRFPhase(duration = config['Snap_Time'], initial_values = ((config['TOP1_MOT'],self.AmplifiedAmp),(config['AOM_MOT'] + config['AOM_Off_Detuning'], self.UnAmplifiedAmp)
                                                                ,(config['AOM_MOT'] + config['AOM_Off_Detuning'],self.AmplifiedAmp),(config['AOM_MOT'], self.UnAmplifiedAmp)))
        # Flash
        self.Flash = QuadRFPhase(duration = config['Flash_duration'], initial_values=((config['TOP1_MOT'], self.AmplifiedAmp), (config['AOM_MOT'], self.UnAmplifiedAmp),
                                                                               (config['AOM_MOT'], self.AmplifiedAmp), (config['AOM_MOT'], self.UnAmplifiedAmp)))

        # Post-Flash
        self.PostFlash = QuadRFPhase(duration=config['Fountain_duration'] - config['Snap_Time'] - config['Flash_duration'],
                                    initial_values=((config['TOP1_MOT'], self.AmplifiedAmp), (config['AOM_MOT'] + config['AOM_Off_Detuning'], self.UnAmplifiedAmp),
                                                    (config['AOM_MOT'] + config['AOM_Off_Detuning'], self.AmplifiedAmp), (config['AOM_MOT'], self.UnAmplifiedAmp)))
    def prepareAllChannelsForTable(self):
        self.prepareChannelForTable(ch = 1,reArm = True, limit = '32dbm')
        self.prepareChannelForTable(ch=2, reArm=True, limit='6dbm')
        self.prepareChannelForTable(ch=3, reArm=True, limit='32dbm')
        self.prepareChannelForTable(ch=4, reArm=True, limit='6dbm')

    def setTopticaLock(self, flag):
        if self.topticaController is None: # if controller does not exist
            self.topticaController = TopticaLockController() # create one
            self.topticaController.setLocking(flag) # set locking on/off
        else: # if controller already existed
            self.topticaController.setLocking(flag)  # set locking on/off
            self.topticaController = None  # destroy the controller, until next change

    def appendMOTSequence(self, ch,f0,a0,f1,a1,f2,a2,f3,a3,f4,a4,f5,a5):
        frqs = [f0, f1, f2, f3, f4, f5]
        amps = [a0, a1, a2, a3, a4, a5]
        ch = int(ch)
        print('Appending MOT table for channel %d.' % ch)

        # Check input...
        if None in frqs + amps:
            print('All values must be provided to appendMOT method...')
            return()
        for i, f in enumerate(frqs):
            if type(f) is not str:
                frqs[i] = str(f) + 'MHz'
        for i, a in enumerate(amps):
            if type(a) is not str:
                amps[i] = str(a) + 'dbm'

        # -- MOT--
        # 2 seconds of @f0 in constant @a0
        self.dev.cmd('TABLE,APPEND,%d, %s, %s, 0x0, %s' % (ch,frqs[0],amps[0],self.MOT.duration))

        # -- PGC --
        # PGC preparation: 3ms of frequency change from f0 to f1, power dropping from a0 to a1
        self.appendLinearEnvelopePulse(ch,frqs[0],frqs[1],amps[0],amps[1],self.PGC_prep.duration) # PGC changing amp & freq stage
        # PGC: Then 9 ms of constant a1 & f1
        self.dev.cmd('TABLE,APPEND,%d, %s, %s, 0x0, %s' %  (ch, frqs[1],amps[1], self.PGC.duration))   # PGC constant stage

        # -- Fountain --
        # Fountain preparation: 1ms of frequency change from f1 to f2, power dropping from a0 to a1
        self.appendLinearEnvelopePulse(ch,frqs[1],frqs[2],amps[1],amps[2],self.Fountain_prep.duration) # Fountain changing amp & freq stage
        # Fountain: 0.2ms of constant f2 & a2
        self.dev.cmd('TABLE,APPEND,%d, %s, %s, 0x0,%s' % (ch, frqs[2],amps[2], self.Fountain.duration))   # Fountain constant stage

        ## -- Free Fall --
        #TODO 0x0 means until further notice
        self.dev.cmd('TABLE,APPEND,%d, %s, %s, 0x0,%s' % (ch, frqs[3], amps[3], self.FreeFall.duration))  # Free-Fall constant stage

        ## -- Flash!  --
        self.dev.cmd('TABLE,APPEND,%d, %s, %s, 0x0,%s' % (ch, frqs[4], amps[4], self.Flash.duration))  # Flash

        ## -- Post-Flash  --
        self.dev.cmd('TABLE,APPEND,%d, %s, %s, 0x0,%s' % (ch, frqs[5], amps[5], self.PostFlash.duration))  # Post Flash constant stage


    def prepareMOTTablesForAllChannels(self):
        # Note: for each channel there exist 4 (constant) stages, hence 4 pairs of freq& amp
        # If numbers are given (instead of strings), dbm & MHz are assumed for units
        # Between constant stages, there are preparation stages: linearEnvelope taking the parameters from one stage to the other.
        # Exception: 3rd to 4th stage - the change between the two is not sloped

        # Top 1 Lock
        self.appendMOTSequence(1, self.MOT.ch_1.freq, self.MOT.ch_1.amp, self.PGC.ch_1.freq, self.PGC.ch_1.amp, self.Fountain.ch_1.freq,
                               self.Fountain.ch_1.amp, self.FreeFall.ch_1.freq, self.FreeFall.ch_1.amp, self.Flash.ch_1.freq, self.Flash.ch_1.amp, self.PostFlash.ch_1.freq, self.PostFlash.ch_1.amp)
        # MOT 0 / -
        self.appendMOTSequence(2, self.MOT.ch_2.freq, self.MOT.ch_2.amp, self.PGC.ch_2.freq, self.PGC.ch_2.amp, self.Fountain.ch_2.freq,
                               self.Fountain.ch_2.amp, self.FreeFall.ch_2.freq, self.FreeFall.ch_2.amp, self.Flash.ch_2.freq, self.Flash.ch_2.amp, self.PostFlash.ch_2.freq, self.PostFlash.ch_2.amp)
         #MOT +
        self.appendMOTSequence(3, self.MOT.ch_3.freq, self.MOT.ch_3.amp, self.PGC.ch_3.freq, self.PGC.ch_3.amp, self.Fountain.ch_3.freq,
                               self.Fountain.ch_3.amp, self.FreeFall.ch_3.freq, self.FreeFall.ch_3.amp, self.Flash.ch_3.freq, self.Flash.ch_3.amp, self.PostFlash.ch_3.freq, self.PostFlash.ch_3.amp)
        # Re-pump / De-pump
        self.appendMOTSequence(4, self.MOT.ch_4.freq, self.MOT.ch_4.amp, self.PGC.ch_4.freq, self.PGC.ch_4.amp, self.Fountain.ch_4.freq,
                               self.Fountain.ch_4.amp, self.FreeFall.ch_4.freq, self.FreeFall.ch_4.amp, self.Flash.ch_4.freq, self.Flash.ch_4.amp, self.PostFlash.ch_4.freq, self.PostFlash.ch_4.amp)

        # AOMAmp = 31  # dbm
        # AOMUnAmp = 5.75  # dbm, unAmplified
        # deTune = 0.3  # MHz ???
        # # Top 1 Lock
        # self.appendMOTSequence(1, '113MHz', AOMAmp, '93MHz', AOMAmp, '93MHz', AOMAmp, '113MHz', AOMAmp)
        # # MOT 0 / -
        # self.appendMOTSequence(2, '110MHz', AOMUnAmp, '110MHz', AOMUnAmp , 110 + deTune, AOMUnAmp , 110 + deTune + 30,AOMUnAmp)
        #  #MOT +
        # self.appendMOTSequence(3, '110MHz', AOMAmp, '110MHz', AOMAmp - 0* 0.08, 110 - deTune, AOMAmp,110 - deTune - 30, AOMAmp)
        # # Re-pump / De-pump
        # self.appendMOTSequence(4, '78.4735MHz', AOMUnAmp, '78.4735MHz', AOMUnAmp, '78.4735MHz', AOMUnAmp, '138.4735MHz',AOMUnAmp)

    def continuousTableForChannel(self, ch, frq, amp):
        self.dev.cmd('TABLE,APPEND,%d, %s, %s, 0x0, 0x0' % (ch,frq,amp))

    def continuousTablesForAllChannels(self):
        self.continuousTableForChannel(1, '113MHz', AOMAmp)
        self.continuousTableForChannel(2, '110MHz', AOMUnAmp)
        self.continuousTableForChannel(3, '110MHz', AOMAmp)
        self.continuousTableForChannel(4, '78.4735MHz', AOMUnAmp)

initTime = time.time()

dev = QuadRFMOTController('169.254.231.196')

#dev.continuousTablesForAllChannels()
dev.prepareMOTTablesForAllChannels()
dev.armChannel(1)
dev.armChannel(2)
dev.armChannel(3)
dev.armChannel(4)
#
dev.setTopticaLock(True)
finalTime = time.time()
print('Load time [sec]: ', (finalTime-initTime))
#### TEST ####

#
# # TODO - This is for debugging!
# dev.dev.cmd('TABLE,RESTART,1,on')
# dev.dev.cmd('TABLE,RESTART,2,on')
# dev.dev.cmd('TABLE,RESTART,3,on')
# dev.dev.cmd('TABLE,RESTART,4,on')

#dev.dev.cmd('TABLE,START,1,2,3,4')
#
#
# #dev.appendLinearEnvelopePulse(1,'10MHz','1MHz',12,-2,'5sec',50)

