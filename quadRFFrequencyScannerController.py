from quadRFController import QuadRFController
from topticaLockController import TopticaLockController
from Config_Table import Frequency_Scan_Config_Values as config
import numpy as np
# from Config_Table import Phases_Names as Phases_Names
# import time


class QuadRFFrequencyScannerController(QuadRFController):
    def __init__(self, MOGdevice = None, devPort = '169.254.231.196', initialValues = None, channel = None, armChannels = True, opticalPowerCalibration = None, topticaLockWhenUpdating = False ,debugging = False):
        if channel is None or channel not in (1,2,3,4):
            self.printError('Did not provide a valid channel for QuadRFFequencyScannerController.')
            return
        else:
            self.channel = channel
        self.initialValues = config if initialValues is None else initialValues
        # Hand calibration of the power of the double-pass AOM (TOP1)
        if opticalPowerCalibration is None and 'opticalPowerCalibrationFunction' in self.initialValues: opticalPowerCalibration = self.initialValues['opticalPowerCalibrationFunction']

        super(QuadRFFrequencyScannerController, self).__init__(MOGdevice, devPort, opticalPowerCalibration = opticalPowerCalibration, debugging=debugging)


        # Amplitude
        self.Amp = 17.25#15.25  # dbm (15.25 -> red AOM connected to 2W amp, +++++with 15db attenuator on input++++.)
        self.zeroAmp = '0x0'

        self.uploadFrequencyScanTable(values = self.initialValues, channel = self.channel, armChannels = armChannels)

    def __getitem__(self, item): # This enables calling attributs as such: self['attr']  = self.attr
        return getattr(self, item)


    def uploadFrequencyScanTable(self, values = None, channel = None, armChannels = True):
        if values is None: values = self.initialValues
        if channel is None: channel = self.channel
        self.prepareChannelsForTable(prepareChannels=[channel])  # Deletes all existing tables on QuadRF. Also sets limits [hard coded!] on each channel

        # ------------- Read values  -------------
        startDelay = values['Start_Delay_time']
        sawTooth = values['Saw_Tooth_Scan'] # if true-> jump in scan; if false -> triangle scan
        freqRange = values['Frequency_Range']
        startFreq = values['Start_Frequency']
        if startFreq is None: startFreq = freqRange[0]
        totalScanDuration = values['Total_Scan_Duration']  # [ms] This should be the same as the TOTAL length od the measurement
        NFreqs = values['N_Different_Frequencies']  # This defines the resolution of the scan
        NRepetitions = values['N_Scan_Repetitions']  # Number of time we go through all frequencies during @Total_Scan_Duration

        #------------- Calculate durations  -------------
        stepDuration = totalScanDuration / (NFreqs * NRepetitions)
        if (stepDuration * 1e3) % 5 != 0:
            self.printError('NOTE! step duration of scan MUST be a whole number of 5us. Current step duration:{d}'.format(d=stepDuration))
            return
        # ------------- Calculate frequencies vector, to be scanned -------------
        freqsInScan = np.linspace(freqRange[0], freqRange[1], NFreqs)
        startFreqIndex = np.abs(freqsInScan - startFreq).argmin() # find closest frequency in scan to startFreq (find it's index)
        freqsInScan = np.roll(freqsInScan, (-1) * startFreqIndex) # Make sure we start scan at @startFreq

        # ------------- Load table to QuadRF -------------
        if startDelay > 0: self.sendCmd('TABLE,APPEND,{channel}, 110000000, 0x0, 0x0, {startDelay}ms'.format(channel=int(channel), startDelay = startDelay))
        for i in range(NRepetitions):
            # For each repetition...
            for freq in freqsInScan:
                self.sendCmd('TABLE,APPEND,{channel}, {freq}Hz, {amp}dbm, 0x0, {duration}ms'.format(channel = int(channel), freq = freq, amp = self.Amp, duration = stepDuration))
            if not sawTooth: freqsInScan = np.flip(freqsInScan) # flip frequencies order, now scan backwards...

        # Turn off channel until next trigger
        self.sendCmd('TABLE,APPEND,{channel}, 110000000, 0x0, 0x0, 0.01ms'.format(channel=int(channel))) # note this is a different command. this is FORCING QuadRF to be off

        # Arm channel and wait for trigger!
        if armChannels: self.armChannels([channel])             # arm channels

    # Note: here limits are set.
    def prepareChannelsForTable(self, prepareChannels = (1,2,3,4),reArm = True):
        if 1 in prepareChannels: self.prepareChannelForTable(ch=1, reArm=reArm, limit='20dbm')
        if 2 in prepareChannels: self.prepareChannelForTable(ch=2, reArm=reArm, limit='20dbm')
        if 3 in prepareChannels: self.prepareChannelForTable(ch=3, reArm=reArm, limit='32dbm')
        if 4 in prepareChannels: self.prepareChannelForTable(ch=4, reArm=reArm, limit='4dbm')



