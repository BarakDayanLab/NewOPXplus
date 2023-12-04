# Written by ngordon1

# pyserial & mogdevice modules for python should be installed.
# pyserial can be installed via pip install
# You can download the file from: https://www.moglabs.com/support/software/connection. Look for python driver; then put it in python-Lib library.

## --------------- Usage example: ----------------------------

# mogRFDevice = QuadRFController('COM3')
# mogRFDevice.armMogRFTable('mogScriptExample.txt')

## --------------- Formatting of commands file: ---------------

# Code ignores empty lines & lines of which 1st char is '#' (could be used for comments)
# Each line should be comprised of 5 data-words:
# channel[1-4] , freq [Hz, kHz, MHz], power [dBm, dB, mW, W], phase [deg, rad], duration [us, ms, s]
# Such that:
#           1,10MHz,10dBm,0x0,3s
# Means:
#           'Operate channel 1 at 10Mhz & 10dBm, 0 phase, for a total of 3 seconds'
# Note: 0x0 means absolute 0.
# So:
#           4,10MHz,4dBm,0x0,0x0 [duration 0x0 = until directed otherwise]
# Means turnning channel 4 indefinitely.
# See MOGLabs QRF manual for further details.


## ------------------------------------------------------------
from Experiments.QuadRF.mogdevice import MOGDevice

import os.path
from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from Utilities.BDLogger import BDLogger
from Utilities.Utils import Utils


class QuadRFController:
    def __init__(self, MOGdevice=None, devPort='169.254.231.196', debugging=False, opticalPowerCalibration=None):

        # Initialize Logger
        self.logger = BDLogger()

        # connect to the device (or use an existing one - passed as argument in 'MOGdevice'
        try:
            self.device = MOGDevice(devPort) if MOGdevice is None else MOGdevice
        except Exception as e:
            self.logger.warn(f'Could not connect to QuadRF device-port {devPort}. Error: {e}. Try different port. Usual values are \'COM3\' or some ip address.')
            return
        # print some information
        self.debugging = debugging
        if self.debugging:
            device_info = self.device.ask('info')
            self.logger.debug(f'Connected to QuadRF. Device info: {device_info}')
        self.lineCount = 0
        self.lines = [{'Amplitudes': [], 'Durations':[], 'Frequencies':[]}, {'Amplitudes': [], 'Durations':[], 'Frequencies':[]}, {'Amplitudes': [], 'Durations':[], 'Frequencies':[]}, {'Amplitudes': [], 'Durations':[], 'Frequencies':[]}]  # 4 channels, each holds: duration, frequency, amplitude
        self.constantOpticalPowerCalibration = opticalPowerCalibration  # Format: c = {'ch_1':{'calibration_func':calibration_func, 'calibrate_RF_Power':'30dbm','calibrate_Optical_Power':0.5}}

        # Previous line means: for channel 1, use calibration function to calibrate 30dbm to be 0.5 power [normalized] (for ANY frequency)
        # Note: calibration_func should be a python function taking two vars: freq & optical power, returning the desired RF power
        if self.constantOpticalPowerCalibration is not None:
            self.logger.info('Note: calibration RF to optical power is ON')

    def __del__(self):
        self.disconnectQuadRF()

    def disconnectQuadRF(self):
        try:
            self.device.close()
        except Exception as e:
            self.logger.error(f'Unable to disconnect from QuadRF. {e}. It was probably already disconnected.')

    def sendCmd(self, line, phase_name=None):
        if str(line).find('TABLE,APPEND') >= 0:
            self.lineCount += 1
            terms = self.processTableAppendLine(line)
            ch = terms['Channel']
            if ch in range(1, 5):
                if self.constantOpticalPowerCalibration and ('ch_%d' % ch) in self.constantOpticalPowerCalibration:
                    cDic = self.constantOpticalPowerCalibration[('ch_%d' % ch)]
                    if ('affected_phases_names' not in cDic or phase_name in cDic['affected_phases_names']) and terms['Amplitude'] != '0x0':
                        if 'calibrated_to_rf_ratio' not in cDic: cDic['calibrated_to_rf_ratio'] = False
                        terms['Amplitude'] = self.calibrateAmplitudeForConstantOpticalPower(cDic['calibration_func'], terms['Frequency'], terms['Amplitude'], float(self.getNumericValueFromString(cDic['calibrate_RF_Power'])), cDic['calibrate_Optical_Power'], calibrated_to_rf_ratio = cDic['calibrated_to_rf_ratio'])
                    line = self.processDictionaryToTableAppendLine(terms)
                self.lines[ch - 1]['Durations'].append(terms['Duration'])
                self.lines[ch - 1]['Frequencies'].append(terms['Frequency'])
                self.lines[ch - 1]['Amplitudes'].append(terms['Amplitude'])
            if self.debugging:
                self.logger.debug(line)
        self.device.cmd(line)

    def getNumericValueFromString(self, s):
        res = ''.join([n for n in str(s) if (n in['-','.'] or n.isdigit())]) if '0x0' not in s else '0x0' # '0x0' IS numeric, it means absolute zero
        return res

    def lineIsTableCommand(self, l):
        if  type(l) is str and len(l.strip()) > 0 and l[0] != '#':
            return True
        else:
            return False

    # @initalAmp, @finalAmp are numbers (dBm are assumed) or strings
    # @ch is channel int, @freq could be either number (MHz assumed) or string.
    # @duration is total duration of pulse [ms assumed]
    # def sendCmdarEnvelopePulse(self, ch, initialFreq, finalFreq, initialAmp, finalAmp, duration, steps = 200):
    def sendCmdarEnvelopePulse(self, ch, initialFreq, finalFreq, initialAmp, finalAmp, duration, linearPowerChange = True, steps = 50):
        # Check: if we have a phase which should be constant, no need to send an envelope.
        if initialFreq == finalFreq and initialAmp == finalAmp:
            self.sendCmd('TABLE,APPEND,%d, %s, %s, 0x0, %s' % (ch, initialFreq, initialAmp, duration))
            return
        iAmp = self.getNumericValueFromString(initialAmp) # string
        fAmp = self.getNumericValueFromString(finalAmp) # string
        ampScale = str(initialAmp).replace(str(iAmp),'')
        if ampScale == '':
            ampScale = 'dBm'

        iFreq = self.getNumericValueFromString(initialFreq)# string
        fFreq = self.getNumericValueFromString(finalFreq)# string
        freqScale = str(initialFreq).replace(str(iFreq),'')
        if freqScale == '':
            freqScale = 'MHz'

        d = self.getNumericValueFromString(duration)
        durationScale = str(duration).replace(str(d), '')
        if durationScale == '':
            durationScale = 'ms'

        if linearPowerChange:
            amps = [Utils.amplitudeMultiplierToDBm(p) for p in np.linspace(Utils.dbmToP(float(iAmp)), Utils.dbmToP(float(fAmp)), num = steps)] # Linear drop in power
        else:
            amps = np.linspace(float(iAmp), float(fAmp), num=steps)  # Linear drop in dbms
        freq = np.linspace(float(iFreq),float(fFreq),steps)
        t = float(d) / steps

        # Phase calculations
        tScale = 0.001 if durationScale == 'ms' else 1  # 1 if [sec]
        fScale = 1e6 if freqScale == 'MHz' else 1  # 1 if [Hz]
        phase = 0  # in degrees (!)
        for index, a in enumerate(amps):
            self.sendCmd('TABLE,APPEND,%d,%.2f%s,%.2f%s,%.2f, %.2f%s' % (ch, float(freq[index]), freqScale, float(a), ampScale, float(phase), float(t), durationScale))
            phase = (phase + freq[index] * fScale * t * tScale) % 360

    def calibrateAmplitudeForConstantOpticalPower(self, calibration_func, freq, RFAmp, calibrate_RF_Power = 32, calibrate_Optical_Power = 0.6, calibrated_to_rf_ratio = False):
        ratio = Utils.dbmToP(RFAmp - calibrate_RF_Power)  # i.e., say I want the equivalent of 22dbm, and am calibrated to 32dbm, then optical power is 10 times weaker
        targetOpticalPower = calibrate_Optical_Power * ratio
        calibratedRF = float(calibration_func(freq, targetOpticalPower))
        if calibratedRF == 0 or np.isnan(calibratedRF):
            self.logger.warn('For f = {f}, RFAmp = {a}, calibrated RF power is 0 (perhaps input outside calibration function limits?)'.format(f=freq, a=RFAmp))
            calibratedRF = calibrate_RF_Power
            self.logger.warn('Setting output RF to maximum ({m}))'.format(m=calibrate_RF_Power))
        if calibrated_to_rf_ratio: calibratedRF = RFAmp + Utils.amplitudeMultiplierToDBm(calibratedRF)  # in this case, calibration function actually gives ratio...
        return calibratedRF

    # ----------------------------------------------------------
    # ----- Arming channels and setting limits ------
    # ----------------------------------------------------------
    # @txtfile: location for table textfile.
    # @reArm: activate rearming option for all chans (see QuadRF manual)
    # @checkCommandsFile: check @txtfile before loading; If it isn't properly constructed, no change will occour
    # @trigger: either 'RISING' or 'FALLING' (see quadRF manual)
    def armMogRFTable(self, txtfile=None, reArm=True, checkCommandsFile=True, trigger='RISING'):
        if not os.path.isfile(txtfile):
            self.logger.warn(f'Could not find {txtfile}. Check for file location.')

        lines = []
        with open(txtfile) as f:
            lines = f.readlines()
            chns = [True] * 4

            # Check commands file
            if checkCommandsFile:
                chns = [False] * 4
                for i, l in enumerate(lines):
                    if self.lineIsTableCommand(l):
                        c = int(l.split(',')[0])
                        try:
                            chns[c - 1] = True
                        except:
                            self.logger.warn('Line %d is corrupted' % (i + 1))
                            self.logger.warn('Table was not loaded.')
                            return

            # Load general commands
            for i, c in enumerate(chns):
                if c is True:
                    self.prepareChannelForTable(i + 1, reArm, trigger)

            # Load table command
            for l in lines:
                if self.lineIsTableCommand(l):
                    self.sendCmd('TABLE,APPEND,%s' % l)
            # ARM
            for i, c in enumerate(chns):
                if c is True:
                    self.armChannel(i + 1)

    def prepareChannelForTable(self, ch, reArm=True, trigger='RISING', limit='4dbm'):
        ch = int(ch)
        self.setLimit(ch, limit)  # Current limit on channel
        self.sendCmd('MODE,%d,TSB' % ch)  # Out in table mode
        self.sendCmd('TABLE,CLEAR,%d' % ch)  # clear all previous commands in table
        self.lineCount = 0
        self.sendCmd('TABLE,EDGE,%d,%s' % (ch, trigger.upper()))  # Set trigger option
        if reArm:
            self.sendCmd('TABLE,REARM,%d,on' % (ch))
        else:
            self.sendCmd('TABLE,REARM,%d,off' % (ch))

    def armChannel(self, ch):
        self.sendCmd('TABLE,ARM,%d' % (ch))

    def armChannels(self, chns=(1, 2, 3, 4)):
        for ch in chns:
            self.armChannel(ch)

    def setLimit(self, ch, limit):
        self.sendCmd('LIMIT, %s, %s' % (str(ch), str(limit)))

    # ----------------------------------------------------------
    # ------ Methods that have to do with debugging and logging
    # ----------------------------------------------------------
    def processDictionaryToTableAppendLine(self, dic): # assume input in Hz, dbm, ms
        ampUnit = 'dbm' if dic['Amplitude'] != '0x0' else '' # This is a patch. there must be a more elegant way.
        return 'TABLE,APPEND,{ch}, {freq}Hz, {amp}{ampUnit}, 0x0, {duration}ms'.format(ch = dic['Channel'], amp = dic['Amplitude'],ampUnit = ampUnit, freq = dic['Frequency'], duration = dic['Duration'])

    def processTableAppendLine(self, l):
        # line example:  'TABLE,APPEND,1, 113MHz, 32db, 0x0, 1ms'
        terms = l.split(',')
        durationNV = self.getNumericValueFromString(str(terms[-1]))
        if durationNV != '0x0':
            duration = float(durationNV) if str(terms[-1]).replace(durationNV,'').strip().lower() == 'ms' else float(durationNV) * 1e3  # ms or sec
        else:
            duration = durationNV        #if duration > 1e3: duration = duration / 1e2
        freqNV = self.getNumericValueFromString(str(terms[3]))
        freq = float(freqNV) / 1e6 if str(terms[-1]).replace(freqNV,'').strip().lower() == 'hz' else float(freqNV)  # either hz or mhz
        ampNV = self.getNumericValueFromString(str(terms[4].replace('dbm', '')))
        amp = float(ampNV) if ampNV != '0x0' else '0x0'                                                         # always dbm
        ch = int(terms[2])
        return {'Duration': duration, 'Frequency': freq, 'Amplitude': amp, 'Channel': ch}

    def get_channel_signals_as_string(self, channel_number):
        """
        Returns all channels signals in a CSV line format (Amplitude, Durations, Frequencies)
        """
        channel_signals = self.lines[channel_number]
        csv_line = ''
        for channel in channel_signals:
            if channel != {'Amplitudes': [], 'Durations': [], 'Frequencies': []}:
                csv_line = f"{channel['Amplitude']},{channel['Durations']},{channel['Frequencies']}\n"
        return csv_line

    def saveLinesAsCSV(self, path):
        for i, channel in enumerate(self.lines):
            if channel != {'Amplitudes': [], 'Durations': [], 'Frequencies': []}:  # That is, if channel isn't empty
                file_name = Path(path).stem
                save_path = path.replace(file_name, f'{file_name}_ch{i + 1}')  # add channel to save path
                try:
                    with open(save_path, "w") as outfile:
                        writer = csv.writer(outfile)
                        writer.writerow(channel.keys())
                        writer.writerows(zip(*channel.values()))
                except IOError as e:
                    self.logger.error(f'I/O error in saveLinesAsCSV {e}')


    def plotTables(self):
        self.logger.debug(self.lines)
        fig, axs = plt.subplots(2,2)
        fig.suptitle('QuadRF Table [MOT not to scale]')
        for ch, terms in enumerate(self.lines):
            channel = ch + 1
            durations = terms['Durations']
            time, t = [], 0
            for i, k in enumerate(durations):
                time.append(t)
                t += durations[i]
            freqs = terms['Frequencies']
            amps = terms['Amplitudes']
            axs[ch // 2, ch % 2].plot(time, freqs, label = 'Frequency [MHz]')
          #  axs[ch // 2, ch % 2].plot(durations, amps, label='Amplitude [dBm]')
            axs[ch // 2, ch % 2].set_title(f'Channel {channel:d}')

        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.show()

    # ----------------------------------------------------------
    # --------------------- Continuous mode methods --------------
    # ----------------------------------------------------------
    def continuousTableForChannel(self, ch, frq, amp):
        self.sendCmd('TABLE,APPEND,%d, %s, %s, 0x0, 0x0' % (ch,frq,amp))

    def continuousTablesForChannels(self, channels=(1,2,3,4), start=True):
        if 1 in channels and self.topticaLockWhenUpdating: self.setTopticaLock(False)
        self.prepareChannelsForTable(prepareChannels=channels,reArm = False)  # Deletes all existing tables on QuadRF. Also sets limits [hard coded!] on each channel
        if 'Continuous' != self.initialValues['Operation_Mode']:
            self.logger.debug('Continuous mode. Parameters taken from quadRFMOTController.py (hard coded)')
            # Hard coding continuous parameters here, if this function was called not due to Operation_Mode = Continuous.
            zeroAmp = '0x0'
            if 1 in channels: self.continuousTableForChannel(1, '113MHz', self.Amp_Ch1)
            if 2 in channels: self.continuousTableForChannel(2, '110MHz', self.Amp_Ch2)
            if 3 in channels: self.continuousTableForChannel(3, '110MHz', self.Amp_Ch3)
            if 4 in channels: self.continuousTableForChannel(4, '133.325MHz', self.Amp_Ch4)  # Depump
        else:
            self.logger.debug('Continuous mode. Parameters taken from Config_Table.py. See Operation_Modes.')
            # Otherwise, get the data from Config_Table
            if 1 in channels: self.continuousTableForChannel(1, self.initialValues['CH1_freq'], self.initialValues['CH1_amp'])
            if 2 in channels: self.continuousTableForChannel(2, self.initialValues['CH2_freq'], self.initialValues['CH2_amp'])
            if 3 in channels: self.continuousTableForChannel(3, self.initialValues['CH3_freq'], self.initialValues['CH3_amp'])
            if 4 in channels: self.continuousTableForChannel(4, self.initialValues['CH4_freq'], self.initialValues['CH4_amp'])

        channelsStr = ''
        for ch in channels:
            channelsStr += ',%d'%ch
        if start:
            self.sendCmd('TABLE,START ' + channelsStr)

        if 1 in channels and self.topticaLockWhenUpdating: self.setTopticaLock(True)

    def turnChannelsOff(self, channels=(1, 2, 3, 4), holdLocking=True):
        if 1 in channels and holdLocking: self.setTopticaLock(False)
        self.prepareChannelsForTable(prepareChannels=channels,reArm=False)  # Deletes all existing tables on QuadRF. Also sets limits [hard coded!] on each channel
        self.logger.debug('Turning all channels off (far detuning them)')
        zeroAmp = '0x0'
        if 1 in channels: self.continuousTableForChannel(1, '113MHz', self.Amp_Ch1)
        if 2 in channels: self.continuousTableForChannel(2, '150MHz', self.Amp_Ch2)
        if 3 in channels: self.continuousTableForChannel(3, '150MHz', self.Amp_Ch3)
        if 4 in channels: self.continuousTableForChannel(4, '150MHz', self.Amp_Ch4)

        channelsStr = ''
        for ch in channels:
            channelsStr += ',%d' % ch
        self.sendCmd('TABLE,START ' + channelsStr)

        if 1 in channels and holdLocking:
            self.setTopticaLock(True)
