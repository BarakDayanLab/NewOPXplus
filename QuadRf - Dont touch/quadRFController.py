# Wrtitten by ngordon1

# pyserial & mogdevice modules for python should be installed.
# pyserial can be installed vias pip install
# You can download the file from: https://www.moglabs.com/support/software/connection. Look for python driver; then put it in python-Lib library.

## --------------- Usage example: ----------------------------

# mogRFDevice = QuadRFController('COM3')
# mogRFDevice.armMogRFTable('mogScriptExample.txt')

## --------------- Formatting of commands file: ---------------

# Code ignores empty lines & lines of which 1st char is '#' (could be used for comments)
# Each line should be compirsed of 5 data-words:
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
#from mogdevice import MOGDevice
import os.path
import numpy as np


class QuadRFController:

    def __init__(self, devPort):
        # connect to the device
        try:
            self.dev = MOGDevice(devPort)
        except:
            print('Could not connect to device-port %s. Try different port. Usual values are \'COM3\'')
            return()  
        # print some information
        print('Connected to QuadRF. \n Device info:', self.dev.ask('info'))

    def __del__(self):
        self.dev.close()
        
    def lineIsTableCommand(self, l):
        if  type(l) is str and len(l.strip()) > 0 and l[0] != '#':
            return(True)
        else:
            return(False)

    # @txtfile: location for table textfile.
    # @reArm: activate rearming option for all chans (see QuadRF manual)
    # @checkCommandsFile: check @txtfile before loading; If it isn't properly constructed, no change will occour
    # @trigger: either 'RISING' or 'FALLING' (see quadRF manual)   
    def armMogRFTable(self, txtfile = None, reArm = True, checkCommandsFile = True, trigger = 'RISING'):
        if not os.path.isfile(txtfile):
            print('Could not find %s; \n check for file location.')
        
        lines = []
        with open(txtfile) as f:
            lines = f.readlines()
            chns = [True] *4

            # Check commands file
            if checkCommandsFile:
                chns = [False] *4
                for i, l in enumerate(lines):
                    if self.lineIsTableCommand(l):
                        c = int(l.split(',')[0])
                        try:
                            chns[c - 1] = True
                        except:
                            print('Line %d is corrupted' % (i+1))
                            print('Table was not loaded.')
                            return

            # Load general commands
            for i,c in enumerate(chns):
                if c is True:
                    self.prepareChannelForTable(i + 1,reArm, trigger)


            # Load table command      
            for l in lines:
                if self.lineIsTableCommand(l):
                    self.dev.cmd('TABLE,APPEND,%s' % l)
            # ARM 
            for i,c in enumerate(chns):
                if c is True:
                    self.armChannel(i+1)
                    

    

    def prepareChannelForTable(self, ch, reArm = True, trigger = 'RISING', limit = '4dbm'):
        ch = int(ch)
        self.setLimit(ch, limit) # Current limit on channel
        self.dev.cmd('MODE,%d,TSB' % ch) # Out in table mode
        self.dev.cmd('TABLE,CLEAR,%d' % ch) # clear all previous commands in table
        self.dev.cmd('TABLE,EDGE,%d,%s' % (ch , trigger.upper())) # Set trigger option
        if reArm:
           self.dev.cmd('TABLE,REARM,%d,on' % (ch))
        else:
            self.dev.cmd('TABLE,REARM,%d,off' % (ch))
        
    def armChannel(self, ch):
        self.dev.cmd('TABLE,ARM,%d' % (ch))

    # @initalAmp, @finalAmp are numbers (dBm are assumed) or strings 
    # @ch is channel int, @freq could be either number (MHz assumed) or string.
    # @duration is total duration of pulse [ms assumed]
    def appendLinearEnvelopePulse(self, ch, initialFreq, finalFreq, initialAmp, finalAmp, duration, steps = 50):
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
        durationScale = str(duration).replace(str(d),'')
        if durationScale == '':
            durationScale = 'ms'   
        
        amps = np.linspace(float(iAmp),float(fAmp),steps)
        freq = np.linspace(float(iFreq),float(fFreq),steps)
        t = float(d) / steps
        for index, a in enumerate(amps):
            #print('TABLE,APPEND,%d,%.2f%s,%.2f%s,0x0,%.2f%s' % (ch, float(f), freqScale,float(a), ampScale, float(t),durationScale))
            self.dev.cmd('TABLE,APPEND,%d,%.2f%s,%.2f%s,0x0,%.2f%s' % (ch, float(freq[index]), freqScale,float(a), ampScale, float(t),durationScale))


        
    def getNumericValueFromString(self, s):
        res = ''.join([n for n in str(s) if (n in['-','.'] or n.isdigit())])
        return(res)
        
    def setLimit(self, ch, limit):
        self.dev.cmd('LIMIT, %s, %s' % (str(ch), str(limit)))

    def amplitudeMultiplierToDBm(self, ratio):
        return 10 * np.log10(ratio)

