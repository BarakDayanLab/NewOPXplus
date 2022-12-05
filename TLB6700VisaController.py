import pyvisa

class TLB6700Controller:
    def __init__(self):

        # Start session
        self.inst = None # Place holder
        self.startSession()

        # Update view...
        self.updateView()
        self.initViewCommands()

    def initViewCommands(self):
        self.view.Chan1VoltageSpinbox['command'] = lambda: self.setVoltage(chan = 1)
        self.view.Chan2VoltageSpinbox['command'] = lambda: self.setVoltage(chan = 2)
        self.view.Chan3VoltageSpinbox['command'] = lambda: self.setVoltage(chan = 3)
        self.view.Chan4VoltageSpinbox['command'] = lambda: self.setVoltage(chan = 4)



    def startSession(self):
        print('Starting session...')
        self.inst = pyvisa.ResourceManager().open_resource(self._instName)
        print(self.inst.query('*IDN?'), ' session started successfully')

    def setVoltage(self, chan, v = None):
        if v is None:
            v = self.getVoltage(chan = chan)
        self.setChannel(chan)
        self.inst.write('VOLT %f' % v)

    def getVoltage(self, chan):
        chanIndex = chan - 1 # index starts from zero, we use interface chan (1-4)
        return(float(self.view.channelsVoltagesSpinboxes[chanIndex].get()))

    def setCurrentLimit(self, chan, v = 10.0):
        self.setChannel(chan)
        self.inst.write('CURR %f' % v)

    def getVoltageReading(self, chan): # From instrument
        self.setChannel(chan)
        return float(self.inst.query('VOLT?'))

    def getCurrent(self, chan):
        self.setChannel(chan)
        return float(self.inst.query('CURR?'))

    def getCurrentLimit(self, chan):
        self.setChannel(chan)
        return self.inst.query('VOLT?')

    def setRemote(self):
        self.inst.write('SYST:REM')

    def setChannel(self, chan):
        if type(chan) is not int and chan <=4 and chan >= 1:
            print('Error: Chan must be an int 1-4')
            return
        self.inst.write('INST OUT%d' % chan)

    # TODO: there should be only one talking-to-inst key, and both query and write should request it before sending mesgs via VISA.
    # (this is for multiple threads)
    def queryInst(self, query):
        return self.inst.quey(query)
    def writeInst(self, msg):
        self.inst.writ(msg)
    _instName : str = 'ASRL7::INSTR'