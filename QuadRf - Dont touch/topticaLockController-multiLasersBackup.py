import telnetlib
import time

locking_host = '132.77.55.156'
top1_port = 60001
top2_port = 60002


class TopticasLockController:
    # Topticas should contain the numbering of the Topticas, i.e. [1,2]; or [1] if only 1 is needed
    def __init__(self, topticas):
        self.lasersIDs = topticas
        self.telnets = [None, telnetlib.Telnet(locking_host, top1_port)]
        # first is None, so as to keep numbering of Toptica 1 & 2 (and avoid "Toptica 0")
        #self.telnets = [None, telnetlib.Telnet(locking_host, top1_port), telnetlib.Telnet(locking_host, top2_port)]
        # self.top1_tn = telnetlib.Telnet(locking_host, top1_port)
        # self.top2_tn = telnetlib.Telnet(locking_host, top2_port)

    def readMessages(self, laserID):
        if laserID not in self.lasersIDs:
            print('Error: wrong laserID; could not read message')
            return
        print('Messages from laserID = ', laserID)
        print(self.telnets[laserID].read_very_eager().decode('ascii'))

    def sendMessage(self, laserID, msg):
        if laserID not in self.lasersIDs:
            print('Error: wrong laserID; could not send message')
            return
        if type(msg) is not str:
            print('Error: message must be a string')
            return
        self.telnets[laserID].write(msg.encode('ascii') + b"\n")
        #time.sleep(3)
        #self.readMessages(laserID)

tc = TopticasLockController([1])
time.sleep(1)
tc.sendMessage(1,'autolock:lock:enable?')
time.sleep(0.15)
tc.readMessages(1)
#time.sleep(2)
#tc.readMessages(1)