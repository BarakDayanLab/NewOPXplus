import telnetlib
import time

#locking_host = '132.77.55.156'
#top1_port = 60001

class TopticaLockController:

    def __init__(self, locking_host = '132.77.55.156', top1_port = 60001):
        self.telnet = telnetlib.Telnet(locking_host, top1_port)
        self.timeout = 1
        # Return a tuple of three items: the index in the list of the first regular expression that matches; the match object returned; and the text read up till and including the match.
        r = self.telnet.expect(['^Welcome to DigiLock110'.encode('ascii')], timeout=self.timeout)
        print('Connected to DigiLock110 successfully [%s].' % (locking_host + ':' + str(top1_port)))

    def readMessages(self,):
        print(self.telnet.read_very_eager().decode('ascii'))

    def sendMessage(self, msg):
        if type(msg) is not str:
            print('Error: message must be a string')
            return
        self.telnet.write(msg.encode('ascii') + b"\n")
        # After sending message, wait @timeout until response.
        # Response contains the message sent + response
        # Return a tuple of three items: the index in the list of the first regular expression that matches; the match object returned; and the text read up till and including the match.
        response = self.telnet.expect([('\^' + msg).encode('ascii')],timeout=self.timeout)
        #print(response[2].decode('ascii').replace(msg, '').strip())
        print(response[2].decode('ascii').strip())

    def setLocking(self, on):
        #msg = 'autolock:lock:enable='
        msg = 'autolock:lock:hold='
        flag = 'true' if on else 'false'
        flag = 'false' if on else 'true'
        self.sendMessage(msg + flag)
    def isLocking(self):
        self.sendMessage('autolock:lock:enable?')

