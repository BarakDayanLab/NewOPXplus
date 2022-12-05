import telnetlib
import time

#locking_host = '132.77.55.156'
#top1_port = 60001

class TopticaLockController:

    def __init__(self, locking_host = '132.77.55.157', top1_port = 60001):
        self.telnet = telnetlib.Telnet(locking_host, top1_port)
        self.timeout = 1
        # Return a tuple of three items: the index in the list of the first regular expression that matches; the match object returned; and the text read up till and including the match.
        r = self.telnet.expect(['^Welcome to DigiLock110'.encode('ascii')], timeout=self.timeout)
        if r[0] == -1 and str(r[2]).find('Welcome to DigiLock110 remote interface') != -1:
            print('\033[92m' + 'Connected to DigiLock110 successfully [%s].' % (locking_host + ':' + str(top1_port)) + '\033[0m') # print green
        else:
            print('\033[93m' + 'Could not connect to DigiLock via Telnet. ' + '\033[0m') # print yellow

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
        print('\033[94m' + response[2].decode('ascii').strip() + '\033[0m') # print blue

    def setLocking(self, on):
        msg = 'autolock:lock:hold='
        flag = 'false' if on else 'true'
       # print(time.ctime() + ' ' + msg + flag)
        self.sendMessage(msg + flag)
    def isLocking(self):
        self.sendMessage('autolock:lock:enable?')

