from UtilityResources.DPO7254DataAcquisition import DPO7254Visa
import matplotlib.pyplot as plt

from quadRFController import QuadRFController
import numpy as np
scope = DPO7254Visa(ip='132.77.54.241')

scope.acquireData(chns=[3,4])
d = scope.wvfm[4]

scope.plotData(ch=3, show=False)
scope.plotData(ch=4)

def appendLinesToQuad(quad, ch, freqs, amps, durations):
    if len(freqs) == len(amps) and len(freqs) == len(durations):
        for i, f in enumerate(freqs):
            l = 'TABLE,APPEND,{ch}, {freq}Hz, {amp}dbm, 0x0, {duration}ms'.format(ch=ch, amp=amps[i], freq=freqs[i], duration=durations[i])
            quad.sendCmd(l)
    else:
        print('QuadRF lines length mismatch!')

q = QuadRFController(debugging=True)  # updates values on QuadRF (uploads table)
q.prepareChannelForTable(1, reArm = True, trigger = 'RISING', limit = '32dbm')

durations = [1000] + [0.8] * 50
freqs = [113] + list(np.linspace(113,90,len(durations) - 1))

amps = [29] + [29] * (len(durations) - 1)

appendLinesToQuad(quad = q, ch = 1, freqs = freqs, amps = amps, durations = durations)
q.armChannel(ch=1)

# plt.plot(scope.timeData, d)
# plt.show()
#
# i06=[]
# for i,f in enumerate(msmntAOMFreq):
#     pForF = normalizedPowerRes[i]
#     idx = (np.abs(pForF - 0.5)).argmin()
#     i06.append(msmntAOMPower[idx])
#
# plt.plot(msmntAOMFreq, i06)