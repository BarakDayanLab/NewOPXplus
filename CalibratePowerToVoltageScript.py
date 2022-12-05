from UtilityResources.PM100DataAxquisition import PM100DVisa
from quadRFMOTController import QuadRFMOTController
import numpy as np
from time import sleep
from datetime import date,datetime
import os

import matplotlib.pyplot as plt
from UtilityResources.scpi import Redpitaya


def saveData(voltage,power):
    now = datetime.now()
    today = date.today()
    datadir = os.path.join("U:\\", "natang", "PM100D calibration")
    todayformated = today.strftime("%B-%d-%Y")
    todaydatadir = os.path.join(datadir, todayformated)
    nowformated = now.strftime("%Hh%Mm%Ss.csv")
    try:
        os.makedirs(todaydatadir)
        print("Created folder DATA/FPO7254/%s" % (todayformated))
        print("Data Saved")
    except FileExistsError:
        # self.saveData()
        print("Data Saved")

    # np.savez_compressed(os.path.join(todaydatadir, nowformated), Voltage=voltage, Power=power)
    np.savetxt(os.path.join(todaydatadir, nowformated), np.asarray([power, voltage]), delimiter=",")


# ---------------- set up measurement ------------------
msmntAOMPower = np.linspace(12,20,15)
powerRes = []
voltageRes = []
sleepTime = 2 # sleep after changing AOM power and msrmnt (in seconds)
nMsrmnts = 3 # avg over nMsrmnts
initialValues = {'Operation_Mode': 'Continuous',
                                    'CH1_freq': '113MHz', #'93MHz',
                                    # In order to turn a channel off, set amp to '0x0'
                                    'CH1_amp': '10.95dbm' }
QuadRFMOTController(initialValues=initialValues, updateChannels=[1], debugging=True)
# ----------------- Measure power ----------------------
input("Place power meter in front of the detector then hit Enter")

PM100D = PM100DVisa()

for power in msmntAOMPower:
    initialValues['CH1_amp'] = '%.3fdbm' % power
    QuadRFMOTController(initialValues=initialValues, updateChannels=[1], debugging=False)
    power = []
    for a in range(nMsrmnts):
        sleep(sleepTime)
        power.append(PM100D.getPower())
    powerRes.append(np.mean(power))
    print(power)

# ----------------- Measure Voltage ----------------------
# input("Move power meter from detector then hit Enter")
#
#
# # rp = Redpitaya(host = '132.77.55.19')
# rp = Redpitaya(host = 'rp-f08c22.local')
# rp.set_triggerSource('NOW')
#
# for power in msmntAOMPower:
#     initialValues['CH3_amp'] = '%.3fdbm' % power
#     QuadRFMOTController(initialValues=initialValues, updateChannels=[3], debugging=False)
#     v = []
#     for a in range(nMsrmnts):
#         sleep(sleepTime)
#         v.append(np.mean(rp.get_trace(1)))
#     voltageRes.append(np.mean(v))
#     print(power)

# ----------------- Plot results ----------------------
print('Done.')
saveData(voltage = voltageRes, power = powerRes)
#fit_coef = np.polyfit(x = powerRes, y = voltageRes, deg = 1)
#fit = np.poly1d(fit_coef)

# plt.plot(voltageRes,powerRes,'x')
# plt.plot(fit(powerRes),'-')

plt.plot(msmntAOMPower, powerRes,'x')
fit_coef = np.polyfit(x = msmntAOMPower, y = powerRes, deg = 2)
fit = np.poly1d(fit_coef)
plt.plot(msmntAOMPower,fit(msmntAOMPower),'-')
plt.show()
