from UtilityResources.PM100DataAxquisition import PM100DVisa
from UtilityResources.DPO7254DataAcquisition import DPO7254Visa

import Config_with_SNSPDs_and_QuadRF as Config
from Archive.CoolingOptimization import CoolingSequenceOptimizer
from quadRFMOTController import QuadRFMOTController
from Config_Table import Initial_Values
import numpy as np
from time import sleep
from datetime import date,datetime
import os

import matplotlib.pyplot as plt

# ----- Create Directory ------
now = datetime.now()
today = date.today()
datadir = os.path.join("U:\\", "Lab_2021-2022", "DATA", "DPO7254")
datadir = os.path.join("U:\\", "natang", "Lock AOM calibration")
todayformated = today.strftime("%B-%d-%Y")
todaydatadir = os.path.join(datadir, todayformated)
nowformated = now.strftime("%Hh%Mm%Ss")
try:
    os.makedirs(todaydatadir)
    print("Created folder DATA/FPO7254/%s" % (todayformated))
    print("Data Saved")
except FileExistsError:
    # self.saveData()
    print("Data Saved")

savePath = os.path.join(todaydatadir, nowformated)

def saveData(voltage,power):
    now = datetime.now()
    today = date.today()
    datadir = os.path.join("U:\\", "natang", "Lock AOM calibration")
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

def measurePower(instr, ch = 1):
    if type(instr) is DPO7254Visa:
        instr.acquireData(chns=[ch])
        msmnt = instr.wvfm[ch]
        start = int(len(msmnt)/2)
        return np.mean(msmnt[0:start])
    elif type(instr) is PM100DVisa:
        return instr.getPower()
    else:
        print('Measurment instrument is not well defined!')
        return None
# ---------------- set up measurement ------------------
msmntAOMFreq = np.linspace(93,135,21) * 1e6  # HZ
msmntAOMPowerMultiplier = np.linspace(1,0.2,num = 20) # multiplier
# note for debugging, it's easier if the sizes of the measured arrays are different.

# We measure the beam power at different frequencies, at different (RF) powers.
# Using this info we can calibrate for the transimition of the AOM

sleepTime = 2.5 # sleep after changing AOM power and msrmnt (in seconds)
nMsrmnts = 1 # avg over nMsrmnts


initialValues = Initial_Values
#initialValues['Operation_Mode'] = 'PrePGC_Fountain' # [ms]
#initialValues['Pulse_1_duration'] = 40 # [ms]

experiment = CoolingSequenceOptimizer(Config.config)
experiment.updateValue('Pulse_1_duration', 40, True)
#QuadRFMOTController(initialValues=initialValues, updateChannels=[1], debugging=True)
# ----------------- Measure power ----------------------
input('Hit ENTER for calbiration start...')
# PM100D = PM100DVisa(port = 'USB0::0x1313::0x8078::P0028892::INSTR')
scope = DPO7254Visa(ip='132.77.54.241')
powerRes = []
rfFreqs, rfPowers_dbm, rfPower_multipliers = [],[],[]
for freq in msmntAOMFreq:
    powerForFreqRes = []
    for powerMult in msmntAOMPowerMultiplier:
        rfPower_dbm = 15.25 + QuadRFMOTController.amplitudeMultiplierToDBm(None, powerMult)
        experiment.updateValue('Pulse_1_amp_f', powerMult, False)
        experiment.updateValue('Pulse_1_CH1_Freq_f', freq, True)
        print(f'Now measuring freq: {freq}HZ; Power: {[powerMult]}')
        power = [] # (for averaging. we don't keep this)
        for a in range(nMsrmnts):
            sleep(sleepTime)
            power.append(measurePower(scope, ch= 1))
        print(power)
        powerForFreqRes.append(np.mean(power))
        rfFreqs.append(freq)
        rfPowers_dbm.append(rfPower_dbm)
        rfPower_multipliers.append(powerMult)
    powerRes.append(powerForFreqRes)


#
print('Done.')
np.savez(savePath, RFFreqs =msmntAOMFreq ,RFPowers_muli=msmntAOMPowerMultiplier, RFPowers_dbm=rfPowers_dbm,OpticalPowers = powerRes)
# file = np.load(r'U:\natang\Lock AOM calibration\November-07-2022\09h23m18s.npz')
# msmntAOMFreq, msmntAOMPowerMultiplier, rfPowers_dbm, powerRes = file['RFFreqs'], file['RFPowers_muli'],file['RFPowers_dbm'],file['OpticalPowers']
# rfPowers_dbm = rfPowers_dbm[:len(msmntAOMPowerMultiplier)]
# # ----------------- Plot results ----------------------

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Plot data.
X = msmntAOMFreq
Y = msmntAOMPowerMultiplier
X, Y = np.meshgrid(Y, X) # note- reverse (shape...)
Z = np.array(powerRes) / np.max(powerRes)
normalizedPowerRes = Z
# Plot the surface.
ax = plt.axes(projection='3d')
#ax.contour3D(Y,X, Z, 50, cmap='binary')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
ax.set_xlabel('RF Power [dbm]')
ax.set_ylabel('Frequency [MHz]')
ax.set_zlabel('Optical Power [arb.]')
plt.show()


# # ----------------- Create calibration function ----------------------

def dataToCalibrationFunction(msmntAOMFreq, msmntAOMPowerMultiplier,powerRes):
    # We would like to reverse data, that is, to get a function:
    # F(aomFreq, normalizedPower) -> aomRFPower

    msmntAOMFreq = msmntAOMFreq * 1e6 # convert MHz to Hz

    normalizedPowerRes = np.array(powerRes) / np.max(powerRes) # Normalize Z
    #x = np.array(list(msmntAOMFreq) * len(msmntAOMPower)) # i.e., [1,2,3,1,2,3...,1,2,3]
    x = np.array([[f] * len(msmntAOMPowerMultiplier) for f in msmntAOMFreq]).flatten() # i.e., [1,1,1,2,2,2...,n,n,n]
    y = normalizedPowerRes.flatten()
    z = np.array(list(msmntAOMPowerMultiplier) * len(msmntAOMFreq)) # i.e., [1,2,3,1,2,3...,1,2,3]
    from scipy import interpolate
    g = interpolate.interp2d(x, y,  z, kind='linear')
    return g

calibration_func = dataToCalibrationFunction(msmntAOMFreq, msmntAOMPowerMultiplier,powerRes)
np.savez(savePath + 'calibrationFunction', data = {'calibration_func':calibration_func, 'calibrate_RF_Power':'32dbm','calibrate_Optical_Power':0.6})


print('Looking at a point...')
# Direct function option 1
h = interpolate.interp2d(msmntAOMFreq, msmntAOMPowerMultiplier,  normalizedPowerRes.flatten(), kind='linear')
h(134,25)
# Direct function option 2
dataPoints = []
for f in msmntAOMFreq:
    for rf in msmntAOMPowerMultiplier:
        dataPoints.append([f, rf])

from scipy.interpolate import griddata
grid_z0 = griddata(dataPoints , normalizedPowerRes.flatten(),([134],[25]), method='nearest')
grid_z0

# Reverse function option 1
g = interpolate.interp2d([[f0] * len(msmntAOMPowerMultiplier) for f0 in msmntAOMFreq], normalizedPowerRes.flatten(), list(msmntAOMPowerMultiplier) * len(msmntAOMFreq),  kind='linear')
# THIS IS BAD!

# Reverse function option 2 -
# Works. cubic interpolation!
dataPoints = []
for i, f in enumerate(msmntAOMFreq):
    for op in normalizedPowerRes[i]:
        dataPoints.append([f, op])
from scipy.interpolate import griddata
grid_z1 = griddata(dataPoints , list(rfPowers_dbm) * len(msmntAOMFreq),(130e6,0.2), method='cubic')
# The follwoing is the function defenition
def foo(f,op):
    return float(griddata(dataPoints , list(rfPowers_dbm) * len(msmntAOMFreq),([f],[op]), method='cubic'))


foo = lambda f,op : float(griddata(dataPoints , list(rfPowers_dbm) * len(msmntAOMFreq),([f],[op]), method='cubic'))
foo1 =  lambda f,op : float(griddata(dataPoints , list(rfPowers_dbm) * len(msmntAOMFreq),([f],[op]), method='linear'))
foo_x = np.linspace(85e6, 140e6, 29)
foo_y = [foo(f, 0.4) for f in foo_x]
fit = np.polyfit(foo_x[:-2], foo_y[:-2],2)
fit_pol = np.poly1d(fit)
plt.figure()
plt.plot(foo_x, foo_y)
np.savez(savePath + 'calibrationFunction_griddata_cubic.npz', data = {'calibration_func':'griddata','calibration_data': [dataPoints , list(msmntAOMPowerMultiplier) * len(msmntAOMFreq)],'calibrate_RF_Power':'15.5dbm','calibrate_Optical_Power':0.45})


o = scipy.interpolate.SmoothBivariateSpline([[f0] * len(msmntAOMPowerMultiplier) for f0 in msmntAOMFreq], normalizedPowerRes.flatten(), list(msmntAOMPowerMultiplier) * len(msmntAOMFreq))

# ------- LinearNDInterpolator ------------
# Also seems to work..?
from scipy.interpolate import LinearNDInterpolator
dataPoints = []
for i, f in enumerate(msmntAOMFreq):
    for op in normalizedPowerRes[i]:
        dataPoints.append([f, op])
u = LinearNDInterpolator(dataPoints, list(msmntAOMPowerMultiplier) * len(msmntAOMFreq))
np.savez(savePath + 'calibrationFunction_u', data = {'calibration_func':u, 'calibrate_RF_Power':'32dbm','calibrate_Optical_Power':0.6})

# checkInterpFunction: returns the difference between measured and interpolated optical power
def checkInterpFunction(func, msmntAOMFreq,msmntAOMPower, normalizedPowerRes, checkForOneOpticalPower = None):
    res = []
    for i, f in enumerate(msmntAOMFreq):
        for j, p in enumerate(msmntAOMPower):
            op = normalizedPowerRes[i][j]
            if checkForOneOpticalPower is not None and checkForOneOpticalPower[0] < op < checkForOneOpticalPower[1]:
                interpP = func(f, op)
                print('F ={f}; OP = {op}, interpP = {interpP}, mark = {r}'.format(f=f,op=op, interpP=interpP, r = interpP - p))
                res.append(interpP - p)
    return(res)
r = checkInterpFunction(foo, msmntAOMFreq,rfPowers_dbm, normalizedPowerRes, checkForOneOpticalPower = [0.40,0.46])
plt.figure()
plt.plot(r)


i06=[]
for i,f in enumerate(msmntAOMFreq):
    pForF = normalizedPowerRes[i]
    idx = (np.abs(pForF - 0.5)).argmin()
    i06.append(msmntAOMPowerMultiplier[idx])

plt.plot(msmntAOMFreq, i06)
plt.show()