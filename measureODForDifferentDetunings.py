from UtilityResources.PM100DataAxquisition import PM100DVisa
from UtilityResources.DPO7254DataAcquisition import DPO7254Visa

import Config_with_SNSPDs_and_QuadRF as Config
from CoolingOptimization import CoolingSequenceOptimizer
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
datadir = os.path.join("U:\\", "Lab_2021-2022", "Experiment_results", "OD")

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
    datadir = os.path.join("U:\\", "Lab_2021-2022", "Experiment_results", "OD")
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


def measureOD(instr, ch = 1, indices = []):
    msmnt = []
    if type(instr) is DPO7254Visa:
        instr.acquireData(chns=[ch])
        msmnt = instr.wvfm[ch]
        darkCountMean = np.mean(msmnt[:indices[0]])
        firstPulseMean = np.mean(msmnt[indices[0]:indices[1]])
        lastPulseMean = np.mean(msmnt[indices[2]:indices[3]])
        OD = np.log(np.abs(lastPulseMean - darkCountMean) / np.abs(firstPulseMean - darkCountMean))
        return OD, msmnt
    else:
        print('Measurment instrument is not well defined!')
        return None
# ---------------- set up measurement ------------------
ODFreqs = np.linspace(113 - 2,113 + 2,30) * 1e6  # HZ
sleepTime = 10 # sleep after changing AOM power and msrmnt (in seconds)

initialValues = Initial_Values
experiment = CoolingSequenceOptimizer(Config.config)

Scope_wvfm_res = []
OD_Res_volts = []
indices = [530,910,2120,2518] # input manually = indices where the od pulses appear on scope:
#scope.acquireData(chns=[2])
#msmnt = scope.wvfm[2]
#plt.plot(msmnt)
#np.where(msmnt > 0.03) # - input here some threshold and find indices where the msmnt complies with it
# ----------------- Measure power ----------------------
input('Hit ENTER for calbiration start...')
scope = DPO7254Visa(ip='132.77.54.241')
for freq in ODFreqs:
    experiment.updateValue('Pulse_1_CH1_Freq_f', freq, True)
    print(f'Now taking scope wavefrom for freq: {freq}HZ')
    sleep(sleepTime)
    ODRes, wvfm = measureOD(scope, ch=2,indices = indices)
    OD_Res_volts.append(ODRes)
    Scope_wvfm_res.append(wvfm)

#
print('Done.')
np.savez(savePath, RFFreqs =ODFreqs ,Scope_wvfm_res=Scope_wvfm_res, OD_Res_volts=OD_Res_volts)
deltaFreqsInMHZ = (np.array(ODFreqs) - 113e6)*1e-6
deltaFreqsInMHZForFit = np.linspace(deltaFreqsInMHZ[0], deltaFreqsInMHZ[-1], 100)

# ---- fit lorentzian ----
from scipy.optimize import curve_fit
def lorentzian( x, x0, a, gam ):
    return a * gam**2 / (gam**2 + ( x - x0 )**2)
popt, pcov = curve_fit(lorentzian, deltaFreqsInMHZ, OD_Res_volts)
x0, a, gam = popt
# ---- plot results ----
plt.figure()
plt.plot(deltaFreqsInMHZ, OD_Res_volts, 'x')
plt.plot(deltaFreqsInMHZForFit, lorentzian(deltaFreqsInMHZForFit, *popt), 'g--')
plt.title(f'OD at different pulse freqs. \n Lorentzian fit: X_0 = {x0}, Amp = {a}, gamma = {gam}')
plt.ylabel('OD')
plt.xlabel('Detuning (of 113MHz) Frequency [MHz]')
plt.show()
# # # ----------------- Plot results ----------------------
#
# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
#
# # Plot data.
# X = msmntAOMFreq
# Y = msmntAOMPowerMultiplier
# X, Y = np.meshgrid(Y, X) # note- reverse (shape...)
# Z = np.array(powerRes) / np.max(powerRes)
# normalizedPowerRes = Z
# # Plot the surface.
# ax = plt.axes(projection='3d')
# #ax.contour3D(Y,X, Z, 50, cmap='binary')
# ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
# ax.set_xlabel('RF Power [dbm]')
# ax.set_ylabel('Frequency [MHz]')
# ax.set_zlabel('Optical Power [arb.]')
# plt.show()

