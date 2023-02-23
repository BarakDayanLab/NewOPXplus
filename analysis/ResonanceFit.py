import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import least_squares
import matplotlib
from scipy.signal import find_peaks
from UtilityResources.DPO7254DataAcquisition import DPO7254Visa,printGreen
import time

class ResonanceFit():
    def __init__(self):
        self.scope = DPO7254Visa(ip='132.77.54.241')
        printGreen(s='please make zero transmission on the scope, to take T=0 trace, and then press Enter')
        input()

        # parameters
        print('taking trace for zero transmission')
        self.trans0,self.refl0,_ = self.get_trace_from_scope()
        fig_T,ax_T = plt.subplots(figsize=(8, 6))
        fig_R,ax_R = plt.subplots(figsize=(8, 6))
        while True:
            time.sleep(30)
            print('taking traces')
            trans, refl, rubidiumLines = self.get_trace_from_scope()


            self.start_index = 0
            self.end_index = -1


            trans100 = np.mean(trans[1:10000]) * np.ones(len(trans))
            self.delta_T = trans - self.trans0
            self.T_abs = self.delta_T / (trans100 - self.trans0)
            self.R_abs = (refl - self.refl0)

            # calibrate range from rubidium lines
            self.t = self.calibrate_x_range(rubidiumLines)

            deep_min, peak_R = self.find_extremum()

            # fit transmission
            print('fit transmission')
            x0_T = [15e6, 10e6, 7e6, deep_min - 50e6, 0]
            self.fit_T,self.kappa_i_T, self.h_fit_T = \
                self.fit_to_resonance(x0=x0_T,bounds=(-np.inf, np.inf),Y=self.T_abs,fit_func=self.func_Lorenzian_T)

            # fit reflection
            print('fit reflection')

            x0_R = [8e6, 20e6, 7e6, peak_R, 0, 2]
            self.fit_R,self.kappa_i_R, self.h_fit_R =\
                self.fit_to_resonance(x0=x0_R,bounds=([0, 0, 0, -np.inf, -np.inf, 0], np.inf),
                                                     Y=self.R_abs,fit_func=self.func_Lorenzian_R)

            # plot
            ax_T.clear()
            ax_R.clear()
            self.plot_fit(fit=self.fit_T,Y=self.T_abs,func=self.func_Lorenzian_T,
                          label='Transmission',ylim=[0, 1.1],res_params=[self.kappa_i_T, self.h_fit_T],fig=ax_T)
            self.plot_fit(fit=self.fit_R,Y=self.R_abs,func=self.func_Lorenzian_R,label='Reflection',
                          ylim=None,res_params=[self.kappa_i_R, self.h_fit_R],fig=ax_R)

    def find_extremum(self):
        # find min and max in resonance
        deep_min = float(min(self.t[np.where(self.T_abs[self.start_index:self.end_index] == self.T_abs[self.start_index:self.end_index].min())]))
        peak_R = float(max(self.t[np.where(self.R_abs[self.start_index:self.end_index] == self.R_abs[self.start_index:self.end_index].max())]))
        return deep_min, peak_R
    def get_trace_from_scope(self):
        self.scope.acquireData(chns=[1, 2, 3, 4])
        trans = np.array(self.scope.wvfm[1])
        refl = np.array(self.scope.wvfm[2])
        rubidiumLines = np.array(self.scope.wvfm[4])
        return trans,refl,rubidiumLines

    def calibrate_x_range(self,rubidiumLines):
        peaks, prop = find_peaks(rubidiumLines, prominence=0.022, distance=1000)  # width=50, rel_height=0.5)
        plt.figure(1)
        plt.clf()
        plt.plot(rubidiumLines)
        plt.plot(peaks, np.array(rubidiumLines)[peaks], "x")
        indx_to_freq = (156.947e6 / 2) / (peaks[-1] - peaks[-2])
        t = np.array(list(range(len(self.trans0)))) * indx_to_freq  # Calibration
        return t

    def fit_to_resonance(self,x0,bounds,Y,fit_func):
        # init guess and fit
        fit = least_squares(fit_func, x0,bounds=bounds, loss='soft_l1', f_scale=0.1,
                                     args=(self.t[self.start_index:self.end_index],
                                           Y[self.start_index:self.end_index]))
        print('fit parameters are:', fit.x)
        kappa_i_h = round(fit.x[1] - fit.x[0], 2)
        h_fit = round(fit.x[2], 2)
        print('k_i=',kappa_i_h)
        print('h=',h_fit)
        return fit,kappa_i_h, h_fit

    def plot_fit(self,fit,Y,func,label,ylim,res_params,fig):
        matplotlib.rcParams.update({'font.size': 22})
        fig.plot((self.t[self.start_index:self.end_index] - fit.x[3]) / 1e6,
                 Y[self.start_index:self.end_index], 'b-', label='data')
        fig.plot((self.t[self.start_index:self.end_index] - fit.x[3]) / 1e6,
                 func(fit.x, self.t[self.start_index:self.end_index], 0), 'r--',
                 label='fit')
        fig.set_ylabel(label, fontsize=28)
        fig.set_xlabel('$\Delta [MHz]$', fontsize=28)
        fig.set_xlim([(self.t[self.start_index] - fit.x[3]) / 1e6, (self.t[self.end_index] - fit.x[3]) / 1e6])
        fig.set_ylim(ylim)
        fig.text(220, 0.05,
                 "$\kappa_i$ = {k:.2f} MHz".format(k=res_params[0] / 1e6) + '\n h = {h:.2f} MHz'.format(h=res_params[1] / 1e6),
                 size=12, rotation=0.1,
                 ha="left", va="bottom",
                 bbox=dict(boxstyle="round",
                           ec=(1., 0, 0),
                           fc=(1, 1, 1),
                           )
                 )
        fig.legend(loc='lower left')
        plt.pause(0.5)


    def func_Lorenzian_T(self,x, t, y):  # with h
        z = x[4] + np.power(np.abs(1 + 2 * 1j * x[0] * (t - x[3] - 1j * x[1]) / \
                                   (np.power((t - x[3] - 1j * x[1]), 2) - np.power(x[2], 2))), 2) - y
        return z

    def func_Lorenzian_R(self,x, t, y):
        z = x[4] + x[5] * np.power(
            np.abs(2 * x[0] * x[2] / (np.power((1j * (t - x[3]) + x[1]), 2) + np.power(x[2], 2))), 2) - y
        return z

if __name__ == '__main__':
    o = ResonanceFit()


#
# Trans0 = Data0['Ch_1_Data']
# Refl0 = Data0['Ch_2_Data']
#
#
# Trans1 = np.mean(Trans[1:10000])*np.ones(len(Trans))
#
# # plt.figure()
# # plt.plot(Trans)
# # plt.plot(Trans0)
# # plt.plot(Trans1)
# # plt.plot(Refl)
# # plt.plot(Refl0)
# # plt.plot(RubidiumLines)
#
#
# peaks, prop = find_peaks(RubidiumLines, prominence=0.022, distance=1000) # width=50, rel_height=0.5)
# plt.figure(1)
# plt.plot(RubidiumLines)
# plt.plot(peaks, np.array(RubidiumLines)[peaks], "x")
# indx_to_freq = (156.947e6 / 2) / (peaks[-1] - peaks[-2])
# t = np.array(list(range(len(Trans)))) * indx_to_freq  # Calibration
# # t = np.array(list(range(len(Trans)))) * (156.947e6/2)/(6.178e4 - 5.862e4)  # Manual Calibration
#
#
# # Gabi's fitting:
#
# T_1 = np.array(Trans1)
# T = np.array(Trans)
# R = np.array(Refl)
# T_0 = np.array(Trans0)
# R_0 = np.array(Refl0)
# delta_T = T - T_0
#
# T_abs = (T - T_0) / (T_1 - T_0)
# R_abs = (R - R_0)
#
# # T_abs = (delta_T)/(np.average(delta_T[3400:3500]))
#
# # %% Choose window manually
#
# start_index = 0
# end_index = -1
#
# # peaks = find_peaks_cwt(1/T_abs[start_index:end_index], np.arange(50, 200))
# # peaks = find_peaks_cwt(-1*T_abs[start_index:end_index], np.arange(50, 200))
# # peaks = np.array(peaks) + start_index
# # print(peaks)
#
# deep_min = float(min(t[np.where(T_abs[start_index:end_index] == T_abs[start_index:end_index].min())]))
# peak_R = float(max(t[np.where(R_abs[start_index:end_index] == R_abs[start_index:end_index].max())]))
#
#
# # %% Plot
#
# # plt.figure(1)
# # plt.plot(t[start_index:end_index], T_abs[start_index:end_index], 'b-', label='data')
# # axes = plt.gca()
# # axes.set_xlim([t[start_index],t[end_index]])
# # axes.set_ylim([0, 1.1])
#
#
# # #%% Fit
# # def func(x, kappa_ex, kappa_t, h, delta, t):
# #     z = t + np.power(np.abs(1+1j*kappa_ex*(x - delta - 1j*kappa_t)/\
# #                            ((x - delta - 1j*kappa_t)**2 - h**2)),2)
# #     return z
# #
# # popt, pcov = curve_fit(func, t[start_index:end_index], T_abs[start_index:end_index],
# #                        bounds=([0., 0., 0., -500, -1], [10., 20., 30, 500, 1]))
# # # plt.plot(t[start_index:end_index], func(t[start_index:end_index], *popt), 'r--', label='fit')
# # # print(popt)
#
# # %% Fit 2
#
# # def func_2(x, t, y):                                                                            # with h
# #     z = x[4] + np.power(np.abs(1+2*1j*x[0]*(t - x[3] - 1j*x[1])/\
# #                            (np.power((t - x[3] - 1j*x[1]),2) - np.power(x[2],2))),2) - y
# #     return z
#
# def func_2(x, t, y):  # with h
#     z = x[4] + np.power(np.abs(1 + 2 * 1j * x[0] * (t - x[3] - 1j * x[1]) / \
#                                (np.power((t - x[3] - 1j * x[1]), 2) - np.power(x[2], 2))), 2) - y
#     return z
#
#
# # def func_3(x, t, y):                                                                            # without h
# #     z = x[3] + np.power(np.abs(1+1j*x[0]/\
# #                            ((t - x[2] - 1j*x[1]) ) ),2) - y
# #     return z
#
#
# def func_reflection(x, t, y):
#     z = x[4] + x[5] * np.power(np.abs(2 * x[0] * x[2] / (np.power((1j * (t - x[3]) + x[1]), 2) + np.power(x[2], 2))), 2) - y
#     return z
#
#
# x0 = [15e6, 10e6, 7e6, deep_min - 50e6, 0]
# res_robust_w = least_squares(func_2, x0, loss='soft_l1', f_scale=0.1,
#                              args=(t[start_index:end_index], T_abs[start_index:end_index]))
# x01 = res_robust_w.x
# # x01[2] = 100
# print(res_robust_w.x)
#
# # x0 = [0.3, 10, deep_min, 1]
# # res_robust_w_1 = least_squares(func_3, x0, loss='soft_l1', f_scale=0.1,
# #                            args=(t[start_index:end_index], T_abs[start_index:end_index]))
#
# # print(res_robust_w_1.x)
# # kappa_i_no_h = round(res_robust_w_1.x[1] - res_robust_w_1.x[0], 2)
# kappa_i_h = round(res_robust_w.x[1] - res_robust_w.x[0], 2)
# h_fit = round(res_robust_w.x[2], 2)
# print(kappa_i_h)
# print(h_fit)
#
# # Plot results
#
# matplotlib.rcParams.update({'font.size': 22})
# plt.figure(2, figsize=(8, 6))
# plt.plot((t[start_index:end_index] - res_robust_w.x[3]) / 1e6, T_abs[start_index:end_index], 'b-', label='data')
# plt.plot((t[start_index:end_index] - res_robust_w.x[3]) / 1e6, func_2(x01, t[start_index:end_index], 0), 'r--',
#          label='fit')
# plt.ylabel('T', fontsize=28)
# plt.xlabel('$\Delta [MHz]$', fontsize=28)
# axes = plt.gca()
# axes.set_xlim([(t[start_index] - res_robust_w.x[3]) / 1e6, (t[end_index] - res_robust_w.x[3]) / 1e6])
# axes.set_ylim([0, 1.1])
# plt.text(220, 0.05, "$\kappa_i$ = {k:.2f} MHz".format(k=kappa_i_h / 1e6) + '\n h = {h:.2f} MHz'.format(h=h_fit / 1e6),
#          size=12, rotation=0.1,
#          ha="left", va="bottom",
#          bbox=dict(boxstyle="round",
#                    ec=(1., 0, 0),
#                    fc=(1, 1, 1),
#                    )
#          )
# plt.legend(loc='lower left')
# # plt.show()
#
# # matplotlib.rcParams.update({'font.size': 22})
# # plt.figure(3, figsize=(8, 6))
# # plt.plot(t[start_index:end_index]-res_robust_w_1.x[2], T_abs[start_index:end_index], 'b-', label='data')
# # plt.plot(t[start_index:end_index]-res_robust_w_1.x[2], func_3(res_robust_w_1.x, t[start_index:end_index], 0), 'r--', label='fit')
# # plt.ylabel('T', fontsize=28)
# # plt.xlabel('$\Delta [MHz]$', fontsize=28)
# # axes = plt.gca()
# # axes.set_xlim([t[start_index]-res_robust_w_1.x[2],t[end_index]-res_robust_w_1.x[2]])
# # axes.set_ylim([0, 1.3])
# # plt.legend(loc='lower left')
# # plt.text(75, 0.05, "$\kappa_i$ = {} MHz".format(kappa_i_no_h), size=12, rotation=0.1,
# #          ha="right", va="bottom",
# #          bbox=dict(boxstyle="round",
# #                    ec=(1., 0, 0),
# #                    fc=(1, 1, 1),
# #                    )
# #          )
# # plt.legend(loc='lower left')
# # plt.show()
#
#
# # For reflection (Dor):
# x0_r = [8e6, 20e6, 7e6, peak_R, 0, 2]
# res_robust_r = least_squares(func_reflection, x0_r, bounds=([0, 0, 0, -np.inf, -np.inf, 0], np.inf), loss='soft_l1',
#                              f_scale=0.1,
#                              args=(t[start_index:end_index], R_abs[start_index:end_index]))
# x01_r = res_robust_r.x
# print(res_robust_r.x)
# kappa_i_h_r = round(res_robust_r.x[1] - res_robust_r.x[0], 2)
# h_fit_r = round(res_robust_r.x[2], 2)
# print(kappa_i_h_r)
# print(h_fit_r)
#
# # Plot results
#
# matplotlib.rcParams.update({'font.size': 22})
# plt.figure(3, figsize=(8, 6))
# plt.plot((t[start_index:end_index] - res_robust_r.x[3]) / 1e6, R_abs[start_index:end_index], 'b-', label='data')
# # plt.plot((t[start_index:end_index] - res_robust_r.x[3]) / 1e6, func_reflection([8.68217686e+06, 1.97951378e+07,
# #                                                                                 7.14531420e+06, 5.56435891e+08,
# #                                                                                 -1.78528038e-03], t[start_index:end_index], 0), 'r--', label='fit')
# plt.plot((t[start_index:end_index] - res_robust_r.x[3]) / 1e6, func_reflection(x01_r, t[start_index:end_index], 0),
#          'r--', label='fit')
# plt.ylabel('R', fontsize=28)
# plt.xlabel('$\Delta [MHz]$', fontsize=28)
# axes = plt.gca()
# axes.set_xlim([(t[start_index] - res_robust_r.x[3]) / 1e6, (t[end_index] - res_robust_r.x[3]) / 1e6])
# # axes.set_ylim([0, 1.1])
# plt.text(220, 0.05,
#          "$\kappa_i$ = {k:.2f} MHz".format(k=kappa_i_h_r / 1e6) + '\n h = {h:.2f} MHz'.format(h=h_fit_r / 1e6), size=12,
#          rotation=0.1,
#          ha="left", va="bottom",
#          bbox=dict(boxstyle="round",
#                    ec=(1., 0, 0),
#                    fc=(1, 1, 1), )
#          )
# plt.legend(loc='lower left')
# plt.show()