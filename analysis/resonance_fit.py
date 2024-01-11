from UtilityResources.DPO7254DataAcquisition import DPO7254Visa
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.interpolate import splrep, BSpline
from scipy.optimize import curve_fit
import time
from threading import Thread

class ResonanceFit:
    def __init__(self, channels_dict: dict, scope_ip='132.77.54.241', k_i=None, h=None):
        self.scope = DPO7254Visa(ip=scope_ip)
        self.channels_dict = channels_dict
        self.default_parameters = {"k_i": k_i, "h": h}

        self.x_axis = None

        self.transmission_0 = 0
        self.transmission_100 = 0
        self.current_transmission_spectrum = None
        self.current_reflection_spectrum = None

        self.rubidium_lines = None
        self.rubidium_peaks = None

        self.optimal_parameters = None

        self.wait_time = 1.15

    def calibrate_x_axis(self):
        self.rubidium_lines = self.read_scope_data("rubidium")
        self.rubidium_peaks, _ = find_peaks(self.rubidium_lines, prominence=0.4, distance=100)
        # plt.plot(np.arange(len(self.rubidium_lines)), self.rubidium_lines)
        # plt.scatter(self.rubidium_peaks, self.rubidium_lines[self.rubidium_peaks])
        # plt.show()
        if len(self.rubidium_peaks) != 4:
            print(ValueError("Could not find 4 peaks in the rubidium spectrum"))
            time.sleep(self.wait_time)
            return self.calibrate_x_axis()
        idx_to_freq = (156.947e6 / 2) / (self.rubidium_peaks[-1] - self.rubidium_peaks[-2])
        x_axis = np.arange(len(self.rubidium_lines)) * idx_to_freq  # Calibration

        self.x_axis = x_axis

    def read_scope_data(self, channel):
        channel_number = self.channels_dict[channel]
        self.scope.acquireData(chns=[channel_number])
        data = np.array(self.scope.wvfm[channel_number])
        # tck = splrep(np.arange(len(data)), data, s=len(data))
        # data = BSpline(*tck)(np.arange(len(tck)))
        data = np.convolve(data, np.ones(50), mode='same')[51:-51]
        return data

    def set_transmission_0(self):
        self.transmission_0 = self.read_scope_data("transmission")

    def set_transmission_100(self):
        self.transmission_100 = self.read_scope_data("transmission")

    @staticmethod
    def transmission_spectrum(x_detuning, k_ex, k_i, h, x_0, y_0, m, b):
        k_total = k_ex + k_i
        k_with_detuning = k_total + 1j * (x_detuning - x_0)
        return (m*x_detuning + b) * np.abs(1 - 2*k_ex*k_with_detuning/(k_with_detuning ** 2 + h ** 2)) ** 2 + y_0

    @staticmethod
    def reflection_spectrum(x_detuning, k_ex, k_i, h, x_0, y_0, m, b):
        k_total = k_ex + k_i
        k_with_detuning = k_total + 1j * (x_detuning - x_0)
        return (m*x_detuning + b) * np.abs(2 * k_ex * h / (k_with_detuning ** 2 + h ** 2)) ** 2 + y_0

    def transmission_spectrum_without_cavity_params(self, x_detuning, k_ex, x_0, y_0, m, b):
        k_i, h = self.default_parameters["k_i"], self.default_parameters["h"]
        return self.transmission_spectrum(x_detuning, k_ex, k_i, h, x_0, y_0, m, b)

    def transmission_spectrum_k_ex(self, x_detuning, k_ex):
        k_i, h = self.default_parameters["k_i"], self.default_parameters["h"]
        m, b = self.default_parameters["m"], self.default_parameters["b"]
        x_0, y_0 = self.default_parameters["x_0"], self.default_parameters["y_0"]
        return self.transmission_spectrum(x_detuning, k_ex, k_i, h, x_0, y_0, m, b)

    def stacked_spectrum(self, x_detuning, k_ex, k_i, h, x_0, y_0, m, b):
        return np.hstack((self.transmission_spectrum(x_detuning, k_ex, k_i, h, x_0, y_0, m, b),
                          self.reflection_spectrum(x_detuning, k_ex, k_i, h, x_0, y_0, m, b)))

    def read_transmission_spectrum(self):
        self.current_transmission_spectrum = self.read_scope_data("transmission")

    def read_reflection_spectrum(self):
        current_reflection = self.read_scope_data("reflection")
        reflection_norm_factor = 1 / (self.transmission_100 - self.transmission_0)
        self.current_reflection_spectrum = reflection_norm_factor * (current_reflection - self.transmission_0)

    def fit_current_spectrum(self):
        y_data = np.hstack((self.current_transmission_spectrum, self.current_reflection_spectrum))
        # noinspection PyTupleAssignmentBalance
        self.optimal_parameters, _ = curve_fit(self.stacked_spectrum, self.x_axis, y_data)

    def fit_transmission_spectrum(self):
        self.calibrate_x_axis()
        self.read_transmission_spectrum()

        # noinspection PyTupleAssignmentBalance
        optimal_parameters, _ = curve_fit(self.transmission_spectrum_k_ex,
                                          self.x_axis,
                                          self.current_transmission_spectrum,
                                          p0=[self.default_parameters["k_ex"]])
        return optimal_parameters

    def initialize_default_parameters(self):
        self.calibrate_x_axis()
        self.read_transmission_spectrum()
        num_points = len(self.x_axis)
        x_0_s = self.x_axis[num_points//2]

        # noinspection PyTupleAssignmentBalance
        optimal_parameters, _ = curve_fit(self.transmission_spectrum_without_cavity_params,
                                          self.x_axis,
                                          self.current_transmission_spectrum,
                                          p0=[7, x_0_s, 1, 1, 1])
        self.default_parameters.update({"k_ex": optimal_parameters[0],
                                        "x_0": optimal_parameters[1],
                                        "y_0": optimal_parameters[2],
                                        "m": optimal_parameters[3],
                                        "b": optimal_parameters[4]})

    def monitor_spectrum(self):
        self.initialize_default_parameters()
        num_changed = 0
        while True:
            new_parameters = self.fit_transmission_spectrum()
            if np.abs(new_parameters[0] - self.default_parameters["k_ex"]) > 0.1:
                print("k_ex changed", new_parameters[0], self.default_parameters["k_ex"])
                num_changed += 1
                self.plot_transmission_fit()
            else:
                num_changed = 0

            if num_changed > 5:
                self.initialize_default_parameters()
            time.sleep(self.wait_time)

    def plot_transmission_fit(self):
        plt.plot(self.x_axis, self.current_transmission_spectrum)
        k_ex = self.default_parameters["k_ex"]
        k_i, h = self.default_parameters["k_i"], self.default_parameters["h"]
        m, b = self.default_parameters["m"], self.default_parameters["b"]
        x_0, y_0 = self.default_parameters["x_0"], self.default_parameters["y_0"]
        plt.plot(self.x_axis, self.transmission_spectrum(self.x_axis, k_ex, k_i, h, x_0, y_0, m, b))
        plt.show()


if __name__ == '__main__':
    channels_dict = {"transmission": 1, "rubidium": 3}
    res_fit = ResonanceFit(channels_dict=channels_dict, k_i=3.9, h=0.6)
    res_fit.monitor_spectrum()
    thread = Thread(target=res_fit.monitor_spectrum)
    thread.start()
    # transmission_channel = input("Enter the transmission signal channel: ")
    # reflection_channel = input("Enter the reflection signal channel: ")
    # rubidium_channel = input("Enter the rubidium channel: ")
    # resonance_fit = ResonanceFit({"transmission": transmission_channel,
    #                               "reflection": reflection_channel,
    #                               "rubidium": rubidium_channel})
    #
    # input("Scan to find the rubidium lines, then press Enter")
    # resonance_fit.calibrate_x_axis()
    #
    # plt.figure()
    # plt.plot(resonance_fit.rubidium_lines)
    # plt.plot(resonance_fit.rubidium_peaks, resonance_fit.rubidium_lines[resonance_fit.rubidium_peaks], "x")
    # plt.show()
    #
    # input("Set reflection and transmission signal to 0, then press Enter")
    # resonance_fit.set_transmission_0()
    # input("Turn on the transmission signal (without the resonance), then press Enter")
    # resonance_fit.set_transmission_100()
    # input("Bring the resonance back")
    # resonance_fit.read_transmission_spectrum()
    #
    # input("Turn the transmission off, and the reflection on, then press Enter")
    # resonance_fit.read_reflection_spectrum()
    #
    # resonance_fit.fit_current_spectrum()
    # resonance_fit.plot_spectrum_fit()
