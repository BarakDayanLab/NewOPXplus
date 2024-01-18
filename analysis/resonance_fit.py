import sys
import time
import logging
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, convolve, find_peaks_cwt
from scipy.optimize import curve_fit
from analysis.scope_connection import Scope, FakeScope

plt.ion()


class BaseResonanceFit:
    def __init__(self, channels_dict: dict, scope_ip='132.77.54.241', rolling_avg=100, wait_time=0.5):
        if scope_ip is None:
            self.scope = FakeScope(channels_dict)
        else:
            self.scope = Scope(ip=scope_ip)
        self.channels_dict = channels_dict
        self.rolling_avg = rolling_avg
        self.wait_time = wait_time

        self.x_axis = None

        format = logging.Formatter('%(levelname)s - %(asctime)s - %(message)s')
        handler = logging.FileHandler("resonance_fit.log")
        handler.setLevel(logging.WARNING)
        handler.setFormatter(format)
        self.logger = logging.getLogger("resonance_fit_logger")
        self.logger.addHandler(handler)

    def read_scope_data(self, channel):
        channel_number = self.channels_dict[channel]
        time_axis, data = self.scope.get_data(channel_number)
        data = self.normalize_data(data)
        return time_axis, data

    def normalize_data(self, data):
        data -= data.min()
        data /= data.max()

        data = np.convolve(data, np.ones(self.rolling_avg) / self.rolling_avg, mode='valid')
        return data

    def calibrate_x_axis(self):
        _, rubidium_lines = self.read_scope_data("rubidium")
        rubidium_peaks, _ = find_peaks(rubidium_lines, prominence=0.03, wlen=2300, distance=500)

        if len(rubidium_peaks) != 6:
            self.logger.error("Could not find 6 peaks in the rubidium spectrum")
            time.sleep(self.wait_time)
            return self.calibrate_x_axis()

        idx_to_freq = (156.947 / 2) / (rubidium_peaks[-1] - rubidium_peaks[-2])
        x_axis = np.arange(len(rubidium_lines)) * idx_to_freq  # Calibration
        self.x_axis = x_axis

    @staticmethod
    def transmission_spectrum(x_detuning, k_ex, k_i, h, x_0, y_0):
        k_total = k_ex + k_i
        k_with_detuning = k_total + 1j * (x_detuning - x_0)
        return np.abs(1 - 2 * k_ex * k_with_detuning / (k_with_detuning ** 2 + h ** 2)) ** 2 + y_0

    @staticmethod
    def reflection_spectrum(x_detuning, k_ex, k_i, h, x_0, y_0):
        k_total = k_ex + k_i
        k_with_detuning = k_total + 1j * (x_detuning - x_0)
        return np.abs(2 * k_ex * h / (k_with_detuning ** 2 + h ** 2)) ** 2 + y_0


class ResonanceFit(BaseResonanceFit):
    def __init__(self, channels_dict: dict, scope_ip='132.77.54.241', rolling_avg=50, wait_time=0.5):
        super().__init__(channels_dict, scope_ip)

    def set_transmission_0(self):
        self.transmission_0 = self.read_scope_data("transmission")

    def set_transmission_100(self):
        self.transmission_100 = self.read_scope_data("transmission")

    def stacked_spectrum(self, x_detuning, k_ex, k_i, h, x_0, y_0):
        return np.hstack((self.transmission_spectrum(x_detuning, k_ex, k_i, h, x_0, y_0),
                          self.reflection_spectrum(x_detuning, k_ex, k_i, h, x_0, y_0)))

    def read_reflection_spectrum(self):
        current_reflection = self.read_scope_data("reflection")
        reflection_norm_factor = 1 / (self.transmission_100 - self.transmission_0)
        reflection_spectrum = reflection_norm_factor * (current_reflection - self.transmission_0)
        return reflection_spectrum

    def fit_current_spectrum(self, transmission_spectrum, reflection_spectrum):
        y_data = np.hstack((transmission_spectrum, reflection_spectrum))
        # noinspection PyTupleAssignmentBalance
        optimal_parameters, _ = curve_fit(self.stacked_spectrum, self.x_axis, y_data)
        return optimal_parameters


class LiveResonanceFit(BaseResonanceFit):
    def __init__(self, channels_dict: dict, k_i: float, h: float, scope_ip='132.77.54.241'):
        super().__init__(channels_dict, scope_ip)
        self.k_i = k_i
        self.h = h

        self.x_0 = 0
        self.y_0 = 0
        self.k_ex = 0

        self.current_k_ex = 0
        self.current_x_0 = 0

        self.current_relevant_area = None
        self.relevant_x_axis = None

        self.fig, self.ax = plt.subplots()

        #Assaf ruined your code
        self.k_ex_error_path = r'C:\temp\refactor_debug\Experiment_results\QRAM\k_ex\k_ex'
        self.k_ex_last_save_time = time.time()

    def start(self):
        self.fig.show()
        self.monitor_spectrum()

    def transmission_spectrum_without_cavity_params(self, x_detuning, k_ex, y_0):
        return self.transmission_spectrum(x_detuning, k_ex, self.k_i, self.h, self.x_0, y_0)

    def transmission_spectrum_k_ex(self, x_detuning, k_ex):
        return self.transmission_spectrum(x_detuning, k_ex, self.k_i, self.h, self.current_x_0, self.y_0)

    def calculate_relevant_area(self):
        self.current_relevant_area = (self.x_0 - 200 < self.x_axis) * (self.x_axis < self.x_0 + 200)
        self.relevant_x_axis = self.x_axis[self.current_relevant_area]

    def initialize_default_parameters(self):
        self.calibrate_x_axis()
        _, transmission_spectrum = self.read_scope_data("transmission")
        self.x_0 = self.x_axis[transmission_spectrum.argmin()]

        self.calculate_relevant_area()
        transmission_spectrum = transmission_spectrum[self.current_relevant_area]

        # noinspection PyTupleAssignmentBalance
        optimal_parameters, _ = curve_fit(self.transmission_spectrum_without_cavity_params,
                                          self.relevant_x_axis,
                                          transmission_spectrum,
                                          p0=[30, 1])

        self.k_ex, self.y_0 = optimal_parameters

    def fit_transmission_spectrum(self, transmission_spectrum):
        # noinspection PyTupleAssignmentBalance
        optimal_parameters, _ = curve_fit(self.transmission_spectrum_k_ex,
                                          self.relevant_x_axis,
                                          transmission_spectrum,
                                          p0=[self.k_ex])

        return optimal_parameters[0]

    def test_x_0(self, transmission_spectrum):
        self.current_x_0 = self.x_axis[transmission_spectrum.argmin()]
        if np.abs(self.current_x_0 - self.x_0) > 5:
            self.logger.warning(f"Measured different x_0, new: {self.current_x_0}, current: {self.x_0}")
            return 1
        return 0

    def test_k_ex(self, transmission_spectrum):
        self.current_k_ex = self.fit_transmission_spectrum(transmission_spectrum)
        if np.abs(self.current_k_ex - self.k_ex) > 5:
            self.logger.warning(f"Measured different k_ex, new: {self.current_k_ex}, current: {self.k_ex}")
            return 1
        return 0

    def monitor_spectrum(self):
        self.initialize_default_parameters()
        logging.warning(f"Initialized k_ex: {self.k_ex}")

        num_k_ex_changed = 0
        num_x_0_changed = 0
        while True:
            _, transmission_spectrum = self.read_scope_data("transmission")

            num_x_0_changed += self.test_x_0(transmission_spectrum)

            self.calculate_relevant_area()
            transmission_spectrum = transmission_spectrum[self.current_relevant_area]

            num_k_ex_changed += self.test_k_ex(transmission_spectrum)

            self.plot_transmission_fit(transmission_spectrum)

            if num_k_ex_changed > 5:
                self.logger.error("k_ex changed too many times, reinitializing")
                self.initialize_default_parameters()
                num_k_ex_changed = 0

            elif num_x_0_changed > 5:
                self.logger.error("x_0 changed too many times, reinitializing")
                self.initialize_default_parameters()
                num_x_0_changed = 0

            now = round(time.time())
            time_passed = now - self.k_ex_last_save_time
            # print('%.3f' %time_passed)
            if time_passed > 10:  # 2 seconds (or more) passed since last write?
                self.k_ex_last_save_time = now
                # try:
                # Save the current error in a file - for the Control code to take
                np.save(self.k_ex_error_path, self.current_k_ex)
                # # Save the current error in a dated file
                # file_name = self.prepare_file_name('Locking_PID_Error', 'locking_err_log', 'txt', True)
                # with open(file_name, 'a') as f:
                #     f.write("%s: %s\n" % (time.strftime("%H:%M:%S"), str(errorSignal)))

            time.sleep(self.wait_time)

    def plot_transmission_fit(self, transmission_spectrum):
        self.ax.clear()
        self.ax.set_title(f"k_ex: {self.current_k_ex:.2f}")
        self.ax.plot(self.relevant_x_axis, transmission_spectrum)
        self.ax.plot(self.relevant_x_axis, self.transmission_spectrum(self.relevant_x_axis, self.current_k_ex,
                                                                      self.k_i, self.h, self.current_x_0, self.y_0))
        plt.pause(0.05)



if __name__ == '__main__':
    channels = {"transmission": 1, "rubidium": 3}
    res_fit = LiveResonanceFit(channels_dict=channels, k_i=3.9, h=0.6)
    res_fit.start()
    # _, transmission_spectrum = res_fit.read_scope_data("transmission")
    #
    # transmission = np.zeros((100, transmission_spectrum.shape[0]))
    # rubidium = np.zeros((100, transmission_spectrum.shape[0]))
    # for i in range(100):
    #     _, transmission_spectrum = res_fit.read_scope_data("transmission")
    #     _, rubidium_spectrum = res_fit.read_scope_data("rubidium")
    #     transmission[i] = transmission_spectrum
    #     rubidium[i] = rubidium_spectrum
    #     time.sleep(1.15)
    #
    # np.save("transmission.npy", transmission)
    # np.save("rubidium.npy", rubidium)


    # res_fit.monitor_spectrum()
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
