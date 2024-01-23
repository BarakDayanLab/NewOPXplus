import os
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from analysis.scope_connection import Scope, FakeScope
from cavity import RubidiumLines, CavityFwhm, CavityKex
import pynput.keyboard as keyboard


class ResonanceFit:
    def __init__(self, calc_k_ex=False, save_folder=None, save_time=60):
        self.cavity = CavityKex(k_i=8, h=0.6) if calc_k_ex else CavityFwhm()
        self.lock_idx = 4
        self.save_folder = save_folder
        self.save_time = save_time
        self.last_save_time = time.time()

        self.rubidium_lines = RubidiumLines()
        self.x_axis = None

        self.current_relevant_area = None
        self.relevant_x_axis = None

        self.fig, self.axes = None, None

        self.prominence = 0.03
        self.rolling_avg = 100
        self.w_len = 2300
        self.distance = 100
        self.width = 1000

    @property
    def lock_error(self):
        return self.cavity.x_0 - self.x_axis[self.rubidium_lines.peaks_idx[self.lock_idx]] + 80

    # ------------------ CALIBRATIONS ------------------ #

    def calibrate_peaks_params(self, num_idx_in_peak):
        self.w_len = num_idx_in_peak * 1.1
        self.width = self.w_len // 5
        self.distance = int(num_idx_in_peak // 10)
        self.rolling_avg = int(num_idx_in_peak // 10)

    def calibrate_peaks_params_gui(self, rubidium_lines):
        self.x_axis = np.arange(len(rubidium_lines))
        fig, ax = plt.subplots(figsize=(10, 3))
        ax.plot(self.x_axis, rubidium_lines)
        points = fig.ginput(2, timeout=0, show_clicks=True)
        plt.close(fig)

        num_idx_in_peak = np.round(np.abs(points[0][0] - points[1][0]))
        self.calibrate_peaks_params(num_idx_in_peak)

    def calibrate_x_axis(self) -> bool:
        self.rubidium_lines.peaks_idx, _ = find_peaks(self.rubidium_lines.data, prominence=self.prominence,
                                                      wlen=self.w_len, distance=self.distance, width=self.width)
        if self.rubidium_lines.num_peaks != 6:
            return False
        self.x_axis = np.arange(self.rubidium_lines.num_points) * self.rubidium_lines.idx_to_freq_factor()
        return True

    def calculate_relevant_area(self):
        self.current_relevant_area = (self.cavity.x_0 - 90 < self.x_axis) * (self.x_axis < self.cavity.x_0 + 90)
        self.relevant_x_axis = self.x_axis[self.current_relevant_area]

    # ------------------ DATA PROCESSING ------------------ #

    def update_rubidium_lines(self, rubidium_lines):
        self.rubidium_lines.data = self.preprocess_data(rubidium_lines)

    def update_transmission_spectrum(self, transmission_spectrum):
        self.cavity.transmission_spectrum = self.preprocess_data(transmission_spectrum)

    def preprocess_data(self, data):
        data -= data.min()
        data /= data.max()
        data = np.convolve(data, np.ones(self.rolling_avg) / self.rolling_avg, mode='valid')
        return data

    # ------------------ FIT ------------------ #
    @staticmethod
    def r2_score(y, f):
        y_bar = y.mean()
        ss_res = ((y - f) ** 2).sum()
        ss_tot = ((y - y_bar) ** 2).sum()
        return 1 - (ss_res / ss_tot)

    def fit_transmission_spectrum(self) -> bool:
        self.cavity.x_0 = self.x_axis[self.cavity.transmission_spectrum.argmin()]
        self.calculate_relevant_area()

        try:
            # noinspection PyTupleAssignmentBalance
            optimal_parameters, covariance = curve_fit(self.cavity.fit,
                                                       self.relevant_x_axis,
                                                       self.cavity.transmission_spectrum[self.current_relevant_area],
                                                       p0=self.cavity.get_fit_parameters(),
                                                       bounds=self.cavity.bounds)
        except Exception as err:
            print(err)
            return False

        self.cavity.set_fit_parameters(*optimal_parameters)
        score = self.r2_score(self.cavity.transmission_spectrum[self.current_relevant_area],
                              self.cavity.transmission_spectrum_func(self.relevant_x_axis))
        if score < 0.6:
            return False
        return True

    # ------------------ PLOT FIT ------------------ #
    def initialize_figure(self):
        plt.ion()
        self.fig, self.axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
        self.axes[0].set_ylabel("Transmission")
        self.axes[1].set_ylabel("Rubidium")
        self.axes[1].set_xlabel("Frequency [MHz]")

    def plot_fit(self, title=None):
        if self.fig is None:
            self.initialize_figure()
        title = title or f"{self.cavity.main_parameter}: {self.cavity.current_fit_value:.2f}, lock error: {self.lock_error:.2f}, $\Delta$x: {self.cavity.x_0 - self.cavity.minimum_of_fit:.2f}"
        self.fig.suptitle(title)
        self.plot_transmission_fit(self.axes[0])
        self.plot_rubidium_lines(self.axes[1])

    def plot_transmission_fit(self, ax):
        ax.clear()
        ax.scatter(self.cavity.x_0, self.cavity.transmission_spectrum_func(self.cavity.x_0), c='g')
        ax.plot(self.x_axis, self.cavity.transmission_spectrum)
        ax.plot(self.relevant_x_axis, self.cavity.transmission_spectrum_func(self.relevant_x_axis))
        plt.pause(0.05)

    def plot_rubidium_lines(self, ax):
        ax.clear()
        ax.scatter(self.x_axis[self.rubidium_lines.peaks_idx],
                   self.rubidium_lines.data[self.rubidium_lines.peaks_idx], c='r')
        ax.scatter(self.x_axis[self.rubidium_lines.peaks_idx[self.lock_idx]],
                   self.rubidium_lines.data[self.rubidium_lines.peaks_idx[self.lock_idx]], c='g')
        ax.plot(self.x_axis, self.rubidium_lines.data)
        plt.pause(0.05)

    # ------------------ SAVE DATA ------------------ #
    def get_save_paths(self):
        date = time.strftime("%Y%m%d")
        hours = time.strftime("%H%M%S")

        transmission_filename = f"{date}-{hours}_cavity_spectrum.npy"
        transmission_path = os.path.join(self.save_folder, date, transmission_filename)

        rubidium_filename = f"{date}-{hours}_rubidium_spectrum.npy"
        rubidium_path = os.path.join(self.save_folder, date, rubidium_filename)
        return transmission_path, rubidium_path

    def save_spectrum(self):
        if self.save_folder is None:
            return

        now = time.time()
        time_since_last_save = now - self.last_save_time
        if time_since_last_save >= self.save_time:
            transmission_path, rubidium_path = self.get_save_paths()
            np.save(transmission_path, self.cavity.transmission_spectrum)
            np.save(rubidium_path, self.rubidium_lines.data)
            self.last_save_time = now


class LiveResonanceFit(ResonanceFit):
    def __init__(self, wait_time=0.5, calc_k_ex=False, save_folder=None, save_time=60, plot=True):
        super().__init__(calc_k_ex=calc_k_ex, save_folder=save_folder, save_time=save_time)
        self.wait_time = wait_time
        self.plot = plot

        self.pause = False
        self.keyboard_listener = keyboard.Listener(on_press=self.keyboard_on_press)
        self.keyboard_listener.start()

        self.stop_condition = False

    # ------------------ KEYBOARD INTERRUPTS ------------------ #
    def keyboard_on_press(self, key):
        if key == keyboard.Key.space:
            self.pause = not self.pause

    def wait(self):
        while self.pause:
            plt.pause(0.1)

    # ------------------ READ DATA ------------------ #
    def read_transmission_spectrum(self):
        raise Exception("read_transmission_spectrum not implemented")

    def read_rubidium_lines(self):
        raise Exception("read_rubidium_lines not implemented")

    # ------------------ MAIN LOOP ------------------ #
    def main_loop_callback(self):
        pass

    def start(self):
        while not self.stop_condition:
            self.wait()
            transmission_spectrum = self.read_transmission_spectrum()
            rubidium_lines = self.read_rubidium_lines()
            if transmission_spectrum is None or rubidium_lines is None:
                time.sleep(self.wait_time)
                continue

            self.update_transmission_spectrum(transmission_spectrum)
            self.update_rubidium_lines(rubidium_lines)

            if not self.calibrate_x_axis():
                continue
            if not self.fit_transmission_spectrum():
                continue

            self.main_loop_callback()

            if self.plot:
                self.plot_fit()
                time.sleep(self.wait_time)


class ScopeResonanceFit(LiveResonanceFit):
    def __init__(self, channels_dict: dict, scope_ip='132.77.54.241', save_data=False):
        save_folder = r'C:\temp\refactor_debug\Experiment_results\QRAM\resonance_params' if save_data else None
        super().__init__(wait_time=0.5, calc_k_ex=False, save_folder=save_folder, save_time=15)
        self.save_data = save_data

        if scope_ip is None:
            self.scope = FakeScope(channels_dict)
        else:
            self.scope = Scope(ip=scope_ip)

        self.channels_dict = channels_dict

        rubidium_lines = self.read_scope_data("rubidium")[1]
        self.calibrate_peaks_params_gui(rubidium_lines)

        # Assaf ruined your code
        self.save_path = r'C:\temp\refactor_debug\Experiment_results\QRAM\resonance_params'
        self.last_save_time_k_ex = time.time()

    # ------------------ READ DATA ------------------ #
    def read_scope_data(self, channel):
        channel_number = self.channels_dict[channel]
        time_axis, data = self.scope.get_data(channel_number)
        return time_axis, data

    def read_transmission_spectrum(self):
        return self.read_scope_data("transmission")[1]

    def read_rubidium_lines(self):
        return self.read_scope_data("rubidium")[1]

    # ------------------ SAVE DATA ------------------ #
    def save_parameter(self):
        if not self.save_data:
            return

        now = round(time.time())
        time_since_last_save = now - self.last_save_time_k_ex
        if time_since_last_save > 3:
            self.last_save_time_k_ex = now
            fwhm_path = os.path.join(self.save_path, "k_ex")
            lock_path = os.path.join(self.save_path, "locking_err")
            print(self.cavity.current_fit_value)
            np.save(fwhm_path, self.cavity.current_fit_value)
            print(self.lock_error)
            np.save(lock_path, self.lock_error)

    def main_loop_callback(self):
        self.save_parameter()


if __name__ == '__main__':
    channels = {"transmission": 1, "rubidium": 3}
    res_fit = ScopeResonanceFit(channels_dict=channels, save_data=False)
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
