import os
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from analysis.scope_connection import Scope, FakeScope
from cavity import RubidiumLines, CavityFwhm, CavityKex
from pynput.keyboard import GlobalHotKeys, Key


class ResonanceFit:
    def __init__(self, calc_k_ex=False, save_folder=None, save_time=60):
        self.cavity = CavityKex(k_i=4.6, h=1.7) if calc_k_ex else CavityFwhm()
        self.lock_idx = 4
        self.save_folder = save_folder
        self.save_time = save_time
        self.last_save_time = time.time()

        self.rubidium_lines = RubidiumLines()
        self.x_axis = None

        self.current_relevant_area = None
        self.relevant_x_axis = None

        self.fit_params_history = None

        self.prominence = 0.03
        self.rolling_avg = 100
        self.w_len = 2300
        self.distance = 100
        self.width = 1000

    @property
    def lock_error(self):
        return self.cavity.x_0 - self.x_axis[self.rubidium_lines.peaks_idx[self.lock_idx]] + 80

    @property
    def rubidium_peaks(self):
        x_vals = self.x_axis[self.rubidium_lines.peaks_idx]
        y_vals = self.rubidium_lines.data[self.rubidium_lines.peaks_idx]
        peaks = np.vstack([x_vals, y_vals]).T
        return peaks

    @property
    def lorentzian_center(self):
        return np.array([self.cavity.x_0, self.cavity.transmission_spectrum_func(self.cavity.x_0)])

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


class ResonanceFitGraphics:
    def __init__(self, resonance_fit: ResonanceFit, show_buttons=True):
        self.fig, self.plot_subfigure, self.buttons_subfigure = None, None, None
        self.show_buttons = show_buttons
        self.buttons = {}
        self.current_active_button = None

        self.resonance_fit = resonance_fit

        self.show_3rd_axis = False
        self.function_3rd_axis = None

    # ------------------ UI ------------------ #
    def update_active_button(self, button_name):
        self.current_active_button = button_name

    def activate_button(self):
        if self.current_active_button == "":
            return

        if self.current_active_button == "choose_line_button":
            self.choose_line()
        elif self.current_active_button == "show_lock_error_button":
            self.toggle_lock_error()
        self.current_active_button = ""

    def choose_line(self):
        prev_title = self.fig.texts[0].get_text()
        self.fig.suptitle("Choose line")

        point = plt.ginput(1, timeout=0, show_clicks=True)[0]
        distances = np.sum((self.resonance_fit.rubidium_peaks - point) ** 2, axis=1)
        self.resonance_fit.lock_idx = np.argmin(distances)

        self.fig.suptitle(prev_title)

    def toggle_lock_error(self):
        if not self.show_3rd_axis:
            self.initialize_plots()
            self.function_3rd_axis = None
        else:
            self.initialize_plots(plot_3rd_axis=("Lock error", "MHz"))
            self.function_3rd_axis = self.plot_lock_error

        self.show_3rd_axis = not self.show_3rd_axis

    # ------------------ PLOT FIT ------------------ #
    def initialize_figure(self):
        plt.ion()
        self.fig = plt.figure(figsize=(16, 9), constrained_layout=True)
        if self.show_buttons:
            self.plot_subfigure, self.buttons_subfigure = self.fig.subfigures(2, 1, height_ratios=[94, 6])
        else:
            self.plot_subfigure = self.fig

        self.initialize_plots()
        if self.show_buttons:
            self.buttons_axis()

    def initialize_plots(self, plot_3rd_axis=None):
        n_rows = 2 if plot_3rd_axis is None else 3
        self.plot_subfigure.clear()
        self.plot_subfigure.add_subplot(n_rows, 1, 1)
        self.plot_subfigure.add_subplot(n_rows, 1, 2, sharex=self.plot_subfigure.axes[0])
        self.plot_subfigure.axes[0].set_ylabel("Transmission")
        self.plot_subfigure.axes[1].set_ylabel("Rubidium")
        self.plot_subfigure.axes[1].set_xlabel("Frequency [MHz]")

        if plot_3rd_axis is not None:
            self.plot_subfigure.add_subplot(3, 1, 3)
            self.plot_subfigure.axes[2].set_title(plot_3rd_axis[0])
            self.plot_subfigure.axes[2].set_ylabel(plot_3rd_axis[1])

    def buttons_axis(self):
        self.buttons_subfigure.subplots(1, 2)

        choose_line_button = plt.Button(self.buttons_subfigure.axes[0], 'choose line')
        self.buttons.update({"choose_line_button": choose_line_button})
        choose_line_button.on_clicked(lambda _: self.update_active_button("choose_line_button"))

        show_lock_error_button = plt.Button(self.buttons_subfigure.axes[1], 'lock error')
        self.buttons.update({"show_lock_error_button": show_lock_error_button})
        show_lock_error_button.on_clicked(lambda _: self.update_active_button("show_lock_error_button"))

    def plot_fit(self, title=None):
        if self.fig is None:
            self.initialize_figure()

        default_title = (f"{self.resonance_fit.cavity.main_parameter}: {self.resonance_fit.cavity.current_fit_value:.2f}, "
                         f"lock error: {self.resonance_fit.lock_error:.2f}, "
                         f"$\Delta$x: {self.resonance_fit.cavity.x_0 - self.resonance_fit.cavity.minimum_of_fit:.2f}")

        title = title or default_title
        self.fig.suptitle(title)
        self.plot_transmission_fit(self.plot_subfigure.axes[0])
        plt.pause(0.05)
        self.plot_rubidium_lines(self.plot_subfigure.axes[1])
        plt.pause(0.05)
        self.show_3rd_axis and self.plot_lock_error(self.plot_subfigure.axes[2])
        plt.pause(0.05)

    def plot_transmission_fit(self, ax):
        ax.clear()

        spectrum = self.resonance_fit.cavity.transmission_spectrum
        fit = self.resonance_fit.cavity.transmission_spectrum_func(self.resonance_fit.relevant_x_axis)
        ax.scatter(*self.resonance_fit.lorentzian_center, c='g')
        ax.plot(self.resonance_fit.x_axis, spectrum)
        ax.plot(self.resonance_fit.relevant_x_axis, fit)
        plt.pause(0.05)

    def plot_rubidium_lines(self, ax):
        ax.clear()
        rubidium_peaks = self.resonance_fit.rubidium_peaks
        ax.scatter(*rubidium_peaks.T, c='r')
        ax.scatter(*rubidium_peaks[self.resonance_fit.lock_idx], c='g')
        ax.plot(self.resonance_fit.x_axis, self.resonance_fit.rubidium_lines.data)

    def plot_lock_error(self, ax):
        ax.clear()
        ax.plot(self.resonance_fit.fit_params_history[-20:, -1])


class LiveResonanceFit(ResonanceFit):
    def __init__(self, wait_time=0.1, calc_k_ex=False, save_folder=None, save_time=60, show_buttons=False):
        super().__init__(calc_k_ex=calc_k_ex, save_folder=save_folder, save_time=save_time)
        self.wait_time = wait_time
        self.fit_params_history = np.empty((0, len(self.cavity.fit_parameters)+1))

        hotkeys = {"<ctrl>+p": self.toggle_pause}
        self.keyboard_listener = GlobalHotKeys(hotkeys)
        self.keyboard_listener.start()
        self.pause = False

        self.stop_condition = False

        self.graphics = ResonanceFitGraphics(self, show_buttons=show_buttons)

    # ------------------ KEYBOARD INTERRUPTS ------------------ #
    def toggle_pause(self):
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

            params = self.cavity.get_fit_parameters() + [self.lock_error]
            self.fit_params_history = np.vstack([self.fit_params_history, params])

            self.main_loop_callback()

            if hasattr(self, "graphics"):
                self.graphics.plot_fit()
                self.graphics.activate_button()
                time.sleep(self.wait_time)


class ScopeResonanceFit(LiveResonanceFit):
    def __init__(self, channels_dict: dict, scope_ip='132.77.54.241', save_data=False):
        save_folder = r'C:\temp\refactor_debug\Experiment_results\QRAM\resonance_params' if save_data else None
        super().__init__(wait_time=0.1, calc_k_ex=False, save_folder=save_folder, save_time=15, show_buttons=False)
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
    res_fit = ScopeResonanceFit(channels_dict=channels, save_data=True)
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
