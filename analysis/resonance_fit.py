import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from analysis.scope_connection import Scope, FakeScope
from resonance_data import RubidiumLines, CavityFwhm, CavityKex


class ResonanceFit:
    def __init__(self, calc_k_ex=False, lock_idx=2):
        self.cavity = CavityKex(k_i=3.9, h=0.6) if calc_k_ex else CavityFwhm()
        self.lock_idx = lock_idx
        self.rubidium_lines = RubidiumLines()
        self.x_axis = None

        self.current_relevant_area = None
        self.relevant_x_axis = None

        plt.ion()
        self.fig, self.axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
        self.axes[0].set_ylabel("Transmission")
        self.axes[1].set_ylabel("Rubidium")
        self.axes[1].set_xlabel("Frequency [MHz]")

        self.prominence = 0.03
        self.rolling_avg = 100
        self.w_len = 2300
        self.distance = 100
        self.width = 1000

    @property
    def lock_error(self):
        return self.cavity.x_0 - self.x_axis[self.rubidium_lines.peaks_idx[self.lock_idx]]

    def calibrate_peaks_params(self, num_idx_in_peak):
        self.w_len = num_idx_in_peak * 1.1
        self.width = self.w_len // 5
        self.distance = int(num_idx_in_peak // 10)
        self.rolling_avg = int(num_idx_in_peak // 5)

    def calibrate_peaks_params_gui(self, rubidium_lines):
        self.x_axis = np.arange(len(rubidium_lines))
        fig, ax = plt.subplots(figsize=(10, 3))
        ax.plot(self.x_axis, rubidium_lines)
        points = fig.ginput(2, timeout=0, show_clicks=True)
        plt.close(fig)

        num_idx_in_peak = np.round(np.abs(points[0][0] - points[1][0]))
        self.calibrate_peaks_params(num_idx_in_peak)

    def preprocess_data(self, data):
        data -= data.min()
        data /= data.max()
        data = np.convolve(data, np.ones(self.rolling_avg) / self.rolling_avg, mode='valid')
        return data

    def update_rubidium_lines(self, rubidium_lines):
        self.rubidium_lines.data = self.preprocess_data(rubidium_lines)

    def update_transmission_spectrum(self, transmission_spectrum):
        self.cavity.transmission_spectrum = self.preprocess_data(transmission_spectrum)

    def calibrate_x_axis(self) -> bool:
        self.rubidium_lines.peaks_idx, _ = find_peaks(self.rubidium_lines.data, prominence=self.prominence,
                                                      wlen=self.w_len, distance=self.distance, width=self.width)
        if self.rubidium_lines.num_peaks != 6:
            return False
        self.x_axis = np.arange(self.rubidium_lines.num_points) * self.rubidium_lines.idx_to_freq_factor()
        return True

    def calculate_relevant_area(self):
        self.current_relevant_area = (self.cavity.x_0 - 100 < self.x_axis) * (self.x_axis < self.cavity.x_0 + 100)
        self.relevant_x_axis = self.x_axis[self.current_relevant_area]

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

    @staticmethod
    def r2_score(y, f):
        y_bar = y.mean()
        ss_res = ((y - f) ** 2).sum()
        ss_tot = ((y - y_bar) ** 2).sum()
        return 1 - (ss_res / ss_tot)

    def plot_fit(self, title=None):
        title = title or f"{self.cavity.main_parameter}: {self.cavity.value:.2f}, lock error: {self.lock_error:.2f}"
        self.fig.suptitle(title)
        self.plot_transmission_fit(self.axes[0])
        self.plot_rubidium_lines(self.axes[1])

    def plot_transmission_fit(self, ax):
        ax.clear()
        ax.plot(self.x_axis, self.cavity.transmission_spectrum)
        if self.current_relevant_area is not None:
            ax.plot(self.relevant_x_axis, self.cavity.transmission_spectrum_func(self.relevant_x_axis))
        plt.pause(0.05)

    def plot_rubidium_lines(self, ax):
        ax.clear()
        ax.plot(self.x_axis, self.rubidium_lines.data)
        ax.scatter(self.x_axis[self.rubidium_lines.peaks_idx],
                   self.rubidium_lines.data[self.rubidium_lines.peaks_idx], c='r')
        plt.pause(0.05)


class ScopeResonanceFit(ResonanceFit):
    def __init__(self, channels_dict: dict, scope_ip='132.77.54.241', calc_k_ex=False):
        super().__init__(calc_k_ex)

        if scope_ip is None:
            self.scope = FakeScope(channels_dict)
        else:
            self.scope = Scope(ip=scope_ip)

        self.channels_dict = channels_dict

        rubidium_lines = self.read_scope_data("rubidium")[1]
        self.calibrate_peaks_params_gui(rubidium_lines)

    def read_scope_data(self, channel):
        channel_number = self.channels_dict[channel]
        time_axis, data = self.scope.get_data(channel_number)
        data = self.preprocess_data(data)
        return time_axis, data

    # def update_rubidium_lines(self):
    #     _, self.rubidium_lines.data = self.read_scope_data("rubidium")
    #     self.calibrate_x_axis()


class LiveResonanceFit(ScopeResonanceFit):
    def __init__(self, channels_dict: dict, scope_ip='132.77.54.241', wait_time=0.5, save_data=False, calc_k_ex=False):
        super().__init__(channels_dict, scope_ip, calc_k_ex)
        self.wait_time = wait_time
        self.save_data = save_data

        # Assaf ruined your code
        self.save_path = r'C:\temp\refactor_debug\Experiment_results\QRAM\k_ex\k_ex'
        self.last_save_time = time.time()

    def save_parameter(self):
        if not self.save_data:
            return

        now = round(time.time())
        time_since_last_save = now - self.last_save_time
        if time_since_last_save > 10:
            self.last_save_time = now
            np.save(self.save_path, self.cavity.value)

    def monitor_spectrum(self):
        while True:
            transmission_spectrum = self.read_scope_data("transmission")[1]
            rubidium_lines = self.read_scope_data("rubidium")[1]

            self.update_transmission_spectrum(transmission_spectrum)
            self.update_rubidium_lines(rubidium_lines)

            if not self.calibrate_x_axis():
                continue
            if not self.fit_transmission_spectrum():
                continue

            self.plot_fit()
            time.sleep(self.wait_time)

    def start(self):
        self.fig.show()
        self.monitor_spectrum()


if __name__ == '__main__':
    channels = {"transmission": 1, "rubidium": 3}
    res_fit = LiveResonanceFit(channels_dict=channels, scope_ip=None)
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
