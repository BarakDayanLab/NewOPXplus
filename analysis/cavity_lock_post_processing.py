from resonance_fit import LiveResonanceFit
import numpy as np
import os
import re
import matplotlib.pyplot as plt


class FolderResonanceFit(LiveResonanceFit):
    cavity_spectrum_regex = re.compile(r"[0-9]+-([0-9]+)_cavity_spectrum.npy")
    rubidium_lines_regex = re.compile(r"[0-9]+-([0-9]+)_rb_lines_spectrum.npy")

    def __init__(self, folder_path, start_time=140000, end_time=180000):
        super().__init__(wait_time=1, plot=True, calc_k_ex=True)
        self.folder_path = folder_path
        self.start_time = start_time
        self.end_time = end_time

        rubidium_lines = next(self.rubidium_lines_generator())
        self.calibrate_peaks_params_gui(rubidium_lines)

        self.transmission_generator = self.transmission_spectrum_generator()
        self.rubidium_lines_generator = self.rubidium_lines_generator()

        self.main_parameter_history = []
        self.lock_error_history = []
        self.diff_history = []

    def transmission_spectrum_generator(self):
        cavity_spectrum_files = [m for f in os.listdir(self.folder_path)
                                 if (m := self.cavity_spectrum_regex.fullmatch(f))]

        cavity_spectrum_files = filter(lambda m: self.start_time <= int(m[1]) <= self.end_time, cavity_spectrum_files)
        cavity_spectrum_files = sorted(cavity_spectrum_files, key=lambda m: int(m[1]))

        for cavity_spectrum_file in cavity_spectrum_files:
            cavity_spectrum = np.load(os.path.join(self.folder_path, cavity_spectrum_file[0]))
            yield -cavity_spectrum

    def rubidium_lines_generator(self):
        rubidium_lines_files = [m for f in os.listdir(self.folder_path)
                                if (m := self.rubidium_lines_regex.fullmatch(f))]

        rubidium_lines_files = filter(lambda m: self.start_time <= int(m[1]) <= self.end_time, rubidium_lines_files)
        rubidium_lines_files = sorted(rubidium_lines_files, key=lambda m: int(m[1]))

        for rubidium_lines_file in rubidium_lines_files:
            rubidium_lines = np.load(os.path.join(self.folder_path, rubidium_lines_file[0]))
            yield rubidium_lines

    def read_transmission_spectrum(self):
        transmission_spectrum = next(self.transmission_generator, None)
        self.stop_condition = transmission_spectrum is None
        return transmission_spectrum

    def read_rubidium_lines(self):
        rubidium_spectrum = next(self.rubidium_lines_generator, None)
        self.stop_condition = rubidium_spectrum is None
        return rubidium_spectrum

    def main_loop_callback(self):
        self.main_parameter_history.append(self.cavity.current_fit_value)
        self.lock_error_history.append(self.lock_error)
        self.diff_history.append(self.cavity.x_0 - self.cavity.minimum_of_fit)

    def start(self):
        super().start()
        plt.plot(self.main_parameter_history)
        plt.show()
        plt.plot(self.lock_error_history)
        plt.show()
        plt.plot(self.diff_history)
        plt.show()


if __name__ == '__main__':
    path = r"U:\Lab_2023\Experiment_results\QRAM\Cavity_Spectrum\20240118"
    folder_resonance_fit = FolderResonanceFit(path, start_time=140000, end_time=180000)
    folder_resonance_fit.start()
