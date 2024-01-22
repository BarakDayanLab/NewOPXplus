from resonance_fit import ResonanceFit
import numpy as np
import os
import re
import time
import pynput.keyboard as keyboard
import matplotlib.pyplot as plt


class FolderResonanceFit(ResonanceFit):
    cavity_spectrum_regex = re.compile(r"[0-9]+-([0-9]+)_cavity_spectrum.npy")
    rubidium_lines_regex = re.compile(r"[0-9]+-([0-9]+)_rb_lines_spectrum.npy")

    def __init__(self, folder_path, start_time=140000, end_time=180000):
        super().__init__()
        self.folder_path = folder_path
        self.start_time = start_time
        self.end_time = end_time
        self.data = self.data_loader()

        rubidium_lines = next(self.data)[1]
        self.calibrate_peaks_params_gui(rubidium_lines)

        self.pause = False
        self.keyboard_listener = keyboard.Listener(on_press=self.keyboard_on_press)
        self.keyboard_listener.start()

    def keyboard_on_press(self, key):
        if key == keyboard.Key.space:
            self.pause = not self.pause

    def data_loader(self):
        cavity_spectrum_files = [m for f in os.listdir(self.folder_path)
                                 if (m := self.cavity_spectrum_regex.fullmatch(f))]
        rubidium_lines_files = [m for f in os.listdir(self.folder_path)
                                if (m := self.rubidium_lines_regex.fullmatch(f))]

        cavity_spectrum_files = filter(lambda m: self.start_time <= int(m[1]) <= self.end_time, cavity_spectrum_files)
        rubidium_lines_files = filter(lambda m: self.start_time <= int(m[1]) <= self.end_time, rubidium_lines_files)

        cavity_spectrum_files = sorted(cavity_spectrum_files, key=lambda m: int(m[1]))
        rubidium_lines_files = sorted(rubidium_lines_files, key=lambda m: int(m[1]))

        for cavity_spectrum_file, rubidium_lines_file in zip(cavity_spectrum_files, rubidium_lines_files):
            cavity_spectrum = np.load(os.path.join(self.folder_path, cavity_spectrum_file[0]))
            rubidium_lines = np.load(os.path.join(self.folder_path, rubidium_lines_file[0]))
            yield -cavity_spectrum, rubidium_lines

    def wait(self):
        while self.pause:
            plt.pause(0.1)

    def run(self):
        for cavity_spectrum, rubidium_lines in self.data:
            self.wait()
            self.update_transmission_spectrum(cavity_spectrum)
            self.update_rubidium_lines(rubidium_lines)
            if not self.calibrate_x_axis():
                continue
            if not self.fit_transmission_spectrum():
                continue

            self.plot_fit()
            self.fig.show()
            time.sleep(1)


if __name__ == '__main__':
    path = r"U:\Lab_2023\Experiment_results\QRAM\Cavity_Spectrum\20240118"
    folder_resonance_fit = FolderResonanceFit(path, start_time=140000, end_time=180000)
    folder_resonance_fit.run()
