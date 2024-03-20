from Experiments.BaseExperiment.BaseExperiment import BaseExperiment
from Experiments.Enums.ExperimentMode import ExperimentMode as ExperimentMode

from Utilities.Utils import Utils

import numpy as np
from PIL import Image
import time
import subprocess
import json
import os
from pathlib import Path
import cv2, glob
import scipy.optimize as opt
import matplotlib
import matplotlib.pyplot as mtl
import pylab as plt
from Utilities.BDMenu import BDMenu
from mpl_toolkits.mplot3d import Axes3D

_Kb = 1.38e-23  #[J/K]
_m = 1.443e-25  #[Kg] of Rb87


class CoolingSequenceOptimizer(BaseExperiment):

    def __init__(self, experiment_mode=ExperimentMode.LIVE):
        # Invoking BaseClass constructor. It will initiate OPX, QuadRF, BDLogger, Camera, BDResults, KeyEvents etc.
        super().__init__(playback_parameters=None, save_raw_data=False, connect_to_camera=True, experiment_mode=experiment_mode)
        pass

    def initialize_experiment_variables(self):

        # Ensure we initialize the basic OPX experiment variables (required for standard MOT)
        super().initialize_experiment_variables()

        # self.SPRINT_Exp_switch(False) # To enable trigger to camera in OPX_control_with_QuadRF_Sprint_Exp
        self.camera = None
        self.NAvg = 1  # Number of photos captured to create an average image
        self.NThrow = 3  # Number of image throwed to garbage we're "skipping" at the begging of each capturing time
        self.NThrow_end = 1 # Number of image throwed to garbage we're "skipping" at the end of each capturing time
        self.imgBounds = (580, 200, 1600, 1450)  # bounds to crop out of the taken pictures
        # self.mm_to_pxl = 8.5/(830-56)  # measured using ruler in focus 13/11/2022
        self.mm_to_pxl = 11/(1228-232)  # measured using ruler in focus 15/01/2024
        self.sigma_bounds = (15, 100)  # This bounds sigma (x & y) of the Gaussian sigma. If value is out of bounds, fit is considered bad and not used in temp-fit

    def connect_disconnect_camera(self):
        if self.camera:
            self.disconnect_camera()
            self.info('Disconnected from MV camera.')
        else:
            success = self.connect_camera()
            if success:
                self.info('Connected to MV camera!')
            else:
                self.info('Attempted to connect to MV camera, but failed. Check if app is in use.')

    def createPathForTemperatureMeasurement(self, path=None, saveConfig=False):
        extraFilesPath = ''
        if path:
            extraFilesPath = os.path.join(path, 'extra_files')
        else:
            # --- create path and save config -----
            path = self.bd_results.get_root()
            extraFilesPath = os.path.join(path, 'extra_files')
            if not os.path.exists(path):
                os.makedirs(path)
            if not os.path.exists(extraFilesPath):
                os.makedirs(extraFilesPath)
            # if saveConfig:
            #     self.save_config_table(path=extraFilesPath)
        return (path, extraFilesPath)

    def fitVy_0FromGaussianFitResults(self, gaussianFitResult, extraFilesPath,fit_for_alpha=False, plotResults = True):
        time_vector = np.array([res[0] for res in gaussianFitResult])
        y_position_vector = np.array([res[-1][2] for res in gaussianFitResult])

        mm_to_pxl = self.mm_to_pxl
        linearFreeFall = lambda t, b, c: 9.8e-6 / 2 / mm_to_pxl * 1e3 * t ** 2 + b * t + c
        quadraticFunc = lambda t, a, b, c: a * t ** 2 + b * t + c
        fitFunc = quadraticFunc if fit_for_alpha else linearFreeFall
        bounds = (-np.inf, np.inf)
        # bounds = (-np.inf, np.inf) if fit_for_alpha else ([9.8e-6/2/self.mm_to_pxl * 1e3, -np.inf,0],[9.8001e-6/2/self.mm_to_pxl * 1e3, np.inf, np.inf])   # fit for the acceleration of the quadratic
        # if fit_for_alpha is FALSE, set acceleration to BE g. not trikim no shtikim

        v_launch_popt, v_launch_cov = opt.curve_fit(fitFunc, time_vector, y_position_vector, bounds=bounds)
        alpha = 9.8e-6 / 2 / v_launch_popt[0] * 1e3 if fit_for_alpha else self.mm_to_pxl  # mm/pixel
        v_launch = v_launch_popt[0] * self.mm_to_pxl  # mm/ms = m/s
        y_position_vector_mm = y_position_vector * self.mm_to_pxl - (y_position_vector[0] - 10) * self.mm_to_pxl  # -10 due to resonator position

        if plotResults:
            titlestr_v = r'Y position fit \n $v_0 = %.2f$' % v_launch + r'[cm/s]; $\alpha = %.3f$' % alpha + '[mm/pixel]'
            resstr_v = r'$v_0 = %.2f$' % v_launch + r'[cm/s]; $\alpha = %.3f$' % alpha
            self.plotDataAndFit(time_vector, y_position_vector, fitFunc=fitFunc, fitParams=v_launch_popt,
            # self.plotDataAndFit(time_vector, y_position_vector_mm, fitFunc=fitFunc, fitParams=v_launch_popt,
                                # title=f'Y position fit \n V_launch = {v_launch}; alpha = {alpha}', ylabel='Y_center [px]',
                                title=titlestr_v, props_str=resstr_v, ylabel='$Y_{center} [px]$',
                                # title=titlestr_v, props_str=resstr_v, ylabel='$Y_{center} [mm]$',
                                saveFilePath=os.path.join(extraFilesPath, 'Y_position.png'), show=False)
        return (v_launch, alpha, v_launch_popt, v_launch_cov)

    def fitTemperaturesFromGaussianFitResults(self, gaussianFitResult, alpha=None, extraFilesPath=None, plotResults=True):
        if alpha is None: alpha = self.mm_to_pxl
        time_vector = np.array([res[0] for res in gaussianFitResult])
        sigma_x_vector = alpha * np.array([res[-1][3] for res in gaussianFitResult])
        sigma_y_vector = alpha * np.array([res[-1][4] for res in gaussianFitResult])
        return self.fitTemperaturesFromSigmasResults(time_vector, sigma_x_vector, sigma_y_vector, alpha = None, extraFilesPath = extraFilesPath, plotResults = plotResults)

    def fitTemperaturesFromSigmasResults(self, time_vector, sigma_x_vector, sigma_y_vector, alpha = None, extraFilesPath = None, plotResults = True):
        if alpha:
            sigma_x_vector = alpha * np.array(sigma_x_vector)
            sigma_y_vector = alpha * np.array(sigma_y_vector)
        tempFromSigmaFunc = lambda t, a, b: np.sqrt(a + b * t ** 2)
        x_temp_popt, x_temp_cov = opt.curve_fit(tempFromSigmaFunc, time_vector, sigma_x_vector, bounds=(0, np.inf))
        y_temp_popt, y_temp_cov = opt.curve_fit(tempFromSigmaFunc, time_vector, sigma_y_vector, bounds=(0, np.inf))
        T_x, T_y = x_temp_popt[1] * _m / _Kb * 1e6, y_temp_popt[1] * _m / _Kb * 1e6  # [uK]
        std_x, std_y = np.sqrt(np.diag(x_temp_cov))[1] * _m / _Kb * 1e6, np.sqrt(np.diag(x_temp_cov))[1] * _m / _Kb * 1e6

        if plotResults:
            titlestr_x = '$T_x$ fit, $\sigma_x$[mm]  vs. Time [ms]\n $T_x$ = %.2f' % T_x  + '$\pm %.2f[uK]$' % std_x
            titlestr_y = '$T_y$ fit, $\sigma_y$[mm]  vs. Time [ms]\n $T_y$ = %.2f' % T_y  + '$\pm %.2f[uK]$' % std_y
            resstr_x = '$T_x$ = %.2f' % T_y + '$\pm %.2f[uK]$' % std_y
            resstr_y = '$T_x$ = %.2f' % T_y + '$\pm %.2f[uK]$' % std_y
            self.plotDataAndFit(time_vector, sigma_x_vector, fitFunc=tempFromSigmaFunc, fitParams=x_temp_popt,
                                # title=f'X Temperature fit, Sigma_x[mm]  vs. Time [ms]\n %T_x% = {T_x} [uK]',
                                # ylabel=f'Sigma_x [mm]', saveFilePath=os.path.join(extraFilesPath, 'X_temp_fit.png'), show=False)
                                title=titlestr_x, props_str=resstr_x,
                                ylabel=f'$\sigma_x$ [mm]', saveFilePath=os.path.join(extraFilesPath, 'X_temp_fit.png'), show=False)
            self.plotDataAndFit(time_vector, sigma_y_vector, fitFunc=tempFromSigmaFunc, fitParams=y_temp_popt,
                                # title=f'Y Temperature fit, Sigma_y[mm]  vs. Time [ms]\n %T_y% = {T_y} [uK]',
                                # ylabel=f'Sigma_y [mm]', saveFilePath=os.path.join(extraFilesPath, 'Y_temp_fit.png'), show=False)
                                title=titlestr_y, props_str=resstr_y,
                                ylabel=f'$\sigma_y$ [mm]', saveFilePath=os.path.join(extraFilesPath, 'Y_temp_fit.png'), show=False)
        return T_x, T_y

    def perform_fit(self, path, fit_for_alpha=False):
        """
        TODO: <TBD>
        """
        # Set the folders we need
        extra_files = os.path.join(path, 'extra_files')
        background_image_path = os.path.join(extra_files, 'background.bmp')

        # Gaussian Fit to results (for each photo)
        gaussianFitResult = self.gaussianFitAllPicturesInPath(path, backgroundPath=background_image_path, saveFitsPath=extra_files, imgBounds=self.imgBounds)

        if len(gaussianFitResult) == 0:
            self.warn('No Gaussian fit on images, not trying any fits. Skipping.')
            return

        # ---- Take Gaussian fit results and get temperature (x,y) and launch speed -----------
        v_launch, alpha, v_launch_popt, v_launch_cov = self.fitVy_0FromGaussianFitResults(gaussianFitResult=gaussianFitResult, extraFilesPath=extra_files, plotResults=True)
        time_vector = np.array([res[0] for res in gaussianFitResult])
        y_position_vector = np.array([res[-1][2] for res in gaussianFitResult])

        # ------ Use v-launch fit to go from px to mm ------
        gaussian_amplitude_vector = alpha * np.array([res[-1][0] for res in gaussianFitResult])
        y_position_vector = alpha * y_position_vector
        x_position_vector = alpha * np.array([res[-1][1] for res in gaussianFitResult])
        sigma_x_vector = alpha * np.array([res[-1][3] for res in gaussianFitResult])
        sigma_y_vector = alpha * np.array([res[-1][4] for res in gaussianFitResult])

        # -------- Fit for x and y temperatures ----
        T_x, T_y = self.fitTemperaturesFromGaussianFitResults(gaussianFitResult, alpha, extra_files, plotResults=True)

        d = {
            'V_y Launch': v_launch,
            'T_x': T_x,
            'T_y': T_y,
            'alpha': alpha
        }
        print(d)
        try:
            with open(os.path.join(extra_files, 'results.json'), 'w') as file:
                json.dump(d, file, indent=4)
        except Exception as err:
            print(f'failed to save results json: {err}')

        plt.close('all')
        return d

    def nothing(self, val):
        pass

    def locate_funnel_in_image(self, image):
        pass

    def locate_chip_in_image(self, image):
        """
        We scan with a moving window from top of image downwards, looking for density change
        """
        chip_y = self.scan_for_density_change(image=image, direction='vertical', offset=0,
                                     frame_width=1300, frame_height=7,
                                     threshold_min=100, threshold_max=130)
        return chip_y

    def locate_fork_in_image(self, image):

        # TODO: Need to complete the work on this. Not working properly yet.
        # TODO: Complete finding of right_bar_2

        left_bar_1 = self.scan_for_density_change(image=image, direction='horizontal', offset=0,
                                     frame_width=7, frame_height=1200,
                                     threshold_min=100, threshold_max=130)

        left_bar_2 = self.scan_for_density_change(image=image, direction='horizontal', offset=left_bar_1,
                                     frame_width=7, frame_height=1200,
                                     threshold_min=100, threshold_max=130)

        right_bar_1 = self.scan_for_density_change(image=image, direction='horizontal', offset=left_bar_2+50,
                                     frame_width=7, frame_height=1200,
                                     threshold_min=100, threshold_max=130)

        right_bar_2 = 99

        print(f'left_bar_start={left_bar_1}, left_bar_end={left_bar_2}')
        print(f'right_bar_start={right_bar_1}, right_bar_end={right_bar_2}')

        pass

    def density_check(self, folder):

        files = Utils.get_files_in_path(folder, opt_in_filter='.bmp', return_full_path=True)
        file_name = files[4]  # Taking 2nd image, but can be any of them

        win_name = 'image'
        cv2.namedWindow(win_name, cv2.WINDOW_NORMAL)  # cv2.WND_PROP_FULLSCREEN)

        #loaded_img = cv2.imread(file_name, cv2.IMREAD_ANYCOLOR) # IMREAD_GRAYSCALE
        #loaded_img = cv2.imread(file_name, cv2.IMREAD_GRAYSCALE)
        loaded_img = cv2.imread(file_name, cv2.IMREAD_UNCHANGED)

        loaded_img = cv2.cvtColor(loaded_img, cv2.COLOR_BGR2RGB)

        line_thickness = 2
        line_color = (250, 150, 50)

        image_height = loaded_img.shape[0]
        image_width = loaded_img.shape[1]

        cv2.createTrackbar('contrast', 'image', 0, 500, self.nothing)
        cv2.createTrackbar('brightness', 'image', 0, 100, self.nothing)


        cv2.createTrackbar('width', 'image', 0, image_height, self.nothing)
        cv2.createTrackbar('height', 'image', 0, image_width, self.nothing)
        cv2.setTrackbarPos('width', 'image', 10)
        cv2.setTrackbarPos('height', 'image', 150)

        cv2.createTrackbar('min', 'image', 0, 255, self.nothing)
        cv2.createTrackbar('max', 'image', 0, 255, self.nothing)
        cv2.setTrackbarPos('min', 'image', 100)
        cv2.setTrackbarPos('max', 'image', 130)

        cv2.createTrackbar('x', 'image', 0, 2000, self.nothing)
        cv2.createTrackbar('y', 'image', 0, 2000, self.nothing)
        cv2.setWindowProperty(win_name, cv2.WND_PROP_FULLSCREEN, cv2.WINDOW_NORMAL)

        while True:

            width = cv2.getTrackbarPos('width', 'image')
            height = cv2.getTrackbarPos('height', 'image')

            min = cv2.getTrackbarPos('min', 'image')
            max = cv2.getTrackbarPos('max', 'image')

            x = cv2.getTrackbarPos('x', 'image')
            y = cv2.getTrackbarPos('y', 'image')

            brightness = cv2.getTrackbarPos('brightness', 'image')
            contrast = cv2.getTrackbarPos('contrast', 'image')
            contrast = 1.0 + (contrast / 100 * 2)

            img = loaded_img.copy()

            if True:
                # alpha = 3.9  # Contrast control (1.0-3.0)
                alpha = contrast
                beta = brightness  # 0  # Brightness control (0-100)
                img = cv2.convertScaleAbs(img, alpha=alpha, beta=beta)

            img = img[y:y + height, x:x + width]

            img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
            ret, thresh = cv2.threshold(img, min, max, cv2.THRESH_BINARY)  # cv2.THRESH_BINARY
            density = cv2.countNonZero(thresh)

            display_img = loaded_img.copy()

            # Change color of pixels in threshold window
            for x0 in range(0, width-1):
                for y0 in range(0, height-1):
                    if thresh[y0, x0] >= max:
                        p = display_img[y0+y, x0+x]
                        p[0] = 0
                        p[1] = 0
                        display_img[y0+y, x0+x] = p
                        pass

            display_img = cv2.rectangle(display_img, (x,y), (x+width, y+height), line_color, line_thickness)

            cv2.putText(display_img, str(density), (150,150), cv2.FONT_HERSHEY_SIMPLEX, 2, (180, 180, 180), 3, 1)

            cv2.imshow(win_name, display_img)
            k = cv2.waitKey(10) & 0xFF
            if k == 27:
                break

        cv2.destroyAllWindows()

        pass

    def compare_clouds_DEP(self):

        win_name = "Success!"

        folder1 = r'C:\temp\playback_data\Temperature\20240314_121527'
        #folder2 = r'C:\temp\playback_data\Temperature\20240314_122553'
        folder2 = r'C:\temp\playback_data\Temperature\20240314_143258'

        files1 = Utils.get_files_in_path(folder1, opt_in_filter='.bmp', return_full_path=True)
        files2 = Utils.get_files_in_path(folder2, opt_in_filter='.bmp', return_full_path=True)
        zipped = zip(files1, files2)

        height = 600  # 300
        width = 500
        line_thickness = 3
        #line_color = (255, 255, 255)
        line_color = (250, 150, 50)
        x = 900
        y = 300  # 550

        for tup in zipped:
            # Load images from folder 1 and folder 2
            image1 = cv2.imread(tup[0], cv2.IMREAD_ANYCOLOR)  #cv2.IMREAD_COLOR)
            image2 = cv2.imread(tup[1], cv2.IMREAD_ANYCOLOR)

            img = cv2.subtract(image1, image2)

            img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)

            alpha = 3.9  # Contrast control (1.0-3.0)
            beta = 0  # Brightness control (0-100)
            img = cv2.convertScaleAbs(img, alpha=alpha, beta=beta)

            top_left_point = (x, y)
            bottom_right = (x+width, y+height)
            img = cv2.rectangle(img, top_left_point, bottom_right, line_color, line_thickness)

            #img = img[y:y + height, x:x + width]

            # xxx = cv2.GaussianBlur(adjusted, (7,7), 0)
            # zzz = cv2.Laplacian(xxx, -1, ksize=3, scale=1, delta=0, borderType=cv2.BORDER_DEFAULT)

            win_name = 'image'  # tup[0]
            cv2.namedWindow(win_name, cv2.WND_PROP_FULLSCREEN)
            cv2.createTrackbar('min', 'image', 0, 255, self.nothing)
            cv2.createTrackbar('max', 'image', 0, 255, self.nothing)
            cv2.setWindowProperty(win_name, cv2.WND_PROP_FULLSCREEN, cv2.WINDOW_NORMAL)

            cv2.imshow(win_name, img)

            cv2.waitKey(0)

            #cv2.destroyAllWindows()
            cv2.destroyWindow(win_name)
        pass

    def measure_temperature(self, pre_pulse_durations=np.arange(1, 9), create_video=False, perform_fit=True):
        """
        Runs a loop over array of PrePulse Duration values.
        For each setting, takes an image and saves in a folder.
        Background image is also taken for reference at the end when there is no MOT
        """

        # Connect to camera. Notify if we fail to connect.
        self.connect_camera()
        if self.is_camera_disconnected():
            self.warn('Cannot connect to camera - maybe app is open?')
            return

        # Create the folders for the images to be saved
        self.bd_results.create_folders()
        path = self.bd_results.get_folder_path('root')
        extra_files = self.bd_results.get_folder_path('extra_files')

        # Iterate over PrePulse Duration parameters and take photos
        try:
            for ppd in pre_pulse_durations:
                # Update prepulse duration
                self.updateValue("PrePulse_duration", float(ppd), update_parameters=True)

                # Capture and average images
                image_name = 'PrePulse_duration={:04.1f}.bmp'.format(ppd)
                image_full_path = os.path.join(path, image_name)
                self.camera.saveAverageImage(image_full_path, NAvg=self.NAvg, NThrow=self.NThrow, NThrow_end=self.NThrow_end, RGB=False)
        except Exception as err:
            self.warn(f'Failed to get images - {err}')

        # Take background picture (after such a long delay, we should have no visible cloud)
        self.updateValue("PrePulse_duration", float(50), update_parameters=True)

        background_image_path = os.path.join(extra_files, 'background.bmp')
        self.camera.saveAverageImage(background_image_path, NAvg=self.NAvg, NThrow=self.NThrow, RGB=False)

        # Update back the pre_pulse_duration to the Config original value
        ppd_at_start = float(self.Exp_Values['PrePulse_duration'])
        self.info(f'Updating PrePulse Duration back to {ppd_at_start}')
        self.updateValue("PrePulse_duration", float(ppd_at_start), update_parameters=True)

        # If required, create video from images
        if create_video:
            self.create_video_from_path(path, save_file_path=os.path.join(extra_files, 'video.avi'))

        if perform_fit:
            self.perform_fit(path)

        pass

    def optimizePGC(self, PrePulseDurations=np.arange(0.5,5,0.5),PGC_final_freq_params = np.linspace(97e6, 90e6, 7),PGC_final_amp_params = np.linspace(0.3, 0.05, 7)):
        """
        This is used to run with different parameters. They are defined within the function

        :param PrePulseDurations:
        :return: <TBD>
        """
        path, extraFilesPath = self.createPathForTemperatureMeasurement()

        # ---- Take background picture ----
        self.updateValue("PrePulse_duration", float(50),
                         update_parameters=True)  # Presumably, after 20ms there's no visible cloud.
        self.camera.saveAverageImage(os.path.join(extraFilesPath, 'background.bmp'), NAvg=self.NAvg, NThrow=self.NThrow,NThrow_end=self.NThrow_end, RGB=False)

        # ---- Define parameters space ------
        xs_key = "PGC_final_freq"
        xs = PGC_final_freq_params  # e.g., PGCFinalFreq
        ys_key = "PGC_final_amp"
        ys = PGC_final_amp_params  # e.g., PGCFinalAmp

        # --- Begin measurement -----
        for ppd in PrePulseDurations:
            self.updateValue("PrePulse_duration", float(ppd))
            for x in xs:
                self.updateValue(xs_key, float(x))
                for y in ys:
                    self.updateValue(ys_key, float(y))
                    self.update_parameters()
                    imgName = 'PrePulse=%s;%s=%s;%s=%s.bmp' %(str(ppd), str(xs_key), str(x), str(ys_key), str(y))
                    self.camera.saveAverageImage(os.path.join(path, imgName), NAvg=self.NAvg, NThrow=self.NThrow, RGB=False)

        print('Finished  optimizePGC. running analyzePGCOptimization...')
        return self.analyzePGCOptimization(path=path)

    def analyzePGCOptimization(self, path, backgroundPath=None):
        # ---- fit gaussians to all photos -----
        gaussianFitsAndParams = self.gaussianFitAllPicturesInPath(path, backgroundPath, saveFitsPath=None, imgBounds=self.imgBounds)
        np.save(os.path.join(path, 'gasuusainFitResultsresults.npy'), gaussianFitsAndParams, allow_pickle=True)
        # ----- extract the parameters we ran over (say, PGC_final_freq, PGC_final_amp) as xs and ys, and the prePulse_duration (ppds) for each photos ----
        xs = np.array([item[1] for item in gaussianFitsAndParams])
        ys = np.array([item[2] for item in gaussianFitsAndParams])
        prePulseD = np.array([item[0] for item in gaussianFitsAndParams])
        # -----
        gaussian_amp = np.array([item[-1][0] for item in gaussianFitsAndParams])
        sigmas_x = np.array([item[-1][3] for item in gaussianFitsAndParams])
        sigmas_y = np.array([item[-1][4] for item in gaussianFitsAndParams])

        uniqueXs = np.unique(xs)
        uniqueYs = np.unique(ys)
        uniqueTemperatures =[]
        xs_for_plot, ys_for_plot = [],[]
        for x in uniqueXs:
            for y in uniqueYs:
                indices = np.argwhere((xs == x) & (ys == y))
                ppds, s_xs_perParam, s_ys_perParam = [],[],[]
                for i in indices:
                    i= i[0] # since np.argwhere returns an array with the index inside
                    s_xs_perParam.append(sigmas_x[i])
                    s_ys_perParam.append(sigmas_y[i])
                    ppds.append(prePulseD[i])
                # ---- sort these lists according to the pre-pulse duration
                s_xs_perParam = [s_xs_perParam for _, s_xs_perParam in sorted(zip(ppds, s_xs_perParam))]
                s_ys_perParam = [s_ys_perParam for _, s_ys_perParam in sorted(zip(ppds, s_ys_perParam))]
                ppds = sorted(ppds)
                # ----- fit for temperature!
                T_x, T_y = self.fitTemperaturesFromSigmasResults(time_vector=ppds, sigma_x_vector=s_xs_perParam, sigma_y_vector=s_ys_perParam, alpha=self.mm_to_pxl, plotResults=False)
                temp = np.sqrt(T_x**2 + T_y**2)
                # ---- Arrange data for plotting
                xs_for_plot.append(x)
                ys_for_plot.append(y)
                uniqueTemperatures.append(temp)
                print('Results for x = {}, y = {} are: {}K'.format(str(x), str(y), str(temp)))

        self.plot3D(xs_for_plot, ys_for_plot, uniqueTemperatures)
        return (xs_for_plot, ys_for_plot, uniqueTemperatures)

    # TODO: we may want to move this to Utils - so everyone can enjoy :-)
    def create_video_from_path(self, path, save_file_path=None):
        Utils.create_video_from_path(path, save_file_path)
        pass

    def open_windows_explorer(self):
        """ Opens Windows Explorer at the folder of the last measurement """

        # Get experiment folder
        experiment_root_path = self.bd_results.get_custom_root('experiment_root')

        # Get all folders in experiment folder and pick the most recent one
        only_folders = [f for f in os.listdir(experiment_root_path) if not os.path.isfile(os.path.join(experiment_root_path, f))]
        most_recent_folder = os.path.join(experiment_root_path, only_folders[-1], 'extra_files')
        Utils.open_windows_explorer(most_recent_folder)
        self.info(f'Opened Windows Explorer at {most_recent_folder}')
        pass

    def gaussian_fit_all_pictures_in_path(self, path, backgroundPath=None, saveFitsPath=None, imgBounds=None):
        if backgroundPath is None:
            backgroundPath = os.path.join(path, 'extra_files', 'background.bmp')

        if imgBounds is None:
            imgBounds = self.imgBounds

        # Get all files in folder (except the 'background' file)
        files = Utils.get_files_in_path(path, opt_out_filter='background', return_full_path=True)

        res = []
        for full_file_name in files:
            try:
                f = Path(full_file_name).stem.replace(',', ';')
                fileRes = []
                ## -- get x and y values from file name ----
                keysAndValues = f.split(';')
                keys = []
                for key in keysAndValues:
                    value = key[key.find('=') + 1:]
                    only_key = key[:key.find('=')]
                    fileRes.append(float(value))
                    keys.append(only_key)
                ## -- get z value from gaussian fit ----
                gaussianFit = self.GaussianFit(full_file_name, background_file=backgroundPath, saveFitsPath=saveFitsPath, imgBounds=imgBounds)
                if gaussianFit is None:
                    continue  # if fit returned None, meaning the fit failed (sigma is out of self.sigma_bounds), discard this results and continue to the next fit
                fileRes.append(keys)
                fileRes.append(gaussianFit)
                res.append(fileRes)
            except Exception as e:
                print(e)

        pass

    def gaussianFitAllPicturesInPath(self, path, backgroundPath=None, saveFitsPath=None, imgBounds=None):
        if backgroundPath is None:
            backgroundPath = os.path.join(path, 'extra_files', 'background.bmp')

        if imgBounds is None:
            imgBounds = self.imgBounds

        files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
        res = []
        for f in files:
            try:
                fullFileName = os.path.join(path, f)
                fileRes = []
                f = f.replace('.jpeg', '').replace('.bmp', '').replace(',', ';')
                if f.find('background'):
                    pass
                ## -- get x and y values from file name ----
                keysAndValues = f.split(';')
                keys = []
                for key in keysAndValues:
                    value = key[key.find('=') + 1:]
                    only_key = key[:key.find('=')]
                    fileRes.append(float(value))
                    keys.append(only_key)
                ## -- get z value from gaussian fit ----
                gaussianFit = self.GaussianFit(fullFileName, background_file=backgroundPath, saveFitsPath=saveFitsPath, imgBounds=imgBounds)
                if gaussianFit is None:
                    continue  # if fit returned None, meaning the fit failed (sigma is out of self.sigma_bounds), discard this results and continue to the next fit
                fileRes.append(keys)
                fileRes.append(gaussianFit)
                res.append(fileRes)
            except Exception as e:
                print(e)
        return (res)


    def GaussianFit(self, file_name_for_fit, background_file, saveFitsPath=None, imgBounds=None, X_PIXEL_LEN=1544, Y_PIXEL_LEN=2064, CROP_IMG_SIZE=180 ,PLOT_IMG=False, PLOT_SLICE=False, SHOW_CROP=False):
        """
        Fit a gaussian to subtracted images
         INPUT:
           Y_PIXEL_LEN - number of vertical pixels
           X_PIXEL_LEN - number of horizontal pixels
           the name of the file for fit with png termination
           CROP_IMG_SIZE - the size of image in each dimension after cropping is 2*CROP_IMG_SIZE
        """
        fileName = Path(file_name_for_fit).stem
        backgroundImg = cv2.imread(background_file, 0)
        ImgToFit = cv2.imread(file_name_for_fit, 0)
        ImgToFit = cv2.subtract(ImgToFit, backgroundImg)

        if not imgBounds:
            h, w = ImgToFit.shape
            imgBounds = (0, 0, w, h)  # that is, don't crop image

        # Show the img bounds on top of the image taken with atoms
        if SHOW_CROP:
            line_thickness = 3
            line_color = (255, 255, 255)  # White
            image_with_bounds = ImgToFit.copy()
            top_left_point = (imgBounds[0], imgBounds[1])
            bottom_right_point = (imgBounds[2], imgBounds[3])
            image_with_bounds = cv2.rectangle(image_with_bounds, top_left_point, bottom_right_point, line_color, line_thickness)
            cv2.imshow("Fit Crop Bounds", image_with_bounds)
            cv2.waitKey(0)

        if imgBounds: ImgToFit = ImgToFit[imgBounds[1]:imgBounds[3], imgBounds[0]:imgBounds[2]]  #
        img_max_index = [np.argmax(np.sum(ImgToFit, axis=0)), np.argmax(np.sum(ImgToFit, axis=1))]

        # img_max_index[1] = np.argmax(np.sum(ImgToFit[10:][img_max_index[0]-CROP_IMG_SIZE:img_max_index[0]+CROP_IMG_SIZE], axis=1))
        print(img_max_index)
        # Parameters
        X_UPPER_BOUND = int(img_max_index[0] + CROP_IMG_SIZE)
        X_LOWER_BOUND = int(img_max_index[0] - CROP_IMG_SIZE)
        if X_LOWER_BOUND < 0: X_LOWER_BOUND = 0
        if X_UPPER_BOUND > X_PIXEL_LEN: X_UPPER_BOUND = X_PIXEL_LEN
        EFFECTIVE_X_PIXEL_LEN = X_UPPER_BOUND - X_LOWER_BOUND

        Y_UPPER_BOUND = int(img_max_index[1] + CROP_IMG_SIZE)
        Y_LOWER_BOUND = int(img_max_index[1] - CROP_IMG_SIZE)
        if Y_LOWER_BOUND < 0: Y_LOWER_BOUND = 0
        if Y_UPPER_BOUND > Y_PIXEL_LEN: Y_UPPER_BOUND = Y_PIXEL_LEN
        EFFECTIVE_Y_PIXEL_LEN = Y_UPPER_BOUND - Y_LOWER_BOUND

        # Create x and y indices
        x = np.linspace(0, EFFECTIVE_Y_PIXEL_LEN - 1, EFFECTIVE_Y_PIXEL_LEN)
        y = np.linspace(0, EFFECTIVE_X_PIXEL_LEN - 1, EFFECTIVE_X_PIXEL_LEN)
        x, y = np.meshgrid(x, y)
        # crop an effective image
        EffectiveImg = ImgToFit[Y_LOWER_BOUND:Y_UPPER_BOUND, X_LOWER_BOUND:X_UPPER_BOUND]
        # plt.imshow(EffectiveImg)
        # plt.show()
        data_noisy = EffectiveImg.ravel()

        # fit the data
        img_max_index = [np.argmax(np.sum(EffectiveImg, axis=0)), np.argmax(np.sum(EffectiveImg, axis=1))]

        #img_max_index = [np.argmax(np.sum(EffectiveImg, axis=1)), np.argmax(np.sum(EffectiveImg, axis=0))]
        initial_guess = (ImgToFit[img_max_index[1]][img_max_index[0]], img_max_index[1], img_max_index[0], EFFECTIVE_X_PIXEL_LEN/10, EFFECTIVE_Y_PIXEL_LEN/10, 0, 10)
        fitBounds =[0, (255, EFFECTIVE_X_PIXEL_LEN, EFFECTIVE_Y_PIXEL_LEN, EFFECTIVE_X_PIXEL_LEN / 2, EFFECTIVE_Y_PIXEL_LEN / 2, 2* np.pi, 255)]
        # print(initial_guess)
        popt, pcov = opt.curve_fit(self.twoD_Gaussian, (x, y), data_noisy, p0=initial_guess, bounds=fitBounds)

        # --- Check sigmas are in-bound ----
        sigma = [popt[3], popt[4]]
        sb = self.sigma_bounds
        if sb and (not sb[0] < sigma[0] < sb[1] or not sb[0] < sigma[1] < sb[1]):  # if either of the sigmas is out of bounds...
            print('In gaussian git, sigma out of bounds. Discarding this result')
            return None

        # ---- plot the results ----
        data_fitted = self.twoD_Gaussian((x, y), *popt)

        if PLOT_IMG or saveFitsPath:
            fig, ax = plt.subplots(1, 1)
            plt.title(fileName)
            plt.text(0.88, 0.95, '\u03C3_x =' + '%.2f' % sigma[0] + '\n' + '\u03C3_y = ' + '%.2f' % sigma[1], color='white',
                 fontsize=16, style='italic', weight='bold', horizontalalignment='center', verticalalignment='center',
                 transform=ax.transAxes, bbox=dict(facecolor='gray', alpha=0.5))
            ax.imshow(data_noisy.reshape(EFFECTIVE_X_PIXEL_LEN, EFFECTIVE_Y_PIXEL_LEN), cmap=plt.cm.jet, origin='lower',
                  extent=(x.min(), x.max(), y.min(), y.max()))
            ax.contour(x, y, data_fitted.reshape(EFFECTIVE_X_PIXEL_LEN, EFFECTIVE_Y_PIXEL_LEN), 8, colors='w')
            plt.savefig(os.path.join(saveFitsPath, fileName + '.png'))
            if PLOT_IMG: plt.show()

        if PLOT_SLICE:
            # plot slice
            mtl.figure()
            data_noisy_mat = data_noisy.reshape(EFFECTIVE_X_PIXEL_LEN, EFFECTIVE_Y_PIXEL_LEN)
            Slice = data_noisy_mat[int(EFFECTIVE_X_PIXEL_LEN / 2) - 10]
            mtl.plot(Slice)
            mtl.plot(data_fitted.reshape(EFFECTIVE_X_PIXEL_LEN, EFFECTIVE_Y_PIXEL_LEN)[int(EFFECTIVE_X_PIXEL_LEN / 2) - 10])
            plt.show()

        # return original x-y
        popt[1] = popt[1] + X_LOWER_BOUND + imgBounds[0]
        popt[2] = popt[2] + Y_LOWER_BOUND + imgBounds[1]

        return popt

    # define model function and pass independent variables x and y as a list
    def twoD_Gaussian(self, x_y, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
        x, y = x_y
        xo = float(xo)
        yo = float(yo)
        a = (np.cos(theta) ** 2) / (2 * sigma_x ** 2) + (np.sin(theta) ** 2) / (2 * sigma_y ** 2)
        b = -(np.sin(2 * theta)) / (4 * sigma_x ** 2) + (np.sin(2 * theta)) / (4 * sigma_y ** 2)
        c = (np.sin(theta) ** 2) / (2 * sigma_x ** 2) + (np.cos(theta) ** 2) / (2 * sigma_y ** 2)
        # g = offset + amplitude * np.exp(- (a * ((x - xo) ** 2) + 2 * b * (x - xo) * (y - yo)
        #                                    + c * ((y - yo) ** 2)))
        g = offset + amplitude * np.exp(- ((x - xo) / (np.sqrt(2) * sigma_x))**2 - ((y - yo) /(np.sqrt(2)*sigma_y))**2)
        return g.ravel()
    def plot3D(self, xs, ys, zs):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # Plot the surface
        ax.plot_trisurf(xs, ys, zs)

    def plotDataAndFit(self, x_data, y_data, fitFunc, fitParams, title='', props_str='', ylabel='', xlabel='Time [ms]',
                       saveFilePath=None, show=False):
        x_vector_for_fit = np.linspace(min(x_data), max(x_data), 100)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        fig, ax = plt.subplots()
        ax.plot(x_data, y_data, 'x', label='Data')
        ax.plot(x_vector_for_fit, fitFunc(x_vector_for_fit, *fitParams), '-', label='Fit')
        ax.set(xlabel=xlabel, ylabel=ylabel)
        ax.set_title(title)
        ax.legend()
        ax.text(0.5, 0.95, props_str, transform=ax.transAxes, fontsize=18, horizontalalignment='center',  verticalalignment='top')
        if saveFilePath:
            plt.savefig(saveFilePath)
        if show:
            plt.show()

    def _scan_vertically_for_density_change(self, image, offset, frame_width, frame_height, threshold_min, threshold_max):
        image_shape = image.shape
        image_height = image_shape[0]
        image_width = image_shape[1]

        t = None
        densities = np.array([])
        center_x = int(image_width / 2)
        half_width = int(frame_width / 2)
        for y in range(0, image_height-frame_height):
            crop = image[y:y + frame_height, (center_x-half_width):(center_x+half_width)]

            ret, thresh = cv2.threshold(crop, threshold_min, threshold_max, cv2.THRESH_BINARY)  # cv2.THRESH_BINARY
            density = cv2.countNonZero(thresh)
            densities = np.append(densities, density)

            print(f'y={y} -> density={density}')

            if density < (np.average(densities) / 5):
                t = y
                break

        # If we couldn't find any density change in image, we return None
        if not t:
            return None

        # We return the y of the change
        y_change = t - int(frame_height/2)
        return y_change

    def _scan_horizontally_for_density_change(self, image, offset, frame_width, frame_height, threshold_min, threshold_max):

        image_shape = image.shape
        image_height = image_shape[0]
        image_width = image_shape[1]

        t = None
        densities = np.array([])
        center_y = int(image_height / 2)
        half_height = int(frame_height / 2)
        for x in range(offset, image_width-frame_width):
            crop = image[(center_y-half_height):(center_y+half_height), x:x + frame_width]

            ret, thresh = cv2.threshold(crop, threshold_min, threshold_max, cv2.THRESH_BINARY)  # cv2.THRESH_BINARY
            density = cv2.countNonZero(thresh)
            densities = np.append(densities, density)

            print(f'x={x} -> density={density}')

            if density != 0 and density < (np.average(densities) / 5):  # Moved from Bright to Dark
                t = x
                break
            elif density > 4000 and density/5 > np.average(densities):  # Moved from Dark to Bright
                t = x
                break

        # If we couldn't find any density change in image, we return None
        if not t:
            return None

        # We return the y of the change
        x_change = t - int(frame_width/2)
        return x_change

    def scan_for_density_change(self, image, direction, offset, frame_width, frame_height, threshold_min, threshold_max):

        # Increase Contrast
        alpha = 3.9  # Contrast control (1.0-3.0)
        beta = 0  # Brightness control (0-100)
        image = cv2.convertScaleAbs(image, alpha=alpha, beta=beta)
        image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

        if direction=='vertical':
            coor = self._scan_vertically_for_density_change(image, offset, frame_width, frame_height, threshold_min, threshold_max)
        elif direction=='horizontal':
            coor = self._scan_horizontally_for_density_change(image, offset,  frame_width, frame_height, threshold_min, threshold_max)
        else:
            print(f'Direction {direction} not supported')
            return None

        return coor

    def update_value_in_opx(self, prepulse_duration):
        """
        This method is invoked directly from the Menu (and not referenced here in the code)
        :param prepulse_duration: The duration before we take the image
        """
        self.updateValue("PrePulse_duration", float(prepulse_duration), update_parameters=True)

if __name__ == "__main__":

    # Initiate the experiment
    # Change to ExperimentMode.OFFLINE if you wish to run outside the lab
    experiment = CoolingSequenceOptimizer()

    # Display menu to get action
    settings = Utils.load_json_from_file(r'./settings.json')
    selection = BDMenu(caller=experiment, menu_file=None, menu_json=settings['menus']).display()

    # TODO: insert this into the cycles/runs schema
    #experiment.measureTemperature()
    #experiment.optimizePGC()

    pass