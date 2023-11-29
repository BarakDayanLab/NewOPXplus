from Experiments.BaseExperiment.BaseExperiment import BaseExperiment

import numpy as np
from PIL import Image
import time
import json
import os
from pathlib import Path
import cv2, glob
import scipy.optimize as opt
import matplotlib.pyplot as mtl
import pylab as plt
from Utilities.BDMenu import BDMenu
from mpl_toolkits.mplot3d import Axes3D

_Kb = 1.38e-23  #[J/K]
_m = 1.443e-25  #[Kg] of Rb87


class CoolingSequenceOptimizer(BaseExperiment):

    def __init__(self):
        # Invoking BaseClass constructor. It will initiate OPX, QuadRF, BDLogger, Camera, BDResults, KeyEvents etc.
        super().__init__()
        pass

    def initialize_experiment_variables(self):
        # self.SPRINT_Exp_switch(False) # To enable trigger to camera in OPX_control_with_QuadRF_Sprint_Exp
        self.camera = None
        self.NAvg = 1  # Number of photos captured to create an average image
        self.NThrow = 3  # Number of throws we're "skipping" until capturing an image
        self.imgBounds = (580, 200, 1600, 1450)  # bounds to crop out of the taken pictures
        self.mm_to_pxl = 8.5/(830-56)  # measured using ruler in focus 13/11/2022
        self.sigma_bounds = (15, 100)  # This bounds sigma (x & y) of the Gaussian sigma. If value is out of bounds, fit is considered bad and not used in temp-fit

    def connect_disconnect_camera(self):
        if self.camera:
            self.disconnect_camera()
        else:
            self.connect_camera()

    def createPathForTemperatureMeasurement(self, path=None, saveConfig=False):
        extraFilesPath = ''
        if path:
            extraFilesPath = os.path.join(path, 'extra_files')
        else:
            # --- create path and save config -----
            timeStamp = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            path = f'U:\\Lab_2021-2022\\Experiment_results\\Temperature\\{timeStamp}'
            extraFilesPath = os.path.join(path, 'extra_files')
            if not os.path.exists(path):
                os.makedirs(path)
            if not os.path.exists(extraFilesPath):
                os.makedirs(extraFilesPath)
            if saveConfig:
                self.saveConfigTable(path=extraFilesPath)

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

        if plotResults:
            self.plotDataAndFit(time_vector, y_position_vector, fitFunc=fitFunc, fitParams=v_launch_popt,
                                title=f'Y position fit \n V_launch = {v_launch}; alpha = {alpha}', ylabel='Y_center [px]',
                                saveFilePath=extraFilesPath + 'Y_position.png', show=False)
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

        if plotResults:
            self.plotDataAndFit(time_vector, sigma_x_vector, fitFunc=tempFromSigmaFunc, fitParams=x_temp_popt,
                                title=f'X Temperature fit, Sigma_x[mm]  vs. Time [ms]\n T_x = {T_x} [uK]',
                                ylabel=f'Sigma_x [mm]', saveFilePath=extraFilesPath + 'X_temp_fit.png', show=False)
            self.plotDataAndFit(time_vector, sigma_y_vector, fitFunc=tempFromSigmaFunc, fitParams=y_temp_popt,
                                title=f'Y Temperature fit, Sigma_y[mm]  vs. Time [ms]\n T_y = {T_y} [uK]',
                                ylabel=f'Sigma_y [mm]', saveFilePath=extraFilesPath + 'Y_temp_fit.png', show=False)
        return (T_x, T_y)

    def measureTemperature(self, path=None, PrePulseDurations=np.arange(1,9), createVideo=False, fit_for_alpha=False):
        extraFilesPath = ''
        if path:
            # if received path, assume there's no need to take new pictures and analyze these photos.
            extraFilesPath = os.path.join(path, 'extra_files')
        else:
            self.connect_camera()
            if self.is_camera_disconnected():
                self.warn('Cannot connect to camera - maybe app is open?')
                return

            # TODO: TEMP - REMOVE
            path = r'C:\temp\dror_debug'

            path, extraFilesPath = self.createPathForTemperatureMeasurement(path, saveConfig=True)

            # ---- Take photos  -----
            try:
                for ppd in PrePulseDurations:
                    # ---- update prepulse duration ----
                    self.updateValue("PrePulse_duration", float(ppd), update_parameters=True)
                    # ---- capture and average pictures ------
                    image_name = 'PrePulse_duration={:04.1f}.bmp'.format(ppd)
                    image_full_path = os.path.join(path, image_name)
                    self.camera.saveAverageImage(image_full_path, NAvg=self.NAvg, NThrow=self.NThrow, RGB=False)
            except Exception as err:
                self.warn(f'Failed to get images - {err}')

            # ---- Take background picture ----
            self.updateValue("PrePulse_duration", float(50), update_parameters=True)  # Presumably, after 20ms there's no visible cloud.

            background_image_path = os.path.join(extraFilesPath, 'background.bmp')
            self.camera.saveAverageImage(background_image_path, NAvg=self.NAvg, NThrow=self.NThrow, RGB=False)

        # --- Gaussian Fit to results (to each photo) ----
        gaussianFitResult = self.gaussianFitAllPicturesInPath(path, backgroundPath=background_image_path, saveFitsPath=extraFilesPath, imgBounds=self.imgBounds)

        # ---- Take Gaussian fit results and get temperature (x,y) and launch speed -----------
        v_launch, alpha, v_launch_popt, v_launch_cov = self.fitVy_0FromGaussianFitResults(gaussianFitResult=gaussianFitResult, extraFilesPath=extraFilesPath, plotResults=True)
        time_vector = np.array([res[0] for res in gaussianFitResult])
        y_position_vector = np.array([res[-1][2] for res in gaussianFitResult])

        # ------ Use v-launch fit to go from px to mm ------
        gaussian_amplitude_vector = alpha * np.array([res[-1][0] for res in gaussianFitResult])
        y_position_vector = alpha * y_position_vector
        x_position_vector = alpha * np.array([res[-1][1] for res in gaussianFitResult])
        sigma_x_vector = alpha * np.array([res[-1][3] for res in gaussianFitResult])
        sigma_y_vector = alpha * np.array([res[-1][4] for res in gaussianFitResult])

        # -------- Fit for x and y temperatures ----
        T_x, T_y = self.fitTemperaturesFromGaussianFitResults(gaussianFitResult, alpha, extraFilesPath, plotResults=True)

        d = {
            'V_y Launch': v_launch,
            'T_x': T_x,
            'T_y': T_y,
            'alpha': alpha
        }
        print(d)
        try:
            with open(os.path.join(extraFilesPath, 'results.json'), 'w') as file:
                json.dump(d, file, indent=4)
        except Exception:
            pass

        if createVideo:
            self.createVideoFromPath(path, saveFilePath=os.path.join(extraFilesPath, 'video.avi'))

        plt.close('all')
        return d

    def optimizePGC(self, PrePulseDurations=np.arange(1,7)):
        """
        This is used to run with different parameters. They are defined within the function

        :param PrePulseDurations:
        :return: <TBD>
        """
        path, extraFilesPath = self.createPathForTemperatureMeasurement()

        # ---- Take background picture ----
        self.updateValue("PrePulse_duration", float(50),
                         update_parameters=True)  # Presumably, after 20ms there's no visible cloud.
        self.camera.saveAverageImage(os.path.join(extraFilesPath, 'background.bmp'), NAvg=self.NAvg, NThrow=self.NThrow, RGB=False)

        # ---- Define parameters space ------
        xs_key = "PGC_final_freq"
        xs = np.linspace(97e6, 90e6, 7)  # e.g., PGCFinalFreq
        ys_key = "PGC_final_amp"
        ys = np.linspace(0.3, 0.05, 7)  # e.g., PGCFinalAmp

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

    def createVideoFromPath(self, path, saveFilePath=None):
        # Taken from: https://theailearner.com/2018/10/15/creating-video-from-images-using-opencv-python/
        if not saveFilePath:
            saveFilePath = os.path.join(path, 'video.avi')
        img_array = []
        fileNames = glob.glob(os.path.join(path,'*.bmp'))
        fileNames.sort()
        for filename in fileNames:
            img = cv2.imread(filename)
            height, width, layers = img.shape
            size = (width, height)
            img_array.append(img)

        out = cv2.VideoWriter(filename=saveFilePath, fourcc=cv2.VideoWriter_fourcc(*'DIVX'), fps=3, frameSize=size)
        for i in range(len(img_array)):
            out.write(img_array[i])
        out.release()

    def gaussianFitAllPicturesInPath(self, path, backgroundPath=None, saveFitsPath=None, imgBounds=None):
        if backgroundPath is None:
            backgroundPath = os.path.join(path, 'extra_files', 'background.bmp')

        if imgBounds is None:
            imgBounds = self.imgBounds

        files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
        res = []
        for f in files:
            try:
                fullFileName = path + f
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


    '''fit a gaussian to subtracted images
    # INPUT:
    #   Y_PIXEL_LEN - number of vertical pixels
    #   X_PIXEL_LEN - number of horizontal pixels
    #   the name of the file for fit with png termination
    #   CROP_IMG_SIZE - the size of image in each dimension after cropping is 2*CROP_IMG_SIZE
    '''
    def GaussianFit(self, file_name_for_fit, background_file, saveFitsPath=None, imgBounds=None, X_PIXEL_LEN=1544, Y_PIXEL_LEN=2064, CROP_IMG_SIZE=180 ,PLOT_IMG=False, PLOT_SLICE=False):
        fileName = Path(file_name_for_fit).stem
        backgroundImg = cv2.imread(background_file, 0)
        ImgToFit = cv2.imread(file_name_for_fit, 0)
        ImgToFit = cv2.subtract(ImgToFit, backgroundImg)

        if not imgBounds:
            h, w = ImgToFit.shape
            imgBounds = (0, 0, w, h)  # that is, don't crop image

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
        popt, pcov = opt.curve_fit(self.twoD_Gaussian, (x, y), data_noisy, p0=initial_guess, bounds = fitBounds)

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
            ax.imshow(data_noisy.reshape(EFFECTIVE_X_PIXEL_LEN, EFFECTIVE_Y_PIXEL_LEN), cmap=plt.cm.jet, origin='bottom',
                  extent=(x.min(), x.max(), y.min(), y.max()))
            ax.contour(x, y, data_fitted.reshape(EFFECTIVE_X_PIXEL_LEN, EFFECTIVE_Y_PIXEL_LEN), 8, colors='w')
            plt.savefig(saveFitsPath + fileName + '.png')
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

    def plotDataAndFit(self, x_data, y_data, fitFunc,fitParams,title='', ylabel='', xlabel='Time [ms]', saveFilePath=None, show=False):
        x_vector_for_fit = np.linspace(min(x_data), max(x_data), 100)
        plt.figure()
        plt.plot(x_data,y_data,'x', label='Data')
        plt.plot(x_vector_for_fit, fitFunc(x_vector_for_fit, *fitParams), '-', label='Fit')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.legend()
        if saveFilePath:
            plt.savefig(saveFilePath)
        if show:
            plt.show()

    def update_value_in_opx(self, prepulse_duration):
        """
        This method is invoked directly from the Menu (and not referenced here in the code)
        :param prepulse_duration: The duration before we take the image
        """
        self.updateValue("PrePulse_duration", float(prepulse_duration), update_parameters=True)

#r = optimizePGC(c)
if __name__ == "__main__":

    # Initiate the experiment
    experiment = CoolingSequenceOptimizer()

    # Display menu to get action
    selection = BDMenu(experiment, r'./menu.json', None).display()

    # TODO: insert this into the cycles/runs schema
    #experiment.measureTemperature()
    #experiment.optimizePGC()