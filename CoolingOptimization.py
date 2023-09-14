from OPX_control_with_QuadRF_TransitPLUS_Exp import OPX
# from OPX_control_with_QuadRF_Sprint_Exp import OPX
import Config_with_SNSPDs_and_QuadRF as Config
# import Config_with_SNSPDs_and_QuadRF_Sprint as Config
import numpy as np
from PIL import Image, ImageDraw
import time, json, os
from pathlib import Path
import cv2, glob
import scipy.optimize as opt
from scipy import fftpack
import matplotlib.pyplot as mtl
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
from pynpunt import keyboard
import numpy as np


_Kb = 1.38e-23; #[J/K]
_m = 1.443e-25; #[Kg] of Rb87

# --------- Camera functionality ------------##
try:
    from UtilityResources import MvCameraController
    # from OPXcontrol.OPX_control_New_v1 import OPX
    from mvIMPACT import acquire
except:
    print('Warning! could not import MvCamera module; in MvCamera.py')

class CoolingSequenceOptimizer(OPX):
    DISTANCE_BETWEEN_FORKS_IN_INVENTOR = 20 # [mm]
    def __init__(self,  config, camera = None):

        super().__init__(config = config)
        # self.SPRINT_Exp_switch(False) # To enable trigger to camera in OPX_control_with_QuadRF_Sprint_Exp
        self.camera = camera
        self.NAvg = 1
        self.NThrow = 3
        self.imgBounds = (580,200,1600,1450) # bounds to crop out of the taken pictures
        self.mm_to_pxl = 8.5/(830-56) # measured using ruler in focus 13/11/2022
        self.sigma_bounds = (15,100) # This bounds sigma (x & y) of the Gaussian sigma. If value is out of bounds, fit is considered bad and not used in temp-fit

    def connectCamera(self):
        try:
            self.camera = MvCameraController.MvCameraController()
            return True
        except Exception as e:
            print('Could not connect to camera. Aborting...', e)
            return False

    def disconnectCamera(self):
        self.camera = None

    def mm_to_pxl_auto_calibration(self, folder_path, backgroundPath=None, print_all_lines=False, apply_lpf=True, show_im=True):
        """
        This function tries to automaticlly calibrate the mm_to_pxl attribute by using image processing to calculate the distance between the two metal bars of the 
        for and devide it by the known distance between them from the inventor file used to create it.
        resulting calibration is stored in self.mm_to_pxl
        :param folder_path: the folder path of the images taken by camera in the current measurement
        :param backgroundPath: a known backgroundPath of an image used to run the function manually
        :param print_all_lines: used for debugging why the algorithm failed and see all the lines identified by the algorithm HoughLinesP
        :param apply_lpf: apply a low pass filter on image using fft and multiplying the image fft with a low pass filter using the function 
                          image_low_pass_filter
        :param show_im: open image with lines to show user the image with the lines in a pop-up fashion

        :return value: returns if calibration was successful or not. 
        """
        # define constants
        AUTO_CALIBRATION_MAXIMUM_LINES = 1000
        OFFSET = 600
        flag = True

        if backgroundPath is None: backgroundPath = os.path.join(folder_path, 'extra_files\\background.bmp')
        img = cv2.imread(backgroundPath)
        gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

        # find brightest pixel to cutout fiber if exists in image
        ind = np.unravel_index(np.argmax(gray, axis=None), gray.shape)

        if apply_lpf:
            gray = self.image_low_pass_filter(gray)

        gray = gray[:][ind[0] + OFFSET:]
        img = img[:][ind[0] + OFFSET:]

        edges = cv2.Canny(gray, 150, 200, apertureSize=3)

        lines = cv2.HoughLinesP(
                    edges, # Input edge image
                        1, # Distance resolution in pixels
                    np.pi/180, # Angle resolution in radians
                    threshold=100, # Min number of votes for valid line
                    minLineLength=5, # Min allowed length of line
                    maxLineGap=10 # Max allowed gap between line for joining them
                    )
        if not lines:
            print("HoughLinesP did not detect any line! make sure that both sides of the forks are in the image")
            flag=False

        if len(lines) > AUTO_CALIBRATION_MAXIMUM_LINES:
            print("HoughLinesP detected to many lines! try using a low pass filter using the flag apply_lpf=True in function mm_to_pxl_auto_calibration \
                  to filter some noise\
                  from the image")
            flag=False

        if flag:
            # abs on delta x is added to verify lines are perpendicular to the x axis
            leftside_lines = np.array([line[0] for line in lines if line[0][0] < len(gray) / 2 and line[0][2] < len(gray) / 2 and abs(line[0][0] - line[0][2]) <= 1])
            rightside_lines = np.array([line[0] for line in lines if line[0][0] > len(gray) / 2 and line[0][2] > len(gray) / 2 and abs(line[0][0] - line[0][2]) <= 1])
            leftmost_line_in_rightside = rightside_lines[np.argmin(rightside_lines[:, 0])]
            rightmost_line_in_leftside = leftside_lines[np.argmax(leftside_lines[:, 0])]
            # calc distance between the two lines
            d = self.distance_between_lines_in_place(leftmost_line_in_rightside, rightmost_line_in_leftside)

            print(f"distance in pixels between lines {d}")
            self.mm_to_pxl = self.DISTANCE_BETWEEN_FORKS_IN_INVENTOR / (d)
            print(f"pixel to mm calibration: { self.mm_to_pxl }")
        
        # for debugging
        if print_all_lines and lines:
            for line in lines:
                x0, y0, x1, y1 = line[0]
                cv2.line(img, (x0, y0), (x1, y1), (0, 0, 255), 2)
        else:
            cv2.line(img, leftmost_line_in_rightside[:2], leftmost_line_in_rightside[2:], (0, 0, 255), 2)
            cv2.line(img, rightmost_line_in_leftside[:2], rightmost_line_in_leftside[2:], (0, 0, 255), 2)

        # for debugging
        if show_im:
            print("showing image with lines. to exit image and continue the run code press '0'")
            cv2.imshow('image', img)
            cv2.waitKey(0)
            cv2.destroyAllWindows()

        cv2.imwrite(backgroundPath.replace(backgroundPath.split(".")[-1], "with_lines_used_in_calibration.bmp"), img)        
        return flag

    def image_low_pass_filter(image1, e_x=500, e_y=500):
        # taken from: https://stackoverflow.com/questions/54641616/low-pass-filter-for-blurring-an-image 
        # e_x, e_y - size of circle
        #convert image to numpy array
        image1_np=np.array(image1)

        #fft of image
        fft1 = fftpack.fftshift(fftpack.fft2(image1_np))

        #Create a low pass filter image
        x,y = image1_np.shape[0],image1_np.shape[1]
        
        #create a box 
        bbox=((x/2)-(e_x/2),(y/2)-(e_y/2),(x/2)+(e_x/2),(y/2)+(e_y/2))

        low_pass=Image.new("L",(image1_np.shape[0],image1_np.shape[1]),color=0)

        draw1=ImageDraw.Draw(low_pass)
        draw1.ellipse(bbox, fill=1)

        low_pass_np=np.array(low_pass)

        #multiply both the images
        filtered=np.multiply(fft1,low_pass_np.T)

        #inverse fft
        ifft2 = np.real(fftpack.ifft2(fftpack.ifftshift(filtered)))
        ifft2 = np.maximum(0, np.minimum(ifft2, 255))
        ifft2_as_int = ifft2.astype(np .uint8)

        #save the image
        return ifft2_as_int

    def createPathForTemperatureMeasurement(self, path = None, saveConfig = False):
        extraFilesPath = ''
        if path:
            extraFilesPath = path + 'extra_files\\'
        else:
            # --- create path and save config -----
            timeStamp = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            path = 'U:\\Lab_2021-2022\\Experiment_results\\Temperature\\%s\\' %timeStamp
            extraFilesPath = path + 'extra_files\\'
            if not os.path.exists(path): os.makedirs(path)
            if not os.path.exists(extraFilesPath): os.makedirs(extraFilesPath)
            if saveConfig: self.saveConfigTable(path = extraFilesPath)
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

    def fitTemperaturesFromGaussianFitResults(self, gaussianFitResult,alpha = None, extraFilesPath = None, plotResults = True):
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
    

    def oneShotFitGuassian(self, path = None, PrePulseDuration= 1, auto_calibrate_mm_to_pxl=False, PLOT_IMG=True):
        self.takePicturesAndGaussianFit(path=path, PrePulseDuration=[PrePulseDuration], auto_calibrate_mm_to_pxl=auto_calibrate_mm_to_pxl, PLOT_IMG=PLOT_IMG)

    def takePicturesAndGaussianFit(self, path = None, PrePulseDurations = np.arange(1,9), createVideo = False, fit_for_alpha = False, auto_calibrate_mm_to_pxl=False, PLOT_IMG=False):
        extraFilesPath = ''
        if path:
            # if received path, assume there's no need to take new pictures and analyze these photos.
            extraFilesPath = path + 'extra_files\\'
        else:
            if not self.connectCamera(): return
            path, extraFilesPath = self.createPathForTemperatureMeasurement(path, saveConfig= True)
            # ---- Take photos  -----
            for ppd in PrePulseDurations:
                # ---- update prepulse duration ----
                self.updateValue("PrePulse_duration", float(ppd), update_parameters = True)
                # ---- capture and average pictures ------
                imgName = 'PrePulse_duration={:04.1f}.bmp'.format(ppd)
                self.camera.saveAverageImage(path + imgName, NAvg = self.NAvg, NThrow = self.NThrow,RGB = False)

            # ---- Take background picture ----
            self.updateValue("PrePulse_duration", float(50), update_parameters = True) # Presumbaly, after 20ms there's no visible cloud.
            bg_path = os.path.join(extraFilesPath, 'background.bmp')
            self.camera.saveAverageImage(bg_path, NAvg=self.NAvg, NThrow=self.NThrow, RGB=False)

        # --- try to auto-calibrate mm_to_pxl using background image --- 
        if auto_calibrate_mm_to_pxl:
            r_val = self.mm_to_pxl_auto_calibration(backgroundPath=bg_path)
            if not r_val:
                print("Failed to finish temprature measurement because mm_to_pxl_auto_calibration failed.")
                return

        # --- Gaussian Fit to results (to each photo) ----
        gaussianFitResult = self.gaussianFitAllPicturesInPath(path, backgroundPath=bg_path, saveFitsPath = extraFilesPath, imgBounds= self.imgBounds, PLOT_IMG=PLOT_IMG)
        return gaussianFitResult


    def measureTemperature(self, path = None, PrePulseDurations = np.arange(1,9), createVideo = False, fit_for_alpha = False, auto_calibrate_mm_to_pxl=False):
        extraFilesPath = ''
        gaussianFitResult = self.takePicturesAndGaussianFit(path, PrePulseDurations, createVideo, fit_for_alpha, auto_calibrate_mm_to_pxl)
        
        # ---- Take Gaussian fit results and get temperature (x,y) and launch speed -----------
        v_launch, alpha, v_launch_popt, v_launch_cov = self.fitVy_0FromGaussianFitResults(gaussianFitResult=gaussianFitResult, extraFilesPath=extraFilesPath, plotResults=True)
        time_vector = np.array([res[0] for res in gaussianFitResult])
        y_position_vector = np.array([res[-1][2] for res in gaussianFitResult])

        # ------ Use v-launch fit to go from pxl to mm ------
        gaussian_amplitude_vector = alpha * np.array([res[-1][0] for res in gaussianFitResult])
        y_position_vector = alpha * y_position_vector
        x_position_vector = alpha * np.array([res[-1][1] for res in gaussianFitResult])
        sigma_x_vector = alpha * np.array([res[-1][3] for res in gaussianFitResult])
        sigma_y_vector = alpha * np.array([res[-1][4] for res in gaussianFitResult])

        # -------- Fit for x and y temperatures ----
        T_x, T_y = self.fitTemperaturesFromGaussianFitResults(gaussianFitResult, alpha, extraFilesPath, plotResults = True)

        d = {'V_y Launch': v_launch, 'T_x': T_x,'T_y': T_y, 'alpha': alpha}
        print(d)
        try:
            with open(f'{extraFilesPath}results.json', 'w') as file: json.dump(d, file, indent=4)
        except Exception: pass

        if createVideo: self.createVideoFromPath(path,saveFilePath = extraFilesPath + 'video.avi')
        plt.close('all')
        return d

    '''optimizePGC is used to run over different parameters.
        Parameters are defined within function.
    '''
    def optimizePGC(self, PrePulseDurations = np.arange(1,7)):
        path, extraFilesPath = self.createPathForTemperatureMeasurement()

        # ---- Take background picture ----
        self.updateValue("PrePulse_duration", float(50),
                         update_parameters=True)  # Presumbaly, after 20ms there's no visible cloud.
        self.camera.saveAverageImage(extraFilesPath + 'background.bmp', NAvg=self.NAvg, NThrow=self.NThrow, RGB=False)

        # ---- Define parameters space ------
        xs_key = "PGC_final_freq"
        xs = np.linspace(97e6, 90e6,7) # e.g., PGCFinalFreq
        ys_key = "PGC_final_amp"
        ys = np.linspace(0.3, 0.05,7) # e.g., PGCFinalAmp

        # --- Begin measurement -----
        for ppd in PrePulseDurations:
            self.updateValue("PrePulse_duration", float(ppd))
            for x in xs:
                self.updateValue(xs_key, float(x))
                for y in ys:
                    self.updateValue(ys_key, float(y))
                    self.update_parameters()
                    imgName = 'PrePulse=%s;%s=%s;%s=%s.bmp' %(str(ppd),str(xs_key),str(x),str(ys_key), str(y))
                    self.camera.saveAverageImage(path + imgName, NAvg=self.NAvg, NThrow=self.NThrow, RGB=False)

        print('Finished  optimizePGC. running analyzePGCOptimization...')
        return self.analyzePGCOptimization(path = path)

    def analyzePGCOptimization(self, path, backgroundPath = None):
        # ---- fit gaussians to all photos -----
        gaussianFitsAndParams = self.gaussianFitAllPicturesInPath(path, backgroundPath, saveFitsPath=None, imgBounds=self.imgBounds)
        np.save(path + '\\gasuusainFitResultsresults.npy', gaussianFitsAndParams, allow_pickle=True)
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
                print('Results for x = {}, y = {} are: {}K'.format(str(x),str(y),str(temp)))

        self.plot3D(xs_for_plot,ys_for_plot,uniqueTemperatures)
        return (xs_for_plot,ys_for_plot,uniqueTemperatures)

    def createVideoFromPath(self, path, saveFilePath = None):
        # Taken from: https://theailearner.com/2018/10/15/creating-video-from-images-using-opencv-python/
        if not saveFilePath: saveFilePath = path + 'video.avi'
        img_array = []
        fileNames = glob.glob(path + '*.bmp')
        fileNames.sort() #
        for filename in fileNames:
            img = cv2.imread(filename)
            height, width, layers = img.shape
            size = (width, height)
            img_array.append(img)

        out = cv2.VideoWriter(filename = saveFilePath, fourcc = cv2.VideoWriter_fourcc(*'DIVX'), fps = 3, frameSize= size)
        for i in range(len(img_array)):
            out.write(img_array[i])
        out.release()

    def gaussianFitAllPicturesInPath(self, path, backgroundPath = None, saveFitsPath = None,imgBounds = None):
        if backgroundPath is None: backgroundPath = path + 'extra_files\\background.bmp'
        if imgBounds is None: imgBounds = self.imgBounds
        files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
        res = []
        for f in files:
            try:
                fullFileName = path + f
                fileRes = []
                f = f.replace('.jpeg','').replace('.bmp','').replace(',',';')
                if f.find('background'): pass
                ## -- get x and y values from file name ----
                keysAndValues = f.split(';')
                keys = []
                for key in keysAndValues:
                    value = key[key.find('=') + 1:]
                    only_key = key[:key.find('=')]
                    fileRes.append(float(value))
                    keys.append(only_key)
                ## -- get z value from gaussian fit ----
                gaussianFit = self.GaussianFit(fullFileName, background_file = backgroundPath,saveFitsPath = saveFitsPath,imgBounds=imgBounds)
                if gaussianFit is None: continue # if fit returned None, meaning the fit failed (sigma is out of self.sigma_bounds), discard this results and continue to the next fit
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
    def GaussianFit(self, file_name_for_fit,background_file, saveFitsPath = None,imgBounds = None,X_PIXEL_LEN=1544, Y_PIXEL_LEN=2064, CROP_IMG_SIZE =180 ,PLOT_IMG=False, PLOT_SLICE =False):
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

    # define model function and pass independant variables x and y as a list
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

    def plotDataAndFit(self, x_data, y_data, fitFunc,fitParams,title = '',ylabel = '',xlabel = 'Time [ms]', saveFilePath = None, show = False):
        x_vector_for_fit = np.linspace(min(x_data), max(x_data), 100)
        plt.figure()
        plt.plot(x_data,y_data,'x', label = 'Data')
        plt.plot(x_vector_for_fit, fitFunc(x_vector_for_fit, *fitParams), '-', label = 'Fit')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.legend()
        if saveFilePath: plt.savefig(saveFilePath)
        if show: plt.show()

#r = optimizePGC(c)
if __name__ == "__main__":
    experiment = CoolingSequenceOptimizer(Config.config)
    # experiment.measureTemperature()
    #experiment.optimizePGC()