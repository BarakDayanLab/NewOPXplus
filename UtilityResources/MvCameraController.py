from PIL import Image, ImageOps
import numpy as np
from os import listdir
from os.path import isfile, join

import cv2
import scipy.optimize as opt
import matplotlib.pyplot as mtl
import pylab as plt
import time

from mvIMPACT import acquire

try:
    from UtilityResources import MvCamera

except:
    import MvCamera
    #print('Warning! could not import MvCamera module; in pgcWidget.py')


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
def printError(s=''): print(f"{bcolors.WARNING}{s}{bcolors.ENDC}")
def printGreen(s=''): print(f"{bcolors.OKGREEN}{s}{bcolors.ENDC}")


class MvCameraController:
    def __init__(self, deviceSerial='FF006583'):
        try:
            self.camera = MvCamera.MvCamera(deviceSerial=deviceSerial)
            printGreen("connected to camera")
        except acquire.EDeviceManager:
            printError("Could not connect to camera. Camera was probably already connected.")
            printError("acquire.EDeviceManager")
            return None

    def captureImage(self):
        img, _ = self.camera.CaptureImage()
        return img

    def captureGrayscaleImage(self):
        img = self.captureImage()
        return ImageOps.grayscale(img)

    def captureImageAsArray(self, RGB=False):
        if RGB:
            return np.array(self.captureImage())
        else:
            return np.array(self.captureGrayscaleImage())

    def saveAverageImage(self, path, NAvg = 3, NThrow = 2,RGB = False):
        for i in range(NThrow):
            self.captureGrayscaleImage()  # Throw images till MOT is stable
        photos = []
        for i in range(NAvg):
            photos.append(np.array(self.captureImageAsArray(RGB=RGB)))
        avgPhoto = np.average(np.array(photos), axis=0)
        img = Image.fromarray(avgPhoto)

        img.convert("L").save(path)

    def gaussianFitAllPicturesInPath(self, path, backgroundPath):
        files = [f for f in listdir(path) if isfile(join(path, f))]
        res = []
        for f in files:
            fullFileName = path + '\\' + f
            fileRes = []
            f = f.replace('.jpeg','')
            ## -- get x and y values from file name ----
            keysAndValues = f.split(';')
            for key in keysAndValues:
                value = key[key.find('=') + 1:]
                fileRes.append(float(value))
            ## -- get z value from gaussian fit ----
            gaussianFit = self.GaussianFit(fullFileName, backgroundPath, PLOT_IMG = True)
            fileRes.append(gaussianFit)
            res.append(fileRes)
        return (res)


    #fit a gaussian to subtracted images
    # INPUT:
    #   Y_PIXEL_LEN - number of vertical pixels
    #   X_PIXEL_LEN - number of horizontal pixels
    #   the name of the file for fit with png termination
    #   CROP_IMG_SIZE - the size of image in each dimension after cropping is 2*CROP_IMG_SIZE
    def GaussianFit(self, file_name_for_fit,background_file, X_PIXEL_LEN=1544, Y_PIXEL_LEN=2064, CROP_IMG_SIZE =180 ,PLOT_IMG=False, PLOT_SLICE =False):
        #ImgToFit = cv2.imread('./Images/SubtractedImages/' + file_name_for_fit, 0)
        backgroundImg = cv2.imread(background_file, 0)
        ImgToFit = cv2.imread(file_name_for_fit, 0)[CROP_IMG_SIZE:,CROP_IMG_SIZE:] #- backgroundImg

        img_max_index = [np.argmax(np.sum(ImgToFit, axis=0)), np.argmax(np.sum(ImgToFit, axis=1))]
        #img_max_index[1] = np.argmax(np.sum(ImgToFit[10:][img_max_index[0]-CROP_IMG_SIZE:img_max_index[0]+CROP_IMG_SIZE], axis=1))
        # Parameters
        X_UPPER_BOUND = int(img_max_index[0] + CROP_IMG_SIZE)
        X_LOWER_BOUND = int(img_max_index[0] - CROP_IMG_SIZE)
        EFFECTIVE_X_PIXEL_LEN = X_UPPER_BOUND - X_LOWER_BOUND

        Y_UPPER_BOUND = int(img_max_index[1] + CROP_IMG_SIZE)
        Y_LOWER_BOUND = int(img_max_index[1] - CROP_IMG_SIZE)
        EFFECTIVE_Y_PIXEL_LEN = Y_UPPER_BOUND - Y_LOWER_BOUND

        # Create x and y indices
        x = np.linspace(0, EFFECTIVE_Y_PIXEL_LEN - 1, EFFECTIVE_Y_PIXEL_LEN)
        y = np.linspace(0, EFFECTIVE_X_PIXEL_LEN - 1, EFFECTIVE_X_PIXEL_LEN)
        x, y = np.meshgrid(x, y)
        # crop an effective image
        EffectiveImg = ImgToFit[X_LOWER_BOUND:X_UPPER_BOUND, Y_LOWER_BOUND:Y_UPPER_BOUND]
        plt.imshow(ImgToFit)
        plt.show()
        data_noisy = EffectiveImg.ravel()

        # fit the data
        img_max_index = [np.argmax(np.sum(EffectiveImg, axis=1)), np.argmax(np.sum(EffectiveImg, axis=0))]
        initial_guess = (ImgToFit[img_max_index[0]][img_max_index[1]], img_max_index[0], img_max_index[1], EFFECTIVE_X_PIXEL_LEN, EFFECTIVE_Y_PIXEL_LEN, 0, 10)
        popt, pcov = opt.curve_fit(self.twoD_Gaussian, (x, y), data_noisy, p0=initial_guess)

        # plot the results
        data_fitted = self.twoD_Gaussian((x, y), *popt)
        sigma = [popt[3], popt[4]]
        if PLOT_IMG:
            fig, ax = plt.subplots(1, 1)

            plt.text(0.88, 0.95, '\u03C3_x =' + '%.2f' % sigma[0] + '\n' + '\u03C3_y = ' + '%.2f' % sigma[1], color='white',
                 fontsize=16, style='italic', weight='bold', horizontalalignment='center', verticalalignment='center',
                 transform=ax.transAxes, bbox=dict(facecolor='gray', alpha=0.5))
            ax.imshow(data_noisy.reshape(EFFECTIVE_X_PIXEL_LEN, EFFECTIVE_Y_PIXEL_LEN), cmap=plt.cm.jet, origin='bottom',
                  extent=(x.min(), x.max(), y.min(), y.max()))
            ax.contour(x, y, data_fitted.reshape(EFFECTIVE_X_PIXEL_LEN, EFFECTIVE_Y_PIXEL_LEN), 8, colors='w')

            plt.show()

        if PLOT_SLICE:
            # plot slice
            mtl.figure()
            data_noisy_mat = data_noisy.reshape(EFFECTIVE_X_PIXEL_LEN, EFFECTIVE_Y_PIXEL_LEN)
            Slice = data_noisy_mat[int(EFFECTIVE_X_PIXEL_LEN / 2) - 10]
            mtl.plot(Slice)
            mtl.plot(data_fitted.reshape(EFFECTIVE_X_PIXEL_LEN, EFFECTIVE_Y_PIXEL_LEN)[int(EFFECTIVE_X_PIXEL_LEN / 2) - 10])
            plt.show()
        return popt

    # define model function and pass independant variables x and y as a list
    def twoD_Gaussian(self, x_y, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
        x, y = x_y
        xo = float(xo)
        yo = float(yo)
        a = (np.cos(theta) ** 2) / (2 * sigma_x ** 2) + (np.sin(theta) ** 2) / (2 * sigma_y ** 2)
        b = -(np.sin(2 * theta)) / (4 * sigma_x ** 2) + (np.sin(2 * theta)) / (4 * sigma_y ** 2)
        c = (np.sin(theta) ** 2) / (2 * sigma_x ** 2) + (np.cos(theta) ** 2) / (2 * sigma_y ** 2)
        g = offset + amplitude * np.exp(- (a * ((x - xo) ** 2) + 2 * b * (x - xo) * (y - yo)
                                           + c * ((y - yo) ** 2)))
        return g.ravel()
#
# c = MvCameraController()#(deviceSerial = 'VD000001') #deviceSerial='FF006524')
# path = 'C:\\Users\\orelb\\Downloads\\NATANPGC\\1'
# p = c.gaussianFitAllPicturesInPath(path, backgroundPath= 'C:\\Users\\orelb\\Downloads\\NATANPGC\\1\\b\\background.png')