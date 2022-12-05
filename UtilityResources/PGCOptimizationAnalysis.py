from os import listdir
from os.path import isfile, join
import cv2
import numpy as np
import pylab as plt
import scipy.optimize as opt
import matplotlib.pyplot as mtl
from pathlib import Path

from mpl_toolkits.mplot3d import Axes3D


class PCGOptimizationAnalyzer:
    def gaussianFitAllPicturesInPath(self, path, backgroundPath, imgBounds = None):
        files = [f for f in listdir(path) if isfile(join(path, f))]
        res = []
        for f in files:
            try:
                fullFileName = path + '\\' + f
                fileRes = []
                f = f.replace('.bmp','').replace(',',';')
                ## -- get x and y values from file name ----
                keysAndValues = f.split(';')
                keys = []
                for key in keysAndValues:
                    value = key[key.find('=') + 1:]
                    only_key = key[:key.find('=') + 1]
                    fileRes.append(float(value))
                    keys.append(only_key)
                ## -- get z value from gaussian fit ----
                gaussianFit = self.GaussianFit(fullFileName, backgroundPath,imgBounds=imgBounds, PLOT_IMG = False)
                fileRes.append(keys)
                fileRes.append(gaussianFit)
                res.append(fileRes)
            except Exception:
                pass
        return (res)


    #fit a gaussian to subtracted images
    # INPUT:
    #   Y_PIXEL_LEN - number of vertical pixels
    #   X_PIXEL_LEN - number of horizontal pixels
    #   the name of the file for fit with png termination
    #   CROP_IMG_SIZE - the size of image in each dimension after cropping is 2*CROP_IMG_SIZE
    def GaussianFit(self, file_name_for_fit,background_file, saveFitsPath = None,imgBounds = None,X_PIXEL_LEN=1544, Y_PIXEL_LEN=2064, CROP_IMG_SIZE =300 ,PLOT_IMG=False, PLOT_SLICE =False):
        #ImgToFit = cv2.imread('./Images/SubtractedImages/' + file_name_for_fit, 0)
        #fileName = file_name_for_fit.split('\\')[-1].split('.')[0]
        fileName = Path(file_name_for_fit).stem
        backgroundImg = cv2.imread(background_file, 0)
        ImgToFit = cv2.imread(file_name_for_fit, 0)
        ImgToFit = cv2.subtract(ImgToFit, backgroundImg)


        if not imgBounds:
            h, w, c = ImgToFit.shape
            imgBounds = (0,0,w,h) # that is, don't crop image

        if imgBounds: ImgToFit = ImgToFit[imgBounds[1]:imgBounds[3],imgBounds[0]:imgBounds[2]] #
        img_max_index = [np.argmax(np.sum(ImgToFit, axis=0)), np.argmax(np.sum(ImgToFit, axis=1))]

        #img_max_index[1] = np.argmax(np.sum(ImgToFit[10:][img_max_index[0]-CROP_IMG_SIZE:img_max_index[0]+CROP_IMG_SIZE], axis=1))
        print(img_max_index)
        # Parameters
        X_UPPER_BOUND = int(img_max_index[0] + CROP_IMG_SIZE)
        X_LOWER_BOUND = int(img_max_index[0] - CROP_IMG_SIZE)
        if X_LOWER_BOUND < 0: X_LOWER_BOUND = 0
        if X_UPPER_BOUND > X_PIXEL_LEN: X_UPPER_BOUND= X_PIXEL_LEN
        EFFECTIVE_X_PIXEL_LEN = X_UPPER_BOUND - X_LOWER_BOUND

        Y_UPPER_BOUND = int(img_max_index[1] + CROP_IMG_SIZE)
        Y_LOWER_BOUND = int(img_max_index[1] - CROP_IMG_SIZE)
        if Y_LOWER_BOUND < 0: Y_LOWER_BOUND = 0
        if Y_UPPER_BOUND > Y_PIXEL_LEN: Y_UPPER_BOUND= Y_PIXEL_LEN
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
        initial_guess = (EffectiveImg[img_max_index[1]][img_max_index[0]], img_max_index[1], img_max_index[0], EFFECTIVE_X_PIXEL_LEN / 10, EFFECTIVE_Y_PIXEL_LEN / 10, 0, 10)
        print(initial_guess)
        popt, pcov = opt.curve_fit(self.twoD_Gaussian, (x, y), data_noisy, p0=initial_guess, bounds = (0,np.inf))


        # plot the results
        data_fitted = self.twoD_Gaussian((x, y), *popt)
        sigma = [popt[3], popt[4]]
        if PLOT_IMG or saveFitsPath:
            fig, ax = plt.subplots(1, 1)
            plt.title(fileName)
            plt.text(0.88, 0.95, '\u03C3_x =' + '%.2f' % sigma[0] + '\n' + '\u03C3_y = ' + '%.2f' % sigma[1], color='white',
                 fontsize=16, style='italic', weight='bold', horizontalalignment='center', verticalalignment='center',
                 transform=ax.transAxes, bbox=dict(facecolor='gray', alpha=0.5))
            ax.imshow(data_noisy.reshape(EFFECTIVE_Y_PIXEL_LEN, EFFECTIVE_X_PIXEL_LEN), cmap=plt.cm.jet, origin='bottom',
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
        g = offset + amplitude * np.exp(- (a * ((x - xo) ** 2) + c * ((y - yo) ** 2)))
        return g.ravel()

    def plot3D(self, xs, ys, zs):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # Plot the surface
        ax.plot_trisurf(xs, ys, zs)


a = PCGOptimizationAnalyzer()
path = 'U:\\Lab_2021-2022\\Experiment_results\\Temperature\\20221030_151451\\'
bgFilePath = path + 'extra_files\\background.bmp'
#p = a.GaussianFit(file_name_for_fit= path + 'PrePulse_duration=00.5.bmp', background_file = bgFilePath,saveFitsPath=path+'extra_files\\', imgBounds= (580,200,1600,1450), PLOT_IMG= True, PLOT_SLICE=True)
p = a.gaussianFitAllPicturesInPath(path, backgroundPath= bgFilePath,imgBounds=(580,200,1600,1450) )


#
xs = np.array([item[1] for item in p])
ys = np.array([item[2] for item in p])
prePulseD = np.array([item[0] for item in p])
gaussian_amp = np.array([item[-1][0] for item in p])
sigmas_x = np.array([item[-1][3] for item in p])
sigmas_y = np.array([item[-1][4] for item in p])
#
#
# def tempFunc(s, a, b):
#     return np.sqrt(a + b * s**2)
#
# alpha = 0.0120 # mm/pixel
# Rb_m = 1.443e-25 # [Kg]
# Kb = 1.38e-23 # [J/K]
#
# uniqueFrqs = np.unique(xs)
# uniqueAmps = np.unique(ys)
# uniqueTemperatures, uniqueArbTemperatures = [], []
# xs_for_plot, ys_for_plot = [],[]
# for f in uniqueFrqs:
#     for amp in uniqueAmps:
#
#         indices = np.argwhere((xs == f) & (ys == amp))
#         ppds, sigmas = [],[]
#         for i in indices:
#             i= i[0] # since np.argwhere returns an array with the index inside
#             ppd = prePulseD[i]
#             s = np.sqrt(sigmas_x[i] ** 2 + sigmas_y[i] **2)
#             #s = s * alpha
#             ppds.append(ppd)
#             sigmas.append(s)
#         plt.plot(ppds, sigmas)
#         plt.show()
#         popt, pcov = opt.curve_fit(tempFunc, ppds, sigmas)
#         arbTemp = popt[1]
#         temp = arbTemp*Rb_m/Kb*1e6 #[uK]
#         # Arrange data for plotting
#         xs_for_plot.append(f)
#         ys_for_plot.append(amp)
#         uniqueArbTemperatures.append(arbTemp)
#         uniqueTemperatures.append(temp)
#
# a.plot3D(append,ys_for_plot,uniqueArbTemperatures)
# np.save(path +'\\gasuusainFitResultsresults.npy', p, allow_pickle=True)
#
