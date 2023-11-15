import os
import pathlib
import numpy as np
from scipy.interpolate import griddata

#-------------------------------------------------------
# Used for AOM Double-Pass Calibration
#
# Calibration data (file) is created by placing power meter in front of the back exit of TOP1 (after double pass!) and
# running CalibrateLockAOM.py
#
#-------------------------------------------------------

file_directory = str(pathlib.Path(__file__).parent.resolve())
calibration_file = os.path.join(file_directory, 'calibrationFunction_griddata_cubic.npz')

calibration_data = np.load(calibration_file, allow_pickle=True)['data']
calibration_data = dict(enumerate(calibration_data.flatten(), 1))[1]  # Converting load data to dictionary

calibrate = lambda f, op: float(griddata(calibration_data['calibration_data'][0], calibration_data['calibration_data'][1], ([f], [op]), method='cubic'))

