import inspect
import re
#from turtle import colorpartial
# import matplotlib
# import matplotlib.axes
import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
# import matplotlib.ticker as mticker
# from matplotlib.patches import FancyArrowPatch
# from mpl_toolkits.mplot3d import Axes3D
# from mpl_toolkits.mplot3d import proj3d
# from astropy.time import Time
# from datetime import datetime
# from tqdm import tqdm

# from functools import partial

# import scipy.stats as stats
# from scipy.stats import rayleigh, uniform, norm, chi2, gamma
# #import numba as nb
# # from numba import jit, vectorize
# from math import sin, cos, sqrt
# import pyautogui

# #from DataAnalysis import LIASignal, SQUID
# import pandas as pd
# from evidently.report import Report
# from evidently.metric_preset import DataDriftPreset



def check(arg):
    """
    Print information of input arg
    
    Example
    ------
    import numpy as np
    
    a = np.zeros((2, 4))
    
    check(a)
    
    a+=1
    
    check(a)
    
    check(len(a))

    TERMINAL OUTPUT: 

    d:\Yu0702\casper-gradient-code\\testofcheckpoint.py @45 a : ndarray(array([[0., 0., 0., 0.], [0., 0., 0., 0.]])) [shape=(2, 4)]
    
    d:\Yu0702\casper-gradient-code\\testofcheckpoint.py @47 a : ndarray(array([[1., 1., 1., 1.], [1., 1., 1., 1.]])) [shape=(2, 4)]
    
    d:\Yu0702\casper-gradient-code\\testofcheckpoint.py @48 len(a) : int(2)
    
    d:\Yu0702\casper-gradient-code\\testofcheckpoint.py @49 a.shape : tuple((2, 4)) [len=2]


    Copyright info: 
    ------
    Adopted from https://gist.github.com/HaleTom/125f0c0b0a1fb4fbf4311e6aa763844b
    
    Author: Tom Hale 
    
    Original comment: Print the line and filename, function call, the class, str representation and some other info
                    Inspired by https://stackoverflow.com/a/8856387/5353461
    
    
    
    """
    frame = inspect.currentframe()
    callerframeinfo = inspect.getframeinfo(frame.f_back)
    try:
        context = inspect.getframeinfo(frame.f_back).code_context
        caller_lines = ''.join([line.strip() for line in context])
        m = re.search(r'check\s*\((.+?)\)$', caller_lines)
        if m:
            caller_lines = m.group(1)
            position = str(callerframeinfo.filename) + " line " + str(callerframeinfo.lineno)

            # Add additional info such as array shape or string length
            additional = ''
            if hasattr(arg, "shape"):
                additional += "[shape={}]".format(arg.shape)
            elif hasattr(arg, "__len__"):  # shape includes length information
                additional += "[len={}]".format(len(arg))

            # Use str() representation if it is printable
            str_arg = str(arg)
            str_arg = str_arg if str_arg.isprintable() else repr(arg)

            print(position, "" + caller_lines + " : ", end='')
            print(arg.__class__.__name__ + "(" + str_arg + ")", additional)
        else:
            print("check: couldn't find caller context")
    finally:
        del frame
        del callerframeinfo



def Lorentzian(x, center, FWHM, area:float=1., offset:float=0.):
    """
    Return the value of the Lorentzian function 
        offset + 0.5*FWHM*area / (np.pi * ( (x-center)**2 + (0.5*FWHM)**2 )      )

                           FWHM A
        offset + ───────────────────────
                  2π ((x-c)^2+(FWHM/2)^2 )
    
    Parameters
    ----------
    
    x : scalar or array_like
        argument of the Lorentzian function 
    center : scalar
        the position of the Lorentzian peak
    FWHM : scalar
        full width of half maximum (FWHM) / linewidth of the Lorentzian peak
    area : scalar
        area under the Lorentzian curve (without taking offset into consideration)
    offset : scalar
        offset for the curve
    
    
    Returns
    -------
    the value of the Lorentzian function : ndarray or scalar
        
    Examples
    --------
    >>> 

    Reference
    ---------- 
    Null

    """
    return offset + 0.5*abs(FWHM)*area / (np.pi * ((x-center)**2 + (0.5*FWHM)**2))



def get_FWHMin(x, y):
    """
    Calculate the Full Width at Half Maximum (FWHM) of a dip.

    Parameters:
        x (array-like): The x-values of the curve.
        y (array-like): The y-values of the curve.

    Returns:
        float: The FWHM of the curve.
    """
    # Ensure inputs are numpy arrays
    x = np.array(x)
    y = np.array(y)
    
    # Find the maximum value of y and its half-maximum
    y_min = np.amin(y)
    Twice_min = y_min * 2.0
    
    # Find indices where y crosses the half-maximum
    # check(np.where(y <= Twice_min))
    indices = np.where(y <= Twice_min)[0]
    if len(indices) < 2:
        raise ValueError("Cannot calculate FWHM: The curve does not have two points crossing the half-maximum.")

    # Extract the first and last indices crossing the half-maximum
    left_index = indices[0]
    right_index = indices[-1]
    
    # # Interpolate to find more precise crossing points
    # x_left = np.interp(Twice_min, [y[left_index - 1], y[left_index]], [x[left_index - 1], x[left_index]])
    # x_right = np.interp(Twice_min, [y[right_index], y[right_index + 1]], [x[right_index], x[right_index + 1]])
    x_left = x[left_index]
    x_right = x[right_index]

    # Calculate FWHM
    FWHMin = abs(x_right - x_left)

    return FWHMin
