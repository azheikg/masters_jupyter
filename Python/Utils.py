#!/usr/bin/python3
"""Utilities and common functions
"""

__author__ = "Kevin Gema"
__copyright__ = "Copyright 2018, Kevin Gema"
__credits__ = ["Kevin Gema"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = ""

import numpy as np
import math

import Constants

# -----------------------------------------------------------------------------
# Helper datatypes
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
def is_main_module():
    """
    Returns whether this notebook is the main module
    i.e. not being run from anotherr notebook
    """
    return __name__ == '__main__' and '__file__' not in globals()

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

def deg2rad(degrees):
    """
    Convert degrees to radians
    """
    return (np.pi/180) * degrees

def rad2deg(radians):
    """
    Convert radians to degrees
    """
    return (180/np.pi) * radians

def GetNoisePSD(T, F):
    """
    Get the noise power spectral density for a given temperature (T), and 
    noise figure (F), N0 = kT0F
    
    T in Kelvin
    F in dB
    
    returns the noise power spectral density, N0, in W/Hz
    """
    return ((math.pow(10, (F/10))) * T * Constants.BOLTZMANS_CONSTANT)