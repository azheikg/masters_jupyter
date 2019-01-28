#!/usr/bin/python3
"""Pulse Repition Frequency
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
import Utils

# -----------------------------------------------------------------------------
# PFR Minimum and maximum
def PrfMin(antennaLength, altitude):
    """
    Get the minimum PRF for a specific altitude and antenna azimuth length.

    antennaLength in meters
    altititude in meters
    
    returns the minimum prf in Hz
    """
     # convert the Re and altitude paramters to km. 
    # Afterwards convert back to m  
    _d = 2/ antennaLength
    _inner = Constants.STANDARD_GRAVITATIONAL_CONSTANT / ((Constants.EARTH_RADIUS + altitude)/1000)
    _prfMin = _d * (np.sqrt(_inner) * 1000)
    return _prfMin

def PrfMax(pulseDuration, sceneDuration):
    """
    Calculates the maximim PRF using the pulse duration along with the scene duration.
    
    pulseDuration in seconds
    sceneDuration in seconds
    
    returns the maximum PRF in Hz
    """
    return (1 / ((2* pulseDuration) + sceneDuration))

def PrfMaxDutyCycle(dutyCycle, sceneDuration):
    """
    Calculates the maximum PRF using the transmitter duty cycle with the scene duration.

    returns the maximum PRF in Hz
    """
    return (1 - 2*dutyCycle) / (sceneDuration)

# -----------------------------------------------------------------------------
# PRF Blind ranges and Nadir Echoes
def PrfBlindRanges(numIntervals, nearRangeDelay, farRangeDelay, pulseDuration):
    """
    Calculates the lower and upper limits for PRF over a number of intervals.
    These limits represent the PRF blind ranges.
    
    numIntervals is a scalar
    nearRangeDelay in seconds
    farRangeDelay in seconds
    pulseDuration in seconds
    
    returns BOTH the lower and upper PRF limits
    """
    _lArr = []
    _uArr = []
    for n in range(1, numIntervals):
        _lower = (n - 1) / (nearRangeDelay - pulseDuration)
        _upper = (n) / (farRangeDelay + pulseDuration)
        _lArr.append(_lower)
        _uArr.append(_upper)
    _lArr = np.asarray(_lArr)
    _uArr = np.asarray(_uArr)
   
    return _lArr, _uArr

def PrfBlindRangesDutyCycle(numIntervals, nearRangeDelay, farRangeDelay, dutyCycle):
    """
    Calculates the lower and upperr limits for PRF over a number of intervals,
    using duty cycle instead of transmitter pulse length.
    These limits represent the PRF blind ranges.

    returns BOTH the lower and upper PRF limits
    """
    _lArr= []
    _uArr = []
    for n in range(1, numIntervals):
        _lArr.append( (n - 1 + dutyCycle) / nearRangeDelay )
        _uArr.append( (n - dutyCycle) / farRangeDelay )

    _lArr = np.asarray(_lArr)
    _uArr = np.asarray(_uArr)
    return _lArr, _uArr


def NadirDelay(altitude):
    """
    Get the nadir echo delay (propagation delay) for a specific altitude
    """
    tnadir = (2*altitude) / Constants.SPEED_OF_LIGHT
    return tnadir

def NadirEchoRanges(numIntervals, nearRangeDelay, farRangeDelay, pulseDuration, nadirDelay):
    """
    Calculates the lower and upper limits for the Nadir echoes over a number of intervals.
    
    numIntervals is a scalar
    nearRangeDelay in seconds
    farRangeDelay in seconds
    pulseDuration in seconds
    nadirDelay in seconds
    
    returns BOTH the lower and upper Nadir echo limits
    """
    _nlArr = []
    _nuArr = []
    for m in range(1, numIntervals):
        _lower = (m - 1) /(nearRangeDelay - pulseDuration - nadirDelay)
        _upper = (m) / (farRangeDelay + pulseDuration - nadirDelay)
        _nlArr.append(_lower)
        _nuArr.append(_upper)
    _nlArr = np.asarray(_nlArr)
    _nuArr = np.asarray(_nuArr)    
    return _nlArr, _nuArr

def NadirEchoRangesDutyCycle(numIntervals, nearRangeDelay, farRangeDelay, dutyCycle, nadirDelay):
    """
    Calculates the lower and upper limits for the Nadir echoes over a number of intervals.

    returns Both the lower and upper nadir echoe limits
    """
    _nlArr = []
    _nuArr = []
    for m in range(1, numIntervals):
        _nlArr.append( (m-1+dutyCycle) / (nearRangeDelay - nadirDelay))
        _nuArr.append( (m - dutyCycle) / (farRangeDelay - nadirDelay))

    _nuArr = np.asarray(_nuArr)
    _nlArr = np.asarray(_nlArr)
    return _nlArr, _nuArr

def NadirEchoRangesCurlander(numIntervals, nearRangeDelay, farRangeDelay, pulseDuration, nadirDelay):
    """
    Calculates the lower and upper limits for the Nadir echoes over a number of intervals.
    
    numIntervals is a scalar
    nearRangeDelay in seconds
    farRangeDelay in seconds
    pulseDuration in seconds
    nadirDelay in seconds
    
    returns BOTH the lower and upper Nadir echo limits
    """
    _nlArr = []
    _nuArr = []
    for m in range(0, numIntervals):
        _lower = (m) /(nearRangeDelay - 2*pulseDuration - nadirDelay)
        _upper = (m) / (farRangeDelay - nadirDelay)
        _nlArr.append(_lower)
        _nuArr.append(_upper)
    _nlArr = np.asarray(_nlArr)
    _nuArr = np.asarray(_nuArr)
    
    return _nlArr, _nuArr