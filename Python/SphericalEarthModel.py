#!/usr/bin/python3
"""Spherical Earth model related functions
"""

from math import sqrt
import numpy as np
import math

import Constants

__author__ = "Kevin Gema"
__copyright__ = "Copyright 2018, Kevin Gema"
__credits__ = ["Kevin Gema"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = ""


# -----------------------------------------------------------------------------
# Look Andgles and Incidence Angles
def IncidenceAngle(Altitude, LookAngle):
    """
    Calculates the incidence angle from the given platform altitude and 
    antenna Look angle
    
    Altitude in meters
    LookAngle in radians
    
    Returns the incidence angle in radians
    """
    RE = Constants.EARTH_RADIUS
    H = Altitude

    return np.arcsin( ((RE + H) / RE) * np.sin(LookAngle) )

def NearFarIncidenceAngleFromBeamwidth(Altitude, MidSwathLookAngle, ElBeamwidth):
    """
    Gets the local near and far incidence angles from the sar orbit altitude,
    the mid-swath look angle and the elevation antenna beamwidth.

    Altitude in meters
    MidSwathLookAngle in radians
    ElBeamwidth in radians

    Returns NearIA and FarIA in radians
    """
    _nia = IncidenceAngle(Altitude, MidSwathLookAngle - (ElBeamwidth/2))
    _fia = IncidenceAngle(Altitude, MidSwathLookAngle + (ElBeamwidth/2))
    return _nia, _fia

def NearFarSlantRangeFromBeamwidth(Altitude, MidSwathLookAngle, ElBeamwidth):
    """
    Get the near and far slant ranges from the sar orbit altitude, the mid-swath look angle and the 
    elevation beamwidth.

    Altitude in meters
    MidSwathLookAngle in radians
    ElBeamwidth in radians

    returns NSR and FSR in meters
    """
    RE = Constants.EARTH_RADIUS
    H = Altitude
    _nia, _fia = NearFarIncidenceAngleFromBeamwidth(H, MidSwathLookAngle, ElBeamwidth)
    _rn = SlantRange(Altitude, _nia, MidSwathLookAngle - (ElBeamwidth/2))
    _rf = SlantRange(Altitude, _fia, MidSwathLookAngle + (ElBeamwidth/2))
    return _rn, _rf

def SwathGroundSwathFromBeamwidth(Altitude, MidSwathLookAngle, ElBeamwidth):
    """
    Get the swath and ground swath widts for a sar orbit altitude, mid-swath look angle and the 
    elevation beamwidth.

    Altitude in meters
    MidSwathLookAngle in radians
    ElBeamwidth in radians

    returns the Swath and GroundSwath in meters
    """
    _rn, _rf = NearFarSlantRangeFromBeamwidth(Altitude, MidSwathLookAngle, ElBeamwidth)
    _wr = _rf - _rn

    # get the midswath incidence angle, 
    # it uses the same formula as the near and far incidence angles, but with beamwidth = 0
    _mia, _ =  NearFarIncidenceAngleFromBeamwidth(Altitude, MidSwathLookAngle, 0)
    
    # ground swath = Swath / sin(mid-swath incidence angle)
    _wgr = _wr / np.sin(_mia)
    return _wr, _wgr


def LookAngle(altitude, incidenceAngle):
    """
    Calculate the look angles corresponding to the specified incidence 
    angles and altitude
    
    altitude in kilometers
    incidenceAngle in radians
    
    returns the look angle in radians
    """
    _num = Constants.EARTH_RADIUS * np.sin(incidenceAngle)
    _denom = Constants.EARTH_RADIUS + altitude
    return np.arcsin(_num / _denom)

def IncidenceAngleFromSlantRange(Altitude, SlantRange):
    """
    Calculate the incidence angle for a given slant range, using the 
    spherical earth model. This equation is useful when working with 
    ambiguities, like when calculating the RASR.
    
    Altitude in meters
    SlantRange in meters
    
    returns the incidence angle in radians
    """
    _num = np.power((1 + (Altitude / Constants.EARTH_RADIUS)), 2) - 1 \
            - np.power((SlantRange / Constants.EARTH_RADIUS), 2)
    _denom = 2* (SlantRange / Constants.EARTH_RADIUS)
    
    return np.arccos(_num / _denom)

# -----------------------------------------------------------------------------
# Slant ranges and round-trip times
def SlantRange(Altitude, IncidenceAngle, LookAngle):
    """
    Calculates the slant range for a given altitude, incidence angle and look 
    angle.
    
    Altitude in meters
    IncidenceAngle in radians
    LookAngle in radians
    
    returns the slant range in meters
    """
    RE = Constants.EARTH_RADIUS
    H = Altitude
    IA = IncidenceAngle
    LA = LookAngle
    _sr = np.sqrt( np.power((RE + H), 2) + (RE**2) - (2*RE*(RE + H)*np.cos(IA - LA)) )

    return _sr

def RoundTripTime(slantRange):
    """
    Calculate the round trip time for a given slant range
    
    slant range in meters.
    
    returns the round trip time in seconds
    """
    return ((slantRange * 2) / Constants.SPEED_OF_LIGHT)

# -----------------------------------------------------------------------------
# Swath Width
def Swath(NearSlantRange, FarSlantRange):
    """
    Calculates the swath width between a given near and far slant range
    
    NearSlantRange in meters
    FarSlantRange in meters
    
    Returns the swath width in meters
    """
    return FarSlantRange - NearSlantRange

def SwathFromGroundSwath(GroundSwath, IncidenceAngle):
    """
    Calculate the slant range swat width using the ground swath width and 
    incidence angle.
    This is just an approximation.
    
    GroundSwath in meters
    IncidenceAngle in radians
    
    returns the swath width in meters
    """
    return (GroundSwath * np.sin(IncidenceAngle))
    
def SwathGroundRange(Swath, IncidenceAngle):
    """
    Calculates a swath ground range corresponding to a swath width and 
    mid-swath incidence angle. This is just an approximation.
    
    Swath in meters
    IncidenceAngle in radians
    
    returns SwathGroundRange in meters
    """    
    return (Swath / np.sin(IncidenceAngle))

def SceneDuration(Swath):
    """
    Gets the scene duration (tau_Wr) for a given swath width.
    Swath in meters
    
    returns the scene duration in seconds.
    """
    return ((2 * Swath) / Constants.SPEED_OF_LIGHT)

# -----------------------------------------------------------------------------
# Swath Width using Interior Core Angles
def NearAndFarSlantRanges(groundSwath, altitude, midSwathIncidenceAngle, midSwathLookAngle):
    """
    Calculates the near and far slant ranges from the mid-swath incidence 
    and look angles.
    
    groundSwath in meters
    altitude in meters
    midSwathIncidenceAngle in radians
    midSwathLlookAngle in radians
    
    returns BOTH the near and far slant ranges in meters
    """
    RE = Constants.EARTH_RADIUS
    H = altitude

    # get alpha_s
    _as = groundSwath / RE
    # get alpha_m
    _am = midSwathIncidenceAngle - midSwathLookAngle
    # get alpha_n and alpha_f
    _an = _am - _as /2
    _af = _am + _as /2

    # get the near and far slant ranges
    _rn = np.sqrt( np.power((RE + H), 2) + (RE**2) - (2*RE*(RE + H)*np.cos(_an)) )
    _rf = np.sqrt( np.power((RE + H), 2) + (RE**2) - (2*RE*(RE + H)*np.cos(_af)) )
    return _rn, _rf

def CoreAngleFromSlantRange(altitude, slantRange):
    """
    Gets the core angle corresponding to a specific slant range
    and platform altitude.
    
    altitude in meters
    slantRange in meters
    
    returns the core angle in radians
    """
    RE = Constants.EARTH_RADIUS
    H = altitude

    _num = math.pow((RE + H), 2) + math.pow(RE, 2) - math.pow(slantRange, 2)
    _denom = 2 * RE*(RE + H)
    return np.arccos(_num / _denom)

def SlantRageFromCoreAngle(altitude, coreAngle):
    """
    Get the slant range corresponding to a specific core angle
    and platform altitude.
    
    altitude in meters
    coreAngle in radians
    
    returns the slant range in meters
    """
    RE = Constants.EARTH_RADIUS
    H = altitude

    _inner = math.pow(RE, 2) + math.pow((RE + H), 2) - 2*RE*(RE + H)*np.cos(coreAngle)
    return np.sqrt(_inner)

def GroundSwathFromCoreAngles(nearCoreAngle, farCoreAngle):
    """
    Gets the ground swath width from the near core angle and the 
    far core angle.
    Wgr = (alpha_f - alpha_n)Re
    
    nearCoreAngle in radians
    farCoreAngle in radians
    
    returns the ground swath in meters
    """
    return (farCoreAngle - nearCoreAngle) * Constants.EARTH_RADIUS

def CoreAngleDeltaFromGroundSwath(groundSwath):
    """
    Gets the difference between the far and near core angles
    for a given ground swath width.
    (alpha_f - alpha_n) = Wgr/Re
    
    groundSwath in meters
    
    returns the delta in radians
    """
    return groundSwath / Constants.EARTH_RADIUS

# -----------------------------------------------------------------------------
# Ground Range
def SlantRangeFromGroundRange(altitude, groundRange):
    """
    Calculates the slant range for a given altitude and ground range
    
    altitude in meters
    groundRange in meters
    
    returns the slant range in meters
    """
    RE = Constants.EARTH_RADIUS
    H = altitude

    _inner = math.pow(RE, 2) + math.pow((RE + H), 2) - 2*RE*(RE + H) * np.cos(groundRange / RE)
    return np.sqrt(_inner)

# -----------------------------------------------------------------------------
# Core Angles
def CoreAngle(altitude, slant_range):
    """
    Calculate the core angle subtended by the given slant range.
    Make use of a spherical Earth model

    return the core angle in radians
    """ 
    Re = Constants.EARTH_RADIUS
    H = altitude
    Rs = slant_range
    _num = math.pow((Re + H), 2) + math.pow(Re, 2) - np.power(Rs, 2)
    _denom = 2*Re*(Re + H)

    return np.arccos(_num / _denom)

def GroundRangeFromCoreAngle(core_angle):
    """
    Calculate the ground range covered by a given core angle.
    Make use of a spherical Earth model.

    returns the ground range in meters
    """
    return core_angle * Constants.EARTH_RADIUS

def CoreAngleFromGroundRange(ground_range):
    """
    Calculate the core angle subtended by the ground range.

    reutrn the core angle in radians
    """
    return ground_range / Constants.EARTH_RADIUS

def IncidenceAngleToGroundRange(altitude, incidence_angle):
    """
    Calculate the ground range from the incidence angle for a given altitude.
    Make use of a spherical Earth model.

    Return the ground range in meters
    """
    # get the look angle for the incidence angle
    _la = LookAngle(altitude, incidence_angle)
    # get the slant range from incidence angle
    _sr = SlantRange(altitude, incidence_angle, _la)
    # get the core angle for this slant range
    _ca = CoreAngle(altitude, _sr)
    # get the ground range for this core angle
    _gr = GroundRangeFromCoreAngle(_ca)
    return _gr

def GroundRangeToIncidenceAngle(altitude, ground_range):
    """
    Calculat the Incidence angle from a given ground range and altitude.

    return the incidencec angle in radians
    """
    # get the core angle from ground range
    _ca = CoreAngleFromGroundRange(ground_range)
    # get the slant range from core angle
    _sr = SlantRageFromCoreAngle(altitude, _ca)
    # get the incidence angle from slant range
    _ia =  IncidenceAngleFromSlantRange(altitude, _sr)
    return _ia

def SlantRangeToGroundRange(altitude, slant_range):
    """
    Calculate the ground range from a given altitude and slant range.

    return the ground range in meters
    """
    # get the core angle for this slant range
    _ca = CoreAngle(altitude, slant_range)
    # get the ground range for this core angle
    _gr = GroundRangeFromCoreAngle(_ca)
    return _gr

# -----------------------------------------------------------------------------
# Orbital Speed
def PlatformVelocity(Altitude):
    """
    calculate the radar platform velocity in m/s

    Altitude in meters
    """
    RE = Constants.EARTH_RADIUS
    H = Altitude
    MU = Constants.STANDARD_GRAVITATIONAL_CONSTANT

    velocity = math.sqrt(MU / ((RE + H)/1000))*1000
    return velocity

# -----------------------------------------------------------------------------
# Doppler
def DopplerFrequencies(Velocity, wavelength, Angles ):
    """
    Gets the range of Doppler frequencies corresponding to a set of given (azimuth)
    angles, for a specified platform velocity and radar wavelength

    Velocity in m/s
    wavelength in meters
    Angles in radians

    returns the Doppler frequencies in Hertz
    """
    return (2*Velocity*np.sin(Angles)) / wavelength

def AzimuthAngleFromDoppler(DopplerFreq, wavelength, Velocity):
    """
    Get the azimuth angles corresponding to specific Doppler frequencies.

    DopplerFreq in Hertz
    wavelength in meters
    Velocity in m/s

    returns the Azimuth angles in radians
    """
    return np.arcsin((wavelength *DopplerFreq )/(2*Velocity))

def ProcessedDopplerBandwidth(Velocity, AzRes):
    """
    Gets the processed Doppler bandwidth for a given platform velocity and azimuth resolution

    Velocity in m/s
    AzRes in meters

    returns the processed Doppler bandwidth in Hertz
    """
    return Velocity / AzRes

# -----------------------------------------------------------------------------
# Resolution
def GroundrangeResToSlantrangeRes(GroundrangeRes, IncidenceAngle):
    """ 
    Calculate the Slant range resolution for a given ground range and incidencec angle

    GroundrangeRes in meters
    IncidenceAngle in radians

    returns the SLant range resolution in meters
    """
    return GroundrangeRes * np.sin(IncidenceAngle)

def SlantrangeResToGroundrangeRes(SlantrangeRes, IncidenceAngle):
    """ 
    Calculate the Ground range resolution for a given slant range and incidencec angle

    SlantrangeRes in meters
    IncidenceAngle in radians

    returns the Ground range resolution in meters
    """
    return SlantrangeRes / np.sin(IncidenceAngle)