#!/usr/bin/python3
"""Circular antenna Class
"""

__author__ = "Kevin Gema"
__copyright__ = "Copyright 2018, Kevin Gema"
__credits__ = ["Kevin Gema"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = ""

import math
import numpy as np
import scipy.special as sp

class CircularAntenna:
    """Model used for a circular aperture antenna.
    A single aperture is defined producing a gain, antenna patterns and 
    beamwidth"""
    def __init__(self, diameter=0, beamwidth=0, efficiency=1):
        self.diam = diameter
        self.bw = beamwidth
        self.efficiency = efficiency

    def set_efficiency(self, efficiency):
        """Sets the antenna efficiency."""
        self.efficiency = efficiency

    def gain(self, wavelength, aperture_efficiency=1):
        """Calculates the antenna gain for a given wavelength and
        aperture efficiency."""
        area = 0.25 * math.pi * math.pow(self.diam, 2)
        return ((4*math.pi * area * aperture_efficiency) / (math.pow(wavelength, 2)))

    def set_diameter(self, diameter):
        """Sets the diameter for the circular aperture."""
        self.diam = diameter

    def beamwidth(self, wavelength):
        """Calculate the beamwidth for the circular aperture are a given 
        wavelength."""
        if self.bw != 0:
            return self.bw
        elif self.diam != 0:
            self.bw = self.efficiency * (wavelength / self.diam)
            return self.bw
        else:
            return 0

    def diameter(self, wavelength):
        """Returns the aperture diameter.
        If the diameter isn't specified, calculate it using the specified 
        wavelength"""
        if self.diam != 0:
            return self.diam
        elif self.bw != 0:
            self.diam = self.efficiency * (wavelength / self.bw)
            return self.diam
        else:
            return 0

    def field_pattern(self, wavelength, angle_range, steering_offset=0):
        """Calculates the antenna elevation field pattern for a given 
        wavelength, over the specified angle range.
        The antenna pattern is steered to the specified offset."""
        u = (self.diam / wavelength) * np.sin(angle_range)
        pat = 2 * sp.j1(2*np.pi*u) / (2*np.pi*u)
        return pat

    def power_pattern(self, wavelength, angle_range, steering_offset=0):
        """Calculates the antenna elevation power pattern for a given 
        wavelength, over the specified angle range.
        The antenna pattern is steered to the specified offset."""
        field_pat = self.field_pattern(wavelength, angle_range, steering_offset)
        return np.power(np.abs(field_pat), 2)