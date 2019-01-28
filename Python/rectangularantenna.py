#!/usr/bin/python3
"""rectangular antenna Class
"""

__author__ = "Kevin Gema"
__copyright__ = "Copyright 2018, Kevin Gema"
__credits__ = ["Kevin Gema"]
__license__ = "GPL"
__version__ = "1.0.2"
__maintainer__ = ""

import math
import numpy as np

class RectangularAntenna:
    """Model used for a Rectangular Aperture Antenna. A single aperture
    is defined, producing a gain, antenna patterns in elevation and 
    azimuth and a beamwidth for a given set of dimensions."""
    def __init__(self, length, height, efficiency=0.886):
        self.az_length  = length
        self.el_height  = height
        self.efficiency = efficiency

    def length(self):
        """Get the rectangular aperture length"""
        return self.az_length

    def height(self):
        """Get the rectangular aperture height"""
        return self.el_height 

    def elevation_beamwidth(self, wavelength):
        """Returns the elevation half-power beamwidth for a given wavelength."""
        return (self.efficiency * wavelength) / self.height()

    def azimuth_beamwidth(self, wavelength):
        """Returns the azimuth half-power beamwidth for a given wavelength."""
        return (self.efficiency * wavelength) / self.length()

    def gain(self, wavelength):
        """Calculates the gain for the antenna using the specified wavelength"""
        area = self.length() * self.height()
        return ((4*math.pi * area * self.efficiency) / (math.pow(wavelength, 2)))

    def elevation_field_pattern(self, wavelength, angle_range, steering_offset=0):
        """Calculates the elevation field pattern for a specified wavelength, 
        over a specified range of angles. The antenna pattern is offset at
        the specified steering offset"""
        angles = np.sin(angle_range - steering_offset)
        pat = np.sinc(angles * self.el_height / wavelength)
        return pat

    def azimuth_field_pattern(self, wavelength, angle_range, steering_offset=0):
        """Calculates the azimuth field pattern for a specified wavelength, 
        over a specified range of angles. The antenna pattern is offset at
        the specified steering offset"""
        angles = np.sin(angle_range - steering_offset)
        pat = np.sinc(angles * self.az_length / wavelength)
        return pat

    def elevation_power_pattern(self, wavelength, angle_range, steering_offset=0):
        """Calculates the elevation power pattern for a aspecified wavevlength, 
        over a specified range of angles. The antenna pattern is offset at
        the specified steering offset."""
        field_pat = self.elevation_field_pattern(wavelength, angle_range, steering_offset)
        return np.power(np.abs(field_pat), 2)

    def azimuth_power_pattern(self, wavelength, angle_range, steering_offset=0):
        """Calculates the azimuth power pattern for a aspecified wavevlength, 
        over a specified range of angles. The antenna pattern is offset at
        the specified steering offset."""
        field_pat = self.azimuth_field_pattern(wavelength, angle_range, steering_offset)
        return np.power(np.abs(field_pat), 2)

    # static method to calculate the length required for a specified beamwidth
    @staticmethod
    def dimension_for_beamwidth(wavelength, beamwidth, efficiency=0.886):
        return (efficiency * wavelength) / beamwidth

    # static method to calculate the beamwidth from a specified dimension
    @staticmethod
    def beamwidth_from_dimension(wavelength, dimension, efficiency):
        return (efficiency * wavelength) / dimension