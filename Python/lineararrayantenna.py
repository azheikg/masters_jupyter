#!/usr/bin/python3
"""Linear array antenna Class
"""

__author__ = "Kevin Gema"
__copyright__ = "Copyright 2018, Kevin Gema"
__credits__ = ["Kevin Gema"]
__license__ = "GPL"
__version__ = "1.0.2"
__maintainer__ = ""

import math
import numpy as np

class LinearArrayAntenna:
    """Model used for a Linear Phased Array Antenna.
    This is a 1-D phased array, producing a gain, antenna pattern and
    beamwidth for a given set of dimensions"""
    def __init__(self, num_elements=1, element_spacing=1, efficiency=0.886):
        self.efficiency   = efficiency
        self.num_elem     = num_elements
        self.elem_spacing = element_spacing

    def length(self):
        """Get the linear array antenna length"""
        return self.elem_spacing * self.num_elem

    def beamwidth(self, wavelength):
        """Get the linear array beamwidth in radians"""
        return (self.efficiency * wavelength) / self.length();

    def num_elements(self):
        """Returns the number of elements"""
        return self.num_elem

    def element_spacing(self):
        """Returns the element spacing"""
        return self.elem_spacing

    def directivity(self, wavelength):
        """Gets the directivity for the linear phased array."""
        _num = math.pow(self.num_elem, 2)
        _denomSum = 0
        for n in range(1, self.num_elem):
            _innerTop = (self.num_elem - n) * np.sin(n * (2*np.pi / wavelength) * self.elem_spacing)
            _innerBottom = n * (2*np.pi / wavelength) * self.elem_spacing
            _denomSum = _denomSum + (_innerTop / _innerBottom)
        _denom = self.num_elem + (2 * _denomSum)
        
        return _num / _denom

    def gain(self, wavelength):
        """Calculates the phased array antenna gain."""
        # check which model to use
        return (self.efficiency * self.directivity(wavelength))

    def field_pattern(self, wavelength, angle_range, steering_offset, elementFactor=1, normalised=True):
        """Get the antenna pattern for a given wavelength, angle range and steering angle offset. 
        A element factor can also be specified, by default the elements are considered to be uniform 
        isotropic radiators. The pattern can be specified to be normalised or un-normalised."""
        _numInner = self.num_elem * (np.pi * self.elem_spacing / wavelength) * (np.sin(angle_range) - np.sin(steering_offset))
        _num = np.sin(_numInner)
        _denomInner = (np.pi * self.elem_spacing / wavelength) * (np.sin(angle_range) - np.sin(steering_offset))
        
        if normalised == True:
            _denom = self.num_elem * np.sin(_denomInner)
        else:
            _denom = np.sin(_denomInner)

        pat = _num / _denom
        return pat

    def field_pattern_alt(self, wavelength, angle_range, steering_offset, elementFactor=1, normalised=True):
        """Get the antenna pattern for a given wavelength, angle range and steering angle offset. 
        A element factor can also be specified, by default the elements are considered to be uniform 
        isotropic radiators. The pattern can be specified to be normalised or un-normalised."""
        _numInner = self.num_elem * (np.pi * self.elem_spacing / wavelength) * (np.cos(angle_range) - np.cos(steering_offset))
        _num = np.sin(_numInner)
        _denomInner = (np.pi * self.elem_spacing / wavelength) * (np.cos(angle_range) - np.cos(steering_offset))
        
        if normalised == True:
            _denom = self.num_elem * np.sin(_denomInner)
        else:
            _denom = np.sin(_denomInner)

        pat = _num / _denom
        return pat

    def power_pattern(self, wavelength, angle_range, steering_offset=0, elementFactor=1, normalised=True):
        """Calculates the power pattern for a aspecified wavevlength, over a specified range of 
        angles. The antenna pattern is offset at the specified steering offset."""
        field_pat = self.field_pattern(wavelength, angle_range, steering_offset, elementFactor, normalised)
        return np.power(np.abs(field_pat), 2)

    def power_pattern_alt(self, wavelength, angle_range, steering_offset=0, elementFactor=1, normalised=True):
        """Calculates the power pattern for a aspecified wavevlength, over a specified range of 
        angles. The antenna pattern is offset at the specified steering offset."""
        field_pat = self.field_pattern_alt(wavelength, angle_range, steering_offset, elementFactor, normalised)
        return np.power(np.abs(field_pat), 2)

    def is_grating_lobes_visible(self, wavelength, steeringAngle):
        """Determines if a grating lobe occurs within the visible range, using 
        wavelength, element spacing and steering angle
        
        Returns True or False
        """
        lhs = self.elem_spacing
        #rhs = wavelength / (1 + np.absolute(np.sin(steeringAngle)))
        rhs = wavelength / (1 + np.absolute(np.cos(steeringAngle)))
        
        test = lhs <= rhs
        return not test

    def grating_lobe_angles(self, wavelength, steeringAngle, numGratingLobes=1):
        """
        Gets the positions of the grating lobes in radians"""
        globes = []
        for r in range(1, numGratingLobes+1):
            gl_pos = np.arcsin(((wavelength*r)/self.elem_spacing) + np.sin(steeringAngle))
            gl_neg = np.arcsin( -((wavelength*r)/self.elem_spacing) + np.sin(steeringAngle))
            globes.append(gl_neg)
            globes.append(gl_pos)
            
        return np.asarray(globes)

    @staticmethod
    def num_elements_for_beamwidth(wavelength, beamwidth, elem_spacing, efficiency=0.886):
        return int(math.floor((efficiency * wavelength) / (elem_spacing * beamwidth)))