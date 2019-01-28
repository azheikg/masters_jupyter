#!/usr/bin/python3
"""Planar array antenna Class
"""

__author__ = "Kevin Gema"
__copyright__ = "Copyright 2018, Kevin Gema"
__credits__ = ["Kevin Gema"]
__license__ = "GPL"
__version__ = "1.0.2"
__maintainer__ = ""

import math
import lineararrayantenna as LAA

class PlanarArrayAntenna:
    """Model used for a planar phased array antenna.
    This is a 2-D phased array, producing a gain, antenna patteren and 
    beamwidth for a given set of dimensions."""

    # Default constructor taking in the number of elements in elevation and azimuth,
    # along with the element spacing in elevation and azimuth.
    # The effieciency can be optionally bespecified
    def __init__(self, num_el_elems, num_az_elems, el_elem_spacing, az_elem_spacing, efficiency=0.886):
        self.efficiency      = efficiency
        self.num_el_elem     = num_el_elems
        self.num_az_elem     = num_az_elems
        self.el_elem_spacing = el_elem_spacing
        self.az_elem_spacing = az_elem_spacing

        # use a linear array antenna for elevation
        self.elevation = LAA.LinearArrayAntenna(num_el_elems, el_elem_spacing, efficiency=1)

        # use linear array antenna for azimuth
        self.azimuth = LAA.LinearArrayAntenna(num_az_elems, az_elem_spacing, efficiency=1)

    @staticmethod
    def panel_to_antenna(panel, num_el_panels, num_az_panels, efficiency):
        efficiency      = efficiency
        num_el_elem     = panel.num_el_elem * num_el_panels
        num_az_elem     = panel.num_az_elem * num_az_panels
        el_elem_spacing = panel.el_elem_spacing
        az_elem_spacing = panel.az_elem_spacing
        return num_el_elem, num_az_elem, el_elem_spacing, az_elem_spacing, efficiency

    @classmethod
    def by_num_panels(cls, panel, num_el_panels, num_az_panels, efficiency=0.886):
        # instantiate the planar array antenna as a combinations of panels (planar arrays)
        return cls(*cls.panel_to_antenna(panel, num_el_panels, num_az_panels, efficiency))

    def length(self):
        """Get the planar array length in azimuth"""
        return (self.num_az_elem * self.az_elem_spacing)

    def height(self):
        """Get the planar array length in elevation"""
        return (self.num_el_elem * self.el_elem_spacing)

    def num_elevation_elements(self):
        """Returns the number of elements in elevation"""
        return self.elevation.num_elements();

    def num_azimuth_elements(self):
        """Returns the number of elements in azimuth"""
        return self.azimuth.num_elements();

    def elevation_element_spacing(self):
        """Returns the element spacing in elevation"""
        return self.elevation.element_spacing();

    def azimuth_element_spacing(self):
        """Returns the element spacing in azimuth"""
        return self.azimuth.element_spacing();

    def elevation_beamwidth(self, wavelength):
        """Get the planar array elevation beamwidth in radians"""
        return (self.efficiency * self.elevation.beamwidth(wavelength))

    def azimuth_beamwidth(self, wavelength):
        """Get the planar array azimuth beamwidth in radians"""
        return (self.efficiency * self.azimuth.beamwidth(wavelength))

    def gain(self, wavelength):
        """Get the planar array antenna gain"""
        area = self.length() * self.height()
        return ((4*math.pi * area * self.efficiency) / (math.pow(wavelength, 2)))

    def gain_alt(self, wavelength):
        """Get the planar array gain using 2 lieanr array approximations for the dimensions"""
        gain1 = self.elevation.gain()
        gain2 = self.azimuth.gain()
        return (self.efficiency * gain1 * gain2)

    def is_grating_lobes_visible(self, wavelength, el_steeringAngle, az_steeringAngle):
        """Determines if a grating lobe occurs within the visible range, using 
        wavelength and steering angle.        
        returns True/False for both elevation and azimuth
        """
        gl_el = self.elevation.is_grating_lobes_visible(wavelength, el_steeringAngle) 
        gl_az = self.azimuth.is_grating_lobes_visible(wavelength, az_steeringAngle) 
        return gl_el, gl_az

    def grating_lobe_angles(self, wavelength, el_steeringAngle, az_steeringAngle, numGratingLobes=1):
        """Gets the positions of the grating lobes.        
        Returns an array with the positions of the grating lobes in radians
        """
        el_lobes = self.elevation.grating_lobe_angles(wavelength, el_steeringAngle, numGratingLobes)
        az_lobes = self.azimuth.grating_lobe_angles(wavelength, az_steeringAngle, numGratingLobes)
        return el_lobes, az_lobes

    def elevation_field_pattern(self, wavelength, angle_range, steering_offset, elementFactor=1, normalised=True):
        """Calculates the elevation field pattern for a given wavelength over
        a specialised angle range. The pattern is steered to the specified offset."""
        return self.elevation.field_pattern(wavelength, angle_range, steering_offset, elementFactor, normalised)

    def azimuth_field_pattern(self, wavelength, angle_range, steering_offset, elementFactor=1, normalised=True):
        """Calculates the azimuth field pattern for a given wavelength over
        a specialised angle range. The pattern is steered to the specified offset."""
        return self.azimuth.field_pattern(wavelength, angle_range, steering_offset, elementFactor, normalised)

    def elevation_power_pattern(self, wavelength, angle_range, steering_offset, elementFactor=1, normalised=True):
        """Calculates the elevation power pattern for a given wavelength over
        a specialised angle range. The pattern is steered to the specified offset."""
        return self.elevation.power_pattern(wavelength, angle_range, steering_offset, elementFactor, normalised)

    def azimuth_power_pattern(self, wavelength, angle_range, steering_offset, elementFactor=1, normalised=True):
        """Calculates the azimuth power pattern for a given wavelength over
        a specialised angle range. The pattern is steered to the specified offset."""
        return self.azimuth.power_pattern(wavelength, angle_range, steering_offset, elementFactor, normalised)

    # Static method to calculate the number of elements required to achieve a specified beamwidth.
    @staticmethod
    def num_elements_for_beamwidth(wavelength, beamwidth, elem_spacing, efficiency=0.886):
        return int(math.ceil((efficiency * wavelength) / (elem_spacing * beamwidth)))
        #return int(math.floor((efficiency * wavelength) / (elem_spacing * beamwidth)))