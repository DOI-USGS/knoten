from knoten import csm, utils
import spiceypy as spice
import numpy as np
import csmapi

class Ellipsoid:
    """
    A biaxial ellipsoid shape model.
    """

    def __init__(self, semi_major, semi_minor = None):
        """
        Create an ellipsoid shape model from a set of radii

        Parameters
        ----------
        semi_major : float
                     The equatorial semi-major radius of the ellipsoid.
        semi_minor : float
                     The polar semi-minor radius of the ellipsoid.
        """
        self.a = semi_major
        self.b = semi_major
        self.c = semi_major

        if semi_minor is not None:
            self.c = semi_minor

    @classmethod
    def from_csm_sensor(cls, sensor):
        semi_major, semi_minor = csm.get_radii(sensor)
        return cls(semi_major, semi_minor)
    
    def get_surface_normal(self, ground_pt):
        """
        Given a ground point, calculate the surface normal.

        Parameters
        ----------
        ground_pt: tuple
                The ground point as an (x, y, z) tuple

        Returns
        -------
        : tuple
        in the form (x, y, z)
        """
        normal = spice.surfnm(self.a, self.b, self.c, np.array([ground_pt.x, ground_pt.y, ground_pt.z]))
        return utils.Point(normal[0], normal[1], normal[2])


    def intersect_surface(self, sensor_pos, look_vec):
        sensor_pos = np.array([sensor_pos.x, sensor_pos.y, sensor_pos.z])
        look_vec = np.array([look_vec.x, look_vec.y, look_vec.z])
        
        ground_pt = spice.surfpt(sensor_pos, look_vec, self.a, self.b, self.c)
        return csmapi.EcefCoord(ground_pt[0], ground_pt[1], ground_pt[2])
        
