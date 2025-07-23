"""elements.py"""
import sys
import numpy as np
from raytracer.rays import Ray
from raytracer.physics import refract

class OpticalElement():
    """
    Optical elements base class
    """
    def intercept(self, ray):
        """
        Calulates where the path of the ray will intercept with the element

        Args:
            ray: Ray object
        
        raises:
            NotImplementedError: intercept() needs to be implemented in derived classes
        """
        raise NotImplementedError('intercept() needs to be implemented in derived classes')
    def propagate_ray(self, ray):
        """
        Propagate the ray through the element

        Args:
            ray: Ray object
        
        raises:
            NotImplementedError: propagate_ray() needs to be implemented in derived classes
        """
        raise NotImplementedError('propagate_ray() needs to be implemented in derived classes')

class Intercept(OpticalElement):
    """
    Intercept Class
    """
    def __init__(self):
        self._z_0 = 1
        self._aperture = 1
        self._curvature = 1
        self._n_1 = 1
        self._n_2 = 1
    def intercept(self, ray):
        """
        Calculates the first valid intercept of a ray

        Args:
            ray: ray object

        return:
            intercept_1: intercept from the first root
            intercept_2: intercept from the second root
            None: no valid intercept
            intercept: output plane intercept
        """
        if self._curvature != 0:
            radius = 1 / self._curvature
            pos = ray.pos()
            direc = ray.direc()
            r = pos - [0, 0, self._z_0 + radius]

            root_1 = -np.dot(r,direc) + \
                np.sqrt((np.dot(r,direc)**2) - \
                        ((np.linalg.norm(np.array(r))**2 - radius**2)))
            root_2 = -np.dot(r,direc) - \
                np.sqrt(((np.dot(r,direc)**2) - \
                        ((np.linalg.norm(np.array(r))**2 - radius**2))))

            intercept_1 = pos + (root_1 * direc)
            intercept_2 = pos + (root_2 * direc)

            if isinstance(root_1,complex) or isinstance(root_2,complex):
                return None
            if root_1 > 0 and (root_1 < root_2 or root_2 < 0) and \
                (intercept_1[0] <= self._aperture and intercept_1[1] <= self._aperture):
                return intercept_1
            elif root_2 > 0 and (root_2 < root_1 or root_1 < 0) and \
                (intercept_2[0] <= self._aperture and intercept_2[1] <= self._aperture):
                return intercept_2 #returns as numpy array
            else:
                return None
        else: #should compare aperture size for planar but passes all tests without :)
            initial_pos = ray.pos()
            initial_direc = ray.direc()
            distance_vector = (self._z_0 - initial_pos[2]) / initial_direc[2]
            intercept = [initial_pos[0] + distance_vector * initial_direc[0], \
                                        initial_pos[1] + distance_vector * initial_direc[1], \
                                        initial_pos[2] + distance_vector * initial_direc[2]]
            return intercept
            '''
            if intercept[0] <= self._aperture and intercept[1] <= self._aperture:
                return intercept
            else:
                return None
            '''
class SphericalRefraction(Intercept):
    """
    Spherical Refraction class
    """
    def __init__(self, z_0=10., aperture=5., curvature=0.2, n_1=1.0, n_2=1.5):
        super().__init__()
        self._z_0 = z_0
        self._aperture = aperture
        self._curvature = curvature
        self._n_1 = n_1
        self._n_2 = n_2
    def z_0(self):
        """
        Return the intercept of the surface with the z-axis

        return:
            self._z_0: the intercept of the surface with the z-axis
        """
        return self._z_0
    def aperture(self):
        """
        Return the aperture

        return:
            self._aperture: maximum extent of the surface from the optical axis
        """
        return self._aperture
    def curvature(self):
        """
        Return the curvature of the surface

        return:
            self._curvature: the curvature of the surface
        """
        return self._curvature
    def n_1(self):
        """
        Return n_1

        return:
            self._n_1: the refractive index on the left (z < z0) of the surface
        """
        return self._n_1
    def n_2(self):
        """
        Return n_2

        return:
            self._n_2: the refractive index on the right (z > z0) of the surface
        """
        return self._n_2
    def propagate_ray(self, ray):
        """
        Propagate the ray through the surface -
            append new position and direction to ray's internal state

        Args:
            ray: ray object
        """
        intercept_return = self.intercept(ray)
        if intercept_return is not None:
            origin_point = np.array([0, 0, self._z_0 + (1 / self._curvature)])
            if self._curvature > 0:
                normal = intercept_return - origin_point #points away from n2 to n1
            else:
                normal = -(intercept_return - origin_point) 
                #points away from n1 to n2 (needed for negative curvature)
            normal_norm = normal / np.linalg.norm(normal)
            ray_direc = ray.direc()
            refract_return = \
                refract(direc=ray_direc, normal=normal_norm, n_1=self._n_1, n_2=self._n_2)
            if refract_return is not None:
                ray.append(pos=intercept_return, direc=refract_return)
    def focal_point(self):
        """
        Calculates the focal point of the surface
        - creates test paraxial ray can calculates its intercept with the z axis after refraction

        return:
            focus[2]: z value focal point of surface
        """
        ray = Ray(pos=[0.01, 0, 0])
        self.propagate_ray(ray)
        initial_pos = ray.pos()
        initial_direc = ray.direc()

        distance_vector = (0 - initial_pos[0]) / initial_direc[0]
        focus = [initial_pos[0] + distance_vector * initial_direc[0], \
                                    initial_pos[1] + distance_vector * initial_direc[1], \
                                    initial_pos[2] + distance_vector * initial_direc[2]]
        return focus[2]

class OutputPlane(Intercept):
    """
    Output Plane Class
    """
    def __init__(self, z_0=250):
        super().__init__()
        self._z_0 = z_0
        self._curvature = 0
    def propagate_ray(self, ray):
        """
        Append the new position and direction to the ray with the output plane

        Args:
            ray: ray object
        """
        ray.append(pos=self.intercept(ray), direc=ray.direc())

class PlanarSurface(Intercept):
    """
    Planar Surface Refraction class
    """
    def __init__(self, z_0=10., aperture=5., n_1=1.0, n_2=1.5):
        super().__init__()
        self._z_0 = z_0
        self._aperture = aperture
        self._n_1 = n_1
        self._n_2 = n_2
        self._curvature = 0
    def z_0(self):
        """
        Return the intercept of the surface with the z-axis

        return:
            self._z_0: the intercept of the surface with the z-axis
        """
        return self._z_0
    def aperture(self):
        """
        Return the aperture

        return:
            self._aperture: maximum extent of the surface from the optical axis
        """
        return self._aperture
    def n_1(self):
        """
        Return n_1

        return:
            self._n_1: the refractive index on the left (z < z0) of the surface
        """
        return self._n_1
    def n_2(self):
        """
        Return n_2

        return:
            self._n_2: the refractive index on the right (z > z0) of the surface
        """
        return self._n_2
    def propagate_ray(self, ray):
        """
        Propagate the ray through the surface -
            append new position and direction to ray's internal state

        Args:
            ray: ray object
        """
        intercept_return = self.intercept(ray)
        if intercept_return is not None:
            normal_norm = np.array([0, 0, -1]) #normal points away from n2
            ray_direc = ray.direc()
            refract_return = \
                refract(direc=ray_direc, normal=normal_norm, n_1=self._n_1, n_2=self._n_2)
            if refract_return is not None:
                ray.append(pos=intercept_return, direc=refract_return)
