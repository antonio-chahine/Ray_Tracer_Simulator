"""lenses.py"""
import sys
from raytracer.elements import SphericalRefraction, PlanarSurface
from raytracer.rays import Ray

class PlanoConvex(SphericalRefraction):
    """
    Plano Convex class
    """
    def __init__(self, z_0=100, curvature=0., curvature1=0., \
                 curvature2=0., n_inside=0., n_outside=0., thickness=0., aperture=0.):
        self.__z_0 = z_0
        self.__curvature = curvature
        self.__curvature1 = curvature1
        self.__curvature2 = curvature2
        self.__n_inside = n_inside
        self.__n_outside = n_outside
        self.__thickness = thickness
        self.__aperture = aperture

    def elements(self):
        """
        Returns a list of the planar and surface elements specified by the inputs

        return:
            pc: plano-convex element list
            cp: convex-planar element list
        """
        if self.__curvature == 0:
            if self.__curvature1 == 0:
                self.__curvature = self.__curvature2
                return[PlanarSurface(z_0=self.__z_0, \
                            aperture=self.__aperture, n_1=self.__n_outside, n_2=self.__n_inside), \
                       SphericalRefraction(z_0=(self.__z_0 + self.__thickness), \
                            aperture=self.__aperture, curvature=self.__curvature, \
                            n_1=self.__n_inside, n_2=self.__n_outside)]
            else:
                self.__curvature = self.__curvature1
                return[SphericalRefraction(z_0=self.__z_0, aperture=self.__aperture, \
                            curvature=self.__curvature, n_1=self.__n_outside, n_2=self.__n_inside),\
                       PlanarSurface(z_0=(self.__z_0+self.__thickness), \
                            aperture=self.__aperture, n_1=self.__n_inside, n_2=self.__n_outside)]
        else:
            if self.__curvature < 0:
                return[PlanarSurface(z_0=self.__z_0, aperture=self.__aperture, \
                            n_1=self.__n_outside, n_2=self.__n_inside), \
                       SphericalRefraction(z_0=(self.__z_0 + self.__thickness), \
                            aperture=self.__aperture, curvature=self.__curvature, \
                            n_1=self.__n_inside, n_2=self.__n_outside)] 
            else:
                return[SphericalRefraction(z_0=self.__z_0, aperture=self.__aperture, \
                            curvature=self.__curvature, n_1=self.__n_outside, n_2=self.__n_inside),\
                       PlanarSurface(z_0=(self.__z_0+self.__thickness), aperture=self.__aperture, \
                            n_1=self.__n_inside, n_2=self.__n_outside)]

    def focal_point(self):
        """
        Returns the focal point given a list of optical elements

        return:
            focus[2]: focal point (z axis)
        """
        ray = Ray(pos=[sys.float_info.min, 0, 0])
        for i, element in enumerate(self.elements()):
            element.propagate_ray(ray)
        initial_pos = ray.pos()
        initial_direc = ray.direc()

        distance_vector = (0 - initial_pos[0]) / initial_direc[0]
        focus = [initial_pos[0] + distance_vector * initial_direc[0], \
                                    initial_pos[1] + distance_vector * initial_direc[1], \
                                    initial_pos[2] + distance_vector * initial_direc[2]] 
        return focus[2]
class ConvexPlano(SphericalRefraction):
    """
    Convex-Plano class
    """
    def __init__(self, z_0=100, curvature=0., curvature1=0., \
                 curvature2=0., n_inside=0., n_outside=0., thickness=0., aperture=0.):
        self.__z_0 = z_0
        self.__curvature = curvature
        self.__curvature1 = curvature1
        self.__curvature2 = curvature2
        self.__n_inside = n_inside
        self.__n_outside = n_outside
        self.__thickness = thickness
        self.__aperture = aperture

    def elements(self):
        """
        Returns a list of the planar and surface elements specified by the inputs

        return:
            pc: plano-convex element list
            cp: convex-planar element list
        """
        if self.__curvature == 0:
            if self.__curvature1 == 0:
                self.__curvature = self.__curvature2
                return[PlanarSurface(z_0=self.__z_0, \
                            aperture=self.__aperture, n_1=self.__n_outside, n_2=self.__n_inside), \
                       SphericalRefraction(z_0=(self.__z_0 + self.__thickness), \
                            aperture=self.__aperture, curvature=self.__curvature, \
                            n_1=self.__n_inside, n_2=self.__n_outside)]
            else:
                self.__curvature = self.__curvature1
                return[SphericalRefraction(z_0=self.__z_0, aperture=self.__aperture, \
                            curvature=self.__curvature, n_1=self.__n_outside, n_2=self.__n_inside),\
                       PlanarSurface(z_0=(self.__z_0+self.__thickness), \
                            aperture=self.__aperture, n_1=self.__n_inside, n_2=self.__n_outside)]
        else:
            if self.__curvature < 0:
                return[PlanarSurface(z_0=self.__z_0, aperture=self.__aperture, \
                            n_1=self.__n_outside, n_2=self.__n_inside), \
                       SphericalRefraction(z_0=(self.__z_0 + self.__thickness), \
                            aperture=self.__aperture, curvature=self.__curvature, \
                            n_1=self.__n_inside, n_2=self.__n_outside)] 
            else:
                return[SphericalRefraction(z_0=self.__z_0, aperture=self.__aperture, \
                            curvature=self.__curvature, n_1=self.__n_outside, n_2=self.__n_inside),\
                       PlanarSurface(z_0=(self.__z_0+self.__thickness), aperture=self.__aperture, \
                            n_1=self.__n_inside, n_2=self.__n_outside)]
