"""rays.py"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, MultipleLocator
from matplotlib.patches import ConnectionPatch, Rectangle
import sys

class Ray:
    """
    Ray Class
    """
    def __init__(self, pos=[0.,0.,0.], direc=[0.,0.,1.]):
        if len(pos) != 3:
            raise Exception("Incorrect position size") 
        if len(direc) != 3:
            raise Exception("Incorrect direction size")
        self.__pos = np.array(pos)
        self.__direc = np.array(direc) / np.linalg.norm(np.array(direc))
        self.__points = [self.__pos.copy()]
    def pos(self):
        """
        Returns the current point of the ray

        return:
            self.__pos: point of ray
        """
        return self.__pos
    def direc(self):
        """
        Returns the current normalised ray direction

        return:
            self.__direc: direction of ray
        """
        return self.__direc
    def append(self, pos=[0,0,0], direc=[0,0,1]):
        """
        Appends a new point and direction to the ray

        Args:
            pos: point of ray
            direc: direction of ray
        
        raises:
            Exception: Incorrect position/direction size
        """
        if len(pos) != 3:
            raise Exception("Incorrect position size")
        if len(direc) != 3:
            raise Exception("Incorrect direction size")
        self.__pos = np.array(pos)
        self.__direc = np.array(direc) / np.linalg.norm(np.array(direc))
        self.__points.append(self.__pos.copy())
    def vertices(self):
        """
        Return all the points along the ray in the form of a list

        return:
            self.__points: points along the ray
        """
        return self.__points

class RayBundle():
    """
    Ray Bundle class

    self.__rays is a list of all the ray objects created in the bundle
    """
    def __init__(self, rmax=5., nrings=5, multi=6):
        self.__rays = []
        for r in range(0, nrings + 1):
            radius = 0
            radius += r * (rmax / nrings)
            n_points = r * multi
            theta = 0
            if n_points == 0:
                x_pos = radius * np.cos(theta)
                y_pos = radius * np.sin(theta)
                self.__rays.append(Ray(pos=[x_pos, y_pos, 0]))
            else:
                for _ in range(0, n_points):
                    x_pos = radius * np.cos(theta)
                    y_pos = radius * np.sin(theta)
                    self.__rays.append(Ray(pos=[x_pos, y_pos, 0]))
                    theta += 2 * (np.pi / n_points)
       
    def propagate_bundle(self, elements=[]):
        """
        Given a list of optical element objects, propagate each ray in the bundle through the object

        Args:
            elements: optical element objects
        """
        for i, ray in enumerate(self.__rays):
            for j, element in enumerate(elements):
                element.propagate_ray(ray)

    def track_plot(self):
        """
        Creates a plot of the path of the rays in the bundle

        return:
            fig: the matplotlib plot
        """
        fig, ax = plt.subplots()
        ray_label = True
        for x in self.__rays:
            vertices = x.vertices()
            x = [v[0] for v in vertices]
            y = [v[1] for v in vertices]
            z = [v[2] for v in vertices]
            if ray_label is True:
                ax.plot(z,y, color="black", label="Rays", linewidth=0.7)
                ray_label = False
            else:
                ax.plot(z,y, color="black", linewidth=0.7)
        ax.legend()

        #Draw Lens
        y = np.linspace(-5, 5, 1000)
        radius = -1/0.03
        z = 100 - np.sqrt(radius**2 - y**2) - radius
        ax.plot(z, y, color='r', label="Lens")
        ax.set_xlabel("z (mm)")
        ax.set_ylabel("y (mm)")

        #Ouput plane
        ax.axvline(x=200, color='blue', linestyle='--', label='Output Plane')
        ax.legend()

        #Zoomed in

        focal_point = 200

        ax_position = 0.65, 0.17, 0.2, 0.2
        ax_zoomed = fig.add_axes(ax_position)
        ax_zoomed.set_xlim(focal_point-0.35, focal_point+0.05)
        ax_zoomed.set_ylim(-0.0004, 0.0004)
        formatter_inset = ScalarFormatter(useOffset=False, useMathText=False)
        ax_zoomed.xaxis.set_major_formatter(formatter_inset)
        ax_zoomed.yaxis.set_major_formatter(formatter_inset)
        ax_zoomed.xaxis.set_major_locator(MultipleLocator(0.1))
        ax_zoomed.tick_params(axis='x', which="both", labelsize=10)
        ax_zoomed.tick_params(axis='y', which="both", labelsize=0, length=0)
        ax_zoomed.yaxis.set_ticks([])

        ray_label = True
        for x in self.__rays:
            vertices = x.vertices()
            x = [v[0] for v in vertices]
            y = [v[1] for v in vertices]
            z = [v[2] for v in vertices]
            if ray_label is True:
                ax_zoomed.plot(z,y, color="black", label="Rays", linewidth=0.5)
                ray_label = False
            else:
                ax_zoomed.plot(z,y, color="black", linewidth=0.5)
            
        
        for spine in ax_zoomed.spines.values():
            spine.set_edgecolor('blue')
            
        
        # Draw a box around the focal point
        box_width = 0.03
        box_height = 0.00005
        rect = Rectangle((focal_point-0.015, -0.000025), box_width, box_height, linewidth=1, edgecolor='b', facecolor='none')
        ax.add_patch(rect)

        # Draw connection lines
        #con1 = ConnectionPatch(xyA=(focal_point, 0), coordsA=ax.transData,
                            #xyB=(focal_point, 0), coordsB=ax_zoomed.transData,
                            #arrowstyle="-", color="black")
        #con2 = ConnectionPatch(xyA=(focal_point-0.03, -0.00005), coordsA=ax.transData,
                            #xyB=(focal_point-0.03, -0.00005), coordsB=ax_zoomed.transData,
                            #arrowstyle="-", color="black")
        con3 = ConnectionPatch(xyA=(focal_point+0.05, 0.0004), coordsA=ax.transData,
                            xyB=(focal_point+0.05, 0.0004), coordsB=ax_zoomed.transData,
                            arrowstyle="-", color="blue")
        con4 = ConnectionPatch(xyA=(focal_point-0.35, 0.0004), coordsA=ax.transData,
                            xyB=(focal_point-0.35, 0.0004), coordsB=ax_zoomed.transData,
                            arrowstyle="-", color="blue")

        for con in [con3, con4]:
            ax.add_patch(con)
        
        ax_zoomed.axvline(x=200, color='blue', linestyle='--', label='Output Plane')


        return fig

    def rms(self):
        """
        Calculates the RMS size of the bundle (spread from optical axis)

        return:
            RMS: root mean square spread from optical axis
        """
        radius = np.array([])
        for i, ray in enumerate(self.__rays):
            position = ray.pos()
            radius = np.append(radius, np.sqrt((position[0])**2 + (position[1])**2))
        return np.sqrt(np.mean(radius**2))

    def spot_plot(self):
        """
        Plot the x-y position of each ray in the bundle at the ray's current positions

        return:
            fig: spot plot figure
        """
        fig, ax = plt.subplots()
        
        ray_label = True
        for i, ray in enumerate(self.__rays):
            point = ray.pos()
            if ray_label is True:
                ax.plot(point[0],point[1],"x",color="black", label="Rays")
                ray_label = False
            else:
                ax.plot(point[0],point[1],"x",color="black", markersize = 4, mew=2)
        ax.legend()
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.set_title("Convex-Plano: z = 198.5 mm")

        return fig
