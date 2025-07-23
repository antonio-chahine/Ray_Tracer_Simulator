"""Analysis module."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, MultipleLocator
from matplotlib.patches import ConnectionPatch, Rectangle
from raytracer.rays import Ray, RayBundle
from raytracer.elements import SphericalRefraction, OutputPlane
from raytracer.lenses import PlanoConvex


params = {
   'axes.labelsize': 18,
   'font.size': 18,
   'font.family': 'sans-serif', 
   'font.serif': 'Arial', 
   'legend.fontsize': 12,
   'xtick.labelsize': 16,
   'ytick.labelsize': 16, 
   'figure.figsize': [8.8, 8.8/1.618] 
} 
plt.rcParams.update(params)

def task8():
    """
    Task 8.

    In this function you should check your propagate_ray function properly
    finds the correct intercept and correctly refracts a ray. Don't forget
    to check that the correct values are appended to your Ray object.
    """
    for i in range(0,3):
        ray = Ray([i, 2., 0.], [0., 0., 1.])
        sr = SphericalRefraction(z_0=10, curvature=0.02, n_1=1., n_2=1.5, aperture=50.)
        sr.propagate_ray(ray)
        print("Points: ", ray.vertices())
        print("Direction: ", ray.direc())

def task10():
    """
    Task 10.

    In this function you should create Ray objects with the given initial positions.
    These rays should be propagated through the surface, up to the output plane.
    You should then plot the tracks of these rays.
    This function should return the matplotlib figure of the ray paths.

    Returns:
        Figure: the ray path plot.
    """
    ray = {
        "ray1": Ray(pos=[0, 4, 0]),
        "ray2": Ray(pos=[0, 1, 0]),
        "ray3": Ray(pos=[0, 0.2, 0]),
        "ray4": Ray(pos=[0, 0, 0]),
        "ray5": Ray(pos=[0, -0.2, 0]),
        "ray6": Ray(pos=[0, -1, 0]),
        "ray7": Ray(pos=[0, -4, 0]),
    }
    sr = SphericalRefraction(z_0=100, curvature=0.03, n_1=1., n_2=1.5, aperture=5.)
    for x in ray.values():
        sr.propagate_ray(x)

    sr_output_plane = OutputPlane(z_0=250)
    for x in ray.values():
        sr_output_plane.propagate_ray(x)
    
    ray_label = True
    fig, ax = plt.subplots()
    for x in ray.values():
        vertices = x.vertices()
        x = [v[0] for v in vertices]
        y = [v[1] for v in vertices]
        z = [v[2] for v in vertices]
        if ray_label is True:
            ax.plot(z,y, color="black", label="Rays")
            ray_label = False
        else:
            ax.plot(z,y, color="black")

    y = np.linspace(-5, 5, 1000)
    radius = -1/0.03
    z = 100 - np.sqrt(radius**2 - y**2) - radius
    ax.plot(z, y, color='r', label="Lens")
    ax.set_xlabel("z (mm)")
    ax.set_ylabel("y (mm)")

    #Zoomed in

    focal_point = sr.focal_point()

    ax_position = 0.62, 0.17, 0.2, 0.2
    ax_zoomed = fig.add_axes(ax_position)
    ax_zoomed.set_xlim(focal_point-0.03, focal_point+0.005)
    ax_zoomed.set_ylim(-0.00005, 0.00005)
    formatter_inset = ScalarFormatter(useOffset=False, useMathText=False)
    ax_zoomed.xaxis.set_major_formatter(formatter_inset)
    ax_zoomed.yaxis.set_major_formatter(formatter_inset)
    ax_zoomed.xaxis.set_major_locator(MultipleLocator(0.02))
    ax_zoomed.tick_params(axis='x', which="both", labelsize=10)
    ax_zoomed.tick_params(axis='y', which="both", labelsize=0, length=0)
    ax_zoomed.yaxis.set_ticks([])

    ray_label = True
    for x in ray.values():
        vertices = x.vertices()
        x = [v[0] for v in vertices]
        y = [v[1] for v in vertices]
        z = [v[2] for v in vertices]
        if ray_label is True:
            ax_zoomed.plot(z,y, color="black", label="Rays")
            ray_label = False
        else:
            ax_zoomed.plot(z,y, color="black")
        
    
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
    con3 = ConnectionPatch(xyA=(focal_point+0.005, 0.00005), coordsA=ax.transData,
                        xyB=(focal_point+0.005, 0.00005), coordsB=ax_zoomed.transData,
                        arrowstyle="-", color="blue")
    con4 = ConnectionPatch(xyA=(focal_point-0.03, 0.00005), coordsA=ax.transData,
                        xyB=(focal_point-0.03, 0.00005), coordsB=ax_zoomed.transData,
                        arrowstyle="-", color="blue")

    for con in [con3, con4]:
        ax.add_patch(con)

    #Ouput plane
    ax.axvline(x=250, color='blue', linestyle='--', label='Output Plane')
    ax.legend()
    
    return fig


def task11():
    """
    Task 11.

    In this function you should propagate the three given paraxial rays through the system
    to the output plane and the tracks of these rays should then be plotted.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for ray paths
    2. the calculated focal point.

    Returns:
        tuple[Figure, float]: the ray path plot and the focal point
    """
    ray = {
        "ray1": Ray(pos=[0.1, 0.1, 0]),
        "ray2": Ray(pos=[0, 0, 0]),
        "ray3": Ray(pos=[-0.1, -0.1, 0]),
    }
    sr = SphericalRefraction(z_0=100, curvature=0.03, n_1=1., n_2=1.5, aperture=5.)
    for x in ray.values():
        sr.propagate_ray(x)
    
    focal_point = sr.focal_point()

    sr_output_plane = OutputPlane(z_0=focal_point)
    for x in ray.values():
        sr_output_plane.propagate_ray(x)

    ray_label = True
    fig, ax = plt.subplots()
    for x in ray.values():
        vertices = x.vertices()
        x = [v[0] for v in vertices]
        y = [v[1] for v in vertices]
        z = [v[2] for v in vertices]
        if ray_label is True:
            ax.plot(z,y, color="black", label="Rays")
            ray_label = False
        else:
            ax.plot(z,y, color="black")
    
    # DRAW LENS
    y = np.linspace(-5, 5, 1000)
    radius = -1/0.03
    z = 100 - np.sqrt(radius**2 - y**2) - radius
    ax.plot(z, y, color='r', label="Lens")
    ax.set_xlabel("z (mm)")
    ax.set_ylabel("y (mm)")



    #Ouput plane
    ax.axvline(x=focal_point, color='blue', linestyle='--', label='Output Plane')
    ax.legend()

    #Zoomed in

  
    print(focal_point)

    ax_position = 0.62, 0.17, 0.2, 0.2
    ax_zoomed = fig.add_axes(ax_position)
    ax_zoomed.set_xlim(focal_point-0.0015, focal_point+0.0005)
    ax_zoomed.set_ylim(-0.0000005, 0.0000005)
    formatter_inset = ScalarFormatter(useOffset=False, useMathText=False)
    ax_zoomed.xaxis.set_major_formatter(formatter_inset)
    ax_zoomed.yaxis.set_major_formatter(formatter_inset)
    ax_zoomed.xaxis.set_major_locator(MultipleLocator(0.001))
    ax_zoomed.tick_params(axis='x', which="both", labelsize=10)
    ax_zoomed.tick_params(axis='y', which="both", labelsize=0, length=0)
    ax_zoomed.yaxis.set_ticks([])

    ray_label = True
    for x in ray.values():
        vertices = x.vertices()
        x = [v[0] for v in vertices]
        y = [v[1] for v in vertices]
        z = [v[2] for v in vertices]
        if ray_label is True:
            ax_zoomed.plot(z,y, color="black", label="Rays")
            ray_label = False
        else:
            ax_zoomed.plot(z,y, color="black")
        
    
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
    con3 = ConnectionPatch(xyA=(focal_point+0.0015, 0.0000005), coordsA=ax.transData,
                        xyB=(focal_point+0.0005, 0.0000005), coordsB=ax_zoomed.transData,
                        arrowstyle="-", color="blue")
    con4 = ConnectionPatch(xyA=(focal_point-0.0005, 0.0000005), coordsA=ax.transData,
                        xyB=(focal_point-0.0015, 0.0000005), coordsB=ax_zoomed.transData,
                        arrowstyle="-", color="blue")

    for con in [con3, con4]:
        ax.add_patch(con)
    
    ax_zoomed.axvline(x=focal_point, color='blue', linestyle='--', label='Output Plane')


    return [fig, focal_point]


def task12():
    """
    Task 12.

    In this function you should create a RayBunble and propagate it to the output plane
    before plotting the tracks of the rays.
    This function should return the matplotlib figure of the track plot.

    Returns:
        Figure: the track plot.
    """
    sr = SphericalRefraction(z_0=100, curvature=0.03, n_1=1., n_2=1.5, aperture=5.)
    rb = RayBundle()
    rb.propagate_bundle([sr])

    focal_point = sr.focal_point()
    sr_output_plane = OutputPlane(z_0=focal_point)
    rb.propagate_bundle([sr_output_plane])

    return rb.track_plot()


def task13():
    """
    Task 13.

    In this function you should again create and propagate a RayBundle to the output plane
    before plotting the spot plot.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the spot plot
    2. the simulation RMS

    Returns:
        tuple[Figure, float]: the spot plot and rms
    """
    sr = SphericalRefraction(z_0=100, curvature=0.03, n_1=1., n_2=1.5, aperture=10.)
    rb = RayBundle()
    rb.propagate_bundle([sr])

    focal_point = sr.focal_point()
    sr_output_plane = OutputPlane(z_0=focal_point)
    rb.propagate_bundle([sr_output_plane])

    return [rb.spot_plot(), rb.rms()]


def task14():
    """
    Task 14.

    In this function you will trace a number of RayBundles through the optical system and
    plot the RMS and diffraction scale dependence on input beam radii.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the diffraction scale plot
    2. the simulation RMS for input beam radius 2.5
    3. the diffraction scale for input beam radius 2.5

    Returns:
        tuple[Figure, float, float]: the plot, the simulation RMS value, the diffraction scale.
    """

    sr = SphericalRefraction(z_0=100, curvature=0.03, n_1=1., n_2=1.5, aperture=50.)
    focal_point = sr.focal_point()
    sr_output_plane = OutputPlane(z_0=focal_point)


    radii = np.arange(0.1, 10, 0.1)
    rms = []
    delta_x = []
    for x in radii:
        rb = RayBundle(rmax=x, multi=10)
        rb.propagate_bundle([sr])
        rb.propagate_bundle([sr_output_plane])
        rms.append(rb.rms())
        delta_x.append((588e-06 * (focal_point-100)) / (2*x))

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(radii,rms, label="RMS", color="blue")
    ax2.plot(radii, delta_x, label="Diffraction Scale", color="red")
    ax1.set_xlabel("radius (mm)")
    ax1.set_ylabel("RMS (mm)")
    ax2.set_ylabel("$\Delta x$ (mm)")
    ax1.set_yticks(np.arange(0, 0.31,0.05))
    ax2.set_yticks(np.arange(0, 0.31,0.05))
    ax1.set_ylim(0,0.31)
    ax2.set_ylim(0,0.31)
    ax1.legend(loc=(0.85,0.9))
    ax2.legend(loc=(0.7,0.8))


    rb_25 = RayBundle(rmax=2.5)
    rb_25.propagate_bundle([sr])
    rb_25.propagate_bundle([sr_output_plane])
    delta_x_25 = (588e-06 * (focal_point-100)) / 5
    return [fig, rb_25.rms(), delta_x_25]


def task15():
    """
    Task 15.

    In this function you will create plano-convex lenses in each orientation and propagate a RayBundle
    through each to their respective focal point. You should then plot the spot plot for each orientation.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the spot plot for the plano-convex system
    2. the focal point for the plano-convex lens
    3. the matplotlib figure object for the spot plot for the convex-plano system
    4  the focal point for the convex-plano lens


    Returns:
        tuple[Figure, float, Figure, float]: the spot plots and rms for plano-convex and convex-plano.
    """

    #plano-convex
    ray_pc = RayBundle()
    pc = PlanoConvex(z_0=100, curvature1=0., curvature2=-0.02, \
                     n_inside=1.5168, n_outside=1., thickness=5., aperture=50.)
    ray_pc.propagate_bundle(pc.elements())
    focal_point_pc = pc.focal_point()
    output_plane = OutputPlane(z_0=focal_point_pc)
    ray_pc.propagate_bundle([output_plane])

    #convex-plano
    ray_cp = RayBundle()
    cp = PlanoConvex(z_0=100, curvature1=0.02, curvature2=0., \
                     n_inside=1.5168, n_outside=1., thickness=5., aperture=50.)
    ray_cp.propagate_bundle(cp.elements())
    focal_point_cp = cp.focal_point()
    output_plane = OutputPlane(z_0=focal_point_cp)
    ray_cp.propagate_bundle([output_plane])

    return [ray_pc.spot_plot(), focal_point_pc , ray_cp.spot_plot(), focal_point_cp]

def task16():
    """
    Task 16.

    In this function you will be again plotting the radial dependence of the RMS and diffraction values
    for each orientation of your lens.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the diffraction scale plot
    2. the RMS for input beam radius 3.5 for the plano-convex system
    3. the RMS for input beam radius 3.5 for the convex-plano system
    4  the diffraction scale for input beam radius 3.5

    Returns:
        tuple[Figure, float, float, float]: the plot, RMS for plano-convex, RMS for convex-plano, diffraction scale.
    """

    #Plano-Convex
    pc = PlanoConvex(z_0=100, curvature1=0., curvature2=-0.02, \
                     n_inside=1.5168, n_outside=1., thickness=5., aperture=50.)
    focal_point_pc = pc.focal_point()
    output_plane_pc = OutputPlane(z_0=focal_point_pc)

    radii = np.arange(0.1, 10, 0.1)
    rms_pc = []
    delta_x_pc = []
    for x in radii:
        rb = RayBundle(rmax=x)
        rb.propagate_bundle(pc.elements())
        rb.propagate_bundle([output_plane_pc])
        rms_pc.append(rb.rms())
        delta_x_pc.append((588e-06 * (focal_point_pc-105)) / (2*x))
    
    fig_pc, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(radii,rms_pc, label="RMS", color="blue")
    ax2.plot(radii, delta_x_pc, label="Diffraction Scale", color="red")
    ax1.set_xlabel("radius (mm)")
    ax1.set_ylabel("RMS (mm)")
    ax2.set_ylabel("$\Delta x$ (mm)")
    ax1.set_yticks(np.arange(0, 0.31,0.05))
    ax2.set_yticks(np.arange(0, 0.31,0.05))
    ax1.set_ylim(0,0.31)
    ax2.set_ylim(0,0.31)
    ax1.legend(loc=(0.1,0.9))
    ax2.legend(loc=(0.1,0.8))
    ax1.set_title("Plano-Convex")

    #Convex-Plano
    cp = PlanoConvex(z_0=100, curvature1=0.02, curvature2=0., \
                     n_inside=1.5168, n_outside=1., thickness=5., aperture=50.)
    focal_point_cp = cp.focal_point()
    output_plane_cp = OutputPlane(z_0=focal_point_cp)

    radii = np.arange(0.1, 10, 0.1)
    rms_cp = []
    delta_x_cp = []
    for x in radii:
        rb = RayBundle(rmax=x)
        rb.propagate_bundle(cp.elements())
        rb.propagate_bundle([output_plane_cp])
        rms_cp.append(rb.rms())
        delta_x_cp.append((588e-06 * (focal_point_cp-105)) / (2*x))

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(radii,rms_cp, label="Convex-Plano RMS", color="blue")
    ax1.plot(radii, rms_pc, label="Plano-Convex RMS", color="green")
    ax2.plot(radii, delta_x_cp, label="Diffraction Scale", color="red")
    ax1.set_xlabel("radius (mm)")
    ax1.set_ylabel("RMS (mm)")
    ax2.set_ylabel("$\Delta x$ (mm)")
    ax1.set_yticks(np.arange(0, 0.31,0.05))
    ax2.set_yticks(np.arange(0, 0.31,0.05))
    ax1.set_ylim(0,0.31)
    ax2.set_ylim(0,0.31)
    ax1.legend(loc=(0.1,0.8))
    ax2.legend(loc=(0.1,0.7))

    #3.5
    ray_35_pc = RayBundle(rmax=3.5)
    ray_35_pc.propagate_bundle(pc.elements())
    ray_35_pc.propagate_bundle([output_plane_pc])

    ray_35_cp = RayBundle(rmax=3.5)
    ray_35_cp.propagate_bundle(cp.elements())
    ray_35_cp.propagate_bundle([output_plane_cp])

    delta_x_pc = (588e-06 * (focal_point_pc-105)) / (2*3.5)

    return [fig, ray_35_pc.rms(), ray_35_cp.rms(), delta_x_pc]


if __name__ == "__main__":

    # Run task 8 function
    #task8()

    # Run task 10 function
    FIG10 = task10()

    # Run task 11 function
    FIG11, FOCAL_POINT = task11()

    # Run task 12 function
    FIG12 = task12()

    # Run task 13 function
    FIG13, TASK13_RMS = task13()

    # Run task 14 function
    FIG14, TASK14_RMS, TASK14_DIFF_SCALE = task14()

    # Run task 15 function
    FIG15_PC, FOCAL_POINT_PC, FIG15_CP, FOCAL_POINT_CP= task15()

    # Run task 16 function
    FIG16, PC_RMS, CP_RMS, TASK16_DIFF_SCALE = task16()

    plt.show()

