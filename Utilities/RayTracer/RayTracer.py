import copy
import numpy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc


class Ray:

    def __init__(self, z, y, theta, n, num_of_beams=None, beams_delta_y=None):
        """
        Ray has few properties: (z, y, theta, n)
        Use num_of_beams and beams_delta_y - to define few beams in the ray (use odd number)
        """
        self.num_of_beams = num_of_beams
        self.beams_delta_y = beams_delta_y

        if num_of_beams is not None:
            self.set_beams(z, y, theta, n, num_of_beams, beams_delta_y)
        else:
            self.z = np.array(z)
            self.y = np.array(y)
            self.theta = np.array(theta)
            self.n = np.array(n)

    def set_beams(self, z, y, theta, n, num_of_beams, beams_delta_y):
        self.z = np.full(self.num_of_beams, z)
        if num_of_beams == 1:
            self.y = np.array([y])
        else:
            half = np.floor(num_of_beams/2) * beams_delta_y
            self.y = np.arange(start=y-half, stop=y+half+self.beams_delta_y, step=self.beams_delta_y)
        self.theta = np.full(self.num_of_beams, theta)
        self.n = np.full(self.num_of_beams, n)

    def __str__(self):
        return f'z={self.z[0]}, y={self.y}, theta={self.theta[0]} n={self.n[0]}'


class Surface:
    """
    This is mainly a superclass. It can be inherited by other classes to implement 'propagate' and 'draw' functions.
    """
    def __init__(self, z, n):
        self.z = z
        self.n = n
        pass

    def propagate(self, ray):
        raise Exception('Surface did not implement propagate function.')

    def draw(self):
        raise Exception('Surface class did not implement draw function.')


class FlatSurface(Surface):

    def __init__(self, z_s, n_s):
        super().__init__(z_s, n_s)
        pass

    def draw(self):
        # Draw a vertical line representing flat surface, and add a label for the index of refraction
        plt.axvline(x=self.z, color='red', linewidth=2.5, linestyle='-')
        plt.text(self.z+0.5, 0, f'n={self.n}', color='red', fontsize=12, verticalalignment='center', horizontalalignment='left')
        pass

    def propagate(self, ray):
        """
        Get a ray and propagate it through the Flat Surface.
        Return the resulting ray - z_out, y_out, theta_out and refracting index after surface
        """

        ray.y = ray.y + (self.z - ray.z) * np.tan(np.radians(ray.theta))

        # Z-out position is the same as the Surface position
        ray.z = np.full(len(ray.y), self.z)

        # Exit angle: Theta_2 = arcsin( n1/n2 * sin(Theta_1)
        ray.theta = np.degrees(np.arcsin(ray.n/self.n * np.sin(np.radians(ray.theta))))

        # Index of refraction at exit is the Surface index of refraction
        ray.n = np.full(len(ray.y), self.n)

        return ray


class CurvedSurface(Surface):

    def __init__(self, R_c, z_c, n_c):
        super().__init__(z_c, n_c)
        self.R_c = R_c

    def draw(self):
        r_x = self.z
        radius = self.R_c
        start_angle = 0
        end_angle = 180

        # Draw the arc representing the curved surface
        arc = Arc(
            (r_x, 0),  # Center of the arc (x=r_x, y=0)
            width=2 * radius,  # Diameter of the arc in x (2 * radius)
            height=2 * radius,  # Diameter of the arc in y (same as x for a circle)
            angle=-90,  # Rotation angle of the arc (0 degrees for no tilt)
            theta1=start_angle,  # Starting angle of the arc
            theta2=end_angle,  # Ending angle of the arc
            color='blue',  # Arc color
            linewidth=2  # Line thickness
        )

        # Add the arc to the plot
        plt.gca().add_patch(arc)

        # Draw text specifying the curved surface properties
        if self.R_c < 0:
            offset = self.z + self.R_c + 0.2
        else:
            offset = self.z + self.R_c - 0.8
        plt.text(offset, 0, f'R={self.R_c}\nn={self.n}', color='blue', fontsize=12, verticalalignment='center', horizontalalignment='left')
        pass

    def propagate(self, ray):

        # Find the z_2 coordinate of the exiting ray
        z1 = ray.z
        theta = np.radians(ray.theta)
        tan_theta = np.tan(theta)

        A = 1 + tan_theta ** 2
        B = 2 * tan_theta * ray.y - 2 * (tan_theta ** 2) * ray.z - 2 * self.z
        C = self.z ** 2 + tan_theta ** 2 * (ray.z ** 2) - 2 * tan_theta * ray.z * ray.y + ray.y ** 2 - self.R_c ** 2

        # Calculate the quadratic equation solution. If value inside sqrt is negative, some values may become nan
        # (this is ok with us as the ray in the display will not appear to propagate)
        sol1 = (-B + np.sqrt(B**2 - 4*A*C)) / (2 * A)
        sol2 = (-B - np.sqrt(B**2 - 4*A*C)) / (2 * A)

        if self.R_c < 0:
            ray.z = sol2
        else:
            ray.z = sol1

        # Find the y_2 coordinate of the exiting ray
        ray.y = ray.y + tan_theta * (ray.z - z1)

        # Find the theta_t angle of the exiting ray
        theta_r = np.arcsin(ray.y / self.R_c)
        n_ratio = ray.n / self.n
        theta_t = np.arcsin(n_ratio * np.sin(theta_r - theta)) - theta_r
        ray.theta = np.degrees(theta_t)
        ray.theta = -ray.theta

        # Ray's refraction index at exit is the Surface's refraction index
        ray.n = np.full(len(ray.y), self.n)

        return ray


class RayTracer:
    """
    This class is like a canvas where elements and rays are placed, and then we can simulate the rays propagation.
    """
    def __init__(self):
        self.elements = []
        pass

    def add_element(self, element):
        self.elements.append(element)
        pass

    def propagate(self, ray):
        """
        This function propagates a ray through all the elements, and saves its state at each segment - this is a journey.
        """
        journey = [ray]
        interim_ray = copy.deepcopy(ray)
        for element in self.elements:
            interim_ray = element.propagate(interim_ray)
            journey.append(copy.deepcopy(interim_ray))
        return journey


class RayStudio:
    """
    This class uses the RayTracer object and allows it to be interactive, by registering to key events and changing
    properties of the ray or elements and then re-plotting them.
    """
    def __init__(self):
        self.lim = 20
        self.ray_tracer = RayTracer()
        self.prepare_display()

    def prepare_display(self):

        # Prepare figure for plot
        self.fig = plt.figure(figsize=(10, 6))

        # Connect the function to the key press event for interactive control
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)

        # Set title, labels, axis, grid
        plt.title("Rays Studio")
        plt.xlabel("Optical Axis (Z)")
        plt.ylabel("Y-coordinate")
        plt.axhline(0, color='gray', linewidth=0.5, linestyle="--")
        plt.axvline(0, color='gray', linewidth=0.5, linestyle="--")
        plt.grid(True)

    def plot_journey(self, journey):
        """
        1) Plot all the elements
        2) Plot the journey segments of the propagating ray
        """
        # Plot elements
        for element in self.ray_tracer.elements:
            element.draw()

        # Plot the rays in journey
        for i in range(0, len(journey)):
            source_ray = journey[i]
            if i == len(journey)-1:
                target_ray = None
            else:
                target_ray = journey[i+1]

            # Extract components
            x = source_ray.z
            y = source_ray.y

            # Calculate end point
            EXTEND = False  # Allows to extend line backward - to find paraxial focal planes
            if target_ray is None:
                inf_ray_len = 80
                x_end = x + np.full(len(y), inf_ray_len)
                y_end = y + inf_ray_len * np.sin(np.radians(source_ray.theta))
                if EXTEND:
                    y = y - x*np.sin(np.radians(source_ray.theta))
                    x = np.zeros(len(y))
            else:
                x_end = target_ray.z
                y_end = target_ray.y
                if EXTEND:
                    y = y - x*np.sin(np.radians(source_ray.theta))
                    x = np.zeros(len(y))

            # Draw the line
            plt.plot([x, x_end], [y, y_end], marker="o")

        plt.xlim(-self.lim, self.lim)
        plt.ylim(-6, 6)

        # If paraxial focal plane was calculated, plot it as a vertical line
        if hasattr(self, 'paraxial_focal_plane') and self.paraxial_focal_plane is not None:
            plt.axvline(x=self.paraxial_focal_plane, color='red', linewidth=1, linestyle='--')

        # Draw text that specifies ray properties
        plt.text(0.01, 0.99, f'{journey[0]}', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top', horizontalalignment='left')

        plt.show(block=True)
        pass

    def on_key_press(self, event):
        """
        Handle various key events that drive the studio's interactivity
        """
        if not hasattr(self, 'original_ray'):
            self.original_ray = copy.deepcopy(self.ray)

        if event.key == 'y':
            self.ray.y -= 1
        elif event.key == 'Y':
            self.ray.y += 1
        elif event.key == 'z':
            self.ray.z -= 1
        elif event.key == 'Z':
            self.ray.z += 1
        elif event.key == 'n':
            self.ray.n = self.ray.n - 0.1
        elif event.key == 'N':
            self.ray.n = self.ray.n + 0.1
        elif event.key == 'm':
            self.ray_tracer.elements[0].n -= 0.1
        elif event.key == 'M':
            self.ray_tracer.elements[0].n += 0.1
        elif event.key == 't':
            self.ray.theta -= 5
        elif event.key == 'T':
            self.ray.theta += 5
        elif event.key == 'R':
            self.ray_tracer.elements[0].R_c += 1
        elif event.key == 'r':
            self.ray_tracer.elements[0].R_c -= 1
        elif event.key == 'D':
            self.ray_tracer.elements[1].z += 1
        elif event.key == 'd':
            self.ray_tracer.elements[1].z -= 1
        elif event.key == 'x':
            self.lim -= 10
        elif event.key == 'X':
            self.lim += 10
        else:
            return

        # Clear display and re-plot the elements and rays
        plt.clf()
        plt.title("Rays Studio")
        plt.xlabel("Z-coordinate")
        plt.ylabel("Y-coordinate")
        plt.axhline(0, color='gray', linewidth=0.5, linestyle="--")
        plt.axvline(0, color='gray', linewidth=0.5, linestyle="--")
        plt.grid(True)

        # Propagate the modified ray through the system and redraw
        new_journey = self.ray_tracer.propagate(self.ray)

        self.plot_journey(new_journey)
        self.fig.canvas.draw()  # Redraw the figure to display the changes
        pass

    def run(self):

        # Mention what is the test you want to run
        test = 'curved surface CONVEX'

        # Test 1 - Flat Surface
        if test == 'flat surface':
            self.ray_tracer.add_element(FlatSurface(z_s=5, n_s=3))
            self.ray = Ray(z=2, y=0, theta=0, n=1, num_of_beams=5, beams_delta_y=0.5)
            journey = self.ray_tracer.propagate(self.ray)

        # Test 2 - Curved Surface Concave
        if test == 'curved surface CONCAVE':
            self.ray_tracer.add_element(CurvedSurface(R_c=5, z_c=10, n_c=1))
            self.ray = Ray(z=2, y=0, theta=0, n=2, num_of_beams=7, beams_delta_y=0.5)
            # self.ray = Ray(z=[2], y=[1], theta=[10], n=[2])
            journey = self.ray_tracer.propagate(self.ray)

        # Test 3 - Curved Surface Convex
        if test == 'curved surface CONVEX':
            self.ray_tracer.add_element(CurvedSurface(R_c=-5, z_c=10, n_c=2))
            self.ray = Ray(z=2, y=0, theta=0, n=1, num_of_beams=7, beams_delta_y=0.5)
            # self.ray = Ray(z=[2], y=[1], theta=[10], n=[1])
            journey = self.ray_tracer.propagate(self.ray)

        # Test 4 - Thin Lens
        if test == 'thin lens':
            self.ray_tracer.add_element(CurvedSurface(R_c=-6, z_c=4, n_c=1.3))
            self.ray_tracer.add_element(CurvedSurface(R_c=6, z_c=-4, n_c=1))
            # self.ray = Ray(z=[2], y=[1], theta=[0], n=[1])
            self.ray = Ray(z=-6, y=0, theta=0, n=1, num_of_beams=7, beams_delta_y=1)
            journey = self.ray_tracer.propagate(self.ray)

        # Test 4 - Plano-convex
        if test == 'plano-convex':
            self.ray_tracer.add_element(CurvedSurface(R_c=-4, z_c=10, n_c=1.3))
            self.ray_tracer.add_element(FlatSurface(z_s=8, n_s=1))
            self.ray = Ray(z=2, y=0, theta=0, n=1, num_of_beams=7, beams_delta_y=0.5)
            journey = self.ray_tracer.propagate(self.ray)

        # Test 5 - 3 rays propagate into CONVEX-PLAN Thorlabs LA4306 Fused Silica Lens
        if test == 'thorlabs-la3406-convex-plano':
            n = 1.457
            _Z_c = 20
            _R_c = -18.4
            _t_c = 7.1
            _Z_s = _Z_c + _R_c + _t_c
            thorlabs_la4306_focal_length = 40.1
            self.ray_tracer.add_element(CurvedSurface(R_c=_R_c, z_c=_Z_c, n_c=1.457))
            self.ray_tracer.add_element(FlatSurface(z_s=_Z_s, n_s=1))
            self.ray = Ray(z=[2, 2, 2], y=[0.1, 1, 10], theta=[0, 0, 0], n=[1, 1, 1])  # 3-beamed Ray, for measuring abberations
            # self.ray = Ray(z=-12, y=0, theta=0, n=1, num_of_beams=9, beams_delta_y=1)  # Multi-Beams Ray, equal spacing!
            # self.ray = Ray(z=[2, 2, 2, 2, 2], y=[0.1, 1, 3, 5, 10], theta=[0, 0, 0, 0, 0], n=[1, 1, 1, 1, 1])  # 3-beamed Ray, for measuring abberations
            # self.ray = Ray(z=[2], y=[5], theta=[0], n=[1])
            journey = self.ray_tracer.propagate(self.ray)
            H_2 = -thorlabs_la4306_focal_length*(n-1)*_t_c/(n*_R_c)
            paraxial_focal_plane = _Z_s - H_2
            self.paraxial_focal_plane = paraxial_focal_plane + thorlabs_la4306_focal_length


        self.plot_journey(journey)

        pass

    @staticmethod
    def run_studio():
        ray_studio = RayStudio()
        ray_studio.run()
        pass


if __name__ == "__main__":

    RayStudio.run_studio()
    pass
