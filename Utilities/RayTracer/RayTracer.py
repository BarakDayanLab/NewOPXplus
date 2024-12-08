import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc

# ------------------
# TODO: Questions
# ------------------
# 1. Do we need to implement it as ABCD formalism?
# 3. Why doesn't my thick len focus to the exact point?
# 4. Why when using R_c=1 and R_c=0.5, the y,z does not hit the curved surface line?
# ------------------

class Ray:

    def __init__(self, z, y, theta, n, num_of_beams=10, beams_delta_y=1):

        self.num_of_beams = num_of_beams
        self.beams_delta_y = beams_delta_y

        self.set_beams(z, y, theta, n, num_of_beams, beams_delta_y)

        return

    def set_beams(self, z, y, theta, n, num_of_beams, beams_delta_y):
        self.z = np.full(self.num_of_beams, z)
        half = np.floor(num_of_beams/2) * self.beams_delta_y
        self.y = np.arange(start=y-half, stop=y+half+self.beams_delta_y, step=self.beams_delta_y)
        # self.y = np.arange(start=y, stop=y+self.num_of_beams*self.beams_delta_y, step=self.beams_delta_y)
        self.theta = np.full(self.num_of_beams, theta)
        self.n = np.full(self.num_of_beams, n)

    def __str__(self):
        return f'z={self.z[0]}, y={self.y[0]}, theta={self.theta[0]} n={self.n[0]} (#beams = {self.num_of_beams})'


class Surface:

    def __init__(self, z, n):
        self.z = z
        self.n = n
        pass

    def propagate(self, ray):
        # TODO: implement
        pass

class FlatSurface(Surface):

    def __init__(self, z_s, n_s):
        """
        z_s - the position of the surface along the optical axis (assuming it's orthogonal to the optical axis)
        n_s - the index of refraction after the surface
        """
        super().__init__(z_s, n_s)
        pass

    def draw(self):
        # Draw a vertical line and add a label near the vertical line with index of refraction
        plt.axvline(x=self.z, color='red', linewidth=2.5, linestyle='-')
        plt.text(self.z+0.5, 0, f'n={self.n}', color='red', fontsize=12, verticalalignment='center', horizontalalignment='left')
        pass

    def propagate(self, ray):
        """
        Get a ray and propagate it through the Flat Surface.
        Return z_out, y_out, theta_out and refracting index after surface
        """

        ray.y = ray.y + (self.z - ray.z) * np.tan(np.radians(ray.theta))

        # Z-out position is the same as the Surface position
        ray.z = np.full(ray.num_of_beams, self.z)

        # Exit angle: Theta_2 = arcsin( n1/n2 * sin(Theta_1)
        ray.theta = np.degrees(np.arcsin(ray.n/self.n * np.sin(np.radians(ray.theta))))

        # Index of refraction at exit is the Surface index of refraction
        ray.n = np.full(ray.num_of_beams, self.n)

        return ray

class CurvedSurface(Surface):

    def __init__(self, R_c, z_c, n_c):
        super().__init__(z_c, n_c)
        self.R_c = R_c
        pass

    def draw(self):
        r_x = self.z
        radius = self.R_c
        start_angle = 0
        end_angle = 180

        # Draw the arc
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

        # Draw a vertical line and add a label near the vertical line with index of refraction
        #plt.axvline(x=self.z, color='blue', linewidth=2.5, linestyle='--')

        if self.R_c < 0:
            offset = self.z + self.R_c + 0.2
        else:
            offset = self.z + self.R_c - 0.8
        plt.text(offset, 0, f'R={self.R_c}\nn={self.n}', color='blue', fontsize=12, verticalalignment='center', horizontalalignment='left')
        pass

    def propagate(self, ray):

        y1 = ray.y
        z1 = ray.z

        theta = np.radians(ray.theta)
        tan_theta = np.tan(theta)

        A = 1 + tan_theta ** 2
        B = 2 * tan_theta * ray.y - 2 * (tan_theta ** 2) * ray.z - 2 * self.z
        C = self.z ** 2 + tan_theta ** 2 * (ray.z ** 2) - 2 * tan_theta * ray.z * ray.y + ray.y ** 2 - self.R_c ** 2

        sol1 = (-B + np.sqrt(B**2 - 4*A*C)) / (2 * A)
        sol2 = (-B - np.sqrt(B**2 - 4*A*C)) / (2 * A)

        if self.R_c < 1:
            ray.z = sol2
        else:
            ray.z = sol1

        ray.y = ray.y + tan_theta * (ray.z - z1)

        theta_r = np.arcsin(ray.y / self.R_c)


        n_ratio = ray.n / self.n
        theta_delta = theta_r - theta

        theta_t = np.arcsin(n_ratio * np.sin(theta_delta)) - theta_r

        ray.theta = np.degrees(theta_t)

        if self.R_c < 1:
            ray.theta = -ray.theta

        # Index of refraction at exit is the Surface index of refraction
        ray.n = np.full(ray.num_of_beams, self.n)

        return ray

class RayTracer:

    def __init__(self):
        self.elements = []
        pass

    def clear_all_elements(self):
        self.elements = []

    def add_element(self, element):
        self.elements.append(element)

        # Add the element before the last
        #self.elements.insert(len(self.elments)-1, element)
        pass

    def propagate(self, ray):

        journey = [ray]
        interim_ray = copy.deepcopy(ray)
        for element in self.elements:
            interim_ray = element.propagate(interim_ray)
            journey.append(copy.deepcopy(interim_ray))
        return journey


class RayStudio:

    def __init__(self):

        self.lim = 20

        self.ray_tracer = RayTracer()

        self.prepare_display()
        pass

    def prepare_display(self):

        # Prepare figure for plot
        self.fig = plt.figure(figsize=(10, 6))

        # Connect the function to the key press event
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)

        # Set title, labels, axis, grid
        plt.title("Rays Studio")
        plt.xlabel("Optical Axis (Z)")
        plt.ylabel("Y-coordinate")
        plt.axhline(0, color='gray', linewidth=0.5, linestyle="--")
        plt.axvline(0, color='gray', linewidth=0.5, linestyle="--")
        plt.grid(True)

        pass

    def plot_journey(self, journey):

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
            if target_ray is None:
                inf_ray_len = 60
                x_end = x + np.full(len(y), inf_ray_len)
                y_end = y + inf_ray_len * np.sin(np.radians(source_ray.theta))
            else:
                x_end = target_ray.z
                y_end = target_ray.y

            # Draw the line
            # plt.plot([x, x_end], [y, y_end], marker="o", label=f"Vector ({x:.1f}, {y:.1f}, {theta:.1f}°)")
            plt.plot([x, x_end], [y, y_end], marker="o")

        plt.xlim(0, self.lim)
        plt.ylim(-4, 4)


        # Show the ray properties
        plt.text(0.01, 0.99, f'{journey[0]}', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top', horizontalalignment='left')

        plt.show(block=True)
        pass

    def on_key_press(self, event):

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
            self.ray.n = self.ray.n - 0.5
        elif event.key == 'N':
            self.ray.n = self.ray.n + 0.5
        elif event.key == 'm':
            self.ray_tracer.elements[0].n -= 0.2
        elif event.key == 'M':
            self.ray_tracer.elements[0].n += 0.2
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
        elif event.key == '+':
            self.ray.num_of_beams += 2
        elif event.key == '-':
            self.ray.num_of_beams -= 2

        elif event.key == 'x':
            self.lim -= 10
        elif event.key == 'X':
            self.lim += 10
        elif event.key == '0':
            self.ray = copy.deepcopy(self.original_ray)
        else:
            return

        # For debug purposes
        print(f'Changing ray to {self.ray}')

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

        # Test 1
        # self.ray_tracer.add_element(FlatSurface(z_s=10, n_s=3))
        # self.ray = Ray(z=2, y=5, theta=0, n=1, num_of_beams=4, beams_delta_y=0.5)

        # Test 2
        # self.ray_tracer.add_element(CurvedSurface(R_c=-5, z_c=10, n_c=3))
        # self.ray = Ray(z=2, y=-1.5, theta=0, n=1, num_of_beams=7, beams_delta_y=0.5)

        # Test 3 - Thick Lens
        # self.ray_tracer.add_element(CurvedSurface(R_c=-7, z_c=10, n_c=3))
        # self.ray_tracer.add_element(CurvedSurface(R_c=7, z_c=14, n_c=1))
        # self.ray = Ray(z=2, y=1, theta=0, n=1, num_of_beams=4, beams_delta_y=0.5)

        # Test 4 - Plano-convex
        # self.ray_tracer.add_element(CurvedSurface(R_c=-4, z_c=10, n_c=3))
        # self.ray_tracer.add_element(FlatSurface(z_s=8, n_s=1))
        # self.ray = Ray(z=2, y=-1, theta=0, n=1, num_of_beams=7, beams_delta_y=0.5)

        # Test 5 - Thorlabs LA4306 Fused Silica Lens
        _Z_c = 30
        _R_c = 18.4
        _t_c = 7.1
        _Z_s = _Z_c - _R_c + _t_c
        self.ray_tracer.add_element(CurvedSurface(R_c=-_R_c, z_c=_Z_c, n_c=1.46))
        self.ray_tracer.add_element(FlatSurface(z_s=_Z_s, n_s=1))
        self.ray = Ray(z=2, y=0, theta=0, n=1, num_of_beams=5, beams_delta_y=0.2)


        # Propagate the ray through the system
        journey = self.ray_tracer.propagate(self.ray)

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
