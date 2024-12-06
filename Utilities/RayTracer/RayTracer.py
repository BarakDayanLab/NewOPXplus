import copy
import numpy as np
import matplotlib.pyplot as plt


# ------------------
# TODO: Questions
# ------------------
# 1. Do we need to implement it as ABCD formalism?
# ------------------

class Ray:

    def __init__(self, z, y, theta, n):
        self.z = z
        self.y = y
        self.theta = theta
        self.n = n

        self.num_of_beams = 10

        self.z = np.full(self.num_of_beams, z)
        self.delta_y = 1
        self.y = np.arange(start=0, stop=self.num_of_beams, step=self.delta_y)
        self.theta = np.full(self.num_of_beams, theta)
        self.n = np.full(self.num_of_beams, n)

        return


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
        plt.axvline(x=self.z, color='red', linewidth=2.5, linestyle='-')

        # Add a label near the vertical line
        plt.text(
            self.z+0.5, 0, f'n={self.n}',  # Position slightly to the right of x=25, and centered at y=0
            color='red',
            fontsize=8,
            verticalalignment='center',
            horizontalalignment='left'
        )
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

    def propagate(self, ray):

        # TODO: implement

        resulting_ray = None
        return resulting_ray

class RayTracer:

    def __init__(self):
        self.elements = []

        # Prepare figure for plot
        self.fig = plt.figure(figsize=(10, 6))

        pass

    def clear_all_elements(self):
        self.elements = []

    def add_element(self, element):
        self.elements.append(element)

        # Add the element before the last
        #self.elements.insert(len(self.elments)-1, element)
        pass

    def propagate(self, ray):
        print('Running Ray Tracer...')

        # Iterate over all elements - from last to first
        intermediate_ray = ray
        for element in self.elements:
            # Plot the element
            element.draw()

            # beams_before_refraction = Ray.clone(intermediate_ray)
            beams_before_refraction = copy.deepcopy(intermediate_ray)
            intermediate_ray = element.propagate(intermediate_ray)
            self.plot_ray(beams_before_refraction, intermediate_ray)

        # Plot the last rays to infinity
        self.plot_ray(intermediate_ray, None)

        return intermediate_ray


    def plot_ray(self, source_ray, target_ray):

        if not hasattr(self, 'fig'):
            self.fig = plt.figure(figsize=(10, 6))

        for i in range(0, source_ray.num_of_beams-1):
            # Extract components
            x = source_ray.z[i]
            y = source_ray.y[i]
            theta = np.radians(source_ray.theta[i])  # Convert to radians

            # Calculate end point
            if target_ray is None:
                delta_x = 10
                x_end = x + delta_x * np.cos(theta)
                y_end = y + delta_x * np.sin(theta)
            else:
                x_end = target_ray.z[i]
                y_end = target_ray.y[i]

            # Draw the line
            # plt.plot([x, x_end], [y, y_end], marker="o", label=f"Vector ({x:.1f}, {y:.1f}, {theta:.1f}°)")
            plt.plot([x, x_end], [y, y_end], label=f"Vector ({x:.1f}, {y:.1f}, {theta:.1f}°)")


        # Configure the plot
        plt.title("Vectors Visualization")
        plt.xlabel("X-coordinate")
        plt.ylabel("Y-coordinate")
        plt.axhline(0, color='gray', linewidth=0.5, linestyle="--")
        plt.axvline(0, color='gray', linewidth=0.5, linestyle="--")
        plt.grid(True)
        #plt.legend(loc="upper left", fontsize="small")
        plt.show()

if __name__ == "__main__":

    ray_tracer = RayTracer()

    # Create the system - add all elements
    ray_tracer.add_element(FlatSurface(z_s=5, n_s=2))
    ray_tracer.add_element(FlatSurface(z_s=10, n_s=1.5))

    # Create a ray
    ray = Ray(z=0, y=0, theta=45, n=1)

    # Propagate the ray through the system
    ray = ray_tracer.propagate(ray)

    pass
