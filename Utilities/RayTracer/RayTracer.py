
class Ray:

    def __init__(self, z, y, theta, n):
        self.z = z
        self.y = y,
        self.theta = theta
        self.n = n
        pass


class Surface:

    def __init__(self, z, n):
        self.z = z
        self.n = n
        pass


class FlatSurface(Surface):

    def __init__(self, z_s, n_s):
        """
        z_s - the position of the surface along the optical axis (assuming it's orthogonal to the optical axis)
        n_s - the index of refraction after the surface
        """
        super().__init__(z_s, n_s)
        pass

    def propagate(self, ray):
        """
        Get a ray and propagate it through the Flat Surface.
        Return z_out, y_out, theta_out and refracting index after surface
        """
        z_out = None
        y_out = None
        theta_out = None
        n_out = None
        return (z_out, y_out, theta_out, n_out)

class CurvedSurface(Surface):

    def __init__(self, R_c, z_c, n_c):
        super().__init__(z_c, n_c)
        self.R_c = R_c
        pass

class RayTracer:

    def __init__(self):
        pass

    def create_rays(self):

        # TODO: make this NP array

        rays_range = [10, 20, 30]
        rays = []
        for ray in rays_range:
            rays.append(Ray(ray, 0))
            pass

        return rays

    def run(self):
        print('Running Ray Tracer...')
        pass


if __name__ == "__main__":

    ray_tracer = RayTracer()
    rays = ray_tracer.create_rays()
    ray_tracer.run(rays)

    pass
