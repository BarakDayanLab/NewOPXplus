import numpy as np
import math
import matplotlib.pyplot as plt


class CalcMagneticField:

    def __init__(self):
        pass

    def rectangle_coil_magnetic_field(self, I, a, b, z, n):
        """
        Calculate the magnetic field B_z at a distance z above the center of a rectangular coil

        Inputs:
        - I: Current in the coil (Amperes)
        - a: Half-length of the rectangle in x direction (Meters)
        - b: Half-length of the rectangle in y direction (Meters)
        - z: distance above the center of the coil (Meters)
        - n: number of loops of the coil

        Output:

        - B_z - Magnetic field strength in the z-direction (Tesla)
        """
        # Permeability of free space
        mu_0 = 4 * np.pi * 1e-7  # H/m

        # Calculate the magnetic field. Note the use of numpy's element-wise operations
        B_z = (n * mu_0 * I) / (4 * np.pi) * (
                (a ** 2) / ((a / 2) ** 2 + z ** 2) * (1 / np.sqrt((a ** 2 / 2) + z ** 2)) +
                (b ** 2) / ((b / 2) ** 2 + z ** 2) * (1 / np.sqrt((b ** 2 / 2) + z ** 2)) )

        # Convert from Tesla to Gauss
        B_z *= 1e4

        return B_z

    def run(self):

        I = 0.200  # 0.2 Ampere = 200 mA
        #I = 1.00  # 0.2 Ampere = 200 mA

        # Create figure
        self.fig = plt.figure(figsize=(10, 6))
        plt.title(f'Magnetic Field vs Z (I={I*1000} mA)')
        plt.xlabel("Z Displacement [mm]")
        plt.ylabel("Magnetic Field [Gauss]")
        plt.grid(True)

        n = 30  # Number of loops
        a = 0.126  # 12.6 cm
        b = 0.114  # 11.4 cm
        z = np.arange(start=0.01, stop=0.3, step=0.01)  # from 0.01 = 1 mm to 30 mm

        z_dist = 0.04  # 4 cm
        x_dist = 0.16  # 16 cm

        Bz_1 = self.rectangle_coil_magnetic_field(I, a, b, z, n)
        Bz_2 = self.rectangle_coil_magnetic_field(I, a, b, z + z_dist, n)
        Bz_2 = Bz_2[::-1]
        Bz_sum = Bz_1 + Bz_2

        # Plot. Translate m back to cm. Translate Tesla to Gauss
        plt.plot(z * 1000, Bz_1, label=f'Top Coil', marker="o")
        plt.plot(z * 1000, Bz_2, label=f'Bottom Coil', marker="o")
        plt.plot(z * 1000, Bz_sum, label=f'Sum', marker="o")

        # Place vertical line at our atom cloud center-of-mass - 30mm:
        plt.axvline(x=30, color='red', linewidth=1, linestyle='--')
        plt.legend()

        plt.show(block=True)

        pass


if __name__ == '__main__':

    cmf = CalcMagneticField()
    cmf.run()
