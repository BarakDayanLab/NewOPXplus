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
        # Create figure
        self.fig = plt.figure(figsize=(10, 6))
        plt.title("Magnetic Field vs Z")
        plt.xlabel("Z Displacement [mm]")
        plt.ylabel("Magnetic Field [Gauss]")
        plt.grid(True)

        I = 0.200  # 0.2 Ampere = 200 mA
        I = 1.00  # 0.2 Ampere = 200 mA

        n = 30  # Number of loops
        # I *= n

        a = 0.126  # 12.6 cm
        b = 0.114  # 11.4 cm
        z = np.arange(start=0.01, stop=0.7, step=0.01)  # from 0.01 = 1 mm to 30 mm

        z_dist = 0.40 # cm  40 mm
        x_dist = 1.6  # 160 mm

        Bz = self.rectangle_coil_magnetic_field(I, a, b, z, n)

        # Plot. Translate m back to cm. Translate Tesla to Gauss
        plt.plot(z * 1000, Bz, label=f'I={I*1000} [ma]', marker="o")
        # plt.plot(z * 100, Bz_n / 1e5, label='I=-200 [mA]', marker="o")


        # Place vertical line at our atom cloud center-of-mass - 30mm:
        plt.axvline(x=30, color='red', linewidth=1, linestyle='--')
        plt.legend()

        plt.show(block=True)

        pass


if __name__ == '__main__':

    cmf = CalcMagneticField()
    cmf.run()
