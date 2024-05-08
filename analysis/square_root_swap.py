import time

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


class SquareRootOfSwap:

    def __init__(self):

        self.r_sigma = 0.18  # %
        self.r_pi = 0.13  # %

        self.gamma_0 = 3*1e6 * 2*np.pi # 3 MHz
        self.kappa_0 = 30*1e6 * 2*np.pi # 30 MHz

        self.sigma_g = 6*1e6

        self.h = 1*1e6 * 2*np.pi  # 1 MHz
        self.kappa_i = 6*1e6 * 2*np.pi  # 6 MHz
        self.kappa_ex = 30*1e6 * 2*np.pi  # 30*10e6
        self.kappa_0 = self.kappa_ex + self.kappa_i

        self.g0 = 16*1e6 * 2*np.pi  # 16 MHz
        self.g1 = 16*1e6 * 2*np.pi  # 16 MHz
        self.g2 = 16*1e6 * 2*np.pi  # 16 MHz

        self.gamma = self.gamma_0
        self.kappa = self.kappa_0

        self.h_squared = self.h**2
        self.h_k_ratio = self.h / self.kappa
        self.h_k_ratio_squared = self.h_k_ratio**2

        self.g0_squared = self.g0**2
        self.g1_squared = self.g1**2
        self.g2_squared = self.g2**2
        self.g0_abs_squared = self.abs_squared(self.g0)
        self.g1_abs_squared = self.abs_squared(self.g1)
        self.g2_abs_squared = self.abs_squared(self.g2)

        self.r_sigma_squared = self.r_sigma**2
        self.r_pi_squared = self.r_pi**2

        self.C_0 = self.g0_squared.real / (2 * self.kappa * self.gamma)
        self.C_tot = (self.g1_squared**2 + self.g2_squared**2).real / (2 * self.kappa * self.gamma)

        self.p = (2 * self.kappa_ex) / (self.kappa + self.h_squared/self.kappa)
        self.D = self.calc_D()

        pass

    def transmission(self):
        alpha0 = self.alpha0()
        alpha1 = self.alpha1()
        alpha2 = self.alpha2()

        alpha0_squared = self.abs_squared(alpha0)
        alpha1_squared = self.abs_squared(alpha1)
        alpha2_squared = self.abs_squared(alpha2)
        t = 1 + 2 * alpha1.real + alpha0_squared + alpha1_squared + alpha2_squared
        return t

    def reflection(self):
        beta0 = self.beta0()
        beta1 = self.beta1()
        beta2 = self.beta2()

        beta0_squared = self.abs_squared(beta0)
        beta1_squared = self.abs_squared(beta1)
        beta2_squared = self.abs_squared(beta2)
        r = beta0_squared + beta1_squared + beta2_squared
        return r

    def calc_D(self):
        term_1 = 1
        term_2 = self.h_k_ratio_squared
        term_3 = 2 * ((1+self.r_sigma_squared) * self.C_tot + 2 * self.r_pi_squared * self.C_0)
        term_4 = -4 * (0+1j) * (self.h/self.kappa) * (self.r_sigma * self.C_tot + self.r_pi_squared * self.C_0)
        term_sum = term_1 + term_2 + term_3 + term_4
        return term_sum

    def u(self, g):
        term_1 = self.r_sigma * (self.h / self.kappa) * g.conjugate()
        term_2 = (0+1j) * g  # TODO: can be: term_2 = complex(0, g)
        term_sum = term_1 + term_2
        return term_sum

    def v(self, g):
        term_1 = -self.r_sigma * g
        term_2 = (0+1j) * (self.h/self.kappa) * g.conjugate()
        term_sum = term_1 + term_2
        return term_sum

    def w(self, g):
        term_1 = (self.h / self.kappa) * g.conjugate()
        term_2 = (0+1j) * g  # TODO: can be: term_2 = complex(0, g)
        product = self.r_pi * (term_1 + term_2)
        return product

    def alpha1(self):
        nom = self.u(self.g1) * self.u(self.g1.conjugate()) * 2 * self.C_tot / self.D
        denom = (self.g1_abs_squared + self.g2_abs_squared)
        product = -self.p * (1 + nom/denom)
        return product

    def alpha2(self):
        nom = self.u(self.g1) * self.v(self.g2.conjugate()) * 2 * self.C_tot / self.D
        denom = (self.g1_abs_squared + self.g2_abs_squared)
        return self.p * (nom/denom)

    def alpha0(self):
        nom = self.u(self.g1) * self.w(self.g0) * self.C_tot / self.D
        denom = (self.g1_abs_squared + self.g2_abs_squared)
        return self.p * (nom/denom)

    def beta1(self):
        nom = self.u(self.g1) * self.v(self.g1) * 2 * self.C_tot / self.D
        denom = (self.g1_abs_squared + self.g2_abs_squared)
        return self.p * (self.h_k_ratio + (nom/denom))

    def beta2(self):
        nom = self.u(self.g1) * self.u(self.g2) * 2 * self.C_tot / self.D
        denom = (self.g1_abs_squared + self.g2_abs_squared)
        return self.p * (nom/denom)

    def beta0(self):
        nom = self.u(self.g1) * self.w(self.g0.conjugate()) * 2 * self.C_tot / self.D
        denom = (self.g1_abs_squared + self.g2_abs_squared)
        return self.p * (nom/denom)

    def abs_squared(self, z):
        #return abs(z) ** 2
        return z.real ** 2 + z.imag ** 2

    def get_t_r(self):

        t = self.transmission()
        r = self.reflection()

        return t, r

    def go(self):

        # self.fig = plt.figure()
        # current_date_time = time.strftime("%H:%M:%S (%d/%m/%Y)")
        # self.fig.canvas.manager.set_window_title(f'{current_date_time} | T-R simulation')

        num_detunings = 100
        start = -50
        end = 50
        atom_detunings = np.linspace(start, end, num=num_detunings) * 1e6 * 2*np.pi
        cavity_detunings = np.linspace(start, end, num=num_detunings) * 1e6 * 2*np.pi

        # Sanity check - no detunings at all
        self.kappa = self.kappa_0 + complex(0, 0)
        self.gamma = self.gamma_0 + complex(0, 0)
        t, r = self.get_t_r()

        n = len(atom_detunings)
        A = np.zeros((n, n))
        for i, d_c in enumerate(cavity_detunings):
            for j, d_a in enumerate(atom_detunings):
                self.kappa = self.kappa_0 + complex(0, d_c)
                self.gamma = self.gamma_0 + complex(0, d_a)

                # Recalc C_0 and C_tot
                self.C_0 = self.g0_squared.real / (2 * self.kappa * self.gamma)
                self.C_tot = (self.g1_squared ** 2 + self.g2_squared ** 2).real / (2 * self.kappa * self.gamma)

                t, r = self.get_t_r()

                z = (t-r)/(t+r)

                A[i][j] = z
            pass

        #data2d = np.sin(t)[:, np.newaxis] * np.cos(t)[np.newaxis, :]

        #data2d = atom_detunings[:, np.newaxis] * cavity_detunings[np.newaxis, :]

        data2d = A
        fig, ax = plt.subplots()
        #im = ax.imshow(data2d)

        extent0 = [-50, 50, -50, 50]
        im = ax.imshow(data2d, extent=extent0)

        ax.set_title('Pan/Zoom colorbar to shift/scale color mapping')

        # Set scales
        SHOW_SCALES = False
        if SHOW_SCALES:
            r = int(n/2)
            plt.xlim(-r, r)
            plt.ylim(-r, r)

        fig.colorbar(im, ax=ax, label='Interactive colorbar')

        # Dashed line
        SHOW_DASHED = False
        if SHOW_DASHED:
            x_points = np.linspace(start, end, num=num_detunings)
            plt.plot(x_points, x_points, linestyle='dashed')

        plt.show(block=True)
        pass

    def test_complex_stuff(self):

        c1 = 3 + 6j
        c2 = complex(2, 5)

        c3 = c1 + c2
        print(c3)

        c_real = c3.real
        c_img = c3.imag
        conj = c3.conjugate()

        abs_squared = self.abs_squared(c1)
        print(abs_squared)

        pass

    def plot(self):
        pass

if __name__ == "__main__":

    sros = SquareRootOfSwap()

    sros.go()

    pass
