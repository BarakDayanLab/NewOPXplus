import time

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


class SquareRootOfSwap:

    def __init__(self):

        self.r_sigma = 0.18  # %
        self.r_pi = 0.13  # %

        self.gamma = 3 * 1e6 * 2*np.pi  # 3 MHz
        self.kappa_total = 27 * 1e6 * 2*np.pi  # 30 MHz

        self.sigma_g = 6*1e6

        self.h = 0.6*1e6 * 2*np.pi  # 1 MHz
        self.kappa_i = 2*1e6 * 2*np.pi  # 6 MHz
        self.kappa_ex = 30*1e6 * 2*np.pi  # 30*1e6
        self.kappa_total = self.kappa_ex + self.kappa_i

        self.avg_g = 10*1e6*2*np.pi  # 16 MHz
        self.g0 = self.avg_g  # 16 MHz
        self.g1 = self.avg_g  # 16 MHz
        self.g2 = self.avg_g  # 16 MHz

        self.g0_squared = self.g0**2
        self.g1_squared = self.g1**2
        self.g2_squared = self.g2**2
        self.g0_abs_squared = self.abs_squared(self.g0)
        self.g1_abs_squared = self.abs_squared(self.g1)
        self.g2_abs_squared = self.abs_squared(self.g2)

        self.gamma_detuned = self.gamma
        self.kappa_total_detuned = self.kappa_total

        self.h_squared = self.h**2

        self.r_sigma_squared = self.r_sigma**2
        self.r_pi_squared = self.r_pi**2

        self.calc_terms()

        pass

    def calc_terms(self):

        self.h_k_ratio = self.h / self.kappa_total_detuned
        self.h_k_ratio_squared = self.h_k_ratio**2

        self.C_0 = self.g0_squared.real / (2 * self.kappa_total_detuned * self.gamma_detuned)

        self.C_tot = (self.g1_squared + self.g2_squared).real / (2 * self.kappa_total_detuned * self.gamma_detuned)

        self.p = (2 * self.kappa_ex) / (self.kappa_total_detuned + (self.h_squared / self.kappa_total_detuned))
        self.D = self.calc_D()

        pass

    def transmission(self):
        alpha0 = self.alpha0()
        alpha1 = self.alpha1()
        alpha2 = self.alpha2()

        alpha0_squared = self.abs_squared(alpha0)
        alpha1_squared = self.abs_squared(alpha1)
        alpha2_squared = self.abs_squared(alpha2)
        t = 1 + (2 * alpha1.real) + alpha1_squared + alpha2_squared + alpha0_squared
        return t

    def reflection(self):
        beta0 = self.beta0()
        beta1 = self.beta1()
        beta2 = self.beta2()

        beta0_squared = self.abs_squared(beta0)
        beta1_squared = self.abs_squared(beta1)
        beta2_squared = self.abs_squared(beta2)
        r = beta1_squared + beta2_squared + beta0_squared
        return r

    def calc_D(self):
        term_1 = 1
        term_2 = self.h_k_ratio_squared
        term_3 = 2 * ((1+self.r_sigma_squared) * self.C_tot + 2 * self.r_pi_squared * self.C_0)
        term_4 = (0-4j) * (self.h / self.kappa_total_detuned) * (self.r_sigma * self.C_tot + self.r_pi_squared * self.C_0)
        term_sum = term_1 + term_2 + term_3 + term_4
        return term_sum

    def u(self, g):
        term_1 = self.r_sigma * (self.h / self.kappa_total_detuned) * g.conjugate()
        term_2 = (0+1j) * g  # TODO: can be: term_2 = complex(0, g)
        term_sum = term_1 + term_2
        return term_sum

    def v(self, g):
        term_1 = -self.r_sigma * g
        term_2 = (0+1j) * (self.h / self.kappa_total_detuned) * g.conjugate()
        term_sum = term_1 + term_2
        return term_sum

    def w(self, g):
        term_1 = (self.h / self.kappa_total_detuned) * g.conjugate()
        term_2 = (0+1j) * g  # TODO: can be: term_2 = complex(0, g)
        product = self.r_pi * (term_1 + term_2)
        return product

    def alpha1(self):
        term1 = self.u(self.g1) * self.u(self.g1.conjugate()) / (self.g1_abs_squared + self.g2_abs_squared)
        term2 = 2 * self.C_tot / self.D
        product = -self.p * (1 + term1 * term2)
        return product

    def alpha2(self):
        term1 = self.u(self.g1) * self.v(self.g2.conjugate()) / (self.g1_abs_squared + self.g2_abs_squared)
        term2 = 2 * self.C_tot / self.D
        product = self.p * term1 * term2
        return product

    def alpha0(self):
        term1 = self.u(self.g1) * self.w(self.g0) / (self.g1_abs_squared + self.g2_abs_squared)
        term2 = 2 * self.C_tot / self.D
        product = self.p * term1 * term2
        return product

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

        n = len(atom_detunings)
        A = np.zeros((n, n))
        for i, d_c in enumerate(cavity_detunings):
            for j, d_a in enumerate(atom_detunings):
                self.kappa_total_detuned = complex(self.kappa_total, d_c)
                self.gamma_detuned = complex(self.gamma, d_a)

                # Recalc C_0 and C_tot
                self.calc_terms()

                # T-R -> ideally we get 50%/50%
                # Get transmissions and reflections
                t, r = self.get_t_r()
                print(f't={t}, r={r}')

                # SQRT-SWAP
                z = (t-r)/(t+r)

                # Sum transaction + reflections
                z = (t+r)

                # Infidelity
                beta1 = self.abs_squared(self.beta1())
                alpha2 = self.abs_squared(self.alpha2())
                z = (beta1 + alpha2) / (t+r)

                A[i][j] = z

            pass

        #data2d = np.sin(t)[:, np.newaxis] * np.cos(t)[np.newaxis, :]

        #data2d = atom_detunings[:, np.newaxis] * cavity_detunings[np.newaxis, :]

        data2d = A
        fig, ax = plt.subplots()
        #im = ax.imshow(data2d)

        extent0 = [-50, 50, -50, 50]
        im = ax.imshow(data2d, origin='lower', extent=extent0)

        ax.set_title('Pan/Zoom colorbar to shift/scale color mapping')

        # Set scales
        SHOW_SCALES = False
        if SHOW_SCALES:
            r = int(n/2)
            plt.xlim(-r, r)
            plt.ylim(-r, r)

        fig.colorbar(im, ax=ax, label='Interactive colorbar')

        # Dashed line
        SHOW_DASHED = True
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
