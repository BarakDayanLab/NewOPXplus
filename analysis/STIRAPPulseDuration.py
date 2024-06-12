import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider


class STIRAPPulseDuration:

    def __init__(self):
        pass

    def calc_photon_output_efficiency(self, k_tot, k_ex, C):
        C2 = 2 * C
        P_s = (k_ex / k_tot) * (C2 / (C2+1))
        return P_s

    def calc_cooperativity(self, g, k_tot, gamma):
        C = (g*g) / (2*k_tot*gamma)
        return C

    def calc_pulse_duration(self, g, k_ex, k_tot, gamma):
        c = self.calc_cooperativity(g=g, k_tot=k_tot, gamma=gamma)
        p_s = self.calc_photon_output_efficiency(k_tot=k_tot, k_ex=k_ex, C=c)
        tau = k_tot * p_s / (g*g)

        # Move to Nano-seconds
        tau = tau*1000
        return tau

    def get_title(self, k_ex, k_i, g):
        return f'Pulse duration minimum (k_i = {int(k_i)}, k_ex = {int(k_ex)}, g={int(g)})'

    def plot_tau(self):

        gamma = 3  # MHz
        k_i = 12  # MHz
        k_ex = 24  # MHz
        k_tot = k_ex + k_i
        g = 20  # MHz

        # Initial value of n
        initial_k_i = 12
        initial_k_ex = 20

        # Generate k_tot values
        k_tot_values = np.linspace(10, 50, 41)

        y_values = self.calc_pulse_duration(g=g, k_ex=initial_k_ex, k_tot=k_tot_values, gamma=gamma)

        # Create the figure and axes
        fig, ax = plt.subplots()
        plt.subplots_adjust(bottom=0.30)  # Adjust the bottom margin for the slider

        # Plot the dots
        dots, = plt.plot(k_tot_values, y_values, 'bo', label=f"Pulse duration [ns]")


        ax.set_xlabel("Kappa (total)")
        ax.set_ylabel("Duration (Tau) [ns]")

        ax.set_title(self.get_title(initial_k_ex, initial_k_i, g))
        ax.set_ylim(0, 700)  # Set y-axis limits dynamically
        ax.grid(True)
        ax.legend()

        # Add a slider for n
        g_ax_slider = plt.axes([0.15, 0.15, 0.65, 0.03], facecolor='lightgoldenrodyellow')  # rec = [left, bottom, width, height]
        g_slider = Slider(g_ax_slider, 'g', 0, 30, valinit=initial_k_i, valstep=1)

        k_ex_ax_slider = plt.axes([0.15, 0.1, 0.65, 0.03], facecolor='lightgoldenrodyellow')
        k_ex_slider = Slider(k_ex_ax_slider, 'K_ex', 0, 30, valinit=initial_k_ex, valstep=1)

        # Update the plot when the slider value changes
        def update(val):
            new_g = g_slider.val
            new_k_ex = k_ex_slider.val

            updated_y_values = self.calc_pulse_duration(g=new_g, k_ex=new_k_ex, k_tot=k_tot_values, gamma=gamma)

            dots.set_ydata(updated_y_values)

            #ax.set_ylim(min(updated_y_values), max(updated_y_values))  # Set y-axis limits dynamically
            ax.set_title(self.get_title(new_k_ex, initial_k_i, new_g))
            fig.canvas.draw_idle()

        g_slider.on_changed(update)
        k_ex_slider.on_changed(update)

        # Show the plot
        plt.show()


if __name__ == "__main__":

    pd = STIRAPPulseDuration()
    pd.plot_tau()

    pass