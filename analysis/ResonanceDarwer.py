import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np

def transmission_func(f, y, k_ex, k_i, h, f_offset, y0):
    z = y0 + np.power(np.abs(1 + 2 * 1j * k_ex * (f - f_offset - 1j * (k_i + k_ex)) /
                             (np.power((f - f_offset - 1j * (k_i + k_ex)), 2) - np.power(h, 2))), 2) - y
    return z

def reflection_func(f, y, k_ex, k_i, h, f_offset, y0):
    z = y0 + np.power(np.abs(2 * k_ex * h /
                             (np.power((1j * (f - f_offset) + (k_i + k_ex)), 2) + np.power(h, 2))), 2) - y
    return z

# # Variables:
max_Detuning = int(100e6)  # in Hz
delta_f = int(0.1e6)  # in Hz
k_ex_0 = 7  # in MHz
k_i_0 = 7  # in MHz
h_0 = 11  # in MHz
f_offset = 0  # in Hz
y0 = 0
f = np.linspace(-max_Detuning, max_Detuning, delta_f)

# General plot parameters
mpl.rcParams['font.family'] = 'Cambria'
mpl.rcParams['font.size'] = 18
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 2

fig = plt.figure(figsize=(10, 8.5))

# Create main axis
ax1 = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.2, top=0.75)

# Create axes for sliders
ax_Kex = fig.add_axes([0.3, 0.85, 0.4, 0.05])
ax_Kex.spines['top'].set_visible(True)
ax_Kex.spines['right'].set_visible(True)

ax_h = fig.add_axes([0.3, 0.90, 0.4, 0.05])
ax_h.spines['top'].set_visible(True)
ax_h.spines['right'].set_visible(True)

ax_Ki = fig.add_axes([0.3, 0.95, 0.4, 0.05])
ax_Ki.spines['top'].set_visible(True)
ax_Ki.spines['right'].set_visible(True)

# Create sliders
s_Kex = Slider(ax=ax_Kex, label='$K_{ex} $', valmin=0, valmax=100.0, valinit=k_ex_0, valfmt=' %1.1f [MHz]',
               facecolor='#cc7000')
s_h = Slider(ax=ax_h, label='h ', valmin=0, valmax=20.0, valinit=h_0, valfmt=' %1.1f [MHz]', facecolor='#cc7000')
s_Ki = Slider(ax=ax_Ki, label='$K_{i} $', valmin=0, valmax=20.0, valinit=k_i_0, valfmt=' %1.1f [MHz]',
              facecolor='#cc7000')

color = 'tab:red'
T, = ax1.plot(np.array(f) / 1e6, transmission_func(f, 0, k_ex_0 * 1e6, k_i_0 * 1e6, h_0 * 1e6, f_offset, y0), label='Transmission', linewidth=2.5)
ax1.set_ylabel('T', fontsize=28)
ax1.set_xlabel('$\Delta [MHz]$', fontsize=28)
ax1.set_ylim([0, 1])

ax2 = ax1.twinx()
R, = ax2.plot(np.array(f) / 1e6, reflection_func(f, 0, k_ex_0 * 1e6, k_i_0 * 1e6, h_0 * 1e6, f_offset, y0), label='Reflection', color=color)
ax2.set_ylabel('R', fontsize=28, color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim([0, 1])

# Update values
def update(val):
    Kex = s_Kex.val
    h = s_h.val
    Ki = s_Ki.val
    T.set_data(np.array(f) / 1e6, transmission_func(f, 0, Kex * 1e6, Ki * 1e6, h * 1e6, f_offset, y0))
    R.set_data(np.array(f) / 1e6, reflection_func(f, 0, Kex * 1e6, Ki * 1e6, h * 1e6, f_offset, y0))
    fig.canvas.draw_idle()


s_Kex.on_changed(update)
s_h.on_changed(update)
s_Ki.on_changed(update)

plt.show()

