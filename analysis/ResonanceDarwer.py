import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np

def find_indx_for_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def transmission_func(f, y, k_ex, k_i, h, f_offset, y0):
    z = y0 + np.power(np.abs(1 + 2 * 1j * k_ex * (f - f_offset - 1j * (k_i + k_ex)) /
                             (np.power((f - f_offset - 1j * (k_i + k_ex)), 2) - np.power(h, 2))), 2) - y
    return z

def reflection_func(f, y, k_ex, k_i, h, f_offset, y0):
    z = y0 + np.power(np.abs(2 * k_ex * h /
                             (np.power((1j * (f - f_offset) + (k_i + k_ex)), 2) + np.power(h, 2))), 2) - y
    return z

# Variables and Values to plot:
max_Detuning = int(150e6)  # in Hz
delta_f = int(0.1e6)  # in Hz
k_ex_0 = 7  # in MHz
k_i_0 = 7  # in MHz
h_0 = 11  # in MHz
f_offset = 0  # in Hz
y0 = 0
f = np.arange(-max_Detuning, max_Detuning + delta_f, delta_f)

x = f / 1e6
y_T0 = transmission_func(f, 0, k_ex_0 * 1e6, k_i_0 * 1e6, h_0 * 1e6, f_offset, y0)
y_R0 = reflection_func(f, 0, k_ex_0 * 1e6, k_i_0 * 1e6, h_0 * 1e6, f_offset, y0)
HM = y_T0[max_Detuning//delta_f] + (1 - y_T0[max_Detuning//delta_f])/2
FWHM = abs(2 * x[find_indx_for_nearest(y_T0, HM)])
minimum_T = y_T0[int(len(y_T0)/2)+1]
max_R = y_R0[int(len(y_R0)/2)+1]


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
T, = ax1.plot(x, y_T0, label='Transmission', linewidth=2.5)
ax1.set_ylabel('T', fontsize=28)
ax1.set_xlabel('$\Delta [MHz]$', fontsize=28)
# ax1.minorticks_on()
# ax1.grid(visible=True, which='both', axis='both', linestyle='-.')
ax1.set_title('FWHM = %1.1f[MHz], ' % FWHM + '$T_{min}$ = %2.1f%%, ' % (100 * minimum_T) + '$R_{max}$ = %2.1f%%' % (100 * max_R))
ax1.set_ylim([0, 1])

ax2 = ax1.twinx()
R, = ax2.plot(x, y_R0, label='Reflection', color=color)
ax2.set_ylabel('R', fontsize=28, color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim([0, 1])

# Update values
def update(val):
    Kex = s_Kex.val
    h = s_h.val
    Ki = s_Ki.val
    y_T = transmission_func(f, 0, Kex * 1e6, Ki * 1e6, h * 1e6, f_offset, y0)
    y_R = reflection_func(f, 0, Kex * 1e6, Ki * 1e6, h * 1e6, f_offset, y0)
    HM = y_T[max_Detuning // delta_f] + (1 - y_T[max_Detuning // delta_f]) / 2
    FWHM = abs(2 * x[find_indx_for_nearest(y_T, HM)])
    minimum_T = y_T[int(len(y_T) / 2) + 1]
    max_R = y_R[int(len(y_R) / 2) + 1]
    T.set_data(x, y_T)
    # ax1.set_title('FWHM = %1.1f[MHz]' % FWHM)
    ax1.set_title('FWHM = %1.1f[MHz], ' % FWHM + '$T_{min}$ = %2.1f%%, ' % (100 * minimum_T) + '$R_{max}$ = %2.1f%%' % (
                100 * max_R))
    R.set_data(x, y_R)
    fig.canvas.draw_idle()


s_Kex.on_changed(update)
s_h.on_changed(update)
s_Ki.on_changed(update)

plt.show()

