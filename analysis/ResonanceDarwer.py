import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np
import scipy.stats as stat


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


def transmission_SPRINT(f, y, g1, g2, gamma, k_ex, k_i, h, f_offset, y0):
    da = f - f_offset
    dc = f - f_offset
    k = k_i + k_ex
    if g1 == 0 and g2 == 0:
        C = 0
        g1g1 = 0
        g1g2 = 0
    else:
        C = (g1 ** 2 + g2 ** 2) / (2 * (k + 1j * dc) * (gamma + 1j * da))
        g1g1 = np.power(np.abs(g1), 2) / (np.power(np.abs(g1), 2) + np.power(np.abs(g2), 2))
        g1g2 = g1 * g2 / (np.power(np.abs(g1), 2) + np.power(np.abs(g2), 2))
        
    t_0 = 1 + 2 * 1j * k_ex * (f - f_offset - 1j * (k_i + k_ex)) / \
          (np.power((f - f_offset - 1j * (k_i + k_ex)), 2) - np.power(h, 2))
    r_0 = 2 * k_ex * h / (np.power((1j * (f - f_offset) + (k_i + k_ex)), 2) + np.power(h, 2))
    
    z = y0 + np.power(np.abs((2 * (k + 1j * dc) * k_ex / (np.power((k + 1j * dc), 2) + np.power(h, 2))) * g1g1 *
                             (2 * C / (1 + np.power(h, 2) / np.power((k + 1j * dc), 2) + 2 * C)) + t_0), 2) + \
        np.power(np.abs(r_0 * (2 * k_ex / (k + 1j * dc)) * g1g2 *
                        (2 * C / (1 + np.power(h, 2) / np.power((k + 1j * dc), 2) + 2 * C))), 2) - y
    return z

def transmission_SPRINT_spread(f, y, g1, g2, gamma, k_ex, k_i, h, f_offset, y0,sig_g,min_detected_g):
    
    gg = np.linspace(max(min_detected_g,g1-4*sig_g), g1+4*sig_g, 100)
    dist_g = stat.norm(g1,sig_g)
    weights = dist_g.pdf(gg)
    weights = weights/np.sum(weights)
    z = np.zeros(np.shape(f))
    for i in range(len(gg)):
        z = z + weights[i] * transmission_SPRINT(f, y, gg[i], g2, gamma, k_ex, k_i, h, f_offset, y0) 
    return z
    
def reflection_SPRINT(f, y, g1, g2, gamma, k_ex, k_i, h, f_offset, y0):
    da = f - f_offset
    dc = f - f_offset
    k = k_i + k_ex
    if g1 == 0 and g2 == 0:
        C = 0
        g1g1 = 0
        g1g2 = 0
    else:
        C = (g1 ** 2 + g2 ** 2) / (2 * (k + 1j * dc) * (gamma + 1j * da))
        g1g1 = np.power(np.abs(g1), 2) / (np.power(np.abs(g1), 2) + np.power(np.abs(g2), 2))
        g1g2 = g1 * g2 / (np.power(np.abs(g1), 2) + np.power(np.abs(g2), 2))
    r_0 = 2 * k_ex * h / (np.power((1j * (f - f_offset) + (k_i + k_ex)), 2) + np.power(h, 2))
    z = y0 + \
        np.power(np.abs(r_0 * (1 - g1g1 * (2 * C / (1 + np.power(h, 2) / np.power((k + 1j * dc), 2) + 2 * C)))), 2) + \
        np.power(np.abs((2 * (k + 1j * dc) * k_ex / (np.power((k + 1j * dc), 2) + np.power(h, 2))) * g1g2 *
                        (2 * C / (1 + np.power(h, 2) / np.power((k + 1j * dc), 2) + 2 * C))), 2) - \
        y
    return z

def reflection_SPRINT_spread(f, y, g1, g2, gamma, k_ex, k_i, h, f_offset, y0,sig_g,min_detected_g):
    
    gg = np.linspace(max(0,g1-4*sig_g), g1+4*sig_g, 100)
    dist_g = stat.norm(g1,sig_g)
    weights = dist_g.pdf(gg)
    weights = weights/np.sum(weights)
    z = np.zeros(np.shape(f))
    for i in range(len(gg)):
        z = z + weights[i] * reflection_SPRINT(f, y, gg[i], g2, gamma, k_ex, k_i, h, f_offset, y0)  
    return z


# Variables and Values to plot:
max_Detuning = int(150e6)  # in Hz
delta_f = int(0.1e6)  # in Hz
k_ex_0 = 7  # in MHz
k_i_0 = 7  # in MHz
h_0 = 11  # in MHz
f_offset = 0  # in Hz
y0 = 0
g1_0 = 25  # in MHz
sig_g_0 = 0.1  # in MHz
min_detected_g = 7 # in MHz
g2_0 = 0  # in MHz
gamma = 6.065 / 2
f = np.arange(-max_Detuning, max_Detuning + delta_f, delta_f)

x = f / 1e6
y_T0 = transmission_func(f, 0, k_ex_0 * 1e6, k_i_0 * 1e6, h_0 * 1e6, f_offset, y0)
y_R0 = reflection_func(f, 0, k_ex_0 * 1e6, k_i_0 * 1e6, h_0 * 1e6, f_offset, y0)
y_T0_SPRINT = transmission_SPRINT_spread(f, 0, g1_0 * 1e6, g2_0 * 1e6, gamma * 1e6, k_ex_0 * 1e6, k_i_0 * 1e6, h_0 * 1e6,
                                  f_offset, y0,sig_g_0*1e6,min_detected_g*1e6)
y_R0_SPRINT = reflection_SPRINT_spread(f, 0, g1_0 * 1e6, g2_0 * 1e6, gamma * 1e6, k_ex_0 * 1e6, k_i_0 * 1e6, h_0 * 1e6,
                                  f_offset, y0,sig_g_0*1e6,min_detected_g*1e6)
HM = y_T0[max_Detuning // delta_f] + (1 - y_T0[max_Detuning // delta_f]) / 2
FWHM = abs(2 * x[find_indx_for_nearest(y_T0, HM)])
f0_indx = find_indx_for_nearest(f, f_offset)
T_f0 = y_T0[f0_indx]
T_SPRINT_f0 = y_T0_SPRINT[f0_indx]
R_f0 = y_R0[f0_indx]
R_SPRINT_f0 = y_R0_SPRINT[f0_indx]

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
fig.subplots_adjust(bottom=0.1, top=0.70)

# Create axes for sliders
ax_sig_g = fig.add_axes([0.3, 0.75, 0.4, 0.05])
ax_sig_g.spines['top'].set_visible(True)
ax_sig_g.spines['right'].set_visible(True)

ax_g1 = fig.add_axes([0.3, 0.80, 0.4, 0.05])
ax_g1.spines['top'].set_visible(True)
ax_g1.spines['right'].set_visible(True)

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
s_sig_g = Slider(ax=ax_sig_g, label='$\sigma_g$', valmin=0.1, valmax=40.0, valinit=sig_g_0, valfmt=' %1.1f [MHz]',
              facecolor='#cc7000')
s_g1 = Slider(ax=ax_g1, label='$g$', valmin=0, valmax=100.0, valinit=g1_0, valfmt=' %1.1f [MHz]',
              facecolor='#cc7000')
s_Kex = Slider(ax=ax_Kex, label='$\kappa_{ex} $', valmin=0, valmax=100.0, valinit=k_ex_0, valfmt=' %1.1f [MHz]',
               facecolor='#cc7000')
s_h = Slider(ax=ax_h, label='h ', valmin=0, valmax=20.0, valinit=h_0, valfmt=' %1.1f [MHz]', facecolor='#cc7000')
s_Ki = Slider(ax=ax_Ki, label='$\kappa_{i} $', valmin=0, valmax=20.0, valinit=k_i_0, valfmt=' %1.1f [MHz]',
              facecolor='#cc7000')

T, = ax1.plot(x, y_T0, label='Transmission of bare cavity', linewidth=2.5)
T_sprint, = ax1.plot(x, y_T0_SPRINT, label='Transmission with atoms', linewidth=2.5, color='tab:green')
ax1.set_ylabel('T \ R', fontsize=28)
ax1.set_xlabel('$\Delta [MHz]$', fontsize=28)
ax1.set_title('FWHM bare cavity = %1.1f[MHz]' % FWHM)
ann_T = ax1.annotate('%2.1f%%' % (100 * T_f0), (x[f0_indx], y_T0[f0_indx]),
                     xytext=(0, 10),  # vertical offset.
                     textcoords='offset points', ha="center", size=15)
T_f0_dot, = ax1.plot(x[f0_indx], y_T0[f0_indx], 'o', color='tab:blue')
ann_T_SPRINT = ax1.annotate('%2.1f%%' % (100 * T_SPRINT_f0), (x[f0_indx], y_T0_SPRINT[f0_indx]),
                            xytext=(0, 10),  # vertical offset.
                            textcoords='offset points', ha="center", size=15)
T_SPRINT_f0_dot, = ax1.plot(x[f0_indx], y_T0_SPRINT[f0_indx], 'o', color='tab:green')

R, = ax1.plot(x, y_R0, label='Reflection of bare cavity', linewidth=2.5, color='tab:red')
R_sprint, = ax1.plot(x, y_R0_SPRINT, label='Reflection with atoms', linewidth=2.5, color='tab:orange')
# ax1.set_ylabel('R', fontsize=28, color='tab:red')
# ax1.tick_params(axis='y', labelcolor='tab:red')
ann_R = ax1.annotate('%2.1f%%' % (100 * R_f0), (x[f0_indx], y_R0[f0_indx]),
                     xytext=(0, 10),  # vertical offset.
                     textcoords='offset points', ha="center", size=15)
R_f0_dot, = ax1.plot(x[f0_indx], y_R0[f0_indx], 'o', color='tab:red')
ann_R_SPRINT = ax1.annotate('%2.1f%%' % (100 * R_SPRINT_f0), (x[f0_indx], y_R0_SPRINT[f0_indx]),
                            xytext=(0, 10),  # vertical offset.
                            textcoords='offset points', ha="center", size=15)
R_SPRINT_f0_dot, = ax1.plot(x[f0_indx], y_R0_SPRINT[f0_indx], 'o', color='tab:orange')
legend1 = ax1.legend(loc='center right')
ax1.set_ylim([0, 1])

list_of_plots = [[T, T_f0_dot, ann_T], [T_sprint, T_SPRINT_f0_dot, ann_T_SPRINT],
                 [R, R_f0_dot, ann_R], [R_sprint, R_SPRINT_f0_dot, ann_R_SPRINT]]
lined = dict()
lined_T = dict()
lined_R = dict()

for legend_line, plot_line in zip(legend1.get_lines(), list_of_plots):
    legend_line.set_picker(5)
    lined[legend_line] = plot_line

# Update values
def update(val):
    sig_g = s_sig_g.val
    g1 = s_g1.val
    Kex = s_Kex.val
    h = s_h.val
    Ki = s_Ki.val
    y_T = transmission_func(f, 0, Kex * 1e6, Ki * 1e6, h * 1e6, f_offset, y0)
    y_R = reflection_func(f, 0, Kex * 1e6, Ki * 1e6, h * 1e6, f_offset, y0)
    y_T_SPRINT = transmission_SPRINT_spread(f, 0, g1 * 1e6, g2_0 * 1e6, gamma * 1e6, Kex * 1e6, Ki * 1e6, h * 1e6,
                                     f_offset, y0,sig_g*1e6,min_detected_g*1e6)
    y_R_SPRINT = reflection_SPRINT_spread(f, 0, g1 * 1e6, g2_0 * 1e6, gamma * 1e6, Kex * 1e6, Ki * 1e6, h * 1e6,
                                   f_offset, y0,sig_g*1e6,min_detected_g*1e6)
    HM = y_T[max_Detuning // delta_f] + (1 - y_T[max_Detuning // delta_f]) / 2
    FWHM = abs(2 * x[find_indx_for_nearest(y_T, HM)])
    T_f0 = y_T[f0_indx]
    T_SPRINT_f0 = y_T_SPRINT[f0_indx]
    R_f0 = y_R[f0_indx]
    R_SPRINT_f0 = y_R_SPRINT[f0_indx]

    T.set_data(x, y_T)
    T_f0_dot.set_data(x[f0_indx], T_f0)
    ann_T.set_text('%2.1f%%' % (100 * T_f0))
    ann_T.xy = (x[f0_indx], T_f0)
    T_sprint.set_data(x, y_T_SPRINT)
    T_SPRINT_f0_dot.set_data(x[f0_indx], T_SPRINT_f0)
    ann_T_SPRINT.set_text('%2.1f%%' % (100 * T_SPRINT_f0))
    ann_T_SPRINT.xy = (x[f0_indx], T_SPRINT_f0)

    ax1.set_title('FWHM = %1.1f[MHz]' % FWHM)

    R.set_data(x, y_R)
    R_f0_dot.set_data(x[f0_indx], R_f0)
    ann_R.set_text('%2.1f%%' % (100 * R_f0))
    ann_R.xy = (x[f0_indx], R_f0)
    R_sprint.set_data(x, y_R_SPRINT)
    R_SPRINT_f0_dot.set_data(x[f0_indx], R_SPRINT_f0)
    ann_R_SPRINT.set_text('%2.1f%%' % (100 * R_SPRINT_f0))
    ann_R_SPRINT.xy = (x[f0_indx], R_SPRINT_f0)

    fig.canvas.draw_idle()

def onpick(event):
    legend_line = event.artist
    plot, dot, ann = lined[legend_line]
    vis = not plot.get_visible()
    plot.set_visible(vis)
    dot.set_visible(vis)
    ann.set(visible=vis)
    legend_line.set_alpha(1.0 if vis else 0.2)
    fig.canvas.draw()

s_sig_g.on_changed(update)
s_g1.on_changed(update)
s_Kex.on_changed(update)
s_h.on_changed(update)
s_Ki.on_changed(update)

fig.canvas.mpl_connect('pick_event', onpick)
plt.show()
