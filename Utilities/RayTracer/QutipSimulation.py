import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from qutip import *
import os
import json


def my_tests():

    creation_op = create(3)
    destruction_op = destroy(3)
    sigma_z_op = sigmaz()
    identity_op = identity(4)
    qzero_op = qzero(3)
    # cnot_gate = cnot()
    # sqrtswap_gate = sqrtswap()

    print("Creation operator matrix:\n", creation_op)
    print("Destruction operator matrix:\n", destruction_op)
    print("Sigma Z operator matrix:\n", sigma_z_op)
    print("Identity operator matrix:\n", identity_op)
    print("QZero operator matrix:\n", qzero_op)

    pass


"""Simulate with Qutip"""
def simulate_with_qutip():
    # Define parameters
    g = 2 * np.pi * 0.004  # Coupling constant (you can adjust this value)
    nth = 0.016
    G = 1 / 42.7  # Decay rate (you can adjust this value)
    Gup = G * nth  # Decay rate (you can adjust this value)
    Gphi = 1 / 1.68  # Dephasing rate (you can adjust this value)
    kappa = 1 / 6.5e3

    # Define the time interval for the simulation
    t_start = 0.0
    t_end = 5e3
    num_points = 200
    t_eval = np.linspace(t_start, t_end, num_points)  # Time points where the solution is evaluated

    # Define the system of ODEs in terms of real variables
    def odes(t, y):
        # Unpack the current values of the variables
        alpha_r, alpha_i, beta_r, beta_i = y

        # Define the ODEs
        dalpha_r_dt = -g * beta_i - kappa / 2 * alpha_r
        dalpha_i_dt = g * beta_r - kappa / 2 * alpha_i
        dbeta_r_dt = -g * alpha_i - G / 2 * beta_r
        dbeta_i_dt = g * alpha_r - G / 2 * beta_i

        return [dalpha_r_dt, dalpha_i_dt, dbeta_r_dt, dbeta_i_dt]

    # Initial conditions
    alpha_r0 = 1.0  # Real part of alpha at t=0
    alpha_i0 = 0.0  # Imaginary part of alpha at t=0
    beta_r0 = 0.0  # Real part of beta at t=0
    beta_i0 = 0.0  # Imaginary part of beta at t=0
    y0 = [alpha_r0, alpha_i0, beta_r0, beta_i0]

    # Numerically solve the ODEs
    sol = solve_ivp(odes, [t_start, t_end], y0, t_eval=t_eval)

    # Extract solutions
    t = sol.t
    alpha_r_num = sol.y[0]
    alpha_i_num = sol.y[1]
    beta_r_num = sol.y[2]
    beta_i_num = sol.y[3]

    # Combine real and imaginary parts into complex numbers
    alpha_num = alpha_r_num + 1j * alpha_i_num
    beta_num = beta_r_num + 1j * beta_i_num

    # Compute the analytical solutions
    # Calculate omega, handling potential complex values
    omega = np.sqrt(g ** 2 - G ** 2 / 16 + 0j)

    # Analytical expressions for alpha(t) and beta(t)
    alpha_analytical = np.exp(-G * t / 4) * (np.cos(omega * t) + (G / (4 * omega)) * np.sin(omega * t))
    beta_analytical = (1j * g / omega) * np.exp(-G * t / 4) * np.sin(omega * t)

    # =========================
    # Define Operators
    # =========================

    # Qubit A operators (first qubit)
    sm_A = tensor(destroy(2), qeye(2))  # lowering operator for A
    sz_A = tensor(sigmaz(), qeye(2))  # Pauli Z for A

    # Qubit B operators (second qubit)
    sm_B = tensor(qeye(2), destroy(2))  # lowering operator for B
    sp_B = tensor(qeye(2), create(2))  # lowering operator for B
    sz_B = tensor(qeye(2), sigmaz())  # Pauli Z for B

    c_ops = []
    c_ops.append(np.sqrt(G) * sm_B)
    c_ops.append(np.sqrt(Gup) * sp_B)
    c_ops.append(np.sqrt(2 * Gphi) * sz_B)
    c_ops.append(np.sqrt(kappa) * sm_A)

    delta_zero = 0  # Zero detuning

    # Define Hamiltonian for Δ = 0
    H_zero = 0.5 * delta_zero * sz_A + g * (sm_A.dag() * sm_B + sm_A * sm_B.dag())

    # Initial state |e,g>
    psig = tensor(basis(2, 1), basis(2, 0))  # psi_g = |e=1,g=0> = |e> X |g>
    psie = tensor(basis(2, 1), basis(2, 1))  # psi_e = |e=1,g=1> = |e> X |e>
    psi0 = (1 - nth) * ket2dm(psig) + nth * ket2dm(psie)

    # Define solver options with increased nsteps
    solver_opts_zero = Options(nsteps=1e6, store_states=False)

    # Solve the master equation with specified options
    print("Simulating dynamics for Δ = 0...")
    result_zero = mesolve(H_zero, psi0, t_eval, c_ops, [(1 - sz_A) / 2, (1 - sz_B) / 2], options=solver_opts_zero)

    # Extract expectation value of sigma_z for qubit A
    expect_A_zero = result_zero.expect[0]
    expect_B_zero = result_zero.expect[1]

    # # Perform fitting
    # try:
    #     popt_zero, pcov_zero = curve_fit(fit_func, t_list_zero, expect_A_zero, p0=p0, maxfev=100000)
    #     fit_zero = fit_func(t_eval, *popt_zero)
    # except RuntimeError:
    #     fit_zero = None
    #     print("Fitting for Δ = 0 failed.")

    # Initial guess curve for plotting
    # initial_guess = fit_func(t_list_zero, *p0)

    # =========================
    # Plot Dynamics for Δ = 0
    # =========================

    Gp = 2 * g ** 2 / G

    return t_eval, expect_B_zero

def plot_avg_data(folder_path, xticks=None, yticks=None, data_label=None, save_fig=None,
              choose_points=None, add_qutip_simulation=None):

    input_json_path = os.path.join(folder_path, 'averaged_result.json')
    with open(input_json_path, 'r') as infile:
        extracted_data = json.load(infile)

    # Separate the extracted data into lists for plotting
    weighted_avg_photon_number = [data['weighted_average_photon_number'] * 1e2 for data in extracted_data]
    weighted_avg_photon_number_error = [data['weighted_average_photon_number_error'] * 1e2 for data in extracted_data]
    beat_duration = [data['beat_duration'] / 10e2 for data in extracted_data]  # Convert to ms

    # Sort data
    # Get sorted indices
    sorted_indices = sorted(range(len(beat_duration)), key=lambda k: beat_duration[k])
    # Sort all lists using the same indices
    avg_photon_number = [weighted_avg_photon_number[i] for i in sorted_indices]
    avg_photon_number_error = [weighted_avg_photon_number_error[i] for i in sorted_indices]
    beat_duration = [beat_duration[i] for i in sorted_indices]

    if choose_points:
        avg_photon_number = avg_photon_number[:16]
        avg_photon_number_error = avg_photon_number_error[:16]
        beat_duration = beat_duration[:16]

    # Plotting
    plt.errorbar(beat_duration, avg_photon_number, avg_photon_number_error, fmt='-o', capsize=5, label=data_label)

    if add_qutip_simulation:
        t_eval, expect_B_zero = simulate_with_qutip()
        plt.plot(t_eval, expect_B_zero*100, 'black', markersize=2, label='Simulation Data')

    plt.xlabel(f'Beat Duration ($\\mu$s)')
    # plt.ylabel(r'T1 ($\mu$s)')
    plt.ylabel(r'Average Photon Number (%)')
    plt.grid(True)

    plt.legend(loc="upper right", frameon=False, fontsize=10)

    title = r'Average Photon Number ($\overline{n}$) vs. Beat Duration'
    # if measurement_time is not None:
    #     measurement_time = measurement_time / 10e5
    #     title += '\n' + f'Measurement time is {measurement_time} ms'
    plt.title(title)

    # Set custom ticks if provided
    if xticks is not None:
        plt.xticks(xticks)
    if yticks is not None:
        plt.yticks(yticks)

    # Save plot to PNG file
    if save_fig:
        output_plot_path = os.path.join(folder_path, 'Temperature_vs_Beat_Duration_with_Simulation.png')

        plt.savefig(output_plot_path, format='png')
        print('Plot was saved in folder', output_plot_path)
        plt.show()
    else:
        print('Finished process, did not save figure')
        plt.show()
    #     if isinstance(folder_path, list):
    #         output_plot_path = os.path.join(folder_path[0], 'Temperature_vs_Beat_Duration.png')
    #     else:
    #         output_plot_path = os.path.join(folder_path, 'Temperature_vs_Beat_Duration.png')
    #
    #     plt.savefig(output_plot_path, format='png')
    #     plt.show()
    #     print('Plot was saved in folder', output_plot_path)
    # else:
    #     print('Finished process, did not save figure')


# For comparing measurements of the same beat duration
if __name__ == "__main__":

    if simulate := True:
        simulate_with_qutip()
    else:
        parent_folder_path = r""  # Replace with your parent folder path

        plot_avg_data(parent_folder_path, data_label='beat = delta', save_fig=True, choose_points=False,
                      add_qutip_simulation=False)


