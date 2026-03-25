import numpy as np
import matplotlib.pyplot as plt
from qiskit_dynamics import Solver, Signal
from qiskit_ibm_runtime import QiskitRuntimeService

from scipy.integrate import quad

import os, sys, argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from qec import config, utils

logpath = config.LOG_DIR
plotpath = config.PLOT_DIR

logger = utils.setup_logging(log_path=logpath + "/leakage_error.txt")

# Initialize parameters
durations = np.linspace(2.0, 30.0, 15)  # Sweep from 2ns to 30ns
beta_values = np.linspace(-3, 3, 50)  # Testing a range of Beta


# --------------------------------------------------------
# ------- Define Hamiltonian + Parameters ----------------
# --------------------------------------------------------
v01 = 5.0            # Qubit frequency (GHz)
alpha = -0.33       # Anharmonicity (GHz) --> v12 = 5.0 - 0.33 = 4.67 GHz
y0 = np.array([1, 0, 0], dtype=complex) # wave function |psi(t)> initiated at |0>
a = np.array([[0, 1, 0], [0, 0, np.sqrt(2)], [0, 0, 0]]) # lowering operator
adag = a.T.conj() # rising operator
h_static = 2 * np.pi * v01 * adag @ a + np.pi * alpha * (adag @ adag @ a @ a) # Full energy levels in the lab frame
h_rot = 2 * np.pi * v01 * adag @ a # ONLY the 5.0 GHz part for the rotating frame reference

'''
Hamiltonian term for driving adjacent energy levels such as |0><1| & |1><2|; 
Dipole operator of the transmon
Multiplied by the MW field later in the Hamiltonian
h_drive: the operator or field that couples to Re[Ω(t) e^{i ω t}]
Qiskit constructs H_drive(t) = Re[Ω(t) e^{i ω t}] * (a + a†)
'''
h_drive = 2 * np.pi * (a + adag)
solver = Solver(static_hamiltonian=h_static, hamiltonian_operators=[h_drive], rotating_frame=h_rot) # time-independent hamiltonian

# --------------------------------------------------------
# ------- Optimize duration + Plot -----------------------
# --------------------------------------------------------
logger.info('Start optimizing MW-Pulse duration')
leakage_results = []
for d in durations:    
    # Update sigma and amp per duration
    sigma = d / 6
    amp = utils.find_physical_amplitude(np.pi, sigma, d) 
    duration_beta_leakage = utils.find_best_beta(beta_values, amp, d, v01, solver) # sweep for best beta value for d
    duration_bestbeta = beta_values[np.argmin(duration_beta_leakage)]
    duration_envelope = utils.sweep_envelope(amp, d, duration_bestbeta)
    sig_obj = Signal(envelope=duration_envelope, carrier_freq=v01)
    sol = solver.solve(t_span=[0, d], y0=np.array([1,0,0], dtype=complex), 
                       signals=[sig_obj]) # Solve only for the final state to save time
    leakage_results.append(np.max(np.abs(sol.y[:, 2])**2)) # save maximum leakage

# Plotting the "Leakage Cliff"
plt.plot(durations, np.array(leakage_results) * 100)
plt.axhline(0.1, color='red', linestyle='--', label='0.1% Threshold') # not a hard limit
plt.yscale('log')
plt.xlabel('Duration (ns)')
plt.ylabel('Max Leakage (%)')
plt.title(f'Finding the Optimum Pulse Duration for IBM-Fez')
plt.legend()
plt.savefig(plotpath + '/maxLeakage_duration.pdf')


# --------------------------------------------------------
# ------- Optimize Best beta value + Plot ----------------
# --------------------------------------------------------
logger.info('Start optimizing Beta-value in H_drive Hamiltonian')
service = QiskitRuntimeService()
backend = service.backend("ibm_fez")
props = backend.properties()
duration = props.gate_length('sx',12) # dummy qubit_id = 12

duration = 24. # ns; Sx gate opertoring time (t_sx); props.gate_length('sx',72)  
sigma = duration/6 # ns; skinnier gaussian operator
t_eval = np.linspace(0, duration, 200)
rotation = np.pi
amp_pi = utils.find_physical_amplitude(rotation, sigma, duration) # For X-gate (theta = pi)

logger.info(f'Required Amplitude for {rotation} rotation = {amp_pi:.2f}, for duration = {duration} ns')

beta_values = np.linspace(-5, 2, 100)  # Testing a range of Beta
beta_leakages = utils.find_best_beta(beta_values, amp_pi, duration, v01, solver)
best_beta = beta_values[np.argmin(beta_leakages)]
minleakage_betaopt = np.min(beta_leakages)

# Beta Calibration Curve
plt.figure(figsize=(8, 5))
plt.plot(beta_values, np.array(beta_leakages) * 100, color='blue', lw=2)
plt.axvline(best_beta, color='red', linestyle='--', label=f'Best Beta: {best_beta:.3f}')
plt.yscale('log') # Log scale helps see the "null" point
plt.xlabel('Beta Value')
plt.ylabel('Max Leakage (%)')
plt.title('DRAG Calibration Sweep (IBM-Fez Simulation)')
plt.legend()
plt.grid(True, which="both", alpha=0.3)
plt.savefig(plotpath + '/leakage_beta.pdf')

logger.info(f"Optimal Beta found: {best_beta:.4e}")
logger.info(f"Minimum Leakage at this Beta: {minleakage_betaopt*100:.3e}%")


# --------------------------------------------------------
# ------- Final leakage computation + Plot ---------------
# --------------------------------------------------------
logger.info('Make Leakage plot, almost Done!')

drag_envelope = utils.sweep_envelope(amp_pi, duration, best_beta) # gate_signal: time-dependent amplitude Ω(t)
gate_signal = Signal(envelope=drag_envelope, carrier_freq=v01) # time-dependent amplitude of the external field
sol = solver.solve(t_span=[0, duration], y0=y0, signals=[gate_signal], t_eval=t_eval) # solves time-dependent Schrodinger equation

pop0 = np.abs(sol.y[:, 0])**2
pop1 = np.abs(sol.y[:, 1])**2
pop2 = np.abs(sol.y[:, 2])**2

plt.figure(figsize=(8, 5))
#plt.plot(t_eval, pop0, label='State |0> (Qubit)', color='blue')
#plt.plot(t_eval, pop1, label='State |1> (Qubit)', color='green')
plt.plot(t_eval, pop2, label='State |2> (Leakage)', color='red', lw=2)
plt.title(f'DRAG Corrected Pulse (Max Leakage: {np.max(pop2)*100:.2f}%)')
plt.legend()
plt.savefig(plotpath + '/leakage_population.pdf')

logger.info(f"New Max Leakage: {np.max(pop2)*100:.2f}%")