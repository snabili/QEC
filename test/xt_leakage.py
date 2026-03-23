
from qiskit_ibm_runtime import QiskitRuntimeService
from qiskit_experiments.library import StandardRB
#from qiskit_experiments.framework import BatchExperiment, ParallelExperiment
from qiskit_experiments.library.randomized_benchmarking.rb_analysis import RBAnalysis

from qiskit_aer.noise import NoiseModel, thermal_relaxation_error, depolarizing_error
from qiskit.quantum_info import average_gate_fidelity, Operator
from qiskit_aer import AerSimulator

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import networkx as nx
import numpy as np
from collections import defaultdict
from scipy.linalg import expm
from scipy.optimize import curve_fit

import os, sys, argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from qec import config, utils

logpath = config.LOG_DIR
plotpath = config.PLOT_DIR

logger = utils.setup_logging(log_path=logpath + "/xtalk_gate-errors.txt")


service = QiskitRuntimeService()
backend = service.backend("ibm_fez")
props = backend.properties()

# ------- Preparing/extracting info from backend -------------------
nqubits = len(props.qubits)
t1_time = np.array([props.t1(i) for i in range(nqubits)])
t2_time = np.array([props.t2(i) for i in range(nqubits)])

cz_errors = {tuple(sorted(g.qubits)): g.parameters[0].value for g in props.gates if g.gate == "cz"}
cz_time   = {tuple(sorted(g.qubits)): g.parameters[1].value for g in props.gates if g.gate == "cz"}

target_pair = min(cz_errors, key=cz_errors.get) # find a quite pair --> xtalk becomes sensible
u, v = target_pair
spectator = utils.find_spectator(target_pair, backend)
logger.info(f'Target: {target_pair}, Spectator: {spectator}')

t1 = props.t1(spectator)
t2 = props.t2(spectator)
tg = props.gate_length("sx", spectator)
t_cz = props.gate_length("cz", target_pair)

# Sanity check to find the validity of thermal_relaxation_error formula
spec_err_manual = (1/6)*(1-np.exp(-tg/t1)) + (1/3)*(1-np.exp(-tg/t2)) # from np.float(N) --> N: spec_err.item()

spec_noise = thermal_relaxation_error(t1, t2, tg)
avg_fid = average_gate_fidelity(spec_noise)
avg_err_calculated = 1 - avg_fid

logger.info(f"Manual Formula result: {spec_err_manual:.8f}")
logger.info(f"Calculated from Object: {avg_err_calculated:.8f}")
logger.info(f"Difference: {abs(spec_err_manual - avg_err_calculated):.8e}")

# --------- NoiseModel: iso & xt modes ---------------------
# --------- N.B.: EPC values here are not ------------------
# ---------     realistic as sx is included redundently ----
bgates = ["id", "rz", "sx", "x", "cz"] # defined instructions in qiskit

'''
IBM tunable coupler Heron machines have extreme
    low values frequency values, in the order of 1 kHz 
    ==> crosstalk effect is negligible
    To be able to see some effect, xi is increased to 100 kHz
'''
xi = 100e03  # Hz
theta = 2 * np.pi * xi * t_cz  # radians
Z = np.array([[1, 0],[0, -1]])
U_z = expm(1j * theta * Z)
z_unitary = Operator(U_z)

logger.info(f'  Tunable coupler frequency = {xi:.2e} Hz')

noise_iso = NoiseModel()
noise_iso.add_quantum_error(spec_noise, "sx", [spectator])
sim_iso = AerSimulator(noise_model=noise_iso, basis_gates=bgates)

# Simultaneous noise model
noise_xt = NoiseModel()
noise_xt.add_quantum_error(spec_noise, "sx", [spectator])
noise_xt.add_quantum_error(z_unitary, "sx", [spectator])
sim_xt = AerSimulator(noise_model=noise_xt, basis_gates=bgates)

# -------- Define RB experiments WITH EXPLICIT ANALYSIS ATTACHED
lengths = [1, 10, 50, 100, 150, 200, 400, 800]
nsamples = 50
nshots = 10000
nseeds = 123

rb_iso = StandardRB((spectator,), lengths=lengths, num_samples=nsamples, seed=nseeds)
rb_iso.analysis = RBAnalysis()   # <-- CRITICAL

rb_xt = StandardRB((spectator,), lengths=lengths, num_samples=nsamples, seed=nseeds)
rb_xt.analysis = RBAnalysis() # <-- CRITICAL

# -------- Run experiments and analysis
data_iso = rb_iso.run(sim_iso, shots=nshots).block_for_results()
data_sim = rb_xt.run(sim_xt, shots=nshots).block_for_results()

# -------- 
'''
Extract decay curves and EPC: Computing using two methods:
    1 - Data extracted methonds: Useful when the error values are not astronomically small :)
    2 - Write a fitting function with scipy: Useful to have control over fitting parameters
    '''
def get_decay_points(expdata):    
    records = expdata.data()       
    groups = defaultdict(list) # group survival values by sequence length    
    for r in records:
        L = r["metadata"]["xval"]
        counts = r["counts"]
        shots = sum(counts.values())
        survival = counts.get("0", 0) / shots
        groups[L].append(survival)    
    # average survival for each length
    x, y = [], []
    for L in sorted(groups.keys()):
        x.append(L)
        y.append(np.mean(groups[L]))    
    return np.array(x), np.array(y)

x_iso, y_iso = get_decay_points(data_iso)
x_sim, y_sim = get_decay_points(data_sim)


def get_epc(expdata):
    res = expdata.analysis_results("EPC")
    return res.value.n, res.quality

epc_iso_data, q_iso_data = get_epc(data_iso)
epc_sim_data, q_sim_data = get_epc(data_sim)
Delta_Data = epc_sim_data - epc_iso_data

logger.info('EPC and Crosstalk values from Data')
logger.info(f'   Isolated EPC = {epc_iso_data:.2e}, Quality: {q_iso_data:s}')
logger.info(f'   Simultaneous EPC = {epc_sim_data:.2e}, Quality: {q_sim_data:s}')
logger.info(f'   ΔEPC Crosstalk = {Delta_Data:.2e}')

plt.figure(figsize=(7,5))
plt.plot(x_iso, y_iso, "o", label="Isolated RB")
plt.plot(x_sim, y_sim, "s", label="Simultaneous RB")

plt.xlabel("Sequence length")
plt.ylabel("Survival probability")
plt.title(f"RB decay on {target_pair} with spectator {spectator}")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.yscale('log')
plt.savefig(plotpath + 'xt_effect.pdf')


# ---- Method2: Manually fitting with scipy ------------
def rb_model(m, A, alpha, B):
    return A * (alpha ** m) + B
def fit_rb_decay(x, y):
    # initial guesses
    A0 = y[0] - y[-1]
    B0 = y[-1]
    alpha0 = 0.999  # shallow decay
    popt, pcov = curve_fit(rb_model, x, y, p0=[A0, alpha0, B0], bounds=([0, 0, 0], [1, 1, 1]))
    A, alpha, B = popt
    return A, alpha, B
def epc_from_alpha(alpha): # for single-qubit: (d-1)/d * (1-alpha)
    return 0.5 * (1 - alpha)

A_iso, alpha_iso, B_iso = fit_rb_decay(x_iso, y_iso)
A_sim, alpha_sim, B_sim = fit_rb_decay(x_sim, y_sim)

epc_iso = epc_from_alpha(alpha_iso)
epc_sim = epc_from_alpha(alpha_sim)
Delta_Manual = epc_sim - epc_iso

logger.info('EPC and Crosstalk values from Manual Computation:')
logger.info(f'  alpha_iso = {alpha_iso:.2e}, alpha_sim = {alpha_sim:.2e}')
logger.info(f'  EPC_iso = {epc_iso:.2e}, EPC_sim = {epc_sim:.2e}')
logger.info(f'  ΔEPC = {Delta_Manual:.2e}')