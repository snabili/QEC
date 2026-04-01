import stim
import pymatching
from qiskit_ibm_runtime import QiskitRuntimeService

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import brentq

import os, sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from qec import config

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
# importing custom codes
from qec import utils, config

logpath = config.LOG_DIR
datapath = config.DATA_DIR
plotpath = config.PLOT_DIR

logger = utils.setup_logging(log_path=logpath + "/logical_error.txt")


# -------------------------------------------------------------------
# -------- Patch Selection and Role Assignment ----------------------
# -------------------------------------------------------------------
service = QiskitRuntimeService()
backend = service.backend("ibm_fez")
dist = 3
error_threshold = 0.045

clean_G = utils.get_clean_subgraph(backend, error_threshold) # to ignore extremely high-error; error_threshold = 0.15 in utils
centers = utils.get_valid_centers(clean_G, dist)
best_c, stats = utils.get_best_patch_center(backend, clean_G, centers, dist) # Use dist to search for a valid "patch" of qubits
found = utils.find_d3_patches(clean_G, dist)
patch_nodes = [i for i in found if i[0]==best_c][0]
# Qubit roles of: data + x-syndrome & z-syndrome
roles = utils.get_true_roles(clean_G, patch_nodes, dist)
result = utils.find_logical_qubit_full(backend, clean_G, dist)

'''
1st step to get a sensible p_phys error with the machine --> ibm_fez
'''

noise_map, gate_map = utils.get_actual_noise_model(result, roles, clean_G, backend)

def find_worst_case_errors(noise_map, gate_map):
    '''
    Finding maximum gate-errors and readout-errors:
        - Find max(readout-error) and the node: max_read_q
        - Find max(gate_error) and the edge (u, v)
        - Strategic Logic:
            - If readout is the bottleneck:
                - find worst gate attached to max_read_q
                - look all edges in gate_errors that touch max_read_q
                - if max_read_q qubit has no gates, fallback to max(gate_error)
            - If gate is the bottleneck:
                - find the worst readout between the two qubits in that gate
    '''
    max_read_q = max(noise_map, key=lambda q: noise_map[q]['readout'])
    p_read_max = noise_map[max_read_q]['readout']
    max_edge = max(gate_map, key=lambda e: gate_map[e])
    p_gate_max = gate_map[max_edge]
    
    if p_read_max >= p_gate_max:
        relevant_gates = [err for edge, err in gate_map.items() if max_read_q in edge]    
        p_gate_final = max(relevant_gates) if relevant_gates else p_gate_max
        p_read_final = p_read_max        
    else:
        u, v = max_edge
        p_read_final = max(noise_map[u]['readout'], noise_map[v]['readout'])
        p_gate_final = p_gate_max        
    return p_gate_final, p_read_final

worst_pg, worst_pr = find_worst_case_errors(noise_map, gate_map)
logger.info(f'gate-error = {worst_pg:.2e}, readout-error = {worst_pr:.2e}')


# code to simulate/compute threshold for logical error
def get_logical_error_custom(d, p_gate, p_read, shots=100000):
    '''
    Alternative way to compute pL (instead of pL = C x (p_phys/p_thresh)^(d+1)/2)
    Count how many success over total number of shots
    '''
    circuit = stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        distance=d,
        rounds=d,
        after_clifford_depolarization=p_gate, # Gates use the fixed high-performance value
        after_reset_flip_probability=p_gate,
        before_round_data_depolarization=p_gate,               
        before_measure_flip_probability=p_read # Measurement uses the variable 'sweep' value
    )
    
    # Standard decoding logic
    model = circuit.detector_error_model(decompose_errors=True)
    matching = pymatching.Matching.from_detector_error_model(model)
    sampler = circuit.compile_detector_sampler()
    
    defects, actual_obs = sampler.sample(shots=shots, separate_observables=True)
    predicted_obs = matching.decode_batch(defects)
    
    # Calculate Logical Error Rate (PL)
    return np.sum(actual_obs.flatten() != predicted_obs.flatten()) / shots

distances = [3, 5, 7]
readout_errors_range = np.linspace(0.001, worst_pr*5, 15)
logical_errors = {}
plt.figure(figsize=(7, 4))
for d in distances:
    logical_errors[d] = [get_logical_error_custom(d, worst_pg, p) for p in readout_errors_range]
    plt.plot(readout_errors_range, logical_errors[d], label=f'd={d}')



# --- Numerical Intersection Finding ---
try:
    log_x = np.log10(readout_errors_range)
    log_y3 = np.log10(logical_errors[3])
    log_y7 = np.log10(logical_errors[7])

    # Create continuous functions for d=3 and d=7
    f3 = interp1d(log_x, log_y3, kind='linear', fill_value="extrapolate")
    f7 = interp1d(log_x, log_y7, kind='linear', fill_value="extrapolate")

    # The intersection is where f3(x) - f7(x) = 0
    def find_crossing(x):
        return f3(x) - f7(x)

    # Solve for the x-axis intersection (Physical Error)
    # Search between the bounds of our sweep
    log_x_intersect = brentq(find_crossing, log_x[0], log_x[-1])
    x_intersect = 10**log_x_intersect

    # Use the solved x to find the y-axis intersection (Logical Error)
    log_y_intersect = f3(log_x_intersect)
    y_intersect = 10**log_y_intersect

    # --- 3. Visualization and Annotation ---
    # Plot the specific point
    plt.plot(x_intersect, y_intersect, 'ro', markersize=12, label='Threshold Point')
    
    # Vertical line to X-axis
    plt.axvline(x=x_intersect, color='red', linestyle='--', alpha=0.5)
    # Horizontal line to Y-axis
    plt.axhline(y=y_intersect, color='red', linestyle='--', alpha=0.5)

    # Text Labels
    plt.text(x_intersect * 0.85, y_intersect * 0.75, 
             f"X (Phys Thresh): {x_intersect:.2%}\nY (Log Thresh): {y_intersect:.4f}",
             bbox=dict(facecolor='white', alpha=0.8), color='red')

    logger.info(f"Intersection Coordinates:")
    logger.info(f"p_threshold (X Intersection): {x_intersect:.5f} ({x_intersect:.2%})")
    logger.info(f"C (Y Intersection):  {y_intersect:.5f}")

except ValueError:
    logger.info("\n[!] No intersection found.")
    logger.info(f"Current Gate Error (worst_pg = {worst_pg:.4f}) may be above the surface code threshold.")
    logger.info("Try lowering worst_pg to 0.001 to see the crossing point.")

# --- 4. Plot Styling ---
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Physical Readout Error ($P_{readout}$)')
plt.ylabel('Logical Error Rate ($P_L$)')
plt.title(f"Distance Sweep: Finding the Threshold Intersection\n(Fixed Gate Noise: {worst_pg:.2e})")
plt.legend()
plt.grid(True, which="both", ls="-", alpha=0.2)
plt.savefig(plotpath + '/p_threshold.pdf')

def logical_error_formula(pThresh, pPhys, C, pL):
    d = 2*(np.log10(pL/C)/np.log10(pPhys/pThresh))-1
    return d
pL = np.linspace(1e-9, 1e-5, 10)
pPhys = worst_pr
pThresh = x_intersect
C = y_intersect

for p in pL:
    lgerror = logical_error_formula(pThresh, pPhys, C, p)
    logger.info(f' logical error: {p:.2e},  distance value = {lgerror:.2f}')