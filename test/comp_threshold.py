import stim
import pymatching
import matplotlib.pyplot as plt
import numpy as np

import os, sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from qec import config

plotpath = config.PLOT_DIR

# code to simulate/compute threshold for logical error

def get_logical_error(d, p_phys, shots=100000):
    '''
    Compute logical error:
        - Simulate 1e+05 stim circuits shots
        - Decode with pymatching
    '''
    p_float = float(p_phys) # p_phys --> std py float; match Stim's C++
    circuit = stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        distance=d,
        rounds=d, 
        after_clifford_depolarization=p_float,
        after_reset_flip_probability=p_float,
        before_measure_flip_probability=p_float,
        before_round_data_depolarization=p_float  # Standard argument name
    )
    model = circuit.detector_error_model(decompose_errors=True)
    matching = pymatching.Matching.from_detector_error_model(model)
    
    sampler = circuit.compile_detector_sampler()
    defects, actual_obs = sampler.sample(shots=shots, separate_observables=True)
    
    predicted_obs = matching.decode_batch(defects)
    errors = np.sum(actual_obs.flatten() != predicted_obs.flatten())
    return errors / shots


distances = [3, 5, 7]
physical_errors = np.linspace(0.005, 0.01, 20) 

plt.figure(figsize=(7, 4))
for d in distances:
# Simulating distance
    logical_errors = [get_logical_error(d, p) for p in physical_errors]
    plt.plot(physical_errors, logical_errors, label=f'd={d}')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Physical Error Rate ($P_{phys}$)')
plt.ylabel('Logical Error Rate ($P_L$)')
plt.grid(True, which="both", ls="-", alpha=0.5)
plt.legend()
plt.title("Surface Code Threshold Crossing")
plt.savefig(plotpath + '/p_threshold.pdf')