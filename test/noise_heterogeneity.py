import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pickle as pkl

from qiskit_ibm_runtime import QiskitRuntimeService

import stim
import pymatching

import os, sys, argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
# importing custom codes
from qec import utils, config

logpath = config.LOG_DIR
datapath = config.DATA_DIR

parser = argparse.ArgumentParser(description="QEC IBM_FEZ")
# Add arguments
parser.add_argument("dist", type=int,  help="logical qubit distance",   default=3)
args = parser.parse_args()

service = QiskitRuntimeService()
backend = service.backend("ibm_fez")

'''with open(datapath + "/graph.pkl", "rb") as f:
    G = pkl.load(f)'''

clean_G = utils.get_clean_subgraph(backend)

centers = utils.get_valid_centers(clean_G, args.dist)
best_c, stats = utils.get_best_patch_center(backend, clean_G, centers, args.dist)

found = utils.find_d3_patches(clean_G, args.dist)
patch_nodes = [i for i in found if i[0]==best_c][0]

logger = utils.setup_logging(log_path=logpath + "/noise_hetero.txt")

# Qubit roles of: data + x-syndrome & z-syndrome
roles = utils.get_true_roles(clean_G, patch_nodes, args.dist)
result = utils.find_logical_qubit_full(backend, clean_G, args.dist)


def get_actual_noise_model(result, roles, clean_G, backend):
    """
    Creates a dictionary of noise parameters mapped to physical qubit IDs.
    """
    props = backend.properties()
    noise_data = {}    
    # Extract specific T1, T2 and Readout for each qubit in the patch
    for q in result['patch']:
        t1 = props.t1(q)
        t2 = props.t2(q)
        readout = min(props.readout_error(q), 0.49)
        noise_data[q] = {'t1': t1, 't2': t2, 'readout': readout}
        
    # Extract specific CX gate errors
    gate_errors = {}
    sub_G = clean_G.subgraph(result['patch'])
    for u, v in sub_G.edges():
        try:
            err = props.gate_error('cz', [u, v])
            if err >= 0.75:
                err = 0.5
        except:
            err = 0.01 # Fallback if data is missing
        gate_errors[(u, v)] = err        
    return noise_data, gate_errors


def build_stim_circuit(result, roles, noise_mode="actual", noise_map=None, gate_map=None, avg_p=0.001, rounds=3):
    c = stim.Circuit()
    all_q = result['patch']
    q_to_idx = {q: i for i, q in enumerate(all_q)}

    num_syndromes = len(roles['x_syndrome'] + roles['z_syndrome'])

    for q in all_q:
        p = noise_map[q]['readout'] if noise_mode == "actual" else avg_p
        c.append("R", [q_to_idx[q]]) # Ideal Reset to |0>
        c.append("X_ERROR", [q_to_idx[q]], p) # Add the probability of a bit-flip during initialization

    for r in range(rounds):
        for (u, v), p in (gate_map.items() if noise_mode == "actual" else []):
            if u in q_to_idx and v in q_to_idx:
                c.append("CZ", [q_to_idx[u], q_to_idx[v]])
                c.append("DEPOLARIZE2", [q_to_idx[u], q_to_idx[v]], p)
        for i, q in enumerate(roles['x_syndrome'] + roles['z_syndrome']):
            p = noise_map[q]['readout'] if noise_mode == "actual" else avg_p
            c.append("M", [q_to_idx[q]], p)            
            if r == 0:
                c.append("DETECTOR", [stim.target_rec(-1)])
            else:
                c.append("DETECTOR", [stim.target_rec(-1), stim.target_rec(-(1 + num_syndromes))])

    for q in roles['data']:
        p = noise_map[q]['readout'] if noise_mode == "actual" else avg_p
        c.append("M", [q_to_idx[q]], p)

    logger.info(f"the data length is = {int(np.sqrt(len(roles['data'])))}")
    d = int(np.sqrt(len(roles['data']))) 
    data_rec_indices = [-(len(roles['data']) - i) for i in range(d)]    
    c.append("OBSERVABLE_INCLUDE", [stim.target_rec(i) for i in data_rec_indices], 0)        
    return c

# --- EXECUTION ---
noise_map, gate_map = get_actual_noise_model(result, roles, clean_G, backend)
avg_p = np.mean([m['readout'] for m in noise_map.values()])
rounds=3

circuit_actual  = build_stim_circuit(result, roles, "actual",   noise_map, gate_map, avg_p=avg_p,rounds=rounds)
circuit_avg     = build_stim_circuit(result, roles, "average",  noise_map, gate_map, avg_p=avg_p,rounds=rounds)

def get_logical_error(circuit,num_shots=100000):
    detector_sampler = circuit.compile_detector_sampler()
    defects, actual_logic_flips = detector_sampler.sample(shots=num_shots, separate_observables=True)
    
    matching = pymatching.Matching.from_stim_circuit(circuit)
    predictions = matching.decode_batch(defects)
    
    num_errors = np.sum(predictions.flatten() != actual_logic_flips.flatten())
    return num_errors / num_shots

p_logical_actual = get_logical_error(circuit_actual)
p_logical_avg = get_logical_error(circuit_avg)

print(f"Logical Error Rate (Actual Chip Noise): {p_logical_actual:.5f}")
print(f"Logical Error Rate (Average Noise):     {p_logical_avg:.5f}")


largest_cc_nodes = max(nx.connected_components(clean_G), key=len)
clean_G_sub = clean_G.subgraph(largest_cc_nodes).copy()


centers = utils.get_valid_centers(clean_G_sub, args.dist)
all_qubits_noise = {}
props = backend.properties()
for i in range(backend.num_qubits):
    try:
        all_qubits_noise[i] = {'readout': props.readout_error(i)}
    except:
        all_qubits_noise[i] = {'readout': 0.5} # Fallback for dead qubits


plotpath = config.PLOT_DIR
filename = plotpath + '/backendhealth_dist.pdf'
utils.plot_backend_health(clean_G, all_qubits_noise, filename)


# Getting statistics of qubits passing health survey:
logger.info(f"Original Qubits: {backend.num_qubits}; Clean Qubits: {clean_G.number_of_nodes()}")
edges = backend.coupling_map.get_edges()
undirected_edges = {tuple(sorted(edge)) for edge in edges}

logger.info(f"Total unique physical connections: {len(undirected_edges)}; Clean Edges: {clean_G.number_of_edges()}")
components = list(nx.connected_components(clean_G))
if components:
    largest_comp = len(max(components, key=len))
    logger.info(f"Largest connected block of qubits: {largest_comp}")
else:
    logger.info("Graph is empty!")