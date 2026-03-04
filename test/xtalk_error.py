from qiskit_ibm_runtime import QiskitRuntimeService
from qiskit_experiments.library import StandardRB
from qiskit_experiments.framework import BatchExperiment

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import networkx as nx

import os, sys, argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from qec import config, utils

logpath = config.LOG_DIR
plotpath = config.PLOT_DIR

logger = utils.setup_logging(log_path=logpath + "/xtalk.txt")

"""
This code intention is to estimate the crosstalk probability.
As the IBM account used for this project is open-instance,
the frequency/qubit is not available ==> used gate_error as a proxy of frequency
Steps to do the study after accessing the backend:
    - Extract CZ Gate Errors, CX (CNOT) is not native for this machine
    - Build a Connectivity Graph
    - Visualize the "Error Map" by create a Log norm object ==> not neccessary, just for fun!
    - Find the best patch with lowest CZ error to let frequency related error to show itself
    - Running Randomize Benchmark in two scenarios:
        - Isolated Experiment (Only running RB on target)
        - Simultaneous Experiment (Running RB on target AND spectator at once)
            to forces the hardware to deal with concurrent pulses, triggering crosstalk
    - Finding target and spectator qubits
    - Shift from math to physical access to IBM machines
    - Compute deltaX and print the results
"""

# 1. Initialize the service and get the backend
service = QiskitRuntimeService()
backend = service.backend("ibm_fez")
props = backend.properties()

# 2. E
cz_errors = {}
for gate in props.gates:
    if gate.gate == 'cz':
        # Extract the error value (usually the first parameter)
        qubits = tuple(gate.qubits)
        error = gate.parameters[0].value
        cz_errors[qubits] = error

# 3. 
G = nx.Graph()
for (u, v), err in cz_errors.items():
    # Use 1/error for weighting so that high error = "longer" or "thinner" edge
    G.add_edge(u, v, weight=err)

# **********************  Plotting  ******************************************** 
pos = nx.spring_layout(G, seed=42) 
edges, weights = zip(*cz_errors.items())


vmin = min(cz_errors.values())
vmax = 0.01 # avoid using max(cz_errors.values()) as the broken a qubit lead to vmax=1 and z-axis is too coarse
norm = colors.LogNorm(vmin=vmin, vmax=vmax)


fig, ax = plt.subplots(figsize=(6, 4))
nodes = nx.draw_networkx_nodes(G, pos, node_color='skyblue', node_size=100)
edges = nx.draw_networkx_edges(
    G, pos, 
    edgelist=list(cz_errors.keys()), 
    edge_color=list(cz_errors.values()),
    edge_cmap=plt.cm.YlOrRd, # Yellow to Red (Higher error = Redder),
    edge_vmin=vmin, # Set min for color mapping
    edge_vmax=vmax, # Set max for color mapping
    width=5,
    ax=ax
)

sm = plt.cm.ScalarMappable(cmap=plt.cm.YlOrRd, norm=norm)
sm.set_array(list(cz_errors.values()))

plt.colorbar(sm, ax=ax, label='CZ Gate Error Rate (Log Scale)')
plt.title(f"Crosstalk Proxy Map for {backend.name}")
plt.savefig(plotpath + '/xtalk.pdf')

# **********************  Find Quiet Patch  ******************************************** 
best_patch = utils.find_best_subgraph(cz_errors)
logger.info(f"Top 5 'Quiet' Qubits for Investigation:, {[q[0] for q in best_patch]}")
quiet_qubits = [q[0] for q in best_patch]

def find_quiet_connected_pair(quiet_list, cz_errors):
    for i, q1 in enumerate(quiet_list):
        for q2 in quiet_list[i+1:]:           
            if (q1, q2) in cz_errors or (q2, q1) in cz_errors: # check if CZ link exists between q1, q2
                return (q1, q2)
    return None

target_pair = find_quiet_connected_pair(quiet_qubits, cz_errors)
logger.info(f"Target Pair for study: {target_pair}")

def find_spectator(target_pair, backend):
    
    coupling_map = backend.coupling_map # get the coupling map, qubits are physically neighbors
    q1, q2 = target_pair        
    for edge in coupling_map: # look for neighbors of q1 or q2 that are NOT q1 or q2
        u, v = edge
        if u == q1 and v != q2:
            return v
        if u == q2 and v != q1:
            return v
    return None

spectator = find_spectator(target_pair, backend)
logger.info(f"Identified Spectator: {spectator}")

# **********************  Run RB on isolated & simultaneous  ******************************************** 
exp_iso = StandardRB(target_pair, lengths=[1, 10, 20, 50, 100, 150])

exp_spectator = StandardRB([spectator], lengths=[1, 10, 20, 50, 100, 150])
exp_sim = BatchExperiment([exp_iso, exp_spectator], flatten_results=True)

# this is the most time consuming part
logger.info("Transition from mathematical definitions to physical execution on the IBM quantum hardware.")
data_iso = exp_iso.run(backend).block_for_results()
data_sim = exp_sim.run(backend).block_for_results()

results = data_sim.analysis_results()

EPC_sim = next(res.value.n for res in results if res.name == 'EPC' and len(res.device_components) == 2)
Quality_sim = next(res.quality for res in results if res.name == 'EPC' and len(res.device_components) == 2)


EPC_iso = data_iso.analysis_results("EPC").value.n
Quality_iso = data_iso.analysis_results("EPC").quality

delta_xt = EPC_sim - EPC_iso

logger.info(f"Isolated Error:     {EPC_iso:.5f}, with quality: {Quality_iso:s}")
logger.info(f"Simultaneous Error: {EPC_sim:.5f}, with quality: {Quality_sim:s}")
logger.info(f"Crosstalk Delta:    {delta_xt:.5f}")
