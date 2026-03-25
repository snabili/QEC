# Iport general packages
import os, os.path as osp, logging, re, time, json
from contextlib import contextmanager
import logging, warnings
import sys, argparse

# plotting packages
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

from scipy.integrate import quad # to compute integral; used in leakage; find_physical_amplitude() function

# IBM Qiskit:
from qiskit.circuit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit_dynamics import Solver, Signal


'''
This code's purpose is to gather functions required to run other codes
'''

def setup_logging(name='qec_logger', log_path=None, level=logging.INFO, console=True):
    """
    Create a reusable logger with both file and console outputs.
    Args:
        name (str): Logger name.
        log_path (str): Path to the log file.
        level (int): Logging level (e.g., logging.INFO, logging.DEBUG).
        console (bool): Whether to log to stdout.
    
    Returns:
        logging.Logger: Configured logger object.
    """
    # Suppress noisy warnings
    warnings.filterwarnings("ignore", category=FutureWarning)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False  # Avoid double logging

    formatter = logging.Formatter(
        "%(asctime)s - [%(levelname)s] - %(name)s - %(message)s", 
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    if logger.hasHandlers():
        logger.handlers.clear()
    # Avoid duplicate handlers
    if not logger.handlers:
        if log_path:
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            file_handler = logging.FileHandler(log_path)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

        if console:
            stream_handler = logging.StreamHandler()
            stream_handler.setFormatter(formatter)
            logger.addHandler(stream_handler)

    logger.info("Logger initialized successfully.")
    return logger

logger = setup_logging()

def pull_arg(*args, **kwargs):
    """
    Pulls specific arguments out of sys.argv.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(*args, **kwargs)
    args, other_args = parser.parse_known_args()
    sys.argv = [sys.argv[0]] + other_args
    return args

def read_arg(*args, **kwargs):
    """
    Reads specific arguments from sys.argv but does not modify sys.argv
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(*args, **kwargs)
    args, _ = parser.parse_known_args()
    return args

#from contextlib import contextmanager
# decorator as command dispatcher
#@contextmanager
class Scripter: # --> scripter decorator
    def __init__(self):
        self.scripts = {}

    def __call__(self, fn):
        self.scripts[fn.__name__] = fn
        return fn

    def run(self):
        script = pull_arg('script', choices=list(self.scripts.keys())).script
        #setup_logging().info('Running %s', script)
        logger.info('Running %s', script)
        self.scripts[script]()

@contextmanager
def time_and_log(begin_msg, end_msg='Done'):
    try:
        t1 = time.time()
        logger.info(begin_msg)
        yield None
    finally:
        t2 = time.time()
        nsecs = t2-t1
        nmins = int(nsecs//60)
        nsecs %= 60
        logger.info(end_msg + f' (took {nmins:02d}m:{nsecs:.2f}s)')

def get_valid_centers(G, dist):
    """
    Finds all qubits that can act as the center for a patch of a given distance.
    In heavy-hex, logical centers should be 'Data' qubits, and they usually have degree 2 or 3.
    By check if the graph can provide enough qubits within the radius
    """
    valid_centers = []
    min_size = 2 * dist**2 - 1
    search_radius = dist + (dist // 2) # Heuristic for heavy-hex stretching

    for node in G.nodes():
        nodes_nearby = nx.single_source_shortest_path_length(G, node, cutoff=search_radius)        
        if len(nodes_nearby) >= min_size:
            valid_centers.append(node)            
    return valid_centers



def get_clean_subgraph(backend, error_threshold=0.15):
    """
    Returns a NetworkX graph of the backend where 'dead' or 
    high-error qubits/gates are removed.
    """
    props = backend.properties()
    G = nx.Graph()
    
    # 1. Filter Qubits (Nodes)
    valid_qubits = []
    for i in range(backend.num_qubits):
        try:
            if props.readout_error(i) < error_threshold:
                valid_qubits.append(i)
        except:
            continue
            
    # 2. Filter Gates (Edges)
    for edge in backend.coupling_map:
        u, v = edge
        if u in valid_qubits and v in valid_qubits:
            try:
                # IBM Fez uses CZ; check the gate error
                err = props.gate_error('cz', [u, v])
                if err < error_threshold:
                    G.add_edge(u, v, weight=err)
            except:
                continue
    return G

def find_d3_patches(graph, dist):
    '''
    In heavy-hex, logical patches are centered on 'Data' qubits
    Look for nodes with degree 3 as potential centers, with 17-qubit (9 data + 8 syndrome) diamond
    '''
    patches = []
    centers = get_valid_centers(graph, dist)
    logger.info(f'Central qubits {centers} with distance={dist}') 
    min_size, search_radius = 2 * dist**2 - 1, 2 * (dist-1) + 2
    for c in centers:
        nodes_at_dist = nx.single_source_shortest_path_length(graph, c, cutoff=search_radius)
        subnodes = list(nodes_at_dist.keys())
        if len(subnodes) >= min_size:
            sorted_nodes = sorted(subnodes, key=lambda x: nodes_at_dist[x])
            patch = sorted_nodes[:min_size]
            patches.append(patch)        
    return patches


def get_best_patch_center(backend, G, centers, dist):
    """
    Evaluates candidate centers and returns the one with the healthiest 
    average calibration metrics (T1, T2, and Gate Error),
    Scoring Criteria:
        high T1 + T2 and LOW gate error. 
        Score = (Avg T1 + Avg T2) / Avg Gate Error
    """
    best_center = None
    best_score = -np.inf
    patch_stats = []

    min_size = 2 * dist**2 - 1
    search_radius = dist + (dist // 2)

    properties = backend.properties()

    for c in centers:
        nodes_at_dist = nx.single_source_shortest_path_length(G, c, cutoff=search_radius)
        patch_nodes = sorted(nodes_at_dist.keys(), key=lambda x: nodes_at_dist[x])[:min_size]
        
        t1s = []
        t2s = []
        for q in patch_nodes:
            t1s.append(properties.t1(q))
            t2s.append(properties.t2(q))
        
        avg_t1 = np.mean(t1s)
        avg_t2 = np.mean(t2s)

        # Gate Metrics (CZ Error)
        gate_errors = []
        subgraph_edges = G.subgraph(patch_nodes).edges()
        for u, v in subgraph_edges:
            error = properties.gate_error('cz', [u, v])
            gate_errors.append(error)
        
        avg_gate_error = np.mean(gate_errors) if gate_errors else 1.0

        score = (avg_t1 + avg_t2) / avg_gate_error

        patch_stats.append({
            "center": c,
            "score": score,
            "avg_t1": avg_t1,
            "avg_t2": avg_t2,
            "avg_gate_error": avg_gate_error
        })
        if score > best_score:
            best_score = score
            best_center = c

    return best_center, patch_stats

def get_true_roles(full_graph, patch_nodes, dist):
    '''
    Heavy-hex layout for d=3: N(data-qubits)=d^2=9, N(ancilla-qubits)=d^2-1=8
    This function find data + ancilla qubits in a subgraph:
        - data-qubits must connect to ancilla (syndrome)
        - X-syndrome, stars, connected to 3 data qubits
        - Z-syndrome, bridges, connected to 2 data qubits
    The function uses Bipartite network to categorize data & syndrome qubits
    It returns sorted:
        data, x_syndrome, z_syndrome
    '''
    sub = full_graph.subgraph(patch_nodes)    
    try:
        coloring = nx.bipartite.color(sub)
    except nx.NetworkXError:
        # If the patch isn't bipartite, it's not a valid QEC patch
        return {"data": [], "x_syndrome": [], "z_syndrome": []}
    
    set_0 = [n for n, color in coloring.items() if color == 0]
    set_1 = [n for n, color in coloring.items() if color == 1]
    
    target_data_count = dist**2
    if len(set_0) == target_data_count:
        data_q, syndrome_q = set_0, set_1
    elif len(set_1) == target_data_count:
        data_q, syndrome_q = set_1, set_0
    else:
        # Fallback: The set with more nodes is usually the data set in these patches
        data_q, syndrome_q = (set_0, set_1) if len(set_0) > len(set_1) else (set_1, set_0)
    x_syn = []
    z_syn = []
    
    for n in syndrome_q:
        # X-syndromes (Stars) connect to more data qubits than Z-syndromes (Bridges)
        # In heavy-hex d=3, X-syndromes usually have degree 3 in the subgraph
        if sub.degree(n) >= dist:
            x_syn.append(n)
        else:
            z_syn.append(n)
    
    return {"data": sorted(data_q), "x_syndrome": sorted(x_syn), "z_syndrome": sorted(z_syn)}


def find_logical_qubit_full(backend, G, dist):
    '''
    Find logical qubit with distance d
        - Identify (2d^2-1)-qubits cluster subgraph around central-node out of all qubits 
        - Assign roles using Bipartite graph function
        - Find logical strings (node=2x(d-1)+2, edges=node-1) candidates 
                that satisfy for example d = 3: D-S-D-S-D (D:data, S:syndrome)
        - Find logical strings that intersect at exactly one node
    Return:
        - the nodes in the patch
        - data nodes
        - logical-x (vertical nodes path)
        - logical-z (horizontal nodes path)
    '''

    min_size        = 2 * dist**2 - 1
    search_radius   = 2 * (dist-1) + 2
    string_len      = 2 * dist - 1
    centers         = get_valid_centers(G, dist)

    logger.info(f'centers: {centers}')

    center_node = get_best_patch_center(backend, G, centers, dist)
    logger.info(f'best patch center: {center_node[0]}')

    nodes_at_dist = nx.single_source_shortest_path_length(G, center_node[0], cutoff = search_radius)
    if len(nodes_at_dist) < min_size:
        return None
    patch_nodes = sorted(nodes_at_dist.keys(), key=lambda x: nodes_at_dist[x])[:min_size]
    sub = G.subgraph(patch_nodes)

    # Assign Roles (Scaling Data Node count)
    try:
        coloring = nx.bipartite.color(sub)
    except nx.NetworkXError:
        return None # Not a valid bipartite QEC patch
    
    set_0 = [n for n, color in coloring.items() if color == 0]
    set_1 = [n for n, color in coloring.items() if color == 1]
    target_data_count = dist**2
    data_nodes = set_0 if len(set_0) == target_data_count else set_1
    
    logical_strings = []
    for start in data_nodes:
        for end in data_nodes:
            if start >= end: continue
            for path in nx.all_simple_paths(sub, start, end, cutoff=string_len - 1):
                if len(path) == string_len:
                    # Every second node must be a data node (indices 0, 2, 4, 6...)
                    if all(path[i] in data_nodes for i in range(0, string_len, 2)):
                        logical_strings.append(path)

    for i, path_a in enumerate(logical_strings):
        for path_b in logical_strings[i+1:]:
            data_a = set(path_a[i] for i in range(0, string_len, 2))
            data_b = set(path_b[i] for i in range(0, string_len, 2))
            intersection = data_a.intersection(data_b)
            
            if len(intersection) == 1:
                return {
                    "patch": patch_nodes,
                    "data": sorted(list(data_nodes)),
                    "logical_x": sorted(list(data_a)),
                    "logical_z": sorted(list(data_b)),
                    "pivot": list(intersection)[0]
                }
    return None

# ================== Stabilizer ===================
def build_heavy_hex_stabilizer(result, roles,G):
    '''
    Building quantum circuit and stabilizer. For the stabilizer to work:
        - Caveat1: Data-Qubits outside the patch must be reset to avoid their noise effect on the patch
        - Caveat2: Data-Qubits inside the patch must NOT be reset to avoid erasing the encoded information
        - Caveat3: Syndrome-Qubits must be reset for the stabilizer algorithm to work
    Steps to build stabilizer:
        - Extract Data, X-Syndrome, Z-syndrome qubits
        - Create circuit with 8 classical bits + reset the qubits to |0> state
        - Z-syndrome measurement to find X-error (Bit-Flip error):
            1 - Loop over Z-syndrome qubits and reset them to avoid noise (to |0>)
            2 - Find data qubits attached to Z-syndromes
            3 - Apply CNOT(D,A) --> control:data target:ancillas(syndrome)
            4 - Measure Syndrome and write it to the corresponding classical bit (bit_idx)
        - X-Syndrome measurement to find Z-error (Phase-Flip error):
            1 - Apply Hadamart gate to data qubits to change basis: |0> (Z-basis) to |+> (X-basis)
            2 - Loop over X-syndrome qubits and reset them to avoid noise (to |0>)
            3 - Apply Hadamart to X-syndrome
            4 - Find data qubits attached to X-syndromes
            5 - Apply CNOT(A,D) --> control:ancillas(syndrome) target:data
            6 - Apply Hadamart again: |+> (X-basis) to |0> (Z-basis) as measurement is in Z-basis
            7 - Measure Syndrome and write it to the corresponding classical bit (bit_idx)
    '''
    data_qubits = roles['data']
    x_syns = roles['x_syndrome']
    z_syns = roles['z_syndrome']
    bit_idx = 0 # assigns the result of the measurement to the corresponding classical bit

    qc = QuantumCircuit(127, 8)  
    for q in range(127): 
        if q not in result['patch']: # should it not reset dataqubits in the patch as well?
            qc.reset(q)

    # ================ Z-Syndrome Measurements to detect X-errors (Bit-Flip Errors) ================
    res = []
    logger.info('Z-syndrom qubits --> neighbors of Z-syndrome: ')
    for s_qubit in z_syns:
        qc.reset(s_qubit)
        neighbors = [n for n in G.neighbors(s_qubit) if n in data_qubits]
        res.append(f"{s_qubit} -> {neighbors}")
        for data_q in neighbors:
            qc.cx(data_q, s_qubit) # CNOT(D,A) control:data, target:ancillas(syndrome)
        qc.measure(s_qubit, bit_idx)
        bit_idx += 1 
    logger.info(res)  
    qc.barrier()

    # ================ X-Syndrome Measurements to detect Z-errors (Phase-Flip Errors) ================
    for d_q in data_qubits:
        qc.h(d_q) # change |0> to |+>
    qc.barrier()
    res = []
    logger.info('X-syndrom qubits --> neighbors of X-syndrome: ')
    for s_qubit in x_syns:
        qc.reset(s_qubit)
        qc.h(s_qubit) # Hadamart gate to change |0> to |+> for syndrome qubits
        neighbors = [n for n in G.neighbors(s_qubit) if n in data_qubits]
        res.append(f"{s_qubit} -> {neighbors}")
        for data_q in neighbors:
            qc.cx(s_qubit, data_q) # CNOT(A,D) control:ancillas(syndrome), target:data --> Kickback Phase
        qc.h(s_qubit) # 2nd Hadamart to change |+> to |0>
        qc.measure(s_qubit, bit_idx)
        bit_idx += 1
    logger.info(res)
    return qc


def check_qubit_connectivity(G, qubit_id, roles):
    '''
    Printout:
        - Connectivity of central-qubit with X and Z syndromes
        - Report the bit-wise location of the flipped (bit & phase) qubits
    '''
    neighbors = list(G.neighbors(qubit_id))
    
    z_checkers = [n for n in neighbors if n in roles['z_syndrome']]
    x_checkers = [n for n in neighbors if n in roles['x_syndrome']]
    
    logger.info(f"--- Pivot Qubit {qubit_id} Analysis ---")
    logger.info(f"Connected Z-Syndromes: {z_checkers}")
    logger.info(f"Connected X-Syndromes: {x_checkers}")
    
    bits = roles['z_syndrome'] + roles['x_syndrome']
    x_bits = [bits.index(s) for s in x_checkers]
    z_bits = [bits.index(s) for s in z_checkers]
    
    logger.info(f"Expected flipped bits in output string: {z_bits + x_bits}")


def run_stabilizer(qc):
    '''
    This function executes the stabilizer quantum circuit on a simulator
    by running it 1024 shots; 
    return syndrome bitstring distribution
    '''
    simulator = AerSimulator(method='statevector') # Initialize a clean simulator without a target backend
    job = simulator.run(qc, shots=1024)
    result = job.result()
    
    counts = result.get_counts()
    print(f"Syndrome Outcomes: {counts}")    
    return counts


def plot_backend_health(G, global_noise_map, filename_path):
    plt.figure(figsize=(12, 8))
    pos = nx.spring_layout(G, seed=42)
    
    # This will now find every node because global_noise_map has all IDs
    readout_errors = [global_noise_map[node]['readout'] for node in G.nodes()]
    
    nodes = nx.draw_networkx_nodes(G, pos, 
                                   node_color=readout_errors,
                                   cmap=plt.cm.YlOrRd, 
                                   node_size=300)
    nx.draw_networkx_edges(G, pos, alpha=0.3)
    nx.draw_networkx_labels(G, pos, font_size=8)
    plt.colorbar(nodes, label='Readout Error Rate')
    
    plt.savefig(filename_path)


def find_best_subgraph(cz_errors, n_qubits=10):
    """
    ranks qubits by their 'Neighborhood Error' (average of all connected CZ links).
    useful function to do the crosstalk error approximation
    as the open-instance ibm account does not provide the frequency map
    """
    qubit_scores = {}
    # Get unique qubits from the links
    all_qubits = set([q for link in cz_errors.keys() for q in link])
    
    for q in all_qubits:
        # Find all CZ links connected to this qubit
        connected_errors = [err for link, err in cz_errors.items() if q in link]
        if connected_errors:
            qubit_scores[q] = np.mean(connected_errors)
    
    # Sort qubits: Lowest average error first
    sorted_qubits = sorted(qubit_scores.items(), key=lambda x: x[1], reverse=False)
    return sorted_qubits[:n_qubits]

# =========================== Extra/Test functions to test the workflow ===========================
def visualize_patch(graph, center_node):
    # 1. Get the 17-qubit patch using BFS (4 steps out)
    nodes_at_dist = nx.single_source_shortest_path_length(graph, center_node, cutoff=4)
    patch_nodes = sorted(nodes_at_dist.keys(), key=lambda x: nodes_at_dist[x])[:17]
    
    # 2. Create subgraph
    subgraph = graph.subgraph(patch_nodes)
    
    # 3. Define layout (spring_layout helps reveal the lattice structure)
    pos = nx.spring_layout(subgraph, seed=42) 
    
    plt.figure(figsize=(5, 3))
    
    # Color nodes: Center is Red, others are Blue
    colors = ['red' if n == center_node else 'skyblue' for n in subgraph.nodes()]
    
    nx.draw(subgraph, pos, with_labels=True, node_color=colors, 
            node_size=300, font_weight='bold', edge_color='gray')
    
    plt.title(f"17-Qubit Logical Patch Centered at {center_node}")
    plt.savefig('files/plots/logicpatch.pdf')

def build_test_stabilizer(result, roles, G):
    '''
    Generates the stabilizer for a specific patch of qubit nodes
    1 - Save dataqubits, X-syndrome and Z-syndrome qubit nodes
    2 - Check Z-syndrome by bit-flipping data qubits next to the syndrome followed by a measurement
    3 - Run the stabilizer test
    '''
    data_qubits = roles['data'] # [0, 2, 4, 6, 16, 22, 24]
    x_syns = roles['x_syndrome']
    z_syns = roles['z_syndrome'] # [1,5]
    
    # Create circuit: 127 qubits, 8 classical bits
    qc = QuantumCircuit(127, 8)    
    
    '''qc.x([2,4]) # bit flip data qubits next to the z-syndromes: 2 and/or 4
    qc.barrier()'''
    # bit_idx tracks which classical bit we are writing to
    bit_idx = 0
    # --- 1. Z-Syndrome Measurements (ZZ) ---
    for s_qubit in z_syns:
        neighbors = [n for n in G.neighbors(s_qubit) if n in data_qubits] # [2,4]
        print(f'Z-syndrom qubits: {s_qubit}, neighbors of Z-syndrome: {neighbors}')
        for data_q in neighbors:
            qc.cx(data_q, s_qubit) # apply CX: data is control, ancillas is target        
        qc.measure(s_qubit, bit_idx)
        bit_idx += 1
        qc.barrier()
    return qc


def find_spectator(pair, backend):
    '''
    to find the neighboor to the target_pair
    useful in comuting crosstalk error
    '''
    q1, q2 = pair
    for a, b in backend.configuration().coupling_map:
        if a == q1 and b != q2: return b
        if b == q1 and a != q2: return a
        if a == q2 and b != q1: return b
        if b == q2 and a != q1: return a
    return None


# ---------------------------------------------------------
# -------- Leakage Functions ------------------------------
# ---------------------------------------------------------

def find_physical_amplitude(target_theta, sigma, duration):
    '''
    To find MW pulse amplitude
    Used in leakage error estimation
    Computes amplitude of MW pulse intended 
        for a specific flip via target_theta
    MW pulse shape used is a gaussian --> could be extended for other MW shapes
    '''    
    unit_area, _ = quad(lambda t: np.exp(-(t - duration/2)**2 / (2 * sigma**2))
                        , 0, duration)  # Calculate the unitless area of the Gaussian curve 
    # Target rotation theta = Integral of (2 * pi * A * gauss)
    # A = theta / (2 * pi * unit_area)
    amp = target_theta / (2 * np.pi * unit_area)
    return amp


def sweep_envelope(amp, d, beta):
    sigma = d / 6
    def envelope(t):
        gauss = amp * np.exp(-(t - d/2)**2 / (2 * sigma**2))
        '''
        Derivative of the Gaussian (The "Q" component)
            d/dt exp(-t^2) = -2t * exp(-t^2)
        '''
        deriv = -(t - d/2) / (sigma**2) * gauss
        return gauss + 1j * beta * deriv
    return envelope

def find_best_beta(beta_values, amp, duration, v01, solver):    
    beta_leakages = []
    for beta in beta_values:
        env = sweep_envelope(amp, duration, beta)
        signal = Signal(envelope=env, carrier_freq=v01)
        sol = solver.solve(t_span=[0, duration], y0=np.array([1,0,0], dtype=complex), signals=[signal])    
        leakage = np.max(np.abs(sol.y[:, 2])**2)
        beta_leakages.append(leakage)
    return beta_leakages


