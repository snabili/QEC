# Iport general packages
import os, os.path as osp, logging, re, time, json
from contextlib import contextmanager
import logging, warnings
import sys, argparse

# plotting packages
import matplotlib.pyplot as plt
import networkx as nx

# IBM Qiskit:
from qiskit.circuit import QuantumCircuit
from qiskit_aer import AerSimulator


'''
This code's purpose is to gather functions required to run other codes
'''

def setup_logging(name='thermo_logger', log_path=None, level=logging.INFO, console=True):
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

def find_d3_patches(graph):
    '''
    In heavy-hex, logical patches are centered on 'Data' qubits
    Look for nodes with degree 3 as potential centers, with 17-qubit (9 data + 8 syndrome) diamond
    '''
    patches = []
    centers = [n for n, d in graph.degree() if d == 3]
    logger.info('Central qubits {centers}')    
    for c in centers:
        nodes_at_dist = nx.single_source_shortest_path_length(graph, c, cutoff=4) # Need to go 4 steps out to encompass the 17-qubit diamond
        subnodes = list(nodes_at_dist.keys())       
        if len(subnodes) >= 17: # satisfy the minimum size.
            # Sort by distance to ensure getting 17-qubits closest to the center 'c'
            sorted_nodes = sorted(subnodes, key=lambda x: nodes_at_dist[x])
            patch = sorted_nodes[:17]
            patches.append(patch)        
    return patches



def get_true_roles(full_graph, patch_nodes):
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
    coloring = nx.bipartite.color(sub) # heavy-hex lattice is Bipartite; color it into two sets.
    
    set_0 = [n for n, color in coloring.items() if color == 0]
    set_1 = [n for n, color in coloring.items() if color == 1]
    
    if len(set_0) == 9:
        data_q, syndrome_q = set_0, set_1
    else:
        data_q, syndrome_q = set_1, set_0
    
    x_syn = [n for n in syndrome_q if full_graph.degree(n) == 3]
    z_syn = [n for n in syndrome_q if full_graph.degree(n) == 2]
    
    return {"data": sorted(data_q), "x_syndrome": sorted(x_syn), "z_syndrome": sorted(z_syn)}


def find_logical_qubit_full(G, center_node):
    '''
    Find logical qubit:
        - Identify 17-qubits cluster (d=3) subgraph around central-node out of all qubits 
        - Assign roles using Bipartite graph function
        - Find logical strings (node=5, edges=4) candidates that satisfy: D-S-D-S-D (D:data, S:syndrome)
        - Find logical strings that intersect at exactly one node
    Return:
        - the nodes in the patch
        - data nodes
        - logical-x (vertical 5 nodes path)
        - logical-z (horizontal 5 nodes path)
    '''
    nodes_at_dist = nx.single_source_shortest_path_length(G, center_node, cutoff=4)
    patch_nodes = sorted(nodes_at_dist.keys(), key=lambda x: nodes_at_dist[x])[:17]
    sub = G.subgraph(patch_nodes)

    coloring = nx.bipartite.color(sub)
    set_0 = [n for n, color in coloring.items() if color == 0]
    set_1 = [n for n, color in coloring.items() if color == 1]
    data_nodes = set_0 if len(set_0) == 9 else set_1
    
    logical_strings = []
    for start in data_nodes:
        for end in data_nodes:
            if start >= end: continue
            for path in nx.all_simple_paths(sub, start, end, cutoff=4):
                if len(path) == 5 and all(path[i] in data_nodes for i in [0, 2, 4]):
                    logical_strings.append(path)

    for i, path_a in enumerate(logical_strings):
        for path_b in logical_strings[i+1:]:
            data_a = {path_a[0], path_a[2], path_a[4]}
            data_b = {path_b[0], path_b[2], path_b[4]}
            intersection = data_a.intersection(data_b)
            
            if len(intersection) == 1:
                return {
                    "patch": patch_nodes,
                    "data": sorted(list(data_nodes)),
                    "logical_x": [path_a[0], path_a[2], path_a[4]],
                    "logical_z": [path_b[0], path_b[2], path_b[4]],
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

