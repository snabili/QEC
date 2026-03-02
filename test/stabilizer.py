import pickle as pkl
import sys, os, argparse 
import numpy as np
from qiskit.circuit import QuantumCircuit
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from qec import utils, config

'''
This code intention is to verify the stabilizers
w/wo injecting any noises
'''

# Add parser argument to build 17 qubits patch based on qubit id 
parser = argparse.ArgumentParser(description="QEC on IBM_FEZ machine")
parser.add_argument("script",type=str,help="Name of script function to run (e.g. stabilizer, stabilizer_noises)")
parser.add_argument("qid", type=int,  help="qubit_id", default=23)
parser.add_argument("zsyn_mon",type=int, help="Z_Syndrome monitoring")
parser.add_argument("xsyn_mon",type=int, help="X_Syndrome monitoring", default=16)

args = parser.parse_args()

# file dirs
datapath = config.DATA_DIR
logpath = config.LOG_DIR
logger = utils.setup_logging(log_path=logpath + "/stabilizer.txt", name="Stabilizer")

# open graph object, to retrieve the entire heavy-hex
with open(datapath + "/graph.pkl", "rb") as f:
    G = pkl.load(f)

result = utils.find_logical_qubit_full(G, args.qid)
roles = utils.get_true_roles(G, result['patch'])

scripter = utils.Scripter()

@scripter
def stabilizer():
    '''
    Check the stabilizer without injecting any noise
    Check connectivity of result['pivot'] on Z + X adjacent Syndromes
    '''
    stabilizer_circuit = utils.build_heavy_hex_stabilizer(result, roles, G)
    counts = utils.run_stabilizer(stabilizer_circuit)
    utils.check_qubit_connectivity(G, result['pivot'], roles)

@scripter
def stabilizer_noise():
    '''
    Check for Z-Errors & X_Errors by adding noise manually
        - Allocate all the qubits in the patch --> might exceed patch size
        - Allocate 8 classical bits for syndrome measurements
        - X-Stabilizer to detect Phase-Flip (Z) Errors; 
            Measure ONLY X-syndromes to avoid Z-measurement noise; map it to classical bits 2-7
        - Z-Stabilizer: to detect Bit-Flip (X) Errors;
            Measure ONLY Z-syndromes (Bits 0-1)
    '''

    max_idx = max(result['patch']) + 1
    qc_z = QuantumCircuit(max_idx, 8)
    qc_x = QuantumCircuit(max_idx, 8)

    # Data-qubits: Reset, Change basis from |0> to |+> (to X-basis)
    for d_q in roles['data']:
        qc_z.reset(d_q)
        qc_z.h(d_q)
    qc_z.barrier()

    qc_z.z(result['pivot']) # Inject Z-error on pivot
    qc_z.barrier()

    for i, s_qubit in enumerate(roles['x_syndrome']):
        qc_z.reset(s_qubit)
        qc_z.h(s_qubit)
        for n in G.neighbors(s_qubit):
            if n in roles['data']:
                qc_z.cx(s_qubit, n) # CNOT Phase-Flip errors: CNOT(A,D)
        qc_z.h(s_qubit)
        qc_z.measure(s_qubit, i + 2) 

    logger.info('X_stabilizer (Should show flips in bits next to X-Syndrome qubits):')
    x_stabilizer = utils.run_stabilizer(qc_z)

    # Initialize data in |0> basis (Z-basis)
    for d_q in roles['data']:
        qc_x.reset(d_q)
    qc_x.barrier()

    # Inject X-error on a qubit that IS monitored by Z-syndromes (e.g., Qubit 4)
    # Note: Pivot might not have Z-syndrome neighbors in the current roles
    qc_x.x(args.zsyn_mon) 
    qc_x.barrier()

    for i, s_qubit in enumerate(roles['z_syndrome']):
        qc_x.reset(s_qubit)
        for n in G.neighbors(s_qubit):
            if n in roles['data']:
                qc_x.cx(n, s_qubit)
        qc_x.measure(s_qubit, i)
    logger.info(f'Z_stabilizer (Should show flip in bit for Qubit {args.zsyn_mon} error):')
    z_stabilizer = utils.run_stabilizer(qc_x)

if __name__ == '__main__':
        scripter.run()