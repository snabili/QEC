import pickle as pkl
import sys, os, argparse 
import numpy as np
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from qec import utils as ut

from qiskit.circuit import QuantumCircuit

# Add parser argument to build 17 qubits patch based on qubit id 
parser = argparse.ArgumentParser(description="QEC on IBM_FEZ machine")
parser.add_argument("qid", type=int,  help="qubit_id", default=23)
parser.add_argument("zsyn_mon",type=int, help="Z_Syndrome monitoring")
parser.add_argument("xsyn_mon",type=int, help="X_Syndrome monitoring", default=16)
args = parser.parse_args()

# open the graph object
with open("files/datafiles/graph.pkl", "rb") as f:
    G = pkl.load(f)

result = ut.find_logical_qubit_full(G, args.qid)
roles = ut.get_true_roles(G, result['patch'])

stabilizer_circuit = ut.build_heavy_hex_stabilizer(result, roles, G)
counts = ut.run_stabilizer(stabilizer_circuit)
# Check connectivity of result['pivot'] to see what are the Z + X Syndromes next to it
ut.check_qubit_connectivity(G, result['pivot'], roles)


# =============== Check for Z-Errors & X_Errors ===============
max_idx = max(result['patch']) + 1
# =============== X-Stabilizer: to detect Phase-Flip (Z) Errors ============
qc_z = QuantumCircuit(max_idx, 8)

# Initialize data states by reseting and changing basis from |0> to |+> (X-Pauli basis)
for d_q in roles['data']:
    qc_z.reset(d_q)
    qc_z.h(d_q)
qc_z.barrier()

# Inject Z-error on pivot
qc_z.z(result['pivot']) 
qc_z.barrier()

# Measure ONLY X-syndromes to avoid Z-measurement noise; map it to classical bits 2-7
for i, s_qubit in enumerate(roles['x_syndrome']):
    qc_z.reset(s_qubit)
    qc_z.h(s_qubit)
    for n in G.neighbors(s_qubit):
        if n in roles['data']:
            qc_z.cx(s_qubit, n) # CNOT for Phase-Flip errors: CNOT(A,D)
    qc_z.h(s_qubit)
    qc_z.measure(s_qubit, i + 2) 

print('X_stabilizer (Should show flips in bits next to X-Syndrome qubits):')
x_stabilizer = ut.run_stabilizer(qc_z)

# =============== Z-Stabilizer: to detect Bit-Flip (X) Errors ===============
qc_x = QuantumCircuit(max_idx, 8)

# Initialize data in |0> basis (Z-basis)
for d_q in roles['data']:
    qc_x.reset(d_q)
qc_x.barrier()

# Inject X-error on a qubit that IS monitored by Z-syndromes (e.g., Qubit 4)
# Note: Pivot might not have Z-syndrome neighbors in the current roles
qc_x.x(args.zsyn_mon) 
qc_x.barrier()

# Measure ONLY Z-syndromes (Bits 0-1)
for i, s_qubit in enumerate(roles['z_syndrome']):
    qc_x.reset(s_qubit)
    for n in G.neighbors(s_qubit):
        if n in roles['data']:
            qc_x.cx(n, s_qubit)
    qc_x.measure(s_qubit, i)

print(f'Z_stabilizer (Should show flip in bit for Qubit {args.zsyn_mon} error):')
z_stabilizer = ut.run_stabilizer(qc_x)