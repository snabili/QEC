import scipy
import stim
import pymatching
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from qiskit_ibm_runtime import QiskitRuntimeService
from qiskit_aer.noise import NoiseModel, thermal_relaxation_error, depolarizing_error, ReadoutError
from qiskit import QuantumCircuit
import math
import networkx as nx

import os, sys, argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
# importing custom codes
from qec import utils as ut


'''
This code's purpose is to find the right patch of qubits,
from an online backend, turn it to graphX object to find qubit's connection,
select logical qubits study
'''

parser = argparse.ArgumentParser(description="QEC IBM_FEZ")
# Add arguments
parser.add_argument("cent", type=int,  help="Central qubit", default=23)
args = parser.parse_args()


# Load a heavy-hex backend from 'ibm_fez'
service = QiskitRuntimeService()
backend = service.backend("ibm_fez")

# Convert the coupling map to a NetworkX graph for easier searching
cmap = backend.coupling_map
# Convert the edge list to a standard Python list, then to a Graph
# Take edges here or use the cmap object
edges = list(cmap.get_edges()) 
G = nx.Graph()
G.add_edges_from(edges)

'''
find 17 qubits centered at args.cent:
    to perform 3x3 logical qubit,
    save it to an object called found and the layout to a pdf file,
    printout qubits id in logical qubits:
        their roles, e.g., data & syndrome
        X and Z syndromes
        logical X and Z 
'''
found = ut.find_d3_patches(G)
patch_nodes = [i for i in found if i[0]==args.cent][0]
print([i[0] for i in found])

print(f"Found {len(found)} candidate patches, centered at {args.cent}: {patch_nodes}.")
if found:
    print(f"Patch centered at {patch_nodes[0]}: {sorted(patch_nodes)}")
    
ut.visualize_patch(G, args.cent)

# Qubit roles of: data + x-syndrome & z-syndrome
roles = ut.get_true_roles(G, patch_nodes)

print(f"Data Qubits: {roles['data']}")
print(f"X-Syndromes: {roles['x_syndrome']}")
print(f"Z-Syndromes: {roles['z_syndrome']}")

# find logical qubits path:
result = ut.find_logical_qubit_full(G, args.cent)

if result:
    print(f"Logical X (Physical X gates): {result['logical_x']}")
    print(f"Logical Z (Physical Z gates): {result['logical_z']}")
    print(f"Intersection Qubit: {result['pivot']}")
