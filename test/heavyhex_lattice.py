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

# importing custom codes
import utils as ut


parser = argparse.ArgumentParser(description="QEC IBM_FEZ")
# Add arguments
parser.add_argument("--cent", type=int,  help="Central qubit", default=23)
args = parser.parse_args()
print(args.cent)

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


# save the 17-qubits to perform 3x3 logical qubit to an object called found
found = ut.find_d3_patches(G)
# find 17 qubits centered at args.cent 
patch_nodes = [i for i in found if i[0]==args.cent][0]
print([i[0] for i in found])# if i[0]==args.cent])

print(f"Found {len(found)} candidate patches, centered at {args.cent}: {patch_nodes}.")
if found:
    print(f"Patch centered at {patch_nodes}: {sorted(patch_nodes)}")
    #print(f"Patch centered at {found[0][0]}: {sorted(found[0])}")

# Try visualizing with node args.cent 
ut.visualize_patch(G, args.cent)

# Get the roles of each qubit: data + x-syndrome & z-syndrome
#patch_nodes = found[0]
#patch_nodes = [i for i in found if i[0]==args.cent]
roles = ut.get_true_roles(G, patch_nodes)

print(f"Data Qubits: {roles['data']}")
print(f"X-Syndromes: {roles['x_syndrome']}")
print(f"Z-Syndromes: {roles['z_syndrome']}")

# find logical qubits path:
result = ut.find_logical_qubit_full(G, args.cent)

if result:
    print(f"Logical X String (Physical X gates): {result['logical_x']}")
    print(f"Logical Z String (Physical Z gates): {result['logical_z']}")
    print(f"Intersection Qubit: {result['pivot']}")

# Drawing logical patch with syndrome identification:
# Execute the drawing
ut.draw_logical_qubit(G, result)