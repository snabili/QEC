from qiskit_ibm_runtime import QiskitRuntimeService

import networkx as nx 
import pickle 


# 1. Load heavy-hex backend from 'ibm_fez'
service = QiskitRuntimeService()
backend = service.backend("ibm_fez")
# 2. Coupling map --> NetworkX graph for easier searching
cmap = backend.coupling_map
# 1. Edge list --> standard Python list --> Graph
edges = list(cmap.get_edges()) 
G = nx.Graph()
G.add_edges_from(edges)
# 3. Saving the graph objects
with open("files/datafiles/graph.pkl", "wb") as f: 
    pickle.dump(G, f)

'''
# Alternatively, use downloaded csv file to sketch the qubits layout in heavy-hex
# Extracted edges from 2-qubit gates
# Caveat: the resulted graph might be different than the online backend results

# 1. Load CSV
csv_path = "files/datafiles/ibm_fez_calibrations_2026-01-06T18_10_35Z.csv"
df = pd.read_csv(csv_path)
# 2. Extract connected qubits
edges_csv = []
for qubit, entry in df['CZ error'].items():
    pairs = entry.split(';')
    for p in pairs:
        target, _ = p.split(':')
        edges_csv.append([qubit, int(target)])
edges_csv = {tuple(sorted(e))for e in edges_csv}
# 3. Tuples --> List
coupling_map = [list(e) for e in edges_csv]
# 4. Map --> Graph
G.add_edges_from(coupling_map)
# 5. Saving graph to pkl file:
with open("files/datafiles/graph.pkl", "wb") as f: 
    pickle.dump(G, f)
'''
