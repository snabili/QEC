import matplotlib.pyplot as plt
import networkx as nx
import pickle as pkl
import sys, os, argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from qec import config, utils

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister

'''
This code collects all plottings in a decorator format
'''

# define path to read/write files
plotpath = config.PLOT_DIR
datapath = config.DATA_DIR

with open(datapath + "/graph.pkl", "rb") as f:
    G = pkl.load(f)


scripter = utils.Scripter()

@scripter
def draw_logical_qubit(G, result):
    # Add parser argument to build 17 qubits patch based on qubit id 
    parser = argparse.ArgumentParser(description="QEC Patch with distance == 3")
    parser.add_argument("qid", type=int,  help="qubit_id", default=23)
    args = parser.parse_args()

    result = utils.find_logical_qubit_full(G, args.qid)

    patch_nodes = result['patch']
    data_nodes = result['data']
    x_path = result['logical_x']
    z_path = result['logical_z']
    syndromes = [n for n in patch_nodes if n not in data_nodes]
    
    sub = G.subgraph(patch_nodes)
    pos = nx.spring_layout(sub, seed=42)  # Maintains consistent shape    
    plt.figure(figsize=(5, 4))   
    
    # Drawing patches with networks
    nx.draw_networkx_edges(sub, pos, alpha=0.3, edge_color='gray') # draw all edges   
    nx.draw_networkx_nodes(sub, pos, nodelist=syndromes, node_color='lightgrey', node_size=300, label='Syndromes') # draw syndrome qubits
    nx.draw_networkx_nodes(sub, pos, nodelist=data_nodes, node_color='skyblue', node_size=600, label='Data Qubits') # draw data qubits
    nx.draw_networkx_nodes(sub, pos, nodelist=x_path, node_color='red', node_size=700, alpha=0.5) # draw logical-x path
    nx.draw_networkx_nodes(sub, pos, nodelist=z_path, node_color='green', node_size=700, alpha=0.5) # draw logical-z path
    nx.draw_networkx_labels(sub, pos, font_size=10, font_family="sans-serif") # add qubits ID
    
    plt.title(f"Heavy-Hex d=3 Logical Qubit\nRed = Logical X | Green = Logical Z")
    plt.legend(scatterpoints=1)
    plt.axis('off')
    plt.savefig(plotpath + '/logicpatch_syndromes.png')

@scripter
def plot_qc():
    parser = argparse.ArgumentParser(description="Quantum Circuit plot")
    parser.add_argument("stb", type=str,  help="stabilizer", default=23)
    args = parser.parse_args()
	# Define named registers
    ancilla = QuantumRegister(1, name='ancilla |0>')
    data_qubit = QuantumRegister(1, name=f'data qubit |\psi>')
    c_reg = ClassicalRegister(1, name='syndrome')
	
	# Create the circuit using the registers
    qc = QuantumCircuit(ancilla, data_qubit, c_reg)
	
    if args.stb == "X":
        # Add gates using register indexing
        qc.h(ancilla[0])
        qc.h(data_qubit[0])
        qc.cx(ancilla[0], data_qubit[0])
        qc.h(ancilla[0])
        qc.h(data_qubit[0])
    else:
        qc.cx(data_qubit[0], ancilla[0])
        qc.measure(ancilla[0], c_reg[0])
	
	# Save as a transparent PNG for GitHub
    filename = '/x_stabilizer.png' if args.stb == "X" else '/z_stabilizer.png'
    qc.draw('mpl', filename=plotpath + filename, style={'backgroundcolor': '#ffffff00'})

if __name__ == '__main__':
    scripter.run()
