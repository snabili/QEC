import matplotlib.pyplot as plt
import networkx as nx
import pickle as pkl
import sys, os, argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from qec import config, utils


# Add parser argument to build 17 qubits patch based on qubit id 
parser = argparse.ArgumentParser(description="QEC on IBM_FEZ backend")
parser.add_argument("qid", type=int,  help="qubit_id", default=23)
args = parser.parse_args()

# define path to read/write files
plotpath = config.PLOT_DIR
datapath = config.DATA_DIR

with open(datapath + "/graph.pkl", "rb") as f:
    G = pkl.load(f)

result = utils.find_logical_qubit_full(G, args.qid)

def draw_logical_qubit(G, result):
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

if __name__ == '__main__':
    draw_logical_qubit(G, result)
