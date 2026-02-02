import matplotlib.pyplot as plt
import networkx as nx


def find_d3_patches(graph):
    patches = []
    # In heavy-hex, logical patches are centered on 'Data' qubits.
    # We look for nodes with degree 3 as potential centers.
    centers = [n for n, d in graph.degree() if d == 3]
    print(centers)
    
    for c in centers:
        # We need to go 4 steps out to encompass the 17-qubit diamond
        nodes_at_dist = nx.single_source_shortest_path_length(graph, c, cutoff=4)
        subnodes = list(nodes_at_dist.keys())
        
        # A d=3 heavy-hex patch requires exactly 17 qubits 
        # (9 data + 8 syndrome). We filter for clusters that 
        # satisfy the minimum size.
        if len(subnodes) >= 17:
            # We sort by distance to ensure we are getting the 
            # 17 qubits closest to our center 'c'
            sorted_nodes = sorted(subnodes, key=lambda x: nodes_at_dist[x])
            patch = sorted_nodes[:17]
            patches.append(patch)
            
    return patches

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
    plt.savefig('plots/logicpatch.pdf')


def get_true_roles(full_graph, patch_nodes):
    # In a d=3 heavy-hex code, there are exactly 9 Data qubits form a 3x3 grid.

    # (In a subgraph, data qubits connect to syndromes, never other data)
    sub = full_graph.subgraph(patch_nodes)
    
    # Heavy-hex lattice is Bipartite; color it into two sets.
    coloring = nx.bipartite.color(sub)
    
    set_0 = [n for n, color in coloring.items() if color == 0]
    set_1 = [n for n, color in coloring.items() if color == 1]
    
    # The set with 9 nodes: Data Qubits
    # The set with 8 nodes: Syndrome Qubits
    if len(set_0) == 9:
        data_q, syndrome_q = set_0, set_1
    else:
        data_q, syndrome_q = set_1, set_0
        
    # Distinguish X and Z syndromes using the FULL graph degree
    # X-syndromes (stars) connect to 3 data qubits in the full lattice
    # Z-syndromes (bridges) connect to 2 data qubits in the full lattice
    x_syn = [n for n in syndrome_q if full_graph.degree(n) == 3]
    z_syn = [n for n in syndrome_q if full_graph.degree(n) == 2]
    
    return {"data": sorted(data_q), "x_syndrome": sorted(x_syn), "z_syndrome": sorted(z_syn)}



def find_logical_qubit_full(G, center_node):
    # 1. Identify the 17-qubit cluster
    nodes_at_dist = nx.single_source_shortest_path_length(G, center_node, cutoff=4)
    patch_nodes = sorted(nodes_at_dist.keys(), key=lambda x: nodes_at_dist[x])[:17]
    sub = G.subgraph(patch_nodes)

    # 2. Assign Roles using Bipartite Coloring
    # This ensures 9 Data qubits and 8 Syndrome qubits
    coloring = nx.bipartite.color(sub)
    set_0 = [n for n, color in coloring.items() if color == 0]
    set_1 = [n for n, color in coloring.items() if color == 1]
    data_nodes = set_0 if len(set_0) == 9 else set_1
    
    # 3. Find all valid D-S-D-S-D strings (Length 4 edges)
    logical_strings = []
    for start in data_nodes:
        for end in data_nodes:
            if start >= end: continue
            for path in nx.all_simple_paths(sub, start, end, cutoff=4):
                if len(path) == 5 and all(path[i] in data_nodes for i in [0, 2, 4]):
                    logical_strings.append(path)

    # 4. Find the X and Z pair (must intersect at exactly ONE data qubit)
    for i, path_a in enumerate(logical_strings):
        for path_b in logical_strings[i+1:]:
            data_a = {path_a[0], path_a[2], path_a[4]}
            data_b = {path_b[0], path_b[2], path_b[4]}
            intersection = data_a.intersection(data_b)
            
            if len(intersection) == 1:
                # We've found our conjugate pair!
                return {
                    "patch": patch_nodes,
                    "data": sorted(list(data_nodes)),
                    "logical_x": [path_a[0], path_a[2], path_a[4]],
                    "logical_z": [path_b[0], path_b[2], path_b[4]],
                    "pivot": list(intersection)[0]
                }
    return None

def draw_logical_qubit(G, result):
    patch_nodes = result['patch']
    data_nodes = result['data']
    x_path = result['logical_x']
    z_path = result['logical_z']
    
    sub = G.subgraph(patch_nodes)
    pos = nx.spring_layout(sub, seed=42)  # Maintains consistent shape
    
    plt.figure(figsize=(5, 4))
    
    # 1. Draw all edges
    nx.draw_networkx_edges(sub, pos, alpha=0.3, edge_color='gray')
    
    # 2. Draw Syndrome Qubits (Small, Grey)
    syndromes = [n for n in patch_nodes if n not in data_nodes]
    nx.draw_networkx_nodes(sub, pos, nodelist=syndromes, node_color='lightgrey', 
                           node_size=300, label='Syndromes')
    
    # 3. Draw Data Qubits (Blue)
    nx.draw_networkx_nodes(sub, pos, nodelist=data_nodes, node_color='skyblue', 
                           node_size=600, label='Data Qubits')
    
    # 4. Highlight Logical X Path (Red Edges/Halo)
    nx.draw_networkx_nodes(sub, pos, nodelist=x_path, node_color='red', 
                           node_size=700, alpha=0.5)
    
    # 5. Highlight Logical Z Path (Green Edges/Halo)
    nx.draw_networkx_nodes(sub, pos, nodelist=z_path, node_color='green', 
                           node_size=700, alpha=0.5)

    # 6. Label the nodes with their IDs
    nx.draw_networkx_labels(sub, pos, font_size=10, font_family="sans-serif")
    
    plt.title(f"Heavy-Hex d=3 Logical Qubit\nRed = Logical X | Green = Logical Z")
    plt.legend(scatterpoints=1)
    plt.axis('off')
    plt.savefig('plots/logicpatch_syndromes.pdf')