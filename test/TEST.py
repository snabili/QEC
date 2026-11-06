import os, sys, argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
# importing custom codes
from qec import utils, config
import pickle as pkl

datapath = getattr(config, 'DATA_DIR', 'data')
graph_path = os.path.join(datapath, "graph.pkl")
with open(graph_path, "rb") as f: G = pkl.load(f)

dist =3 

def test_get_valid_centers():
    centers = utils.get_valid_centers(G, dist)
    assert isinstance(centers, list)

def test_find_d3_patches():
    patches = utils.find_d3_patches(G,dist)
    assert isinstance(patches, list)