from utils import build_heavy_hex_stabilizer, find_logical_qubit_full, get_true_role

def test_build_heavy_hex_stabilizer():
    datapath = files/datafiles
    with open(datapath + "/graph.pkl", "rb") as f: G = pkl.load(f)
    result = utils.find_logical_qubit_full(G, args.qid)
    roles = utils.get_true_roles(G, result['patch'])

