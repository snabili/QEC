import os
import numpy as np
import matplotlib.pyplot as plt

import sys, argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
# importing custom codes
from qec import utils, config

plotpath = config.PLOT_DIR
datapath = config.DATA_DIR
logpath = config.LOG_DIR


data = np.load(datapath + '/x.npz')
X00, X11, X01, X10 = data['X00'], data['X11'], data['X01'], data['X10']

def snr_double_qubits(X_prim, X_sec, nq):
    if nq == 1:
        cluster0 = X_prim[:, :2]
        cluster1 = X_sec[:, :2]
    else:
        cluster0 = X_prim[:, 2:]
        cluster1 = X_sec[:, 2:]

    mu0 = np.mean(cluster0, axis=0)
    mu1 = np.mean(cluster1, axis=0)

    center_distance = np.linalg.norm(mu1 - mu0)
    direction_vector = (mu1 - mu0) / center_distance
    cov0 = np.cov(cluster0, rowvar=False)
    cov1 = np.cov(cluster1, rowvar=False)

    var0_projected = direction_vector.T @ cov0 @ direction_vector
    var1_projected = direction_vector.T @ cov1 @ direction_vector

    sigma0 = np.sqrt(var0_projected)
    sigma1 = np.sqrt(var1_projected)

    average_noise = 0.5 * (sigma0 + sigma1)
    snr = center_distance / average_noise
    return snr


grid_config = [
    (0, 0, 127, 0, 1, X00, '|00>', X10, '|10>', 0.2), 
    (0, 1, 128, 2, 3, X00, '|00>', X01, '|01>', 0.2),
    (1, 0, 127, 0, 1, X01, '|01>', X11, '|11>', 0.2),
    (1, 1, 128, 2, 3, X10, '|10>', X11, '|11>', 0.2)
]

fig, axes = plt.subplots(2, 2, figsize=(11, 9))

for row, col, i, i_col, q_col, X_prim, lbl_prim, X_sec, lbl_sec, a_sec in grid_config:
    if i == 127: snr = snr_double_qubits(X_prim, X_sec, 1)
    else: snr = snr_double_qubits(X_prim, X_sec, 2)
    ax = axes[row, col]
    ax.scatter(X_prim[:, i_col], X_prim[:, q_col], s=20, alpha=1.0, label=lbl_prim)
    ax.scatter(X_sec[:, i_col], X_sec[:, q_col], s=20, alpha=a_sec, label=lbl_sec)
    ax.set_xlabel(fr"$I_{{{i}}}$ (Volts)", fontsize=14)
    ax.set_ylabel(fr"$Q_{{{i}}}$ (Volts)", fontsize=14)
    ax.set_title(fr'$Qubit_{{{i}}}$: {lbl_prim} vs {lbl_sec}, SNR = {snr:.2f}', fontsize=13, fontweight='bold')
    ax.legend(loc='upper right', fontsize=11)
    ax.grid(True, linestyle=':', alpha=0.6)

plt.tight_layout()
plt.savefig(plotpath + '/iq_plot.png')
