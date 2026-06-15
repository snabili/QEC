import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

import sys, argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
# importing custom codes
from qec import utils, config

plotpath = config.PLOT_DIR
datapath = config.DATA_DIR
logpath = config.LOG_DIR
logger = utils.setup_logging(log_path=logpath + "/IQ_plots.txt")



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
plt.savefig(plotpath + '/iq_plot.pdf')


# To make 1D I&Q distribution along with gaussian fit
parser = argparse.ArgumentParser(description="QEC IBM_FEZ")
# Add arguments
parser.add_argument("dist", type=str,  help="1D Distribution",   default='I')
parser.add_argument("qu",   type=int,  help="Qubit-ID",   default=128)

args = parser.parse_args()
logger.info(f'Distribution: {args.dist}, Qubit-ID: {args.qu}')

grid_config = [
    (0, 0, args.qu, X00, '|00>'),
    (0, 1, args.qu, X01, '|01>'),
    (1, 0, args.qu, X10, '|10>'),
    (1, 1, args.qu, X11, '|11>')
]
iq_dist = {'I': [0,2], 'Q': [1,3]}
iq_index = {127: iq_dist[args.dist][0],
            128: iq_dist[args.dist][1]}

fig, axes = plt.subplots(2, 2, figsize=(11, 9))
for row, col, i, X, lbl in grid_config:
    ax = axes[row, col]
    ind = iq_index[args.qu]
    IQ = X[:,ind]
    mu, std = norm.fit(IQ)
    x = np.linspace(min(IQ), max(IQ), 60)
    stats_text = f"State: {lbl}\n$\mu$: {mu:.2e}\n$\sigma$: {std:.2e}"
    ax.text(
        0.05, 0.95, stats_text, transform=ax.transAxes, verticalalignment='top', horizontalalignment='left', fontsize=10, 
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.7, edgecolor='gray')
    )

    ax.hist(IQ, bins=60, density=True, alpha=0.6)
    ax.plot(x, norm.pdf(x, mu, std), 'k', linewidth=2)
    ax.set_xlabel(args.dist)
    ax.set_ylabel("Density")

fig.suptitle(f"{args.dist} Distribution\nQubit-ID = {args.qu}", fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(plotpath + '/'+args.dist+'_dist_qu'+str(args.qu) +'.pdf')