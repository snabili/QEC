import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from sklearn.mixture import GaussianMixture
from sklearn.model_selection import train_test_split
from sklearn.covariance import MinCovDet, EmpiricalCovariance

import os, sys, argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
# importing custom codes
from qec import utils, config

logpath = config.LOG_DIR
datapath = config.DATA_DIR
plotpath = config.PLOT_DIR

logger = utils.setup_logging(log_path=logpath + "/state_2qu_maha.txt")



parser = argparse.ArgumentParser(description="QEC IBM_FEZ")
# Add arguments
parser.add_argument("nst", type=int,  help="States: 0:|00>, 1:|01>, 2:|10>, 3:|11>",   default=0)
parser.add_argument("nqu", type=int,  help="Qubits: 1:127, 2:128",   default=1)
args = parser.parse_args()

logger.info(f'argument state: {args.nst}, qubits: {args.nqu}')

data = np.load(datapath + '/x.npz')

X00 = data['X00']
X01 = data['X01']
X10 = data['X10']
X11 = data['X11']

X = np.vstack([X00, X01, X10, X11])

y = np.concatenate([
    np.zeros(len(X00)),
    np.ones(len(X01)),
    np.full(len(X10), 2),
    np.full(len(X11), 3)
])

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=42)

states = {0:'00', 1:'01', 2:'10', 3:'11'}
qubit_title = {0:f"$Qubit_{{{127}}}$", 1:f"$Qubit_{{{128}}}$"}


def train_mcd_discriminator(X0, X1):
    mcd0 = MinCovDet().fit(X0)
    mu0, cov0 = mcd0.location_, mcd0.covariance_
    
    mcd1 = MinCovDet().fit(X1)
    mu1, cov1 = mcd1.location_, mcd1.covariance_
    
    return (mu0, cov0), (mu1, cov1)

def predict_mcd(X_test, params0, params1):
    mu0, cov0 = params0
    mu1, cov1 = params1
    
    inv_cov0 = np.linalg.inv(cov0)
    diff0 = X_test - mu0
    dist0 = np.sqrt(np.sum(np.dot(diff0, inv_cov0) * diff0, axis=1))
    
    inv_cov1 = np.linalg.inv(cov1)
    diff1 = X_test - mu1
    dist1 = np.sqrt(np.sum(np.dot(diff1, inv_cov1) * diff1, axis=1))
    
    return (dist1 < dist0).astype(int)


def pars_qubit(nqu, ini_state):
    error_q = []  
    tar_state = []  
    test_data_com = {}  
    params1_dict = {} 
    for i in (0,1,2,3):
        if i == ini_state: continue
        if nqu == 1:
            X0_train = X_train[y_train == ini_state, :2]
            X1_train = X_train[y_train == i,         :2]
            X1_test  = X_test [y_test  == i,         :2]
            test_data_ref = X_test[y_test == ini_state, :2]
        else:
            X0_train = X_train[y_train == ini_state, 2:]
            X1_train = X_train[y_train == i,         2:]
            X1_test  = X_test [y_test  == i,         2:]
            test_data_ref = X_test[y_test == ini_state, 2:]

        params0, params1 = train_mcd_discriminator(X0_train, X1_train)
        y_pred = predict_mcd(X1_test, params0, params1)

        error_q.append(np.mean(y_pred).item())
        tar_state.append(i)
        test_data_com[i] = X1_test
        params1_dict[i] = params1
    return error_q, tar_state, params0, params1_dict, test_data_ref, test_data_com



def draw_mcd_ellipse(mu, cov, ax, color, label):
    """Draws a 2-sigma (~95% of 2D dist) MCD ellipse"""
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    
    # 2-sigma width and height
    width, height = 2 * 2 * np.sqrt(vals)
    
    ell = Ellipse(xy=mu, width=width, height=height, angle=theta,
                  edgecolor=color, facecolor=color, alpha=0.2, lw=2, label=label)
    ax.add_patch(ell)
    # Add a border line for clarity
    ell_border = Ellipse(xy=mu, width=width, height=height, angle=theta,
                         edgecolor=color, facecolor='none', lw=2)
    ax.add_patch(ell_border)

nqu, ini_state = args.nqu, args.nst
error_q, tar_state, params0, params1_dict, test_data_ref, test_data_com = pars_qubit(nqu, ini_state)
fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True, sharey=True)

for ax, n in zip(axes, tar_state):
    mu0, cov0 = params0
    mui, covi = params1_dict[n]

    title_str = f"|{states[ini_state]}> vs |{states[n]}>"
    err_idx = tar_state.index(n)

    ax.scatter(test_data_ref[:, 0], test_data_ref[:, 1], s=5, alpha=0.3,
               color='gray', label='Test ' + states[ini_state])

    ax.scatter(test_data_com[n][:, 0], test_data_com[n][:, 1], s=5, alpha=0.3,
               color=f"C{n+1}", label='Test ' + states[n])

    draw_mcd_ellipse(mu0, cov0, ax, 'black', f'MCD {states[ini_state]}')
    draw_mcd_ellipse(mui, covi, ax, f"C{n+1}", f'MCD {states[n]}')

    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    xx, yy = np.meshgrid(np.linspace(x_min, x_max, 100), np.linspace(y_min, y_max, 100))
    grid_points = np.c_[xx.ravel(), yy.ravel()]

    inv_cov0 = np.linalg.inv(cov0)
    inv_covi = np.linalg.inv(covi)

    d0 = np.sqrt(np.sum((grid_points - mu0) @ inv_cov0 * (grid_points - mu0), axis=1))
    di = np.sqrt(np.sum((grid_points - mui) @ inv_covi * (grid_points - mui), axis=1))

    boundary = (d0 - di).reshape(xx.shape)
    ax.contour(xx, yy, boundary, levels=[0], colors='red', linestyles='--')

    ax.set_title(f"{qubit_title[nqu-1]}  {title_str}  Fidelity: {error_q[err_idx]:.2f}")

    ax.set_xlabel("I (Volts)")
    ax.set_ylabel("Q (Volts)")
    ax.legend(loc='upper right', fontsize='small')

plt.tight_layout()

qname = {1: '127', 2: '128'}

plt.savefig(plotpath+ '/mahala_'+states[ini_state]+'_qubit'+qname[nqu]+'.pdf')

def draw_ellipse(mcd, ax, color, label, linestyle='-'):
    """Draws a 3-sigma Mahalanobis contour."""
    cov = mcd.covariance_
    center = mcd.location_
        
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
        
    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))        
    width, height = 2 * 3 * np.sqrt(vals)
        
    ell = Ellipse(xy=center, width=width, height=height, angle=theta,
                      edgecolor=color, facecolor='none', lw=2, 
                      linestyle=linestyle, label=label)
    ax.add_patch(ell)

fig, ax = plt.subplots(1, 2, figsize=(12, 5))
for i,nqu in enumerate([127, 128]):
    if nqu == 127:
        ax[i].scatter(X10[:,0], X10[:,1], color='navy', alpha=0.2, s=15, label=f'IQ State: |10>')
        standard_cov = EmpiricalCovariance().fit(X10[:,:2])
        robust_cov = MinCovDet().fit(X10[:,:2]) 
    else:
        ax[i].scatter(X01[:,2], X01[:,3], color='navy', alpha=0.2, s=15, label=f'IQ State: |01>')
        standard_cov = EmpiricalCovariance().fit(X01[:,2:])
        robust_cov = MinCovDet().fit(X01[:,2:])
    draw_ellipse(standard_cov, ax[i], 'black', 'Standard Mahalanobis (Skewed)', '--')
    draw_ellipse(robust_cov, ax[i], 'green', 'Robust MCD Mahalanobis (Correct)')
    ax[i].set_xlabel(fr"$I_{{{nqu}}}$ (Volts)", fontsize=16)
    ax[i].set_ylabel(fr"$Q_{{{nqu}}}$ (Volts)", fontsize=16)
    ax[i].legend()
    ax[i].grid(True, alpha=0.3)

plt.savefig(plotpath + '/maha_qu1-10_qu2-01.png')
