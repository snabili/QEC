import numpy as np
import os
import joblib

import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split

import sys, argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
# importing custom codes
from qec import utils, config

plotpath = config.PLOT_DIR
datapath = config.DATA_DIR
logpath = config.LOG_DIR

parser = argparse.ArgumentParser(description="QEC IBM_FEZ")
# Add arguments
parser.add_argument("mlalgo", type=str,  help="ml: svc, bdt, gmm, lr, mah",   default='bdt')
args = parser.parse_args()

logger = utils.setup_logging(log_path=logpath + "/ml_readout.txt")
logger.info(f'ML algo: {args.mlalgo}')

# Data prep: download level=1 I&Q data for 15,000 shots for double-qubits four states: |00>, |11>, |01>, |10>; train-test split 
data = np.load(datapath +'/x.npz')
X00, X11, X01, X10 = data['X00'], data['X11'], data['X01'], data['X10']

X = np.vstack([X00, X01, X10, X11])
y = np.concatenate([np.zeros(len(X00)),np.ones(len(X01)),np.full(len(X10), 2),np.full(len(X11), 3)])

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=42)


if args.mlalgo != 'mah': mdl = joblib.load(datapath + '/' + args.mlalgo + '.joblib')
else: mdl = np.load(datapath + '/' + args.mlalgo + '.npz')

y_pred = mdl.predict(X_test)
cm = confusion_matrix(y_test, y_pred)

'''
# to get rid of scientific notation
np.set_printoptions(suppress=True, precision=5)
'''

def run_2qubit_ibu(measured_counts, matrix_A, max_iterations=15):
    total_shots = np.sum(measured_counts)
    P_measured = measured_counts / total_shots
    num_states = matrix_A.shape[0]
    P_true_guess = np.ones(num_states) / num_states    
    for iteration in range(max_iterations):
        P_true_next = np.zeros(num_states)        
        for m in range(num_states):
            if P_measured[m] == 0:
                continue                
            P_m_total = np.sum(matrix_A[m, :] * P_true_guess)
            for t in range(num_states):
                P_t_given_m = (matrix_A[m, t] * P_true_guess[t]) / P_m_total
                P_true_next[t] += P_measured[m] * P_t_given_m
        P_true_guess = P_true_next        
    return P_true_guess * total_shots

# scikit-learn matrix_A: Row = True, Columns = Measured ===> matrix [i,j]: P(Predicted=j | Actual=i)
cm_norm_rows = confusion_matrix(y_test, y_pred, normalize='true') 

# IBU matrix_A: Row = Measured, Column = True ===> matrix [i,j]: P(Predicted=i | Actual=j)
matrix_A_crosstalk = cm_norm_rows.T
assert np.allclose(np.sum(matrix_A_crosstalk, axis=0), 1.), "Columns must sum to 1.0!"

state_labels = ['|00>', '|01>', '|10>', '|11>']

logger.info("IBU effect on CM Diagonal elements Fidelity")
logger.info("State  | Raw Counts | xtalk-Mit Counts | Truth Counts")
logger.info("-" * 50)
for i in range(4):
    raw_hardware_counts = cm[i, :] 
    mitigated_counts = run_2qubit_ibu(raw_hardware_counts, matrix_A_crosstalk, max_iterations=15)
    logger.info(f"|{state_labels[i]}>   |    {raw_hardware_counts[i]:4d}   |     {round(mitigated_counts[i]):4d}        |   {round(sum(raw_hardware_counts))}")


matrix_A_mitigated = []
for i in range(4):
    raw_hardware_counts = cm[i,:]
    matrix_A_mitigated.append(run_2qubit_ibu(raw_hardware_counts, matrix_A_crosstalk, max_iterations=15))

matrix_A_final = np.array(matrix_A_mitigated).round().astype(int)

fig, axes = plt.subplots(1, 2, figsize=(9, 4))

# Raw Confusion Matrix
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=state_labels, yticklabels=state_labels, ax=axes[0],cbar_kws={'label': 'Counts'})
axes[0].set_title('Raw Hardware Confusion Matrix\n(BDT Classifier Predictions)', fontsize=11, pad=10)
axes[0].set_xlabel('Predicted State', fontsize=8)
axes[0].set_ylabel('Prepared State (True)', fontsize=8)

# Mitigated Matrix
sns.heatmap(matrix_A_final, annot=True, fmt='d', cmap='Greens', xticklabels=state_labels, yticklabels=state_labels, ax=axes[1],cbar_kws={'label': 'Counts'})
axes[1].set_title('Crosstalk-Mitigated Matrix\n(After 15 Iterations of IBU)', fontsize=11, pad=10)
axes[1].set_xlabel('Mitigated Target State', fontsize=8)
axes[1].set_ylabel('Prepared State (True)', fontsize=8)

plt.tight_layout()
plt.savefig(plotpath + '/readout_cm.png')