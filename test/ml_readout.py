import numpy as np
import matplotlib.pyplot as plt
import joblib
import os

import xgboost as xgb
from sklearn.covariance import MinCovDet
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC # --> SVC: classification class of SVM (State Vector Machine)
from sklearn.mixture import GaussianMixture
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split

import sys, argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
# importing custom codes
from qec import utils, config

plotpath = config.PLOT_DIR
logpath = config.LOG_DIR
datapath = config.DATA_DIR

logger = utils.setup_logging(log_path=logpath + "/ml_readout.txt")

# Data prep: download level=1 I&Q data for 15,000 shots for double-qubits four states: |00>, |11>, |01>, |10>; train-test split 
data = np.load(os.getcwd()+'/files/datafiles/x.npz')
X00, X11, X01, X10 = data['X00'], data['X11'], data['X01'], data['X10']

X = np.vstack([X00, X01, X10, X11])
y = np.concatenate([np.zeros(len(X00)),np.ones(len(X01)),np.full(len(X10), 2),np.full(len(X11), 3)])

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=42)

# MinimumCovarianceDeterminant (MCD) on Mahalanobis distance to reduce outliers due to coherence decay time
def train_mcd_discriminator_multi(X_list):
    params = []
    for X in X_list:
        mcd = MinCovDet().fit(X)
        mu = mcd.location_
        inv_cov = np.linalg.inv(mcd.covariance_)
        params.append((mu, inv_cov))    
    return params

def predict_mcd_multi(X_test, all_params):
    distances = []
    for mu, inv_cov in all_params:
        diff = X_test - mu
        dist = np.sqrt(np.sum(np.dot(diff, inv_cov) * diff, axis=1))
        distances.append(dist)
    distances = np.array(distances) 
    return np.argmin(distances, axis=0) 

XT = [X_train[y_train==i] for i in range(4)]
mcd_params = train_mcd_discriminator_multi(XT)
predictions = predict_mcd_multi(X_test, mcd_params)

accuracy = np.mean(predictions == y_test)
logger.info(f"MCD Multi-class Accuracy: {accuracy * 100:.2f}%")

state_map = {0: "|00>", 1: "|01>", 2: "|10>", 3: "|11>"}
mah_confmatrix = {}
for i in range(4):
    X_test_i = X_test[y_test == i]
    pred_01 = predict_mcd_multi(X_test_i, mcd_params)
    unique, counts = np.unique(pred_01, return_counts=True)
    results = dict(zip(unique, counts))
    logger.info(f"Prediction Breakdown for True State: {state_map[i]}")
    per = []
    for label, count in results.items():
        p = (count / len(pred_01))
        per.append(np.round(p,3).item())
        logger.info(f"  Predicted as {state_map[label]}: {count} times ({p*100:.2f}%)")
    mah_confmatrix[i] = per
    #logger.info('\n','--'*100)

state_labels = ["|00>", "|01>", "|10>", "|11>"]

# LR Model
clf = LogisticRegression()
clf.fit(X_train, y_train)
y_pred_lr = clf.predict(X_test)
cm_lr = confusion_matrix(y_test, y_pred_lr)

# GMM used MCD initial means and precision to exclude outliers
initial_means = np.array([p[0] for p in mcd_params])
initial_precisions = np.array([p[1] for p in mcd_params])
gmm = GaussianMixture(n_components=4, covariance_type='full', means_init=initial_means, random_state=42) # means_init forces the mapping!
gmm.fit(X_train)
y_pred_gmm = gmm.predict(X_test)
cm_gmm = confusion_matrix(y_test, y_pred_gmm)

# SVC Model
svc = SVC(kernel='rbf', gamma='scale', C=0.1) # used grid search method to tune svc parameters
svc.fit(X_train, y_train)
y_pred_svc = svc.predict(X_test)
cm_svc = confusion_matrix(y_test, y_pred_svc)

# BDT Model
bdt = xgb.XGBClassifier(n_estimators=150, max_depth=4, learning_rate=0.05, objective='binary:logistic',eval_metric='logloss',random_state=42,n_jobs=-1)
bdt.fit(X_train, y_train)
y_pred_bdt = bdt.predict(X_test)
cm_bdt = confusion_matrix(y_test, y_pred_bdt)


# Save the trained model to a file
joblib.dump(svc, os.getcwd()+'/files/datafiles/svc.joblib')
joblib.dump(bdt, os.getcwd()+'/files/datafiles/bdt.joblib')
joblib.dump(gmm, os.getcwd()+'/files/datafiles/gmm.joblib')
joblib.dump(clf, os.getcwd()+'/files/datafiles/lr.joblib')

np.savez(datapath + '/mah.npz', mah_confmatrix)


logger.info("Logistic Regression CM:\n", cm_lr)
logger.info("\n SVM CM:\n", cm_svc)
logger.info("\n GMM CM:\n", cm_gmm)
logger.info("\n Mah CM:\n", np.array([mah_confmatrix[i] for i in mah_confmatrix.keys()]))
logger.info("\n BDT CM:\n", cm_bdt)
logger.info("\n", "*"*150, "\n")

def fidelity(cm):
    P_i_given_i = []
    for i in range(cm.shape[0]):
        Nii = cm[i,i]
        Nij = sum([cm[i,j] for j in range(cm.shape[1]) if j!= i])
        P_i_given_i.append(Nii/(Nii+ Nij))
    return P_i_given_i
    
F_lr = fidelity(cm_lr)
F_svc = fidelity(cm_svc)
F_gmm = fidelity(cm_gmm)
F_mah = [mah_confmatrix[i][i] for i in mah_confmatrix.keys()]
F_bdt = fidelity(cm_bdt)


logger.info(f"{'State':<10} | {'LR  Fidelity':<17} | {'SVC Fidelity':<13} | {'GMM Fidelity':<13} |  {'MAH Fidelity':<15} | {'BDT Fidelity':<15}")
logger.info("-" * 100)

for i in range(4):
    state_label = f"|{bin(i)[2:].zfill(2)}>"
    logger.info(f"{state_label:<10} | "
          f"{F_lr[i]:>10.3f}        |  "
          f"{F_svc[i]:>10.3f}   |"
          f"{F_gmm[i]:>10.3f}     |"
          f"{F_mah[i]:>10.3f}     |"
          f"{F_bdt[i]:>10.3f}"
          )
logger.info(f'\n mean:\n   LR={np.mean(F_lr):.3f},\n   SVC={np.mean(F_svc):.3f},\n   GMM={np.mean(F_gmm):.3f},\n   MAH={np.mean(F_mah):.3f},\n   BDT={np.mean(F_bdt):.3}')

'''
Now plotting I&Q with decision boundaries

'''
parser = argparse.ArgumentParser(description="QEC IBM_FEZ")
# Add arguments
parser.add_argument("mlalgo", type=str,  help="ml: svc, bdt, gmm, lr, mah",   default='bdt')
args = parser.parse_args()
logger.info(f'Plot ML boundaries: {args.mlalgo}')

if args.mlalgo != 'mah': mdl = joblib.load(datapath + '/' + args.mlalgo + '.joblib')
else: mdl = np.load(datapath + '/' + args.mlalgo + '.npz')

y_pred = mdl.predict(X_test)
cm = confusion_matrix(y_test, y_pred)


def plot_iq_decision_boundaries(X_test, y_test, model, qubit_index=1):
    """
    Plots the IQ scatter data and RBF SVC decision boundaries.
    qubit_index: 1 for (I1, Q1), 2 for (I2, Q2)
    """
    # Select indices based on which qubit we are plotting
    if qubit_index == 1:
        idx_i, idx_q = 0, 1
        other_i, other_q = 2, 3
        title = "Qubit 1 Decision Boundaries (holding Q2 at mean)"
    else:
        idx_i, idx_q = 2, 3
        other_i, other_q = 0, 1
        title = "Qubit 2 Decision Boundaries (holding Q1 at mean)"

    state_labels = ["|00>", "|01>", "|10>", "|11>"]

    # Create a mesh grid to plot boundaries
    i_min, i_max = X_test[:, idx_i].min() - 0.1, X_test[:, idx_i].max() + 0.1
    q_min, q_max = X_test[:, idx_q].min() - 0.1, X_test[:, idx_q].max() + 0.1
    ii, qq = np.meshgrid(np.linspace(i_min, i_max, 100),
                         np.linspace(q_min, q_max, 100))

    # To predict on the grid, we need 4D input. 
    # We hold the "other" qubit at its average value.
    other_i_val = np.mean(X_test[:, other_i])
    other_q_val = np.mean(X_test[:, other_q])
    
    grid_points = np.c_[ii.ravel(), qq.ravel()] # --> what is this doing??
    # Construct 4D grid based on which qubit is being visualized
    if qubit_index == 1:
        full_grid = np.c_[grid_points, 
                          np.full(len(grid_points), other_i_val), 
                          np.full(len(grid_points), other_q_val)]
    else:
        full_grid = np.c_[np.full(len(grid_points), other_i_val), 
                          np.full(len(grid_points), other_q_val),
                          grid_points]

    # Predict across the grid
    Z = model.predict(full_grid)
    Z = Z.reshape(ii.shape)

    # Plotting
    plt.figure(figsize=(8, 6))

    # Draw decision regions
    plt.contourf(ii, qq, Z, alpha=0.2, cmap='viridis')
    
    # Scatter plot of actual test points
    scatter = plt.scatter(X_test[:, idx_i], X_test[:, idx_q], c=y_test, 
                          edgecolors='k', s=20, cmap='viridis')
    
    #CREATE LEGEND:
    # legend_elements returns handles and labels for each unique value in c=y_test
    handles, _ = scatter.legend_elements()
    plt.legend(handles, state_labels, title="Quantum States", loc="upper right")
    
    plt.xlabel(f'I{qubit_index} (V)')
    plt.ylabel(f'Q{qubit_index} (V)')
    plt.title(title)
    plt.colorbar(scatter, ticks=[0, 1, 2, 3], label='States: |00>, |01>, |10>, |11>')
    plt.grid(True, alpha=0.3)
    plt.ylim(min(X_test[:,1])*1.2, max(X_test[:,1])*2.0,)
    plt.savefig(plotpath + '/iqmap_ml-' + args.mlalgo + '_boundaries.png')

# Usage:
plot_iq_decision_boundaries(X_test, y_test, mdl, qubit_index=1)