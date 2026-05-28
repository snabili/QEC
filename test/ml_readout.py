import numpy as np
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


data = np.load(os.getcwd()+'/files/datafiles/x.npz')
X00, X11, X01, X10 = data['X00'], data['X11'], data['X01'], data['X10']

X = np.vstack([X00, X01, X10, X11])

y = np.concatenate([
    np.zeros(len(X00)),
    np.ones(len(X01)),
    np.full(len(X10), 2),
    np.full(len(X11), 3)
])

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=42)

# MCD
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
debug = False
print(f"MCD Multi-class Accuracy: {accuracy * 100:.2f}%")

state_map = {0: "|00>", 1: "|01>", 2: "|10>", 3: "|11>"}
mah_confmatrix = {}
for i in range(4):
    X_test_i = X_test[y_test == i]
    pred_01 = predict_mcd_multi(X_test_i, mcd_params)
    unique, counts = np.unique(pred_01, return_counts=True)
    results = dict(zip(unique, counts))
    if debug == True: print(f"Prediction Breakdown for True State: {state_map[i]}")
    per = []
    for label, count in results.items():
        p = (count / len(pred_01))
        per.append(np.round(p,3).item())
        if debug == True: print(f"  Predicted as {state_map[label]}: {count} times ({p*100:.2f}%)")
    mah_confmatrix[i] = per
    if debug == True: print('\n','--'*100)

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


if debug == True:
    print("Logistic Regression CM:\n", cm_lr)
    print("\n SVM CM:\n", cm_svc)
    print("\n GMM CM:\n", cm_gmm)
    print("\n Mah CM:\n", np.array([mah_confmatrix[i] for i in mah_confmatrix.keys()]))
    #print("\n BDT CM:\n", cm_bdt)
    print("\n", "*"*150, "\n")

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


print(f"{'State':<10} | {'LR  Fidelity':<17} | {'SVC Fidelity':<13} | {'GMM Fidelity':<13} |  {'MAH Fidelity':<15} | {'BDT Fidelity':<15}")
print("-" * 100)

for i in range(4):
    state_label = f"|{bin(i)[2:].zfill(2)}>"
    print(f"{state_label:<10} | "
          f"{F_lr[i]:>10.3f}        |  "
          f"{F_svc[i]:>10.3f}   |"
          f"{F_gmm[i]:>10.3f}     |"
          f"{F_mah[i]:>10.3f}     |"
          f"{F_bdt[i]:>10.3f}"
          )
print(f'\n mean:\n   LR={np.mean(F_lr):.3f},\n   SVC={np.mean(F_svc):.3f},\n   GMM={np.mean(F_gmm):.3f},\n   MAH={np.mean(F_mah):.3f},\n   BDT={np.mean(F_bdt):.3}')