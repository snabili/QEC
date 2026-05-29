import numpy as np
import pandas as pd
import os

import xgboost as xgb

from sklearn.covariance import MinCovDet, EmpiricalCovariance
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.svm import SVC

import sys, argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
# importing custom codes
from qec import utils, config

plotpath = config.PLOT_DIR
logpath = config.LOG_DIR
datapath = config.DATA_DIR

logger = utils.setup_logging(log_path=logpath + "/ml_gridsearch.txt")

parser = argparse.ArgumentParser(description="QEC IBM_FEZ")
# Add arguments
parser.add_argument("mlalgo", type=str,  help="ml: svc, bdt",   default='bdt')
args = parser.parse_args()
logger.info(f'GridSearch ML: {args.mlalgo}')


data = np.load(datapath + '/x.npz')
X00, X11, X01, X10 = data['X00'], data['X11'], data['X01'], data['X10']
X = np.vstack([X00, X01, X10, X11])
y = np.concatenate([np.zeros(len(X00)),np.ones(len(X01)),np.full(len(X10), 2),np.full(len(X11), 3)])
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=42)


if args.mlalgo == 'svc':
    param_grid = {'svc__C': [0.1, 1, 10, 100],'svc__gamma': ['scale', 'auto', 0.01, 0.001]}
    model = SVC(kernel='rbf')
elif args.mlalgo == 'bdt':
    param_grid = {'max_depth': [3, 5], 'learning_rate': [0.01, 0.1],'n_estimators': [100, 200]}
    model = xgb.XGBClassifier(objective='multi:softprob', eval_metric='mlogloss')

grid_search = GridSearchCV(estimator=model, param_grid=param_grid, cv=5, scoring='accuracy',n_jobs=-1) # 5 crossvalidation
grid_search.fit(X_train, y_train)

best_model = grid_search.best_estimator_
print(f"Best parameters found: {grid_search.best_params_}")

results_df = pd.DataFrame(grid_search.cv_results_)
#logger.info(results_df.sort_values('rank_test_score'))
if args.mlalgo == 'svc': 
    logger.info(results_df[['param_svc__C', 'param_svc__gamma', 'mean_test_score', 'rank_test_score']].sort_values('rank_test_score'))
elif args.mlalgo == 'bdt': 
    logger.info(results_df[['param_max_depth', 'param_learning_rate', 'param_n_estimators', 'mean_test_score', 'rank_test_score']].sort_values('rank_test_score'))



# save the results as it takes ~20 minutes to run
results_df.to_csv(datapath + '/' + args.mlalgo+ '_gridsearch.csv', index=True)

'''# To read it later
q = pd.read_csv(datapath + '/' + args.mlalgo+ '_gridsearch.csv')
logger.info(results_df.sort_values('rank_test_score'))'''
