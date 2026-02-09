import os

# Base dirs
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(BASE_DIR, ".."))

# Key project dirs
FILE_DIR       = os.path.join(ROOT_DIR, "files")
LOG_DIR        = os.path.join(FILE_DIR, "logs")
PLOT_DIR       = os.path.join(FILE_DIR, "plots")
DATA_DIR      = os.path.join(FILE_DIR, "datafiles")

# Create dirs if not exist
for path in [FILE_DIR, LOG_DIR, PLOT_DIR, DATA_DIR]:
    os.makedirs(path, exist_ok=True)


STRUCTURES_JSON = os.path.join(FILE_DIR, 'structures_dict.json')
LOG_FILE        = os.path.join(LOG_DIR, 'qec.txt')

# Processing settings
MAX_JOBS = 4
BATCH_COUNT = 50  # Number of batches to divide the dataset into
CUTOFF_DISTANCE = 5  # For MinimumDistanceNN

# Debugging
DEBUG = False
