#!/bin/bash

# run from repo root dir

python3 -m venv venv || exit

source ./venv/bin/activate

mkdir -p data/pdb_files
# download data automatically if possible

pip install umap-learn
pip install pandas
pip install h5py
pip install openpyxl
pip install xgboost
pip install torch
pip install matplotlib
pip install seaborn
pip install plotly
pip install faiss-cpu
pip install faiss-gpu
pip install scikit-learn

