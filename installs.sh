#!/bin/bash

# run from repo root dir

mkdir data
# download data automatically if possible

python3 -m venv venv

pip install umap-learn
pip install pandas
pip install h5py
pip install openpyxl
pip install xgboost
pip install torch
pip install matplotlib
pip install seaborn
pip install plotly
pip install faiss
pip install --upgrade scikit-learn xgboost
