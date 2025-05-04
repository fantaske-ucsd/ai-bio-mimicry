# Function library that may or maynot be useful

import h5py
from tqdm import tqdm
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
import os
import joblib
from sklearn.model_selection import train_test_split
import numpy as np
import random
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import precision_score, recall_score, f1_score
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import cross_val_score, KFold
from sklearn.model_selection import LeaveOneOut
import pickle
import seaborn as sns
import torch
import torch.nn.functional as F
import warnings

# cleans up output when not debugging
#warnings.filterwarnings("ignore")

embed_data = {
    'label': "",
    'keys': [],
    'data': [],
    'clusters': [],
}

def new_embed_data():
    return embed_data.deepcopy()

def cluster_embeddings(embeddings, distance_threshold):
    """
    Cluster embeddings using a CD-HIT like approach with cosine distance.
    Assumes we want all pairs of embeddings in different clusters to have
    distance > threshold.


    Args:
        embeddings: numpy array of shape (n_samples, n_features)
        distance_threshold: float, maximum cosine distance for clusters
    """
    n_samples = len(embeddings)

    # Normalize embeddings for cosine distance - do this once
    print("\nNormalizing embeddings...")
    norms = np.linalg.norm(embeddings, axis=1)
    normalized_embeddings = embeddings / norms[:, np.newaxis]

    # Initialize clusters
    cluster_labels = np.full(n_samples, -1)  # -1 means unassigned
    current_cluster = 0
    cluster_representatives = []  # Store indices of cluster representatives

    print("\nStarting clustering...")

    # First sequence becomes first cluster representative
    cluster_representatives.append(0)
    cluster_labels[0] = current_cluster
    current_cluster += 1

    # For each remaining sequence
    for i in tqdm(range(0, n_samples), desc="Clustering"):
        # Calculate distances to all cluster representatives
        # Since vectors are normalized, dot product gives cosine similarity directly
        similarities = np.dot(normalized_embeddings[i],
                            normalized_embeddings[cluster_representatives].T)
        distances = 1 - similarities

        # If all distances are greater than threshold, start new cluster
        if np.all(distances > distance_threshold):
            cluster_representatives.append(i)
            cluster_labels[i] = current_cluster
            current_cluster += 1
        else:
            # Assign to closest cluster within threshold
            closest_cluster = cluster_labels[cluster_representatives[np.argmin(distances)]]
            cluster_labels[i] = closest_cluster

    print(f"\nFinished clustering. Found {current_cluster} clusters")

    return cluster_labels


def sample_clusters_for_split(cluster_labels, rando, target_train_fraction=0.8):
    """
    Randomly sample clusters until we get approximately target_fraction of proteins


    Args:
        cluster_labels: array where each element is a cluster ID for that protein
        target_train_fraction: what fraction of total proteins we want in training
    """
    total_proteins = len(cluster_labels)
    target_train_size = int(total_proteins * (1 - target_train_fraction))

    # Get unique cluster IDs and shuffle them
    cluster_ids = np.unique(cluster_labels)
    np.random.shuffle(cluster_ids)

    # Keep track of how many proteins we've accumulated
    current_size = 0
    test_clusters = []

    # Keep adding clusters until we reach target size
    for c in cluster_ids:
        if rando :
            cluster_id = random.choice(cluster_ids)
            cluster_ids.remove(cluster_id)
        else :
            cluster_id = c
        cluster_size = np.sum(cluster_labels == cluster_id)
        if current_size + cluster_size <= target_train_size:
            test_clusters.append(cluster_id)
            current_size += cluster_size

    # Get the protein indices for each set
    test_mask = np.isin(cluster_labels, test_clusters)
    train_sum = total_proteins - np.sum(test_mask)

    print(f"Split results:")
    print(f"Training set: {train_sum} proteins ({train_sum/total_proteins:.2%})")
    print(f"Test set: {np.sum(test_mask)} proteins ({np.sum(test_mask)/total_proteins:.2%})")

    return test_mask, (len(test_clusters), len(cluster_ids) - len(test_clusters))

def get_positive_pairs(df, embed, N, col):
    res = []
    for _ in range(0, N) :
        # choose a random protein and make sure it interacts
        while True :
            val1 = random.choice(df.index.tolist())
            if df.at[val1, col] == 1 :
                one = embed['keys'].index(val1)
                break
        # choose a second unique random protein and make sure it interacts
        while True :
            val2 = random.choice(df.index.tolist())
            if df.at[val2, col] == 1 and val2 != val1 :
                two = embed['keys'].index(val2)
                res.append((one, two, 1))
                break

    return res

def get_clustered_hard_negatives(df, embed, N, col):
    res = []
    for _ in range(0, N) :
        # pick a binding protein
        while True :
            val1 = random.choice(df.index.tolist())
            if df.at[val1, col] == 1 :
                break

        # get it's cluster and the indicies of the rest of the cluster
        cl = embed['clusters'][embed['keys'].index(val1)]
        ind = np.where(embed['clusters'] == cl)[0]
        one = embed['keys'].index(val1)

        # get a non binder
        for i in ind :
            two = random.choice(ind)
            if df.at[embed['keys'][two], col] == 0 :
                res.append((one, two, -1))
                break

    return res

def get_random_negative_pairs(df, embed, N, col):
    res = []
    for _ in range(0, N) :
        # choose a random protein and make sure it interacts
        while True :
            val1 = random.choice(df.index.tolist())
            if df.at[val1, col] == 1 :
                one = embed['keys'].index(val1)
                break
        # choose a second unique random protein and make sure it doesn't bind
        while True :
            val2 = random.choice(df.index.tolist())
            if df.at[val2, col] == 0 :
                two = embed['keys'].index(val2)
                res.append((one, two, -1))
                break

    return res

def h5_to_np(h5path):
    with h5py.File(h5path, "r") as h5f:
        keys = list(h5f.keys())
        datasets = [np.array(h5f[key]) for key in keys]
        return  (keys, np.vstack(datasets))

def h5_to_dict(h5file):
    def recursively_load_h5_group(h5_group):
        result = {}
        for key, item in h5_group.items():
            if isinstance(item, h5py.Group):
                result[key] = recursively_load_h5_group(item)  # Recursive call for groups
            elif isinstance(item, h5py.Dataset):
                result[key] = item[()]  # Load datasets as NumPy arrays
        return result

    with h5py.File(h5file, 'r') as file:
        return recursively_load_h5_group(file)

def get_cos_sim(embed1, embed2):
    t1 = torch.tensor(embed1, dtype=torch.float32).unsqueeze(0)
    t2 = torch.tensor(embed2, dtype=torch.float32).unsqueeze(0)
    return F.cosine_similarity(t1, t2)

def get_files_in_dir(fdir):
    return [file for file in os.listdir(fdir) if os.path.isfile(os.path.join(fdir, file))]

def dir_h5_to_np(fpath):
    """
    returns a list of embed_data dictionaries
    """

    r=[]
    files = get_files_in_dir(fpath)
    for f in files:
        if f.endswith('.h5') :
            new = new_embed_data()
            ret = h5_to_np(fpath + '/' + f)
            new['keys'] = list(ret[0])
            new['data'] = list(ret[1])
            new['label'] = f.split(".")[0]
            r.append(new)

    return r

def xlsx_to_df(xlsx_file):
    return None
