import h5py
import numpy as np
import faiss
import pandas as pd
from tqdm import tqdm
import os
import torch
from bio_functions import *

print("Torch sees GPU:", torch.cuda.is_available())
print("Device count:", torch.cuda.device_count())
print("CUDA_VISIBLE_DEVICES:", os.environ.get("CUDA_VISIBLE_DEVICES"))

shard = False

# File paths
bacteria_h5_list = ["data/salmonella-unreviewed-5k.h5", "data/listeria-unreviewed-3k.h5", "data/tuberculosis-unreviewed-4k.h5"]
human_h5 = "data/human-unreviewed-83k.h5"

# Load human embeddings
with h5py.File(human_h5, 'r') as f:
    human_ids = list(f.keys())
    human_embeddings = np.stack([f[seq_id][:] for seq_id in human_ids])

human_df = xlsx_to_df("data/human.xlsx")

human_embeddings = human_embeddings.astype(np.float32)
faiss.normalize_L2(human_embeddings)

gpu_count = faiss.get_num_gpus()
if gpu_count > 0 :
    gpu_resources = [faiss.StandardGpuResources() for _ in range(gpu_count)]
    gpu_options = faiss.GpuMultipleClonerOptions()
    gpu_options.shard = shard  # depending on your goal (sharding vs replication)

for bacteria_h5 in bacteria_h5_list :

    # Load embeddings from .h5 files}
    with h5py.File(bacteria_h5, 'r') as f:
        bacteria_ids = list(f.keys())
        bacteria_embeddings = np.stack([f[seq_id][:] for seq_id in bacteria_ids])

    # Convert embeddings to float32 before normalization
    bacteria_embeddings = bacteria_embeddings.astype(np.float32)

    # Normalize for cosine similarity (inner product in FAISS)
    faiss.normalize_L2(bacteria_embeddings)

    # Create FAISS GPU index
    # Check how many GPUs are visible (respects CUDA_VISIBLE_DEVICES)
    try:
        if gpu_count == 0:
            raise RuntimeError("No GPUs available, falling back to CPU mode.")
        
        print(f"Using {gpu_count} GPU(s) for FAISS")
        cpu_index = faiss.IndexFlatIP(bacteria_embeddings.shape[1])
        cpu_index.add(bacteria_embeddings.astype(np.float32))
        gpu_index = faiss.index_cpu_to_gpu_multiple_py(gpu_resources, cpu_index, gpu_options)
        index = gpu_index

    except Exception as e:
        print(f"[WARNING] GPU initialization failed: {e}")
        print("Using CPU FAISS instead.")

        index = faiss.IndexFlatIP(bacteria_embeddings.shape[1])
        index.add(bacteria_embeddings.astype(np.float32))

    # Search: top_k hits for each human protein
    top_k = 50
    D, I = index.search(human_embeddings.astype(np.float32), top_k)

    patho_name = bacteria_h5.split("/")[1].split("-")[0]
    if "tub" in patho_name :
        patho_df = xlsx_to_df("data/tub.xlsx")
    elif  "salm" in patho_name :
        patho_df = xlsx_to_df("data/salm.xlsx")
    elif "list" in patho_name :
        patho_df = xlsx_to_df("data/list.xlsx")

    results = []
    low = 2.0
    low_a = 10000.0
    high_a = 0.000000001

    for i, human_id in enumerate(tqdm(human_ids, desc="Formatting results")):
        for j in range(top_k):
            if low > D[i, j] :
                low = D[i, j]

            if D[i, j] >= 0.75 :
                patho_seq = patho_df.at[bacteria_ids[I[i, j]], "Sequence"]
                human_seq = human_df.at[human_id, "Sequence"]
                ident, ident2, ta = align_sequences(human_seq, patho_seq)

                if ta.score > high_a :
                    high_a = ta.score

                if ta.score < low_a :
                    low_a = ta.score

                results.append({
                    "human_id": human_id,
                    "bacteria_id": bacteria_ids[I[i, j]],
                    "human_seq": human_seq,
                    "bacteria_seq": patho_seq,
                    "cosine_similarity": D[i, j],
                    "seq_align_score": ta.score,
                    "seq_identity_aligned": ident,
                    "seq_identity_min": ident2,
                })

    print(f"Found {len(results)} hits for {patho_name}")
    print(f"Low value for {patho_name} is {low}")
    print(f"Lowest align: {low_a}, Highest: {high_a}")
    fn = "data/" + patho_name + "-human-results.xlsx"
    df_hits = pd.DataFrame(results)
    df_hits.to_excel(fn, index=False)
