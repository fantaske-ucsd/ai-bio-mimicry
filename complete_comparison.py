import h5py
import numpy as np
import faiss
import pandas as pd
from tqdm import tqdm
import os
import torch
print("Torch sees GPU:", torch.cuda.is_available())
print("Device count:", torch.cuda.device_count())
print("CUDA_VISIBLE_DEVICES:", os.environ.get("CUDA_VISIBLE_DEVICES"))

# File paths
bacteria_h5 = "data/salmonella-unreviewed-5k.h5"
human_h5 = "data/human-unreviewed-83k.h5"

# Load embeddings from .h5 files}
with h5py.File(bacteria_h5, 'r') as f:
    bacteria_ids = list(f.keys())
    bacteria_embeddings = np.stack([f[seq_id][:] for seq_id in bacteria_ids])

# Load human embeddings
with h5py.File(human_h5, 'r') as f:
    human_ids = list(f.keys())
    human_embeddings = np.stack([f[seq_id][:] for seq_id in human_ids])

# Convert embeddings to float32 before normalization
bacteria_embeddings = bacteria_embeddings.astype(np.float32)
human_embeddings = human_embeddings.astype(np.float32)

# Normalize for cosine similarity (inner product in FAISS)
faiss.normalize_L2(bacteria_embeddings)
faiss.normalize_L2(human_embeddings)

# Create FAISS GPU index
# Check how many GPUs are visible (respects CUDA_VISIBLE_DEVICES)
try:
    gpu_count = faiss.get_num_gpus()
    if gpu_count == 0:
        raise RuntimeError("No GPUs available, falling back to CPU mode.")
    
    print(f"Using {gpu_count} GPU(s) for FAISS")
    gpu_resources = [faiss.StandardGpuResources() for _ in range(gpu_count)]
    gpu_options = faiss.GpuMultipleClonerOptions()
    gpu_options.shard = True  # or False, depending on your goal (sharding vs replication)

    cpu_index = faiss.IndexFlatIP(bacteria_embeddings.shape[1])
    cpu_index.add(bacteria_embeddings.astype(np.float32))
    gpu_index = faiss.index_cpu_to_gpu_multiple_py(gpu_resources, cpu_index, gpu_options)
    index = gpu_index

except Exception as e:
    print(f"[WARNING] GPU initialization failed: {e}")
    print("Using CPU FAISS instead.")

    index = faiss.IndexFlatIP(bacteria_embeddings.shape[1])
    index.add(bacteria_embeddings.astype(np.float32))

# Search: top 100 hits for each human protein
top_k = 500
D, I = index.search(human_embeddings.astype(np.float32), top_k)

# Create results DataFrame
results = []
for i, human_id in enumerate(tqdm(human_ids, desc="Formatting results")):
    for j in range(top_k):
        if D[i, j] >= 0.5 :
            results.append({
                "human_id": human_id,
                "bacteria_id": bacteria_ids[I[i, j]],
                "cosine_similarity": D[i, j]
            })

df_hits = pd.DataFrame(results)
df_hits.to_csv("top_100_hits_per_human_protein.csv", index=False)
