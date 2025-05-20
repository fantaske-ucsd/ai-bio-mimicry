import pandas as pd
import requests
import os
from tqdm import tqdm

csv_path = "data/tuberculosis-human-results.csv"     # Path to your CSV file
tsv_path = "data/pdb_chain_uniprot.tsv"
output_dir = "data/pdb_files"      # Folder to save downloaded PDBs
pdb_columns = [0, 1]          # Index of the first two columns with PDB IDs

uniprot_df = pd.read_csv(csv_path)
pdb_df = pd.read_csv(tsv_path, sep='\t', low_memory=False).set_index("SP_PRIMARY")

for i, _ in tqdm(enumerate(uniprot_df.index)) :
    human_id = uniprot_df.at[i, 'human_id']
    if human_id not in pdb_df.index :
        continue


    bac_id = uniprot_df.at[i, 'bacteria_id']
    if bac_id not in pdb_df.index :
        continue

    t = pdb_df.at[human_id, "PDB"]
    p = pdb_df.at[bac_id, "PDB"]

    print(f"h: {human_id} b: {bac_id} cycle {i}")

    print(f"\nHuman len {len(t)}")
    print(type(t))

    print(f"\nPathogen len {len(p)}")
    print(type(p))

    break

#t = pdb_df.at[bac_id, "PDB"]

"""
# === EXTRACT UNIQUE PDB IDs ===
pdb_ids = set()
for col in pdb_columns:
    pdb_ids.update(df.iloc[:, col].dropna().astype(str).str.strip().str.lower())

# === DOWNLOAD EACH PDB FILE ===
for pdb_id in pdb_ids:
    if len(pdb_id) != 4:
        print(f"Skipping invalid PDB ID: {pdb_id}")
        continue
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    out_path = os.path.join(output_dir, f"{pdb_id}.pdb")

    try:
        response = requests.get(url)
        response.raise_for_status()
        with open(out_path, "w") as f:
            f.write(response.text)
        print(f"Downloaded: {pdb_id}")
    except Exception as e:
        print(f"Failed to download {pdb_id}: {e}")
"""
