import pandas as pd
import requests
import os
from tqdm import tqdm

csv_path = "data/tuberculosis-human-results.csv"     # Path to your CSV file
tsv_path = "data/pdb_chain_uniprot.tsv"
output_dir = "data/pdb_files"      # Folder to save downloaded PDBs

match = {
    'uniprot_human_id' : '',
    'uniprot_bac_id' : '',
    'pdb_human_ids' : [],
    'pdb_bac_ids' : [],
}

matches = []

uniprot_df = pd.read_csv(csv_path)
pdb_df = pd.read_csv(tsv_path, sep='\t', low_memory=False).set_index("SP_PRIMARY")

human_ids = uniprot_df['human_id'].unique()
human_groups = uniprot_df.groupby('human_id')

for human_id in tqdm(enumerate(human_ids)) :
    if human_id[1] not in pdb_df.index :
        continue

    bac_ids = human_groups.get_group(human_id[1])['bacteria_id'].unique()
    for bac_id in enumerate(bac_ids) :
        if bac_id[1] not in pdb_df.index :
            continue

        t = pdb_df.at[human_id[1], "PDB"]
        p = pdb_df.at[bac_id[1], "PDB"]

#        print(f"h: {human_id} b: {bac_id}")

        if not isinstance(t, str):
            tu = t.dropna().unique()
        else:
            tu = []
            tu.append(t)
#        print(f"Unique human len {len(tu)}")
#        print(tu)

        if not isinstance(p, str):
            pu = p.dropna().unique()
        else:
            pu = []
            pu.append(p)
#        print(f"Unique pathogen len {len(pu)}")
#        print(pu)

        new_match = dict(match)
        new_match['uniprot_human_id'] = human_id[1]
        new_match['uniprot_bac_id'] = bac_id[1]
        new_match['pdb_human_ids'] = list(tu)
        new_match['pdb_bac_ids'] = list(pu)
        matches.append(new_match)

print(f"Found {len(matches)} matches")

#print(t.iloc[0])
#print(p.iloc[0])

"""
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
