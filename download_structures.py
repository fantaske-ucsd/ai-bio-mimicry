from bio_functions import *

pkl_path = ["data/tuberculosis-human-pretmalign.pkl", "data/listeria-human-pretmalign.pkl", "data/salmonella-human-pretmalign.pkl"]
output_dir = "data/pdb_files/" # pdb file output

match = {
    'uniprot_human_id' : '',
    'uniprot_bac_id' : '',
    'pdb_human_ids' : [],
    'pdb_bac_ids' : [],
    "human_pdb_exists": True,
    "bac_pdb_exists": True,
}


for m in matches :
    print(f"Checking h:{m['uniprot_human_id']} vs b:{m['uniprot_bac_id']}")
    for ph in m['pdb_human_ids'] :
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

        for pb in m['pdb_bac_ids']:
            print(f"Fetching {ph.upper()} and {pb.upper()}")
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
