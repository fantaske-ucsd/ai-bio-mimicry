from bio_functions import *
import multiprocessing as mp

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

def mk_file(fn) :
    try:
        with open(fn, "w") as file:
            file.write("")
    except FileExistsError :
        print("Skipping {fn}")

def get_alpha(alpha, output_dir) :
    success = download_alphafold_pdb(alpha, output_dir)
    if not success :
        mk_file(output_dir + alpha + ".fail")

def get_all_pdbs(m) :
    print(f"Checking h:{m['uniprot_human_id']} vs b:{m['uniprot_bac_id']}")
    if m["human_pdb_exists"] == True :
        retry = False
        for pdb_id in m['pdb_human_ids'] :
            retry |= download_pdb(pdb_id, output_dir)

        if retry :
#            print(f"All H PDB downloads failed, trying alphafold for {m['uniprot_human_id']}")
            get_alpha(m['uniprot_human_id'], output_dir)
    else :
#        print(f"Downloading H {m['uniprot_human_id']} from alphafold")
        get_alpha(m['uniprot_human_id'], output_dir)

    if m["bac_pdb_exists"] == True :
        retry = False
        for pdb_id in m['pdb_bac_ids']:
            retry |= download_pdb(pdb_id, output_dir)

        if retry :
#            print(f"All B PDB downloads failed, trying alphafold for {m['uniprot_bac_id']}")
            get_alpha(m['uniprot_bac_id'], output_dir)
    else :
#        print(f"Downloading B {m['uniprot_bac_id']} from alphafold")
        get_alpha(m['uniprot_bac_id'], output_dir)

for p in pkl_path :
    matches = get_dict_from_pkl(p)
    if __name__ == "__main__":
        with mp.Pool(processes=16) as pool:  # adjust based on your machine
            pool.map(get_all_pdbs, matches)
