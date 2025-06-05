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
def get_all_match_pdb(org, m, output_dir) :
    if m[f"{org}_pdb_exists"] == True :
        retry = False
        for pdb_id in m[f'pdb_{org}_ids'] :
            retry |= download_pdb(pdb_id, output_dir, False)

        if retry :
            download_pdb(m[f'uniprot_{org}_id'], output_dir, True)
    else :
        download_pdb(m[f'uniprot_{org}_id'], output_dir, True)

def get_all_pdbs(m) :
    print(f"Checking h:{m['uniprot_human_id']} vs b:{m['uniprot_bac_id']}")
    get_all_match_pdb("human", m, output_dir)
    get_all_match_pdb("bac", m, output_dir)

for p in pkl_path :
    matches = get_dict_from_pkl(p)
    if __name__ == "__main__":
        with mp.Pool(processes=24) as pool:  # adjust based on your machine
            pool.map(get_all_pdbs, matches)
