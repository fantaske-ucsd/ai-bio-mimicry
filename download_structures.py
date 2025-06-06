from bio_functions import *
import multiprocessing as mp
import time
import os

pkl_path = ["data/tuberculosis-human-pretmalign.pkl", "data/listeria-human-pretmalign.pkl", "data/salmonella-human-pretmalign.pkl"]
output_dir = "data/pdb_files/" # pdb file output

def get_all_match_pdb(org, m, output_dir) :
    if m[f"{org}_pdb_exists"] == True :
        for pdb_id in m[f'pdb_{org}_ids'] :
            download_pdb(pdb_id, output_dir, False)

        download_pdb(m[f'uniprot_{org}_id'], output_dir, True)
    else :
        download_pdb(m[f'uniprot_{org}_id'], output_dir, True)

def get_all_pdbs(m) :
#    print(f"Checking h:{m['uniprot_human_id']} vs b:{m['uniprot_bac_id']}")
    get_all_match_pdb("human", m, output_dir)
    get_all_match_pdb("bac", m, output_dir)

for p in pkl_path :
    matches = get_dict_from_pkl(p)
    print(f"Starting {p}")
    print(f"Length {len(matches)}")
    if __name__ == "__main__":
        with mp.Pool(processes=24) as pool:  # adjust based on your machine
            result = pool.map_async(get_all_pdbs, matches)
            ready = result.ready()
            while ready == False :
                num_files = len([f for f in os.listdir(output_dir) if os.path.isfile(os.path.join(output_dir, f))])
                print(f"\rFiles Downloaded: {num_files}", end='', flush=True)
                time.sleep(1)
                ready = result.ready()
            print() 
