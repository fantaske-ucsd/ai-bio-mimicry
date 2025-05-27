import pandas as pd
import requests
import os
from tqdm import tqdm
from bio_functions import *

xlsx_path = ["data/tuberculosis-human-results.xlsx", "data/listeria-human-results.xlsx", "data/salmonella-human-results.xlsx"]
tsv_path = "data/pdb_chain_uniprot.tsv"

match = {
    'uniprot_human_id' : '',
    'uniprot_bac_id' : '',
    'pdb_human_ids' : [],
    'pdb_bac_ids' : [],
    "human_pdb_exists": True,
    "bac_pdb_exists": True,
}

matches = []
results = pd.DataFrame()

pdb_df = pd.read_csv(tsv_path, sep='\t', low_memory=False).set_index("SP_PRIMARY")

for x in xlsx_path :
    uniprot_df = pd.read_excel(x)
    human_ids = uniprot_df['human_id'].unique()
    human_groups = uniprot_df.groupby('human_id')

    patho_name = x.split("/")[1].split("-")[0]
    print(f"Processing hits for {patho_name}")

    for human_id in tqdm(human_ids) :
        human_pdb_exists = True
        if human_id not in pdb_df.index :
            human_pdb_exists = False

        group_df = human_groups.get_group(human_id)
        bac_ids = group_df['bacteria_id'].unique()
        for bac_id in bac_ids :
            patho_pdb_exists = True
            if bac_id not in pdb_df.index :
                patho_pdb_exists = False

            zscores = group_df[ (group_df['bacteria_id'] == bac_id) & (group_df['human_id'] == human_id) ]
            scores = zscores['seq_identity_aligned']
            for sc in scores :
                if sc < 0.25 :
#                    print(f"Found very disimilar protein, score {sc}")
#                    print(f"h: {human_id} b: {bac_id}")
#                    print(zscores)

                    if results.empty :
                        results = zscores
                    else :
                        results = pd.concat([results, zscores])

                    if human_pdb_exists :
                        t = pdb_df.at[human_id, "PDB"]
                        if not isinstance(t, str):
                            tu = t.dropna().unique()
                        else:
                            tu = []
                            tu.append(t)
                    else :
                        tu = [] # TODO fill the list with what we need for alphafold

                    if patho_pdb_exists :
                        p = pdb_df.at[bac_id, "PDB"]
                        if not isinstance(p, str):
                            pu = p.dropna().unique()
                        else:
                            pu = []
                            pu.append(p)
                    else :
                        pu = [] # TODO fill the list with what we need for alphafold

                    new_match = dict(match)
                    new_match['uniprot_human_id'] = human_id
                    new_match['uniprot_bac_id'] = bac_id
                    new_match['pdb_human_ids'] = list(tu)
                    new_match['pdb_bac_ids'] = list(pu)
                    new_match['human_pdb_exists'] = human_pdb_exists
                    new_match['bac_pdb_exists'] = patho_pdb_exists
                    matches.append(new_match)

    print(f"Found {len(matches)} matches")

    fn = "data/" + patho_name + "-human-pretmalign.xlsx"
    df_hits = pd.DataFrame(results)
    df_hits.to_excel(fn, index=False)

    fn = fn.split('.')[0] + '.pkl'
    save_dict_to_pkl(matches, fn)
