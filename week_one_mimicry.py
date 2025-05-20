from bio_functions import *
import sys

print(sys.argv)

data_dir = "data"

if len(sys.argv) > 1 and not os.path.isdir(sys.argv[1]):
    print(f"Directory not found! Using default")
else :
    data_dir = sys.argv[1]

data = dir_h5_to_np(data_dir)
df = dir_xlsx_to_df(data_dir)

combined_data = new_embed_data()

for d in tqdm(data):
    combined_data['data'].extend(d['data'])
    combined_data['keys'].extend(d['keys'])
    combined_data['label'].extend(d['label'])

combined_df = pd.DataFrame()
for d in df:
    combined_df = pd.concat([combined_df, d], ignore_index=False)

print(list(combined_df.columns))

missing_entries = []
if len(combined_df) != len(combined_data['data']) :
    print(f"MISMATCH: df length {len(combined_df)} embed length {len(combined_data['data'])}")
    if len(combined_df) > len(combined_data['data']) :
        missing_entries = list(combined_df.index)
        entries = list(combined_data['keys'])
    else :
        missing_entries = list(combined_data['keys'])
        entries = list(combined_df.index)

    for k in tqdm(entries):
        if k in missing_entries :
            missing_entries.remove(k)

    if len(missing_entries) > 0 :
        print(missing_entries)

combined_data['clusters'] = cluster_embeddings(combined_data['data'], 0.7)
#print(f"{len(combined_data['data'])} {len(combined_data['clusters'])} {len(combined_data['label'])}")

plot_embed_data(combined_data, "umap_data")
