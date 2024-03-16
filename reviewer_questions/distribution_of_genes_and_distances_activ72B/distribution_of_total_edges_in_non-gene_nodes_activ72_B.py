# srun --nodes=1 --ntasks=1 --cpus-per-task=8 --gpus=1 --mem=100g --time=06:00:00 --pty bash -i
# Sconda conda activate ghmc1; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/reviewer_questions/distribution_of_genes_and_distances_activ72B; python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import torch
import torch_geometric
from   scipy.sparse import load_npz
import os

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

data_dir='/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/src/data/B72_CMS_10000bp/'
df_genes = pd.read_pickle(os.path.join(data_dir,'df_genes_reg0.pkl'))
extract = torch_geometric.utils.from_scipy_sparse_matrix(
    load_npz(
        os.path.join(
            data_dir,
            'hic_sparse.npz'
            )
        )
    )
edge_value = extract[1]
edge_ends = extract[0]

df_node_coord = pd.read_pickle(os.path.join(data_dir,'df_node_coord.pkl')).dropna().drop(columns = ['gene_catalog_name','Abbrev', 'expression_lvl']).astype(int)
df_node_coord = (
    df_node_coord
    .merge(df_genes, how='inner', on=['hic_node_id'])
    .drop(columns = ['gene_catalog_name'])
    .rename(columns = {
        'abbrev':'gene_id',
        'expression_lvl':'expression_class',
        'chr_num':'chr',
        'chr_start':'start',
        'chr_finish':'end',
    }))
widow_size = 10_000

node_id_with_gene = torch.tensor(df_node_coord['hic_node_id'].tolist())

mask_0 = torch.isin(edge_ends[0], node_id_with_gene)
mask_1 = torch.isin(edge_ends[1], node_id_with_gene)

combined_mask = mask_0 | mask_1

# Apply the mask to subset extract
edge_ends = extract[0][:, ~combined_mask]
edge_value = extract[1][~combined_mask]

# unique after removing
# filtered_indices[0].unique().shape
# filtered_indices[1].unique().shape

def count_ids_in_tensor(edge_ends, ids_of_interest):
    counts_per_id = []
    for node_id in ids_of_interest:
        mask = (edge_ends[1] == node_id) | (edge_ends[0] == node_id)
        edges = torch.stack([edge_ends[0][mask], edge_ends[1][mask]], dim=1)
        edges_sorted, _ = torch.sort(edges, dim=1)
        edges_unique, _ = torch.unique(edges_sorted, dim=0, return_inverse=True)
        counts_per_id.append(edges_unique.shape[0])
    return counts_per_id


counts_per_id = count_ids_in_tensor(edge_ends, edge_ends[0].unique())

# Merge and save node ids with gene promoter and n interactions
n_edges = pd.DataFrame({
    'n_edges' : counts_per_id,
    })

font_name = 'Nimbus Roman'
plt.figure(figsize=(10, 6))
ax = sns.histplot(data=n_edges, x='n_edges', bins=250, kde=True)
ax.set_title('GhmCN - Distribution of Edges between Nodes w/ No Promoter', fontsize=24, fontname=font_name)
ax.set_xlabel('Edges', fontsize=22, fontname=font_name)
ax.set_ylabel('Frequency', fontsize=22, fontname=font_name)
plt.xticks(fontsize=22, fontname=font_name)
plt.yticks(fontsize=22, fontname=font_name)
plt.xlim(0,50)
plt.savefig('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/reviewer_questions/distribution_of_genes_and_distances_activ72B/distribution_of_edges_between_non-gene-interacting_nodes_torch_GhmCN.pdf')
plt.close()
