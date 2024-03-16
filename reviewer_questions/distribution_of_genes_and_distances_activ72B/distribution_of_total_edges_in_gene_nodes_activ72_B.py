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
edge_value = torch.tensor(extract[1], device=device)
edge_ends = torch.tensor(extract[0], device=device)

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

node_id_with_gene = torch.tensor(df_node_coord['hic_node_id'].tolist(), device=device)

def count_ids_in_tensor(edge_ends, ids_of_interest):
    counts_per_id = []
    for node_id in ids_of_interest:
        mask = (edge_ends[1] == node_id) | (edge_ends[0] == node_id)
        edges = torch.stack([edge_ends[0][mask], edge_ends[1][mask]], dim=1)
        edges_sorted, _ = torch.sort(edges, dim=1)
        edges_unique, _ = torch.unique(edges_sorted, dim=0, return_inverse=True)
        counts_per_id.append(edges_unique.shape[0])
    return counts_per_id

counts_per_id = count_ids_in_tensor(edge_ends, node_id_with_gene)

# Merge and save node ids with gene promoter and n interactions
n_edges = pd.DataFrame({
    'hic_node_id': node_id_with_gene.tolist(),
    'n_edges' : counts_per_id,
    })
df_node_coord = df_node_coord.merge(n_edges, how='inner', on=['hic_node_id'])
df_node_coord['connected'] = df_node_coord['connected'].astype(int)
reorder = ['chr', 'start', 'end', 'hic_node_id', 'gene_id', 'expression_class', 'connected', 'n_edges']
df_node_coord = df_node_coord[reorder]
df_node_coord.to_csv('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/reviewer_questions/distribution_of_genes_and_distances_activ72B/gene_node_ids_with_n_edges.bed',sep='\t', index=False)

font_name = 'Nimbus Roman'
plt.figure(figsize=(10, 6))
ax = sns.histplot(data=df_node_coord, x='n_edges', bins=250, kde=True)
ax.set_title('GhmCN - Distribution of Edges in Nodes w/Promoter', fontsize=24, fontname=font_name)
ax.set_xlabel('Edges', fontsize=22, fontname=font_name)
ax.set_ylabel('Frequency', fontsize=22, fontname=font_name)
plt.xticks(fontsize=22, fontname=font_name)
plt.yticks(fontsize=22, fontname=font_name)
plt.xlim(5,50)
plt.savefig('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/reviewer_questions/distribution_of_genes_and_distances_activ72B/distribution_of_edges_in_promoter_bearing_nodes_torch_GhmCN.pdf')
plt.close()
