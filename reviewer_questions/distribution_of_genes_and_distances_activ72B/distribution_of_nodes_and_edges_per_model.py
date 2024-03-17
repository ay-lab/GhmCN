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

def analyze_graph(extract):
    # Unpack the edge index tensor and the edge weight tensor
    edge_indices, edge_weights = extract
    
    # Flatten the edge index tensor and find the unique nodes
    unique_nodes = torch.unique(edge_indices.flatten())
    total_nodes = unique_nodes.size(0)
    
    # The total number of edges is the size of the edge weights tensor
    total_edges = edge_weights.size(0)
    
    # Average number of edges per node is total edges divided by total nodes
    average_edges_per_node = total_edges / total_nodes
    
    print(f"Total number of nodes: {total_nodes}")
    print(f"Total number of edges: {total_edges}")
    print(f"Average number of edges (interactions) per node: {round(average_edges_per_node,3)}\n")



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
analyze_graph(extract)
# Whole GhmCN graph B72
# Total number of nodes: 238616
# Total number of edges: 3284006
# Average number of edges (interactions) per node: 13.762723371441982

hic_nodes_mask = pd.read_csv('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/reviewer_questions/top_interactions_fithic/fithic_results/hicpro_to_fithic/hic_node_id.csv')
all_genes_mask = pd.read_pickle(os.path.join(data_dir,'df_node_coord.pkl')).dropna().drop(columns = ['gene_catalog_name','Abbrev', 'expression_lvl']).astype(int)

hic_nodes_mask = torch.tensor(hic_nodes_mask['hic_node_id'].values)
all_genes_mask = torch.tensor(all_genes_mask['hic_node_id'].values)
mask_0 = torch.isin(extract[0][0], hic_nodes_mask)
mask_1 = torch.isin(extract[0][1], hic_nodes_mask)
combined_mask = mask_0 | mask_1
gene_mask_0 = torch.isin(extract[0][0], all_genes_mask)
gene_mask_1 = torch.isin(extract[0][1], all_genes_mask)
combined_gene_mask = gene_mask_0 | gene_mask_1
extract_fithic = (extract[0][:, combined_mask], extract[1][combined_mask])
extract_with_gene = (extract[0][:, combined_gene_mask], extract[1][combined_gene_mask])
extract_without_gene = (extract[0][:, ~combined_gene_mask], extract[1][~combined_gene_mask])

analyze_graph(extract_fithic)
# Total number of nodes: 228969
# Total number of edges: 2219048
# Average number of edges (interactions) per node: 9.691477885652644

analyze_graph(extract_with_gene)
# Total number of nodes: 129935
# Total number of edges: 490280
# Average number of edges (interactions) per node: 3.773

analyze_graph(extract_without_gene)
# Total number of nodes: 114592
# Total number of edges: 1064958
# Average number of edges (interactions) per node: 9.29347598436191

# Nod do it iteratively for the networks used in the main text:

cells=[
  'B00',
  'B72',
  'DP',
  'Naive_CD4T',
  'Naive_CD8T',
  'Th2',
  ]

cell_dir='/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/src/data/'
for cell in cells:
    extract = torch_geometric.utils.from_scipy_sparse_matrix(
    load_npz(
        os.path.join(
            cell_dir,
            f'{cell}_CMS_10000bp/hic_sparse.npz'
            )
        )
    )
    print(f'{cell} Stats:')
    analyze_graph(extract)



exit()
"""
B00 Stats:
Total number of nodes: 238668
Total number of edges: 3570126
Average number of edges (interactions) per node: 14.959

B72 Stats:
Total number of nodes: 238616
Total number of edges: 3284006
Average number of edges (interactions) per node: 13.763

DP Stats:
Total number of nodes: 236528
Total number of edges: 3349702
Average number of edges (interactions) per node: 14.162

Naive_CD4T Stats:
Total number of nodes: 238522
Total number of edges: 3438270
Average number of edges (interactions) per node: 14.415

Naive_CD8T Stats:
Total number of nodes: 236014
Total number of edges: 3235358
Average number of edges (interactions) per node: 13.708

Th2 Stats:
Total number of nodes: 235384
Total number of edges: 1389000
Average number of edges (interactions) per node: 5.901
"""