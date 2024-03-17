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

# 
# 
# sparse matrix with node-node (edges) and value arrays
data_dir='/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/src/data/B72_CMS_10000bp/'
extract = torch_geometric.utils.from_scipy_sparse_matrix(
    load_npz(
        os.path.join(
            data_dir,
            'hic_sparse.npz'
            )
        )
    )

# 
# 
# Get node IDs with genes

# # Data frame with all node's info
df_node_coord = pd.read_pickle(os.path.join(data_dir,'df_node_coord.pkl')).drop(columns = ['gene_catalog_name','Abbrev', 'expression_lvl']).astype(int)

# # Data frame indicating which nodes have promoter(s)
df_genes = pd.read_pickle(os.path.join(data_dir,'df_genes_reg0.pkl'))

# # Subset further with `df_genes` data frame
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
    })
    )

# # tensor of node IDs with genes
node_id_with_gene = torch.tensor(df_node_coord['hic_node_id'].tolist())

# 
# 
# Generate the FitHiC starting graph

fithic_nodes_mask = pd.read_csv('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/reviewer_questions/top_interactions_fithic/fithic_results/hicpro_to_fithic/hic_node_id.csv')
# # Get nodes forming the matrix
node_id_fithic = torch.tensor(fithic_nodes_mask['hic_node_id'].tolist())
# # Generate masks to subset complete GhmCN matrix
fithic_mask_0 = torch.isin(extract[0][0], node_id_fithic)
fithic_mask_1 = torch.isin(extract[0][1], node_id_fithic)
combined_fithic_mask = fithic_mask_0 | fithic_mask_1
# # Subset GhmCN matrix to generate FitHiC 
extract_fithic = (extract[0][:, combined_fithic_mask], extract[1][combined_fithic_mask])


# 
# 
# Function to leverage torch (GPUs) and count edges per node for all nodes in graph
def count_edges_per_node(edge_index):
    # Combine both directions since the graph might be undirected
    all_nodes = torch.cat((edge_index[0], edge_index[1]))
    # Count occurrences of each node to get the degree (total interactions per node)
    unique_nodes, counts_per_node = torch.unique(all_nodes, return_counts=True)
    return unique_nodes, counts_per_node//2


# 
# 
# DRY: Generate custom annotation to be used in plt histogram
def summary_edges(counts):
    annotation_str = f'N={counts.size()[0]}\nE={counts.sum().item()}\navg(E)={round(counts.float().mean().item(),2)}'
    print(annotation_str)
    return annotation_str


# 
# 
# DRY: Wrapper function to plt histogram and annotation of graph
def plot_wrap(counts,title,file_name):
    annotation_str = summary_edges(counts)
    font_name = 'Nimbus Roman'
    plt.figure(figsize=(10, 6))
    ax = sns.histplot(data=counts.numpy(), bins=250, kde=True)
    ax.set_title(title, fontsize=24, fontname=font_name)
    ax.set_xlabel('Edges', fontsize=22, fontname=font_name)
    ax.set_ylabel('Frequency', fontsize=22, fontname=font_name)
    plt.xticks(fontsize=22, fontname=font_name)
    plt.yticks(fontsize=22, fontname=font_name)
    plt.xlim(0,50)
    plt.text(0.9, 0.9, annotation_str,fontsize=18, horizontalalignment='right', verticalalignment='top', transform=plt.gca().transAxes, fontname=font_name)
    plt.savefig(file_name)
    plt.close()


# 
# 
# Count edges for both GhmCN and FitHiC complete matrices
nodes_ghmcn, counts_ghmcn = count_edges_per_node(extract[0])
nodes_fithic, counts_fithic = count_edges_per_node(extract_fithic[0])

# Plot histogram from the whole GhmCN matrix
plot_wrap(counts_ghmcn,
          'GhmCN - Distribution of Edges in All Nodes',
          '/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/reviewer_questions/distribution_of_genes_and_distances_activ72B/distribution_of_edges_in_all_nodes_ghmcn.pdf',
)

# Plot histogram from the whole FitHiC matrix
plot_wrap(counts_fithic,
          'FitHiC - Distribution of Edges in All Nodes',
          '/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/reviewer_questions/distribution_of_genes_and_distances_activ72B/distribution_of_edges_in_all_nodes_fithic.pdf',
)


# Plot subset of GhmCN of nodes with genes:
mask_nodes_with_genes = torch.isin(nodes_ghmcn, node_id_with_gene)
counts_per_node_w_gene = counts_ghmcn[mask_nodes_with_genes]
plot_wrap(counts_per_node_w_gene,
          'GhmCN - Distribution of Edges in Nodes w/Promoter',
          '/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/reviewer_questions/distribution_of_genes_and_distances_activ72B/distribution_of_edges_in_nodes_with_promoter_ghmcn.pdf',
)

# Plot subset of GhmCN of nodes without genes (logical NOT of mask of nodes with genes):
counts_per_node_wo_gene = counts_ghmcn[~mask_nodes_with_genes]
plot_wrap(counts_per_node_wo_gene,
          'GhmCN - Distribution of Edges in Nodes w/o Promoter',
          '/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/reviewer_questions/distribution_of_genes_and_distances_activ72B/distribution_of_edges_in_nodes_without_promoter_ghmcn.pdf',
)

# Plot subset of FitHiC of nodes with genes:
mask_nodes_fithic_with_genes = torch.isin(nodes_fithic, node_id_with_gene)
counts_fithic_per_node_w_gene = counts_fithic[mask_nodes_fithic_with_genes]
plot_wrap(counts_fithic_per_node_w_gene,
          'FitHiC - Distribution of Edges in Nodes w/Promoter',
          '/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/reviewer_questions/distribution_of_genes_and_distances_activ72B/distribution_of_edges_in_nodes_with_promoter_fithic.pdf',
)

# Plot subset of FitHiC of nodes without genes (logical NOT of mask of nodes with genes):
counts_fithic_per_node_wo_gene = counts_fithic[~mask_nodes_fithic_with_genes]
plot_wrap(counts_fithic_per_node_wo_gene,
          'FitHiC - Distribution of Edges in Nodes w/Promoter',
          '/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/reviewer_questions/distribution_of_genes_and_distances_activ72B/distribution_of_edges_in_nodes_without_promoter_fithic.pdf',
)
