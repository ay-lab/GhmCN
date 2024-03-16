# Sconda conda activate ghmc1; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/reviewer_questions/distribution_of_genes_and_distances_activ72B; python
import pandas as pd
import numpy as np
import torch
import torch_geometric
from   scipy.sparse import load_npz
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

df_genes = pd.read_pickle('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/src/data/B72_CMS_10000bp/df_genes_reg0.pkl')
mat = load_npz('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/src/data/B72_CMS_10000bp/hic_sparse.npz')
extract = torch_geometric.utils.from_scipy_sparse_matrix(mat)
df_node_coord = pd.read_pickle('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/src/data/B72_CMS_10000bp/df_node_coord.pkl')

df_node_coord['hic_node_id'] = df_node_coord['hic_node_id'].astype(int)
df_node_coord['chr_num'] = df_node_coord['chr_num'].astype(int)
df_node_coord['chr_start'] = df_node_coord['chr_start'].astype(int)
df_node_coord['chr_finish'] = df_node_coord['chr_finish'].astype(int)
df_node_coord.shape
df_node_coord_w_gene = df_node_coord.merge(df_genes,how='inner',on=['hic_node_id'])
widow_size = 10000

left = extract[0][0].cpu().detach().numpy()
right = extract[0][1].cpu().detach().numpy()
interaction_distances = []
max_per_chr = []
for chr in range(1,19):
    subset_nodes = df_node_coord_w_gene['hic_node_id'][df_node_coord_w_gene['chr_num'] == chr].tolist()
    mask_chr = np.isin(left, subset_nodes) | np.isin(right, subset_nodes)
    chr_distances = (abs(extract[0][1][mask_chr] - extract[0][0][mask_chr])*widow_size).cpu().detach().tolist()
    interaction_distances = interaction_distances + chr_distances
    max_per_chr.append(max(chr_distances))
    print(max_per_chr[-1])

len(interaction_distances)
df_node_coord_w_gene.shape[0]*24

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Set the style for seaborn
sns.set_theme(style="whitegrid", palette="muted")

# Generate dummy data: 1000 data points that are multiples of 10000
data = np.array(interaction_distances)
# Plotting the histogram
plt.figure(figsize=(10, 6))
sns.histplot(data, bins=30, kde=False, color='skyblue')
plt.title('Distribution of Distance of Interactions in Nodes w/Genes', fontsize=16)
plt.xlabel('Distance of Interacting Node', fontsize=14)
plt.ylabel('Frequency', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Enhancing the aesthetics for publication quality
plt.tight_layout()

# Display the plot
plt.show()
plt.savefig('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/reviewer_questions/distribution_of_distances.pdf')

threshold = 0
data = np.array([x for x in interaction_distances if x >= threshold])
# Plotting the histogram
plt.figure(figsize=(10, 6))
sns.histplot(data, bins=300, kde=False, color='skyblue')
plt.title('Distribution of Distance of Interactions in Nodes w/Genes', fontsize=16)
plt.xlabel('Distance of Interacting Node', fontsize=14)
plt.ylabel('Frequency', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylim(0,100000)
plt.xlim(0,20_000_000)

# Enhancing the aesthetics for publication quality
plt.tight_layout()

# Display the plot
plt.show()
plt.savefig('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/reviewer_questions/distribution_of_distances_above_50k_bp.pdf')

np.percentile(interaction_distances, [25, 50, 75, 95])
38_720_000
4_880_000

