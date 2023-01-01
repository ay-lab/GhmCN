"""
python
GNNExplainNode.py

Purpose: Load trained networks, mount in device and explain node of interest.

"""

import argparse
import matplotlib.cm as mcm
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import os
import pandas as pd
import re
import random
import sys
import time
import torch
import torch.nn.functional as F
import torch_geometric
import torch_geometric.transforms as T

from datetime import datetime, date
from inspect import signature
from math import sqrt
from matplotlib.colors import to_hex as ToHex
from model_classes_ import GCN_classification, GCN_regression
from scipy.sparse import load_npz
from sklearn.metrics import roc_auc_score, precision_recall_curve, roc_curve, auc 

from torch_geometric.data  import Data
from torch_geometric.nn    import GNNExplainer
from torch_geometric.utils import k_hop_subgraph, to_networkx

def to_cpu_npy(x):
    '''
    Simple helper function to transfer GPU tensors to CPU numpy matrices
    
    Parameters
    ----------
    x [tensor]: PyTorch tensor stored on GPU
    
    Returns
    -------
    new_x [array]: Numpy array stored on CPU
    
    '''
    new_x = x.cpu().detach().numpy()
    return new_x

def PlotCoords(coords,save_file_coords):
    fig, axs =plt.subplots(1,1)
    collabel=("Node", "chr:str-end")
    axs.axis('tight')
    axs.axis('off')
    the_table = axs.table(cellText=coords,colLabels=collabel,loc='center')
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(10)
    plt.savefig(fname = save_file_coords, transparent=True, format="pdf")


def ReScale(data,new_min,new_max,old_min=0,old_max=1):
    NewValue = (((data - old_min) * (new_max - new_min)) / (old_max - old_min)) + new_min
    return NewValue

def NiceCoors(df_node_coord,i):
    NewValue = '-'.join(map(str, df_node_coord[['chr_num','chr_start', 'chr_finish']].iloc[i].tolist())).replace("-",":",1)
    return NewValue

def Midpoint(p1, p2):
    return [(p1[0]+p2[0])/2, (p1[1]+p2[1])/2]

def MidpointY(p1, p2):
    return [(p1[0]+p2[0])/2, ((p1[1]+p2[1])/2)-0.04]

def MidpointX(p1, p2):
    return [((p1[0]+p2[0])/2)-0.06, (p1[1]+p2[1])/2]

def EGA_visualize_subgraph(node_idx, edge_index, edge_mask, explainer,df_genes, y, df_node_coord,
                            threshold=None, node_alpha=None, seed=10, Top10=True, DrawWeight=False,DrawCoords=True):
    assert edge_mask.size(0) == edge_index.size(1)
    subset, edge_index, _, hard_edge_mask = k_hop_subgraph(
        node_idx, explainer.num_hops, edge_index, relabel_nodes=True,
        num_nodes=None, flow='target_to_source')
    
    edge_mask = edge_mask[hard_edge_mask]
    if threshold is not None:
        edge_mask = (edge_mask >= threshold).to(torch.float)
    
    if y is None:
        y = torch.zeros(edge_index.max().item() + 1, device=edge_index.device)
    else:
        y = y[subset].to(torch.float) / y.max().item()
    
    # >>> Gradient color for edges
    EdgeGradientColor = plt.get_cmap('gray',12 if Top10 else 24) 
    
    # >>> Rank top 10 interactions if more than 11 nodes:
    if Top10:
        new_node_idx = [i for i, x in enumerate(subset == node_idx) if x][0]
        edge_mask    = edge_mask[   edge_index[1,:]    == new_node_idx]
        edge_index   = edge_index[:,edge_index[1,:]    == new_node_idx]
        subset_edge = np.argsort(-1*edge_mask.cpu())[0:10].tolist()
        subset_edge.sort()
        subset_node = edge_index[0,subset_edge].tolist() + [new_node_idx]
        subset_node.sort()
        subset      = subset[subset_node]
        edge_index  = edge_index[:,subset_edge]
        edge_mask   = edge_mask[subset_edge]
        y           = y[subset_node]
        subs        = {i: k for k, i in enumerate(edge_index.unique().tolist())}
        edge_index  = torch.LongTensor([[subs[i] for i in edge_index[0,:].tolist() ], [subs[i] for i in edge_index[1,:].tolist() ]]).to(device)
    
    # >>> Renaming node labes if they bear a gene:
    rename_subset = [ df_genes['abbrev'][df_genes['hic_node_id']==i.tolist()].to_list()[0] if any(df_genes['hic_node_id']==i.tolist()) else i.tolist() for i in subset.cpu() ] 
    coordsname    = [ str(i).replace("\n", "+") for i in rename_subset ]
    coords        = [ NiceCoors(df_node_coord,i) for i in subset.cpu().tolist()   ]
    Rename_subset = [ re.sub(":","\n",re.sub("-\d*", "",coords[i])) if str(rename_subset[i]).isnumeric() else rename_subset[i] for i in range(len(rename_subset)) ]
    coords        = np.array([item for pair in zip(coordsname,coords) for item in pair]).reshape(len(subset),2)
    
    # Get edge colors based on the significance gradient
    gnne_rank = np.argsort(np.argsort(-1*edge_mask.cpu())).tolist()
    edge_color = [ToHex(EdgeGradientColor(i)) for i in gnne_rank]
    
    # Generate
    edge_labels = [ str(round(edge_mask[i].tolist(),4)) for i in range(edge_index.size(1)) ]
    data    = Data(edge_index=edge_index, att=edge_mask, edge_color=edge_color, y=y, num_nodes=y.size(0), edge_label = edge_labels, gnne_rank = gnne_rank).to('cpu')
    G       = to_networkx(data, node_attrs=['y'], edge_attrs=['att', 'edge_color', 'edge_label','gnne_rank'])
    mapping = {k: i for k, i in enumerate(Rename_subset)} if DrawCoords else {k: i for k, i in enumerate(rename_subset)}
    G       = nx.relabel_nodes(G, mapping)
    
    node_args                 = set(signature(nx.draw_networkx_nodes).parameters.keys())
    node_kwargs               = {}
    node_kwargs['node_size']  = 800
    node_kwargs['alpha']      = node_alpha
    node_kwargs['cmap']       = 'brg'
    label_args                = set(signature(nx.draw_networkx_labels).parameters.keys())
    label_kwargs              = {}
    label_kwargs['font_size'] = 8
    
    pos = nx.spiral_layout(G)
    ax  = plt.gca()
    for source, target, data in G.edges(data=True):
        ax.annotate(
            '', xycoords='data', xy = pos[target], xytext=pos[source],
            textcoords = 'data', arrowprops=dict(
                arrowstyle      = "->",
                color           = data['edge_color'],
                shrinkA         = sqrt(node_kwargs['node_size']) / 2.0, # Adjust where the arrow beggins
                shrinkB         = sqrt(node_kwargs['node_size']) / 2.0, # Adjust where the arrow ends
                connectionstyle = "arc3,rad=0.05",
            ))
    
    if DrawWeight:
        for source, target, data in G.edges(data=True):
            ax.annotate(data['edge_label'], 
                        xy = pos[target], 
                        xycoords='data', 
                        xytext=Midpoint(pos[source], pos[target]),
                        fontsize=6,
                        color= data['edge_color']
                        )
            ax.annotate(data['gnne_rank'], 
                        xy = pos[target], 
                        xycoords='data', 
                        xytext=MidpointX(pos[source], pos[target]),
                        fontsize=6,
                        color= data['edge_color']
                        )
    
    nx.draw_networkx_nodes(G, pos, node_color=y.tolist(), **node_kwargs)
    nx.draw_networkx_labels(G, pos, **label_kwargs)
    
    return ax, G, coords

# >>> Argparser example
# # python GNNExplainNode.py -ep 800 -ni 3443 -c Naive_CD4T -rf 0 -hm 2 -cr 10000 -df data -mt ./data/Naive_CD4T_CMS_10000bp/saved_runs/model_2021-10-19-at-14-28-28.pt
parser = argparse.ArgumentParser()

parser.add_argument('-c',  '--cell_line',                    default='Naive_CD4T',  type=str)
parser.add_argument('-df', '--data_folder',                  default='data',        type=str)
parser.add_argument('-rf', '--regression_flag',              default=0,             type=int)
parser.add_argument('-rs', '--random_seed',                  default=0,             type=int)
parser.add_argument('-cr', '--chip_resolution',              default=10000,         type=int)
parser.add_argument('-hr', '--hic_resolution',               default=10000,         type=int)
parser.add_argument('-hm', '--histone_modifications',        default=2,             type=int)
parser.add_argument('-mt', '--model_trained',                default='./data/Naive_CD4T_CMS_10000bp/saved_runs/model_2021-10-19-at-14-28-28.pt' ,               type=str)
parser.add_argument('-mo', '--model_origin',                 default='Unknown',     type=str)
parser.add_argument('-ni', '--node_idx',                     default=3443,          type=int)
parser.add_argument('-gn', '--gene_name',                    default=None,          type=str)
parser.add_argument('-ep', '--epochs',                       default=800,           type=int)
parser.add_argument('-od', '--save_dir',                     default=None,          type=str)
parser.add_argument(       '--noTop10',                      action='store_false')
parser.add_argument(       '--noPlotTable',                  action='store_false')
parser.add_argument(       '--DrawWeight',                   action='store_true')
parser.add_argument(       '--DrawCoords',                   action='store_true')
parser.add_argument(       '--UsePredExp',                   action='store_true')
parser.add_argument('-cn', '--num_graph_conv_layers',        default=2,      type=int)
parser.add_argument('-gs', '--graph_conv_layer_size',        default=256,    type=int)
parser.add_argument('-ln', '--num_lin_layers',               default=3,      type=int)
parser.add_argument('-ls', '--lin_hidden_layer_size',        default=256,    type=int)


args                  = parser.parse_args()
cell_line             = args.cell_line
data_folder           = args.data_folder
regression_flag       = args.regression_flag
random_seed           = args.random_seed

chip_res              = args.chip_resolution
hic_res               = args.hic_resolution
num_hm                = args.histone_modifications
path_to_model         = args.model_trained
model_origin          = args.model_origin
node_idx              = args.node_idx
epochs                = args.epochs
gene_name             = args.gene_name
Top10                 = args.noTop10
PlotTable             = args.noPlotTable
DrawWeight            = args.DrawWeight
DrawCoords            = args.DrawCoords
UsePredExp            = args.UsePredExp
num_graph_conv_layers = args.num_graph_conv_layers
graph_conv_embed_size = args.graph_conv_layer_size
num_lin_layers        = args.num_lin_layers
lin_hidden_size       = args.lin_hidden_layer_size

num_feat = int((hic_res/chip_res)*num_hm)

if regression_flag == 0:
    num_classes = 2
    task = 'Classification'
else:
    num_classes = 1
    task = 'Regression'

# random_seed = random.randint(0,10000)
random.seed(random_seed)
np.random.seed(random_seed)
torch.manual_seed(random_seed)


# >>> Set Device
dev       = "cuda" if torch.cuda.is_available() else "cpu"
device    = torch.device(dev)  

# >>> Load input files
base_path               = os.path.join(os.getcwd(), data_folder,cell_line)
save_dir                = os.path.join(base_path, 'nodes_explained') if args.save_dir == None else args.save_dir
hic_sparse_mat_file     = os.path.join(base_path, 'hic_sparse.npz')
np_nodes_lab_genes_file = os.path.join(base_path, 'np_nodes_lab_genes_reg' + str(regression_flag) + '.npy')
np_hmods_norm_all_file  = os.path.join(base_path, 'np_hmods_norm_chip_' + str(chip_res) + 'bp.npy')
df_genes_file           = os.path.join(base_path, 'df_genes_reg' + str(regression_flag) + '.pkl')
dict_gene_dups_reg_file = os.path.join(base_path, 'dict_gene_dups_reg' + str(regression_flag) + '.pkl')
df_node_coord_file      = os.path.join(base_path, 'df_node_coord.pkl')
predictions_file        = os.path.join(os.path.dirname(args.model_trained) , os.path.basename(args.model_trained).replace("_","_predictions_").replace(".pt",".csv"))
df_genes                = pd.read_pickle(df_genes_file)
dict_gene_dups_reg      = pd.read_pickle(dict_gene_dups_reg_file)
df_node_coord           = pd.read_pickle(df_node_coord_file)
predictions             = pd.read_csv(predictions_file)

# >>> Update genes dictionary to include all gene names:
pd.set_option('mode.chained_assignment',None)
for i in dict_gene_dups_reg.keys():
    if any(df_genes['hic_node_id']==i):
        df_genes['abbrev'][df_genes['hic_node_id']==i ] = '\n'.join(dict_gene_dups_reg[i][0])

df_genes = df_genes.drop(columns = ['gene_catalog_name' , 'expression_lvl', 'connected'])

# >>> Drop unused columns from df_node_coord:
df_node_coord = df_node_coord.drop(columns=['gene_catalog_name','Abbrev','expression_lvl'])
df_node_coord = df_node_coord.astype(int)

# >>> Attempt finding queried gene name:
#     WARNING: if more than one, only the first occurence will be used.
if gene_name is not None and any(df_genes['abbrev'].str.contains(gene_name)):
    query    = df_genes[df_genes['abbrev'].str.contains(gene_name)]['hic_node_id'].to_list()
    node_idx = query[0]
    if len(query)>1:
        print('\n\nWARNING: More than one matching _NodeID_ was found using queried _GeneName_ "' + gene_name + '":')
        print(f'\tNode\tGene(s)\n\t----\t-------')
        for n, g in zip(query , [ df_genes['abbrev'][df_genes['hic_node_id']== i ].iloc[0].replace("\n", ", ") for i in query ]):
            print(f'\t{n}\t{g}')
        
        print('\n\tWhen more than one _NodeID_ is matched, the first occurence will be used: ' + str(query[0]))
        print('\tIf not the Gene of Interest, Try using one of the suggested Node IDs.\n')
elif gene_name is not None:
    sys.exit('Could not find gene name: "' + gene_name + '" -- Check Spelling')

# >>> Attempt finding selected node_idx gene name:
NodeName = df_genes['abbrev'][df_genes['hic_node_id']== node_idx ].iloc[0] if any(df_genes['hic_node_id']==node_idx) else "Enhancer"
NodeName = NodeName.replace("\n", "+")

# >>> Set PDF file output
save_file_graph  = os.path.join(save_dir, 
                                str(cell_line) + 
                                '.Reg_' + str(regression_flag) + 
                                '.Top10_' + str(Top10) + 
                                '.PredY_' + str(UsePredExp) + 
                                '.e' + str(epochs).zfill(4) + 
                                '.Model_' + str(model_origin) + 
                                '.nID_' + str(node_idx) + '.' + 
                                NodeName + '_graph.pdf')
save_file_coords = os.path.join(save_dir, 
                                str(cell_line) + 
                                '.Reg_' + str(regression_flag) + 
                                '.Top10_' + str(Top10) + 
                                '.PredY_' + str(UsePredExp) + 
                                '.e' + str(epochs).zfill(4) + 
                                '.Model_' + str(model_origin) + 
                                '.nID_' + str(node_idx) + '.' + 
                                NodeName + '_coords.pdf')

# >>> Define model inputs
# Sparse Matrix
mat            = load_npz(hic_sparse_mat_file)
# Histone Modifications
allNodes_hms   = np.load(np_hmods_norm_all_file)
if cell_line == "RAVCrep1" or cell_line == "RAVCrep2":
    row = np.array([allNodes_hms[-1][0]+1, 0, 0])
    allNodes_hms = np.append([row],allNodes_hms,axis= 0)
hms            = allNodes_hms[:, 1:] #only includes features, not node ids
X              = torch.tensor(hms).float().reshape(-1, num_feat) 
allNodes       = allNodes_hms[:, 0].astype(int)
# Label genes
geneNodes_labs = np.load(np_nodes_lab_genes_file)
geneNodes      = geneNodes_labs[:, -2].astype(int)
allLabs        = -1*np.ones(np.shape(allNodes))

targetNode_mask = torch.tensor(geneNodes).long()

# >>> Define Predicted Ys:
if UsePredExp:
    allLabs = -1*np.ones(np.shape(allNodes))
    allLabs[np.array(predictions['Node ID'].tolist())] = np.array(predictions['Predicted Label'].tolist())
    Y = torch.tensor(allLabs).long()
else:
    # >>> Define Real Y
    if regression_flag == 0:
        geneLabs           = geneNodes_labs[:, -1].astype(int)
        allLabs[geneNodes] = geneLabs
        Y                  = torch.tensor(allLabs).long()
    else:
        geneLabs           = geneNodes_labs[:, -1].astype(float)
        allLabs[geneNodes] = geneLabs
        Y                  = torch.tensor(allLabs).float()


# >>> Build network
extract = torch_geometric.utils.from_scipy_sparse_matrix(mat)
data    = torch_geometric.data.Data(edge_index = extract[0], edge_attr = extract[1], x = X, y = Y)
G       = data

# >>> Define convolutional and linear layer input/output sizes
graph_conv_layer_sizes = [num_feat] + \
    [int(max(graph_conv_embed_size, lin_hidden_size)) \
          for i in np.arange(1, num_graph_conv_layers, 1)] + [lin_hidden_size]

lin_hidden_sizes = [graph_conv_layer_sizes[-1]] + \
    [int(max(lin_hidden_size, num_classes)) \
          for i in np.arange(1, num_lin_layers, 1)] + [num_classes]

# >>> Instantiate neural network model, choose optimizer, and print model parameters
if regression_flag == 0:
    model = GCN_classification(num_feat, num_graph_conv_layers, graph_conv_layer_sizes, num_lin_layers, lin_hidden_sizes, num_classes)
else:
    model = GCN_regression(num_feat, num_graph_conv_layers, graph_conv_layer_sizes, num_lin_layers, lin_hidden_sizes, num_classes)

# >>> Mount model into device
model.load_state_dict(torch.load(path_to_model, map_location=device) )
model.to(device)
graph  = G.to(device)

# >>> Run Explainer
explainer = GNNExplainer(model, epochs=epochs, return_type='log_prob',num_hops=1)
node_feat_mask, edge_mask = explainer.explain_node(node_idx, graph.x, graph.edge_index)

# >>> Visualize Node:
ax, G, coords = EGA_visualize_subgraph(node_idx, graph.edge_index,edge_mask,explainer,df_genes,data.y, df_node_coord, node_alpha=0.5, Top10=Top10, DrawWeight = DrawWeight, DrawCoords = DrawCoords)
plt.savefig(fname = save_file_graph, transparent=True, format="pdf")
plt.close()
if PlotTable:
    PlotCoords(coords,save_file_coords)
    plt.close()

exit()

