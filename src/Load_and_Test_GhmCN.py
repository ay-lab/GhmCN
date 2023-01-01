"""
python
Load_and_Test_GhmCN.py

Purpose: Load trained networks and test how given datasets predicts gene expression.

Usage: python ./Load_etTest_EGA.py [-c <str>] [-rf <int>] [-mo <str>]

Arguments:
    '-c', '--cell_line', default='E116', type=str
    '-rf', '--regression_flag', default=1 (1 = regression; 0 = classification), type=int
    '-cr', '--chip_resolution', default=10000, type=int
    '-hr', '--hic_resolution',  default=10000, type=int
    '-hm', '--histone_modifications',  default=2, type=int
    '-hc', '--highest_chromosome_excluding',  default=23, type=int
    '-mt', '--model_trained', type=str
    '-mo', '--model_origin',  type=str

# # python Load_etTest_EGA.py -c Naive_CD4T -rf 1 -hm 2 -cr 10000 -df data -mt /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/01.Downloaded/GC-MERGE/src/data/Naive_CD4T_CMS_10000bp/saved_runs/model_2021-10-19-at-14-49-45.pt

Processed inputs: 
    In ./data/cell_line subdirectory:
        ./hic_sparse.npz: Concatenated Hi-C matrix in sparse CSR format
        ./np_nodes_lab_genes_reg[rf].npy: Numpy array stored in binary format, where 
            rf denotes regression flag (rf = 1 for regression, 0 for classification);
            2-column array that stores IDs of nodes corresponding to genes
            and the node label (expression level)
        ./np_hmods_norm_chip_10000bp.npy: Numpy array stored in binary format;
            (F+1)-column array where the 0th column contains node IDs
            and columns 1..F contain feature values, where F = total number of features
        ./df_genes_reg[rf].pkl: Pandas dataframe stored in .pkl format, where 
            rf denotes regression flag (rf = 1 for regression, 0 for classification);
            5-column dataframe, where columns = [ENSEMBL ID, 
            gene name abbreviation, node ID, expression level, connected status]
    *Note: Users can prepare these files or use process_inputs.py script provided

Outputs:  
    In ./data/cell_line/saved_runs subdirectory:
        model_[date_and_time].pt: Model state dictionary stored in .pt (PyTorch) format
        model_predictions_[date_and_time].csv: Predictions for each gene with the following columns:
            Classification: [Dataset, Node ID, ENSEMBL ID, 
                gene name abbreviation, true label, predicted label, classification [TP/TN/FP/FN]]
            Regression: [Dataset, Node ID, ENSEMBL ID, 
                gene name abbreviation, true expression, predicted expression]
                *Note: Expression values are obtained by taking the base-10 logarithm
                    of the RNA-seq counts and adding a pseudocount of 1 prior to taking the logarithm
        model_[date_and_time]_info.txt: Text file containing summary of model
            evaluation metrics as well as hyperparameter settings

"""

import argparse
import csv
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import random
import time
import torch
import torch_geometric

from   datetime import datetime, date
from   scipy.sparse import load_npz
from   scipy.stats import pearsonr
from   sklearn.metrics import roc_auc_score, precision_recall_curve, roc_curve, auc, f1_score 

from    model_classes_ import GCN_classification, GCN_regression


def eval_model_classification(model, graph, targetNode_mask,test_idx, path_to_model):
    '''
    Runs fully trained classification model and compute evaluation statistics

    Parameters
    ----------
    model  [GCN_classification]: Instantiation of model class
    graph      [PyG Data class]: PyTorch Geometric Data object representing the graph
    targetNode_mask    [tensor]: Mask ensuring model only trains on nodes with genes
    train_idx           [array]: Node IDs

    Returns
    -------
    test_AUROC  [float]: Test set AUROC score;
    test_AUPR   [float]: Test set AUPR score
    test_pred   [array]: Test set predictions;
    test_labels [array]: Test set labels;
    '''
    
    # Import model to device:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.load_state_dict(torch.load(path_to_model, map_location=device) )
    model.to(device)
    graph  = graph.to(device)
    
    # Establish labels to test on:
    test_labels = to_cpu_npy(graph.y[targetNode_mask[test_idx]])
    
    # Evaluate model:
    model.eval()
    train_status=False
    forward_scores = model(graph.x.float(), graph.edge_index, train_status)[targetNode_mask]
    
    test_scores = forward_scores[test_idx]
    test_softmax, test_pred = model.calc_softmax_pred(test_scores) 
    
    test_softmax = to_cpu_npy(test_softmax)
    test_pred    = to_cpu_npy(test_pred)
    test_AUROC   = roc_auc_score(test_labels, test_softmax[:,1], average="micro")
    test_precision, test_recall, thresholds = precision_recall_curve(test_labels, test_softmax[:,1])
    test_AUPR    = auc(test_recall, test_precision)
    test_F1      = f1_score(test_labels, test_pred, average="micro")
    
    return test_AUROC, test_AUPR, test_pred, test_labels, test_softmax, test_F1

def eval_model_regression(model, graph, targetNode_mask, test_idx, path_to_model):
    '''
    Runs fully trained regression model and compute evaluation statistics

    Parameters
    ----------
    model  [GCN_classification]: Instantiation of model class
    graph      [PyG Data class]: PyTorch Geometric Data object representing the graph
    targetNode_mask    [tensor]: Mask ensuring model only trains on nodes with genes
    train_idx           [array]: Node IDs

    Returns
    -------
    test_pearson  [float]: PCC for test set;
    test_pred     [array]: Test set predictions;
    test_labels   [array]: Test set labels (expression values);

    '''
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.load_state_dict(torch.load(path_to_model, map_location=device) )
    model.to(device)
    graph  = graph.to(device)
    
    model.eval()
    train_status=False
    
    forward_scores = model(graph.x.float(), graph.edge_index, train_status)[targetNode_mask]
    
    test_scores = forward_scores[test_idx]
    test_pred = to_cpu_npy(test_scores)
    test_labels = to_cpu_npy(graph.y[targetNode_mask[test_idx]])
    test_pearson = calc_pearson(test_pred, test_labels)
    
    return test_pearson, test_pred, test_labels


def calc_pearson(scores, targets):
    '''
    Calculates Pearson correlation coefficient (PCC) between predicted \
        expression levels and true expression levels

    Parameters
    ----------
    scores [array]: Predicted expression levels
    targets [array]: True expression levels

    Returns
    -------
    pcc [float]: Pearson correlation coefficient

    '''
    pcc, _ = pearsonr(scores, targets)
    
    return pcc


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


###Set hyperparameters and miscellaneous options
#  -c B00 -rf 1 -hm 2 -cr 10000 -df data -mo B00 -mt /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/01.Downloaded/GC-MERGE/src/data/B00_CMS_10000bp/saved_runs/model_2021-05-24-at-13-11-12.pt
parser = argparse.ArgumentParser()

parser.add_argument('-c',  '--cell_line',                    default='Naive_CD4T',  type=str)
parser.add_argument('-df', '--data_folder',                  default='data',        type=str)
parser.add_argument('-rf', '--regression_flag',              default=0,             type=int)
parser.add_argument('-rs', '--random_seed',                  default=0,             type=int)
parser.add_argument('-cr', '--chip_resolution',              default=10000,         type=int)
parser.add_argument('-hr', '--hic_resolution',               default=10000,         type=int)
parser.add_argument('-hm', '--histone_modifications',        default=2,             type=int)
parser.add_argument('-mt', '--model_trained',                default='./data/Naive_CD4T_CMS_10000bp/saved_runs/model_2021-11-29-at-21-23-15.pt' ,               type=str)
parser.add_argument('-mo', '--model_origin',                 default='Naive_CD4T',  type=str)
parser.add_argument('-sv', '--save_dir',                     default=None,          type=str)

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
cell_from_model       = args.model_origin

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


###Initialize start time
start_time    = time.time()
today         = date.today()
mdy           = today.strftime("%Y-%m-%d")
clock         = datetime.now()
hms           = clock.strftime("%H-%M-%S")
hm            = clock.strftime("%Hh-%Mm")
hm_colon      = clock.strftime("%H:%M")
date_and_time = mdy + '-at-' + hms

###Test for GPU availability
dev       = "cuda" if torch.cuda.is_available() else "cpu"
device    = torch.device(dev)  

###Load input files
base_path               = os.getcwd()
save_dir                = os.path.join(base_path, data_folder, cell_line, 'saved_predictions')    if args.save_dir is None else args.save_dir
hic_sparse_mat_file     = os.path.join(base_path, data_folder, cell_line, 'hic_sparse.npz')
np_nodes_lab_genes_file = os.path.join(base_path, data_folder, cell_line, 'np_nodes_lab_genes_reg' + str(regression_flag) + '.npy')
np_hmods_norm_all_file  = os.path.join(base_path, data_folder, cell_line, 'np_hmods_norm_chip_'    + str(chip_res)        + 'bp.npy')
df_genes_file           = os.path.join(base_path, data_folder, cell_line, 'df_genes_reg'           + str(regression_flag) + '.pkl')
df_genes                = pd.read_pickle(df_genes_file)

# ###Print model specifications
# model_specs_info = f"""
# # {' '*14}Model Specifications Info
# # {'-'*50}
# #             Program name: run_models #os.path.basename(__file__)
# #              Date & time: {date_and_time}
# #                Cell line: {cell_line}
# #                     Task: {task}
# #          ChIP resolution: {str(chip_res)}
# # {'-'*50}

# """
# print(model_specs_info)


###Define model inputs
# print("\n# >>>>>> Loading Files:\n-Sparse Matrix")
mat            = load_npz(hic_sparse_mat_file)
# print("-Features (Histone Modifications)")
allNodes_hms   = np.load(np_hmods_norm_all_file)
hms            = allNodes_hms[:, 1:] #only includes features, not node ids
X              = torch.tensor(hms).float().reshape(-1, num_feat) 
allNodes       = allNodes_hms[:, 0].astype(int)
# print("-Nodes label genes")
geneNodes_labs = np.load(np_nodes_lab_genes_file)
geneNodes      = geneNodes_labs[:, -2].astype(int)
allLabs        = -1*np.ones(np.shape(allNodes))

targetNode_mask = torch.tensor(geneNodes).long()

if regression_flag == 0:
    geneLabs           = geneNodes_labs[:, -1].astype(int)
    allLabs[geneNodes] = geneLabs
    Y                  = torch.tensor(allLabs).long()
else:
    geneLabs           = geneNodes_labs[:, -1].astype(float)
    allLabs[geneNodes] = geneLabs
    Y                  = torch.tensor(allLabs).float()

extract = torch_geometric.utils.from_scipy_sparse_matrix(mat)
data    = torch_geometric.data.Data(edge_index = extract[0], edge_attr = extract[1], x = X, y = Y)
G       = data

###Define convolutional and linear layer input/output sizes
graph_conv_layer_sizes = [num_feat] + \
    [int(max(graph_conv_embed_size, lin_hidden_size)) \
          for i in np.arange(1, num_graph_conv_layers, 1)] + [lin_hidden_size]

lin_hidden_sizes = [graph_conv_layer_sizes[-1]] + \
    [int(max(lin_hidden_size, num_classes)) \
          for i in np.arange(1, num_lin_layers, 1)] + [num_classes]

###Instantiate neural network model, choose optimizer, and print model parameters
if regression_flag == 0:
    model = GCN_classification(num_feat, num_graph_conv_layers, graph_conv_layer_sizes, num_lin_layers, lin_hidden_sizes, num_classes)
else:
    model = GCN_regression(num_feat, num_graph_conv_layers, graph_conv_layer_sizes, num_lin_layers, lin_hidden_sizes, num_classes)

###Randomize node order and split into 70%/15%/15% training/validation/test sets
pred_idx_shuff = torch.randperm(targetNode_mask.shape[0])

# >>> Randomize if sample tested is the same cell-type as the one used to train the model
#   > Otherwise use all the dataset
if cell_from_model == cell_line:
    fin_valid    = np.floor(0.85*pred_idx_shuff.shape[0]).astype(int)
    test_idx     = pred_idx_shuff[fin_valid:]
    test_gene_ID = targetNode_mask[test_idx].numpy()
else:
    test_idx     = pred_idx_shuff[:]
    test_gene_ID = targetNode_mask[test_idx].numpy()

# ###Instantiate neural network model, choose optimizer, and print model parameters

### For classification:
if regression_flag == 0:
    ### Evaluate model
    test_AUROC, test_AUPR, test_pred, test_labels, test_softmax, test_F1 =  eval_model_classification(model, G, targetNode_mask, test_idx, path_to_model)
    
    ### Plot AU-ROC
    fpr, tpr, thresholds = roc_curve(test_labels, test_softmax[:,1])
    plt.clf()
    fig = plt.figure()
    fig.suptitle('AU-ROC', fontweight='bold')
    ax  = fig.add_subplot(1, 1, 1)
    fig.subplots_adjust(top=0.85)
    ax.plot(fpr,    tpr,                    color='tab:orange', label="Model")
    ax.legend(loc='upper left', frameon=False)
    plt.ylabel('TPR')
    plt.xlabel('FPR')
    Title = "Area Under the ROC: " + str(test_AUROC.round(4)) + "\nCelltype for (Training,Testing) = ("+ str(cell_from_model)+ " , " + str(cell_line) + ")"
    plt.title(Title)
    plt.savefig(os.path.join(save_dir, 'Load_etTested_model_prediction_rf' + str(regression_flag) + '_fromCell_'+ cell_from_model + '_toCell_' + cell_line + '_Time_' + date_and_time + '_AUROC.png'))
    
    ### Save Metrics
    # AUC-FPR,TPR,CellName
    AUC_plot = np.array([item for pair in zip(fpr,tpr,[cell_line] * len(fpr)) for item in pair]).reshape(len(fpr),3)
    with open(os.path.join(save_dir, 'Load_etTested_model_prediction_rf' + str(regression_flag) + '_fromCell_'+ cell_from_model + '_toCell_' + cell_line + '_Time_' + date_and_time + '_AUC_FPR_TPR.csv'), "wt") as f:
        fileWriter = csv.writer(f)
        fileWriter.writerows(AUC_plot)
    
    # General Metrics
    test_metrics = [test_gene_ID, test_pred, test_labels, test_AUROC, test_AUPR, ['na']]
    np.save(os.path.join(save_dir, 'Load_etTested_model_prediction_rf' + str(regression_flag) + '_fromCell_'+ cell_from_model + '_toCell_' + cell_line + '_Time_' + date_and_time + '_metrics.npy'), test_metrics)

    # Full Metrics
    dataset_list = [test_metrics]
    df_full_metrics = pd.DataFrame(columns=['Dataset','Node ID','True Label','Predicted Label','Classification'])
    dataset_metrics = dataset_list[0]
    partial_metrics = pd.DataFrame()
    partial_metrics['Node ID'] = dataset_metrics[0]
    partial_metrics['True Label'] = dataset_metrics[2]
    partial_metrics['Predicted Label'] = dataset_metrics[1]
    partial_metrics['Classification'] = dataset_metrics[1]*1 + dataset_metrics[2]*2
    partial_metrics['Classification'].replace(to_replace=0, value='TN', inplace=True)
    partial_metrics['Classification'].replace(to_replace=1, value='FP', inplace=True)
    partial_metrics['Classification'].replace(to_replace=2, value='FN', inplace=True)
    partial_metrics['Classification'].replace(to_replace=3, value='TP', inplace=True)
    partial_metrics['Dataset'] = 'modelTesting'
    df_full_metrics = df_full_metrics.append(partial_metrics)
    df_gene_names = df_genes.iloc[:,:3]
    df_gene_names = df_gene_names.rename(columns={"gene_catalog_name": "ENSEMBL_ID", "abbrev": "Abbreviation", "hic_node_id" : 'Node ID'})
    df_full_metrics = pd.merge(df_full_metrics, df_gene_names, how='inner', on='Node ID')
    df_full_metrics = df_full_metrics[df_full_metrics.columns[[0,1,5,6,2,3,4]]]
    
### For regression:
elif regression_flag == 1:
    ### Evaluate model
    test_pearson, test_pred, test_labels = eval_model_regression(model, G, targetNode_mask, test_idx, path_to_model)
    test_metrics = [test_gene_ID, test_pred, test_labels, test_pearson, ['na']]
    np.save(os.path.join(save_dir, 'Load_etTested_model_prediction_rf' + str(regression_flag) + '_fromCell_'+ cell_from_model + '_toCell_' + cell_line + '_Time_' + date_and_time + '_metrics.npy'), test_metrics)
    
    ### Plot Pearson Correlaiton
    m, b = np.polyfit(test_labels,    test_pred, 1) #obtain m (slope) and b(intercept) of linear regression line
    fig = plt.figure()
    fig.suptitle('PCC', fontweight='bold')
    ax = fig.add_subplot(1, 1, 1)
    fig.subplots_adjust(top=0.85)
    ax.plot(test_labels,    test_pred,  '.',                  color='tab:orange')
    ax.plot(test_labels,    m*test_labels+b)        #Add slope
    ax.legend(loc='upper left', frameon=False)
    plt.ylabel('Predicted Expression')
    plt.xlabel('Observed Expression')
    plt.ylim((-0.3,5.2))
    plt.xlim((-0.3,5.2))
    Title = "Pearson Corr Coef:  " + str(test_pearson.round(4)) + "\nCelltype for (Training,Testing) = ("+ str(cell_from_model)+ " , " + str(cell_line) + ")"
    plt.title(Title)
    plt.savefig(os.path.join(save_dir, 'Load_etTested_model_prediction_rf' + str(regression_flag) + '_fromCell_'+ cell_from_model + '_toCell_' + cell_line + '_Time_' + date_and_time + '_PCC.png'))

    dataset_list = [test_metrics]
    df_full_metrics = pd.DataFrame(columns=['Dataset','Node ID','True Label','Predicted Label'])
    
    dataset_metrics = dataset_list[0]
    partial_metrics = pd.DataFrame()
    partial_metrics['Node ID'] = dataset_metrics[0]
    partial_metrics['True Label'] = dataset_metrics[2]
    partial_metrics['Predicted Label'] = dataset_metrics[1]
    partial_metrics['Dataset'] = 'modelTesting'
    
    df_full_metrics = df_full_metrics.append(partial_metrics)
    
    df_gene_names = df_genes.iloc[:,:3]
    df_gene_names = df_gene_names.rename(columns={"gene_catalog_name": "ENSEMBL_ID", "abbrev": "Abbreviation",
                                  "hic_node_id" : 'Node ID'})
    df_full_metrics = pd.merge(df_full_metrics, df_gene_names, how='inner', on='Node ID')
    df_full_metrics = df_full_metrics[df_full_metrics.columns[[0,1,4,5,2,3]]]

print('\nModel Performance:\nTained on: ' + str(cell_from_model) + '\nTested on: ' + str(cell_line))
if regression_flag == 0:
    print('Test AUROC:', test_AUROC)
    print('Test  AUPR:', test_AUPR)
    print('Test    F1:', test_F1)
elif regression_flag == 1:
    print('Test   PCC:', test_pearson)


### Save model predictions CSV file
df_full_metrics_filename = os.path.join(save_dir, 'Load_etTested_model_prediction_rf' + str(regression_flag) + '_fromCell_'+ cell_from_model + '_toCell_' + cell_line + '_Time_' + date_and_time + '_metrics.csv')
df_full_metrics.to_csv(df_full_metrics_filename, index=False)

model_info_filename = os.path.join(save_dir,'Load_etTested_model_prediction_rf' + str(regression_flag) + '_fromCell_'+ cell_from_model + '_toCell_' + cell_line + '_Time_' + date_and_time + '_info.txt')
f = open(model_info_filename, 'w')
f.write('         Model Path: ' + str(path_to_model) + '\n')
f.write('               Task: ' + task + '\n')
f.write('Tested on Cell line: ' + cell_line + '\n')
f.write('        Performance:\n')
if regression_flag == 0:
    f.write('           Test AUROC: ' + str(test_AUROC) + '\n')
    f.write('           Test  AUPR: ' + str(test_AUPR)  + '\n')
    f.write('           Test    F1: ' + str(test_F1)    + '\n')
elif regression_flag == 1:
    f.write('           Test   PCC: ' + str(test_pearson) + '\n\n')
f.close()


