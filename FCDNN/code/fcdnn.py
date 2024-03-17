# Sconda conda activate TensorFlow_CPU_01; cd /mnt/bioadhoc-temp/Groups/RaoLab/Edahi/ForFerhatGit/FCDNN/code
#!/mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/01.Downloaded/anaconda3/envs/TensorFlow_CPU_01/bin/python
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
import tensorflow as tf
from pathlib import Path

from typing import List
from sklearn.metrics import precision_score, recall_score, f1_score, roc_curve, auc

# data_path = '/mnt/bioadhoc-temp/Groups/RaoLab/Edahi/ForFerhatGit/FCDNN/data/C2P'
# data_path = '/mnt/bioadhoc-temp/Groups/RaoLab/Edahi/ForFerhatGit/FCDNN/data/C2B'
def load_dataset(data_prefix: str, subset: str,categories: int = 1, print_shape:bool=False,is_npz:bool=False):
    if subset not in ["Train", "Dev", "Test"]:
        raise ValueError("Invalid value for 'subset'. Allowed values are 'Train', 'Dev', and 'Test'.")
    if is_npz:
        suffix_x = "_"+subset+"_X.npz"
        suffix_y = "_"+subset+"_Y.npz"
        x_data = np.load(data_prefix + suffix_x)['arr_0']
        y_data = np.load(data_prefix + suffix_y)['arr_0']
        y_data_hot = one_hot_matrix(y_data, categories) if categories > 1 else y_data
    else:
        suffix_x = "_"+subset+"_X.txt"
        suffix_y = "_"+subset+"_Y.txt"
        x_data = np.loadtxt(data_prefix + suffix_x,delimiter='\t',dtype = float)
        y_data = np.loadtxt(data_prefix + suffix_y, delimiter='\t', dtype = float)
        y_data_hot = one_hot_matrix(y_data.reshape(y_data.shape[0]), categories) if categories > 1 else y_data.reshape((1, y_data.shape[0]))
    if print_shape:
        print ('X_%s shape:\t%s\nY_%s shape:\t%s' % ( subset, x_data.shape, subset, y_data_hot.shape))
    return x_data, y_data, y_data_hot

# calculate accuracy
def run_predictions(x_data, y_data, y_data_hot, parameters, output_layer, output_file_pred):
    prob, pred, acc = predict(x_data, y_data_hot, parameters, output_layer)
    highest_prob = np.max(prob,axis = 0)
    pred = pred.flatten() if output_layer == 'sigmoid' else pred
    obs = np.squeeze(y_data.astype(int))
    side2prob = pd.DataFrame({'Observed': obs,'Predicted': pred, 'Probability': highest_prob})
    side2prob.to_csv(output_file_pred, index=False)
    return obs, pred, highest_prob, acc, prob

# calculate unbiased metrics
def calculate_metrics(obs, pred, prob, acc, output_layer, output_file_metrics = None, name = None, save_metrics:bool=False, print_metrics:bool=True):
    if output_layer == "sigmoid":
        precision = precision_score(obs, pred)
        recall = recall_score(obs, pred)
        f1 = f1_score(obs, pred)
        fpr, tpr, _ = roc_curve(obs, prob)
        auc_score = auc(fpr, tpr)
        if print_metrics:
            print("Precision:\t%.4f\nRecall: \t%.4f\nF1 score:\t%.4f\nAUC score:\t%.4f"  %  (precision, recall, f1, auc_score) )
    else:
        precision = precision_score(obs, pred, average = "weighted")
        recall = recall_score(obs, pred, average = "weighted")
        f1 = f1_score(obs, pred, average = "weighted")
        fpr, tpr, _ = np.nan, np.nan, np.nan
        auc_score = np.nan
        if print_metrics:
            print("Precision:\t%.4f\nRecall: \t%.4f\nF1 score:\t%.4f\nAUC score:\t%.4f"  %  (precision, recall, f1, auc_score) )
    if save_metrics:
        metrics = pd.DataFrame({'Accuracy': acc,
                                'Precision' : precision,
                                'Recall' : recall,
                                'F1_score' : f1,
                                'AUC_score' : auc_score}, index=[name]).round(4)
        metrics.to_csv(output_file_metrics)
    return fpr, tpr, auc_score, f1

# ptotting
def plot_roc(fpr, tpr, auc_score, f1, name, output_file_roc):
    plt.figure()
    plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.4f)' % auc_score)
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f"Receiver Operating Characteristic (ROC) Curve\n{name}")
    plt.legend(loc='lower right')
    plt.text(0.825, 0.125, f'F1: {f1:.4f}', transform=plt.gca().transAxes)
    plt.savefig(output_file_roc)
    plt.close()

def main(args):
    # define preffix
    if args.name is None:
        # Run the whole dataset!
        # all_files_path = os.path.join(Path(__file__).parents[1], 'all_files_surnames.csv')
        all_files_path = os.path.join(Path(__file__).parents[1], 'filtered_files_surnames.csv')
        names = pd.read_csv(all_files_path, header = None)[7].tolist()
    else:
        names = [args.name]
    for name in names:
        print(f'\n# {name}')
        #define data prefix
        data_prefix = os.path.join(args.data_path,name)
        # define output layer
        output_layer = "sigmoid" if args.output_neurons == 1 else "softmax"
        # declate figures outputs
        output_figure_file = os.path.join(args.output_fig,'accuracy_plot_'+name+'.png')
        output_file_roc = os.path.join(args.output_fig,'roc_curve_'+name+'.png')
        output_file_pred = os.path.join(args.output_pred,'predictions_'+name+'.csv')

        if args.train:
            Parameters=None
            if args.transfer_parameters:
                with open(args.transfer_parameters, 'rb') as f:
                    Parameters = pickle.load(f)
            # declate model outputs (train)
            output_file_model = os.path.join(args.output_model,'fcdnn_model_parameters_'+name+'.pkl')
            # load data for train or test
            x_train, y_train, y_train_hot = load_dataset(data_prefix, 'Train', categories = args.output_neurons, print_shape = args.print_shape, is_npz = args.is_npz)
            x_dev, y_dev, y_dev_hot = load_dataset(data_prefix, 'Dev', categories = args.output_neurons, print_shape = args.print_shape, is_npz = args.is_npz)
            # layers configuration (train)
            layers = [x_train.shape[0]] + [int(x) for x in args.hidden_layers.split(',')] + [int(args.output_neurons)]
            # Train model and save parameters
            parameters, _, _, _ = FCDNN(x_train, y_train_hot, x_dev, y_dev_hot, layers_dim = layers,
                                        Parameters=Parameters, n_layers = args.transfer_n_layers, learning_rate = args.learning_rate,
                                        num_epochs = args.epoch,  minibatch_size = args.mini_batch, print_acc = args.print_acc,
                                        generate_plots = True, AppendEach = args.append_each, keepProb=args.dropout_keep_prob, betaReg=args.beta_reg, seed = args.seed,
                                        decay = args.decay_rate, decay_schedule = args.epoch_decay, OutputFigure=output_figure_file)
            # save parameters
            with open(output_file_model, 'wb') as f:
                pickle.dump(parameters, f, protocol=4)
            # calculate accuracy
            obs, pred, prob, acc, _ = run_predictions(x_dev, y_dev, y_dev_hot, parameters, output_layer, output_file_pred)
            # calculate unbiased metrics
            fpr, tpr, auc_score, f1 = calculate_metrics(obs, pred, prob, acc, output_layer)
        elif args.test:
            if args.cross_test is not None:
                print(f'# --> Tested with parameters from: {args.cross_test}')
                parameters_file_name = 'fcdnn_model_parameters_'+args.cross_test+'.pkl'
                metrics_file_name = 'metrics_'+args.cross_test+'_tested_on_'+name+'.csv'
                predictions_file_name = 'predictions_'+args.cross_test+'_tested_on_'+name+'.csv'
                output_file_roc = os.path.join(args.output_fig,'roc_curve_'+args.cross_test+'_tested_on_'+name+'.png')
            else:
                parameters_file_name = 'fcdnn_model_parameters_'+name+'.pkl'
                metrics_file_name = 'metrics_'+name+'.csv'
                predictions_file_name = 'predictions_'+name+'.csv'
            # declare parameters to load
            parameters_file_model = os.path.join(args.parameters,parameters_file_name)
            # declate metrics outputs
            output_file_metrics = os.path.join(args.output_metrics,metrics_file_name)
            # declate predictions output
            output_file_pred = os.path.join(args.output_pred,predictions_file_name)
            # load data for train or test
            x_test, y_test, y_test_hot = load_dataset(data_prefix, 'Test', categories = args.output_neurons, print_shape = args.print_shape, is_npz = args.is_npz)
            # load parameters
            with open(parameters_file_model, 'rb') as f:
                parameters = pickle.load(f)
            # calculate accuracy
            obs, pred, prob, acc, _ = run_predictions(x_test, y_test, y_test_hot, parameters, output_layer, output_file_pred)
            # calculate unbiased metrics
            fpr, tpr, auc_score, f1 = calculate_metrics(obs, pred, prob, acc, output_layer, output_file_metrics, name, True)
        
        # ptotting
        if output_layer == 'sigmoid':
            plot_roc(fpr, tpr, auc_score, f1, name, output_file_roc)

if __name__ == '__main__':
    from tf_utils import *
    parser = argparse.ArgumentParser(description="FCDNN TensorFlow Train")
    parser.add_argument("--data_path", '-d', type=str, required=True, help="<Path> where Training/Test/Dev files are located")
    parser.add_argument("--name", '-n', type=str, help="Sample name to load", default=None)
    parser.add_argument("--output_fig", '-f', type=str, required=True, help="Output directory for generated figures")
    parser.add_argument("--output_pred", '-p', type=str, required=True, help="Output directory for DEV/TRAIN datasets predictions")
    parser.add_argument("--output_neurons", '-o', type=int, help="Number of output neurons", default=1)
    parser.add_argument("--hidden_layers", '-l', type=str, help="Comma-separated list of hidden layers dimensions", default='200,100,50')
    parser.add_argument("--seed", '-s', type=int, help="Seed to be used for repeatability", default = 1921)
    parser.add_argument("--mini_batch", '-b', type=int, help="Mini batch size", default = 128)
    parser.add_argument("--epoch", '-e', type=int, help="Total training epochs", default = 40)
    parser.add_argument("--epoch_decay", '-ed', type=int, help="Decay rate starts at this epoch", default = 30)
    parser.add_argument("--learning_rate", '-lr', type=float, help="Learning rate", default = 0.0001)
    parser.add_argument("--decay_rate", '-dr', type=float, help="Learning decay rate", default = 0.975)
    parser.add_argument("--beta_reg", '-l2', type=float, help="Beta regularizer value", default = 0.01)
    parser.add_argument("--dropout_keep_prob", '-dp', type=float, help="Dropout probability", default = 0.85)
    parser.add_argument("--train", action='store_true', help="Select if running training. REQ: output_model")
    parser.add_argument("--print_acc", action='store_true', help="Print Accuracy")
    parser.add_argument("--print_shape", action='store_true', help="Print Loaded datasets shape")
    parser.add_argument("--output_model", '-m', type=str, help="Output directory for generated model", default = None)
    parser.add_argument("--test", action='store_true', help="Select if running Test. REQ: parameters & output_metrics")
    parser.add_argument("--parameters", '-pr', type=str, help="Output directory of generated model (Test exclusive)", default=None)
    parser.add_argument("--output_metrics", '-r', type=str, help="Output directory for metrics results (Test exclusive)", default=None)
    parser.add_argument("--cross_test", '-n2', type=str, help="Sample name to test on (if different than the used for model training)", default=None)
    parser.add_argument("--transfer_parameters", '-tp', type=str, help="Path to learnt pkl parameters to transfer from", default=None)
    parser.add_argument("--transfer_n_layers", '-tl', type=int, help="Total layers from transfer_parameters to transfer", default=None)
    parser.add_argument("--append_each", '-ae', type=int, help="Append training info each nth iteration", default=5)
    parser.add_argument("--is_npz", action='store_true', help="Flag to indicate data is stored as numpy binary npz")
    args = parser.parse_args()
    if args.train and args.test:
        raise ValueError("Select only Training (--train) or Testing (--test), not both.")
    elif args.train and args.output_model is None:
        raise ValueError("Training (--train) requires --output_model <str>")
    elif args.test and (args.parameters is None or args.output_metrics is None):
        raise ValueError("Testing (--test) requires both --parameters <str> AND --output_metrics <str>")
    main(args)
else:
    from tf_utils import *
