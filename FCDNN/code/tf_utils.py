import math
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
from tensorflow.python.framework import ops

# >>>>>> Formatter functions:
def create_placeholders(n_x, n_y):
    """
    Creates the placeholders for the tensorflow session.
    Arguments:
        n_x -- scalar, size of vector
        n_y -- scalar, number of classes (if using OneHot from 0 to 5, then -> 6)
    Returns:
        X   -- placeholder for the data input, of shape [n_x, None] and dtype "float"
        Y   -- placeholder for the input labels, of shape [n_y, None] and dtype "float"
    """
    X = tf.placeholder(tf.float32, shape=(n_x,None), name="X")
    Y = tf.placeholder(tf.float32, shape=(n_y,None), name="Y")
    E = tf.placeholder(tf.float32, name="E") # <<< Epochs for exponential decay
    return X, Y, E


def convert_to_one_hot(Y, C):
    Y = np.eye(C)[Y.reshape(-1)].T
    return Y


def one_hot_matrix(labels, C):
    """
    Creates a matrix where the i-th row corresponds to the ith class number and the jth column corresponds to 
        the jth training example. So if example j had a label i. Then entry (i,j) will be 1. 
    Arguments:
        label s -- vector containing the labels 
        C       -- number of classes, the depth of the one hot dimension
    Returns: 
        one_hot -- one hot matrix
    """
    C = tf.constant(C, name = "C")
    one_hot_matrix = tf.one_hot(labels, depth=C, axis=0)
    sess = tf.Session()
    one_hot = sess.run(one_hot_matrix)
    sess.close()
    return one_hot


def ones(shape):
    """
    Creates an array of ones of dimension shape
    Arguments:
        shape   -- shape of the array you want to create
    Returns:    
        ones    -- array containing only ones
    """
    ones = tf.ones(shape=shape)
    sess = tf.Session()
    ones = sess.run(ones)
    sess.close()
    return ones


def random_mini_batches(X, Y, mini_batch_size = 64, seed = 0):
    """
    Creates a list of random minibatches from (X, Y)
    Arguments:
        X               -- input data, of shape (input size, number of examples)
        Y               -- true "label" vector (containing 0 if cat, 1 if non-cat), of shape (1, number of examples)
        mini_batch_size -- size of the mini-batches, integer
        seed            -- Different mini batches each time
    Returns:
        mini_batches    -- list of synchronous (mini_batch_X, mini_batch_Y)
    """
    
    m = X.shape[1] # number of training examples
    mini_batches = []
    np.random.seed(seed)
    
    # Step 1: Shuffle (X, Y)
    permutation = list(np.random.permutation(m))
    shuffled_X = X[:, permutation]
    shuffled_Y = Y[:, permutation].reshape((Y.shape[0],m))

    # Step 2: Partition (shuffled_X, shuffled_Y). Minus the end case.
    num_complete_minibatches = math.floor(m/mini_batch_size) # number of mini batches of size mini_batch_size in your partitionning
    for k in range(0, num_complete_minibatches):
        mini_batch_X = shuffled_X[:, k * mini_batch_size : k * mini_batch_size + mini_batch_size]
        mini_batch_Y = shuffled_Y[:, k * mini_batch_size : k * mini_batch_size + mini_batch_size]
        mini_batch = (mini_batch_X, mini_batch_Y)
        mini_batches.append(mini_batch)
    
    # Handling the end case (last mini-batch < mini_batch_size)
    if m % mini_batch_size != 0:
        mini_batch_X = shuffled_X[:, num_complete_minibatches * mini_batch_size : m]
        mini_batch_Y = shuffled_Y[:, num_complete_minibatches * mini_batch_size : m]
        mini_batch = (mini_batch_X, mini_batch_Y)
        mini_batches.append(mini_batch)
    
    return mini_batches


def exp_learn_rate_decay(epoch, learning_rate, decay = 0.99, decay_schedule = 30.0):
    return tf.cond(decay_schedule > epoch, lambda: learning_rate, lambda: learning_rate * tf.pow(decay, (epoch - decay_schedule)))


# >>>>>> Costs functions:
def cost(logits, labels):
    """
    Computes the cost using the sigmoid cross entropy
    Arguments:
        logits -- vector containing z, output of the last linear unit (before the final sigmoid activation)
        labels -- vector of labels y (1 or 0)
    Returns:
        cost   -- runs the session of the cost (formula (2))
    """
    z = tf.placeholder(tf.float32, name="z")
    y = tf.placeholder(tf.float32, name="y")
    
    cost = tf.nn.sigmoid_cross_entropy_with_logits(logits=z, labels=y)
    sess = tf.Session()
    cost = sess.run(cost, feed_dict={z : logits, y : labels})
    sess.close()
    
    return cost


def compute_cost(ZL, Y):
    """
    Computes the cost
    Arguments:
        ZL        -- output of forward propagation (output of the last linear unit), of shape (n_features, n_examples)
        Y         -- "true" labels vector placeholder, same shape as Z3
    Returns:
        cost      -- Tensor of the cost function
    """
    logits = tf.transpose(ZL)
    labels = tf.transpose(Y)
    cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=logits, labels=labels))
    return cost


def compute_cost_L2_Regularizer(ZL, Y, parameters, betaReg = 0.01):
    """
    Computes the cost
    Arguments:
        ZL        -- output of forward propagation (output of the last linear unit), of shape (n_features, n_examples)
        Y         -- "true" labels vector placeholder, same shape as Z3
    Returns:
        cost      -- Tensor of the cost function
    """
    L = int(len(parameters)/2 + 1)                                    
    regularizers = 0
    logits = tf.transpose(ZL)
    labels = tf.transpose(Y)
    cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=logits, labels=labels))
    for l in range(1,L):
        regularizers += tf.nn.l2_loss(parameters["W"+str(l)])    
    cost = tf.reduce_mean(cost + betaReg * regularizers)
    return cost


def Regular_Cost(ZL, Y, parameters, betaReg = 0.01, OutputLayer="softmax"):
    """
    Computes the cost
    Arguments:
        ZL           -- Output of forward propagation (output of the last linear unit)
        Y            -- "true" labels vector placeholder, same shape as ZL
        betaReg      -- magnitude of beta regularization, set betaReg=0 to nullify its effects.
        OutputLayer  -- Nature of output layer ("softmax" or "sigmoid")
    Returns:
        cost         -- Tensor of the cost function
    """
    regularizers = 0
    logits = tf.transpose(ZL)
    labels = tf.transpose(Y)
    if   OutputLayer == "softmax" :
        cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=logits, labels=labels))
    elif OutputLayer == "sigmoid" :
        cost = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(logits=logits, labels=labels))
    if betaReg > 0:
        L = int(len(parameters)/2 + 1)
        for l in range(1,L):
            regularizers += tf.nn.l2_loss(parameters["W"+str(l)])
    cost = tf.reduce_mean(cost + betaReg * regularizers)
    return cost


# >>>>>> Parameter functions:
def initialize_parameters(layers_dim):
    """
    Randomly initializes parameters to met the given network structure:
    Arguments:
        layers_dim -- an array whose lenght indicate NN's depth and content indicate neurons in each layer
    Returns:
        parameters -- a dictionary of tensors containing W1, b1, W2, b2, ... , WN, bN
    """
    parameters = {}
    L = len(layers_dim) # number of layers in the network
    for l in range(1, L):
        parameters['W' + str(l)] = tf.get_variable("W"+str(l), shape = [layers_dim[l],layers_dim[l-1]], initializer=tf.contrib.layers.xavier_initializer()) # Removed `seed=1` from initializer EGA 2023/08/23
        parameters['b' + str(l)] = tf.get_variable("b"+str(l), shape = [layers_dim[l],1], initializer=tf.zeros_initializer())
    return parameters


def initialize_given_parameters(Parameters):
    """
    Initializes parameters to Given Values:
    Arguments:
        Parameters -- A dictionary that cointains the weight and the bias for the whole network's neurons
    Returns:
        parameters -- a dictionary of tensors containing W1, b1, W2, b2, ... , WN, bN
    """
    parameters = {}
    L = int(len(Parameters)/2 + 1) # number of layers in the network
    for l in range(1, L):
        parameters['W' + str(l)] = tf.Variable(Parameters['W' + str(l)], trainable=True)
        parameters['b' + str(l)] = tf.Variable(Parameters['b' + str(l)], trainable=True)
    return parameters

def initialize_transfer_parameters(Parameters,n_layers,layers_dim):
    """
    Transfer parameters up to `n_layers`:
    Arguments:
        Parameters -- A dictionary that cointains the weight and the bias for the whole network's neurons
        n_layers   -- an int indicating up to n'th layer to transfer
        layers_dim -- an array whose lenght indicate NN's depth and content indicate neurons in each layer
    Returns:
        parameters -- a dictionary of tensors containing W1, b1, W2, b2, ... , WN, bN
    """
    parameters = {}
    L = len(layers_dim) # number of layers in the network
    # Transfer parameters
    for l in range(1, n_layers):
        parameters['W' + str(l)] = tf.Variable(Parameters['W' + str(l)], trainable=True)
        parameters['b' + str(l)] = tf.Variable(Parameters['b' + str(l)], trainable=True)
    # Initialize new parameters
    for l in range(n_layers, L):
        parameters['W' + str(l)] = tf.get_variable("W"+str(l), shape = [layers_dim[l],layers_dim[l-1]], initializer=tf.contrib.layers.xavier_initializer()) # Removed `seed=1` from initializer EGA 2023/08/23
        parameters['b' + str(l)] = tf.get_variable("b"+str(l), shape = [layers_dim[l],1], initializer=tf.zeros_initializer())
    return parameters




# >>>>>> Propagation functions:
def forward_propagation_DropOut(X, parameters, keepProb=0.7):
    """
    Implements the forward propagation for the model: LINEAR -> ReLU -> ReLU_DropOut -> ... -> LINEAR -> ReLU -> ReLU_DropOut -> LINEAR -> SOFTMAX
    Arguments:
        X          -- input dataset placeholder, of shape (input size, number of examples)
        parameters -- python dictionary containing your parameters "W1", "b1", "W2", ..., "W_N", "b_N" the shapes are given in initialize_parameters
        layers_dim -- Network sizes
        keepProb   -- Probability of keeping a node, if not DropOut, set keepProb=1
    Returns:
        Z_N -- the output of the last LINEAR unit
    """
    L = int(len(parameters)/2 + 1) # number of layers in the network
    for l in range(1, L):
        ZL = tf.add(tf.matmul(parameters["W"+str(l)],X),parameters["b"+str(l)]) # Z1 = np.dot(W1, X) + b1
        relu_layer = tf.nn.relu(ZL)
        X = relu_layer if keepProb == 1 else tf.nn.dropout(relu_layer, keepProb)
    return ZL


# >>>>>> Prediction functions:
def predict(X_dev, Y_dev, Parameters, OutputLayer):
    # Create Placeholders of shape (n_x, n_y)
    (n_x, m) = X_dev.shape
    n_y = Y_dev.shape[0]  
    X, Y, _ = create_placeholders(n_x, n_y)
    # Initialize parameters
    parameters = initialize_given_parameters(Parameters)
    ZLt = forward_propagation_DropOut(X=X, parameters=parameters, keepProb=1.0)
    init = tf.global_variables_initializer() # Initialize all the variables
    sSconf = tf.ConfigProto( intra_op_parallelism_threads=8, inter_op_parallelism_threads=2, allow_soft_placement=True, device_count = {'CPU': 8})
    with tf.Session(config=sSconf) as sess: # Start the session to compute the tensorflow graph
        sess.run(init) # Run the initialization
        if OutputLayer == "softmax":
            Probabilities = tf.nn.softmax(ZLt, axis = 0)
            Prediction = tf.argmax(ZLt)
            Precision = tf.equal(Prediction, tf.argmax(Y))             
        elif OutputLayer == "sigmoid":
            Probabilities = tf.nn.sigmoid(ZLt)
            Prediction = tf.round(Probabilities)
            Precision = tf.equal(Prediction, Y) 
        else:
            raise ValueError(f"Output Layer '{OutputLayer}' not implemented")
        accuracy = tf.reduce_mean(tf.cast(Precision, "float"))
        accuracy = accuracy.eval({X: X_dev, Y: Y_dev})
        print ("Accuracy:\t%.4f"  %  (accuracy), flush=True)
        return Probabilities.eval({X: X_dev}), Prediction.eval({X: X_dev}), accuracy


# >>>>>> Models 
def FCDNN(X_train, Y_train, X_dev, Y_dev, layers_dim, Parameters:dict=None, n_layers:int=None, 
          learning_rate = 0.001, num_epochs = 1500, minibatch_size = 32, print_acc = True,
          generate_plots = True, AppendEach = 5, keepProb=0.7,betaReg=0.01, seed=3, 
          decay = 0.975, decay_schedule = 5.0, OutputFigure='./modelFigs.png'):
    """
    v1.0
    Implements a Multilayered Fully Connected tensorflow Deep Neural Network:
        Linear -> ReLU -> DropOutReLU -> Linear ->  ...  -> Linear -> OutputActivationNeuron(s).
    Arguments:
        X_train         -- Training set X features.
        Y_train         -- Training set Y label.
        X_dev           -- Development set X features.
        Y_dev           -- Development set Y Labels.
        layers_dim      -- list, length = model's depth (total layers), content=neurons per layer
                           Total Neurons in output layer will determine if sigmoid or softmax activation function to use.
        Parameters      -- If Follow Up in Parameters, take them into the new parameters
        n_layers        -- If Transfering parameters, how many to conserve minus 1 (max eq Parameters)
        learning_rate   -- learning rate of the optimization
        num_epochs      -- number of epochs of the optimization loop
        minibatch_size  -- size of a minibatch
        print_acc       -- True to print the Accuracy every 'AppendEach' (def:5)
        generate_plots  -- True to Generate output Accuracy figure and output cost figure
        AppendEach      -- Will recollect models info each N-th epoch iteration (def:5)
        keepProb        -- Probability of keeping a node
        betaReg         -- Loss function with L2 Regularization with beta=betaReg
        seed            -- "Random" seed
        decay           -- Strength of exponential learning rate adjustment (set to 1 to opt out).
        decay_schedule  -- Number of epochs before learning rate adjustment (set to num_epochs to opt out).
    Returns:
        parameters -- Dict: parameters learnt by the model. They can then be used to predict.
        costs      -- List: The costs computed each N-th epoch (default each 5, change with 'AppendEach')
        TrAcc      -- List: The accuracy in the Training    Set computed each N-th epoch (default each 5, change with 'AppendEach')
        DevAcc     -- List: The accuracy in the Development Set computed each N-th epoch (default each 5, change with 'AppendEach')
    """
    ops.reset_default_graph() # to be able to rerun the model without overwriting tf variables
    sSconf  = tf.ConfigProto( intra_op_parallelism_threads=8, inter_op_parallelism_threads=2, allow_soft_placement=True, device_count = {'CPU': 8})
    tf.set_random_seed(seed) # to keep consistent results
    (n_x, m) = X_train.shape # (n_x: n_features, m : n_examples train data set)
    n_y = Y_train.shape[0] # n_y : n_classes
    num_minibatches = int(m / minibatch_size) # number of minibatches of size minibatch_size in the train set
    costs = [] # Cost tracker
    TrAcc = [] # Accuracy in training tracker
    DevAcc = [] # Accuracy in development tracker
    X, Y, E = create_placeholders(n_x, n_y) # Create Placeholders of shape (n_x, n_y)
    # Initialize parameters following setting
    if Parameters == None:
        parameters = initialize_parameters(layers_dim = layers_dim)
    elif n_layers == None:
        parameters = initialize_given_parameters(Parameters)
    else:
        parameters = initialize_transfer_parameters(Parameters,n_layers,layers_dim)
    OutputLayer = "sigmoid" if layers_dim[-1] == 1 else "softmax"
    ZL = forward_propagation_DropOut(X=X, parameters=parameters, keepProb=keepProb)
    ZLt = forward_propagation_DropOut(X=X, parameters=parameters, keepProb=1.0)
    cost = Regular_Cost(ZL=ZL, Y=Y, parameters=parameters, betaReg=betaReg, OutputLayer=OutputLayer)
    optimizer = tf.train.AdamOptimizer(learning_rate = exp_learn_rate_decay( epoch = E, learning_rate = learning_rate, decay = decay, decay_schedule = decay_schedule) ).minimize(cost) # Backpropagation: Define the tensorflow optimizer. Use an AdamOptimizer.
    predicted = tf.equal(tf.argmax(ZLt), tf.argmax(Y)) if OutputLayer == "softmax" else tf.equal(tf.round(tf.nn.sigmoid(ZLt)), Y) # Calculate predictions (1 or 0 for accurate or not)
    accuracy = tf.reduce_mean(tf.cast(predicted, "float")) # mean value of 1 and 0s [0, 1] accuracy
    init = tf.global_variables_initializer() # Initialize all the variables
    with tf.Session(config=sSconf) as sess: # Start the session to compute the tensorflow graph
    #    with tf.Session() as sess: # Start the session to compute the tensorflow graph
        sess.run(init) # Run the initialization
        for epoch in range(num_epochs): # Do the training loop
            epoch_cost = 0. # Defines a cost related to an epoch
            seed += 1
            minibatches = random_mini_batches(X_train, Y_train, minibatch_size, seed)
            for minibatch in minibatches:
                (minibatch_X, minibatch_Y) = minibatch # Select a minibatch
                _ , minibatch_cost = sess.run([optimizer, cost], feed_dict={X: minibatch_X, Y: minibatch_Y, E: epoch}) # Run the session to execute the "optimizer" and the "cost", the feedict should contain a minibatch for (X,Y) and current Epochs E.
                epoch_cost += minibatch_cost / num_minibatches
            if (epoch+1) % AppendEach == 0:
                costs.append(epoch_cost)
                TrAcc.append(accuracy.eval({X: X_train, Y: Y_train}))
                DevAcc.append(accuracy.eval({X: X_dev, Y: Y_dev}))
                if print_acc:
                    print ("epoch" + str(epoch+1).zfill(3) + ": Train Acc:", TrAcc[-1] , flush=True)
                    print ("epoch" + str(epoch+1).zfill(3) + ":   Dev Acc:", DevAcc[-1], flush=True)
        if generate_plots:
            Title=f'lr={str(learning_rate)} L2={str(betaReg)} p(keep)={str(keepProb)}\nminiB={str(minibatch_size)} e={str(num_epochs)}'
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            x_vals = range(AppendEach, num_epochs+1, AppendEach)
            ax.plot(x_vals, TrAcc, color='tab:blue',   label="Training")
            ax.plot(x_vals, DevAcc, color='tab:orange', label="Dev")
            ax.legend(loc='upper left', frameon=False)
            plt.title(Title)
            plt.ylim(0.5,1)
            plt.ylabel('Accuracy')
            plt.xlabel('Epochs')
            plt.savefig(OutputFigure)
            plt.close()
        parameters = sess.run(parameters)
        print ("Train Acc:\t%.4f"  %  (accuracy.eval({X: X_train, Y: Y_train})), flush=True)

        return parameters, costs, TrAcc, DevAcc
