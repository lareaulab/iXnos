#from __future__ import print_function

import sys
import os
import glob
import time
import copy
import random
import pickle
import numpy as np
import theano
import theano.tensor as T

import lasagne
from rp_predict import process as proc
from rp_predict import plot

class LSTM(object):
    def __init__(
            self, X_tr, y_tr, X_te, y_te, X_val=False, y_val=False, 
            name="my_nn", out_dir="lasagne_nn", reloaded=False):
        self.name = name
        self.out_dir = out_dir
        self.num_epochs = 0
        self.train_err_by_epoch = []
        self.test_err_by_epoch = []

class RegressionMLP(object):
    def __init__(
            self, X_tr, y_tr, X_te, y_te, X_val=False, y_val=False, 
            name="my_nn", out_dir="lasagne_nn", reloaded=False,
            nonnegative=False):
        self.name = name
        self.out_dir = out_dir
        self.num_epochs = 0
        self.train_err_by_epoch = []
        self.test_err_by_epoch = []
        self.nonnegative = nonnegative
        #Make file structure for nn storage
        if not reloaded:
            self.make_out_dir(self.out_dir)
            self.make_nn_dir(self.out_dir, self.name)
            self.make_init_data_dir(self.out_dir, self.name)
        #Store data matrices
        self.X_tr = X_tr
        self.y_tr = y_tr
        self.X_te = X_te
        self.y_te = y_te
        self.X_val = X_val if X_val else self.X_te
        self.y_val = y_val if y_val else self.y_te
        #Check conformity of data matrices
        self.check_match_dtype([self.X_tr, self.X_te, self.X_val])
        self.check_match_ndim([self.X_tr, self.X_te, self.X_val])
        self.check_match_dims([self.X_tr, self.X_te, self.X_val])
        #Store properties of data matrices
        self.input_dtype = self.X_tr.dtype
        self.input_ndim = self.X_tr.ndim
        self.input_dims = [None] + list(self.X_tr.shape[1:])
        #Check conformity of output matrices
        self.check_match_dtype([self.y_tr, self.y_te, self.y_val])
        self.check_match_ndim([self.y_tr, self.y_te, self.y_val])
        self.check_match_dims([self.y_tr, self.y_te, self.y_val])
        #Store properties of output matrices
        self.output_dtype = self.y_tr.dtype
        self.output_ndim = self.y_tr.ndim
        self.output_dims = [None] + list(self.y_tr.shape[1:])
        #Check that output matrices are 2D
        self.check_ndarray_dim(self.y_tr, 2)
        self.check_ndarray_dim(self.y_te, 2)
        self.check_ndarray_dim(self.y_val, 2)
        #Store more matrix parameters
        self.num_outputs = self.output_dims[1] 

    def get_nonlin_fn(self, nonlinearity):
        if nonlinearity == "tanh":
            return lasagne.nonlinearities.tanh
        elif nonlinearity == "rectify":
            return lasagne.nonlinearities.rectify
        elif nonlinearity == "linear":
            return lasagne.nonlinearities.linear
        else:
            print "ERROR! Nonlinearity must be in [tanh, rectify, linear]"

    def check_match_dtype(self, nparrays):
        dtype = nparrays[0].dtype
        for array in nparrays:
            if array.dtype != dtype:
                print "ERROR! Arrays with mismatching data type"
                print nparrays

    def check_match_ndim(self, nparrays):
        ndim = len(nparrays[0].shape)
        for array in nparrays:
            if len(array.shape) != ndim:
                print "ERROR! Arrays with mismatching ndim"
                print nparrays

    def check_match_dims(self, nparrays):
        dims = nparrays[0].shape
        for array in nparrays:
            for i in range(len(dims)):
                if i == 0: continue
                else:
                    if array.shape[i] != dims[i]:
                        print "ERROR! Arrays with mismatching dimension >= 1"
                        print nparrays

    def check_ndarray_dim(self, matrix, dim):
        if matrix.ndim != dim:
            print "ERROR! ndarray is not of dim {}".format(dim)
            print matrix

    def get_theano_var(self, dtype, ndim, var_name=None):
        if dtype == "float64":
            if ndim == 1:
                var = T.dvector(var_name)
            elif ndim == 2:
                var = T.dmatrix(var_name)
            elif ndim == 3:
                var = T.dtensor3(var_name)
            elif ndim == 4:
                var = T.dtensor4(var_name)
            else:
                print "Unsupported theano var dimension {0}".format(ndim)
                raise ValueError
        elif dtype == "int64":
            if ndim == 1:
                var = T.lvector(var_name)
            elif ndim == 2:
                var = T.lmatrix(var_name)
            elif ndim == 3:
                var = T.ltensor3(var_name)
            elif ndim == 4:
                var = T.ltensor4(var_name)
            else:
                print "Unsupported theano var dimension {0}".format(ndim)
                raise ValueError
        else:
            print "Unsupported theano var dtype {0}".format(dtype)
            raise ValueError
        return var
  
    def get_tr_pred_var(self, network):
        return lasagne.layers.get_output(network)

    def get_te_pred_var(self, network):
        return lasagne.layers.get_output(network, deterministic=True)

    def get_tr_loss_var(self, tr_pred_var, target_var, loss_fn):
        if loss_fn == "L2":
            return lasagne.objectives.squared_error(
                tr_pred_var, target_var).mean()
        elif loss_fn == "L1":
            return abs(tr_pred_var - target_var).mean()

    def get_te_loss_var(self, te_pred_var, target_var, loss_fn):
        if loss_fn == "L2":
            return lasagne.objectives.squared_error(
                te_pred_var, target_var).mean()
        elif loss_fn == "L1":
            return abs(te_pred_var - target_var).mean()

    def get_all_updates(self, network, tr_loss_var, update_method, 
                        learning_rate, momentum=0.9):
        params = lasagne.layers.get_all_params(network, trainable=True)
        if update_method == "sgd":
            updates = lasagne.updates.sgd(
                tr_loss_var, params, learning_rate=learning_rate)
        elif update_method == "momentum":
            updates = lasagne.updates.momentum(
                tr_loss_var, params, learning_rate=learning_rate,
                momentum=momentum)
        elif update_method == "nesterov":
            updates = lasagne.updates.nesterov_momentum(
                tr_loss_var, params, learning_rate=learning_rate, 
                momentum=momentum)
        else:
            print "Error: Method must be in [sgd, momentum, nesterov]"
            raise ValueError
        return updates     

    def iterate_minibatches(self, inputs, targets, batchsize, shuffle=False):
        #print inputs.shape
        #print targets.shape
        assert len(inputs) == len(targets)
        if shuffle:
            indices = np.arange(len(inputs))
            np.random.shuffle(indices)
        for start_idx in range(0, len(inputs) - batchsize + 1, batchsize):
            if shuffle:
                excerpt = indices[start_idx:start_idx + batchsize]
            else:
                excerpt = slice(start_idx, start_idx + batchsize)
            yield inputs[excerpt], targets[excerpt]

    def make_input_layer(self, input_dims, input_var):
        input_layer = lasagne.layers.InputLayer(
            shape=tuple(input_dims), input_var=input_var)
        return input_layer

    def make_dense_layers(self, input_layer, widths, nonlin_fn, 
                          input_drop_rate=0, hidden_drop_rate=0):
        print "Making dense mlp"
        network = input_layer
        if input_drop_rate:
            network = lasagne.layers.dropout(network, p=input_drop_rate)
        for width in widths:
            network = lasagne.layers.DenseLayer(
                network, width, nonlinearity=nonlin_fn)
            if hidden_drop_rate:
                network = lasagne.layers.dropout(network, p=hidden_drop_rate)
        return network

    def make_regression_output_layer(self, network, num_outputs):
        linear = lasagne.nonlinearities.linear
        rectify = lasagne.nonlinearities.rectify
        if not self.nonnegative:
            print "regular model"
            network = lasagne.layers.DenseLayer(
                network, num_outputs, nonlinearity=linear)
        if self.nonnegative:
            print "nonnegative model"
            network = lasagne.layers.DenseLayer(
                network, num_outputs, nonlinearity=rectify)
        return network

    def run_epoch(self):
        #Run full epoch of training, and compute test error
        start = time.time()
        train_err = self.train_epoch(
            self.X_tr, self.y_tr, self.batch_size, self.train_fn)
        test_err = self.test_epoch(
            self.X_te, self.y_te, self.batch_size, self.test_fn)
        end = time.time()
        total = end - start
        self.num_epochs += 1
        #Store epoch data
        self.y_tr_hat = self.pred_fn(self.X_tr)
        self.y_te_hat = self.pred_fn(self.X_te)
        self.train_err_by_epoch.append(train_err)
        self.test_err_by_epoch.append(test_err)
        self.weights = [param.get_value() for param in self.params]
        self.pickle_epoch(
            self.weights, self.train_err_by_epoch, self.test_err_by_epoch, 
            self.y_tr_hat, self.y_te_hat, self.num_epochs)
        # Then we print the results for this epoch:
        print "Epoch {} took {:.3f}s".format(
            self.num_epochs, total)
        print "  training loss:\t\t{:.6f}".format(train_err)
        print "  test loss:\t\t{:.6f}".format(test_err)

    def run_epochs(self, num_epochs):
        while self.num_epochs < num_epochs:
            self.run_epoch()

    def split_matrix_by_col(self, matrix, sizes):
        start_idx = 0
        splits = []
        for size in sizes:
            splits.append(matrix[:,start_idx:start_idx+size])
            start_idx += size
        if start_idx != matrix.shape[1]:
            print "ERROR! Sum of sizes {0} ".format(start_idx) +\
                  "not equal to width of matrix {1}".format(matrix.shape[1])
            raise ValueError
        return splits

    def pickle_epoch(
            self, weights, tr_costs, te_costs, y_tr_hat, y_te_hat, epoch):
        self.make_epoch_dir(self.out_dir, self.name, epoch)
        epoch_dir = "{0}/{1}/epoch{2}".format(self.out_dir, self.name, epoch)
        self.pickle_object(weights, epoch_dir + "/weights.pkl")
        self.pickle_object(tr_costs, epoch_dir + "/tr_cost_by_epoch.pkl")
        self.pickle_object(te_costs, epoch_dir + "/te_cost_by_epoch.pkl")
        self.pickle_object(y_tr_hat, epoch_dir + "/y_tr_hat.pkl")
        self.pickle_object(y_te_hat, epoch_dir + "/y_te_hat.pkl")

    def unpickle_epoch(self, epoch):
        epoch_dir = self.out_dir + "/{0}/epoch{1}".format(self.name, epoch)
        train_errs = self.unpickle_object(epoch_dir + "/tr_cost_by_epoch.pkl")
        test_errs = self.unpickle_object(epoch_dir + "/te_cost_by_epoch.pkl")
        weights = self.unpickle_object(epoch_dir + "/weights.pkl")
        y_tr_hat = self.unpickle_object(epoch_dir + "/y_tr_hat.pkl")
        y_te_hat = self.unpickle_object(epoch_dir + "/y_te_hat.pkl")

        self.train_err_by_epoch = train_errs
        self.test_err_by_epoch = test_errs
        self.weights = weights
        self.y_tr_hat = y_tr_hat
        self.y_te_hat = y_te_hat
        self.num_epochs = epoch
      
        for idx, param in enumerate(self.params):
            param.set_value(weights[idx])

    def make_out_dir(self, out_dir):
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

    def make_nn_dir(self, out_dir, name):
        if not os.path.exists(out_dir):
            print "ERROR: out_dir {0} does not exist".format(out_dir)
            raise NameError
        dir_name = out_dir + "/" + name
        epoch_dirs = glob.glob(dir_name + "/" + "epoch*")
        if os.path.exists(dir_name):
            if len(epoch_dirs) > 0:
                print "BAD ERROR: dir with name {0} already exists".format(
                    dir_name)
                raise NameError
            else: return
        else:
            os.makedirs(dir_name)
  
    def make_init_data_dir(self, out_dir, name):
        dir_name = "{0}/{1}/init_data".format(out_dir, name)
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

    def make_epoch_dir(self, out_dir, name, epoch):
        dir_name = "{0}/{1}/epoch{2}".format(out_dir, name, epoch)
        if os.path.exists(dir_name):
            print "BAD ERROR: dir with name {0} already exists".format(dir_name)
            raise NameError
        else:
            os.makedirs(dir_name)

    def pickle_object(self, obj, pkl_fname):
        with open(pkl_fname, "w") as f:
            pickle.dump(obj, f)

    def unpickle_object(self, pkl_fname):
        with open(pkl_fname, "r") as f:
            obj = pickle.load(f)
        return obj

    def get_network(self, *args, **kwargs):
        raise NotImplementedError

    def train_fn(self, *args, **kwargs):
        raise NotImplementedError

    def test_fn(self, *args, **kwargs):
        raise NotImplementedError

    def pred_fn(self, *args, **kwargs):
        raise NotImplementedError

    def train_epoch(self, *args, **kwargs):
        raise NotImplementedError

    def test_epoch(self, *args, **kwargs):
        raise NotImplementedError

class FeedforwardMLP(RegressionMLP):    
    def __init__(
            self, X_tr, y_tr, X_te, y_te, widths=[200], nonlinearity="tanh", 
            X_val=False, y_val=False, input_drop_rate=0, hidden_drop_rate=0, 
            num_outputs=1, update_method="sgd", learning_rate=0.01, 
            loss_fn = "L2", momentum=0.9, batch_size=500, name="my_nn", 
            out_dir="lasagne_nn", reloaded=False, nonnegative=False):
       
        super(FeedforwardMLP, self).__init__(
            X_tr, y_tr, X_te, y_te, X_val, y_val, name=name, out_dir=out_dir, 
            reloaded=reloaded, nonnegative=nonnegative)
        #Store network dimensions
        self.widths = widths
        self.num_hidden_layers = len(widths)
        #Create theano variables 
        self.input_var = self.get_theano_var(
            self.input_dtype, self.input_ndim, var_name='inputs')
        self.output_var = self.get_theano_var(
            self.output_dtype, self.output_ndim, var_name='outputs')
        #Store NN parameters
        self.nonlinearity = nonlinearity
        self.input_drop_rate = input_drop_rate
        self.hidden_drop_rate = hidden_drop_rate
        self.update_method = update_method
        self.learning_rate = learning_rate
        self.momentum = momentum
        self.batch_size = batch_size
        self.loss_fn = loss_fn
        #Make lasagne objects
        self.nonlin_fn = self.get_nonlin_fn(self.nonlinearity)
        self.network = self.get_network(
            self.input_dims, self.input_var, self.widths, self.nonlin_fn, 
            self.input_drop_rate, self.hidden_drop_rate, self.num_outputs) 
        self.params = lasagne.layers.get_all_params(self.network)
        #print self.params
        self.weights = [param.get_value() for param in self.params]
        #for w in self.weights:
        #    print w.shape
        self.tr_pred_var = self.get_tr_pred_var(self.network)
        self.te_pred_var = self.get_te_pred_var(self.network)
        self.tr_loss_var = self.get_tr_loss_var(
            self.tr_pred_var, self.output_var, self.loss_fn)
        self.te_loss_var = self.get_te_loss_var(
            self.te_pred_var, self.output_var, self.loss_fn)
        self.updates_var = self.get_all_updates(
            self.network, self.tr_loss_var, self.update_method, 
            self.learning_rate, self.momentum)
        self.train_fn = self.get_train_fn(
            self.input_var, self.output_var, self.tr_loss_var, 
            self.updates_var)
        self.test_fn = self.get_test_fn(
            self.input_var, self.output_var, self.te_loss_var)
        self.pred_fn = self.get_pred_fn(self.input_var, self.te_pred_var)
        #Initialize batch training and test errors
        #try:
        train_err = self.test_epoch(self.X_tr, self.y_tr, 
                                     self.batch_size, self.test_fn)
        #except AssertionError:
        #    print self.X_tr
        #    print self.y_tr
        #    print self.X_tr.shape
        #    print self.y_tr.shape
        test_err = self.test_epoch(self.X_te, self.y_te, 
                                   self.batch_size, self.test_fn)
        self.train_err_by_epoch.append(train_err)
        self.test_err_by_epoch.append(test_err)
        self.y_tr_hat = self.pred_fn(self.X_tr)
        self.y_te_hat = self.pred_fn(self.X_te)
        #Store initial epoch
        if not reloaded:
            proc.pickle_obj(
                self.y_tr, "{0}/{1}/init_data/y_tr.pkl".format(out_dir, name))
            proc.pickle_obj(
                self.y_te, "{0}/{1}/init_data/y_te.pkl".format(out_dir, name))
            self.pickle_epoch(
                self.weights, self.train_err_by_epoch, self.test_err_by_epoch, 
                self.y_tr_hat, self.y_te_hat, self.num_epochs)

    def get_network(self, input_dims, input_var, widths, nonlin_fn, 
                    input_drop_rate, hidden_drop_rate, num_outputs):
        input_layer = self.make_input_layer(self.input_dims, self.input_var)
        network = self.make_dense_layers(
            input_layer, self.widths, self.nonlin_fn, self.input_drop_rate, 
            self.hidden_drop_rate)
        network = self.make_regression_output_layer(network, num_outputs)
        return network
       
    def change_training_vars(self, update_method=None, learning_rate=None, 
                             momentum=None):
        if update_method is not None:
            self.update_method = update_method
        if learning_rate is not None:
            self.learning_rate = learning_rate
        if momentum is not None:
            self.momentum = momentum
        self.updates_var = self.get_all_updates(
            self.network, self.tr_loss_var, self.update_method,
            self.learning_rate, self.momentum)
        self.train_fn = self.get_train_fn(
            self.input_var, self.output_var, self.tr_loss_var, self.updates_var)
       
    def get_train_fn(self, input_var, output_var, loss, updates):
        train_fn = theano.function([input_var, output_var], loss, 
                                   updates=updates)
        return train_fn
    
    def get_test_fn(self, input_var, output_var, te_loss_var):
    # Compile a second function computing the validation loss and accuracy:
        test_fn = theano.function([input_var, output_var], te_loss_var)
        return test_fn

    def get_pred_fn(self, input_var, te_pred_var):
        pred_fn = theano.function([input_var], te_pred_var)
        return pred_fn

    def train_epoch(self, X_tr, y_tr, batch_size, train_fn):
        train_err = 0
        train_batches = 0
        for batch in self.iterate_minibatches(
                X_tr, y_tr, batch_size, shuffle=True):
            inputs, targets = batch
            train_err += train_fn(*[inputs, targets])
            train_batches += 1
        return train_err / train_batches

    def test_epoch(self, X_te, y_te, batch_size, test_fn):
        test_err = 0
        test_batches = 0
        for batch in self.iterate_minibatches(
                X_te, y_te, batch_size, shuffle=False):
            inputs, targets = batch
            err = test_fn(*[inputs, targets])
            test_err += err
            test_batches += 1
        return test_err / test_batches

class SplitMLP(RegressionMLP):    
    def __init__(
            self, X_tr, y_tr, X_te, y_te, rel_cod_idxs, rel_nt_idxs, 
            cod_widths=[200], nt_widths=[50], nonlinearity="tanh", 
            X_val=False, y_val=False, input_drop_rate=0, hidden_drop_rate=0, 
            num_outputs=1, update_method="sgd", learning_rate=0.01, 
            momentum=0.9, batch_size=500, name="my_nn", out_dir="lasagne_nn", 
            reloaded=False):
        super(SplitMLP, self).__init__(
            X_tr, y_tr, X_te, y_te, X_val, y_val, name=name, out_dir=out_dir, 
            reloaded=reloaded)
        #Store network dimensions
        self.rel_cod_idxs = rel_cod_idxs
        self.rel_nt_idxs = rel_nt_idxs
        self.cod_widths = cod_widths
        self.nt_widths = nt_widths
        self.split_sizes = [64 * len(self.rel_cod_idxs), 4 * len(self.rel_nt_idxs)]
        self.X_tr_split = self.split_matrix_by_col(self.X_tr, self.split_sizes)
        self.X_val_split = self.split_matrix_by_col(self.X_val, self.split_sizes)
        self.X_te_split = self.split_matrix_by_col(self.X_te, self.split_sizes)
        #Create theano variables 
        self.cod_input_var = self.get_theano_var(
            self.input_dtype, self.input_ndim, var_name='cod_inputs')
        self.nt_input_var = self.get_theano_var(
            self.input_dtype, self.input_ndim, var_name='nt_inputs')
        self.input_vars = [self.cod_input_var, self.nt_input_var]
        self.output_var = self.get_theano_var(
            self.output_dtype, self.output_ndim, var_name='outputs')
        #Store NN parameters
        self.nonlinearity = nonlinearity
        self.input_drop_rate = input_drop_rate
        self.hidden_drop_rate = hidden_drop_rate
        self.update_method = update_method
        self.learning_rate = learning_rate
        self.momentum = momentum
        self.batch_size = batch_size
        #Make lasagne objects
        self.nonlin_fn = self.get_nonlin_fn(self.nonlinearity)
        self.network = self.get_network(
            self.split_sizes, self.input_vars, self.cod_widths, self.nt_widths, 
            self.nonlin_fn, self.input_drop_rate, self.hidden_drop_rate, 
            self.num_outputs) 
        self.params = lasagne.layers.get_all_params(self.network)
        self.weights = [param.get_value() for param in self.params]
        self.tr_pred_var = self.get_tr_pred_var(self.network)
        self.te_pred_var = self.get_te_pred_var(self.network)
        self.tr_loss_var = self.get_tr_loss_var(
            self.tr_pred_var, self.output_var)
        self.te_loss_var = self.get_tr_loss_var(
            self.te_pred_var, self.output_var)
        self.updates_var = self.get_all_updates(
            self.network, self.tr_loss_var, self.update_method, 
            self.learning_rate, self.momentum)
        self.train_fn = self.get_train_fn(
            self.input_vars, self.output_var, self.tr_loss_var, 
            self.updates_var)
        self.test_fn = self.get_test_fn(
            self.input_vars, self.output_var, self.te_loss_var)
        self.pred_fn = self.get_pred_fn(self.input_vars, self.te_pred_var)
        #Initialize batch training and test errors
        train_err = self.test_epoch(self.X_tr, self.y_tr, 
                                     self.batch_size, self.test_fn)
        test_err = self.test_epoch(self.X_te, self.y_te, 
                                   self.batch_size, self.test_fn)
        self.train_err_by_epoch.append(train_err)
        self.test_err_by_epoch.append(test_err)
        self.y_tr_hat = self.pred_fn(*self.X_tr_split)
        self.y_te_hat = self.pred_fn(*self.X_te_split)
        #Store initial epoch
        if not reloaded:
            proc.pickle_obj(
                self.y_tr, "{0}/{1}/init_data/y_tr.pkl".format(out_dir, name))
            proc.pickle_obj(
                self.y_te, "{0}/{1}/init_data/y_te.pkl".format(out_dir, name))
            self.pickle_epoch(
                self.weights, self.train_err_by_epoch, self.test_err_by_epoch, 
                self.y_tr_hat, self.y_te_hat, self.num_epochs)

    def get_network(
            self, split_sizes, input_vars, cod_widths, nt_widths, nonlin_fn, 
            input_drop_rate, hidden_drop_rate, num_outputs):
        cod_input_layer = self.make_input_layer(
            (None, split_sizes[0]), self.input_vars[0])
        nt_input_layer = self.make_input_layer(
            (None, split_sizes[1]), self.input_vars[1])
        cod_network = self.make_dense_layers(
            cod_input_layer, self.cod_widths, self.nonlin_fn, self.input_drop_rate, 
            self.hidden_drop_rate)
        nt_network = self.make_dense_layers(
            nt_input_layer, self.nt_widths, self.nonlin_fn, self.input_drop_rate, 
            self.hidden_drop_rate)
        cod_output = self.make_regression_output_layer(cod_network, 1)
        nt_output = self.make_regression_output_layer(nt_network, 1)
        merge_outputs = lasagne.layers.ConcatLayer([cod_output, nt_output]) 
        network = self.make_regression_output_layer(merge_outputs, num_outputs)
        return network
       
    def change_training_vars(self, update_method=None, learning_rate=None, 
                             momentum=None):
        if update_method is not None:
            self.update_method = update_method
        if learning_rate is not None:
            self.learning_rate = learning_rate
        if momentum is not None:
            self.momentum = momentum
        self.updates_var = self.get_all_updates(
            self.network, self.tr_loss_var, self.update_method,
            self.learning_rate, self.momentum)
        self.train_fn = self.get_train_fn(
            self.input_var, self.output_var, self.tr_loss_var, self.updates_var)
       
    def get_train_fn(self, input_vars, output_var, loss, updates):
        all_vars = input_vars + [output_var]
        train_fn = theano.function(all_vars, loss, updates=updates)
        return train_fn
    
    def get_test_fn(self, input_vars, output_var, te_loss_var):
    # Compile a second function computing the validation loss and accuracy:
        all_vars = input_vars + [output_var]
        test_fn = theano.function(all_vars, te_loss_var)
        return test_fn

    def get_pred_fn(self, input_vars, te_pred_var):
        pred_fn = theano.function(input_vars, te_pred_var)
        return pred_fn

    def train_epoch(self, X_tr, y_tr, batch_size, train_fn):
        train_err = 0
        train_batches = 0
        for batch in self.iterate_minibatches(
                X_tr, y_tr, batch_size, shuffle=True):
            inputs, targets = batch
            inputs = self.split_matrix_by_col(inputs, self.split_sizes)
            all_inputs = inputs + [targets]
            train_err += train_fn(*all_inputs)
            train_batches += 1
        return train_err / train_batches

    def test_epoch(self, X_te, y_te, batch_size, test_fn):
        test_err = 0
        test_batches = 0
        for batch in self.iterate_minibatches(
                X_te, y_te, batch_size, shuffle=False):
            inputs, targets = batch
            inputs = self.split_matrix_by_col(inputs, self.split_sizes)
            all_inputs = inputs + [targets]
            err = test_fn(*all_inputs)
            test_err += err
            test_batches += 1
        return test_err / test_batches

    def run_epoch(self):
        #Run full epoch of training, and compute test error
        start = time.time()
        train_err = self.train_epoch(
            self.X_tr, self.y_tr, self.batch_size, self.train_fn)
        test_err = self.test_epoch(
            self.X_te, self.y_te, self.batch_size, self.test_fn)
        end = time.time()
        total = end - start
        self.num_epochs += 1
        #Store epoch data
        self.y_tr_hat = self.pred_fn(*self.X_tr_split)
        self.y_te_hat = self.pred_fn(*self.X_te_split)
        self.train_err_by_epoch.append(train_err)
        self.test_err_by_epoch.append(test_err)
        self.weights = [param.get_value() for param in self.params]
        self.pickle_epoch(
            self.weights, self.train_err_by_epoch, self.test_err_by_epoch, 
            self.y_tr_hat, self.y_te_hat, self.num_epochs)
        # Then we print the results for this epoch:
        print "Epoch {} took {:.3f}s".format(
            self.num_epochs, total)
        print "  training loss:\t\t{:.6f}".format(train_err)
        print "  test loss:\t\t{:.6f}".format(test_err)

class AdjacencyMLP(RegressionMLP):    
    def __init__(
            self, X_tr, y_tr, X_te, y_te, rel_cod_idxs=[], cod_adj_idxs=[], 
            rel_nt_idxs=[], nt_adj_idxs=[], widths=[200], nonlinearity="tanh", 
            X_val=False, y_val=False, input_drop_rate=0, hidden_drop_rate=0, 
            num_outputs=1, update_method="sgd", learning_rate=0.01, 
            momentum=0.9, batch_size=500, name="my_nn", out_dir="lasagne_nn", 
            reloaded=False, loss_fn = "L2"):
        super(AdjacencyMLP, self).__init__(
            X_tr, y_tr, X_te, y_te, X_val, y_val, name=name, out_dir=out_dir, 
            reloaded=reloaded)
        #Store network dimensions
        #NOTE: Make rules for cod/nt adj idxs
        self.rel_cod_idxs = rel_cod_idxs
        self.cod_adj_idxs = cod_adj_idxs
        self.rel_nt_idxs = rel_nt_idxs
        self.nt_adj_idxs = nt_adj_idxs
        self.split_sizes = [64 for elt in self.rel_cod_idxs] + \
                           [4 for elt in self.rel_nt_idxs]
        self.X_tr_split = self.split_matrix_by_col(self.X_tr, self.split_sizes)
        self.X_val_split = self.split_matrix_by_col(
            self.X_val, self.split_sizes)
        self.X_te_split = self.split_matrix_by_col(self.X_te, self.split_sizes)
        self.widths = widths
        self.num_hidden_layers = len(widths)
        #Create theano variables 
        self.cod_input_vars = [
            self.get_theano_var(
                self.input_dtype, self.input_ndim, 
                var_name='input_cod_' + str(cod)) 
            for cod in rel_cod_idxs]
        self.nt_input_vars = [
            self.get_theano_var(
                self.input_dtype, self.input_ndim, 
                var_name='input_nt_' + str(nt)) 
            for nt in rel_nt_idxs]
        self.input_vars = self.cod_input_vars + self.nt_input_vars
        self.output_var = self.get_theano_var(
            self.output_dtype, self.output_ndim, var_name='outputs')
        #Store NN parameters
        self.nonlinearity = nonlinearity
        self.input_drop_rate = input_drop_rate
        self.hidden_drop_rate = hidden_drop_rate
        self.update_method = update_method
        self.learning_rate = learning_rate
        self.momentum = momentum
        self.batch_size = batch_size
        self.loss_fn = loss_fn
        #Make lasagne objects
        self.nonlin_fn = self.get_nonlin_fn(self.nonlinearity)
        self.network = self.get_network(
            self.cod_input_vars, self.rel_cod_idxs, self.cod_adj_idxs, 
            self.nt_input_vars, self.rel_nt_idxs, self.nt_adj_idxs, 
            self.widths, self.nonlin_fn, self.input_drop_rate, 
            self.hidden_drop_rate, self.num_outputs) 
        self.params = lasagne.layers.get_all_params(self.network)
        self.weights = [param.get_value() for param in self.params]
        self.tr_pred_var = self.get_tr_pred_var(self.network)
        self.te_pred_var = self.get_te_pred_var(self.network)
        self.tr_loss_var = self.get_tr_loss_var(
            self.tr_pred_var, self.output_var, self.loss_fn)
        self.te_loss_var = self.get_tr_loss_var(
            self.te_pred_var, self.output_var, self.loss_fn)
        self.updates_var = self.get_all_updates(
            self.network, self.tr_loss_var, self.update_method, 
            self.learning_rate, self.momentum)
        self.train_fn = self.get_train_fn(
            self.input_vars, self.output_var, self.tr_loss_var, 
            self.updates_var)
        self.test_fn = self.get_test_fn(
            self.input_vars, self.output_var, self.te_loss_var)
        self.pred_fn = self.get_pred_fn(
            self.input_vars, self.te_pred_var, self.split_sizes)
        #Initialize batch training and test errors
        train_err = self.test_epoch(
            self.X_tr, self.y_tr, self.batch_size, self.test_fn)
        test_err = self.test_epoch(
            self.X_te, self.y_te, self.batch_size, self.test_fn)
        self.train_err_by_epoch.append(train_err)
        self.test_err_by_epoch.append(test_err)
        self.y_tr_hat = self.pred_fn(self.X_tr)
        self.y_te_hat = self.pred_fn(self.X_te)
        #Store initial epoch
        if not reloaded:
            proc.pickle_obj(
                self.y_tr, "{0}/{1}/init_data/y_tr.pkl".format(out_dir, name))
            proc.pickle_obj(
                self.y_te, "{0}/{1}/init_data/y_te.pkl".format(out_dir, name))
            self.pickle_epoch(
                self.weights, self.train_err_by_epoch, self.test_err_by_epoch, 
                self.y_tr_hat, self.y_te_hat, self.num_epochs)

    def get_network(
            self, cod_input_vars, rel_cod_idxs, cod_adj_idxs, nt_input_vars, 
            rel_nt_idxs, nt_adj_idxs, widths, nonlin_fn, input_drop_rate, 
            hidden_drop_rate, num_outputs):
        print "execute 1"
        cod_input_layers = [
            self.make_input_layer((None, 64), cod_input_var) 
            for cod_input_var in cod_input_vars] 
        nt_input_layers = [
            self.make_input_layer((None, 4), nt_input_var) 
            for nt_input_var in nt_input_vars] 
        input_layers = []
        for idx, layer in zip(rel_cod_idxs, cod_input_layers):
            if idx not in cod_adj_idxs:
                print "execute 2"
                input_layers.append(layer)
            elif idx in cod_adj_idxs and idx - 1 not in rel_cod_idxs:
                print "execute 3"
                input_layers.append(layer)
            else:
                print "execute 4"
                prev_layer = cod_input_layers[rel_cod_idxs.index(idx)]
                adj_layer = lasagne.layers.AdjacencyLayer(prev_layer, layer)
                print lasagne.layers.get_all_params(adj_layer)
                input_layers.append(layer)
        for idx, layer in zip(rel_nt_idxs, nt_input_layers):
            if idx not in nt_adj_idxs:
                input_layers.append(layer)
            elif idx in nt_adj_idxs and idx - 1 not in rel_nt_idxs:
                input_layers.append(layer)
            else:
                prev_layer = nt_input_layers[rel_nt_idxs.index(idx)]
                adj_layer = lasagne.layers.AdjacencyLayer(prev_layer, layer)
                input_layers.append(layer)
        print "execute"
        merge_layer = lasagne.layers.ConcatLayer(input_layers)
        network = self.make_dense_layers(
            merge_layer, self.widths, self.nonlin_fn, self.input_drop_rate, 
            self.hidden_drop_rate)
        network = self.make_regression_output_layer(network, num_outputs)
        return network
       
    def change_training_vars(
            self, update_method=None, learning_rate=None, momentum=None):
        if update_method is not None:
            self.update_method = update_method
        if learning_rate is not None:
            self.learning_rate = learning_rate
        if momentum is not None:
            self.momentum = momentum
        self.updates_var = self.get_all_updates(
            self.network, self.tr_loss_var, self.update_method,
            self.learning_rate, self.momentum)
        self.train_fn = self.get_train_fn(
            self.cod_input_vars + self.nt_input_vars, self.output_var, 
            self.tr_loss_var, self.updates_var)
       
    def get_train_fn(self, input_vars, output_var, loss, updates):
        all_vars = input_vars + [output_var]
        train_fn = theano.function(all_vars, loss, updates=updates)
        return train_fn
    
    def get_test_fn(self, input_vars, output_var, te_loss_var):
    # Compile a second function computing the validation loss and accuracy:
        all_vars = input_vars + [output_var]
        test_fn = theano.function(all_vars, te_loss_var)
        return test_fn

    def get_pred_fn(self, input_vars, te_pred_var, split_sizes):
        theano_pred_fn = theano.function(input_vars, te_pred_var)
        def pred_fn(X):
            X_split = self.split_matrix_by_col(X, split_sizes)
            return theano_pred_fn(*X_split)
        return pred_fn

    def train_epoch(self, X_tr, y_tr, batch_size, train_fn):
        train_err = 0
        train_batches = 0
        for batch in self.iterate_minibatches(
                X_tr, y_tr, batch_size, shuffle=True):
            inputs, targets = batch
            inputs = self.split_matrix_by_col(inputs, self.split_sizes)
            all_inputs = inputs + [targets]
            train_err += train_fn(*all_inputs)
            train_batches += 1
        return train_err / train_batches

    def test_epoch(self, X_te, y_te, batch_size, test_fn):
        split_sizes = [64 for elt in self.rel_cod_idxs] +\
                      [4 for elt in self.rel_nt_idxs]
        test_err = 0
        test_batches = 0
        for batch in self.iterate_minibatches(
                X_te, y_te, batch_size, shuffle=False):
            inputs, targets = batch
            inputs = self.split_matrix_by_col(inputs, self.split_sizes)
            all_inputs = inputs + [targets]
            err = test_fn(*all_inputs)
            test_err += err
            test_batches += 1
        return test_err / test_batches

    def run_epoch(self):
        #Run full epoch of training, and compute test error
        start = time.time()
        train_err = self.train_epoch(
            self.X_tr, self.y_tr, self.batch_size, self.train_fn)
        test_err = self.test_epoch(
            self.X_te, self.y_te, self.batch_size, self.test_fn)
        end = time.time()
        total = end - start
        self.num_epochs += 1
        #Store epoch data
        self.y_tr_hat = self.pred_fn(self.X_tr)
        self.y_te_hat = self.pred_fn(self.X_te)
        self.train_err_by_epoch.append(train_err)
        self.test_err_by_epoch.append(test_err)
        self.weights = [param.get_value() for param in self.params]
        self.pickle_epoch(
            self.weights, self.train_err_by_epoch, self.test_err_by_epoch, 
            self.y_tr_hat, self.y_te_hat, self.num_epochs)
        # Then we print the results for this epoch:
        print "Epoch {} took {:.3f}s".format(
            self.num_epochs, total)
        print "  training loss:\t\t{:.6f}".format(train_err)
        print "  test loss:\t\t{:.6f}".format(test_err)

def load_dataset(
        project_dir, expt_dir, sam_fname, gene_seq_fname, gene_len_fname, 
        tr_codons_fname, te_codons_fname, outputs_fname, rel_cod_idxs, 
        rel_nt_idxs):

    tr_codon_bounds = proc.load_codon_set_bounds(tr_codons_fname)
    te_codon_bounds = proc.load_codon_set_bounds(te_codons_fname)
    tr_codon_set = proc.expand_codon_set(tr_codon_bounds)
    te_codon_set = proc.expand_codon_set(te_codon_bounds)

    len_dict = proc.get_len_dict(gene_len_fname)
    cds_dict = proc.get_cds_dict(gene_seq_fname, len_dict)
    rel_struc_idxs=False
    struc_dict=False
    X_tr = proc.get_X(tr_codon_set, cds_dict, rel_cod_idxs=rel_cod_idxs,
             rel_nt_idxs=rel_nt_idxs, rel_struc_idxs=rel_struc_idxs,
             struc_dict=struc_dict)
    X_te = proc.get_X(te_codon_set, cds_dict, rel_cod_idxs=rel_cod_idxs,
             rel_nt_idxs=rel_nt_idxs, rel_struc_idxs=rel_struc_idxs,
             struc_dict=struc_dict)

    outputs = proc.load_outputs(outputs_fname)
    y_tr = proc.get_y(tr_codon_set, outputs)
    y_te = proc.get_y(te_codon_set, outputs)

    X_train = X_tr.transpose(1, 0)
    X_val = X_te.transpose(1, 0)
    X_test = X_val
    y_train = y_tr.transpose(1, 0)
    y_val = y_te.transpose(1, 0)
    y_test = y_val
    return X_train, y_train, X_val, y_val, X_test, y_test

def load_weinberg_dataset():    
    project_dir = "/mnt/lareaulab/rtunney/Regression"
    expt_dir = project_dir + "/expts/weinberg"
    sam_fname = expt_dir + "/process/weinberg.transcript.sam"
    gene_seq_fname = project_dir + "/genome_data/yeast_13cds10.fa"
    gene_len_fname = project_dir + "/genome_data/yeast_13cds10_lengths.txt"
    tr_codons_fname = expt_dir + "/process/tr_set_bounds.trunc.20.20.min_cts.200.min_cod.100.top.500.txt"
    te_codons_fname = expt_dir + "/process/te_set_bounds.trunc.20.20.min_cts.200.min_cod.100.top.500.txt"
    outputs_fname = expt_dir + "/process/outputs.txt"
    rel_cod_idxs = range(-7, 6)
    rel_nt_idxs = range(-20, 18)
    X_train, y_train, X_val, y_val, X_test, y_test = load_dataset(
        project_dir, expt_dir, sam_fname, gene_seq_fname, gene_len_fname,
        tr_codons_fname, te_codons_fname, outputs_fname, rel_cod_idxs,
        rel_nt_idxs)
    return X_train, y_train, X_val, y_val, X_test, y_test

def load_human_unt_large_dataset():
    project_dir = "/mnt/lareaulab/rtunney/Regression"
    expt_dir = project_dir + "/expts/human_unt_large"
    sam_fname = expt_dir + "/process/human_unt.transcript.sam"
    gene_seq_fname = project_dir + "/genome_data/gencode.v22.transcript.13cds10.fa"
    gene_len_fname = project_dir + "/genome_data/gencode.v22.transcript.13cds10.lengths.txt"
    tr_codons_fname = expt_dir + "/process/tr_set_bounds.trunc.20.20.min_cts.200.min_cod.100.top.300.txt"
    te_codons_fname = expt_dir + "/process/te_set_bounds.trunc.20.20.min_cts.200.min_cod.100.top.300.txt"
    outputs_fname = expt_dir + "/process/outputs.txt"
    rel_cod_idxs = range(-7, 6)
    rel_nt_idxs = range(-20, 18)
    X_train, y_train, X_val, y_val, X_test, y_test = load_dataset(
        project_dir, expt_dir, sam_fname, gene_seq_fname, gene_len_fname,
        tr_codons_fname, te_codons_fname, outputs_fname, rel_cod_idxs,
        rel_nt_idxs)
    return X_train, y_train, X_val, y_val, X_test, y_test

def load_human_unt_small_dataset():
    project_dir = "/mnt/lareaulab/rtunney/Regression"
    expt_dir = project_dir + "/expts/human_unt_small"
    sam_fname = expt_dir + "/process/human_unt.transcript.sam"
    gene_seq_fname = project_dir + "/genome_data/gencode.v22.transcript.13cds10.fa"
    gene_len_fname = project_dir + "/genome_data/gencode.v22.transcript.13cds10.lengths.txt"
    tr_codons_fname = expt_dir + "/process/tr_set_bounds.trunc.20.20.min_cts.200.min_cod.100.top.300.txt"
    te_codons_fname = expt_dir + "/process/te_set_bounds.trunc.20.20.min_cts.200.min_cod.100.top.300.txt"
    outputs_fname = expt_dir + "/process/outputs.txt"
    rel_cod_idxs = range(-7, 6)
    rel_nt_idxs = range(-20, 18)
    X_train, y_train, X_val, y_val, X_test, y_test = load_dataset(
        project_dir, expt_dir, sam_fname, gene_seq_fname, gene_len_fname,
        tr_codons_fname, te_codons_fname, outputs_fname, rel_cod_idxs,
        rel_nt_idxs)
    return X_train, y_train, X_val, y_val, X_test, y_test

def load_stanford12_large_dataset():
    project_dir = "/mnt/lareaulab/rtunney/Regression"
    expt_dir = project_dir + "/expts/stanford12_large"
    sam_fname = expt_dir + "/process/yeast12.transcript.sam"
    gene_seq_fname = project_dir + "/genome_data/gencode.v22.transcript.13cds10.fa"
    gene_len_fname = project_dir + "/genome_data/gencode.v22.transcript.13cds10.lengths.txt"
    tr_codons_fname = expt_dir + "/process/tr_set_bounds.trunc.20.20.min_cts.200.min_cod.100.top.400.txt"
    te_codons_fname = expt_dir + "/process/te_set_bounds.trunc.20.20.min_cts.200.min_cod.100.top.400.txt"
    outputs_fname = expt_dir + "/process/outputs.txt"
    rel_cod_idxs = range(-7, 6)
    rel_nt_idxs = range(-20, 18)
    X_train, y_train, X_val, y_val, X_test, y_test = load_dataset(
        project_dir, expt_dir, sam_fname, gene_seq_fname, gene_len_fname,
        tr_codons_fname, te_codons_fname, outputs_fname, rel_cod_idxs,
        rel_nt_idxs)
    return X_train, y_train, X_val, y_val, X_test, y_test

def load_stanford12_small_dataset():
    project_dir = "/mnt/lareaulab/rtunney/Regression"
    expt_dir = project_dir + "/expts/stanford12_small"
    sam_fname = expt_dir + "/process/yeast12.transcript.sam"
    gene_seq_fname = project_dir + "/genome_data/gencode.v22.transcript.13cds10.fa"
    gene_len_fname = project_dir + "/genome_data/gencode.v22.transcript.13cds10.lengths.txt"
    tr_codons_fname = expt_dir + "/process/tr_set_bounds.trunc.20.20.min_cts.200.min_cod.100.top.300.txt"
    te_codons_fname = expt_dir + "/process/te_set_bounds.trunc.20.20.min_cts.200.min_cod.100.top.300.txt"
    outputs_fname = expt_dir + "/process/outputs.txt"
    rel_cod_idxs = range(-7, 6)
    rel_nt_idxs = range(-20, 18)
    X_train, y_train, X_val, y_val, X_test, y_test = load_dataset(
        project_dir, expt_dir, sam_fname, gene_seq_fname, gene_len_fname,
        tr_codons_fname, te_codons_fname, outputs_fname, rel_cod_idxs,
        rel_nt_idxs)
    return X_train, y_train, X_val, y_val, X_test, y_test

def run_all_models():
    X_tr, y_tr, X_val, y_val, X_te, y_te = load_weinberg_dataset()
    my_nn = FeedforwardMLP(X_tr, y_tr, X_te, y_te, learning_rate=0.01, update_method="nesterov", widths=[200])
    my_nn.run_epochs(10)
    print my_nn.y_te.transpose()
    print my_nn.y_te_hat.transpose()
    params = lasagne.layers.get_all_params(my_nn.network, trainable=True)
    wts_dir = "/mnt/lareaulab/rtunney/Regression/expts/weinberg/lasagne_wts/L1W200"
    proc.pickle_obj(params[0].get_value(), wts_dir + "/W_0.pkl")
    proc.pickle_obj(params[1].get_value(), wts_dir + "/b_0.pkl")
    proc.pickle_obj(params[2].get_value(), wts_dir + "/W_1.pkl")
    proc.pickle_obj(params[3].get_value(), wts_dir + "/b_1.pkl")

    X_tr, y_tr, X_val, y_val, X_te, y_te = load_human_unt_large_dataset()
    my_nn1 = FeedforwardMLP(X_tr, y_tr, X_te, y_te, learning_rate=0.01, update_method="nesterov", widths=[200])
    my_nn1.run_epochs(10)
    print my_nn1.y_te.transpose()
    print my_nn1.y_te_hat.transpose()
    params = lasagne.layers.get_all_params(my_nn1.network, trainable=True)
    wts_dir = "/mnt/lareaulab/rtunney/Regression/expts/human_unt_large/lasagne_wts/L1W200"
    proc.pickle_obj(params[0].get_value(), wts_dir + "/W_0.pkl")
    proc.pickle_obj(params[1].get_value(), wts_dir + "/b_0.pkl")
    proc.pickle_obj(params[2].get_value(), wts_dir + "/W_1.pkl")
    proc.pickle_obj(params[3].get_value(), wts_dir + "/b_1.pkl")

    X_tr, y_tr, X_val, y_val, X_te, y_te = load_human_unt_small_dataset()
    my_nn2 = lasagnenn.FeedforwardMLP(X_tr, y_tr, X_te, y_te, learning_rate=0.01, update_method="nesterov", widths=[200])
    my_nn2.run_epochs(5)
    print my_nn2.y_te.transpose()
    print my_nn2.y_te_hat.transpose()
    params = lasagne.layers.get_all_params(my_nn2.network, trainable=True)
    wts_dir = "/mnt/lareaulab/rtunney/Regression/expts/human_unt_small/lasagne_wts/L1W200"
    proc.pickle_obj(params[0].get_value(), wts_dir + "/W_0.pkl")
    proc.pickle_obj(params[1].get_value(), wts_dir + "/b_0.pkl")
    proc.pickle_obj(params[2].get_value(), wts_dir + "/W_1.pkl")
    proc.pickle_obj(params[3].get_value(), wts_dir + "/b_1.pkl")

    X_tr, y_tr, X_val, y_val, X_te, y_te = lasagnenn.load_human_unt_large_dataset()
    my_nn3 = lasagnenn.FeedforwardMLP(X_tr, y_tr, X_te, y_te, learning_rate=0.01, update_method="nesterov", widths=[200])
    my_nn3.run_epochs(10)
    print my_nn3.y_te.transpose()
    print my_nn3.y_te_hat.transpose()
    params = lasagne.layers.get_all_params(my_nn3.network, trainable=True)
    wts_dir = "/mnt/lareaulab/rtunney/Regression/expts/stanford12_large/lasagne_wts/L1W200"
    proc.pickle_obj(params[0].get_value(), wts_dir + "/W_0.pkl")
    proc.pickle_obj(params[1].get_value(), wts_dir + "/b_0.pkl")
    proc.pickle_obj(params[2].get_value(), wts_dir + "/W_1.pkl")
    proc.pickle_obj(params[3].get_value(), wts_dir + "/b_1.pkl")

    X_tr, y_tr, X_val, y_val, X_te, y_te = lasagnenn.load_human_unt_small_dataset()
    my_nn4 = lasagnenn.FeedforwardMLP(X_tr, y_tr, X_te, y_te, learning_rate=0.01, update_method="nesterov", widths=[200])
    my_nn4.run_epochs(5)
    print my_nn4.y_te.transpose()
    print my_nn4.y_te_hat.transpose()
    params = lasagne.layers.get_all_params(my_nn4.network, trainable=True)
    wts_dir = "/mnt/lareaulab/rtunney/Regression/expts/stanford12_small/lasagne_wts/L1W200"
    proc.pickle_obj(params[0].get_value(), wts_dir + "/W_0.pkl")
    proc.pickle_obj(params[1].get_value(), wts_dir + "/b_0.pkl")
    proc.pickle_obj(params[2].get_value(), wts_dir + "/W_1.pkl")
    proc.pickle_obj(params[3].get_value(), wts_dir + "/b_1.pkl")

if __name__ == "__main__":
    X_tr, y_tr, X_val, y_val, X_te, y_te = load_weinberg_dataset()
    rel_cod_idxs = range(-7, 6)
    cod_adj_idxs = range(-3, 3)
    rel_nt_idxs = range(-20, 18)
    nt_adj_idxs = range(-20, -10) + range(8, 18)
    #print X_tr.shape
    out_dir = "/mnt/lareaulab/rtunney/Regression/expts/weinberg/lasagne_nn"
    name = "file_struc_test1"
    my_nn = FeedforwardMLP(
        X_tr, y_tr, X_te, y_te, learning_rate=0.01, update_method="nesterov", 
        widths=[200], out_dir=out_dir, name=name)
#    cod_adj_idxs = range(-2, 3)
#    rel_nt_idxs = []
#    nt_adj_idxs = [] 
#    my_nn = AdjacencyMLP(X_tr, y_tr, X_te, y_te, rel_cod_idxs=rel_cod_idxs, rel_nt_idxs=rel_nt_idxs, cod_adj_idxs=cod_adj_idxs, nt_adj_idxs=nt_adj_idxs, learning_rate=0.01, update_method="nesterov", widths=[200], nonlinearity="tanh", name=name, out_dir=out_dir)
#    my_nn = SplitMLP(X_tr, y_tr, X_te, y_te, rel_cod_idxs=rel_cod_idxs, rel_nt_idxs=rel_nt_idxs, learning_rate=0.03, update_method="sgd", cod_widths=[200], nt_widths=[100], nonlinearity="tanh", name=name, out_dir=out_dir)
    my_nn.run_epochs(30)
#    print my_nn.y_te.transpose()
#    print my_nn.y_te_hat.transpose()
