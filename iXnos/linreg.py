import numpy as np
import scipy
import os
import pickle

def train(X_tr, y_tr):
    XXt = np.dot(X_tr, X_tr.transpose())
    n = XXt.shape[0]
    XXt_plus1 = XXt + np.identity(n)
    L = np.linalg.cholesky(XXt_plus1)
    wts = scipy.linalg.cho_solve((L, True), np.dot(X_tr, y_tr.transpose()))
    return wts

def predict(X, wts):
    return np.dot(X.transpose(), wts).reshape((1, -1))

def make_linreg_data_file(
        linreg_dir, gene_seq_fname, gene_len_fname, tr_codons_fname, 
        te_codons_fname, outputs_fname, rel_cod_idxs=False, 
        rel_nt_idxs=False):
    f = open("{0}/linreg_data.txt".format(linreg_dir), "w")
    f.write("Transcripts_fasta_file: {0}".format(gene_seq_fname) + "\n")
    f.write("Transcripts_length_file: {0}".format(gene_len_fname) + "\n")
    f.write("Training_codon_bounds: {0}".format(tr_codons_fname) + "\n")
    f.write("Test_codon_bounds: {0}".format(te_codons_fname) + "\n")
    f.write("Outputs: {0}".format(outputs_fname) + "\n")
    if rel_cod_idxs:
        f.write(
            "Relative_codon_indices: " +\
            "\t".join([str(elt) for elt in rel_cod_idxs])) + "\n")
    if rel_nt_idxs:
        f.write(
            "Relative_nt_indices: " +\
            "\t".join([str(elt) for elt in rel_nt_idxs])) + "\n")
    f.close()

def make_linreg_parent_dir(expt_dir):
    parent_dir = "{0}/linreg".format(expt_dir)
    if not os.path.exists(parent_dir):
        os.makedirs(parent_dir)

def make_linreg_dir(expt_dir, name):
    dir_name = "{0}/linreg/{1}".format(expt_dir, name)
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

def linreg_init(
        expt_dir, name, gene_seq_fname, gene_len_fname, tr_codons_fname, 
        te_codons_fname, outputs_fname, rel_cod_idxs=False, 
        rel_nt_idxs=False):
    make_linreg_parent_dir(expt_dir)
    make_linreg_dir(expt_dir, name)
    linreg_dir = "{0}/linreg/{1}".format(expt_dir, name)
    make_linreg_data_file(
        linreg_dir, gene_seq_fname, gene_len_fname, tr_codons_fname, 
        te_codons_fname, outputs_fname, rel_cod_idxs=rel_cod_idxs, 
        rel_nt_idxs=rel_nt_idxs)

def train_and_predict(X_tr, X_te, y_tr, y_te):
    wts = train(X_tr, y_tr)
    y_tr_hat = predict(X_tr, wts)
    y_te_hat = predict(X_te, wts)
    return wts, y_tr_hat, y_te_hat

def get_linreg_MSE(linreg_dir, test=True, train=False):
    assert test != train, "Must choose linreg MSE on training xor test set"
    if test: 
        y = pickle.load(open(linreg_dir + "/y_te.pkl", "r"))    
        y_hat = pickle.load(open(linreg_dir + "/y_te_hat.pkl", "r"))    
    if train: 
        y = pickle.load(open(linreg_dir + "/y_tr.pkl", "r"))    
        y_hat = pickle.load(open(linreg_dir + "/y_tr_hat.pkl", "r"))    
    return ((y.ravel() - y_hat.ravel())**2).mean()

def aggregate_linreg_MSEs(
        expt_linreg_dir, series_names, out_fname):
    # Open out file, write header
    f = open(out_fname, "w")
    header = "series\ttrain_MSE\ttest_MSE\n"
    f.write(header)
    # Populate out file
    for series_name in series_names:
        linreg_dir = expt_linreg_dir + "/" + series_name
        train_mse = get_linreg_MSE(linreg_dir, train=True, test=False)
        test_mse = get_linreg_MSE(linreg_dir, train=False, test=True)
        line = "{0}\t{1}\t{2}\n".format(series_name, train_mse, test_mse)
        f.write(line)
    f.close()

