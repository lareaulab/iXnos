import os
import math
import copy
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
import pickle

alpha="ACGT"
nts = ["A", "C", "G", "T"]
codons = [x+y+z for x in alpha for y in alpha for z in alpha]
cod2id = {codon:idx for idx, codon in enumerate(codons)}
id2cod = {idx:codon for codon, idx in cod2id.items()}
nt2id = {nt:idx for idx, nt in enumerate(alpha)}
id2nt = {idx:nt for nt, idx in nt2id.items()}

def get_lasagne_codon_featimp1(nn, rel_cod_idxs, rel_nt_idxs=False):
    print "Running old version of feat imps!"
    # 2D list for feature importance matrix (positions x codons)
    feat_imps = []
    num_cod_feats = len(rel_cod_idxs) * 64
    y_te_hat = nn.pred_fn(nn.X_te)
    X_te_copy = copy.deepcopy(nn.X_te)
    N_te = nn.X_te.shape[0]
    #Iterate over codon positions
    for cod_idx in rel_cod_idxs:
        print "running position " + str(cod_idx)
        #Feature importance list for that positions
        cod_feat_imps = []
        #Set codons in position to their expected value in X_te_copy
         #Save idxs of altered features in data matrices
        alt_feat_idxs = [] 
         #Get indexes of codon features at given cod_idx in data matrices
        cod_feat_idx = rel_cod_idxs.index(cod_idx)
        pos_cod_idxs = range((cod_feat_idx) * 64, (cod_feat_idx + 1) * 64)
        alt_feat_idxs.extend(pos_cod_idxs)
         #Get distribution of codons at cod_idx
        pos_cod_cts = nn.X_te[:,pos_cod_idxs]
        pos_cod_dist = pos_cod_cts.mean(axis=0)
        X_te_copy[:,pos_cod_idxs] = np.zeros((N_te, 64)) + pos_cod_dist
        #Get nt_idxs in rel_nt_idxs that overlap with cod_idx
        nt_idxs = []
        if rel_nt_idxs:
            for nt_idx in rel_nt_idxs:
                if cod_idx * 3 <= nt_idx < (cod_idx + 1) * 3:
                    nt_idxs.append(nt_idx)
        #Set overlapping nts to their expected value in X_te_copy
        for nt_idx in nt_idxs:
            nt_feat_idx = rel_nt_idxs.index(nt_idx)
            pos_nt_idxs = range(num_cod_feats + (nt_feat_idx * 4), 
                           num_cod_feats + ((nt_feat_idx + 1) * 4))
            alt_feat_idxs.extend(pos_nt_idxs)
            pos_nt_cts = nn.X_te[:,pos_nt_idxs]
            pos_nt_dist = pos_nt_cts.mean(axis=0)
            X_te_copy[:,pos_nt_idxs] = np.zeros((N_te, 4)) + pos_nt_dist
        y_te_tilde = nn.pred_fn(X_te_copy)
        y_te_delta = y_te_hat - y_te_tilde
        #Iterate over codons
        for cod_id in range(64):
            #Get index of that codon in nn weight matrix
            cod_var_idx = rel_cod_idxs.index(cod_idx) * 64 + cod_id
            pt_bools = nn.X_te[:,cod_var_idx]
            num_pts = pt_bools.sum()
            mean_imp = np.dot(pt_bools, y_te_delta) / float(num_pts)
            cod_feat_imps.append(mean_imp[0])
        #Return altered features to original state
        X_te_copy[:,alt_feat_idxs] = nn.X_te[:,alt_feat_idxs]
        feat_imps.append(cod_feat_imps)
    return np.array(feat_imps).transpose()

def get_lasagne_codon_featimp2(nn, rel_cod_idxs, rel_nt_idxs=False):
    print "Running new version of feat imps!"
    # 2D list for feature importance matrix (positions x codons)
    feat_imps = []
    num_cod_feats = len(rel_cod_idxs) * 64
    y_te_hat = nn.pred_fn(nn.X_te)
    X_te_copy = copy.deepcopy(nn.X_te)
    N_te = nn.X_te.shape[0]
    #Iterate over codon positions
    for cod_idx in rel_cod_idxs:
        print "running position " + str(cod_idx)
        #Feature importance list for that positions
        cod_feat_imps = []
        #Set codons in position to their expected value in X_te_copy
         #Save idxs of altered features in data matrices
        alt_feat_idxs = [] 
         #Get indexes of codon features at given cod_idx in data matrices
        cod_feat_idx = rel_cod_idxs.index(cod_idx)
        pos_cod_idxs = range((cod_feat_idx) * 64, (cod_feat_idx + 1) * 64)
        alt_feat_idxs.extend(pos_cod_idxs)
         #Get distribution of codons at cod_idx
        pos_cod_cts = nn.X_te[:,pos_cod_idxs]
        pos_cod_dist = pos_cod_cts.mean(axis=0)
        X_te_copy[:,pos_cod_idxs] = np.zeros((N_te, 64))
        #Get nt_idxs in rel_nt_idxs that overlap with cod_idx
        nt_idxs = []
        if rel_nt_idxs:
            for nt_idx in rel_nt_idxs:
                if cod_idx * 3 <= nt_idx < (cod_idx + 1) * 3:
                    nt_idxs.append(nt_idx)
        #Set overlapping nts to their expected value in X_te_copy
        pos_nts_idxs = []
        for nt_idx in nt_idxs:
            nt_feat_idx = rel_nt_idxs.index(nt_idx)
            pos_nt_idxs = range(num_cod_feats + (nt_feat_idx * 4), 
                           num_cod_feats + ((nt_feat_idx + 1) * 4))
            pos_nts_idxs.extend(pos_nt_idxs)
            alt_feat_idxs.extend(pos_nt_idxs)
            pos_nt_cts = nn.X_te[:,pos_nt_idxs]
            pos_nt_dist = pos_nt_cts.mean(axis=0)
            X_te_copy[:,pos_nt_idxs] = np.zeros((N_te, 4))
        #y_te_tilde = nn.pred_fn(X_te_copy)
        #y_te_delta = y_te_hat - y_te_tilde
        #Iterate over codons
        y_te_tildes = []
        for cod_id in range(64):
            pos_cod_idx = cod_feat_idx * 64 + cod_id
            X_te_copy[:,pos_cod_idx] = 1
            for nt_idx in nt_idxs: 
                frame = nt_idx % 3
                cod = id2cod[cod_id]
                nt = cod[frame]
                nt_id = nt2id[nt]
                pos_nt_idx = num_cod_feats + (nt_feat_idx * 4) + nt_id
                X_te_copy[:,pos_nt_idx] = 1
            y_te_tilde = nn.pred_fn(X_te_copy)
            y_te_tildes.append(y_te_tilde)
            X_te_copy[:,pos_cod_idxs] = np.zeros((N_te, 64))
            X_te_copy[:,pos_nts_idxs] = np.zeros((N_te, len(pos_nts_idxs)))
        y_te_tilde_mat = np.hstack(y_te_tildes)
        E_y_te_tilde = np.dot(y_te_tilde_mat, pos_cod_dist).reshape(N_te, 1)
        y_te_delta = y_te_hat - E_y_te_tilde
        for cod_id in range(64):
            #Get index of that codon in nn weight matrix
            cod_var_idx = rel_cod_idxs.index(cod_idx) * 64 + cod_id
            pt_bools = nn.X_te[:,cod_var_idx]
            num_pts = pt_bools.sum()
            mean_imp = np.dot(pt_bools, y_te_delta) / float(num_pts)
            cod_feat_imps.append(mean_imp[0])
        #Return altered features to original state
        X_te_copy[:,alt_feat_idxs] = nn.X_te[:,alt_feat_idxs]
        feat_imps.append(cod_feat_imps)
    return np.array(feat_imps).transpose()

def plot_cod_inst_feat_imp_scatters(fi_mat1, fi_mat2, xlabel, ylabel, rel_cod_idxs, parent_dir, dir_name):
    num_cods = fi_mat1.shape[1]
    if fi_mat2.shape[1] != num_cods: 
        print "Error! Feature importance matrices should have same width"
    if len(rel_cod_idxs) != num_cods:
        print "Error! Feature importance matrices don't match size of rel_cod_idxs"
    if not os.path.exists(parent_dir):
        print "Error! Parent dir {0} does not exist".format(parent_dir)
    if not os.path.exists("{0}/{1}".format(parent_dir, dir_name)):
        os.makedirs("{0}/{1}".format(parent_dir, dir_name))
    for idx, cod_idx in enumerate(rel_cod_idxs):
        fig = plt.scatter(fi_mat1[:,idx], fi_mat2[:,idx])
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title("Feature importance comparison, codon {0}".format(cod_idx))
        xlims = fig.axes.get_xlim() 
        ylims = fig.axes.get_ylim()
        print xlims, ylims

        #lims = [min(xlims[0], ylims[0]), max(xlims[0], ylims[0])]
        #fig.axes.set_xlim(lims) 
        #fig.axes.set_ylim(lims)
        text_x = xlims[0] + (xlims[1] - xlims[0]) * 0.1     
        text_y = ylims[0] + (ylims[1] - ylims[0]) * 0.9
        not_stop_idxs = np.logical_not(
            np.isnan(np.sum(np.hstack([fi_mat1, fi_mat2]), axis=1)))
        pear_r, p = pearsonr(fi_mat1[:,idx][[not_stop_idxs]].flatten(), 
                             fi_mat2[:,idx][[not_stop_idxs]].flatten())
        plt.text(text_x, text_y, "r = {0}".format(round(pear_r, 4))) 
        out_fname = "{0}/{1}/featimp_cod{2}_scatter.pdf".format(parent_dir, dir_name, cod_idx)  
        plt.savefig(out_fname)
        plt.close()

def plot_2d_fi_colormap(fi_mat, rel_cod_idxs, title, out_fname):
    # Make a 2d heatmap for feature importance scores
    # This has been edited but not tested
    alpha = "AGCT"
    alpha_idxs = [0, 2, 1, 3]
    yticks = np.array(range(0, 52, 4) + [50, 53, 57]) 
    yticklabels = [x + y for x in alpha for y in alpha]
    cod_sort_idxs = [x * 16 + y * 4 + z for x in alpha_idxs 
        for y in alpha_idxs for z in alpha_idxs]
    nonstop_idxs = range(64)
    nonstop_idxs.remove(48)
    nonstop_idxs.remove(49)
    nonstop_idxs.remove(52)
    fi_mat = fi_mat[cod_sort_idxs]
    fi_mat = fi_mat[nonstop_idxs]
    fig = plt.imshow(fi_mat, interpolation="nearest")
    ax = fig.get_axes()
    ax.set_xticks(range(0, len(rel_cod_idxs), 3))
    ax.set_xticklabels(rel_cod_idxs[::3])
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)
    plt.xlabel("Position")
    plt.ylabel("Codon")
    plt.colorbar()
    plt.title(title)
    plt.savefig(out_fname)
    plt.close()
    
def plot_weight_matrix_heatmaps(
        weights, rel_cod_idxs, title_prefix, outfile_prefix, 
        rel_nt_idxs=[], rel_struc_idxs=[], cmap_name="coolwarm"):
    # This has been edited but not tested
    # Make a set of 2d heatmaps to visualize NN weight matrices
    num_cods = len(rel_cod_idxs)
    num_cod_feats = num_cods * 64
    num_nts = len(rel_nt_idxs)
    num_nt_feats = num_nts * 4
    num_struc_feats = len(rel_struc_idxs)
    # Store parameters for first layer plots
    alpha = "ACGT"
    yticks = np.array(range(0, 52, 4) + [50, 54, 57]) 
    yticklabels = [x + y for x in alpha for y in alpha]
    nonstop_idxs = range(64)
    nonstop_idxs.remove(48)
    nonstop_idxs.remove(50)
    nonstop_idxs.remove(56)
    # Inputs -> HL1 weight matrix and b vector
    W0 = weights[0]
    b0 = weights[1]
    b0 = b0.reshape(b0.size, 1)
    layer_size1 = W0.shape[1]
    # Make heatmap for bias vector
    bias_out_fname = outfile_prefix + ".bias.pdf"
    fig = plt.imshow(b0, interpolation="nearest", cmap=plt.get_cmap(cmap_name))
    ax = fig.get_axes()
    #ax.set_xticks(range(len(rel_struc_idxs)))
    ax.set_yticks(np.arange(b0.size))
    ax.set_yticklabels(np.arange(b0.size))
    #ax.set_yticks([])
    #ax.set_yticklabels([])
    plt.ylabel("Hidden Unit")
    plt.colorbar()
    plt.title(title_prefix + " bias weights")
    plt.savefig(bias_out_fname)
    plt.close()
    # Make heatmap for each hidden unit
    for unit_idx in xrange(layer_size1):
        cod_mat = W0[:num_cod_feats,unit_idx].reshape(num_cods, 64).transpose()
        start_row = num_cod_feats
        if num_nt_feats > 0:
            end_row = start_row + num_nt_feats
            nt_mat = W0[start_row:end_row,unit_idx].reshape(num_nts,4).transpose()
            start_row = end_row
        if num_struc_feats > 0:
            end_row = start_row + num_struc_feats
            struc_mat = W0[start_row:end_row,unit_idx].transpose().reshape(num_struc_feats, 1)
        # Make codon heatmap
        cod_mat = cod_mat[nonstop_idxs]
        cod_out_fname = outfile_prefix + ".cod.HU{0}.pdf".format(unit_idx)
        fig = plt.imshow(cod_mat, interpolation="nearest", cmap=plt.get_cmap(cmap_name))
        ax = fig.get_axes()
        ax.set_xticks(range(0, len(rel_cod_idxs), 3))
        ax.set_xticklabels(rel_cod_idxs[::3])
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)
        plt.xlabel("Position")
        plt.ylabel("Codon")
        plt.colorbar()
        plt.title(title_prefix + " codon feats, HU{0}".format(unit_idx))
        plt.savefig(cod_out_fname)
        plt.close()
        # Make nt heatmap
        if num_nt_feats > 0:
            nt_out_fname = outfile_prefix +\
                ".nt.HU{0}.pdf".format(unit_idx)
            fig = plt.imshow(nt_mat, interpolation="nearest", cmap=plt.get_cmap(cmap_name))
            ax = fig.get_axes()
            ax.set_xticks(range(0, len(rel_nt_idxs), 3))
            ax.set_xticklabels(rel_nt_idxs[::3])
            ax.set_yticks(np.arange(4))
            ax.set_yticklabels(['A', 'C', 'G', 'T'])
            plt.xlabel("Position")
            plt.ylabel("Nucleotide")
            plt.colorbar()
            plt.title(title_prefix + " nt feats, HU{0}".format(unit_idx))
            plt.savefig(nt_out_fname)
            plt.close()
        # Make struc heatmap
        if num_struc_feats > 0:
            struc_out_fname = outfile_prefix +\
                ".struc.HU{0}.pdf".format(unit_idx)
            fig = plt.imshow(struc_mat, interpolation="nearest", cmap=plt.get_cmap(cmap_name))
            ax = fig.get_axes()
            ax.set_yticks(range(len(rel_struc_idxs)))
            ax.set_yticklabels(rel_struc_idxs)
            #ax.set_yticks([])
            #ax.set_yticklabels([])
            plt.ylabel("Position")
            plt.colorbar()
            plt.title(title_prefix +\
                " struc feats, HU{0}".format(unit_idx))
            plt.savefig(struc_out_fname)
            plt.close()
    
def get_test_err(epoch_dir):
    test_errs_by_epoch = pickle.load(open(
        epoch_dir + "/te_cost_by_epoch.pkl", "r"))
    test_err = test_errs_by_epoch[-1]
    return test_err

def get_series_test_errs(epoch_dirs):
    test_errs = []
    for epoch_dir in epoch_dirs: 
        test_err = get_test_err(epoch_dir)
        test_errs.append(test_err)
    return test_errs

def get_test_corrs(nn_dir, epoch, method="pearson"):
    if method not in ["pearson", "spearman"]:
        print "Error! Correlation method must be 'pearson' or 'spearman'"
    y_te_fname = nn_dir + "/init_data/y_te.pkl"
    y_te_hat_fname = nn_dir + "/epoch{0}/y_te_hat.pkl".format(epoch)
    y_te = pickle.load(open(y_te_fname, "r"))
    y_te_hat = pickle.load(open(y_te_hat_fname, "r"))
    if method == "pearson":
        return pearsonr(y_te, y_te_hat)[0]
    elif method == "spearman":
        return spearmanr(y_te, y_te_hat)[0]

def save_nocod_series_test_MSE_diffs(ref_epoch_dir, epoch_dirs, leaveout_idxs, out_fname):
    if len(epoch_dirs) != len(leaveout_idxs):
        print "ERROR! should input same # epoch_dirs and leaveout_idxs"
    ref_errs_by_epoch = pickle.load(open(
        ref_epoch_dir + "/te_cost_by_epoch.pkl", "r"))
    ref_err = ref_errs_by_epoch[-1]
    test_errs = get_series_test_errs(epoch_dirs)
    outfile = open(out_fname, "w")
    outfile.write("codon\ttest_MSE\tdifference\n")
    outfile.write("full\t{0}\t0\n".format(round(ref_err, 4)))
    for idx, cod_idx in enumerate(leaveout_idxs):
        outfile.write("{0}\t{1}\t{2}\n".format(
            cod_idx, round(test_errs[idx], 4), round(test_errs[idx] - ref_err, 4)))
    outfile.close()

def save_rep_series_MSEs(nn_parent_dir, series_names, epochs, num_reps, out_fname, round_digits=8):
    """
        series_names: list of base nn names (with _rep# appended, reps must be 
                     zero indexed and consecutive)
        epochs: list of epochs to pull data from for each nn name_model
        num_reps: number of reps for each name_model (1 int, same for all modls)
        out_fname: file to write to
    """
    # If only one epoch entered, use that as the epoch for all series
    if type(epochs) == int:
        epochs = [epochs for i in range(len(series_names))]
    # Open out file, write header
    outfile = open(out_fname, "w")
    outfile_colstring = "rep_series"
    for rep in range(num_reps):
        outfile_colstring += "\trep_{0}_MSE".format(rep)
    outfile_colstring += "\tmean_MSE\tstd_err_mean"
    outfile.write(outfile_colstring + "\n")
    # For each series
    for idx, name in enumerate(series_names):
        epoch_dirs = [nn_parent_dir + "/{0}_rep{1}/epoch{2}".format(
                      name, rep, epochs[idx]) 
                      for rep in range(num_reps)]
        test_errs = get_series_test_errs(epoch_dirs)
        mean_err = sum(test_errs) / len(test_errs)
        std_err = np.std(np.array(test_errs))/math.sqrt(len(test_errs))
        line = name 
        for rep in range(num_reps):
            line += "\t{0}".format(round(test_errs[rep], round_digits))
        line += "\t{0}".format(round(mean_err, round_digits))
        line += "\t{0}".format(round(std_err, round_digits))  
        outfile.write(line + "\n")
    outfile.close()


def save_rep_series_corrs(
        nn_parent_dir, series_names, epochs, num_reps, out_fname, 
        round_digits=8, method="pearson"):
    """
    Inputs: 
        nn_parent_dir (str): Directory above nn directories
        series_names (list of strs): base nn series names (with _rep# 
            appended, reps must be zero indexed and consecutive)
        epochs: list of epochs to pull data from for each nn model name
        num_reps (int): number of reps for each nn model name (1 int, same
            for all models)
        out_fname (str): output file name
    """
    if method not in ["pearson", "spearman"]:
        print "Error! Correlation method must be 'pearson' or 'spearman'"
    # If only one epoch entered, use that as the epoch for all series
    if type(epochs) == int:
        epochs = [epochs for i in range(len(series_names))]
    # Open out file, write header
    outfile = open(out_fname, "w")
    outfile_colstring = "rep_series"
    for rep in range(num_reps):
        outfile_colstring += "\trep_{0}_corr".format(rep)
    outfile.write(outfile_colstring + "\n")
    # For each series
    for idx, series_name in enumerate(series_names):
        epoch = epochs[idx]
        line = series_name
        for rep in range(num_reps):
            nn_dir = nn_parent_dir + "/{0}_rep{1}".format(series_name, rep)
            corr = get_test_corrs(
                nn_dir, epoch, method=method)
            line += "\t{0}".format(round(corr, round_digits))
        outfile.write(line + "\n")
    outfile.close()

def add_MSE_diff_column(in_fname, ref_mse, out_fname):
    in_file = open(in_fname, "r")
    out_file = open(out_fname, "w")
    colstring = in_file.readline().strip()
    colstring += "\tMSE diff\n"
    out_file.write(colstring)
    for line in in_file:
        line = line.strip().split()
        mean_err = float(line[-2])
        mse_diff = mean_err - ref_mse
        mse_diff = str(round(mse_diff, 8))
        line.append(mse_diff)
        line = "\t".join(line)
        out_file.write(line + "\n")
    out_file.close()
        

def add_rep_series_MSE_diff_stats( 
        nn_parent_dir, ref_series, epoch, num_reps, fname, round_digits=8):
    #get data for ref series
    ref_epoch_dirs = [nn_parent_dir +\
        "/{0}_rep{1}/epoch{2}".format(ref_series, rep, epoch) 
        for rep in range(num_reps)]
    ref_errs = get_series_test_errs(ref_epoch_dirs)
    #compute statistics for ref series
    ref_mean = sum(ref_errs) / len(ref_errs)
    ref_std_err = np.std(np.array(ref_errs)) / math.sqrt(len(ref_errs))
    #open temp. outfile
    outfile = open(fname, "r")
    tmp_fname = fname + ".tmp"
    tmp_outfile = open(tmp_fname, "w")
    #make header 
    header = outfile.readline().strip() + "\tMSE_diff\tstd_err_diff"
    tmp_outfile.write("#Reference series for differences: {0}\n".format(ref_series))
    tmp_outfile.write(header + "\n")
    #modify data for each series
    for line in outfile: 
        line = line.strip().split()
        mean_MSE = float(line[-2])
        std_err_mean = float(line[-1])
        diff_MSE = mean_MSE - ref_mean
        if line[0] == ref_series:
            std_err_diff = 0
        else: 
            std_err_diff = math.sqrt(std_err_mean**2 + ref_std_err**2)
        line.append(str(round(diff_MSE, round_digits)))
        line.append(str(round(std_err_diff, round_digits)))
        line = "\t".join(line) + "\n"
        tmp_outfile.write(line)
    tmp_outfile.close()
    #move tmp file to where old file was
    os.rename(tmp_fname, fname)
