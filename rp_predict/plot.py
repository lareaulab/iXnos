import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import collections as mc
import pylab as pl
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.stats import gaussian_kde
from scipy.interpolate import LSQUnivariateSpline
import pickle
import os
from sklearn.metrics import roc_curve, auc
from rp_predict.process import load_codon_set_bounds, expand_codon_set, \
    load_outputs, get_y, load_cts_by_codon, has_enough_cts, \
    sort_genes_by_density, get_y_data, get_gene2symbol_dict, \
    get_mat_idxs_by_gene, save_gene_list, sort_genes_by_density, \
    get_genes_by_density2

def make_roc_plot(y, y_hat, title, fname, y_cutoff=0):
    
    fpr, tpr, thresh = roc_curve(y, y_hat)
    roc_auc = auc(fpr, tpr)

    plt.figure()
    plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(title)
    plt.legend(loc="lower right")
    plt.savefig(fname)
    plt.close()
   
def make_cts_by_size_plot(cts_by_size_and_frame, title, fname, sizes=False):
    if not sizes:
        sizes = cts_by_size_and_frame.keys()
    sorted_sizes = sorted(sizes)
    tickpos = [size + 0.45 for size in sorted_sizes]
    cts = []
    #Possibly kludgy replacement for line below
    #for size in sorted_sizes:
    #    try: cts.append(sum(cts_by_size_and_frame[size]))
    #    except KeyError: cts.append(0)
    cts = [sum(cts_by_size_and_frame[size]) for size in sorted_sizes]
    plt.bar(sorted_sizes, cts)
    plt.xlabel("FP size")
    plt.ylabel("Counts")
    plt.title(title)
    plt.xticks(tickpos, sorted_sizes)
    plt.savefig(fname)
    plt.close()

def make_frame_by_size_plot(size, data, title, fname, cts=False, freq=False):
    if cts == freq:
        print "ERROR: must choose cts or freq for y axis units"
        raise AssertionError
    frames = [0, 1, 2]
    tickpos = [frame + 0.45 for frame in frames]
    colors = ["red", "green", "blue"]
    plt.bar(frames, data, color=colors)
    plt.xticks(tickpos, frames)
    plt.xlabel("Frame")
    if cts:
        ylab = "Counts"
    elif freq:
        ylab = "Frequency"
    plt.ylabel(ylab)
    plt.title(title)
    plt.savefig(fname)
    plt.close()

def get_te_spearmanr_by_gene(
        y_te, y_te_hat, te_genes, outputs, cod_trunc_5p, cod_trunc_3p):
    start_idx = 0
    spearmanr_vals = []
    p_vals = []
    for gene in te_genes:
        num_cods = len(outputs[gene]) - cod_trunc_5p - cod_trunc_3p
        end_idx = start_idx + num_cods
        y_te_hat_gene = y_te_hat.flatten()[start_idx:end_idx]
        y_te_gene = y_te.flatten()[start_idx:end_idx]
        spearmanr_val, p_val = spearmanr(y_te_hat_gene, y_te_gene)
        spearmanr_vals.append(spearmanr_val)
        p_vals.append(p_val)
    return spearmanr_vals, p_vals

def make_density_scatter_plot(
        y, y_hat, title, xlab, ylab, xlims, ylims, textpos, fname, bins=0, 
        no_title=False):
    spear_r, spear_p = spearmanr(y.flatten(), y_hat.flatten())
    pear_r, pear_p = pearsonr(y.flatten(), y_hat.flatten())
    spear_r = round(spear_r, 4)
    pear_r = round(pear_r, 4)
    num_pts = y.size
    mse = float(((y - y_hat)**2).sum()) / num_pts
    mse = round(mse, 4)
    xy = np.vstack([y,y_hat])
    z = gaussian_kde(xy)(xy)
    idx = z.argsort()
    y, y_hat, z = y[idx], y_hat[idx], z[idx]
    fig, ax = plt.subplots()
    plt.axhline(y=0, linewidth=1, color = '#d3d3d3', zorder=1)
    plt.axvline(x=0, linewidth=1, color = '#d3d3d3', zorder=2)
    plt.scatter(y, y_hat, c=z, s=10, edgecolor='', zorder=3)
    if not no_title:
        plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.xlim(xlims)
    plt.ylim(ylims)
    textstring = "Pear. {0} Spear. {1}\nMSE {2}".format(pear_r, spear_r, mse)
    #PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
    plt.colorbar() 
    if bins > 0:
        sort_idxs = np.argsort(y.flatten())
        y_sorted = y.flatten()[sort_idxs]
        y_hat_sorted = y_hat.flatten()[sort_idxs]
        bin_bounds = []
        bin_spear_rs = []
        bin_pear_rs = []
        bin_mses = []
        bin_size = num_pts / bins
        for i in range(bins):
            start_idx = i * bin_size
            end_idx = (i + 1) * bin_size
            bin_y = y_sorted[start_idx:end_idx]
            bin_y_hat = y_hat_sorted[start_idx:end_idx]
            bin_spear_r, p = spearmanr(bin_y, bin_y_hat) 
            bin_pear_r, p = spearmanr(bin_y, bin_y_hat)
            bin_spear_r = round(bin_spear_r, 4)
            bin_pear_r = round(bin_pear_r, 4)
            bin_mse = float(((bin_y - bin_y_hat)**2).sum()) / bin_size
            bin_mse = round(bin_mse, 4)
            bin_spear_rs.append(bin_spear_r)
            bin_pear_rs.append(bin_pear_r)
            bin_mses.append(bin_mse)
            start_bound = y_sorted[start_idx]
            bin_bounds.append(start_bound)
            if i != 0: 
                plt.plot([start_bound, start_bound], [ylims[0], ylims[0] + 1], 
                         color="red")
            textstring += "\nBin {0}: Pear. {1} MSE {2}".format(
                i, bin_pear_r, bin_mse)    
        textpos = (textpos[0], textpos[1] - 0.4 * bins)
    plt.text(textpos[0], textpos[1], textstring)
    plt.savefig(fname)
    plt.close()

def make_scatter_plot_with_xticks(
        xticks, y, title, xlab, ylab, fname, tickskip=0):
    num_xticks = len(xticks)
    #x = range(num_xticks)
    fig, ax = plt.subplots()
    ax.scatter(xticks, y)
    if tickskip:
        skip_ticks = range(min(xticks), max(xticks) + 2, tickskip+1)
        plt.xticks(skip_ticks, skip_ticks)
    else:
        plt.xticks(xticks, xticks)
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.savefig(fname)
    plt.close()

def make_smoothed_err_plot(y, y_hat, title, fname, order=3, num_knots=10):
    sort_idxs = np.argsort(y.flatten())
    y_sorted = y.flatten()[sort_idxs]
    y_hat_sorted = y_hat.flatten()[sort_idxs]
    sq_errs = (y_sorted - y_hat_sorted)**2
    cum_sq_errs = np.cumsum(sq_errs)

    winrad = 100
    avg_errs = np.array([sq_errs[i-winrad:i+winrad].sum()/(2 * winrad) 
                        for i in range(winrad, y.size - winrad)])
    log_avg_errs = np.log10(avg_errs)

    fig, ax = plt.subplots()
    x = np.arange(avg_errs.size)
    ax.plot(x, log_avg_errs, color="red")
    knots = np.linspace(x.min() + 1, x.max() - 1, num_knots)
    spline = LSQUnivariateSpline(x, log_avg_errs, k=order, t=knots)
    log_avg_errs_fit = spline(x)
    
    ax.plot(x, log_avg_errs_fit, "b-")
    plt.title(title)
    plt.xlabel("Points ordered by true outputs")
    plt.ylabel("Log_10 of squared error")
    plt.savefig(fname)
    plt.close()

def make_binned_err_plot(y, y_hat, title, fname, num_bins, no_title=False):
    sort_idxs = np.argsort(y.flatten())
    y_sorted = y.flatten()[sort_idxs]
    y_hat_sorted = y_hat.flatten()[sort_idxs]
    sq_errs = (y_sorted - y_hat_sorted)**2
    num_pts = y.size
    bin_size = num_pts/num_bins
    avg_errs = np.array([sq_errs[i*bin_size : (i+1)*bin_size].sum()/(bin_size)
                        for i in range(num_bins)])
    log_avg_errs = np.log10(avg_errs)

    fig, ax = plt.subplots()
    x = np.arange(avg_errs.size) #* bin_size
    #ax.plot(x, log_avg_errs, "b.")
    ax.plot(x, avg_errs, "b.")
    plt.xlabel("Bin no. (codons ranked by true scaled counts)")
    plt.ylabel("Bin MSE")
    ax.set_yscale('log')
    ax2 = ax.twiny()
    x2tick_location = ax.xaxis.get_ticklocs() #Get the tick locations in data coordinates as a numpy array
    ax2.set_xticks(x2tick_location)
    x2ticklabels = np.array([round(y_sorted[int(loc * bin_size)], 2) for loc in x2tick_location if loc < y.size])
    #x2ticklabels = np.array([round(tick, 2) for tick in y_sorted[x2tick_idxs]])
    ax2.set_xticklabels(x2ticklabels)
    #Remove plot title, change to upper x axis label
    #plt.title(title, y = 1.06)
    plt.title("True scaled counts", y = 1.06)
    plt.legend("topleft")
    plt.savefig(fname)
    plt.close()

def make_binned_err_diff_plot(
        y1, y1_hat, y2, y2_hat, title, fname, num_bins, no_title=False):
    sort_idxs1 = np.argsort(y1.flatten())
    y1_sorted = y1.flatten()[sort_idxs1]
    y1_hat_sorted = y1_hat.flatten()[sort_idxs1]
    sq_errs1 = (y1_sorted - y1_hat_sorted)**2
    sort_idxs2 = np.argsort(y2.flatten())
    y2_sorted = y2.flatten()[sort_idxs2]
    y2_hat_sorted = y2_hat.flatten()[sort_idxs2]
    sq_errs2 = (y2_sorted - y2_hat_sorted)**2
    num_pts = y1.size
    bin_size = num_pts/num_bins
    avg_errs1 = np.array([sq_errs1[i*bin_size : (i+1)*bin_size].sum()/(bin_size)
                        for i in range(num_bins)])
    avg_errs2 = np.array([sq_errs2[i*bin_size : (i+1)*bin_size].sum()/(bin_size)
                        for i in range(num_bins)])
    diff_errs = avg_errs1 - avg_errs2
    fig, ax = plt.subplots()
    x = np.arange(diff_errs.size) #* bin_size
    ax.plot(x, diff_errs, "b.")
    plt.axhline(y=0, color='r')
    plt.xlabel("Bin no. (codons ranked by true scaled counts)")
    plt.ylabel("MSE Difference")
    ax2 = ax.twiny()
    x2tick_location = ax.xaxis.get_ticklocs() #Get the tick locations in data coordinates as a numpy array
    ax2.set_xticks(x2tick_location)
    x2ticklabels = np.array([round(y1_sorted[int(loc * bin_size)], 2) for loc in x2tick_location if loc < y1.size])
    #x2ticklabels = np.array([round(tick, 2) for tick in y_sorted[x2tick_idxs]])
    ax2.set_xticklabels(x2ticklabels)
    #Remove plot title, change to upper x axis label
    #plt.title(title, y = 1.06)
    plt.title("True scaled counts", y = 1.06)
    plt.legend("topleft")
    plt.savefig(fname)
    plt.close()

    
def make_cum_err_plot(y, y_hat, title, fname):
    sort_idxs = np.argsort(y.flatten())
    y_sorted = y.flatten()[sort_idxs]
    y_hat_sorted = y_hat.flatten()[sort_idxs]
    sq_errs = (y_sorted - y_hat_sorted)**2
    cum_sq_errs = np.cumsum(sq_errs)
    cum_mses = cum_sq_errs / np.arange(1, y.size + 1, dtype="float")

    fig, ax1 = plt.subplots()
    ax1.plot(np.arange(y_sorted.size), cum_sq_errs, 'b-')
    ax1.set_xlabel("Codon rank by true output")
    ax1.set_ylabel("Sum of squared error", color='b')
    #for tl in ax1.get_yticklabels():
    #    tl.set_color('b')

    ax2 = ax1.twiny()
    x2tick_location = ax1.xaxis.get_ticklocs()
    ax2.set_xticks(x2tick_location)
    x2ticklabels = np.array([round(y_sorted[int(loc)], 2) for loc in x2tick_location if loc < y.size])
    ax2.set_xticklabels(x2ticklabels)

    #ax3 = ax1.twinx()
    #ax3.plot(np.arange(y.size), cum_mses, 'r-')
    #ax3.set_ylabel('Cumulative MSEs', color='r')
    #for tl in ax3.get_yticklabels():
    #    tl.set_color('r')
    
    plt.title(title, y = 1.06)
    plt.savefig(fname)
    plt.close()

def make_err_thresh_plot(y, y_hat, title, fname):
    num_codons = y.size
    sort_idxs = np.argsort(y.flatten())
    y_sorted = y.flatten()[sort_idxs]
    y_hat_sorted = y_hat.flatten()[sort_idxs]
    sq_errs = (y_sorted - y_hat_sorted)**2

    sq_errs_sorted = np.sort(sq_errs)
    log10_sq_errs_sorted = np.log10(sq_errs_sorted)
    frac_below_err = np.arange(y_sorted.size)/float(y_sorted.size)   

    plt.semilogx(sq_errs_sorted, frac_below_err)
    #ax.set_yscale('log')
    plt.xlabel("Squared error threshold")
    plt.ylabel("Frac. codons below error threshold")    

    #xmin = np.min(log10_sq_errs_sorted)
    #frac = 0.9
    #plt.axhline(frac, color="red")
    #threshold = log10_sq_errs_sorted[int(frac * num_codons)]
    #plt.text(xmin, frac, "Error threshold = {0}".format(round(10**threshold, 4)))
    #frac = 0.8
    #plt.axhline(frac, color="red")
    #threshold = log10_sq_errs_sorted[int(frac * num_codons)]
    #plt.text(xmin, frac, "{0}".format(round(10**threshold, 4)))
    #frac = 0.7
    #plt.axhline(frac, color="red")
    #threshold = log10_sq_errs_sorted[int(frac * num_codons)]
    #plt.text(xmin, frac, "{0}".format(round(10**threshold, 4)))
    #frac = 0.6
    #plt.axhline(frac, color="red")
    #threshold = log10_sq_errs_sorted[int(frac * num_codons)]
    #plt.text(xmin, frac, "{0}".format(round(10**threshold, 4)))

    plt.title(title)
    plt.savefig(fname)
    plt.close()

def make_mse_v_quant_cutoff_plot(y, y_hat, title, fname):
    sort_idxs = np.argsort(y.flatten())
    y_sorted = y.flatten()[sort_idxs]
    y_hat_sorted = y_hat.flatten()[sort_idxs]
    sq_errs = (y_sorted - y_hat_sorted)**2
    cum_sq_errs = np.cumsum(sq_errs)

    fig, ax1 = plt.subplots()
    ax1.plot(np.arange(y_sorted.size), cum_sq_errs, 'b-')
    ax1.set_xlabel("Codon rank by true output")
    ax1.set_ylabel("Sum of squared error", color='b')
    ax2 = ax1.twiny()
    x2tick_location = ax1.xaxis.get_ticklocs()
    ax2.set_xticks(x2tick_location)
    x2ticklabels = np.array([round(y_sorted[int(loc)], 2) for loc in x2tick_location if loc < y.size])
    ax2.set_xticklabels(x2ticklabels)

    plt.title(title, y = 1.06)
    plt.savefig(fname)
    plt.close()

def get_submatrix_frobenius_norms(weight_mtx, num_sub, width_sub, start_idx):
    fnorms = []
    for i in range(num_sub):
        sub_mtx_i = weight_mtx[:,start_idx + i * width_sub: start_idx + (i+1) * width_sub]
        fnorms.append((sub_mtx_i**2).sum())
    return fnorms

def make_weight_norms_plots(
        weight_mtx, title_prefix, fname_prefix, rel_cod_idxs=False, 
        rel_nt_idxs=False, rel_struc_idxs=False):
    feat_start_idx = 0
    if rel_cod_idxs:
        num_cods = len(rel_cod_idxs)
        cod_ss_weights = get_submatrix_frobenius_norms(
            weight_mtx, num_cods, 64, feat_start_idx)
        cod_tickskip =  num_cods / 15
        make_scatter_plot_with_xticks(
            rel_cod_idxs, cod_ss_weights, title_prefix + " Codon", 
            "Codon relative to A site", "Weight submatrix Frobenius norm", 
            fname_prefix + ".codon.pdf", tickskip=cod_tickskip)
        feat_start_idx += num_cods * 64
    if rel_nt_idxs:
        num_nts = len(rel_nt_idxs)
        nt_ss_weights = get_submatrix_frobenius_norms(
            weight_mtx, num_nts, 4, feat_start_idx)
        nt_tickskip = num_nts / 15
        make_scatter_plot_with_xticks(
            rel_nt_idxs, nt_ss_weights, title_prefix + " Nucleotide", 
            "Nucleotide relative to start of A site", 
            "Weight submatrix Frobenius norm", fname_prefix + ".nt.pdf", 
            tickskip=nt_tickskip)
        feat_start_idx += num_nts * 4
    if rel_struc_idxs:
        num_pos = len(rel_struc_idxs)
        struc_ss_weights = get_submatrix_frobenius_norms(
            weight_mtx, num_pos, 1, feat_start_idx)
        struc_tickskip = num_pos / 15
        make_scatter_plot_with_xticks(
            rel_struc_idxs, struc_ss_weights, title_prefix + "Structure", 
            "Nucleotide relative to start of A site", 
            "Weight submatrix Frobenius norm", fname_prefix + ".struc.pdf", 
            tickskip=struc_tickskip)

def make_cost_by_epoch_plot(cost_by_epoch, title, fname):
    num_epochs = len(cost_by_epoch)
    epoch_idxs = range(num_epochs)
    plt.plot(epoch_idxs, cost_by_epoch)
    xlab = "Epoch number"
    ylab = "Batch squared error cost"
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    plt.text(num_epochs * 0.6, cost_by_epoch[0] * 0.9, 
             "Final cost: {0}".format(round(cost_by_epoch[num_epochs-1], 4)))
    plt.savefig(fname)
    plt.close()

def unpickle_object(pkl_fname):
    with open(pkl_fname, "r") as f:
        obj = pickle.load(f)
    return obj

def unpickle_epoch(self, epoch):
    epoch_dir = "{0}/{1}/epoch{2}".format(self.out_dir, self.name, epoch)
    weights = self.unpickle_object(epoch_dir + "/weights.pkl")
    tr_cost_by_epoch = self.unpickle_object(epoch_dir + "/tr_cost_by_epoch.pkl")
    te_cost_by_epoch = self.unpickle_object(epoch_dir + "/te_cost_by_epoch.pkl")
    return weights, tr_cost_by_epoch, te_cost_by_epoch

def load_nn_data(out_dir, name, epoch):
    nn_dir = "{0}/{1}".format(out_dir, name) 
    init_dir = nn_dir + "/init_data"
    y_tr = unpickle_object(init_dir + "/y_tr.pkl")
    y_te = unpickle_object(init_dir + "/y_te.pkl")
    
    epoch_dir = "{0}/epoch{1}".format(nn_dir, epoch)
    weights = unpickle_object(epoch_dir + "/weights.pkl")
    tr_cost_by_epoch = unpickle_object(epoch_dir + "/tr_cost_by_epoch.pkl")
    te_cost_by_epoch = unpickle_object(epoch_dir + "/te_cost_by_epoch.pkl")
    y_tr_hat = unpickle_object(epoch_dir + "/y_tr_hat.pkl")
    y_te_hat = unpickle_object(epoch_dir + "/y_te_hat.pkl")
    
    return y_tr, y_te, y_tr_hat, y_te_hat, weights, tr_cost_by_epoch, te_cost_by_epoch
    
def make_plot_dir(out_dir, name, overwrite=False):
    dir_name = "{0}/{1}/plots".format(out_dir, name)
    if os.path.exists(dir_name) and not overwrite:
        print "WARNING: dir with name {0} already exists".format(dir_name)
        print "Run with keyword parameter overwrite=True to overwrite"
        raise NameError
    if os.path.exists(dir_name) and overwrite:
        return 
    else:
        os.makedirs(dir_name)

def make_all_plots(
        out_dir, name, epoch, scat_xlim, scat_ylim, scat_textpos, 
        rel_cod_idxs=False, rel_nt_idxs=False, rel_struc_idxs=False, 
        cost_by_epoch_plots=True, weight_plots=True, scat_plots=True, 
        overwrite=False, ext="pdf"):
    (y_tr, y_te, y_tr_hat, y_te_hat, weights, tr_cost_by_epoch, 
        te_cost_by_epoch) = load_nn_data(out_dir, name, epoch)
    make_plot_dir(out_dir, name, overwrite)
    plot_dir = "{0}/{1}/plots".format(out_dir, name)
    if cost_by_epoch_plots:
        make_cost_by_epoch_plot(
            tr_cost_by_epoch, name + " Training Cost by Epoch", 
            "{0}/tr_cost_by_epoch.pdf".format(plot_dir))
        make_cost_by_epoch_plot(
            te_cost_by_epoch, name + " Test Cost by Epoch", 
            "{0}/te_cost_by_epoch.pdf".format(plot_dir))
    if weight_plots:
        make_weight_norms_plots(
            weights[0], 
            name + " Input Weight Sq. Norm by", "{0}/ssweights".format(plot_dir), 
            rel_cod_idxs=rel_cod_idxs, rel_nt_idxs=rel_nt_idxs, 
            rel_struc_idxs=rel_struc_idxs)
    if scat_plots:
        make_density_scatter_plot(
            y_tr.flatten(), y_tr_hat.flatten(), 
            name + " Training Predicted vs. True Outputs", "True Outputs", 
            "Predicted Outputs", scat_xlim, scat_ylim, scat_textpos, 
            "{0}/tr_scatter.{1}".format(plot_dir, ext))
        make_density_scatter_plot(
            y_te.flatten(), y_te_hat.flatten(), 
            name + " Test Predicted vs. True Outputs", "True Outputs", 
            "Predicted Outputs", scat_xlim, scat_ylim, scat_textpos, 
            "{0}/te_scatter.{1}".format(plot_dir, ext))

def lasagne_make_all_plots(
        out_dir, name, epoch, scat_xlim, scat_ylim, scat_textpos, scat_bins=0, 
        rel_cod_idxs=False, rel_nt_idxs=False, rel_struc_idxs=False, 
        cost_by_epoch_plots=True, weight_plots=True, scat_plots=True, 
        cum_err_plots=True, overwrite=False, ext="pdf"):
    (y_tr, y_te, y_tr_hat, y_te_hat, weights, tr_cost_by_epoch,
        te_cost_by_epoch) = load_nn_data(out_dir, name, epoch)
    make_plot_dir(out_dir, name, overwrite)
    plot_dir = "{0}/{1}/plots".format(out_dir, name)
    if cost_by_epoch_plots:
        make_cost_by_epoch_plot(
            tr_cost_by_epoch, name + " Training Cost by Epoch", 
            "{0}/tr_cost_by_epoch.pdf".format(plot_dir))
        make_cost_by_epoch_plot(
            te_cost_by_epoch, name + " Test Cost by Epoch", 
            "{0}/te_cost_by_epoch.pdf".format(plot_dir))
    #if weight_plots:
    #    make_weight_norms_plots(
    #        weights[0].transpose(), 
    #        name + " Input Weight Sq. Norm by", "{0}/ssweights".format(plot_dir), 
    #        rel_cod_idxs=rel_cod_idxs, rel_nt_idxs=rel_nt_idxs, 
    #        rel_struc_idxs=rel_struc_idxs)
    if scat_plots:
        make_density_scatter_plot(
            y_tr.flatten(), y_tr_hat.flatten(), 
            name + " Training Predicted vs. True Scaled Cts.", "True Outputs", 
            "Predicted Outputs", scat_xlim, scat_ylim, scat_textpos,
            "{0}/tr_scatter.{1}".format(plot_dir, ext), bins=scat_bins)
        make_density_scatter_plot(
            y_te.flatten(), y_te_hat.flatten(), 
            name + " Test Predicted vs. True Outputs", "True Scaled Counts", 
            "Predicted Scaled Counts", scat_xlim, scat_ylim, scat_textpos,
            "{0}/te_scatter.{1}".format(plot_dir, ext), bins=scat_bins,
            no_title=True)
    if cum_err_plots:
        make_binned_err_plot(
            y_tr.flatten(), y_tr_hat.flatten(), 
            name + " Training Sq. Err. v. Sorted Outputs", 
            "{0}/tr_bin_err.{1}".format(plot_dir, ext), 700)
        #make_smoothed_err_plot(
        #    y_tr.flatten(), y_tr_hat.flatten(), 
        #    name + " Training Sq. Err. v. Sorted Outputs", 
        #    "{0}/tr_smoothed_err.{1}".format(plot_dir, ext))
        make_cum_err_plot(
            y_tr.flatten(), y_tr_hat.flatten(), 
            name + " Training Cumulative Sq. Err. v. Sorted Outputs", 
            "{0}/tr_cum_err.{1}".format(plot_dir, ext))
        #make_err_thresh_plot(
        #    y_tr.flatten(), y_tr_hat.flatten(), 
        #    name + " Training fraction of points below error threshold", 
        #    "{0}/tr_err_thresh.{1}".format(plot_dir, ext))
        make_binned_err_plot(
            y_te.flatten(), y_te_hat.flatten(), 
            name + " Test Sq. Err. v. Sorted Outputs", 
            "{0}/te_bin_err.{1}".format(plot_dir, ext), 700)
        #make_smoothed_err_plot(
        #    y_te.flatten(), y_te_hat.flatten(), 
        #    name + " Test Sq. Err. v. Sorted Outputs", 
        #    "{0}/te_smoothed_err.{1}".format(plot_dir, ext))
        make_cum_err_plot(
            y_te.flatten(), y_te_hat.flatten(), 
            name + " Test Cumulative Sq. Err. v. Sorted Outputs", 
            "{0}/te_cum_err.{1}".format(plot_dir, ext))
        #make_err_thresh_plot(
        #    y_te.flatten(), y_te_hat.flatten(), 
        #    name + " Test fraction of points below error threshold", 
        #    "{0}/te_err_thresh.{1}".format(plot_dir, ext))

def plot_preds_by_gene(
        tr_codons_fname, te_codons_fname, cts_by_codon_fname, outputs_fname, 
        cod_trunc_5p, cod_trunc_3p, nn, plot_dirname, save_gene_lists=False, 
        tr_genes_out_fname="tr_genes_density_sorted.txt", 
        te_genes_out_fname="te_genes_density_sorted.txt", 
        all_genes_out_fname="all_genes_density_sorted.txt", 
        min_tot_cts=200, min_cod_w_cts=100):

    gene2sym_fname = "/mnt/lareaulab/lareau/StanfordFootprints/GenomeEtc/" +\
                     "S_cer_id_symbol.txt"
    gene2symbol = proc.get_gene2symbol_dict(
        gene2sym_fname, lambda x: x + "_13cds10")
    
    #Load codon sets, cts_by_codon, outputs, y_tr/te
    cts_by_codon = load_cts_by_codon(cts_by_codon_fname)
    y_tr, y_te, tr_codon_set, te_codon_set, outputs = get_y_data(
        tr_codons_fname, te_codons_fname, outputs_fname)

    tr_sorted_genes = sorted(tr_codon_set.keys())
    te_sorted_genes = sorted(te_codon_set.keys())

    tr_gene_data = get_mat_idxs_by_gene(tr_sorted_genes, tr_codon_set)
    te_gene_data = get_mat_idxs_by_gene(te_sorted_genes, te_codon_set)

    #Sort all genes by data sensity, save list if indicated
    gene_set = cts_by_codon.keys()
    gene_set = filter(lambda gene: len(cts_by_codon[gene]) > (cod_trunc_5p +
                               cod_trunc_3p), gene_set)
    gene_set = filter(lambda gene: has_enough_cts(cts_by_codon[gene][cod_trunc_5p:-cod_trunc_3p], min_tot_cts, min_cod_w_cts), gene_set)
    all_genes_by_density = sort_genes_by_density(gene_set, cts_by_codon,
                            cod_trunc_5p, cod_trunc_3p, descend=True)

    #Sort training genes by data density, save list if indicated
    tr_genes_by_density = sort_genes_by_density(
        tr_sorted_genes, cts_by_codon, cod_trunc_5p, cod_trunc_3p)
    te_genes_by_density = sort_genes_by_density(
        te_sorted_genes, cts_by_codon, cod_trunc_5p, cod_trunc_3p)

    #Sort genes by data density, save list if indicated
    if save_gene_lists:
        save_gene_list(all_genes_by_density, all_genes_out_fname)
        save_gene_list(tr_genes_by_density, tr_genes_out_fname)
        save_gene_list(te_genes_by_density, te_genes_out_fname)

    #Make dir to save gene plots, if needed
    if not os.path.exists(plot_dirname):
        os.makedirs(plot_dirname)

    #Make gene plots
    for idx, gene in enumerate(te_genes_by_density):
        num_codons = te_gene_data[gene]["length"]
        start_idx = te_gene_data[gene]["start_idx"]
        y_gene = nn.y_te.flatten()[start_idx: start_idx + num_codons]
        y_hat_gene = nn.y_te_hat.flatten()[start_idx: start_idx + num_codons]

        fig, ax = plt.subplots()
        bars = ax.bar(range(num_codons), y_gene)
        line = ax.plot(np.arange(num_codons) + .45 , y_hat_gene, color="r")
        out_fname = "{0}/{1}_{2}_pred_v_true.pdf".format(plot_dirname, idx, gene)
        plt.title(gene2symbol.get(gene, gene[:-8]))
        plt.savefig(out_fname)
        plt.close()

    return

def plot_preds_by_gene_2(
        nn, te_set_data_table, gene2sym, plot_dir):
    #Make dir to save gene plots, if needed
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    #Get te_genes sorted by density in descending order
    te_genes_by_density = get_genes_by_density2(te_set_data_table)
    te_genes_by_pos = np.array([row[0] for row in te_set_data_table])
    #outputs_by_pos = np.array([row[4] for row in te_set_data_table])
    for gene_idx, gene in enumerate(te_genes_by_density):
        gene_row_bools = (te_genes_by_pos == gene)
        y_idxs = np.where(gene_row_bools)
        #for i in range(len(te_set_row_labels)):
        #    if te_set_row_labels[i][0] == gene:
        #        y_idxs.append(i) 
        out_fname = "{0}/{1}_{2}_pred_v_true.pdf".format(
            plot_dir, gene_idx, gene)
        y = nn.y_te.take(y_idxs, axis=0)
        y_hat = nn.y_te_hat.take(y_idxs, axis=0)
        gene_sym = gene2sym.get(gene, gene)
        make_preds_by_gene_plot(y, y_hat, gene, gene_sym, out_fname)
        
def make_preds_by_gene_plot(y, y_hat, gene_name, gene_sym, out_fname):
    num_codons = y.size
    fig, ax = plt.subplots()
    #print np.arange(num_codons).shape
    #print y.ravel().shape
    bars = ax.bar(np.arange(num_codons), y.ravel())
    line = ax.plot(np.arange(num_codons) + .45 , y_hat.ravel(), color="r")
    plt.title(gene_sym)
    plt.savefig(out_fname)
    plt.close()

def make_cts_by_codon_freq_plot(
        cts, title, xlab, ylab, out_fname, xlims=False):
    """
        Inputs: 
            cts - numpy array of cts for each codon in a set (training?)
            xlims - optional list of 2 integers for max/min cts
            rest of arguments are strings
    """
    # Make cts_freq_dict
    min_cts = cts.min()
    max_cts = cts.max()
    cts_vals = np.unique(cts)
    cts_freq_dict = {}
    for i in cts_vals:
        cts_freq_dict[i] = (cts == i).sum()

    # Make plot
    items = cts_freq_dict.items()
    keys = [item[0] for item in items]
    values = [item[1] for item in items]
    plt.scatter(keys, values)
    if xlims:
        plt.xlim(xlims)
        title = title + ", truncated"
        out_fname = out_fname[:-4] + ".{0}-{1}.pdf".format(xlims[0], xlims[1])
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.savefig(out_fname)
    plt.close()

def make_outputs_freq_plot(
        outputs, title, xlab, ylab, out_fname, xlims=False, num_bins=500):
    """
        Inputs: 
            outputs - numpy array of outputs for each codon in a set (training?)
            xlims - optional list of 2 integers for max/min cts
            num_bins - # of bins in which to put outputs values
            rest of arguments are strings
    """
    # Make cts_freq_dict
    min_output = outputs.min()
    max_output = outputs.max()
    bins = np.linspace(min_output, max_output, num_bins)
    width = bins[1] - bins[0]
    output_freq_dict = {}
    for i in bins:
        cts = (np.vstack([outputs > i, outputs <= i + width]).all(axis=0)).sum()
        output_freq_dict[i] = cts

    # Make plot
    items = output_freq_dict.items()
    keys = [item[0] for item in items]
    values = [item[1] for item in items]
    plt.scatter(keys, values)
    if xlims:
        plt.xlim(xlims)
        title = title + ", truncated"
        out_fname = out_fname[:-4] + ".{0}-{1}.pdf".format(xlims[0], xlims[1])
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.savefig(out_fname)
    plt.close()

def make_log_outputs_freq_plot(
        outputs, title, xlab, ylab, out_fname, num_bins=500):
    """
        Inputs: 
            cts - numpy array of cts for each codon in a set (training?)
            num_bins - # of bins in which to put outputs values
            rest of arguments are strings
    """
    # Make output_freq_dict
    log_outputs = np.log(outputs)
    min_output = log_outputs.min()
    max_output = log_outputs.max()
    bins = np.linspace(min_output, max_output, num_bins)
    width = bins[1] - bins[0]
    output_freq_dict = {}
    for i in bins:
        cts = (np.vstack([log_outputs > i, log_outputs <= i + width]).all(
            axis=0)).sum()
        output_freq_dict[i] = cts

    # Make plot
    items = output_freq_dict.items()
    keys = [item[0] for item in items]
    values = [item[1] for item in items]
    plt.scatter(keys, values)
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    #NOTE: Remove this later
    plt.xlim([-12, 6])

    mu = log_outputs.mean()
    sigma = log_outputs.std()
    var = sigma * sigma

    x = np.linspace(mu - 5 * sigma, mu + 5 * sigma, 1000)
    y = mlab.normpdf(x, mu, sigma)
    scale_factor = max(output_freq_dict.values())/max(y)
    
    plt.plot(x, scale_factor * y, "red")
    plt.savefig(out_fname)
    plt.close()

