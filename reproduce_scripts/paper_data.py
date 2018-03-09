import sys
import pickle
import numpy as np
from scipy.stats import pearsonr, spearmanr

if __name__ == "__main__":
    nn_dir = sys.argv[1]
    nn_epoch = int(sys.argv[2])
    feat_nb_mses_fname = sys.argv[3]
    # NOTE: New
    feat_nb_corrs_fname = sys.argv[4]
    leaveout_mses_fname = sys.argv[5]
    # NOTE: New
    leaveout_corrs_fname = sys.argv[6]
    struc_mses_fname = sys.argv[7]
    # NOTE: New
    struc_corrs_fname = sys.argv[8]
    weinberg_cod_scores_fname = sys.argv[9] #.tsv file, not .pkl
    green_cod_scores_fname = sys.argv[10] #.tsv file, not .pkl
    weinberg_28_cod_scores_fname = sys.argv[11] #.tsv file, not .pkl
    green_28_cod_scores_fname = sys.argv[12] #.tsv file, not .pkl
    lareau_cod_scores_fname = sys.argv[13] #.tsv file, not .pkl
    iwasaki_cod_scores_fname = sys.argv[14] #.tsv file, not .pkl
    out_fname = sys.argv[15]

    # Open out_file for writing results
    out_file = open(out_fname, "w")
    out_file.write("# Collected results included in text of Tunney et al. 2017\n\n")

    # Get full model w. struc R and MSE values
    #NOTE: Change this to average MSE in the future?
    y_te_fname = nn_dir + "/init_data/y_te.pkl"
    y_te_hat_fname = nn_dir + "/epoch{0}/y_te_hat.pkl".format(nn_epoch)
    y_te = pickle.load(open(y_te_fname)).ravel()
    y_te_hat = pickle.load(open(y_te_hat_fname)).ravel()
    pear_r, pear_p = pearsonr(y_te, y_te_hat)
    spear_r, spear_p = spearmanr(y_te, y_te_hat)
    mse = ((y_te - y_te_hat)**2).mean()

    out_file.write("Complete model with structure features:\n")
    out_file.write("MSE: {0}\n".format(mse))
    out_file.write("Pearson r: {0}\t Spearman r: {1}\n".format(pear_r, spear_r))
    out_file.write("\n")

    struc_mses_file = open(struc_mses_fname, "r")
    avg_mse = "MISSING"
    std_mse = "MISSING"
    for line in struc_mses_file:
        line = line.strip().split()
        if line[0] == "str_n17n15_cod_n5p4_nt_n15p14":
            avg_mse = float(line[11])
            std_mse = float(line[12])
        if line[0] == "max_str_p13p42_cod_n5p4_nt_n15p14":
            max_str_avg_mse = float(line[11])
   
    struc_corrs_file = open(struc_corrs_fname, "r")
    struc_corrs_file.readline()
    for line in struc_corrs_file:
        line = line.strip().split()
        if line[0] == "str_n17n15_cod_n5p4_nt_n15p14":
            med_corr = np.median([float(elt) for elt in line[1:]])
        if line[0] == "max_str_p13p42_cod_n5p4_nt_n15p14":
            max_str_med_corr = np.median([float(elt) for elt in line[1:]])

    out_file.write("Complete model with structure features rep_series:\n")
    out_file.write("Avg. MSE: {0}\n".format(avg_mse))
    out_file.write("Std. MSE: {0}\n".format(std_mse))
    out_file.write("Median pearson's r: {0}\n".format(med_corr))
    out_file.write("\n")

    out_file.write("Complete model with downstream struc windows (p13p42) rep_series:\n")
    out_file.write("Avg. MSE: {0}\n".format(max_str_avg_mse))
    out_file.write("Median pearson's r: {0}\n".format(max_str_med_corr))
    out_file.write("\n")

    # MSE for points with true scaled counts < 2:
    idxs_y_leq_2 = np.where(y_te <= 2)[0]
    y_te_leq_2_cts = idxs_y_leq_2.size
    y_te_cts = y_te.size
    y_te_leq_2_frac = float(y_te_leq_2_cts)/y_te_cts

    y_te_leq_2 = y_te[idxs_y_leq_2]
    y_te_hat_leq_2 = y_te_hat[idxs_y_leq_2]

    y_te_leq_2_mse = ((y_te_leq_2 - y_te_hat_leq_2)**2).mean()
    y_te_leq_2_pcorr = pearsonr(y_te_leq_2, y_te_hat_leq_2)[0]

    out_file.write("Test set analysis for true scaled cts leq 2\n")
    out_file.write("Pts with true s.c. leq 2: {0}\n".format(y_te_leq_2_cts))
    out_file.write("Total test set pts: {0}\n".format(y_te_cts))
    out_file.write("Fraction pts with true s.c. leq 2: {0}\n".format(y_te_leq_2_frac))
    out_file.write("MSE pts with true s.c. leq 2: {0}\n".format(y_te_leq_2_mse))
    out_file.write("Pearson r pts with true s.c. leq 2: {0}\n".format(y_te_leq_2_pcorr))
    out_file.write("\n")

    # Feature neighborhood MSES, and + structure
    feat_nb_mses_file = open(feat_nb_mses_fname, "r")
    feat_nb_corrs_file = open(feat_nb_corrs_fname, "r")
    model_names = []
    avg_mses = []
    med_corrs = []

    # Skip header lines
    feat_nb_mses_file.readline()
    feat_nb_mses_file.readline()
    feat_nb_corrs_file.readline()

    for line in feat_nb_mses_file:
        line = line.strip().split()
        model_name = line[0]
        avg_mse = float(line[11])
        model_names.append(model_name)
        avg_mses.append(avg_mse)

    for idx, line in enumerate(feat_nb_corrs_file):
        line = line.strip().split()
        model_name = line[0]
        assert model_names[idx] == model_name, "model names in mses and corrs file don't match!"
        med_corr = np.median([float(elt) for elt in line[1:]])
        med_corrs.append(med_corr)
            
    out_file.write("Average MSEs and median correlations by model:\n")
    out_file.write("model_name\tMSE\tpearson_r\n")
    for i in range(len(model_names)):
        out_file.write("{0}\t{1}\t{2}\n".format(model_names[i], avg_mses[i], med_corrs[i]))
    out_file.write("\n")

    # Write difference between correlations for -5 to +4 model with and without nt features

    n5p4_n15p14_idx = model_names.index('full_cod_n5p4_nt_n15p14')
    n5p4_idx = model_names.index('full_cod_n5p4')
    mse_diff = avg_mses[n5p4_idx] - avg_mses[n5p4_n15p14_idx]
    corr_diff = med_corrs[n5p4_n15p14_idx] - med_corrs[n5p4_idx]
    out_file.write("Improvement in MSE and pearson r from full_cod_n5p4 to full_cod_n5p4_nt_n15p14\n")
    out_file.write("Avg. MSE difference: {0}\n".format(mse_diff))
    out_file.write("Median pearson r difference: {0}\n".format(corr_diff))
    out_file.write("\n")

    # Store for differences in leave-one-out. 
    # NOTE: Change this if base model changes
    n7p5_n21p17_idx = model_names.index('full_cod_n7p5_nt_n21p17')
    n7p5_n21p17_corr = med_corrs[n7p5_n21p17_idx]

    #NOTE: need to add in structure series MSEs

    # Leaveout series MSEs
    leaveout_mses_file = open(leaveout_mses_fname, "r")
    leaveout_corrs_file = open(leaveout_corrs_fname, "r")
    #NOTE: Assumes -7 to +5 codon feature neighborhood
    cod_idxs = range(-7, 6)
    model_names = []
    delta_mses = []
    med_corrs = []
    delta_corrs = []

    leaveout_mses_file.readline()
    leaveout_mses_file.readline()
    leaveout_corrs_file.readline()
    for line in leaveout_mses_file:
        line = line.strip().split()
        model_name = line[0]
        model_names.append(model_name)
        delta_mse = float(line[13])
        delta_mses.append(delta_mse)

    for idx, line in enumerate(leaveout_corrs_file):
        line = line.strip().split()
        model_name = line[0]
        assert model_names[idx] == model_name, "model names in mses and corrs file don't match!"
        med_corr = np.median([float(elt) for elt in line[1:]])
        med_corrs.append(med_corr)

    out_file.write("Weinberg leaveout series delta MSE and pearson r by codon\n")
    out_file.write("codon_idx\tdelta_MSE\tpearson_r\tdelta_r\n")
    for i in range(len(cod_idxs)):
        out_file.write("{0}\t{1}\t{2}\t{3}\n".format(cod_idxs[i], delta_mses[i], med_corrs[i], n7p5_n21p17_corr - med_corrs[i]))
    out_file.write("\n")
    
    ### Codon scores analysis
    weinberg_cod_scores = np.loadtxt(weinberg_cod_scores_fname)
    green_cod_scores = np.loadtxt(green_cod_scores_fname)
    weinberg_28_cod_scores = np.loadtxt(weinberg_28_cod_scores_fname)
    green_28_cod_scores = np.loadtxt(green_28_cod_scores_fname)
    print lareau_cod_scores_fname
    lareau_cod_scores = np.loadtxt(lareau_cod_scores_fname)
    iwasaki_cod_scores = np.loadtxt(iwasaki_cod_scores_fname)

    alpha = "ACGT"
    cods = [x + y + z for x in alpha for y in alpha for z in alpha]
    cods = np.array(cods)
    #NOTE: Assumes -5 to +4 codon feature neighborhood
    weinberg_A_site = weinberg_cod_scores[:,5].ravel()
    weinberg_A_site_sort_idxs = np.argsort(weinberg_A_site)
    weinberg_A_site_sorted = weinberg_A_site[weinberg_A_site_sort_idxs]
    weinberg_cods_sorted = cods[weinberg_A_site_sort_idxs]
    
    out_file.write("Weinberg codons sorted by A site nn codon scores\n")
    out_file.write("codon\tcodon_score\n")
    for i in range(len(weinberg_cods_sorted)):
        out_file.write("{0}\t{1}\n".format(weinberg_cods_sorted[i], weinberg_A_site_sorted[i]))
    out_file.write("\n")

    num_cods = lareau_cod_scores.shape[1]
    lareau_iwasaki_corrs_by_cod = []
    for i in range(num_cods):
        keep_idxs = np.logical_not(np.isnan(lareau_cod_scores[:,i]))
        lareau_scores_i = lareau_cod_scores[:,i][keep_idxs].ravel()
        iwasaki_scores_i = iwasaki_cod_scores[:,i][keep_idxs].ravel()
        #print lareau_scores_i
        #print iwasaki_scores_i
        #print lareau_scores_i.size
        #print iwasaki_scores_i.size
        pear_r, pear_p = pearsonr(lareau_scores_i, iwasaki_scores_i)
        lareau_iwasaki_corrs_by_cod.append(pear_r)
    
    #NOTE: Assumes -5 to +4 codon feature neighborhood
    cod_idxs = range(-5, 4)
    out_file.write("Lareau/Iwasaki s28_cod_n5p4_nt_n15p14 codon score correlations by codon position\n")
    out_file.write("codon_idx\tpear_r\n")
    for i, cod_idx in enumerate(cod_idxs):
        out_file.write(
            "{0}\t{1}\n".format(cod_idx, lareau_iwasaki_corrs_by_cod[i]))
    out_file.write("\n")

    num_cods = weinberg_28_cod_scores.shape[1]
    weinberg_green_28_corrs_by_cod = []
    for i in range(num_cods):
        keep_idxs = np.logical_not(np.isnan(weinberg_28_cod_scores[:,i]))
        weinberg_scores_i = weinberg_cod_scores[:,i][keep_idxs].ravel()
        green_scores_i = green_28_cod_scores[:,i][keep_idxs].ravel()
        pear_r, pear_p = pearsonr(weinberg_scores_i, green_scores_i)
        weinberg_green_28_corrs_by_cod.append(pear_r)
    
    out_file.write("Weinberg/Green s28_cod_n5p4_nt_n15p14 codon score correlations by codon position\n")
    out_file.write("codon_idx\tpear_r\n")
    for i, cod_idx in enumerate(cod_idxs):
        out_file.write(
            "{0}\t{1}\n".format(cod_idx, weinberg_green_28_corrs_by_cod[i]))
    out_file.write("\n")

    weinberg_green_corrs_by_cod = []
    for i in range(num_cods):
        keep_idxs = np.logical_not(np.isnan(weinberg_cod_scores[:,i]))
        weinberg_scores_i = weinberg_cod_scores[:,i][keep_idxs].ravel()
        green_scores_i = green_cod_scores[:,i][keep_idxs].ravel()
        pear_r, pear_p = pearsonr(weinberg_scores_i, green_scores_i)
        weinberg_green_corrs_by_cod.append(pear_r)
    
    out_file.write("Weinberg/Green full_cod_n5p4_nt_n15p14 codon score correlations by codon position\n")
    out_file.write("codon_idx\tpear_r\n")
    for i, cod_idx in enumerate(cod_idxs):
        out_file.write(
            "{0}\t{1}\n".format(cod_idx, weinberg_green_corrs_by_cod[i]))

    out_file.close()

