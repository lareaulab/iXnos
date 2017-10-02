import sys
import pickle
import numpy as np
from scipy.stats import pearsonr, spearmanr

if __name__ == "__main__":
    nn_dir = sys.argv[1]
    nn_epoch = int(sys.argv[2])
    feat_nb_mses_fname = sys.argv[3]
    leaveout_mses_fname = sys.argv[4]
    struc_mses_fname = sys.argv[5]
    weinberg_cod_scores_fname = sys.argv[6] #.tsv file, not .pkl
    lareau_cod_scores_fname = sys.argv[7] #.tsv file, not .pkl
    iwasaki_cod_scores_fname = sys.argv[8] #.tsv file, not .pkl
    out_fname = sys.argv[9]

    # Open out_file for writing results
    out_file = open(out_fname, "w")
    out_file.write("# Collected results included in text of Tunney et al. 2017\n\n")

    # Get full model w. struc R and MSE values
    #NOTE: Change this to average MSE in the future?
    #NOTE: need to add in structure series MSEs
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
        if line[0] == "str_n17n15_cod_n7p5_nt_n21p17":
            avg_mse = float(line[11])
            std_mse = float(line[12])
    
    out_file.write("Complete model with structure features rep_series:\n")
    out_file.write("Avg. MSE: {0}\n".format(avg_mse))
    out_file.write("Std. MSE: {0}\n".format(std_mse))
    out_file.write("\n")

    # MSE for points with true scaled counts < 2:
    idxs_y_leq_2 = np.where(y_te <= 2)[0]
    y_te_leq_2_cts = idxs_y_leq_2.size
    y_te_cts = y_te.size
    y_te_leq_2_frac = float(y_te_leq_2_cts)/y_te_cts

    y_te_leq_2 = y_te[idxs_y_leq_2]
    y_te_hat_leq_2 = y_te_hat[idxs_y_leq_2]

    y_te_leq_2_mse = ((y_te_leq_2 - y_te_hat_leq_2)**2).mean()

    out_file.write("Test set analysis for true scaled cts leq 2\n")
    out_file.write("Pts with true s.c. leq 2: {0}\n".format(y_te_leq_2_cts))
    out_file.write("Total test set pts: {0}\n".format(y_te_cts))
    out_file.write("Fraction pts with true s.c. leq 2: {0}\n".format(y_te_leq_2_frac))
    out_file.write("MSE pts with true s.c. leq 2: {0}\n".format(y_te_leq_2_mse))
    out_file.write("\n")

    # Feature neighborhood MSES, and + structure
    feat_nb_mses_file = open(feat_nb_mses_fname, "r")
    model_names = []
    avg_mses = []
    
    feat_nb_mses_file.readline()
    feat_nb_mses_file.readline()
    for line in feat_nb_mses_file:
        line = line.strip().split()
        model_name = line[0]
        avg_mse = float(line[11])
        model_names.append(model_name)
        avg_mses.append(avg_mse)

    out_file.write("Average MSEs by model:\n")
    out_file.write("model_name\tMSE\n")
    for i in range(len(model_names)):
        out_file.write("{0}\t{1}\n".format(model_names[i], avg_mses[i]))
    out_file.write("\n")

    #NOTE: need to add in structure series MSEs

    # Leaveout series MSEs
    leaveout_mses_file = open(leaveout_mses_fname, "r")
    #NOTE: Assumes -7 to +5 codon feature neighborhood
    cod_idxs = range(-7, 6)
    delta_mses = []

    leaveout_mses_file.readline()
    leaveout_mses_file.readline()
    for line in leaveout_mses_file:
        line = line.strip().split()
        delta_mse = float(line[13])
        delta_mses.append(delta_mse)

    out_file.write("Weinberg leaveout series delta MSE by codon\n")
    out_file.write("codon_idx\tdelta_MSE\n")
    for i in range(len(cod_idxs)):
        out_file.write("{0}\t{1}\n".format(cod_idxs[i], delta_mses[i]))
    out_file.write("\n")
    
    ### Codon scores analysis
    weinberg_cod_scores = np.loadtxt(weinberg_cod_scores_fname)
    lareau_cod_scores = np.loadtxt(lareau_cod_scores_fname)
    iwasaki_cod_scores = np.loadtxt(iwasaki_cod_scores_fname)

    alpha = "ACGT"
    cods = [x + y + z for x in alpha for y in alpha for z in alpha]
    cods = np.array(cods)
    #NOTE: Assumes -7 to +5 codon feature neighborhood
    weinberg_A_site = weinberg_cod_scores[:,7].ravel()
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
        spear_r, spear_p = spearmanr(lareau_cod_scores[:,i].ravel(), iwasaki_cod_scores[:,i].ravel())
        lareau_iwasaki_corrs_by_cod.append(spear_r)
    
    out_file.write("Lareau/Iwasaki codon score correlations by codon position\n")
    out_file.write("codon_idx\tspear_r\n")
    for i, cod_idx in enumerate(cod_idxs):
        out_file.write("{0}\t{1}\n".format(cod_idx, lareau_iwasaki_corrs_by_cod[i]))
   
    out_file.close()

