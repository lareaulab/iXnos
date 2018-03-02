import random
import copy
import sys
import numpy as np
from scipy.stats import spearmanr, pearsonr
from matplotlib import pyplot as plt

import iXnos.interface as inter
import iXnos.process as proc

if __name__ == "__main__":

    model_names = ["str_n17n15_cod_n5p4_nt_n15p14"]

    model_name = sys.argv[1]
    model_rep = sys.argv[2]
    epoch = int(sys.argv[3])
    nn_parent_dir = sys.argv[4]
    gene_len_fname = sys.argv[5]
    gene_seq_fname = sys.argv[6]
    cts_by_codon_fname = sys.argv[7]
    outputs_fname = sys.argv[8]
    struc_fname = sys.argv[9]
    tr_bounds_fname = sys.argv[10]
    paralog_groups_fname = sys.argv[11]
    out_fname = sys.argv[12]

#    gene_len_fname = genome_dir + "/scer.transcripts.13cds10.lengths.txt"
#    gene_seq_fname = genome_dir + "/scer.transcripts.13cds10.fa"
#    struc_fname = struc_dir + "/scer.13cds10.windows.30len.fold"

#    cts_by_codon_fname = wein_proc_dir + "/cts_by_codon.size.27.31.txt"
#    outputs_fname = wein_proc_dir + "/outputs.size.27.31.txt"
#    paralog_groups_fname = genome_dir + "/scer_100_80.paralogs.id.txt"


    assert model_name in model_names, "model name {0} not in " \
        "list of models".format(model_name)

    if model_name == "str_n17n15_cod_n5p4_nt_n15p14":
        # Init vars
        rel_cod_idxs = range(-5, 5)
        rel_nt_idxs = range(-15, 15)
        rel_struc_idxs = range(-17, -14)

        nn_dir = nn_parent_dir + "/" + model_name +\
            "_rep{0}".format(model_rep)

        cod_trunc_5p = 20
        cod_trunc_3p = 20
        min_tot_cts = 200
        min_cod_w_cts = 100

        # Load transcriptome dicts
        len_dict = proc.get_len_dict(gene_len_fname)
        cds_dict = proc.get_cds_dict(gene_seq_fname, len_dict)
        struc_dict = proc.get_struc_dict(struc_fname)

        # Compute cts_by_codon and outputs
        cts_by_codon = proc.load_cts_by_codon(cts_by_codon_fname)
        outputs = proc.load_outputs(outputs_fname)
        paralog_groups = proc.load_paralog_groups(paralog_groups_fname)

        # Get gene list to compute performance metrics
        gene_set = cts_by_codon.keys()
        # Filter genes that are shorter than truncation regions
        gene_set = filter(lambda gene: len(cts_by_codon[gene]) > (cod_trunc_5p +
            cod_trunc_3p), gene_set)
        # Filter genes that don't have enough counts to meet cutoffs
        gene_set = filter(lambda gene: proc.has_enough_cts( \
            cts_by_codon[gene][cod_trunc_5p:-cod_trunc_3p], min_tot_cts, \
            min_cod_w_cts), gene_set)
        # Sort genes by footprint density
        genes_by_density = proc.sort_genes_by_density(
            gene_set, cts_by_codon, cod_trunc_5p, cod_trunc_3p, descend=True)
        # NOTE: I'm not sure if I want to filter out paralogs here. 
        #       We're interested in testing performance on all genes
        #       Or could this muddle what is/isn't a training gene?
        ## Filter out all but first example from paralog groups
        genes_by_density = proc.filter_paralogs(
            genes_by_density, paralog_groups) 

        # Load NN model
        my_nn = inter.load_lasagne_feedforward_nn(nn_dir, epoch)

        # Make dictionary of training genes:
        tr_genes_dict = {}
        with open(tr_bounds_fname, "r") as f:
            for line in f:
                gene_name = line.strip().split()[0]
                tr_genes_dict[gene_name] = True

        # Init vectors to store performance metrics
        # All vectors are for genes in order of fp density
        trunc_fp_density_by_density = []
        training_gene_bools = []
        mses_by_density = []
        pearson_r_by_density = []
        pearson_r_subsample_1000 = []
        pearson_r_subsample_2000 = []
        pearson_r_subsample_3000 = []
        pearson_r_subsample_4000 = []
        
        # For each gene
        for idx, gene_name in enumerate(genes_by_density):
            # Mark gene as training or test
            training_gene_bools.append(tr_genes_dict.get(gene_name, False))
            # Get number of codons in full gene
            num_cods = len(cts_by_codon[gene_name])
            # Get codon indices for truncated gene
            codon_set = {gene_name: [cod for cod in range(cod_trunc_5p, num_cods - cod_trunc_3p)]}
            # Get fp density on truncated gene, and save
            trunc_fp_density = proc.get_fp_density(
                cts_by_codon[gene_name], cod_trunc_5p, cod_trunc_3p)
            trunc_fp_density_by_density.append(trunc_fp_density)
            # Compute feature matrix
            X_te = proc.get_X(
                codon_set, cds_dict, rel_cod_idxs=rel_cod_idxs,
                rel_nt_idxs=rel_nt_idxs, rel_struc_idxs=rel_struc_idxs,
                struc_dict=struc_dict).transpose()
            # Get true and predicted scaled counts
            y_te = proc.get_y(codon_set, outputs).transpose()
            y_te_hat = my_nn.pred_fn(X_te)
            # Compute MSE for gene, append to mses_by_density
            mse = ((np.array(y_te) - np.array(y_te_hat))**2).sum() / y_te.size
            mses_by_density.append(mse)
            # Compute pearson corr. for gene, append to pearson_r_by_density
            pcorr = pearsonr(y_te, y_te_hat)[0][0]
            pearson_r_by_density.append(pcorr)
            
        for benchmark_gene_rank in [1000, 2000, 3000, 4000]:
            benchmark_gene = genes_by_density[benchmark_gene_rank]
            benchmark_gene_density = proc.get_fp_density(
                cts_by_codon[benchmark_gene], cod_trunc_5p, cod_trunc_3p)
            # Make dictionary for subsampled counts by codon
            subsampled_cbc = copy.deepcopy(cts_by_codon)
            # Subsample gene cts by codon and put in subsampled_cbc
            for rank, gene_name in enumerate(genes_by_density):
                # Skip genes with less dense data than benchmark gene
                if rank >= benchmark_gene_rank:
                    continue
                # Name, cts_by_codon of gene being subsampled
                gene_cbc = subsampled_cbc[gene_name]
                # Length of truncated region, target cts for that region
                gene_len = len(gene_cbc)
                trunc_len = gene_len - cod_trunc_5p - cod_trunc_3p
                target_cts = benchmark_gene_density * trunc_len
                # Vector for subsampled cts
                sample = [0 for i in range(gene_len)]
                # Expanded cts vector for sampling
                expanded_cbc = []
                for i in range(cod_trunc_5p, cod_trunc_5p + trunc_len):
                    cod_i_cts = gene_cbc[i]
                    while cod_i_cts > 0:
                        if cod_i_cts >= 1:
                            expanded_cbc.append((i, 1))
                        else:
                            expanded_cbc.append((i, cod_i_cts))
                        cod_i_cts -= 1
                gene_cts = 0
                random.shuffle(expanded_cbc)
                # Idx to iterate over cod_idx, count tuples
                j = 0
                while gene_cts < target_cts:
                    cod_idx, ct = expanded_cbc[j]
                    sample[cod_idx] += ct
                    gene_cts += ct
                    j += 1
                subsampled_cbc[gene_name] = sample
            # Get subsampled outputs after subsampling cts_by_codon
            subsampled_outputs = proc.get_outputs(
                subsampled_cbc, cod_trunc_5p, cod_trunc_3p)
            # Predict cts on subsampled genes
            for rank, gene_name in enumerate(genes_by_density):
                # Skip genes with less dense data than benchmark gene
                if rank >= benchmark_gene_rank:
                    continue
                # Predict scaled counts on subsampled gene
                #     and compute correlation with real subsampled data
                num_cods = len(subsampled_cbc[gene_name])
                codon_set = {gene_name: [cod for cod in range(cod_trunc_5p, num_cods - cod_trunc_3p)]}
                X_te = proc.get_X(
                    codon_set, cds_dict, rel_cod_idxs=rel_cod_idxs,
                    rel_nt_idxs=rel_nt_idxs, rel_struc_idxs=rel_struc_idxs,
                    struc_dict=struc_dict).transpose()
                y_te = proc.get_y(codon_set, subsampled_outputs).transpose()
                y_te_hat = my_nn.pred_fn(X_te)
                pcorr = pearsonr(y_te, y_te_hat)[0][0]
                if benchmark_gene_rank == 1000:
                    pearson_r_subsample_1000.append(pcorr)
                elif benchmark_gene_rank == 2000:
                    pearson_r_subsample_2000.append(pcorr)
                elif benchmark_gene_rank == 3000:
                    pearson_r_subsample_3000.append(pcorr)
                elif benchmark_gene_rank == 4000:
                    pearson_r_subsample_4000.append(pcorr)
        # Fill in missing entries in subsampled corrs vectors
        num_genes = len(genes_by_density)
        subsample_vecs = [pearson_r_subsample_1000, 
            pearson_r_subsample_2000, 
            pearson_r_subsample_3000, 
            pearson_r_subsample_4000]
        for vec in subsample_vecs:
            while len(vec) < num_genes:
                vec.append("NA")
        
    # Open file to save model predictive performance data per gene
    with open(out_fname, "w") as out_file:
    # Record parameters used to process data for file
        out_file.write("#cod_trunc_5p: {0}\tcod_trunc_3p: {1}\n".format(
        cod_trunc_5p, cod_trunc_3p))
        out_file.write("#nn model: str_n17n15_cod_n5p4_nt_n15p14_rep0\n")
        # Make column names header
        out_file.write(
            "gene_name\tfp_density\ttr_gene\tmse\tpearson_r\t"
            "pearson_r_sub1000\tpearson_r_sub2000\tpearson_r_sub3000\n")
        # Record predictive performance for each gene in subsampled data set, 
        #    ordered by original fp density in genes
        for i, gene_name in enumerate(genes_by_density):
            format_string = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n"
            out_file.write(format_string.format(
                gene_name, trunc_fp_density_by_density[i], 
                training_gene_bools[i], mses_by_density[i], 
                pearson_r_by_density[i], pearson_r_subsample_1000[i], 
                pearson_r_subsample_2000[i], pearson_r_subsample_3000[i],
                pearson_r_subsample_4000[i]))

