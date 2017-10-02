import iXnos.interface as inter
import sys

if __name__ == "__main__":
    model_names = [
        "cod_p0", "cod_p0_nt_p0p2", "cod_n3p2", "cod_n3p2_nt_n9p8", 
        "cod_n5p4", "cod_n5p4_nt_n15p14", "cod_n7p5", "cod_n7p5_nt_n21p17"]
    model_name = sys.argv[1]
    expt_dir = sys.argv[2]
    gene_len_fname = sys.argv[3]
    gene_seq_fname = sys.argv[4]
    tr_codons_fname = sys.argv[5]
    te_codons_fname = sys.argv[6]
    outputs_fname = sys.argv[7]

    assert model_name in model_names, "model name {0} not in " \
        "list of models".format(model_name)

    #Parse codon positions, assumes 1 continuous codon neighborhood
    model_params = model_name.split("_")
    cod_idxs = model_params[1]
    cod_start_sign = 1 if cod_idxs[0] == "p" else -1
    cod_start_idx = cod_start_sign * int(cod_idxs[1])
    cod_end_sign = 1 if cod_idxs[-2] == "p" else -1
    cod_end_idx = cod_end_sign * int(cod_idxs[-1])
    rel_cod_idxs = range(cod_start_idx, cod_end_idx + 1)

    #Parse nt positions, assumes 1 continuous nt neighborhood
    if len(model_params) > 2:
        #Get start sign
        nt_idxs = model_params[3]
        nt_start_sign = 1 if nt_idxs[0] == "p" else -1
        nt_idxs = nt_idxs[1:]
        nt_idxs = nt_idxs.split("n")
        if len(nt_idxs) == 2:
            nt_end_sign = -1
        else:
            nt_idxs = nt_idxs[0]
            nt_idxs = nt_idxs.split("p")
            nt_end_sign = 1
        nt_start_idx = nt_start_sign * int(nt_idxs[0])
        nt_end_idx = nt_end_sign * int(nt_idxs[1])
        rel_nt_idxs = range(nt_start_idx, nt_end_idx + 1)
    else:
        rel_nt_idxs = []

    print model_name
    print rel_cod_idxs
    print rel_nt_idxs

    name = "lr_" + model_name
    wts, y_tr_hat, y_te_hat, y_tr, y_te = inter.make_linreg(
        expt_dir, name, gene_seq_fname, gene_len_fname, tr_codons_fname, 
        te_codons_fname, outputs_fname, rel_cod_idxs=rel_cod_idxs, 
        rel_nt_idxs=rel_nt_idxs)
    #skip plotting
