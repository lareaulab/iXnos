import iXnos.interface as inter
import sys

if __name__ == "__main__":
    model_names = [
        "cod_p0", "cod_p0_nt_p0p2", "cod_n3p2", "cod_n3p2_nt_n9p8", 
        "cod_n5p4", "cod_n5p4_nt_n15p14", "cod_n7p5", "cod_n7p5_nt_n21p17"]
    model_name = sys.argv[1]
    model_rep = sys.argv[2]
    expt_dir = sys.argv[3]
    sam_fname = sys.argv[4]
    gene_len_fname = sys.argv[5]
    gene_seq_fname = sys.argv[6]
    tr_codons_fname = sys.argv[7]
    te_codons_fname = sys.argv[8]
    outputs_fname = sys.argv[9]
    num_epochs = int(sys.argv[10])
    lr_decay = float(sys.argv[11])

    assert model_name in model_names, "model name {0} not in " \
        "list of models".format(model_name)

    model_params = model_name.split("_")
    cod_idxs = model_params[1]
    cod_start_sign = 1 if cod_idxs[0] == "p" else -1
    cod_start_idx = cod_start_sign * int(cod_idxs[1])
    cod_end_sign = 1 if cod_idxs[-2] == "p" else -1
    cod_end_idx = cod_end_sign * int(cod_idxs[-1])
    rel_cod_idxs = range(cod_start_idx, cod_end_idx + 1)

    #NOTE: Still needs testing on next rerun - check if this works! 
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

    name = "full_" + model_name + "_rep{0}".format(model_rep)
    neural_net = inter.make_lasagne_feedforward_nn(
        name, expt_dir, gene_seq_fname, gene_len_fname, tr_codons_fname,
        te_codons_fname, outputs_fname, rel_cod_idxs=rel_cod_idxs,
        rel_nt_idxs=rel_nt_idxs, lr_decay=lr_decay,
        nonlinearity="tanh", widths=[200], update_method="nesterov")
    neural_net.run_epochs(num_epochs)
    #skip plotting
