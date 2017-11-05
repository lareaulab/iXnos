import iXnos.interface as inter
import sys

if __name__ == "__main__":
    model_names = ["nocod{0}_cod_n7p5_nt_n21p17".format(cod) 
                   for cod in range(-7, 6)]

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

    leaveout_idx = int(model_params[0][5:])
    #leaveout_str = model_params[0][5:]
    #sign = 1 if leaveout_str[0] == "p" else -1
    #leaveout_idx = sign * int(leaveout_str[1:])

    rel_cod_idxs = range(-7, leaveout_idx) + range(leaveout_idx + 1, 6)
    rel_nt_idxs = range(-21, 3*leaveout_idx) + range(3*(leaveout_idx + 1), 18)

    print model_name
    print rel_cod_idxs
    print rel_nt_idxs

    name = model_name + "_rep{0}".format(model_rep)
    neural_net = inter.make_lasagne_feedforward_nn(
        name, expt_dir, gene_seq_fname, gene_len_fname, tr_codons_fname,
        te_codons_fname, outputs_fname, rel_cod_idxs=rel_cod_idxs,
        rel_nt_idxs=rel_nt_idxs, lr_decay=lr_decay,
        nonlinearity="tanh", widths=[200], update_method="nesterov")
    neural_net.run_epochs(num_epochs)
    #skip plotting
