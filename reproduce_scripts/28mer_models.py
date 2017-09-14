import rp_predict.interface as inter
import sys

if __name__ == "__main__":
    model_names = ["s28_cod_n7p5_nt_n21p17"]

    model_name = sys.argv[1]
    expt_dir = sys.argv[2]
    sam_fname = sys.argv[3]
    gene_len_fname = sys.argv[4]
    gene_seq_fname = sys.argv[5]
    tr_codons_fname = sys.argv[6]
    te_codons_fname = sys.argv[7]
    outputs_fname = sys.argv[8]
    num_epochs = int(sys.argv[9])

    assert model_name in model_names, "model name {0} not in " \
        "list of models".format(model_name)

    rel_cod_idxs = range(-7, 6)
    rel_nt_idxs = range(-21, 18)

    print model_name
    print rel_cod_idxs
    print rel_nt_idxs

    name = model_name 
    neural_net = inter.make_lasagne_feedforward_nn(
        name, expt_dir, gene_seq_fname, gene_len_fname, tr_codons_fname,
        te_codons_fname, outputs_fname, rel_cod_idxs=rel_cod_idxs,
        rel_nt_idxs=rel_nt_idxs,
        nonlinearity="tanh", widths=[200], update_method="nesterov")
    neural_net.run_epochs(num_epochs)
    #skip plotting
