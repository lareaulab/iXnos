import iXnos.interface as inter
import shutil
import numpy as np
import sys

if __name__ == "__main__":
    # NOTE: If changing these, change idx parameters below
    model_names = ["n17n15_cod_n5p4_nt_n15p14", "p13p42_cod_n5p4_nt_n15p14",]
    model_name = sys.argv[1]
    model_rep = sys.argv[2]
    expt_dir = sys.argv[3]
    sam_fname = sys.argv[4]
    gene_len_fname = sys.argv[5]
    gene_seq_fname = sys.argv[6]
    struc_fname = sys.argv[7]
        #project_dir + "/structure_data/yeast_13cds10.windows.30len.fold"
    tr_codons_fname = sys.argv[8]
    te_codons_fname = sys.argv[9]
    outputs_fname = sys.argv[10]
    num_epochs = int(sys.argv[11])
    lr_decay = float(sys.argv[12])

    assert model_name in model_names, "model name {0} not in " \
        "list of models".format(model_name)

    rel_cod_idxs = range(-5, 5)
    rel_nt_idxs = range(-15, 15)
    learning_rate = 0.005

    print model_name
    print rel_cod_idxs
    print rel_nt_idxs

    if model_name == "n17n15_cod_n5p4_nt_n15p14":
        rel_struc_idxs = range(-17, -14)
        name = "str_" + model_name + "_rep{0}".format(model_rep)
        success = False
        while not success:
            my_nn = inter.make_lasagne_feedforward_nn(
                name, expt_dir, gene_seq_fname, gene_len_fname, tr_codons_fname,
                te_codons_fname, outputs_fname, rel_cod_idxs=rel_cod_idxs,
                rel_nt_idxs=rel_nt_idxs, rel_struc_idxs=rel_struc_idxs,
                struc_fname=struc_fname, lr_decay=lr_decay,
                nonlinearity="tanh", widths=[200], update_method="nesterov")
            failed = False
            for i in range(num_epochs + 1):
                my_nn.run_epochs(i)
                if np.isnan(my_nn.test_err_by_epoch[-1]):
                    failed = True
                    shutil.rmtree(my_nn.out_dir + "/" + my_nn.name)
                    break
            if not failed:
                success = True
    
    if model_name == "p13p42_cod_n5p4_nt_n15p14":
        max_struc_start_idx = 13
        max_struc_width = 30
        name = "max_str_" + model_name + "_rep{0}".format(model_rep)
        success = False
        while not success:
            learning_rate = 0.005
            my_nn = inter.make_lasagne_feedforward_nn(
                name, expt_dir, gene_seq_fname, gene_len_fname, tr_codons_fname,
                te_codons_fname, outputs_fname, rel_cod_idxs=rel_cod_idxs,
                rel_nt_idxs=rel_nt_idxs, lr_decay=lr_decay,
                max_struc_start_idx=max_struc_start_idx,
                struc_fname=struc_fname, max_struc_width=max_struc_width,
                nonlinearity="tanh", widths=[200], update_method="nesterov")
            failed = False
            for i in range(num_epochs + 1):
                my_nn.run_epochs(i)
                if np.isnan(my_nn.test_err_by_epoch[-1]):
                    failed = True
                    shutil.rmtree(my_nn.out_dir + "/" + my_nn.name)
                    break
            if not failed:
                success = True
