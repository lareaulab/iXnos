#NOTE: 
import sys
import pickle
import numpy as np
import iXnos.interface as inter

if __name__ == "__main__":
    nn_dir = sys.argv[1]
    epoch = int(sys.argv[2])
    nn = inter.load_lasagne_feedforward_nn(nn_dir, epoch)

    init_data_fname = nn_dir + "/init_data/init_data.pkl"
    init_data = pickle.load(open(init_data_fname, "r"))
    rel_cod_idxs = init_data["rel_cod_idxs"]
    rel_nt_idxs = init_data["rel_nt_idxs"]

    cod_score_mat = inter.feat.get_lasagne_codon_featimp2(
        nn, rel_cod_idxs, rel_nt_idxs)

    epoch_dir = nn_dir + "/epoch{0}".format(epoch)
    np.savetxt(epoch_dir + "/codon_scores.tsv", cod_score_mat, delimiter="\t")
    inter.proc.pickle_obj(cod_score_mat, epoch_dir + "/codon_scores.pkl")

    model_name = init_data["name"]
    colormap_title = "{0} Epoch {1} Codon Scores".format(model_name, epoch)
    colormap_fname = epoch_dir + "/codon_scores_colormap.pdf"
    inter.feat.plot_2d_fi_colormap(
        cod_score_mat, rel_cod_idxs, colormap_title, colormap_fname)
