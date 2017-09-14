import rp_predict.interface as inter
import sys

if __name__ == "__main__":
    #Load parameters
    nn_dir = sys.argv[1]
    epoch = int(sys.argv[2])
    te_data_table_fname = sys.argv[3]
    gene2sym_fname = sys.argv[4]
    plot_dir = sys.argv[5]
    #Load nn
    my_nn = inter.load_lasagne_feedforward_nn(nn_dir, epoch)
    #Plot predictions vs. true scaled cts for each gene
    inter.plot_nn_preds_by_gene(
        my_nn, te_data_table_fname, gene2sym_fname, plot_dir,
        gene2sym_mod_fn=False)
