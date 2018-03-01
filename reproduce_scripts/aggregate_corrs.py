import sys
import iXnos.interface as inter

if __name__ == "__main__":
    methods = ["pearson", "spearman"]
    method = sys.argv[1]
    assert method in methods, \
        "method must be in [pearson, spearman]"
    expt_nn_dir = sys.argv[2]
    epoch = int(sys.argv[3])
    num_reps = int(sys.argv[4])
    out_fname = sys.argv[5]
    series_names = sys.argv[6:]

    #Aggregate corrs into a file
    inter.feat.save_rep_series_corrs(
        expt_nn_dir, series_names, epoch, num_reps, out_fname, 
        method=method)
