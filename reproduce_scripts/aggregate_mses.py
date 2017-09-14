import sys
import rp_predict.interface as inter

if __name__ == "__main__":
    procedures = ["aggregate", "aggregate_w_diffs"]
    #aggregate just collects MSEs for series and computes statistics
    #aggregate_w_diffs takes a reference series as well, and computes
    #   additional difference statistics
    procedure = sys.argv[1]
    assert procedure in procedures, \
        "procedure must be in [aggregate, aggregate_w_diffs]"
    expt_nn_dir = sys.argv[2]
    epoch = int(sys.argv[3])
    num_reps = int(sys.argv[4])
    out_fname = sys.argv[5]
    if procedure == "aggregate":
        series_names = sys.argv[6:]
    elif procedure == "aggregate_w_diffs":
        ref_series = sys.argv[6]
        series_names = sys.argv[7:]
        if ref_series not in series_names: 
            series_names += [ref_series]

    #Aggregate MSEs and stats into a file
    inter.feat.save_rep_series_MSEs(
        expt_nn_dir, series_names, epoch, num_reps, out_fname)
    if procedure == "aggregate_w_diffs":
        inter.feat.add_rep_series_MSE_diff_stats(
            expt_nn_dir, ref_series, epoch, num_reps, out_fname)

