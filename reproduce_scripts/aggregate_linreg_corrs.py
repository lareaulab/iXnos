import sys
import iXnos.interface as inter

if __name__ == "__main__":
    method = sys.argv[1]
    expt_linreg_dir = sys.argv[2]
    out_fname = sys.argv[3]
    model_names = sys.argv[4:]
    print
    print model_names
    print

    #Aggregate MSEs and stats into a file
    inter.linreg.aggregate_linreg_corrs(expt_linreg_dir, model_names, out_fname)
