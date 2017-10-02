import sys
import iXnos.interface as inter

if __name__ == "__main__":
    expt_linreg_dir = sys.argv[1]
    out_fname = sys.argv[2]
    series_names = sys.argv[3:]
    print
    print series_names
    print

    #Aggregate MSEs and stats into a file
    inter.linreg.aggregate_linreg_MSEs(expt_linreg_dir, series_names, out_fname)
