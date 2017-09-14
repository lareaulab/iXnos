import sys
import numpy as np
import pickle

if __name__ == "__main__":
    in_fname = sys.argv[1]
    out_fname = sys.argv[2]

    in_file = open(in_fname, "r")
    matrix = pickle.load(in_file)
    np.savetxt(out_fname, matrix, delimiter="\t")
    in_file.close()
