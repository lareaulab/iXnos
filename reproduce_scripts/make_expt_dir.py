import sys
import iXnos.interface as inter

if __name__ == "__main__":
    expts_dir = sys.argv[1]
    expt_name = sys.argv[2]
    inter.make_expt_dirs(expts_dir, expt_name)

