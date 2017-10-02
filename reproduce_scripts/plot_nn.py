import iXnos.interface as inter
import sys

if __name__ == "__main__":
    expt_dir = sys.argv[1]
    name = sys.argv[2]
    epoch = int(sys.argv[3])

    scat_xlim = (-1, 10)
    scat_ylim = (-1, 10)
    scat_textpos = (1, 8.5)

    inter.plot_lasagne_nn(
        expt_dir, name, epoch, scat_xlim, scat_ylim, scat_textpos)
