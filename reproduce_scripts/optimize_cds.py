import sys
import iXnos.interface as inter

if __name__ == "__main__":
    nn_dir = sys.argv[1]
    epoch = int(sys.argv[2])
    nt_feats = bool(sys.argv[3])
    aa_seq = sys.argv[4]
    out_fname = sys.argv[5]
    num_samples = 100000

    min_seq, min_score = inter.get_lasagne_optimal_codons(
        nn_dir, epoch, aa_seq, maximum=False, nt_feats=nt_feats)
    rand_scores, rand_seqs = inter.get_protein_score_dist(
        nn_dir, epoch, aa_seq, num_samples, nt_feats=nt_feats)
    max_seq, max_score = inter.get_lasagne_optimal_codons(
        nn_dir, epoch, aa_seq, maximum=True, nt_feats=nt_feats)

    out_file = open(out_fname, "w")
    model_name = nn_dir.split("/")[-1]
    header = "CDS series under model {0}\n".format(model_name)
    header += "Protein sequence {1}\n".format(aa_seq)
    out_file.write(header)
    out_file.write(">Minimum score: {0}\n".format(min_score))
    out_file.write(min_seq + "\n")
    for i in [0, 33333, 66666, 99999]:
        out_file.write(">{0} score: {1}\n".format(i, rand_scores[i]))
        out_file.write(rand_seqs[i] + "\n")
    out_file.write(">Maximum score: {0}\n".format(max_score))
    out_file.write(max_seq + "\n")
    out_file.close()
