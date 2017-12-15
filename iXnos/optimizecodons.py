import os
import numpy as np
import copy
import random
import math

gfp_seq = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"

cit_seq = 'MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLGYGLMCFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK'

cit_dna = "ATGTCTAAAGGTGAAGAATTATTCACTGGTGTTGTCCCAATTTTGGTTGAATTAGATGGTGATGTTAATGGTCACAAATTTTCTGTCTCCGGTGAAGGTGAAGGTGATGCTACTTACGGTAAATTGACCTTAAAATTTATTTGTACTACTGGTAAATTGCCAGTTCCATGGCCAACCTTAGTCACTACTTTAGGTTATGGTTTGATGTGTTTTGCTAGATACCCAGATCATATGAAACAACATGACTTTTTCAAGTCTGCCATGCCAGAAGGTTATGTTCAAGAAAGAACTATTTTTTTCAAAGATGACGGTAACTACAAGACCAGAGCTGAAGTCAAGTTTGAAGGTGATACCTTAGTTAATAGAATCGAATTAAAAGGTATTGATTTTAAAGAAGATGGTAACATTTTAGGTCACAAATTGGAATACAACTATAACTCTCACAATGTTTACATCATGGCTGACAAACAAAAGAATGGTATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGTTCTGTTCAATTAGCTGACCATTATCAACAAAATACTCCAATTGGTGATGGTCCAGTCTTGTTACCAGACAACCATTACTTATCCTATCAATCTGCCTTATCCAAAGATCCAAACGAAAAGAGAGACCACATGGTCTTGTTAGAATTTGTTACTGCTGCTGGTATTACCCATGGTATGGATGAATTGTACAAATAA"

alpha="ACGT"
codons = [x+y+z for x in alpha for y in alpha for z in alpha]
cod2id = {codon:idx for idx, codon in enumerate(codons)}
id2cod = {idx:codon for codon, idx in cod2id.items()}
nt2id = {nt:idx for idx, nt in enumerate(alpha)}
id2nt = {idx:nt for nt, idx in nt2id.items()}

def get_letter_to_codons():
    let2cod = {}
    let2cod["F"] = ["TTT", "TTC"]
    let2cod["L"] = ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]
    let2cod["I"] = ["ATT", "ATC", "ATA"]
    let2cod["M"] = ["ATG"]
    let2cod["V"] = ["GTT", "GTC", "GTA", "GTG"]
    let2cod["S"] = ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"]
    let2cod["P"] = ["CCT", "CCC", "CCA", "CCG"]
    let2cod["T"] = ["ACT", "ACC", "ACA", "ACG"]
    let2cod["A"] = ["GCT", "GCC", "GCA", "GCG"]
    let2cod["Y"] = ["TAT", "TAC"]
    let2cod["H"] = ["CAT", "CAC"]
    let2cod["Q"] = ["CAA", "CAG"]
    let2cod["N"] = ["AAT", "AAC"]
    let2cod["K"] = ["AAA", "AAG"]
    let2cod["D"] = ["GAT", "GAC"]
    let2cod["E"] = ["GAA", "GAG"]
    let2cod["C"] = ["TGT", "TGC"]
    let2cod["W"] = ["TGG"]
    let2cod["R"] = ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]
    let2cod["G"] = ["GGT", "GGC", "GGA", "GGG"]
    let2cod["_"] = ["___"]
    return let2cod

def get_codon_to_letters():
    let2cod = get_letter_to_codons()
    cod2let = {codon:letter for letter in let2cod for codon in let2cod[letter]}
    return cod2let

def get_abbrev_to_codons():
    abb2cod = {}
    abb2cod["Phe"] = ["TTT", "TTC"]
    abb2cod["Leu"] = ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]
    abb2cod["Ile"] = ["ATT", "ATC", "ATA"]
    abb2cod["Met"] = ["ATG"]
    abb2cod["Val"] = ["GTT", "GTC", "GTA", "GTG"]
    abb2cod["Ser"] = ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"]
    abb2cod["Pro"] = ["CCT", "CCC", "CCA", "CCG"]
    abb2cod["Thr"] = ["ACT", "ACC", "ACA", "ACG"]
    abb2cod["Ala"] = ["GCT", "GCC", "GCA", "GCG"]
    abb2cod["Tyr"] = ["TAT", "TAC"]
    abb2cod["His"] = ["CAT", "CAC"]
    abb2cod["Gln"] = ["CAA", "CAG"]
    abb2cod["Asn"] = ["AAT", "AAC"]
    abb2cod["Lys"] = ["AAA", "AAG"]
    abb2cod["Asp"] = ["GAT", "GAC"]
    abb2cod["Glu"] = ["GAA", "GAG"]
    abb2cod["Cys"] = ["TGT", "TGC"]
    abb2cod["Trp"] = ["TGG"]
    abb2cod["Arg"] = ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]
    abb2cod["Gly"] = ["GGT", "GGC", "GGA", "GGG"]
    return abb2cod

def get_abbrev_to_letter():
    abb2let = {"Phe":"F", "Leu":"L", "Ile":"I", "Met":"M", "Val":"V", 
               "Ser":"S", "Pro":"P", "Thr":"T", "Ala":"A", "Tyr":"Y", 
               "His":"H", "Gln":"Q", "Asn":"N", "Lys":"K", "Asp":"D", 
               "Glu":"E", "Cys":"C", "Trp":"W", "Arg":"R", "Gly":"G"}
    return abb2let

def get_letter_to_abbrev():
    abb2let = get_abbrev_to_letter()
    let2abb = {value:key for key, value in abb2let.iteritems()}
    return let2abb

def get_cod_feat_seqs(aa_feat_seq, let2cod):
    if len(aa_feat_seq) == 0: 
        return []
    elif len(aa_feat_seq) == 1: 
        cod_feat_seqs = let2cod[aa_feat_seq]
        return cod_feat_seqs
    else: 
        first_aa = aa_feat_seq[0]
        prefixes = let2cod[first_aa]
        suffixes = get_cod_feat_seqs(aa_feat_seq[1:], let2cod)
        cod_feat_seqs = [pref + suff for pref in prefixes for suff in suffixes]
        return cod_feat_seqs

def get_cod_feat_vec(cod_feat_seq):
    if len(cod_feat_seq) % 3 != 0: 
        print "Error: codon sequence {0} length is not multiple of 3".format(
            cod_feat_seq) 
    num_cods = len(cod_feat_seq) / 3
    cod_list = [cod_feat_seq[i*3:(i+1)*3] for i in range(num_cods)]
    num_features = 64*len(cod_list)
    features = [0 for i in range(num_features)]
    for i, cod in enumerate(cod_list):
        if cod == "___": 
            pass
        else: 
            features[i*64 + cod2id[cod]] = 1
    return np.array(features).reshape(-1,1)
    
def get_nt_feat_vec(cod_feat_seq):
    if len(cod_feat_seq) % 3 != 0: 
        print "Error: codon sequence {0} length is not multiple of 3".format(
            cod_feat_seq) 
    num_cods = len(cod_feat_seq) / 3
    nt_list = list(cod_feat_seq)
    num_features = 4*len(nt_list)
    features = [0 for i in range(num_features)]
    for i, nt in enumerate(nt_list):
        if nt == "_": 
            pass
        else: 
            features[i*4 + nt2id[nt]] = 1
    return np.array(features).reshape(-1,1)
    
def get_cod_feat_matrix(cod_feat_seqs):
    return np.hstack([get_cod_feat_vec(cod_feat_seq) 
                      for cod_feat_seq in cod_feat_seqs]) 

def get_nt_feat_matrix(nt_feat_seqs):
    return np.hstack([get_nt_feat_vec(nt_feat_seq) 
                      for nt_feat_seq in nt_feat_seqs]) 

def get_pred_outputs_linreg(wts, cod_feat_matrix):
    return np.dot(cod_feat_matrix.transpose(), wts).reshape(1, -1)

def get_seqs_to_linreg_outputs(aa_seq, Asite_idx, wts, rel_cod_idxs):
    num_aa = len(aa_seq)
    let2cod = get_letter_to_codons()
    aa_idxs = [Asite_idx + elt for elt in rel_cod_idxs]
    filt_aa_idxs = filter(lambda x: x >= 0 and x < num_aa, aa_idxs)
    aa_feats = "".join([aa_seq[idx] if idx in filt_aa_idxs else "_" for idx in aa_idxs])
    cod_feat_seqs = get_cod_feat_seqs(aa_feats, let2cod)
    cod_feat_matrix = get_cod_feat_matrix(cod_feat_seqs)
    pred_outputs = get_pred_outputs_linreg(wts, cod_feat_matrix)
    seqs_to_outputs = {seq:out for seq, out in 
        zip(cod_feat_seqs, pred_outputs.ravel())}
    return seqs_to_outputs

def get_seqs_to_nn_outputs(aa_seq, Asite_idx, my_nn, rel_cod_idxs):
    num_aa = len(aa_seq)
    let2cod = get_letter_to_codons()
    aa_idxs = [Asite_idx + elt for elt in rel_cod_idxs]
    filt_aa_idxs = filter(lambda x: x >= 0 and x < num_aa, aa_idxs)
    filt_aa_idx_pos = [i for i in range(len(aa_idxs)) if aa_idxs[i] in filt_aa_idxs]
    wts = copy.deepcopy(my_nn.weights)
    wts_0 = wts[0]
    bias_idx = wts_0.shape[1] - 1
    filt_wts_pos = [i*64 + j for i in filt_aa_idx_pos for j in range(64)]
    filt_wts_pos.append(bias_idx)
    filt_wts_0 = wts_0.take(filt_wts_pos, axis=1)
    wts[0] = filt_wts_0
    aa_feats = "".join([aa_seq[idx] for idx in filt_aa_idxs])
    cod_feat_seqs = get_cod_feat_seqs(aa_feats, let2cod)
    cod_feat_matrix = get_cod_feat_matrix(cod_feat_seqs)
    pred_outputs = my_nn.forward(cod_feat_matrix, wts)
    seqs_to_outputs = {seq:out for seq, out in 
        zip(cod_feat_seqs, pred_outputs.ravel())}
    return seqs_to_outputs

def get_seqs_to_lasagne_outputs(
        aa_seq, Asite_idx, my_nn, rel_cod_idxs, nt_feats=False):
    """
        nt_feats: bool to include coterminous nt features with codon features
    """
    num_aa = len(aa_seq)
    let2cod = get_letter_to_codons()
    aa_idxs = [Asite_idx + elt for elt in rel_cod_idxs]
    filt_aa_idxs = filter(lambda x: x >= 0 and x < num_aa, aa_idxs)
    aa_feats = "".join([aa_seq[idx] if idx in filt_aa_idxs else "_" for idx in aa_idxs])
    cod_feat_seqs = get_cod_feat_seqs(aa_feats, let2cod)
    cod_feat_matrix = get_cod_feat_matrix(cod_feat_seqs)
    cod_feat_matrix = cod_feat_matrix.transpose()
    feat_matrix = cod_feat_matrix
    if nt_feats:
        nt_feat_matrix = get_nt_feat_matrix(cod_feat_seqs)
        nt_feat_matrix = nt_feat_matrix.transpose()
        feat_matrix = np.hstack([cod_feat_matrix, nt_feat_matrix])
    pred_outputs = my_nn.pred_fn(feat_matrix)
    seqs_to_outputs = {seq:out for seq, out in zip(cod_feat_seqs, pred_outputs.ravel())}
    return seqs_to_outputs

def get_suffs_to_opt_prefs(seqs_to_vit_scores, maximum=True):
    """
        Input: dictionary of sequence neighborhoods to viterbi scores at one 
        codon position

        Output: dictionary of n-1 codon suffixes to 1 codon prefix that 
        maximizes the viterbi score of the suffix
    """
    suffs_to_opt_prefs = {}
    cod2let = get_codon_to_letters()
    let2cod = get_letter_to_codons()
    some_seq = seqs_to_vit_scores.keys()[0]
    first_cod = some_seq[:3]
    first_aa = cod2let[first_cod]
    pref_cods = let2cod[first_aa]
    for seq in seqs_to_vit_scores:
        pref = seq[:3]
        suff = seq[3:]
        if not suffs_to_opt_prefs.get(suff, False):
            suffs_to_opt_prefs[suff] = pref
        else:
            curr_opt_pref = suffs_to_opt_prefs[suff]
            curr_opt_seq = curr_opt_pref + suff
            if maximum == True:
                if seqs_to_vit_scores[seq] > seqs_to_vit_scores[curr_opt_seq]:
                    suffs_to_opt_prefs[suff] = pref
            else:
                if seqs_to_vit_scores[seq] < seqs_to_vit_scores[curr_opt_seq]:
                    suffs_to_opt_prefs[suff] = pref
    return suffs_to_opt_prefs

def get_seqs_to_vit_scores(seqs_to_outputs, prev_seqs_to_vit_scores, prev_suff_to_opt_pref, incl_start=False, past_end=False):
    seqs_to_vit_scores = {}
    for seq in seqs_to_outputs:
        prev_suff = seq
        #if not past_end:
        prev_suff = seq[:-3]
        #if incl_start: 
        #    prev_suff = prev_suff[3:]
        prev_pref = prev_suff_to_opt_pref[prev_suff]
        prev_seq = prev_pref + prev_suff
        prev_vit_score = prev_seqs_to_vit_scores[prev_seq]
        vit_score = prev_vit_score + seqs_to_outputs[seq]
        seqs_to_vit_scores[seq] = vit_score
    return seqs_to_vit_scores

def viterbi_backtrack(
        seqs_to_vit_scores_by_pos, suff_to_opt_pref_by_pos, rel_cod_idxs, 
        maximum=True):
    num_aa = len(seqs_to_vit_scores_by_pos)
    if maximum == True:
        opt_end_seq = max(seqs_to_vit_scores_by_pos[-1].keys(), 
                          key=lambda seq: seqs_to_vit_scores_by_pos[-1][seq])
    else:
        opt_end_seq = min(seqs_to_vit_scores_by_pos[-1].keys(), 
                          key=lambda seq: seqs_to_vit_scores_by_pos[-1][seq])
    opt_vit_score = seqs_to_vit_scores_by_pos[-1][opt_end_seq]
    opt_total_seq = opt_end_seq
    curr_idx = len(seqs_to_vit_scores_by_pos) - 1
    curr_opt_seq = opt_end_seq
    while curr_idx > 0:
        prev_opt_suff = curr_opt_seq
        prev_opt_suff = prev_opt_suff[:-3]
        if curr_idx + min(rel_cod_idxs) > 0:
            try:
                prev_opt_pref = suff_to_opt_pref_by_pos[curr_idx-1][prev_opt_suff]
            except KeyError:
                print suff_to_opt_pref_by_pos[curr_idx-1]
        else:
            prev_opt_pref = ""
        #print prev_opt_pref
        opt_total_seq = prev_opt_pref + opt_total_seq
        #print opt_total_seq
        prev_opt_seq = prev_opt_pref + prev_opt_suff
        curr_opt_seq = prev_opt_seq
        curr_idx -= 1
    return opt_total_seq, opt_vit_score

def get_seqs_to_linreg_outputs_by_pos(aa_seq, wts, rel_cod_idxs):
    seqs_to_outputs_by_pos = []
    num_aa = len(aa_seq) 
    for i in range(num_aa):
        seqs_to_outputs = get_seqs_to_linreg_outputs(
                              aa_seq, i, wts, rel_cod_idxs)
        seqs_to_outputs_by_pos.append(seqs_to_outputs)
    return seqs_to_outputs_by_pos

def get_seqs_to_nn_outputs_by_pos(aa_seq, nn, rel_cod_idxs):
    seqs_to_outputs_by_pos = []
    num_aa = len(aa_seq) 
    for i in range(num_aa):
        seqs_to_outputs = get_seqs_to_nn_outputs(
                              aa_seq, i, nn, rel_cod_idxs)
        seqs_to_outputs_by_pos.append(seqs_to_outputs)
    return seqs_to_outputs_by_pos

def get_seqs_to_lasagne_outputs_by_pos(
        aa_seq, nn, rel_cod_idxs, nt_feats=False):
    """
        Inputs: 
            rel_cod_idxs: a list of consecutive codon indices for seq nbhood
            nt_feats: flag for seq nb that includes redundant nt. feats
                NOTE: seq nb in model for nts must cover same region as codons
        Returns a list of length len(aa_seq)
            Each entry is a dictionary mapping all synonymous sequence
            neighborhoods around an a site to their score under nn
    """
    seqs_to_outputs_by_pos = []
    num_aa = len(aa_seq) 
    for i in range(num_aa):
        seqs_to_outputs = get_seqs_to_lasagne_outputs(
            aa_seq, i, nn, rel_cod_idxs, nt_feats=nt_feats)
        seqs_to_outputs_by_pos.append(seqs_to_outputs)
    return seqs_to_outputs_by_pos

def do_viterbi_propagate(
        num_aa, seqs_to_outputs_by_pos, rel_cod_idxs, maximum=True):
    seqs_to_vit_scores_by_pos = []
    suffs_to_opt_prefs_by_pos = []
    for i in range(num_aa):
        if i == 0: 
            seqs_to_vit_scores = seqs_to_outputs_by_pos[0]
            suffs_to_opt_prefs = get_suffs_to_opt_prefs(seqs_to_vit_scores, 
                                   maximum=maximum)
        else:
            incl_start = True if min(rel_cod_idxs) + i <= 0 else False
            past_end = True if max(rel_cod_idxs) + i > num_aa - 1 else False
            seqs_to_vit_scores = get_seqs_to_vit_scores(
                seqs_to_outputs_by_pos[i], seqs_to_vit_scores_by_pos[i-1], 
                suffs_to_opt_prefs_by_pos[i-1], incl_start=incl_start, 
                past_end=past_end)
            suffs_to_opt_prefs = get_suffs_to_opt_prefs(
                seqs_to_vit_scores, maximum=maximum)
        seqs_to_vit_scores_by_pos.append(seqs_to_vit_scores)
        suffs_to_opt_prefs_by_pos.append(suffs_to_opt_prefs)
    return seqs_to_vit_scores_by_pos, suffs_to_opt_prefs_by_pos

def get_optimal_codons_linreg(aa_seq, wts, rel_cod_idxs, maximum=True):
    #NOTE: Expects rel_cod_idxs to be a consecutive series of integers
    #NOTE: Expects wts from linreg trained on only rel_cod_idxs as features
    num_aa = len(aa_seq)
    seqs_to_outputs_by_pos = get_seqs_to_linreg_outputs_by_pos(
        aa_seq, wts, rel_cod_idxs)
    seqs_to_vit_scores_by_pos, suff_to_opt_pref_by_pos = do_viterbi_propagate(
        num_aa, seqs_to_outputs_by_pos, rel_cod_idxs, maximum=maximum)
    #print seqs_to_vit_scores_by_pos
    opt_total_seq, opt_vit_score = viterbi_backtrack(
        seqs_to_vit_scores_by_pos, suff_to_opt_pref_by_pos, rel_cod_idxs, 
        maximum=maximum)
    return opt_total_seq, opt_vit_score

def get_optimal_codons_nn(aa_seq, nn, rel_cod_idxs, maximum=True):
    #NOTE: Expects rel_cod_idxs to be a consecutive series of integers
    #NOTE: Expects nn trained on only rel_cod_idxs as features
    num_aa = len(aa_seq)
    seqs_to_outputs_by_pos = get_seqs_to_nn_outputs_by_pos(
        aa_seq, nn, rel_cod_idxs)
    seqs_to_vit_scores_by_pos, suff_to_opt_pref_by_pos = do_viterbi_propagate(
        num_aa, seqs_to_outputs_by_pos, rel_cod_idxs, maximum=maximum)
    #print seqs_to_vit_scores_by_pos
    opt_total_seq, opt_vit_score = viterbi_backtrack(seqs_to_vit_scores_by_pos, suff_to_opt_pref_by_pos, rel_cod_idxs, maximum=maximum)
    return opt_total_seq, opt_vit_score

def get_optimal_codons_lasagne(
        aa_seq, nn, rel_cod_idxs, maximum=True, nt_feats=False, unlog=False):
    #NOTE: Expects rel_cod_idxs to be a consecutive series of integers
    #NOTE: Expects nn trained on only rel_cod_idxs as features
    num_aa = len(aa_seq)
    seqs_to_outputs_by_pos = get_seqs_to_lasagne_outputs_by_pos(
        aa_seq, nn, rel_cod_idxs, nt_feats=nt_feats)
    #if unlog:
    #    for pos_dict in seqs_to_outputs_by_pos:
    #        for seq in pos_dict:
    #            pos_dict[seq] = math.exp(pos_dict[seq])
    seqs_to_vit_scores_by_pos, suff_to_opt_pref_by_pos = do_viterbi_propagate(
        num_aa, seqs_to_outputs_by_pos, rel_cod_idxs, maximum=maximum)
    opt_total_seq, opt_vit_score = viterbi_backtrack(
        seqs_to_vit_scores_by_pos, suff_to_opt_pref_by_pos, 
        rel_cod_idxs, maximum=maximum)
    return opt_total_seq.strip("_"), opt_vit_score

def get_codons_by_aa(aa_seq, cod_seq):
    aas = "ACDEFGHIKLMNPQRSTVWY"
    codons_by_aa = {aa:[] for aa in aas}
    for i in range(len(aa_seq)):
        aa = aa_seq[i]
        cod = cod_seq[3*i : 3*(i+1)]
        codons_by_aa[aa].append(cod)
    return codons_by_aa

def get_random_cod_seq(aa_seq):
    let2cod = get_letter_to_codons()
    codons = [random.choice(let2cod[aa]) for aa in aa_seq]
    cod_seq = ''.join(codons)
    return cod_seq

def score_cod_seq_full(
        cod_seq, my_nn, rel_cod_idxs, rel_nt_idxs=False, rel_struc_idxs=False, 
        struc_dict=False, unlog=False):
    """
        Scores a coding sequence by getting a predicted score at each codon
        in the CDS, and inputs sequence outside of the CDS as zeroes
    """
    # Get number of codons in coding sequence
    num_cod = len(cod_seq)/3
    # Split coding sequence into list of codons
    cod_seq_split = [cod_seq[i:i+3] for i in range(0, len(cod_seq), 3)]
    cod_feat_seqs_split =\
        [[cod_seq_split[i + rel_idx] if 0 <= i + rel_idx < num_cod else "___" 
            for rel_idx in rel_cod_idxs] for i in range(num_cod)]
    cod_feat_seqs = ["".join(cods) for cods in cod_feat_seqs_split]
    cod_feat_matrix = get_cod_feat_matrix(cod_feat_seqs)
    cod_feat_matrix = cod_feat_matrix.transpose()
    feat_matrix = cod_feat_matrix
    if rel_nt_idxs and rel_nt_idxs != []:
        nt_feat_matrix = get_nt_feat_matrix(cod_feat_seqs)
        nt_feat_matrix = nt_feat_matrix.transpose()
        feat_matrix = np.hstack([cod_feat_matrix, nt_feat_matrix])
    outputs = my_nn.pred_fn(feat_matrix)
    if unlog:
        return np.exp(outputs).sum()
    return outputs.sum()

def score_cod_seq_trunc(
        cds, my_nn, rel_cod_idxs, rel_nt_idxs=False, unlog=False):
    """
        Scores a coding sequence by getting a predicted score at each codon
        in the CDS where the sequence neighborhood is contained in the CDS
    """
    #if bool(rel_struc_idxs) != bool(struc_dict):
    #    print "Error: must set both or neither of rel_struc_idxs and struc_dict"
    #    raise ValueError
    # Get number of codons in coding sequence
    num_cod = len(cds)/3
    num_nt = len(cds)
    # Get codon indices to score in CDS
    cds_valid_idxs = range(num_cod)
    # Filter invalid codon neighborhoods
    def is_valid_cod_nb(cds, Asite_idx, rel_cod_idxs):
        num_cod = len(cds)/3
        is_valid = (0 <= Asite_idx + min(rel_cod_idxs) and\
                    Asite_idx + max(rel_cod_idxs) < num_cod)
        return is_valid
    cds_valid_idxs = filter(
        lambda x: is_valid_cod_nb(cds, x, rel_cod_idxs), cds_valid_idxs)
    # If we use nt features, filter invalid nt neighborhoods
    if rel_nt_idxs:
        def is_valid_nt_nb(cds, Asite_idx, rel_nt_idxs):
            num_nts = len(cds)
            is_valid = (0 <= Asite_idx*3 + min(rel_nt_idxs) and\
                        Asite_idx*3 + max(rel_nt_idxs) < num_nt)
            return is_valid
        cds_valid_idxs = filter(
            lambda x: is_valid_nt_nb(cds, x, rel_nt_idxs), cds_valid_idxs)
    # Split coding sequence into list of codons
    cds_cod_list = [cds[i:i+3] for i in range(0, len(cds), 3)]
    # Get sequences of codon neighborhoods at each valid codon position
    cod_nbs = []
    for i in cds_valid_idxs:
        cods = [cds_cod_list[i + rel_idx] for rel_idx in rel_cod_idxs]
        cod_nb = "".join(cods)
        cod_nbs.append(cod_nb)
    # Get codon feature matrix, initialize feature matrix
    cod_feat_matrix = get_cod_feat_matrix(cod_nbs)
    cod_feat_matrix = cod_feat_matrix.transpose()
    feat_matrix = cod_feat_matrix
    # If we use nt features, 
    #   get sequences of nt neighborhoods at each valid codon position
    if rel_nt_idxs:
        nt_nbs = []
        for i in cds_valid_idxs:
            nts = [cds[i + rel_idx] for rel_idx in rel_nt_idxs]
            nt_nb = "".join(nts)
            nt_nbs.append(nt_nb)
    #   Get nt feature matrix, append to feature matrix
        nt_feat_matrix = get_nt_feat_matrix(nt_nbs)
        nt_feat_matrix = nt_feat_matrix.transpose()
        feat_matrix = np.hstack([cod_feat_matrix, nt_feat_matrix])
    outputs = my_nn.pred_fn(feat_matrix)
    if unlog:
        return np.exp(outputs).sum()
    return outputs.sum()

def get_random_cod_seqs(aa_seq, num_samples):
    return [get_random_cod_seq(aa_seq) for i in range(num_samples)]

def get_score_dist(aa_seq, my_nn, rel_cod_idxs, num_samples, nt_feats=False):
    if nt_feats:
        rel_nt_idxs = [cod_idx * 3 + i for cod_idx in rel_cod_idxs for i in range(3)]
    else:
        rel_nt_idxs = False
    cod_seqs = get_random_cod_seqs(aa_seq, num_samples)
    print "computing score dist"
    scores = [score_cod_seq_full(cod_seq, my_nn, rel_cod_idxs, rel_nt_idxs) 
        for cod_seq in cod_seqs]
    print "done computing score dist"
    idxs = sorted(range(len(scores)), key=lambda x: scores[x])
    scores_sorted = [scores[idx] for idx in idxs]
    cod_seqs_sorted = [cod_seqs[idx] for idx in idxs]
    return scores_sorted, cod_seqs_sorted

def make_opt_dir(opt_dir):
    if not os.path.exists(opt_dir):
        os.makedirs(opt_dir)
