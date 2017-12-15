import numpy as np
import math
import random
import copy
import matplotlib.pyplot as plt
import pickle
import os
import imp
import glob

alpha="ACGT"
codons = [x+y+z for x in alpha for y in alpha for z in alpha]
cod2id = {codon:idx for idx, codon in enumerate(codons)}
id2cod = {idx:codon for codon, idx in cod2id.items()}
nt2id = {nt:idx for idx, nt in enumerate(alpha)}
id2nt = {idx:nt for nt, idx in nt2id.items()}

aas = list("ACDEFGHIKLMNPQRSTVWY")
cod2aa = {
        'AAA':'K', 'AAG':'K', 'AAC':'N', 'AAT':'N', 'ACA':'T', 'ACC':'T', 
        'ACG':'T', 'ACT':'T', 'AGA':'R', 'AGG':'R', 'CGA':'R', 'CGC':'R', 
        'CGG':'R', 'CGT':'R', 'AGC':'S', 'AGT':'S', 'TCA':'S', 'TCC':'S', 
        'TCG':'S', 'TCT':'S', 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'CAA':'Q', 'CAG':'Q', 'CAC':'H', 'CAT':'H', 'CCA':'P', 'CCC':'P', 
        'CCG':'P', 'CCT':'P', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'TTA':'L', 'TTG':'L', 'TTC':'F', 'TTT':'F', 'GAA':'E', 'GAG':'E',
        'GAC':'D', 'GAT':'D', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'GTA':'V', 'GTC':'V',
        'GTG':'V', 'GTT':'V', 'TAA':'STOP', 'TAG':'STOP', 'TGA':'STOP',
        'TAC':'Y', 'TAT':'Y', 'TGG':'W', 'TGC':'C', 'TGT':'C'}
aa2id = {aa:idx for idx, aa in enumerate(aas)}
aa2id["STOP"] = 20
id2aa = {idx:aa for aa, idx in aa2id.items()}

def sam_add_simple_map_wts(in_fname, out_fname):
    """
     in_fname: name of input sam file
     out_fname: name of output sam file, identical to original with extra 
                terminal field for number of mappings of read
     NOTE: Expects multiple mappings of a read to be on consecutive lines
    """
    in_file = open(in_fname, "r")
    out_file = open(out_fname, "w")
    first = True
    prev_read = "thisisnotareadname"
    read_ct = 0
    read_maps = []
    for line in in_file:
        if line[0] == "@":
            out_file.write(line)
            continue
        line = line.strip().split()
        read_name = line[0]
        #Sam positions are 1 based, make 0 based for mod arithmetic
        if read_name != prev_read:
            if not first:
                simple_wt = float(1) / read_ct
                for read_map in read_maps:
                    out_file.write(read_map + "\t{0}\n".format(simple_wt))
                read_ct = 0
                read_maps = []
        read_maps.append("\t".join(line))
        read_ct += 1
        prev_read = read_name
        first = False
    #Get last read
    for read_map in read_maps:
        out_file.write(read_map + "\t{0}\n".format(read_ct))
    in_file.close()
    out_file.close()

def sam_flag_set(bit, flag):
    #NOTE: use 0-indexed bits, starting from lowest contribution to flag
    bin_flag = bin(int(flag)).zfill(12)
    return int(bin_flag[11 - bit])
    
def sam_filter_unmapped(in_fname, out_fname):
    """
     in_fname: name of input sam file
     out_fname: name of output sam file, identical to original with
                unmapped lines removed from file
    """
    in_file = open(in_fname, "r")
    out_file = open(out_fname, "w")
    for line in in_file:
        if line[0] == "@":
            out_file.write(line)
            continue
        line = line.strip().split()
        flag = int(line[1])
        if sam_flag_set(2, flag):
            continue
        else: 
            out_file.write("\t".join(line) + "\n")
    in_file.close()
    out_file.close()

def RSEM_add_weights(in_sam_fname, out_sam_fname):
    """
     in_sam_fname: name of input RSEM sam file
     out_sam_fname: name of output sam file, identical to original with extra 
                terminal field for each line equal to RSEM mapping weight
    """
    #NOTE: Assumes each non-header line ends with map weight field of format ZW:f:{weight}
    in_file = open(in_sam_fname, "r")
    out_file = open(out_sam_fname, "w")
    for line in in_file:
        if line[0] == "@":
            out_file.write(line)
            continue
        line = line.strip().split()
        weight = line[-1][5:]
        line.append(weight)
        out_file.write("\t".join(line) + "\n")
    in_file.close()
    out_file.close()

def sam_subsample(in_sam_fname, out_sam_fname, pct):
    in_file = open(in_sam_fname, "r")
    out_file = open(out_sam_fname, "w")
    line_buffer = []
    current_line = ''
    first_line = True
    for line in in_file:
        if line[0] == "@":
            out_file.write(line)
            continue
        line = line.strip().split()
        if line[0] != current_line: 
            if not first_line:
                randunif = random.uniform(0, 1)
                if randunif < pct:
                    for saved_line in line_buffer:
                        out_file.write("\t".join(saved_line) + "\n")
                line_buffer = []           
                current_line = line[0]
            else: 
                current_line = line[0] 
                first_line = False
            line_buffer.append(line)
        else: 
            line_buffer.append(line)
    in_file.close()
    out_file.close()

def get_cts_by_size_and_frame(sam_fname, len_dict, verbose=False):
    """
     sam_fname: name of input sam file
     verbose: bool, set to True to output fp size range of sam file
    """
    cts_by_size_and_frame = {}
    in_file = open(sam_fname, "r")
    for line in in_file:
        if line[0] == "@":
            continue
        line = line.strip().split()
        gene = line[2]
        utr5_len = len_dict[gene]["utr5"]
        cds_len = len_dict[gene]["cds"]
        map_pos = int(line[3]) - 1 - utr5_len
        if map_pos < -15 or map_pos > cds_len - 9:
            continue
        frame = map_pos % 3
        read = line[9]
        size = len(read)
        ct = float(line[-1])
        if not cts_by_size_and_frame.get(size, False):
            cts_by_size_and_frame[size] = [0, 0, 0]
        cts_by_size_and_frame[size][frame] += ct
    in_file.close()
    #make [0, 0, 0] entries for intermediate sizes not in dict
    min_size = min(cts_by_size_and_frame.keys())
    max_size = max(cts_by_size_and_frame.keys())
    if verbose:
        print "Minimum fp size: {0}".format(min_size)
        print "Maximum fp size: {0}".format(max_size)
    for size in range(min_size, max_size + 1):
        if not cts_by_size_and_frame.get(size, False):
            cts_by_size_and_frame[size] = [0, 0, 0]
    #for size in cts_by_size_and_frame: 
    #    print size, "\t".join([str(round(elt, 4)) for elt in cts_by_size_and_frame[size]])
    return cts_by_size_and_frame

def fill_cts_by_size_and_frame(cts_by_size_and_frame, min_size, max_size):
    for size in range(min_size, max_size + 1):
        if not cts_by_size_and_frame.get(size, False):
            cts_by_size_and_frame[size] = [0, 0, 0]
    return cts_by_size_and_frame

def normalize_cts_by_size_and_frame(cts_by_size_and_frame):
    norm_cts_by_size_and_frame = {}
    for size in cts_by_size_and_frame: 
        tot_cts = sum(cts_by_size_and_frame[size])        
        if tot_cts != 0:
            norm_cts_by_size_and_frame[size] = []
            for i in range(3):
                norm_cts_by_size_and_frame[size].append(cts_by_size_and_frame[size][i] / float(tot_cts))
    return norm_cts_by_size_and_frame

def load_fasta(in_fname):
  #Create dictionary to store {label:sequence} for fasta file                  i
  d = {}
  f = open(in_fname, "r")            #Open fasta file
  first = True                       #Bool for "Am I in the first seq?"

  for line in f:
    line = line.strip()
    if line[0] == ">":               #If line is a label:
      if not first:                      #If we finished an old seq
        d[name] = "".join(seq)               #Put old label:seq in dictionary
      name = line.split()[0][1:]         #Save new name
      seq = []                           #Start a new sequence
    else:
      seq.append(line)               #Else: add line to currentsequence
    first = False                    #Change start bool 
  d[name] = "".join(seq)             #Add final label:seq to dict
  return d

def get_cds_dict(gene_seq_fname, len_dict):
  #Create dictionary to store {label:sequence} for fasta file
  d = {}
  f = open(gene_seq_fname, "r")            #Open fasta file
  first = True                       #Bool for "Am I in the first seq?"

  for line in f:
    line = line.strip()
    if line == "": continue
    if line[0] == ">":               #If line is a label:
      if not first:                      #If we finished an old seq
        gene_seq = "".join(seq)                #get gene sequence
        cds_seq = gene_seq[cds_start_idx:cds_end_idx] #truncate to cds sequence
        d[name] = cds_seq                      #Put old label:seq in dictionary
      name = line.split()[0][1:]         #Save new name
      utr5_len = len_dict[name]["utr5"]
      cds_len = len_dict[name]["cds"]
      utr3_len = len_dict[name]["utr3"]
      cds_start_idx = utr5_len
      cds_end_idx = utr5_len + cds_len
      seq = []                           #Start a new sequence
    else:
      seq.append(line)               #Else: add line to currentsequence
    first = False                    #Change start bool 
  gene_seq = "".join(seq)
  cds_seq = gene_seq[cds_start_idx:cds_end_idx]
  d[name] = cds_seq                      #Add final label:seq to dict
  return d

def get_seq_dict(gene_seq_fname):
  #Create dictionary to store {label:sequence} for fasta file
  d = {}
  f = open(gene_seq_fname, "r")            #Open fasta file
  first = True                       #Bool for "Am I in the first seq?"

  for line in f:
    line = line.strip()
    if line == "": continue
    if line[0] == ">":               #If line is a label:
      if not first:                      #If we finished an old seq
        gene_seq = "".join(seq)                #get gene sequence
        d[name] = gene_seq                      #Put old label:seq in dictionary
      name = line.split()[0][1:]         #Save new name
      seq = []                           #Start a new sequence
    else:
      seq.append(line)               #Else: add line to currentsequence
    first = False                    #Change start bool 
  d[name] = "".join(seq)             #Add final label:seq to dict
  return d

def write_fasta(out_fname, d):
    f = open(out_fname, "w")
    for seq_name in d:
        f.write(">{0}\n".format(seq_name))
        seq = d[seq_name]
        seq_len = len(seq)
        num_segments = seq_len / 60
        segments = []
        for i in range(num_segments):
            segments.append(seq[i*60:(i+1)*60])
        if seq_len % 60 != 0:
            segments.append(seq[num_segments*60:])
        for segment in segments:
            f.write(segment + "\n")
    f.close()

def make_structure_fasta(out_fname, seq_dict, len_dict, window_len):
    f = open(out_fname, "w")
    for gene in seq_dict:
        utr5_len = len_dict[gene]["utr5"]
        gene_seq = seq_dict[gene]
        gene_len = len(gene_seq)
        num_windows = gene_len - window_len + 1
        for i in range(num_windows):
            pos = i - utr5_len
            window = gene_seq[i:i+window_len]
            f.write(">{0}\t{1}\n".format(gene, pos))
            f.write(window + "\n")
    f.close()

def get_len_dict(len_dict_fname):
    len_dict = {}
    f = open(len_dict_fname, "r")
    for line in f:
        line = line.strip().split()
        gene = line[0]
        try:
            utr5_len = int(line[1])
        except:
            print len_dict_fname
            print line
        cds_len = int(line[2])
        utr3_len = int(line[3])
        len_dict[gene] = {"utr5":utr5_len, "cds":cds_len, "utr3":utr3_len}
    f.close()
    return len_dict

def get_noutr_len_dict(cds_dict):
    len_dict = {}
    for gene in cds_dict:
        cds_len = len(cds_dict[gene])
        utr5_len = 0
        utr3_len = 0
        len_dict[gene] = {"utr5":utr5_len, "cds":cds_len, "utr3":utr3_len}
    return len_dict

def get_struc_dict(struc_fname):
    struc_dict = {}
    f = open(struc_fname, "r")
    line_idx = 0
    for line in f: 
        line = line.strip().split()
        if line_idx % 3 == 0:
            gene = line[0][1:]
            pos = int(line[1])
        if line_idx % 3 == 2:
            score = float(line[-1].strip("(").strip(")"))
            if not struc_dict.get(gene, False):
                struc_dict[gene] = {}
            struc_dict[gene][pos] = score
        line_idx += 1
    return struc_dict
    
def get_cts_by_codon(in_sam_fname, cds_dict, len_dict, shift_dict, min_fp_len, max_fp_len):
    #NOTE: Expects that mappings in in_sam are 1-indexed and done with 
    #      respect to the cds of the transcript to which a read is mapped.
    #      Expects a mapping count or weight (alone, no text) in the last field of each line
    #Coding note: need extended shift dict for variable [min|max]_len
    #Coding note: need to handle ambiguous mapping in future
    cts_by_codon = {gene:[0 for i in range(len(cds_dict[gene])/3)] for gene in cds_dict}
    in_sam = open(in_sam_fname, "r")
    for line in in_sam:
        if line[0] == "@":
            continue
        line = line.strip().split()
        read = line[9]
        if len(read) < min_fp_len or len(read) > max_fp_len:
            continue
        gene = line[2]
        utr5_len = len_dict[gene]["utr5"]
        map_nt = int(line[3]) - 1 - utr5_len
        frame = map_nt % 3
        shift = shift_dict[len(read)][frame]
        if not shift:
            continue
        gene = line[2]
        codon = (map_nt + shift) / 3
        if codon < 0 or codon >= len(cds_dict[gene]) / 3:
            continue
        map_wt = float(line[-1])
        cts_by_codon[gene][codon] += map_wt
    in_sam.close()
    return cts_by_codon

def write_cts_by_codon(out_fname, cts_by_codon):
    out_file = open(out_fname, "w")
    for gene in cts_by_codon:
        out_file.write("\t".join([gene] + [str(ct) for ct in cts_by_codon[gene]]) + "\n")
    out_file.close()

def load_cts_by_codon(in_fname):
    in_file = open(in_fname, "r")
    cts_by_codon = {}
    for line in in_file:
        line = line.strip().split()
        gene = line[0]
        cts = [float(elt) for elt in line[1:]]
        cts_by_codon[gene] = cts
    in_file.close()
    return cts_by_codon

def get_rel_cod_feats(gene_cds, A_site, rel_cod_idxs):
    cod_idxs = [A_site + idx for idx in rel_cod_idxs]
    cods = [gene_cds[idx*3:(idx+1)*3] for idx in cod_idxs]
    num_features = 64*len(cod_idxs)
    features = [0 for i in range(num_features)]
    for i, cod in enumerate(cods):
        features[i*64 + cod2id[cod]] = 1
    return np.array(features).reshape(-1, 1)

def get_rel_nt_feats(gene_cds, A_site, rel_nt_idxs):
    nt_idxs = [3*A_site + idx for idx in rel_nt_idxs]
    nts = [gene_cds[idx] for idx in nt_idxs]
    num_features = 4*len(nt_idxs)
    features = [0 for i in range(num_features)]
    for i, nt in enumerate(nts):
        features[i*4 + nt2id[nt]] = 1
    return np.array(features).reshape(-1, 1)

def get_rel_struc_feats(struc_dict, gene, A_site, rel_struc_idxs):
    struc_idxs = [3*A_site + idx for idx in rel_struc_idxs]
    features = [struc_dict[gene][idx] for idx in struc_idxs]
    return np.array(features).reshape(-1, 1)

def get_max_struc_feat(
        struc_dict, gene, A_site, max_struc_start_idx, max_struc_width):
    start_struc_idx = max_struc_start_idx
    end_struc_idx = max_struc_start_idx + max_struc_width #note, 1 too big!
    rel_struc_idxs = range(start_struc_idx, end_struc_idx) #fixes ^ !
    struc_idxs = [3*A_site + idx for idx in rel_struc_idxs]
    print max(struc_idxs), max(struc_dict[gene].keys())
    try:
        features = [struc_dict[gene][idx] for idx in struc_idxs]
    except KeyError: 
        print gene, idx, len(struc_dict[gene])
        #print struc_dict[gene]
    return min(features)

def get_outputs(cts_by_codon, cod_trunc_5p, cod_trunc_3p, raw_psct=0):
    outputs = {}
    for gene in cts_by_codon:
        cts_by_gene = copy.copy(cts_by_codon[gene])
        num_codons = len(cts_by_gene)
        num_codons_in_range = num_codons - cod_trunc_5p - cod_trunc_3p
        if num_codons_in_range <= 0:
            continue
        for i in range(cod_trunc_5p):
            cts_by_gene[i] = 0.0
        for i in range(num_codons - cod_trunc_3p, num_codons):
            cts_by_gene[i] = 0.0
        if raw_psct:
            for i in range(cod_trunc_5p, num_codons - cod_trunc_3p):
                cts_by_gene[i] += raw_psct
        tot_cts = sum(cts_by_gene)
        if tot_cts == 0:
            #print "{0} has no counts".format(gene)
            continue
        #Recompute tot_cs after zeroing out truncation regions
        avg_cts = float(tot_cts / num_codons_in_range)
        normed_cts_by_gene = [ct / avg_cts for ct in cts_by_gene]
        outputs[gene] = normed_cts_by_gene
    return outputs

def write_outputs(out_fname, outputs):
    out_file = open(out_fname, "w")
    for gene in outputs:
        out_file.write("\t".join([gene] + [str(elt) for elt in outputs[gene]]) + "\n")
    out_file.close()

def load_outputs(in_fname):
    in_file = open(in_fname, "r")
    outputs = {}
    for line in in_file:
        line = line.strip().split()
        gene = line[0]
        cts = [float(elt) for elt in line[1:]]
        outputs[gene] = cts
    in_file.close()
    return outputs

def log_transform_outputs(outputs, prior_ct):
    new_outputs = {}
    for gene in outputs:
        new_outputs[gene] = [math.log(elt + prior_ct) for elt in outputs[gene]]
    return new_outputs

def bin_transform_outputs(outputs, cutoff):
    new_outputs = {}
    for gene in outputs:
        new_outputs[gene] = [1 if elt >= cutoff else 0 for elt in outputs[gene]]
    return new_outputs

def shuffle_outputs(outputs, cdtr_5p, cdtr_3p):
    shuffled_outputs = {}
    for gene in outputs:
        num_cods = len(outputs[gene])
        num_in_rng = num_cods - cdtr_5p - cdtr_3p
        shuffled = random.sample(
            outputs[gene][cdtr_5p:cdtr_5p + num_in_rng], num_in_rng)
        shuffled_outputs[gene] = [0 for i in range(cdtr_5p)] + shuffled + [0 for i in range(cdtr_3p)]
    return shuffled_outputs

def shuffle_y(y):
    num_samples = y.size
    shuff_idxs = random.sample(range(num_samples), num_samples)
    y_shuff = y.take(shuff_idxs, axis=1)
    return y_shuff
 
def load_codon_set_bounds(in_fname):
    f = open(in_fname, "r")
    codon_set_bounds = {}
    for line in f:
        line = line.strip().split()
        if len(line) % 2 != 1:
            print "Error: codon bounds file has line with even # of entries"
            print "Line for gene {0}".format(line[0])
        gene = line[0]
        line = line[1:]
        num_intervals = len(line) / 2
        gene_bounds = [(int(line[i*2]), int(line[i*2 + 1])) for i in range(num_intervals)]
        codon_set_bounds[gene] = gene_bounds
    f.close()
    return codon_set_bounds

def write_codon_set_bounds(out_fname, codon_set_bounds):
    f = open(out_fname, "w")
    for gene in codon_set_bounds:
        intervals = codon_set_bounds[gene]
        bounds_list = [str(elt) for interval in intervals for elt in interval]
        fields = [gene] + bounds_list
        f.write("\t".join(fields) + "\n")
    f.close()

def make_codon_set_bounds(len_dict, genes, cod_trunc_5p, cod_trunc_3p):
    codon_set_bounds = {}
    for gene in genes:
        codon_len = len_dict[gene]
        start_idx = cod_trunc_5p
        end_idx = codon_len - cod_trunc_3p - 1
        if start_idx < codon_len and end_idx > 0 and end_idx >= start_idx:
            codon_set_bounds[gene] = [(start_idx, end_idx)]
    return codon_set_bounds
    
def expand_codon_set(codon_set_bounds):
    """
        Input: codon set bounds 
               {gene:[(int_start, int_end) for gene interval of codons in set]}
        Output: codon set {gene:[0-indexed codons in set for gene]}
    """
    codon_set = {}
    for gene in codon_set_bounds:
        gene_bounds = codon_set_bounds[gene]
        exp_intervals = [range(bound[0], bound[1] + 1) for bound in gene_bounds]
        gene_codons = [elt for interval in exp_intervals for elt in interval]
        codon_set[gene] = gene_codons
    return codon_set

def compress_codon_set(codon_set):
    """
        Input: codon set {gene:[0-indexed codons in set for gene]}
        Output: codon set bounds 
                {gene:[(int_start, int_end) for gene interval of codons in set]}
    """
    codon_set_bounds = {}
    for gene in codon_set:
        gene_codons = sorted(codon_set)
        gene_bounds = []
        int_start_idx = 0
        curr_codon = gene_codons[0]
        i = 1
        while i < len(gene_codons):
            prev_codon = curr_codon
            curr_codon = gene_codons[i]
            if curr_codon - prev_codon > 1:
                int_end_idx = i - 1
                gene_bounds.append((int_start_idx, int_end_idx))
                int_start_idx = i
            elif curr_codon - prev_codon < 1:
                print "Error: codons are closer than consecutive integers"
        gene_bounds.append((int_start_idx, gene_codons[-1]))
        codon_set_bounds[gene] = gene_bounds
    return codon_set_bounds

def check_codon_set(codon_set, rel_cod_idxs, rel_nt_idxs, cds_dict):
    """
        Checks that all codon and nt features fall within CDS for 
        all codons in codon_set
    """
    min_rel_cod = min(rel_cod_idxs)
    max_rel_cod = max(rel_cod_idxs)
    min_rel_nt = min(rel_nt_idxs)
    max_rel_nt = max(rel_nt_idxs)
    for gene in codon_set:
        cds_len = len(cds_dict[gene])
        sorted_codons = sorted(codon_set[gene])
        min_cod = sorted_codons[0]
        if min_cod + min_rel_cod < 0:
            print "Codon {0} of {1}\t{2} below 0".format(
                               min_rel_cod, gene, min_cod)
            return False
        max_cod = sorted_codons[-1]
        if max_cod + max_rel_cod >= cds_len / 3:
            print "Codon {0} of {1}\t{2} above {3}".format(
                    max_rel_cod, gene, max_cod, cds_len / 3)
            return False
        if min_cod * 3 + min_rel_nt < 0:
            print "nt {0} of {1}\t{2} below 0".format(
                         min_rel_nt, gene, min_cod)
            return False
        if max_cod * 3 + max_rel_nt >= cds_len:
            print "nt {0} of {1}\t{2} above {3}".format(
                    max_rel_nt, gene, max_cod, cds_len)  
            return False
    return True

def get_fp_density(cts_vector):
    return float(sum(cts_vector))/len(cts_vector)

def sort_genes_by_density(
        genes, cts_by_codon, cod_trunc_5p, cod_trunc_3p, descend=True):
    long_enough_genes = filter(
        lambda x: len(cts_by_codon[x]) > (cod_trunc_5p + cod_trunc_3p), genes)  
    genes_by_density = sorted(long_enough_genes, key=lambda x: get_fp_density(
        cts_by_codon[x][cod_trunc_5p:-cod_trunc_3p]), reverse=descend)
    return genes_by_density

def get_top_n_genes(cts_by_codon, cod_trunc_5p, cod_trunc_3p, n):
    genes_by_density = sort_genes_by_density(cts_by_codon, cod_trunc_5p, 
                                              cod_trunc_3p, descend=True)
    return genes_by_density[:n]

def split_gene_set(gene_set, n):
    num_genes = len(gene_set)
    sample_idxs = random.sample(range(num_genes), n)
    comp_idxs = np.setdiff1d(range(num_genes), sample_idxs)
    gene_set_1 = [gene_set[i] for i in sample_idxs]
    gene_set_2 = [gene_set[j] for j in comp_idxs]
    return gene_set_1, gene_set_2
                
def get_X(
        codon_set, cds_dict, rel_cod_idxs=False, rel_nt_idxs=False, 
        rel_struc_idxs=False, struc_dict=False, max_struc_start_idx=None, 
        max_struc_width=None, aa_feats=False):
    if not rel_cod_idxs and not rel_nt_idxs and not rel_struc_idxs:
        print "Error: must input some valid feature type"
        raise ValueError
    if bool(rel_struc_idxs) and not bool(struc_dict):
        print "Error: must input rel_struc_idxs " +\
            "and struc_dict if using structure"
        print struc_dict
        print rel_struc_idxs
        raise ValueError
    if (max_struc_start_idx is None) != (max_struc_width is None):
        print "Error: must input max_struc_start_idx " +\
            "and max_struc_width if using a max structure feature"
        print max_struc_start_idx
        print max_struc_width
        raise ValueError
    if not (max_struc_start_idx is None) and (struc_dict == False):
        print "Error: must input struc_dict if using a max structure feature"
        raise ValueError
    X_submatrices = []
    sorted_genes = sorted(codon_set.keys())
    if rel_cod_idxs:
        cod_X = np.hstack(
            [get_rel_cod_feats(cds_dict[gene], A_site, rel_cod_idxs)
            for gene in sorted_genes for A_site in codon_set[gene]])
        X_submatrices.append(cod_X)
    if aa_feats:
        aa_X = np.hstack(
            [get_rel_aa_feats(cds_dict[gene], A_site, rel_cod_idxs)
            for gene in sorted_genes for A_site in codon_set[gene]])
        X_submatrices.append(aa_X)
    if rel_nt_idxs:
        nt_X = np.hstack(
            [get_rel_nt_feats(cds_dict[gene], A_site, rel_nt_idxs)
            for gene in sorted_genes for A_site in codon_set[gene]])
        X_submatrices.append(nt_X)
    if rel_struc_idxs:
        struc_X = np.hstack(
            [get_rel_struc_feats(struc_dict, gene, A_site, rel_struc_idxs) 
            for gene in sorted_genes for A_site in codon_set[gene]])
        X_submatrices.append(struc_X)
    #print "max_struc_width: {0}".format(max_struc_width)
    if max_struc_width:
        print "woohoo! executing max struc!"
        max_struc_X = np.hstack(
            [get_max_struc_feat(
                struc_dict, gene, A_site, max_struc_start_idx, max_struc_width)
            for gene in sorted_genes for A_site in codon_set[gene]])
        X_submatrices.append(max_struc_X)
    #for matrix in X_submatrices: 
    #    print matrix.shape
    return np.vstack(X_submatrices)

def filter_outliers_by_pct(X_tr, y_tr, pct):
    num_pts = y_tr.size()
    cutoff_pctile = 1 - float(pct)/100
    cutoff_idx = int(num_pts * cutoff_pctile)
    sorted_y = sorted(y_tr.ravel())
    cutoff_val = sorted_y[cutoff_idx]
    i_vals, j_vals = np.where(y_tr < cutoff_val)
    y_tr_filtered = y_tr.take(j_vals, axis=1)
    X_tr_filtered = X_tr.take(j_vals, axis=1)
    return X_tr_filtered, y_tr_filtered

def filter_by_max_output(X_tr, y_tr, filter_max):
    pt_idxs, zeros = np.where(y_tr < filter_max)
    y_tr_filtered = y_tr.take(pt_idxs, axis=0)
    X_tr_filtered = X_tr.take(pt_idxs, axis=0)
    return X_tr_filtered, y_tr_filtered

def pickle_obj(X, pkl_fname):
    f = open(pkl_fname, "w")
    pickle.dump(X, f)
    f.close()

def load_obj(pkl_fname):
    f = open(pkl_fname, "r")
    X = pickle.load(f)
    return X

def get_y(codon_set, outputs):
    sorted_genes = sorted(codon_set.keys())
    return np.array([outputs[gene][i] for gene in sorted_genes
              for i in codon_set[gene]]).reshape(1, -1)

def get_default_vars():
    rel_cod_idxs = range(-5, 6)
    rel_nt_idxs = range(-20,-10) + range(8,18)
    cod_trunc_5p = 20
    cod_trunc_3p = 20
    return rel_cod_idxs, rel_nt_idxs, cod_trunc_5p, cod_trunc_3p

def write_genes(genes, f_name):
    f = open(f_name, "w")
    for gene in genes:
        f.write(gene + "\n")
    f.close()

def load_genes(f_name):
    genes = []
    f = open(f_name, "r")
    for line in f:
        genes.append(line.strip())
    f.close()
    return genes

def enforce_dir_exists(dir_name):
    if not os.path.exists(dir_name):
        print "ERROR: dir {0} does not exist".format(dir_name)
        raise NameError

def enforce_dir_doesnt_exist(dir_name):
    if os.path.exists(dir_name):
        print "ERROR: dir {0} already exists".format(dir_name)
        raise NameError

def make_expt_dirs(parent_dir, name):
    """
    In parent_dir, makes dir for expt 'name' with subdirs plots, process, 
    lasagne_nn, and linreg 

    Args: 
        parent_dir (str): parent directory to make subdir for expt 'name'
        name (str): name of experiment

    Returns: 
        void, just makes dirs
    """
    enforce_dir_exists(parent_dir)
    expt_dir = parent_dir + "/" + name
    enforce_dir_doesnt_exist(expt_dir)
    os.makedirs(expt_dir)
    os.makedirs(expt_dir + "/plots")
    #Deprecated, now we have lasagne_nn
    #os.makedirs(expt_dir + "/nn_data")
    os.makedirs(expt_dir + "/process")
    os.makedirs(expt_dir + "/lasagne_nn")
    os.makedirs(expt_dir + "/linreg")

def make_out_dir(out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

def make_nn_dir(out_dir, name):
    enforce_dir_exists(out_dir)
    dir_name = out_dir + "/" + name
    epoch_dirs = glob.glob(dir_name + "/" + "epoch*")
    if os.path.exists(dir_name):
        if len(epoch_dirs) > 0:
            print "BAD ERROR: dir with name {0} already exists".format(dir_name)
            raise NameError
        else: return
    else:
        os.makedirs(dir_name)

def make_init_data_dir(out_dir, name):
    dir_name = "{0}/{1}/init_data".format(out_dir, name)
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

def make_init_data_file(
        out_dir, name, gene_seq_fname, gene_len_fname, tr_codons_fname, 
        te_codons_fname, outputs_fname, cost_fn, act_fn, num_hidden,
        num_outputs, lam, rel_cod_idxs=False, rel_nt_idxs=False):
    f = open("{0}/{1}/init_data/init_data.txt".format(out_dir, name), "w")
    #f.write("**************\nInit data\n**************\n")
    f.write("Transcripts_fasta_file: {0}".format(gene_seq_fname) + "\n")
    f.write("Transcripts_length_file: {0}".format(gene_len_fname) + "\n")
    f.write("Training_codon_bounds: {0}".format(tr_codons_fname) + "\n")
    f.write("Test_codon_bounds: {0}".format(te_codons_fname) + "\n")
    f.write("Outputs: {0}".format(outputs_fname) + "\n")
    f.write("Cost_fn: {0}".format(cost_fn) + "\n")
    f.write("Activation_fn: {0}".format(act_fn) + "\n")
    f.write("Hidden_units_by_layer: {0}".format("\t".join([str(elt) for elt in num_hidden]) + "\n"))
    f.write("Number_outputs: {0}".format(num_outputs) + "\n")
    f.write("L2_reg_parameter: {0}".format(lam) + "\n")
    #f.write("\n")
    #f.write("**************\nInput Feature Neighborhoods\n**************\n")
    if rel_cod_idxs:
        f.write("Relative_codon_indices: {0}".format("\t".join([str(elt) for elt in rel_cod_idxs])) + "\n")
    if rel_nt_idxs:
        f.write("Relative_nt_indices: {0}".format("\t".join([str(elt) for elt in rel_nt_idxs])) + "\n")
    f.close()

def load_init_data_file(init_data_fname):
    d = {}
    f = open(init_data_fname, "r")
    for line in f:
        line = line.strip().split()
        key = line[0]
        value = "\t".join(line[1:])
        d[key] = value
    f.close()
    if d.get("Transcripts_fasta_file:", False): d["gene_seq_fname"] = d["Transcripts_fasta_file:"]
    if d.get("Transcripts_length_file:", False): d["gene_len_fname"] = d["Transcripts_length_file:"]
    if d.get("Training_codon_bounds:", False): d["tr_codons_fname"] = d["Training_codon_bounds:"]
    if d.get("Test_codon_bounds:", False): d["te_codons_fname"] = d["Test_codon_bounds:"]
    if d.get("Outputs:", False): d["outputs_fname"] = d["Outputs:"]
    if d.get("Cost_fn:", False): d["cost_fn"] = d["Cost_fn:"]
    if d.get("Activation_fn:", False): d["act_fn"] = d["Activation_fn:"]
    if d.get("Hidden_units_by_layer:", False): d["num_hidden"] = [int(elt) for elt in d["Hidden_units_by_layer:"].split()]
    if d.get("Number_outputs:", False): d["num_outputs"] = int(d["Number_outputs:"])
    if d.get("L2_reg_parameter:", False): d["lam"] = float(d["L2_reg_parameter:"])
    if d.get("Relative_codon_indices:", False): d["rel_cod_idxs"] = [int(elt) for elt in d["Relative_codon_indices:"].split()]
    if d.get("Relative_nt_indices:", False): d["rel_nt_idxs"] = [int(elt) for elt in d["Relative_nt_indices:"].split()]
    return d

def load_linreg_data_file(linreg_data_fname):
    d = {}
    f = open(linreg_data_fname, "r")
    for line in f:
        line = line.strip().split()
        key = line[0]
        value = "\t".join(line[1:])
        d[key] = value
    f.close()
    if d.get("Transcripts_fasta_file:", False): 
        d["gene_seq_fname"] = d["Transcripts_fasta_file:"]
    if d.get("Transcripts_length_file:", False): 
        d["gene_len_fname"] = d["Transcripts_length_file:"]
    if d.get("Training_codon_bounds:", False): 
        d["tr_codons_fname"] = d["Training_codon_bounds:"]
    if d.get("Test_codon_bounds:", False): 
        d["te_codons_fname"] = d["Test_codon_bounds:"]
    if d.get("Outputs:", False): 
        d["outputs_fname"] = d["Outputs:"]
    if d.get("Relative_codon_indices:", False): 
        d["rel_cod_idxs"] = \
            [int(elt) for elt in d["Relative_codon_indices:"].split()]
    if d.get("Relative_nt_indices:", False): 
        d["rel_nt_idxs"] = \
            [int(elt) for elt in d["Relative_nt_indices:"].split()]
    return d

def has_enough_cts(gene_cts_by_codon, min_tot_cts, min_cod_w_cts):
    num_codons = len(gene_cts_by_codon)
    tot_cts = sum(gene_cts_by_codon)
    cod_w_cts = sum([elt > 0 for elt in gene_cts_by_codon])
    return (tot_cts >= min_tot_cts) and (cod_w_cts >= min_cod_w_cts)    

def make_codon_set_files(
        cds_dict, cts_by_codon_fname, outputs_fname, tr_out_fname, 
        te_out_fname, num_tr_genes, num_te_genes, cod_trunc_5p, cod_trunc_3p, 
        min_tot_cts, min_cod_w_cts, skip_genes=0, overwrite=False, verbose=True,
        paralog_groups_fname=False):
    cts_by_codon = load_cts_by_codon(cts_by_codon_fname)
    tot_set_sizes = num_tr_genes + num_te_genes
    gene_set = cts_by_codon.keys()
    #Filter genes that are shorter than truncation regions
    gene_set = filter(lambda gene: len(cts_by_codon[gene]) > (cod_trunc_5p +
                               cod_trunc_3p), gene_set)
    #Filter genes that don't have enough counts to meet cutoffs
    gene_set = filter(lambda gene: has_enough_cts( \
        cts_by_codon[gene][cod_trunc_5p:-cod_trunc_3p], min_tot_cts, \
        min_cod_w_cts), gene_set)
    genes_by_density = sort_genes_by_density(
        gene_set, cts_by_codon, cod_trunc_5p, cod_trunc_3p, descend=True)
    if paralog_groups_fname:
        paralog_groups = load_paralog_groups(paralog_groups_fname)
        genes_by_density = filter_paralogs(genes_by_density, paralog_groups)
    if verbose:
        print "{0} genes meeting data cutoff".format(len(genes_by_density))
    if tot_set_sizes + skip_genes > len(genes_by_density):
        print "Error: Skip {0} + Get {1} > {2} Post-Filter Genes".format(
            skip_genes, tot_set_sizes, len(genes_by_density))
        raise ValueError
    gene_set = genes_by_density[skip_genes:skip_genes + tot_set_sizes]
    len_dict = {gene:len(cts_by_codon[gene]) for gene in gene_set}
    tr_set, te_set = split_gene_set(gene_set, num_tr_genes)
    tr_set_bounds = make_codon_set_bounds(
        len_dict, tr_set, cod_trunc_5p, cod_trunc_3p)
    te_set_bounds = make_codon_set_bounds(
        len_dict, te_set, cod_trunc_5p, cod_trunc_3p)
    if not overwrite:
        if not os.path.isfile(tr_out_fname):
            write_codon_set_bounds(tr_out_fname, tr_set_bounds)
        else:
            print "training set file {0} already exists".format(tr_out_fname) \
            + "\nUse parameter overwrite=True to overwrite"
        if not os.path.isfile(te_out_fname):
            write_codon_set_bounds(te_out_fname, te_set_bounds)
        else:
            print "test set file {0} already exists".format(te_out_fname) \
            + "\nUse parameter overwrite=True to overwrite"
    else:
        write_codon_set_bounds(tr_out_fname, tr_set_bounds)
        write_codon_set_bounds(te_out_fname, te_set_bounds)
    tr_data_table_fname = tr_out_fname[:-4] + ".data_table.txt"
    te_data_table_fname = te_out_fname[:-4] + ".data_table.txt"
    make_set_data_table(
        tr_out_fname, cds_dict, cts_by_codon_fname, outputs_fname, 
        tr_data_table_fname)
    make_set_data_table(
        te_out_fname, cds_dict, cts_by_codon_fname, outputs_fname, 
        te_data_table_fname)

def load_paralog_groups(paralog_groups_fname):
    paralog_groups = []
    f = open(paralog_groups_fname)
    for line in f: 
        p_group = line.strip().split()
        paralog_groups.append(p_group)
    f.close()
    return paralog_groups
   
def add_unique_paralogs(
        paralog_groups_fname, gene_seq_fname, out_paralog_groups_fname):
    """
        Input: -file with paralog groups (possibly lacking singletons)
               -file with gene names
               -output filename
        Output: -None, remakes paralog groups file at output filename, but 
                    adds singleton paralogs (i.e. genes w. no paralogs) on 
                    their own lines. 
    """
    # Step 1: get full list of genes
    f = open(gene_seq_fname, "r")
    all_genes = []
    for line in f:
        line = line.strip().split()
        if line[0][0] == ">":
            gene_name = line[0][1:]
            all_genes.append(gene_name)
    f.close()
    # Step 2: get list of paralog groups
    paralog_groups = load_paralog_groups(paralog_groups_fname)
    paralog_genes = [gene for p_group in paralog_groups for gene in p_group]
    is_paralog_gene = {gene:True for gene in paralog_genes}
    # Step 3: Add paralog singletons to paralog groups
    for gene in all_genes: 
        if not is_paralog_gene.get(gene, False):
            #print gene
            paralog_groups.append([gene])
    # Step 4: Write paralog groups file at out_fname
    g = open(out_paralog_groups_fname, "w")
    for group in paralog_groups: 
        g.write("\t".join(group) + "\n")
    g.close()
    # Step 5: Some random analyses: 
    print "Total genes: " + str(len(all_genes))
    print "Singleton genes: " +\
        str(len(filter(lambda x: len(x) == 1, paralog_groups)))
    print "Number paralog groups: " +\
        str(len(filter(lambda x: len(x) > 1, paralog_groups)))
    print "Number paralog genes: " +\
        str(sum([len(group) for group in filter(lambda x: len(x) > 1, paralog_groups)]))
    for i in range(115):
        print "Genes in group > {0}: ".format(i) +\
            str(sum([len(group) for group in filter(lambda x: len(x) > i, paralog_groups)]))

def filter_paralogs(genes_by_density, paralog_groups):
    """
        Input: -list of genes sorted by data density IN DESCENDING ORDER
               -2D list, outer list is groups of paralogs, inner lists are
                   sets of genes that are mutually paralogs
        Output: -Filtered list of genes sorted by data density, taking only 
                   the top paralog in a group by data density
    """
    print "filtering paralogs"
    accepted_genes = 0
    rejected_genes = 0
    para_group_used = [False for p_group in paralog_groups]
    gene2groupidx = {}
    #Remove after testing
    groupidx2reprgene = {}
    #End remove after testing
    for idx, p_group in enumerate(paralog_groups): 
        for gene in p_group: 
            gene2groupidx[gene] = idx
    filtered_genes_by_density = []
    for gene in genes_by_density:
        group_idx = gene2groupidx[gene]
        if not para_group_used[group_idx]:
            filtered_genes_by_density.append(gene)
            print "accept gene: " + gene
            para_group_used[group_idx] = True
            groupidx2reprgene[group_idx] = gene
            accepted_genes += 1
        else: 
            print "reject gene: " + gene
            print "    already have paralog: " + groupidx2reprgene[group_idx]
            rejected_genes += 1
    print "finished filtering paralogs"
    return filtered_genes_by_density

def get_data_matrices(
        cds_dict, tr_codons_fname, te_codons_fname, outputs_fname, 
        rel_cod_idxs=False, rel_nt_idxs=False, rel_struc_idxs=False, 
        struc_dict=False):
    tr_codon_bounds = load_codon_set_bounds(tr_codons_fname)
    te_codon_bounds = load_codon_set_bounds(te_codons_fname)
    tr_codon_set = expand_codon_set(tr_codon_bounds)
    te_codon_set = expand_codon_set(te_codon_bounds)
    outputs = load_outputs(outputs_fname)
    X_tr = get_X(tr_codon_set, cds_dict, rel_cod_idxs=rel_cod_idxs,
                 rel_nt_idxs=rel_nt_idxs, rel_struc_idxs=rel_struc_idxs, 
                 struc_dict=struc_dict)
    X_te = get_X(te_codon_set, cds_dict, rel_cod_idxs=rel_cod_idxs,
                 rel_nt_idxs=rel_nt_idxs, rel_struc_idxs=rel_struc_idxs, 
                 struc_dict=struc_dict)
    y_tr = get_y(tr_codon_set, outputs)
    y_te = get_y(te_codon_set, outputs)
    return X_tr, X_te, y_tr, y_te

def get_y_data(tr_codons_fname, te_codons_fname, outputs_fname):
    #Load codon sets, cts_by_codon, outputs, y_tr/te
    tr_codon_bounds = load_codon_set_bounds(tr_codons_fname)
    te_codon_bounds = load_codon_set_bounds(te_codons_fname)
    tr_codon_set = expand_codon_set(tr_codon_bounds)
    te_codon_set = expand_codon_set(te_codon_bounds)
    outputs = load_outputs(outputs_fname)
    y_tr = get_y(tr_codon_set, outputs)
    y_te = get_y(te_codon_set, outputs)
    return y_tr, y_te, tr_codon_set, te_codon_set, outputs

def process_data(
        out_dir, name, gene_seq_fname, gene_len_fname, tr_codons_fname, 
        te_codons_fname, outputs_fname, cost_fn, act_fn, num_hidden, 
        num_outputs, lam, rel_cod_idxs=False, rel_nt_idxs=False, 
        filter_max=False, filter_pct=False, rel_struc_idxs=False, 
        struc_fname=False):
    if filter_max and filter_pct:
        print "Error! Must pick one filter [max, pct] for codon data"
        return
    make_out_dir(out_dir)
    make_nn_dir(out_dir, name)
    make_init_data_dir(out_dir, name)
    make_init_data_file(out_dir, name, gene_seq_fname, gene_len_fname, 
        tr_codons_fname, te_codons_fname, outputs_fname, 
        cost_fn, act_fn, num_hidden, num_outputs, lam, 
        rel_cod_idxs=rel_cod_idxs, rel_nt_idxs=rel_nt_idxs)
    len_dict = get_len_dict(gene_len_fname)
    cds_dict = get_cds_dict(gene_seq_fname, len_dict)
    if struc_fname:
        struc_dict = get_struc_dict(struc_fname)
    else:
        struc_dict = False
    X_tr, X_te, y_tr, y_te = get_data_matrices(
        cds_dict, tr_codons_fname, te_codons_fname, outputs_fname, 
        rel_cod_idxs=rel_cod_idxs, rel_nt_idxs=rel_nt_idxs,
        rel_struc_idxs=rel_struc_idxs, struc_dict=struc_dict)
    if filter_max: 
        X_tr, y_tr = filter_by_max_output(X_tr, y_tr, filter_max)
    if filter_pct:
        X_tr, y_tr = filter_outliers_by_pct(X_tr, y_tr, filter_pct)
    pickle_obj(y_tr, "{0}/{1}/init_data/y_tr.pkl".format(out_dir, name))
    pickle_obj(y_te, "{0}/{1}/init_data/y_te.pkl".format(out_dir, name))
    return X_tr, X_te, y_tr, y_te
        
def load_lasagne_data(
        gene_len_fname, gene_seq_fname, tr_codons_fname, te_codons_fname, 
        outputs_fname, rel_cod_idxs=False, rel_nt_idxs=False, 
        rel_struc_idxs=False, struc_fname=False, max_struc_start_idx=None,
        max_struc_width=None, aa_feats=False,
        filter_max = False, filter_test=False, filter_pct=False):
    if filter_max and filter_pct:
        print "Error! Must pick one filter [max, pct] for codon data"
        return
    len_dict = get_len_dict(gene_len_fname)
    cds_dict = get_cds_dict(gene_seq_fname, len_dict)
    if struc_fname:
        struc_dict = get_struc_dict(struc_fname)
    else:
        struc_dict = False
    X_tr, y_tr, X_te, y_te = get_data_matrices_lasagne(
        cds_dict, tr_codons_fname, te_codons_fname, outputs_fname,
        rel_cod_idxs=rel_cod_idxs, rel_nt_idxs=rel_nt_idxs,
        rel_struc_idxs=rel_struc_idxs, struc_dict=struc_dict,
        max_struc_start_idx=max_struc_start_idx, 
        max_struc_width=max_struc_width, aa_feats=aa_feats)
    if filter_max: 
        X_tr, y_tr = filter_by_max_output(X_tr, y_tr, filter_max)
        if filter_test: 
            X_te, y_te = filter_by_max_output(X_te, y_te, filter_max)
    if filter_pct:
        X_tr, y_tr = filter_outliers_by_pct(X_tr, y_tr, filter_pct)
    return X_tr, y_tr, X_te, y_te

def get_data_matrices_lasagne(
        cds_dict, tr_codons_fname, te_codons_fname, outputs_fname, 
        rel_cod_idxs, rel_nt_idxs, rel_struc_idxs, struc_dict, 
        max_struc_start_idx, max_struc_width, aa_feats):

    # Load tr_codon_bounds
    if type(tr_codons_fname) == str:
        tr_codon_bounds = load_codon_set_bounds(tr_codons_fname)
    elif type(tr_codons_fname) == list:
        tr_bounds_subsets = \
            [load_codon_set_bounds(tr_subset_fname) 
                for tr_subset_fname in tr_codons_fname]
        tr_codon_bounds = {k:v for d in tr_bounds_subsets for k,v in d.items()}
    
    # Load te_codon_bounds
    if type(te_codons_fname) == str:
        te_codon_bounds = load_codon_set_bounds(te_codons_fname)
    elif type(te_codons_fname) == list:
        te_bounds_subsets = \
            [load_codon_set_bounds(te_subset_fname) 
                for te_subset_fname in te_codons_fname]
        te_codon_bounds = {k:v for d in te_bounds_subsets for k,v in d.items()}

    tr_codon_set = expand_codon_set(tr_codon_bounds)
    te_codon_set = expand_codon_set(te_codon_bounds)

    #rel_struc_idxs=False
    #struc_dict=False
    X_tr = get_X(tr_codon_set, cds_dict, rel_cod_idxs=rel_cod_idxs,
             rel_nt_idxs=rel_nt_idxs, rel_struc_idxs=rel_struc_idxs,
             struc_dict=struc_dict, max_struc_start_idx=max_struc_start_idx,
             max_struc_width=max_struc_width, aa_feats=aa_feats)
    X_te = get_X(te_codon_set, cds_dict, rel_cod_idxs=rel_cod_idxs,
             rel_nt_idxs=rel_nt_idxs, rel_struc_idxs=rel_struc_idxs,
             struc_dict=struc_dict, max_struc_start_idx=max_struc_start_idx,
             max_struc_width=max_struc_width, aa_feats=aa_feats)

    outputs = load_outputs(outputs_fname)
    y_tr = get_y(tr_codon_set, outputs)
    y_te = get_y(te_codon_set, outputs)

    X_train = X_tr.transpose(1, 0)
    X_test = X_te.transpose(1, 0)
    y_train = y_tr.transpose(1, 0)
    y_test = y_te.transpose(1, 0)
    return X_train, y_train, X_test, y_test

def make_linreg_parent_dir(expt_dir):
    parent_dir = "{0}/linreg".format(expt_dir)
    if not os.path.exists(parent_dir):
        os.makedirs(parent_dir)

def make_linreg_dir(expt_dir, name):
    dir_name = "{0}/linreg/{1}".format(expt_dir, name)
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

####################
### LSTM Functions
####################
def get_ragged_y_matrix(codon_set, outputs):
    sorted_genes = sorted(codon_set.keys())
    num_genes = len(sorted_genes)
    max_gene_len = max([len(codon_list) for codon_list in codon_set.values()])
    y_mat = np.zeros((num_genes, max_gene_len))
    for i, gene in enumerate(sorted_genes):
        for j, cod_idx in enumerate(codon_set[gene]):
            y_mat[i, j] = outputs[gene][cod_idx]
    return y_mat

def get_ragged_mask(codon_set):
    sorted_genes = sorted(codon_set.keys())
    num_genes = len(sorted_genes)
    max_gene_len = max([len(codon_list) for codon_list in codon_set.values()])
    y_mask = np.zeros((num_genes, max_gene_len), dtype="int64")
    for i, gene in enumerate(sorted_genes):
        for j, cod_idx in enumerate(codon_set[gene]):
            y_mask[i, j] = 1
    return y_mask

def get_cod_feat_vec(codon):
    vec = [0 for i in range(64)]
    vec[cod2id[codon]] = 1
    return vec

def get_nt_feat_vec(nt):
    vec = [0 for i in range(4)]
    vec[nt2id[nt]] = 1
    return vec
 
def get_ragged_X_matrix(codon_set, cds_dict, nt=False):
    sorted_genes = sorted(codon_set.keys())
    num_genes = len(sorted_genes)
    max_gene_len = max([len(codon_list) for codon_list in codon_set.values()])
    num_feats = 64
    if nt: 
        num_feats += 12
    X_mat = []
    for gene in sorted_genes:
        gene_mat = []
        for cod_idx in codon_set[gene]:
            codon = cds_dict[gene][cod_idx * 3 : (cod_idx + 1) * 3]
            feat_vec = get_cod_feat_vec(codon)
            if nt:
                #print codon
                for char in codon:
                    #print char
                    feat_vec += get_nt_feat_vec(char)
            gene_mat.append(feat_vec)
        for i in range(max_gene_len - len(codon_set[gene])):
            gene_mat.append([0 for j in range(num_feats)])
        X_mat.append(gene_mat)
    X_mat = np.array(X_mat)
    return X_mat

def get_rel_cod_feats(gene_cds, A_site, rel_cod_idxs):
    cod_idxs = [A_site + idx for idx in rel_cod_idxs]
    cods = [gene_cds[idx*3:(idx+1)*3] for idx in cod_idxs]
    num_features = 64*len(cod_idxs)
    features = [0 for i in range(num_features)]
    for i, cod in enumerate(cods):
        features[i*64 + cod2id[cod]] = 1
    return np.array(features).reshape(-1, 1)

def get_rel_nt_feats(gene_cds, A_site, rel_nt_idxs):
    nt_idxs = [3*A_site + idx for idx in rel_nt_idxs]
    nts = [gene_cds[idx] for idx in nt_idxs]
    num_features = 4*len(nt_idxs)
    features = [0 for i in range(num_features)]
    for i, nt in enumerate(nts):
        features[i*4 + nt2id[nt]] = 1
    return np.array(features).reshape(-1, 1)

def get_rel_aa_feats(gene_cds, A_site, rel_cod_idxs):
    cod_idxs = [A_site + idx for idx in rel_cod_idxs]
    cods = [gene_cds[idx*3:(idx+1)*3] for idx in cod_idxs]
    aas = [cod2aa[cod] for cod in cods]
    num_features = 21*len(cod_idxs)
    features = [0 for i in range(num_features)]
    for i, aa in enumerate(aas):
        features[i*21 + aa2id[aa]] = 1
    return np.array(features).reshape(-1, 1)

def get_gene2symbol_dict(gene2sym_fname, sysmod_fn=False):
    #Make dictionary to convert systematic gene name to common symbol
    f = open(gene2sym_fname, "r")
    gene2symbol = {}
    for line in f:
        systematic, gene_symbol = line.strip().split("\t")
        if sysmod_fn:
            systematic = sysmod_fn(systematic)
        gene2symbol[systematic] = gene_symbol
    f.close()
    return gene2symbol
    
def get_mat_idxs_by_gene(sorted_gene_list, codon_set):
    #Get slicing indices for genes in y vector (or x matrix)
    start_idx = 0
    gene_data = {}
    for gene in sorted_gene_list:
        num_codons = len(codon_set[gene])
        gene_data[gene] = { "length":num_codons, "start_idx":start_idx }
        start_idx += num_codons
    return gene_data

def save_gene_list(gene_list, out_fname, gene2symbol):
    f = open(out_fname, "w")
    for gene in gene_list:
        f.write("{0}\t{1}\n".format(gene, gene2symbol.get(gene, "")))
    f.close()

def get_gene_subset_test_idxs(
        te_codons_fname, partial_gene_set, cts_by_codon_fname, outputs_fname):

    subset_test_idxs = []

    #Load codon sets, cts_by_codon, outputs, y_tr/te
    cts_by_codon = load_cts_by_codon(cts_by_codon_fname)
    te_codon_bounds = load_codon_set_bounds(te_codons_fname)
    te_codon_set = expand_codon_set(te_codon_bounds)
    outputs = load_outputs(outputs_fname)
    y_te = get_y(te_codon_set, outputs)

    te_sorted_genes = sorted(te_codon_set.keys())
    te_gene_data = get_mat_idxs_by_gene(te_sorted_genes, te_codon_set)
    
    for gene in sorted(partial_gene_set):
        start_idx = te_gene_data[gene]["start_idx"]
        num_codons = te_gene_data[gene]["length"]
        subset_test_idxs.extend(range(start_idx, start_idx + num_codons))

    return subset_test_idxs
       
def get_set_data_table(
        set_bounds_fname, cds_dict, cts_by_codon_fname, outputs_fname):
    cts_by_codon = load_cts_by_codon(cts_by_codon_fname)
    outputs = load_outputs(outputs_fname)
    f = open(set_bounds_fname, "r")
    set_bounds_data = []
    #Parse gene names and bounds from set_bounds file
    for line in f:
        gene, start_idx, end_idx = line.strip().split() 
        start_idx = int(start_idx)
        end_idx = int(end_idx)
        set_bounds_data.append((gene, start_idx, end_idx))
    set_bounds_data = sorted(set_bounds_data, key=lambda x: x[0])
    f.close()
    #Populate gene_names, codon_idxs vectors
    gene_names = []
    codon_idxs = []
    for gene, start_idx, end_idx in set_bounds_data:
        for idx in range(start_idx, end_idx + 1):
            gene_names.append(gene)
            codon_idxs.append(idx)
    #Make raw_cts, scaled_cts, and cod_seq vectors
    raw_cts = []
    scaled_cts = []
    cod_seqs = []
    for i in range(len(gene_names)):
        gene = gene_names[i]
        cod_idx = codon_idxs[i]
        raw_cts.append(cts_by_codon[gene][cod_idx])
        scaled_cts.append(outputs[gene][cod_idx])
        cod_seqs.append(cds_dict[gene][cod_idx*3:(cod_idx+1)*3])
    #Return set row labels
    return zip(gene_names, codon_idxs, cod_seqs, raw_cts, scaled_cts)

def save_set_data_table(set_data_table, out_fname):
    f = open(out_fname, "w")
    f.write("gene\tcod_idx\tcod_seq\traw_cts\tscaled_cts\n")
    for row in set_data_table:
        f.write("\t".join([str(elt) for elt in row]) + "\n")
    f.close()

def load_set_data_table(in_fname):
    data_table = []
    f = open(in_fname, "r")
    f.readline()
    for line in f:
        #Data table has columns 1) gene 2) cod_idx 3) cod_seq 4) raw_cts 
            #5) scaled_cts
        gene, cod_idx, cod_seq, raw_cts, scaled_cts = line.strip().split() 
        cod_idx = int(cod_idx)
        raw_cts = float(raw_cts)
        scaled_cts = float(scaled_cts)
        data_table.append((gene, cod_idx, cod_seq, raw_cts, scaled_cts))
    f.close()
    return data_table

def make_set_data_table(
        set_bounds_fname, cds_dict, cts_by_codon_fname, outputs_fname, 
        out_fname):
    set_data_table = get_set_data_table(
        set_bounds_fname, cds_dict, cts_by_codon_fname, outputs_fname)
    save_set_data_table(set_data_table, out_fname)

def get_genes_by_density(
        set_bounds_fname, cts_by_codon_fname, cod_trunc_5p, cod_trunc_3p):
    cts_by_codon = load_cts_by_codon(cts_by_codon_fname)
    genes = []
    f = open(set_bounds_fname)
    for line in f:
        genes.append(line.strip().split()[0]) 
    genes_by_density = sort_genes_by_density(
        genes, cts_by_codon, cod_trunc_5p, cod_trunc_3p)    
    return genes_by_density

def get_genes_by_density2(te_set_data_table):
    #Compute te_genes_by_density
    te_genes_by_pos = np.array([row[0] for row in te_set_data_table])
    cts_by_pos = np.array([row[3] for row in te_set_data_table])
    te_genes = np.unique(te_genes_by_pos)

    def get_cts_density(gene):
        gene_row_bools = (te_genes_by_pos == gene)
        tot_cts = sum(cts_by_pos[gene_row_bools])
        num_cods = sum(gene_row_bools)
        return float(tot_cts)/num_cods

    te_genes_by_density = sorted(
        list(te_genes), key=lambda gene: get_cts_density(gene), reverse=True)
    return te_genes_by_density

