#print "importing libraries"
import os
import inspect
import iXnos.process as proc
import iXnos.plot as plot
import iXnos.linreg as linreg
import iXnos.optimizecodons as opt
import iXnos.lasagnenn as lasagnenn
import iXnos.featureanalysis as feat
import numpy as np
from shutil import copyfile
import math
#print "...finished importing libraries"

def fp_size_analysis(
        cts_by_size_and_frame, expt_dir, min_size=13, max_size=36):
        #sam_fname, gene_len_fname, plot_fname, min_size=13, max_size=36, 
        #verbose=False):
    """
    Makes cts_by_size plot in expt_dir/plots

    Args:
        cts_by_size_and_frame (dict):
            {fp_size (int): 
                [cts_frame_0 (int), cts_frame_1 (int), cts_frame_2 (int)]
                for fp_size in range(min_fp_size, max_fp_size in sam file)}
        expt_dir (str) - name of experiment directory
        min_size (int) - minimum footprint size to plot
        max_size (int) - maximum footprint size to plot

    Returns: 
        void, just makes cts by size plot
    """
    # cts_by_size plot file name
    size_plot_fname = expt_dir + "/plots/cts_by_size.pdf"
    # make cts_by_size plot
    plot.make_cts_by_size_plot(
        cts_by_size_and_frame, "Counts by FP Size",
        size_plot_fname, sizes=range(min_size, max_size + 1))

def fp_frame_by_size_analysis(
        cts_by_size_and_frame, expt_dir, min_size=13, max_size=36):
    """
    Makes cts_by_size_and_frame plots in expt_dir/plots

    Args:
        cts_by_size_and_frame (dict):
            {fp_size (int): 
                [cts_frame_0 (int), cts_frame_1 (int), cts_frame_2 (int)]
                for fp_size in range(min_fp_size, max_fp_size in sam file)}
        expt_dir (str) - name of experiment directory
        min_size (int) - minimum footprint size to plot
        max_size (int) - maximum footprint size to plot

    Returns: 
        void, just makes cts by size and frame plots
    """
    # for each footprint size
    for size in range(min_size, max_size+1):
        # cts by size and frame plot file name
        plot_fname = expt_dir + "/plots/size_{0}_by_frame.pdf".format(size)
        # make cts by size and frame plot
        plot.make_frame_by_size_plot(
            size, cts_by_size_and_frame[size],
            "Size {0} by Frame".format(size), plot_fname, cts=True)

def do_frame_and_size_analysis(
        expt_dir, sam_fname, gene_len_fname, min_size=13, max_size=36, 
        verbose=False):
    """
    Makes cts_by_size and cts_by_size_and_frame plots in expt_dir/plots

    Args:
        expt_dir (str) - name of experiment directory
        sam_fname (str) - name of input sam file
        gene_len_fname (str) - name of gene lengths file
        plot_fname (str) - output directory for plotting
        min_size (int) - minimum footprint size to plot
        max_size (int) - maximum footprint size to plot
        verbose (bool) - flag to print max/min fp size in sam file

    Returns: 
        void, just makes cts by size and cts by size and frame plots
    """
    #NOTE: Set verbose=True if you want to see min and max fp size in sam_fname
    # Get plot directory for experiment
    plot_dir = expt_dir + "/plots"
    # load cts_by_size_and_frame dictionary
    len_dict = proc.get_len_dict(gene_len_fname)
    cts_by_size_and_frame = proc.get_cts_by_size_and_frame(
        sam_fname, len_dict, verbose=verbose)
    cts_by_size_and_frame = proc.fill_cts_by_size_and_frame(
        cts_by_size_and_frame, min_size, max_size)
    #Make cts by size plot
    fp_size_analysis(
        cts_by_size_and_frame, expt_dir, min_size=min_size, max_size=max_size)
    #Make cts by size and frame plots
    fp_frame_by_size_analysis(
        cts_by_size_and_frame, expt_dir, min_size=min_size, max_size=max_size)

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
    proc.make_expt_dirs(parent_dir, name)

def get_fname_nodir(fname):
    """
    Removes directory prefix (i.e. anything before and including '/') from a 
    file name

    Args: 
        fname (str): string, name of a file

    Returns: 
        fname_nodir: fname with any leading directories stripped
    """
    last_slash_idx = fname.rfind("/")
    fname_nodir = fname[last_slash_idx + 1:]
    return fname_nodir

def edit_sam_file(
        expt_dir, sam_fname, filter_unmapped=False, 
        sam_add_simple_map_wts=False, RSEM_add_map_wts=False, 
        sample_frac=False):
    """
    Does edit operations on sam file of mapped ribosome footprints

    Args:
        expt_dir (str) - name of experiment directory
        sam_fname (str) - name of input sam file
            NOTE: sam_fname must end with '.sam'
        filter_unmapped (bool): 
            filter out unmapped reads from sam file
        sam_add_simple_map_wts (bool): 
            add field at end of sam file with uniform mapping weights for 
            multimapped reads
        RSEM_add_map_wts (bool):
            add field at end of sam file with rsem mapping weight
            NOTE: REQUIRES sam_fname TO BE RSEM OUTPUT
        sample_frac (float): 
            subsample reads from sam file with probability sam_frac
            NOTE: Includes all mappings for a read, or none
            NOTE: Assumes multimappings for a read are consecutive in sam file

    Returns: 
       sam_fname (str): name of edited sam file with infixes before '.sam'
    """

    #NOTE sam_fname should end with .sam extension
    if sam_add_simple_map_wts and RSEM_add_map_wts:
        print "Error: Cannot add both simple_map_wts and RSEM_map_wts"
        raise ValueError
    # Get name of sam file with no directory
    sam_fname_nodir = get_fname_nodir(sam_fname)
    # Filter unmapped reads from sam file
    if filter_unmapped: 
        out_fname = expt_dir + "/process/" + sam_fname_nodir[:-4] +\
            ".mapped.sam"
        proc.sam_filter_unmapped(sam_fname, out_fname)
        sam_fname = out_fname
        sam_fname_nodir = get_fname_nodir(sam_fname)
    # add field at end of sam file with uniform mapping weights for 
    # multimapped reads
    if sam_add_simple_map_wts:
        out_fname = expt_dir + "/process/" + sam_fname_nodir[:-4] +\
            ".simple_wts.sam"
        proc.sam_add_simple_map_wts(sam_fname, out_fname)
        sam_fname = out_fname
        sam_fname_nodir = get_fname_nodir(sam_fname)
    # add field at end of sam file with RSEM mapping weight
    if RSEM_add_map_wts: 
        out_fname = expt_dir + "/process/" + sam_fname_nodir[:-4] + ".wts.sam"
        proc.RSEM_add_weights(sam_fname, out_fname)
        sam_fname = out_fname
        sam_fname_nodir = get_fname_nodir(sam_fname)
    # subsample reads from sam file with probability sam_frac
    if sample_frac:
        out_fname = expt_dir + "/process/" + sam_fname_nodir[:-4] +\
            "{0}pct.sam".format(sample_frac * 100)
        proc.sam_subsample(sam_fname, out_fname, sample_frac)
        sam_fname = out_fname
        sam_fname_nodir = get_fname_nodir(sam_fname)
    return sam_fname

def process_sam_file(
        expt_dir, sam_fname, gene_seq_fname, gene_len_fname, shift_dict, 
        cod_trunc_5p, cod_trunc_3p, min_fp_size, max_fp_size, num_tr_genes, 
        num_te_genes, min_cts_per_gene, min_cod_w_data, raw_psct=0, 
        paralog_groups_fname=False, overwrite=False, folds=False):
    """
    Processes an RP sam file for an experiment.
    Makes in expt_dir/process:
        cts_by_codon file   (sum sam map wts per codon)
        outputs file        (scaled cts_by_codon, each gene mean centered at 1)
        te_bounds file      (first and last codon idxs per test set gene)
        te_data_table file  (data table for test set codons)
        tr_bounds file      (first and last codon idxs per training set gene)
        tr_data_table file  (data table for training set codons)

    Args: 
        expt_dir (str) - name of experiment directory
        sam_fname (str) - name of input sam file
        gene_seq_fname (str) - name of transcriptome fasta file
        gene_len_fname (str) - name of gene lengths file
        shift_dict (dict): 
            {fp_size (int): 
                {frame (int): shift (int, or False) for frame in range(2)}
                for fp_size in range(min_fp_size, max_fp_size + 1) }
        cod_trunc_5p (int): number of codons to exclude at start of each CDS
        cod_trunc_3p (int): number of codons to exclude at end of each CDS
        min_fp_size (int): minimum size footprint to accept in sam file
        max_fp_size (int): maximum size footprint to accept in sam file
        num_tr_genes (int): number of genes to sort into training set
        num_te_genes (int): number of genes to sort into test set
        min_cts_per_gene (int): 
            cutoff for total cts on gene to include gene in tr/te sets
        min_cod_w_data (int):
            cutoff for codons with data to include gene in tr/te sets
        raw_psct (float): psct to add to raw cts_by_codon values
        paralog_groups_fname (str):
            file containing one group of mutually paralogous genes on each line
            results in filtering out all but one paralog (NOTE which one?)
        overwrite (bool): flag to overwrite processed files, default False
        folds (bool): 
            might implement this later to divide genes into folds
            rather than explicit training/test sets

    Returns: 
        void, makes files listed above
    """
    # Load CDS dict, len dict
    len_dict = proc.get_len_dict(gene_len_fname)
    cds_dict = proc.get_cds_dict(gene_seq_fname, len_dict)

    # Process and write cts by codon file
    cts_by_codon = proc.get_cts_by_codon(
        sam_fname, cds_dict, len_dict, shift_dict, min_fp_size, max_fp_size)
    cts_by_codon_fname = expt_dir + \
        "/process/cts_by_codon.size.{0}.{1}.txt".format(
            min_fp_size, max_fp_size)
    if not os.path.isfile(cts_by_codon_fname):
        print "making file " + cts_by_codon_fname
        proc.write_cts_by_codon(cts_by_codon_fname, cts_by_codon)
    else: 
        print "file " + cts_by_codon_fname + " already exists"

    #Process and write outputs file
    outputs = proc.get_outputs(
        cts_by_codon, cod_trunc_5p, cod_trunc_3p, raw_psct=raw_psct)
    outputs_fname = expt_dir + \
        "/process/outputs.size.{0}.{1}.txt".format(min_fp_size, max_fp_size)
    if raw_psct:
        outputs_fname = outputs_fname[:-4] +\
            ".raw_psct.{0}.txt".format(raw_psct)
    if not os.path.isfile(outputs_fname):
        print "making file " + outputs_fname
        proc.write_outputs(outputs_fname, outputs)
    else: 
        print "file " + outputs_fname + " already exists"

    #Make training and test set codon files
    tr_set_fname = expt_dir + "/process/tr_set_bounds.size." + \
        "{0}.{1}.trunc.{2}.{3}.min_cts.{4}.min_cod.{5}.top.{6}.txt".format(
            min_fp_size, max_fp_size, cod_trunc_5p, cod_trunc_3p, 
            min_cts_per_gene, min_cod_w_data, num_tr_genes + num_te_genes)
    te_set_fname = expt_dir + "/process/te_set_bounds.size." + \
        "{0}.{1}.trunc.{2}.{3}.min_cts.{4}.min_cod.{5}.top.{6}.txt".format(
            min_fp_size, max_fp_size, cod_trunc_5p, cod_trunc_3p, 
            min_cts_per_gene, min_cod_w_data, num_tr_genes + num_te_genes)
    print "making file " + tr_set_fname
    print "making file " + te_set_fname
    proc.make_codon_set_files(
        cds_dict, cts_by_codon_fname, outputs_fname, tr_set_fname, 
        te_set_fname, num_tr_genes, num_te_genes, cod_trunc_5p, cod_trunc_3p, 
        min_cts_per_gene, min_cod_w_data, 
        paralog_groups_fname=paralog_groups_fname, overwrite=overwrite)

def make_log_outputs(expt_dir, outputs_fname, scaled_psct):
    """
    Makes log transformed outputs file with a scaled pseudocount

    Args: 
        expt_dir (str) - name of experiment directory
        outputs_fname (str) - name of original outputs file
        scaled_psct (float) - scaled pseudocount to add before log transform

    Returns: 
        void, just makes log transformed outputs file w. scaled psct.
    """
    # load outputs
    outputs = proc.load_outputs(outputs_fname)
    # add scaled pseudocount, log transform outputs
    log_outputs = proc.log_transform_outputs(outputs, scaled_psct)
    # write log transformed pseudocounts
    log_out_fname = expt_dir +\
        "/process/log_outputs.scaled_psct_{0}.txt".format(scaled_psct)
    proc.write_outputs(log_outputs, log_out_fname)

def make_bin_outputs(expt_dir, outputs_fname, cutoff):
    """
    Makes binary transformed outputs file
        output = 1 if output >= cutoff else 0

    Args: 
        expt_dir (str) - name of experiment directory
        outputs_fname (str) - name of original outputs file
        cutoff (float) - scaled counts cutoff to transform to 1, else 0

    Returns: 
        void, just makes binary transformed outputs file w. scaled cts cutoff
    """
    outputs = proc.load_outputs(outputs_fname)
    bin_outputs = proc.bin_transform_outputs(outputs, cutoff)
    bin_out_fname = expt_dir + "/process/bin_outputs.cutoff_{0}.txt".format(cutoff)
    proc.write_outputs(bin_outputs, bin_out_fname)

def setup_lasagne_nn(name, expt_dir):
        #name, expt_dir, gene_seq_fname, gene_len_fname, tr_codons_fname, 
        #te_codons_fname, outputs_fname, rel_cod_idxs=False, rel_nt_idxs=False, 
        #nonlinearity="tanh", widths=[200], input_drop_rate=0, 
        #hidden_drop_rate=0, num_outputs=1, update_method="sgd", 
        #filter_max=False, filter_test=False,
        #filter_pct=False, rel_struc_idxs=False, struc_fname=False, 
        #max_struc_start_idx=None, max_struc_width=None, aa_feats=False, 
        #learning_rate=0.01, momentum=0.09, batch_size=500, raw_psct=0, 
        #lr_decay=16):
    """
    Makes initial directories for a neural network model

    Args:
        name (str) - name of neural network model
        expt_dir (str) - name of experiment directory

    Returns: 
        void, just makes directories for neural network model
    """
    # Get names of initial directories
    expt_nn_dir = expt_dir + "/lasagne_nn"
    nn_dir = expt_nn_dir + "/" + name
    init_data_dir = nn_dir + "/init_data"
    # Make initial directories
    proc.make_out_dir(expt_nn_dir)
    proc.make_out_dir(nn_dir)
    proc.make_out_dir(init_data_dir)
     
def make_lasagne_feedforward_nn(
        name, expt_dir, gene_seq_fname, gene_len_fname, tr_codons_fname,
        te_codons_fname, outputs_fname, rel_cod_idxs=[], rel_nt_idxs=[],
        nonlinearity="tanh", widths=[200], input_drop_rate=0, 
        hidden_drop_rate=0, num_outputs=1, update_method="sgd",
        filter_max = False, filter_test=False,
        filter_pct=False, rel_struc_idxs=False, struc_fname=False, 
        max_struc_start_idx=None, max_struc_width=None, aa_feats=False,
        learning_rate=0.01, lr_decay=16, 
        momentum=0.9, batch_size=500,
        log_y=False, scaled_psct=0, raw_psct=False, loss_fn="L2", 
        drop_zeros=False, nonnegative=True):
    """
    Sets up neural network model directory, 
        initializes neural network model, 
        saves initial parameters, 
        and returns neural network model

    Args:
        name (str) - name of neural network model
        expt_dir (str) - name of experiment directory
        gene_seq_fname (str) - name of transcriptome fasta file
        gene_len_fname (str) - name of gene lengths file
        tr_codons_fname (str) - name of training set codons file
        te_codons_fname (str) - name of test set codons file
        outputs_fname (str) - name of outputs file
        rel_cod_idxs (list of ints) - indices of codon features in model
        rel_nt_idxs (list of ints) - indices of nucleotide features in model
        nonlinearity (str) - name of nonlinearity fn [tanh|rectify|linear]
        widths (list of ints) - # of units in each hidden layer, in order
        input_drop_rate (float) - dropout rate for inputs
        hidden_drop_rate (float) - dropout rate for hidden unit
        num_outputs (int) - number of units in output layer
        update_method (str) - name of update method [sgd|momentum|nesterov]
        NOTE: more arg descriptions here

    Returns: 
        my_nn (lasagnenn.FeedforwardMLP) - neural network object
    """
    # Initialize neural network directories
    setup_lasagne_nn(name, expt_dir)
    
    # Load neural network data matrices
    X_tr, y_tr, X_te, y_te = proc.load_lasagne_data(
        gene_len_fname, gene_seq_fname, tr_codons_fname, te_codons_fname,
        outputs_fname, rel_cod_idxs=rel_cod_idxs,
        rel_nt_idxs=rel_nt_idxs, rel_struc_idxs=rel_struc_idxs, 
        struc_fname=struc_fname, max_struc_start_idx=max_struc_start_idx, 
        max_struc_width=max_struc_width, aa_feats=aa_feats, 
        filter_max=filter_max, filter_pct=filter_pct, filter_test=filter_test)

    # NOTE: Should I remove this?
    # Handle log transformation of y values
    if log_y:
        #Must have either a scaled psct to add, or a raw psct that has already
        #been put in the counts_by_codon when making the outputs file
        # Maybe change this scheme for raw pscts in the future?
        if scaled_psct <= 0 and not raw_psct and not drop_zeros:
            raise ValueError("Pseudocount must be >= 0 for log y")
        if scaled_psct > 0:
            y_tr = np.log(y_tr + scaled_psct)
            y_te = np.log(y_te + scaled_psct)
        if (not scaled_psct > 0) and raw_psct:
            y_tr = np.log(y_tr)
            y_te = np.log(y_te)
        if ((not scaled_psct > 0) and not raw_psct) and drop_zeros:
            positive = (y_tr > 0).ravel()
            y_tr = y_tr[positive]
            X_tr = X_tr[positive]

    # Save initial parameters
    out_dir = expt_dir + "/lasagne_nn"
    _, _, _, params = inspect.getargvalues(inspect.currentframe())
    del params["X_tr"]
    del params["X_te"]
    proc.pickle_obj(params, 
        out_dir + "/{0}/init_data/init_data.pkl".format(name))

    # Make neural network object
    my_nn = lasagnenn.FeedforwardMLP(
        X_tr, y_tr, X_te, y_te, name=name, out_dir=out_dir,
        learning_rate=learning_rate, lr_decay=lr_decay,
        update_method=update_method, widths=widths,
        nonlinearity=nonlinearity, input_drop_rate=input_drop_rate, 
        hidden_drop_rate=hidden_drop_rate, num_outputs=num_outputs, 
        momentum=momentum, batch_size=batch_size, loss_fn=loss_fn, 
        nonnegative=nonnegative)

    # Return neural network object
    return my_nn

def load_lasagne_feedforward_nn(nn_dir, epoch):
    init_data_pkl = nn_dir + "/init_data/init_data.pkl"
    params = proc.load_obj(init_data_pkl)
    #Kludge, take this out in the future
    if not params.get("max_struc_start_idx", False):
        params["max_struc_start_idx"] = None
    if not params.get("max_struc_width", False):
        params["max_struc_width"] = None
    if not params.get("aa_feats", False):
        params["aa_feats"] = False
    if not params.get("nonnegative", False):
        params["nonnegative"] = False
    X_tr, _, X_te, _ = proc.load_lasagne_data(
        params["gene_len_fname"], params["gene_seq_fname"],
        params["tr_codons_fname"], params["te_codons_fname"],
        params["outputs_fname"], rel_cod_idxs=params["rel_cod_idxs"],
        rel_nt_idxs=params["rel_nt_idxs"], 
        rel_struc_idxs=params["rel_struc_idxs"],
        struc_fname=params["struc_fname"], 
        max_struc_start_idx=params["max_struc_start_idx"], 
        max_struc_width=params["max_struc_width"], 
        aa_feats=params["aa_feats"],
        filter_pct=params["filter_pct"])

    my_nn = lasagnenn.FeedforwardMLP(
        X_tr, params["y_tr"], X_te, params["y_te"],  
        name=params["name"], out_dir=params["out_dir"],
        learning_rate=params["learning_rate"], 
        update_method=params["update_method"], widths=params["widths"],
        nonlinearity=params["nonlinearity"], 
        input_drop_rate=params["input_drop_rate"],
        hidden_drop_rate=params["hidden_drop_rate"], 
        num_outputs=params["num_outputs"],
        momentum=params["momentum"], batch_size=params["batch_size"], 
        nonnegative=params["nonnegative"],
        reloaded=True)

    my_nn.unpickle_epoch(epoch)
  
    return my_nn

def make_lasagne_adjacency_nn(
        name, expt_dir, gene_seq_fname, gene_len_fname, tr_codons_fname,
        te_codons_fname, outputs_fname, rel_cod_idxs=[], cod_adj_idxs=[], 
        rel_nt_idxs=[], nt_adj_idxs=[], nonlinearity="tanh", widths=[200], 
        input_drop_rate=0, hidden_drop_rate=0, num_outputs=1, 
        update_method="sgd", filter_pct=False, rel_struc_idxs=False, 
        struc_fname=False, max_struc_start_idx=None, 
        max_struc_width=None, aa_feats=False, learning_rate=0.01, 
        momentum=0.9, batch_size=500):

    X_tr, y_tr, X_te, y_te = setup_lasagne_nn(
        name, expt_dir, gene_seq_fname, gene_len_fname, tr_codons_fname,
        te_codons_fname, outputs_fname, rel_cod_idxs=rel_cod_idxs,
        rel_nt_idxs=rel_nt_idxs, nonlinearity=nonlinearity, widths=widths,
        input_drop_rate=input_drop_rate, hidden_drop_rate=hidden_drop_rate,
        num_outputs=num_outputs, update_method=update_method,
        filter_max=filter_max,
        filter_pct=filter_pct, rel_struc_idxs=rel_struc_idxs,
        struc_fname=struc_fname, max_struc_start_idx=max_struc_start_idx, 
        max_struc_width=max_struc_width, aa_feats=aa_feats, 
        learning_rate=learning_rate,
        momentum=momentum, batch_size=batch_size)

    out_dir = expt_dir + "/lasagne_nn"
    _, _, _, params = inspect.getargvalues(inspect.currentframe())
    proc.pickle_obj(params, 
        out_dir + "/{0}/init_data/init_data.pkl".format(name))

    my_nn = lasagnenn.AdjacencyMLP(
        X_tr, y_tr, X_te, y_te, name=name, out_dir=out_dir, 
        rel_cod_idxs=rel_cod_idxs, cod_adj_idxs=cod_adj_idxs, 
        rel_nt_idxs=rel_nt_idxs, 
        nt_adj_idxs=nt_adj_idxs, learning_rate=learning_rate,
        update_method=update_method, widths=widths, nonlinearity=nonlinearity,
        input_drop_rate=input_drop_rate, hidden_drop_rate=hidden_drop_rate,
        num_outputs=num_outputs, momentum=momentum, batch_size=batch_size)

    return my_nn

def load_lasagne_adjacency_nn(nn_dir, epoch):
    init_data_pkl = nn_dir + "/init_data/init_data.pkl"
    params = proc.load_obj(init_data_pkl)
    try: 
        max_struc_start_idx=params["max_struc_start_idx"], 
    except KeyError: 
        max_struc_start_idx=None
    try: 
        max_struc_width=params["max_struc_width"], 
    except KeyError: 
        max_struc_width=None
    X_tr, y_tr, X_te, y_te = proc.load_lasagne_data(
        params["gene_len_fname"], params["gene_seq_fname"],
        params["tr_codons_fname"], params["te_codons_fname"],
        params["outputs_fname"], rel_cod_idxs=params["rel_cod_idxs"],
        rel_nt_idxs=params["rel_nt_idxs"],
        rel_struc_idxs=params["rel_struc_idxs"],
        struc_fname=params["struc_fname"], 
        max_struc_start_idx=max_struc_start_idx,
        max_struc_width=max_struc_width, 
        filter_pct=params["filter_pct"])

    my_nn = lasagnenn.AdjacencyMLP(
        X_tr, y_tr, X_te, y_te, name=params["name"], out_dir=params["out_dir"],
        rel_cod_idxs=params["rel_cod_idxs"], 
        cod_adj_idxs=params["cod_adj_idxs"],
        rel_nt_idxs=params["rel_nt_idxs"],
        nt_adj_idxs=params["nt_adj_idxs"],
        learning_rate=params["learning_rate"],
        update_method=params["update_method"], widths=params["widths"],
        nonlinearity=params["nonlinearity"],
        input_drop_rate=params["input_drop_rate"],
        hidden_drop_rate=params["hidden_drop_rate"],
        num_outputs=params["num_outputs"],
        momentum=params["momentum"], batch_size=params["batch_size"],
        reloaded=True)

    my_nn.unpickle_epoch(epoch)

    return my_nn

def make_lasagne_split_nn(
        name, expt_dir, gene_seq_fname, gene_len_fname, tr_codons_fname,
        te_codons_fname, outputs_fname, rel_cod_idxs, rel_nt_idxs,
        nonlinearity="tanh", widths=[200], input_drop_rate=0,
        hidden_drop_rate=0, num_outputs=1, update_method="sgd",
        filter_pct=False, rel_struc_idxs=False, struc_fname=False,
        max_struc_start_idx=None, max_struc_width=None,
        learning_rate=0.01, momentum=0.9, batch_size=500):

    X_tr, y_tr, X_te, y_te = setup_lasagne_nn(
        name, expt_dir, gene_seq_fname, gene_len_fname, tr_codons_fname,
        te_codons_fname, outputs_fname, rel_cod_idxs=rel_cod_idxs,
        rel_nt_idxs=rel_nt_idxs, nonlinearity=nonlinearity, widths=widths,
        input_drop_rate=input_drop_rate, hidden_drop_rate=hidden_drop_rate,
        num_outputs=num_outputs, update_method=update_method,
        filter_max=filter_max,
        filter_pct=filter_pct, rel_struc_idxs=rel_struc_idxs,
        struc_fname=struc_fname, max_struc_start_idx=max_struc_start_idx, 
        max_struc_width=max_struc_width, learning_rate=learning_rate,
        momentum=momentum, batch_size=batch_size)

    out_dir = expt_dir + "/lasagne_nn"
    _, _, _, params = inspect.getargvalues(inspect.currentframe())
    proc.pickle_obj(params, 
        out_dir + "/{0}/init_data/init_data.pkl".format(name))

    my_nn = lasagnenn.SplitMLP(
        X_tr, y_tr, X_te, y_te,  
        rel_cod_idxs, rel_nt_idxs, name=name, out_dir=out_dir, 
        learning_rate=learning_rate, update_method=update_method, 
        widths=widths, nonlinearity=nonlinearity,
        input_drop_rate=input_drop_rate, hidden_drop_rate=hidden_drop_rate,
        num_outputs=num_outputs, momentum=momentum, batch_size=batch_size)

    return my_nn

def load_lasagne_split_nn(nn_dir, epoch):
    init_data_pkl = nn_dir + "/init_data/init_data.pkl"
    params = proc.load_obj(init_data_pkl)
    X_tr, y_tr, X_te, y_te = proc.load_lasagne_data(
        params["gene_len_fname"], params["gene_seq_fname"],
        params["tr_codons_fname"], params["te_codons_fname"],
        params["outputs_fname"], rel_cod_idxs=params["rel_cod_idxs"],
        rel_nt_idxs=params["rel_nt_idxs"],
        rel_struc_idxs=params["rel_struc_idxs"],
        struc_fname=params["struc_fname"], 
        max_struc_start_idx=params["max_struc_start_idx"], 
        max_struc_width=params["max_struc_width"], 
        filter_pct=params["filter_pct"])

    my_nn = lasagnenn.AdjacencyMLP(
        X_tr, y_tr, X_te, y_te, name=params["name"], out_dir=params["out_dir"],
        rel_cod_idxs=params["rel_cod_idxs"], 
        rel_nt_idxs=params["rel_nt_idxs"],
        learning_rate=params["learning_rate"],
        update_method=params["update_method"], widths=params["widths"],
        nonlinearity=params["nonlinearity"],
        input_drop_rate=params["input_drop_rate"],
        hidden_drop_rate=params["hidden_drop_rate"],
        num_outputs=params["num_outputs"],
        momentum=params["momentum"], batch_size=params["batch_size"],
        reloaded=True)

    my_nn.unpickle_epoch(epoch)

    return my_nn

def plot_lasagne_nn(
        expt_dir, name, epoch, xlims, ylims, textpos, rel_cod_idxs=False,
        rel_nt_idxs=False, rel_struc_idxs=False, cost_by_epoch_plots=True,
        weight_plots=True, scat_plots=True, cum_err_plots=True, 
        overwrite=False, ext="pdf", scat_bins=0):
    out_dir = expt_dir + "/lasagne_nn"
    plot.lasagne_make_all_plots(
        out_dir, name, epoch, xlims, ylims, textpos, rel_cod_idxs=rel_cod_idxs,
        rel_nt_idxs=rel_nt_idxs, rel_struc_idxs=rel_struc_idxs,
        cost_by_epoch_plots=cost_by_epoch_plots, weight_plots=weight_plots,
        scat_plots=scat_plots, cum_err_plots=cum_err_plots, 
        overwrite=overwrite, ext=ext, scat_bins=scat_bins
        )

def make_linreg(
        expt_dir, name, gene_seq_fname, gene_len_fname, tr_codons_fname, 
        te_codons_fname, outputs_fname, rel_cod_idxs=False, rel_nt_idxs=False, 
        rel_struc_idxs=False, struc_fname=False):
    linreg.linreg_init(
        expt_dir, name, gene_seq_fname, gene_len_fname, tr_codons_fname, 
        te_codons_fname, outputs_fname, rel_cod_idxs=rel_cod_idxs, 
        rel_nt_idxs=rel_nt_idxs)
    len_dict = proc.get_len_dict(gene_len_fname)
    cds_dict = proc.get_cds_dict(gene_seq_fname, len_dict)
    if struc_fname:
        struc_dict = proc.get_struc_dict(struc_fname)
    else:
        struc_dict = False
    X_tr, X_te, y_tr, y_te = proc.get_data_matrices(
        cds_dict, tr_codons_fname, te_codons_fname, outputs_fname, 
        rel_cod_idxs=rel_cod_idxs, rel_nt_idxs=rel_nt_idxs, 
        rel_struc_idxs=rel_struc_idxs, struc_dict=struc_dict)
    wts, y_tr_hat, y_te_hat = linreg.train_and_predict(X_tr, X_te, y_tr, y_te)
    linreg_fname = "{0}/linreg/{1}".format(expt_dir, name)
    proc.pickle_obj(wts, "{0}/wts.pkl".format(linreg_fname))
    proc.pickle_obj(y_tr, "{0}/y_tr.pkl".format(linreg_fname))
    proc.pickle_obj(y_tr_hat, "{0}/y_tr_hat.pkl".format(linreg_fname))
    proc.pickle_obj(y_te, "{0}/y_te.pkl".format(linreg_fname))
    proc.pickle_obj(y_te_hat, "{0}/y_te_hat.pkl".format(linreg_fname))
    return wts, y_tr_hat, y_te_hat, y_tr, y_te

def plot_linreg(
        expt_dir, name, y_tr, y_tr_hat, y_te, y_te_hat, scat_xlim, scat_ylim, 
        scat_textpos, wts=False, rel_cod_idxs=False, rel_nt_idxs=False, 
        rel_struc_idxs=False, weight_plots=True, scat_plots=True, 
        cum_err_plots=True, ext="pdf"):
    linreg_dir = "{0}/linreg/{1}".format(expt_dir, name)
    plot_dir = linreg_dir
    if scat_plots:
        plot.make_density_scatter_plot(
            y_tr.flatten(), y_tr_hat.flatten(), 
            name + " Training Predicted vs. True Outputs", "True Outputs", 
            "Predicted Outputs", scat_xlim, scat_ylim, scat_textpos, 
            "{0}/tr_scatter.{1}".format(linreg_dir, ext))
        plot.make_density_scatter_plot(
            y_te.flatten(), y_te_hat.flatten(), 
            name + " Test Predicted vs. True Outputs", "True Outputs", 
            "Predicted Outputs", scat_xlim, scat_ylim, scat_textpos, 
            "{0}/te_scatter.{1}".format(linreg_dir, ext))
    if weight_plots:
        if type(wts) != bool:
            title_prefix = name + " Input Weight Sq. Norm by"
            fname_prefix = "{0}/ssweights".format(linreg_dir)
            plot.make_weight_norms_plots(
                wts.transpose(), title_prefix, fname_prefix, 
                rel_cod_idxs=rel_cod_idxs, rel_nt_idxs=rel_nt_idxs, 
                rel_struc_idxs=rel_struc_idxs)
    if cum_err_plots:
        plot.make_binned_err_plot(
            y_tr.flatten(), y_tr_hat.flatten(),
            name + "Training Sq. Err. v. Sorted Outputs",
            "{0}/tr_bin_err.{1}".format(plot_dir, ext), 700)
        plot.make_smoothed_err_plot(
            y_tr.flatten(), y_tr_hat.flatten(),
            name + "Training Sq. Err. v. Sorted Outputs",
            "{0}/tr_smoothed_err.{1}".format(plot_dir, ext))
        plot.make_cum_err_plot(
            y_tr.flatten(), y_tr_hat.flatten(),
            name + "Training Cumulative Sq. Err. v. Sorted Outputs",
            "{0}/tr_cum_err.{1}".format(plot_dir, ext))
        plot.make_err_thresh_plot(
            y_tr.flatten(), y_tr_hat.flatten(),
            name + " Training fraction of points below error threshold",
            "{0}/tr_err_thresh.{1}".format(plot_dir, ext))
        plot.make_binned_err_plot(
            y_te.flatten(), y_te_hat.flatten(),
            name + "Test Sq. Err. v. Sorted Outputs",
            "{0}/te_bin_err.{1}".format(plot_dir, ext), 700)
        plot.make_smoothed_err_plot(
            y_te.flatten(), y_te_hat.flatten(),
            name + "Test Sq. Err. v. Sorted Outputs",
            "{0}/te_smoothed_err.{1}".format(plot_dir, ext))
        plot.make_cum_err_plot(
            y_te.flatten(), y_te_hat.flatten(),
            name + "Test Cumulative Sq. Err. v. Sorted Outputs",
            "{0}/te_cum_err.{1}".format(plot_dir, ext))
        plot.make_err_thresh_plot(
            y_te.flatten(), y_te_hat.flatten(),
            name + " Test fraction of points below error threshold",
            "{0}/te_err_thresh.{1}".format(plot_dir, ext))

def get_linreg_optimal_codons(linreg_dir, aa_seq, maximum=False):
    linreg_data_fname = linreg_dir + "/linreg_data.txt"
    d = proc.load_linreg_data_file(linreg_data_fname)
    rel_cod_idxs = d["rel_cod_idxs"]
    wts = proc.load_obj(linreg_dir + "/wts.pkl")
    opt_total_seq, opt_vit_score = opt.get_optimal_codons_linreg(aa_seq, wts, rel_cod_idxs, maximum=maximum)
    return opt_total_seq, opt_vit_score

def get_lasagne_optimal_codons(
        nn_dir, epoch, aa_seq, nn_type="Feedforward", maximum=False, 
        nt_feats=False):
    init_data_pkl = nn_dir + "/init_data/init_data.pkl"
    params = proc.load_obj(init_data_pkl)
    if nn_type=="Feedforward":
        my_nn = load_lasagne_feedforward_nn(nn_dir, epoch)
    elif nn_type=="Adjacency":
        my_nn = load_lasagne_adjacency_nn(nn_dir, epoch)
    elif nn_type=="Split":
        my_nn = load_lasagne_split_nn(nn_dir, epoch)
    else: print "nn_type must be in [Feedforward, Adjacency, Split]"
    rel_cod_idxs = params["rel_cod_idxs"]
    opt_total_seq, opt_vit_score = opt.get_optimal_codons_lasagne(
        aa_seq, my_nn, rel_cod_idxs, maximum=maximum, nt_feats=nt_feats)
    return opt_total_seq, opt_vit_score
    
def save_lasagne_optimal_codons(
        nn_dir, epoch, aa_seq, aa_seq_name, nn_type="Feedforward", 
        nt_feats=False):
    max_total_seq, max_vit_score = get_lasagne_optimal_codons(
        nn_dir, epoch, aa_seq, nn_type=nn_type, maximum=True, nt_feats=nt_feats)
    min_total_seq, min_vit_score = get_lasagne_optimal_codons(
        nn_dir, epoch, aa_seq, nn_type=nn_type, maximum=False, 
        nt_feats=nt_feats)
    opt_dir = nn_dir + "/optimize"
    opt.make_opt_dir(opt_dir)
    f = open(opt_dir + "/{0}.txt".format(aa_seq_name), "w")
    f.write("max seq: {0}\n".format(max_total_seq))
    f.write("max score: {0}\n".format(max_vit_score))
    f.write("min seq: {0}\n".format(min_total_seq))
    f.write("min score: {0}\n".format(min_vit_score))
    f.close()

def get_protein_score_dist(nn_dir, epoch, aa_seq, num_samples, nt_feats=False):
    init_data_pkl = nn_dir + "/init_data/init_data.pkl"
    params = proc.load_obj(init_data_pkl)
    rel_cod_idxs = params["rel_cod_idxs"]
    my_nn = load_lasagne_feedforward_nn(nn_dir, epoch)
    scores_sorted, cod_seqs_sorted = opt.get_score_dist(
        aa_seq, my_nn, rel_cod_idxs, num_samples, nt_feats=nt_feats)
    return scores_sorted, cod_seqs_sorted

def make_structure_fasta(out_fname, gene_seq_fname, gene_len_fname, window_len):
    seq_dict = proc.get_seq_dict(gene_seq_fname)
    len_dict = proc.get_len_dict(gene_len_fname)
    proc.make_structure_fasta(out_fname, seq_dict, len_dict, window_len)
  
def plot_nn_preds_by_gene(
        nn, te_set_data_table_fname, gene2sym_fname, plot_dir, 
        gene2sym_mod_fn=False):
    #In case you need to modify the names of your genes in gene2sym_fname
        #e.g. adding _13cds10 to gene names
    if gene2sym_mod_fn:
        gene2sym = proc.get_gene2symbol_dict(gene2sym_fname, gene2sym_mod_fn)
    #And if not...
    else:
        gene2sym = proc.get_gene2symbol_dict(gene2sym_fname)
    #Get data table for test set
    te_set_data_table = proc.load_set_data_table(te_set_data_table_fname)
    # And now plot predictions by gene!
    plot.plot_preds_by_gene_2(
        nn, te_set_data_table, gene2sym, plot_dir)

def plot_cts_and_outputs(
        expt_name, expt_dir, cts_by_codon_fname, outputs_fname, 
        process_label, tr_codons_fname, cts_by_codon_xlims=False, 
        outputs_xlims=False, log_plot=False):
    """
        process_label: string that denotes parameters of the processing used
            for tr_codons_fname, including fp size, # genes, and type/# of psct
    """

    #Load cts_by_codon and outputs
    cts_by_codon = proc.load_cts_by_codon(cts_by_codon_fname)
    outputs = proc.load_outputs(outputs_fname)
    #Load tr_codon_set
    tr_codon_bounds = proc.load_codon_set_bounds(tr_codons_fname)
    tr_codon_set = proc.expand_codon_set(tr_codon_bounds)
    #load training counts and outputs
    cts = proc.get_y(tr_codon_set, cts_by_codon).astype("int")
    outputs = proc.get_y(tr_codon_set, outputs)
    #Make plots
    xlab = "FP counts"
    ylab = "# Codons"
    title = "{0} Counts by Codon Frequencies".format(expt_name)
    xlims = cts_by_codon_xlims
    out_fname = expt_dir +\
        "/plots/cts_by_codon_freqs.{0}.pdf".format(process_label)
    plot.make_cts_by_codon_freq_plot(
        cts, title, xlab, ylab, out_fname, xlims=xlims)
    xlab = "Scaled counts"
    ylab = "# Codons"
    title = "{0} Scaled Counts Frequencies".format(expt_name)
    xlims = outputs_xlims
    out_fname = expt_dir + "/plots/output_freqs.{0}.pdf".format(process_label)
    plot.make_outputs_freq_plot(
        outputs, title, xlab, ylab, out_fname, xlims=xlims, num_bins=500)
    if log_plot:
        xlab = "Log scaled counts"
        ylab = "# Codons"
        title = "{0} Log Scaled Counts Frequencies".format(expt_name)
        out_fname = expt_dir +\
            "/plots/log_output_freqs.{0}.pdf".format(process_label)
        plot.make_log_outputs_freq_plot(
            outputs, title, xlab, ylab, out_fname, num_bins=500)
