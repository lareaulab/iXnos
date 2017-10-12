# iXnos

iXnos (pronounced ick-noss) is a software package with two main functions: 
* Train a neural network model to predict local mRNA to protein translation rates, across a set of transcripts. 
* Using a trained model, design the fastest or slowest predicted coding sequence for a gene.

## Installation

iXnos is supported for Linux/OSX environments

### Dependencies

iXnos is written in python version 2.7, and requires the following python packages: 
* lasagne
* theano
* matplotlib
* sklearn

These packages can be installed with: 
```
pip install [options] <package>
```

The current released version of Lasagne (v. 0.1) is not compatible with the current released version of theano (v. 0.9.0), so please install the development version of lasagne with: 
```
pip install [options] https://github.com/Lasagne/Lasagne/archive/master.zip
```

iXnos also depends on the following software. The binaries for these packages should be accessible on your PATH:
* Bowtie 1
* RSEM
* samtools
* SRA toolkit

### Installing iXnos

You can install iXnos by cloning our git repository: 
```
git clone https://github.com/lareaulab/iXnos.git
```

And then add iXnos to your PYTHONPATH by appending the following line to your config file of choice (e.g. .profile, .bash_rc, etc.)
```
export PYTHONPATH="${PYTHONPATH}:/path/above/repo/iXnos"
```

And update your PYTHONPATH in your shell
```
source <config_filename>
```

## Getting Started

### Creating a New Experiment

To use iXnos with a new ribosome profiling data set, first create experiment directories with: 
```
python /path/above/repo/iXnos/reproduce_scripts/make_expt_dir.py /path/above/repo/iXnos/expts <expt_name>
```

Do any preprocessing steps necessary for your data set (e.g., trimming adapters). 

### Read mapping

We provide genome and mapping files for S. cerevisiae, and humans, in iXnos/genome_data. We map fastq data to a transcriptome using bowtie 1, and obtain mapping weights for multimapping reads with RSEM (recommended, but not required).

For yeast data, an example workflow is provided in iXnos/expts/weinberg/process/weinberg.sh. 
* filter out rRNA reads by mapping against ScerRRNA bowtie1 index with options bowtie -v 2 -p 36 -S 
* filter out noncoding RNA reads by mapping against rna_coding bowtie1 index with options bowtie -v 2 -p 36 -S
* map to yeast transcriptome with bowtie1 index scer.transcripts.13cds10 and options bowtie -a --norc -v 2 -p 36 -S
* (recommended) obtain mapping weights with rsem-calculate-expression on our bowtie output and an RSEM index built on scer.transcripts.13cds10, then convert bam output to sam format with samtools view -h

A similar example is provided for human data in iXnos/expts/iwasaki/process/iwasaki.sh

You can use your own rRNA and noncoding RNA files to filter unwanted reads, but you need to make your own bowtie indices. 
You can use your own transcriptome fasta file, but you need to create your own bowtie and RSEM indices. In addition, you need to create your own gene_lengths_file in iXnos/genome_data with one line for each gene in your transcriptome, in the format: 

<fasta_gene_name> <5' UTR length (nt)> <CDS length (nt)> <3' UTR length (nt)>

Our transcriptomes are padded with 13 nt on the 5' end of each CDS, and 10 nt on the 3' end. 
Length files are provided for our transcriptomes in: 
* iXnos/genome_data/scer.transcripts.13cds10.lengths.txt (S.cer)
* iXnos/genome_data/gencode.v22.transcript.13cds10.lengths.txt (human)

Your final sam file should go in iXnos/expts/<expt_name>/process/your_sam_file.sam

### Processing mapped data

First, we append a mapping weights field to the end of our sam file:
```
import iXnos.interface as inter
expt_dir = "/path/above/repo/iXnos/expts/<expt_name>"
sam_fname = expt_dir + "/process/your_sam_file.sam"
inter.edit_sam_file(
    expt_dir, sam_fname, filter_unmapped=True, 
    sam_add_simple_map_wts=False, RSEM_add_map_wts=True)
```
If you did not use RSEM to get mapping weights, change sam_add_simple_map_wts to True, and RSEM_add_map_wts to False. This will assign uniform weights to all mappings for multimapping reads. By default, RSEM includes unmapped reads in the output file. If you have already filtered these out, you can set filter_unmapped to False. 

After this step, your new sam file will be located in .../<expt_dir>/process, and will end in ".wts.sam". Use this sam file for future commands.

Next, we choose A site assignment rules for our reads. We start by plotting the frames of the 5' ends of our footprints, for each footprint size. 
```
# gene lengths file for yeast, or your own file with a custom transcriptome
gene_len_fname = /path/to/repo/iXnos/genome_data/scer.transcripts.13cds10.lengths.txt
inter.do_frame_and_size_analysis(
    expt_dir, sam_fname, gene_len_fname, min_size=13, 
    max_size=36, verbose=False)
```
You can find the output plots in .../<expt_dir>/plots. 

Use these plots to define A site offset rules. Generally, we find that 28mers map consistently to the 0th frame, and their A site offset is 15 nt from the 5' end. You can assign additional offsets by looking for adjacent footprint lengths that map consistently to 1-2 frames, and reasoning about the over/underdigestion that led to this class of reads. For example, a 29mer mapping to the 0th frame is underdigested by 1 nt at the 3' end relative to a 28mer in the 0th frame, so its A site offset is still 15 nt from the 5' end. A 29mer in the 2nd frame is underdigested by 1 nt at its 5' end, so its A site offset is 16 nt from the 5' end. We do this for all footprint sizes with large amounts of data in our experiment. A histogram of overall counts by size is also provided in <expt_dir>/plots.

Finally, we process our sam file for use in model training. A site offset rules are used to define a shift_dict. 

You can pass an optional paralog_groups_file if you wish to only allow one member from each group of annotated paralogs. This file should be tab delimited, with one set of paralogous gene names on each line. See an example in .../iXnos/genome_data/scer_100_80.paralogs.id.txt
```
# A site offset rules
shift_dict = {
    27:{0:False, 1:14, 2:False}, 28:{0:15, 1:False, 2:False},
    29:{0:15, 1:False, 2:16}, 30:{0:15, 1:False, 2:16},
    31:{0:15, 1:False, 2:16}}
# Number of codons to ignore at either end of CDS in model training
cod_trunc_5p = 20
cod_trunc_3p = 20
# Minimum and maximum allowed footprint sizes
min_fp_size = 27
max_fp_size = 31
# Choose number of genes in training and test sets
num_tr_genes = 333
num_te_genes = 167
# Set coverage cutoffs for genes in training/test sets
# Minimum footprint counts per gene
min_cts_per_gene = 200
# Minimum A site codons with nonzero counts
min_cod_w_data = 100
# Transcriptome file, e.g. for yeast
/path/to/repo/iXnos/genome_data/scer.transcripts.13cds10.fa
# Optional paralog groups file, else pass False
paralog_groups_fname = /path/to/repo/iXnos/genome_data/scer_100_80.paralogs.id.txt
inter.process_sam_file(
    expt_dir, sam_fname, gene_seq_fname, gene_len_fname,
    shift_dict, cod_trunc_5p, cod_trunc_3p, min_fp_size,
    max_fp_size, num_tr_genes, num_te_genes, min_cts_per_gene,
    min_cod_w_data, paralog_groups_fname=paralog_groups_fname,
    overwrite=False)
```
This will populate all the necessary files in .../<expt_dir>/process to run a neural network model.

### Training a Neural Network Model

An example for training a neural network model on your processed data is provided below:
```
import iXnos.interface as inter
name = <your_nn_name>
expt_dir = "/path/to/repo/iXnos/expts/<expt_name>"
# Transcriptome file
gene_seq_fname = "/path/to/repo/iXnos/genome_data/<transcriptome_fname>"
# Gene lengths file
gene_len_fname = "/path/to/repo/iXnos/genome_data/<gene_len_fname>"
# Training/test set genes files, generated in last step
# Look up in <expt_dir>/process
tr_codons_fname = expt_dir + "/process/tr_set_bounds.size.<min_fp_size>.<max_fp_size>.trunc.<cod_trunc_5p>.<cod_trunc_3p>.min_cts.<min_cts_per_gene>.min_cod.<min_cod_w_data>.top.<num_genes_te+tr>.txt
te_codons_fname = expt_dir + "/process/te_set_bounds.size.<min_fp_size>.<max_fp_size>.trunc.<cod_trunc_5p>.<cod_trunc_3p>.min_cts.<min_cts_per_gene>.min_cod.<min_cod_w_data>.top.<num_genes_te+tr>.txt
# Indices of codons and nts in desired sequence neighborhood
rel_cod_idxs = range(-5,5)
rel_nt_idxs = range(-15,15)
# Number of hidden units per layer [1st, 2nd, ...]
# Number of elements specifies model depth
widths = [200]
# Specify activation function on hidden layers, [tanh|rectify|linear]
nonlinearity = "tanh"
# Specify lasagne update method, [sgd|momentum|nesterov]
update_method = "nesterov"
# Create neural network
neural_net = inter.make_lasagne_feedforward_nn(
    name, expt_dir, gene_seq_fname, gene_len_fname, tr_codons_fname,
    te_codons_fname, outputs_fname, rel_cod_idxs=rel_cod_idxs,
    rel_nt_idxs=rel_nt_idxs, nonlinearity=nonlinearity, widths=widths, 
    update_method=update_method)
# Run neural network for a desired Number of epochs
num_epochs = 30
neural_net.run_epochs(num_epochs)
``` 
The model will print out training and test errors for each epoch. We recommend terminating model training when the test error plateaus. 

The state of the model at each epoch is saved in .../iXnos/expts/<expt_dir>/lasagne_nn/<model_name>/epoch<#> To reload a model at a given epoch, call:
```
import iXnos.interface as inter
nn_dir = "/path/to/repo/iXnos/expts/<expt_dir>/lasagne_nn/<model_name>"
epoch = <epoch_number>
neural_net = inter.load_lasagne_feedforward_nn(nn_dir, epoch)
```

iXnos will not allow you to overwrite a model or an epoch, so if you want to retrain a model with the same name, or overwrite an epoch, you will have to delete the corresponding directories in .../<expt_dir>/lasagne_nn.

### Neural Networks With mRNA Structure Features

Instructions will be added later

### Optimizing Gene Sequences Under a Neural Network Model

After you have trained a neural network model, you can use iXnos to find the minimum and maximum translation speed coding sequence for a given protein, under that model. 

We currently support optimization for models that are only trained with sequence input features (i.e. no structure features). The sequence input features must satisfy one of these conditions: 
* Codon features only (set nt_feats=False)
* Nucleotide features spanning the same sequence neighborhood as the codon features (set nt_feats=True)
```
import iXnos.interface as inter
# Specify the directory of your trained neural network
nn_dir = /path/to/repo/iXnos/expts/<expt_name>/lasagne_nn/<model_name>
# Specify the final epoch of your model training
epoch = <int>
# The amino acid sequence of your protein for optimization, as a string
aa_seq = "MACDEFGHIKLNPQRSTVWY"
# Say we trained a neural network model with only codon features
nt_feats = False
# By default, maximum=False, yielding the minimum translation time coding sequence. 
# If you want the maximum translation time coding sequence, set maximum=True
min_seq, min_score = inter.get_lasagne_optimal_codons(
        nn_dir, epoch, aa_seq, maximum=False, nt_feats=nt_feats)
max_seq, max_score = inter.get_lasagne_optimal_codons(
        nn_dir, epoch, aa_seq, maximum=True, nt_feats=nt_feats)
```

min_seq and max_seq are the sequences of the minimum and maximum translation time coding sequence of aa_seq, under your model. 

min_score and max_score are the total translation time score of these sequences under our model, which is equal to the sum of the predicted scaled counts under the model at all codons in the sequence.

## Reproducing Research

You can re-run most of the analyses in our research paper with the makefile provided in the root iXnos directory. There is some variability in the results, due to randomness in parameter initialization and model training.

We train a large number of models for each experiment, so we recommend running make for the experiments separately, with multithreading (make -j #threads) enabled. 
```
make -j 20 weinberg_expt
make -j 10 lareau_expt
make -j 10 iwasaki_expt
make -j 10 green_expt
```
After all experiment analyses have been run, you can recreate the figures and values cited in the publication with:
```
make figures
make paper_data
```
The results, figures, and values generated can be found in iXnos/results/[weinberg|lareau|iwasaki|green|figures|paper_data]
