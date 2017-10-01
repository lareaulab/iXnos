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
export PYTHONPATH="/path/above/repo/iXnos/"
```

And update your PYTHONPATH in your shell
```
source <config_filename>
```

## Reproducible Research

You can re-run most of the analyses in our research paper with the makefile provided in the root iXnos directory. There is some variability in the results, due to randomness in parameter initialization and model training.

We train a large number of models for each experiment, so we recommend running make for the experiments separately, with multithreading (make -j #threads) enabled. 
```
make -j 20 weinberg_expt
make -j 10 lareau_expt
make -j 10 iwasaki_expt
make -j 10 green_expt
make figures
make paper_data
```

The results, figures, and values presented in the paper can be found in iXnos/results/[weinberg|lareau|iwasaki|green|figures|paper_data]

Currently only weinberg_expt, iwasaki_expt, and green_expt can be run in the makefile, because we have not released our lareau_expt data yet.
