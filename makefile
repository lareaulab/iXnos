
### Top level installation variable definitions

parent_dir = $(shell pwd)

top_dir = $(parent_dir)

lib_dir_name = iXnos
lib_dir = $(top_dir)/$(lib_dir_name)
lib_files = \
	$(lib_dir)/interface.py \
	$(lib_dir)/lasagnenn.py \
	$(lib_dir)/linreg.py \
	$(lib_dir)/optimizecodons.py \
	$(lib_dir)/plot.py \
	$(lib_dir)/process.py \
	$(lib_dir)/featureanalysis.py

genome_dir_name = genome_data
genome_dir = $(top_dir)/$(genome_dir_name)
yeast_gene_len_file = $(genome_dir)/scer.transcripts.13cds10.lengths.txt
yeast_gene_seq_file = $(genome_dir)/scer.transcripts.13cds10.fa
yeast_paralogs_file = $(genome_dir)/scer_100_80.paralogs.id.txt
yeast_gene_symbol_file = $(genome_dir)/scer_id_symbol.txt
yeast_codon_props_file = $(genome_dir)/yeast_codon_properties.txt
yeast_codon_anticodon_file = $(genome_dir)/yeast_codon_anticodon.csv
human_gene_len_file = $(genome_dir)/gencode.v22.transcript.13cds10.lengths.txt
human_gene_seq_file = $(genome_dir)/gencode.v22.transcript.13cds10.fa

genome_files = \
	$(yeast_gene_len_file) \
	$(yeast_gene_seq_file) \
	$(yeast_paralogs_file) \
	$(yeast_gene_symbol_file) \
	$(yeast_codon_props_file) \
	$(human_gene_len_file) \
	$(human_gene_seq_file) 

struc_dir_name = structure_data
struc_dir = $(top_dir)/$(struc_dir_name)
yeast_30_windows_seqs_file_zip = $(struc_dir)/scer.13cds10.windows.30len.fa.gz
yeast_30_windows_str_scores_file_zip = $(struc_dir)/scer.13cds10.windows.30len.fold.gz
yeast_30_windows_seqs_file = $(struc_dir)/scer.13cds10.windows.30len.fa
yeast_30_windows_str_scores_file = $(struc_dir)/scer.13cds10.windows.30len.fold
struc_files = \
	$(yeast_30_windows_seqs_file) \
	$(yeast_30_windows_str_scores_file) 

wetlab_dir_name = wetlab_data
wetlab_dir = $(top_dir)/$(wetlab_dir_name)
circligase_qpcr_file = $(wetlab_dir)/circligase_qpcr.csv
mrna_data_file = $(wetlab_dir)/mrna_cy0_summary_20170906.csv
facs_data_zip_file = $(wetlab_dir)/gated-facs-data-20170829.csv.gz
facs_data_file = $(wetlab_dir)/gated-facs-data-20170829.csv

wetlab_files = \
	$(circligase_qpcr_file) \
	$(mrna_data_file) \
	$(facs_data_zip_file) \
	$(facs_data_file) 

repro_dir_name = reproduce_scripts
repro_dir = $(top_dir)/$(repro_dir_name)
repro_files = \
	$(repro_dir)/process_data.py \
	$(repro_dir)/feat_neighborhood_nn_series.py \
	$(repro_dir)/feat_neighborhood_linreg_series.py \
	$(repro_dir)/leaveout_series.py \
	$(repro_dir)/struc_series.py \
	$(repro_dir)/28mer_models.py \
	$(repro_dir)/codon_scores.py \
	$(repro_dir)/aggregate_mses.py \
	$(repro_dir)/aggregate_linreg_mses.py \
	$(repro_dir)/plot_nn.py \
	$(repro_dir)/plot_genes.py \
	$(repro_dir)/paper_data.py \
	$(repro_dir)/optimize_cds.py \
	$(repro_dir)/pkl2txt.py \
	$(repro_dir)/figure_1B_scaledcts.R \
	$(repro_dir)/figure_1C_mse.R \
	$(repro_dir)/figure_1D_scatter.R \
	$(repro_dir)/figure_1E_indiv_gene.R \
	$(repro_dir)/figure_1F_binned_error.R \
	$(repro_dir)/figure_2A_codonmse.R \
	$(repro_dir)/figure_2B_heatmap.R \
	$(repro_dir)/figure_2C_tai.R \
	$(repro_dir)/figure_2D_wobble.R \
	$(repro_dir)/figure_2E_lareaucodonmse.R \
	$(repro_dir)/figure_2F_5prime.R \
	$(repro_dir)/figure_2G_asite.R \
	$(repro_dir)/figure_2H_cl2.R \
	$(repro_dir)/figure_2I_cl1.R \
	$(repro_dir)/figure_3B_citrine_dist.R \
	$(repro_dir)/figure_3C_facs.R \
	$(repro_dir)/figure_3D_te.R \
	$(repro_dir)/supp_table_codon_scores.py \
	$(repro_dir)/supp_figure_facs.R \
	$(repro_dir)/supp_figure_mrna_qpcr.R \
	$(repro_dir)/supp_figure_greenmse.R \
	$(repro_dir)/supp_figure_iwasakimse.R \
	$(repro_dir)/supp_figure_cl1v2.R

expts_dir_name = expts
expts_dir = $(top_dir)/$(expts_dir_name)

expt_subdirs = process plots lasagne_nn linreg

results_dir_name = results
results_dir = $(top_dir)/$(results_dir_name)

fig_dir = $(results_dir)/figures

citrine_aa_seq = MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLGYGLMCFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK

### General experiment variable definitions

fp_sizes = $(shell seq -s ' ' 13 36)
10_epochs = $(shell seq -s ' ' 0 10)
15_epochs = $(shell seq -s ' ' 0 15)
20_epochs = $(shell seq -s ' ' 0 20)
25_epochs = $(shell seq -s ' ' 0 25)
30_epochs = $(shell seq -s ' ' 0 30)
35_epochs = $(shell seq -s ' ' 0 35)
nn_init_data_file_names = init_data.pkl y_te.pkl y_tr.pkl
nn_init_data_files = $(addprefix init_data/,$(nn_init_data_file_names))
nn_epoch_file_names = te_cost_by_epoch.pkl tr_cost_by_epoch.pkl weights.pkl \
                      y_te_hat.pkl  y_tr_hat.pkl
nn_35_epoch_files = \
	$(nn_init_data_files) \
	$(foreach num,$(35_epochs),\
		$(addprefix epoch$(num)/,$(nn_epoch_file_names)))
nn_30_epoch_files = \
	$(nn_init_data_files) \
	$(foreach num,$(30_epochs),\
		$(addprefix epoch$(num)/,$(nn_epoch_file_names)))
nn_25_epoch_files = \
	$(nn_init_data_files) \
	$(foreach num,$(25_epochs),\
		$(addprefix epoch$(num)/,$(nn_epoch_file_names)))
nn_20_epoch_files = \
	$(nn_init_data_files) \
	$(foreach num,$(20_epochs),\
		$(addprefix epoch$(num)/,$(nn_epoch_file_names)))
nn_15_epoch_files = \
	$(nn_init_data_files) \
	$(foreach num,$(15_epochs),\
		$(addprefix epoch$(num)/,$(nn_epoch_file_names)))
nn_10_epoch_files = \
	$(nn_init_data_files) \
	$(foreach num,$(10_epochs),\
		$(addprefix epoch$(num)/,$(nn_epoch_file_names)))

plot_files = cts_by_size.pdf \
	${addprefix size_,${addsuffix _by_frame.pdf,${fp_sizes}}}

plot_files_pattern = cts_by_size%pdf \
	${addprefix size_,${addsuffix _by_frame%pdf,${fp_sizes}}}

10_reps = $(shell seq -s ' ' 0 9)
feat_neighborhoods = \
	cod_p0 cod_p0_nt_p0p2 cod_n3p2 cod_n3p2_nt_n9p8 \
	cod_n5p4 cod_n5p4_nt_n15p14 cod_n7p5 cod_n7p5_nt_n21p17
leaveout_cods = $(shell seq -s ' ' -7 5)

linreg_data_files = \
	linreg_data.txt wts.pkl \
	y_te.pkl y_te_hat.pkl \
	y_tr.pkl y_tr_hat.pkl

linreg_plot_files = \
	tr_scatter.pdf te_scatter.pdf \
	tr_bin_err.pdf te_bin_err.pdf \
	tr_cum_err.pdf te_cum_err.pdf

nn_plot_files = \
	tr_cost_by_epoch.pdf te_cost_by_epoch.pdf \
	tr_scatter.pdf te_scatter.pdf \
	tr_bin_err.pdf te_bin_err.pdf \
	tr_cum_err.pdf te_cum_err.pdf

### Weinberg variable definitions

weinberg_expt_dir_name = weinberg
weinberg_expt_dir = $(expts_dir)/$(weinberg_expt_dir_name)
weinberg_proc_dir = $(weinberg_expt_dir)/process
weinberg_plot_dir = $(weinberg_expt_dir)/plots
weinberg_nn_dir = $(weinberg_expt_dir)/lasagne_nn
weinberg_lr_dir = $(weinberg_expt_dir)/linreg
weinberg_subdirs  = $(addprefix $(weinberg_expt_dir)/,$(expt_subdirs))
weinberg_subdirs_pattern = $(addprefix $(weinberg_expt_dir)%,$(expt_subdirs))

weinberg_gene_len_file = $(genome_dir)/yeast_13cds10_lengths.txt
weinberg_gene_seq_file = $(genome_dir)/yeast_13cds10.fa

weinberg_raw_fastq_file = $(weinberg_proc_dir)/SRR1049521.fastq 
weinberg_trimmer_script = $(weinberg_proc_dir)/trim_linker_weinberg.pl
weinberg_mapping_script = $(weinberg_proc_dir)/weinberg.sh

weinberg_raw_sam_file = $(weinberg_expt_dir)/process/weinberg.transcript.sam
weinberg_mapped_sam_file = $(weinberg_expt_dir)/process/weinberg.transcript.mapped.sam
weinberg_sam_file = $(weinberg_expt_dir)/process/weinberg.transcript.mapped.wts.sam


weinberg_plot_files = ${addprefix ${weinberg_plot_dir}/,${plot_files}}

weinberg_plot_files_pattern = \
	${addprefix ${weinberg_plot_dir}/,${plot_files_pattern}}

weinberg_27_31_cts_by_codon = $(weinberg_proc_dir)/cts_by_codon.size.27.31.txt
weinberg_27_31_outputs = $(weinberg_proc_dir)/outputs.size.27.31.txt
weinberg_27_31_te_bounds = $(weinberg_proc_dir)/te_set_bounds.size.27.31.trunc.20.20.min_cts.200.min_cod.100.top.500.txt
weinberg_27_31_te_data_table = $(weinberg_proc_dir)/te_set_bounds.size.27.31.trunc.20.20.min_cts.200.min_cod.100.top.500.data_table.txt
weinberg_27_31_tr_bounds = $(weinberg_proc_dir)/tr_set_bounds.size.27.31.trunc.20.20.min_cts.200.min_cod.100.top.500.txt
weinberg_27_31_tr_data_table = $(weinberg_proc_dir)/tr_set_bounds.size.27.31.trunc.20.20.min_cts.200.min_cod.100.top.500.data_table.txt

weinberg_27_31_proc_sam_files = \
	$(weinberg_27_31_cts_by_codon) \
	$(weinberg_27_31_outputs) \
	$(weinberg_27_31_te_bounds) \
	$(weinberg_27_31_te_data_table) \
	$(weinberg_27_31_tr_bounds) \
	$(weinberg_27_31_tr_data_table)

weinberg_27_31_proc_sam_pattern = \
	$(subst .txt,%txt,$(weinberg_27_31_cts_by_codon)) \
	$(subst .txt,%txt,$(weinberg_27_31_outputs)) \
	$(subst .txt,%txt,$(weinberg_27_31_te_bounds)) \
	$(subst .txt,%txt,$(weinberg_27_31_te_data_table)) \
	$(subst .txt,%txt,$(weinberg_27_31_tr_bounds)) \
	$(subst .txt,%txt,$(weinberg_27_31_tr_data_table))

weinberg_28_cts_by_codon = $(weinberg_proc_dir)/cts_by_codon.size.28.28.txt
weinberg_28_outputs = $(weinberg_proc_dir)/outputs.size.28.28.txt
weinberg_28_te_bounds = $(weinberg_proc_dir)/te_set_bounds.size.28.28.trunc.20.20.min_cts.200.min_cod.100.top.500.txt
weinberg_28_te_data_table = $(weinberg_proc_dir)/te_set_bounds.size.28.28.trunc.20.20.min_cts.200.min_cod.100.top.500.data_table.txt
weinberg_28_tr_bounds = $(weinberg_proc_dir)/tr_set_bounds.size.28.28.trunc.20.20.min_cts.200.min_cod.100.top.500.txt
weinberg_28_tr_data_table = $(weinberg_proc_dir)/tr_set_bounds.size.28.28.trunc.20.20.min_cts.200.min_cod.100.top.500.data_table.txt

weinberg_28_proc_sam_files = \
	$(weinberg_28_cts_by_codon) \
	$(weinberg_28_outputs) \
	$(weinberg_28_te_bounds) \
	$(weinberg_28_te_data_table) \
	$(weinberg_28_tr_bounds) \
	$(weinberg_28_tr_data_table)

weinberg_28_proc_sam_pattern = \
	$(subst .txt,%txt,$(weinberg_28_cts_by_codon)) \
	$(subst .txt,%txt,$(weinberg_28_outputs)) \
	$(subst .txt,%txt,$(weinberg_28_te_bounds)) \
	$(subst .txt,%txt,$(weinberg_28_te_data_table)) \
	$(subst .txt,%txt,$(weinberg_28_tr_bounds)) \
	$(subst .txt,%txt,$(weinberg_28_tr_data_table))

weinberg_feat_nb_series_files = $(foreach feat_nb,$(feat_neighborhoods),$(foreach rep,$(10_reps),$(addprefix $(weinberg_nn_dir)/full_$(feat_nb)_rep$(rep)/,$(nn_30_epoch_files))))
weinberg_feat_nb_series_pattern = $(addprefix $(weinberg_nn_dir)/full_%/,$(nn_30_epoch_files))

weinberg_leaveout_series_files = $(foreach lo_cod,$(leaveout_cods),$(foreach rep,$(10_reps),$(addprefix $(weinberg_nn_dir)/nocod$(lo_cod)_cod_n7p5_nt_n21p17_rep$(rep)/,$(nn_30_epoch_files))))
weinberg_leaveout_series_pattern = $(addprefix $(weinberg_nn_dir)/nocod%/,$(nn_30_epoch_files))

weinberg_28mer_files = $(addprefix $(weinberg_nn_dir)/s28_cod_n7p5_nt_n21p17/,$(nn_25_epoch_files))
weinberg_28mer_pattern = $(addprefix $(weinberg_nn_dir)/s28_%/,$(nn_25_epoch_files))

weinberg_fp_struc_files = $(foreach rep,$(10_reps),$(addprefix $(weinberg_nn_dir)/str_n17n15_cod_n7p5_nt_n21p17_rep$(rep)/,$(nn_30_epoch_files)))
weinberg_fp_struc_pattern = $(addprefix $(weinberg_nn_dir)/str_%/,$(nn_30_epoch_files))

weinberg_max_struc_files = $(foreach rep,$(10_reps),$(addprefix $(weinberg_nn_dir)/max_str_p13p42_cod_n7p5_nt_n21p17_rep$(rep)/,$(nn_35_epoch_files)))
weinberg_max_struc_pattern = $(addprefix $(weinberg_nn_dir)/max_str_%/,$(nn_35_epoch_files))

weinberg_linreg_series_files = $(foreach feat_nb,$(feat_neighborhoods),$(addprefix $(weinberg_lr_dir)/lr_$(feat_nb)/,$(linreg_data_files)))
weinberg_linreg_series_pattern = $(addprefix $(weinberg_lr_dir)/lr_%/,$(linreg_data_files)) 

weinberg_struc_plot_files = $(addprefix $(weinberg_nn_dir)/str_n17n15_cod_n7p5_nt_n21p17_rep0/plots/,$(nn_plot_files))
weinberg_struc_plot_pattern = $(addprefix $(weinberg_nn_dir)/str_n17n15_cod_n7p5_nt_n21p17_rep0/plot%/,$(nn_plot_files))

weinberg_full_plot_files = $(addprefix $(weinberg_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/plots/,$(nn_plot_files))
weinberg_full_plot_pattern = $(addprefix $(weinberg_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/plot%/,$(nn_plot_files))

weinberg_28_plot_files = $(addprefix $(weinberg_nn_dir)/s28_cod_n7p5_nt_n21p17/plots/,$(nn_plot_files))
weinberg_28_plot_pattern = $(addprefix $(weinberg_nn_dir)/s28_cod_n7p5_nt_n21p17/plot%/,$(nn_plot_files))

weinberg_struc_gene_plot_dir = $(weinberg_nn_dir)/str_n17n15_cod_n7p5_nt_n21p17_rep0/gene_plots

# Weinberg results files variable definitions

weinberg_results_dir = $(results_dir)/weinberg
weinberg_results_feat_neighborhood_dir = \
	$(weinberg_results_dir)/feat_neighborhood_series
weinberg_results_leaveout_dir = $(weinberg_results_dir)/leaveout_series
weinberg_results_struc_series_dir = $(weinberg_results_dir)/structure_series
weinberg_results_struc_dir = $(weinberg_results_dir)/str_n17n15_cod_n7p5_nt_n21p17_rep0
weinberg_results_struc_epoch_dir = $(weinberg_results_struc_dir)/epoch30
weinberg_results_full_dir = $(weinberg_results_dir)/full_cod_n7p5_nt_n21p17_rep0
weinberg_results_full_epoch_dir = $(weinberg_results_full_dir)/epoch30
weinberg_results_28_dir = $(weinberg_results_dir)/s28_cod_n7p5_nt_n21p17
weinberg_results_28_epoch_dir = $(weinberg_results_28_dir)/epoch25
weinberg_results_opt_model_dir = $(weinberg_results_dir)/full_cod_n3p2_nt_n9p8_rep0
weinberg_results_opt_model_epoch_dir = $(weinberg_results_opt_model_dir)/epoch30

weinberg_full_analysis_epoch_dir = $(weinberg_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/epoch30
weinberg_28_analysis_epoch_dir = $(weinberg_nn_dir)/s28_cod_n7p5_nt_n21p17/epoch25

weinberg_full_codon_scores_results_files = \
	$(weinberg_results_full_epoch_dir)/codon_scores.tsv \
	$(weinberg_results_full_epoch_dir)/codon_scores_colormap.pdf

weinberg_full_codon_scores_results_files_pattern = \
	$(weinberg_results_full_epoch_dir)/codon_scores%tsv \
	$(weinberg_results_full_epoch_dir)/codon_scores_colormap%pdf

weinberg_28_codon_scores_results_files = \
	$(weinberg_results_28_epoch_dir)/codon_scores.tsv \
	$(weinberg_results_28_epoch_dir)/codon_scores_colormap.pdf

weinberg_28_codon_scores_results_files_pattern = \
	$(weinberg_results_28_epoch_dir)/codon_scores%tsv \
	$(weinberg_results_28_epoch_dir)/codon_scores_colormap%pdf

weinberg_feat_nb_mses_file = \
	$(weinberg_results_feat_neighborhood_dir)/feat_neighborhood_mses.txt
weinberg_feat_nb_linreg_mses_file = \
	$(weinberg_results_feat_neighborhood_dir)/linreg_mses.txt
weinberg_leaveout_mses_file = \
	$(weinberg_results_leaveout_dir)/leaveout_mses.txt
weinberg_struc_mses_file = \
	$(weinberg_results_struc_series_dir)/struc_mses.txt

weinberg_final_model_y_te = $(weinberg_nn_dir)/str_n17n15_cod_n7p5_nt_n21p17_rep0/init_data/y_te.txt
weinberg_final_model_y_te_hat = $(weinberg_nn_dir)/str_n17n15_cod_n7p5_nt_n21p17_rep0/epoch30/y_te_hat.txt \

weinberg_results_final_model_y_te = $(weinberg_results_struc_epoch_dir)/y_te.txt
weinberg_results_final_model_y_te_hat = $(weinberg_results_struc_epoch_dir)/y_te_hat.txt

weinberg_opt_series_file = $(weinberg_results_opt_model_epoch_dir)/citrine_opt_series.txt

### Lareau variable definitions

lareau_expt_dir_name = lareau
lareau_expt_dir = $(expts_dir)/$(lareau_expt_dir_name)
lareau_proc_dir = $(lareau_expt_dir)/process
lareau_plot_dir = $(lareau_expt_dir)/plots
lareau_nn_dir = $(lareau_expt_dir)/lasagne_nn
lareau_lr_dir = $(lareau_expt_dir)/linreg
lareau_subdirs  = $(addprefix $(lareau_expt_dir)/,$(expt_subdirs))
lareau_subdirs_pattern = $(addprefix $(lareau_expt_dir)%,$(expt_subdirs))

lareau_gene_len_file = $(genome_dir)/yeast_13cds10_lengths.txt
lareau_gene_seq_file = $(genome_dir)/yeast_13cds10.fa

lareau_raw_fastq_file_1 = $(lareau_proc_dir)/LLMG004_S31_L007_R1_001.fastq
lareau_raw_fastq_file_2 = $(lareau_proc_dir)/LLMG005_S1_L001_R1_001.fastq

lareau_raw_fastq_files = \
	$(lareau_proc_dir)/LLMG004_S31_L007_R1_001.fastq \
	$(lareau_proc_dir)/LLMG005_S1_L001_R1_001.fastq 
lareau_trimmer_script = $(lareau_proc_dir)/trim_linker_bc.pl
lareau_mapping_script = $(lareau_proc_dir)/lareau.sh

lareau_raw_sam_file = $(lareau_expt_dir)/process/lareau.transcript.sam
lareau_mapped_sam_file = $(lareau_expt_dir)/process/lareau.transcript.mapped.sam
lareau_sam_file = $(lareau_expt_dir)/process/lareau.transcript.mapped.wts.sam

lareau_plot_files = ${addprefix ${lareau_plot_dir}/,${plot_files}}

lareau_plot_files_pattern = \
	${addprefix ${lareau_plot_dir}/,${plot_files_pattern}}

lareau_27_29_cts_by_codon = $(lareau_proc_dir)/cts_by_codon.size.27.29.txt
lareau_27_29_outputs = $(lareau_proc_dir)/outputs.size.27.29.txt
lareau_27_29_te_bounds = $(lareau_proc_dir)/te_set_bounds.size.27.29.trunc.20.20.min_cts.200.min_cod.100.top.500.txt
lareau_27_29_te_data_table = $(lareau_proc_dir)/te_set_bounds.size.27.29.trunc.20.20.min_cts.200.min_cod.100.top.500.data_table.txt
lareau_27_29_tr_bounds = $(lareau_proc_dir)/tr_set_bounds.size.27.29.trunc.20.20.min_cts.200.min_cod.100.top.500.txt
lareau_27_29_tr_data_table = $(lareau_proc_dir)/tr_set_bounds.size.27.29.trunc.20.20.min_cts.200.min_cod.100.top.500.data_table.txt

lareau_27_29_proc_sam_files = \
        $(lareau_27_29_cts_by_codon) \
        $(lareau_27_29_outputs) \
        $(lareau_27_29_te_bounds) \
        $(lareau_27_29_te_data_table) \
        $(lareau_27_29_tr_bounds) \
        $(lareau_27_29_tr_data_table)

lareau_27_29_proc_sam_pattern = \
        $(subst .txt,%txt,$(lareau_27_29_cts_by_codon)) \
        $(subst .txt,%txt,$(lareau_27_29_outputs)) \
        $(subst .txt,%txt,$(lareau_27_29_te_bounds)) \
        $(subst .txt,%txt,$(lareau_27_29_te_data_table)) \
        $(subst .txt,%txt,$(lareau_27_29_tr_bounds)) \
        $(subst .txt,%txt,$(lareau_27_29_tr_data_table))

lareau_28_cts_by_codon = $(lareau_proc_dir)/cts_by_codon.size.28.28.txt
lareau_28_outputs = $(lareau_proc_dir)/outputs.size.28.28.txt
lareau_28_te_bounds = $(lareau_proc_dir)/te_set_bounds.size.28.28.trunc.20.20.min_cts.200.min_cod.100.top.500.txt
lareau_28_te_data_table = $(lareau_proc_dir)/te_set_bounds.size.28.28.trunc.20.20.min_cts.200.min_cod.100.top.500.data_table.txt
lareau_28_tr_bounds = $(lareau_proc_dir)/tr_set_bounds.size.28.28.trunc.20.20.min_cts.200.min_cod.100.top.500.txt
lareau_28_tr_data_table = $(lareau_proc_dir)/tr_set_bounds.size.28.28.trunc.20.20.min_cts.200.min_cod.100.top.500.data_table.txt

lareau_28_proc_sam_files = \
	$(lareau_28_cts_by_codon) \
	$(lareau_28_outputs) \
	$(lareau_28_te_bounds) \
	$(lareau_28_te_data_table) \
	$(lareau_28_tr_bounds) \
	$(lareau_28_tr_data_table)

lareau_28_proc_sam_pattern = \
	$(subst .txt,%txt,$(lareau_28_cts_by_codon)) \
	$(subst .txt,%txt,$(lareau_28_outputs)) \
	$(subst .txt,%txt,$(lareau_28_te_bounds)) \
	$(subst .txt,%txt,$(lareau_28_te_data_table)) \
	$(subst .txt,%txt,$(lareau_28_tr_bounds)) \
	$(subst .txt,%txt,$(lareau_28_tr_data_table))

lareau_full_model_files = \
	$(foreach rep,$(10_reps),$(addprefix \
	$(lareau_nn_dir)/full_cod_n7p5_nt_n21p17_rep$(rep)/,\
	$(nn_15_epoch_files)))
lareau_full_model_pattern = $(addprefix $(lareau_nn_dir)/full_%/,\
	$(nn_15_epoch_files))

lareau_leaveout_series_files = \
	$(foreach lo_cod,$(leaveout_cods),\
	$(foreach rep,$(10_reps), $(addprefix \
	$(lareau_nn_dir)/nocod$(lo_cod)_cod_n7p5_nt_n21p17_rep$(rep)/,\
	$(nn_15_epoch_files))))
lareau_leaveout_series_pattern = $(addprefix $(lareau_nn_dir)/nocod%/,\
	$(nn_15_epoch_files))

lareau_28mer_files = $(addprefix $(lareau_nn_dir)/s28_cod_n7p5_nt_n21p17/,\
	$(nn_10_epoch_files))
lareau_28mer_pattern = $(addprefix $(lareau_nn_dir)/s28_%/,$(nn_10_epoch_files))

lareau_full_plot_files = $(addprefix $(lareau_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/plots/,$(nn_plot_files))
lareau_full_plot_pattern = $(addprefix $(lareau_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/plot%/,$(nn_plot_files))

lareau_28_plot_files = $(addprefix $(lareau_nn_dir)/s28_cod_n7p5_nt_n21p17/plots/,$(nn_plot_files))
lareau_28_plot_pattern = $(addprefix $(lareau_nn_dir)/s28_cod_n7p5_nt_n21p17/plot%/,$(nn_plot_files))

# Lareau results files variable definitions

lareau_results_dir = $(results_dir)/lareau
lareau_results_leaveout_dir = $(lareau_results_dir)/leaveout_series
lareau_results_full_dir = $(lareau_results_dir)/full_cod_n7p5_nt_n21p17_rep0
lareau_results_full_epoch_dir = $(lareau_results_full_dir)/epoch15
lareau_results_28_dir = $(lareau_results_dir)/s28_cod_n7p5_nt_n21p17
lareau_results_28_epoch_dir = $(lareau_results_28_dir)/epoch10

lareau_full_analysis_epoch_dir = $(lareau_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/epoch15
lareau_28_analysis_epoch_dir = $(lareau_nn_dir)/s28_cod_n7p5_nt_n21p17/epoch10

lareau_full_codon_scores_results_files = \
	$(lareau_results_full_epoch_dir)/codon_scores.tsv \
	$(lareau_results_full_epoch_dir)/codon_scores_colormap.pdf

lareau_full_codon_scores_results_files_pattern = \
	$(lareau_results_full_epoch_dir)/codon_scores%tsv \
	$(lareau_results_full_epoch_dir)/codon_scores_colormap%pdf

lareau_28_codon_scores_results_files = \
        $(lareau_results_28_epoch_dir)/codon_scores.tsv \
        $(lareau_results_28_epoch_dir)/codon_scores_colormap.pdf

lareau_28_codon_scores_results_files_pattern = \
        $(lareau_results_28_epoch_dir)/codon_scores%tsv \
        $(lareau_results_28_epoch_dir)/codon_scores_colormap%pdf

lareau_leaveout_mses_file = \
	$(lareau_results_leaveout_dir)/leaveout_mses.txt

### Iwasaki variable definitions

iwasaki_expt_dir_name = iwasaki
iwasaki_expt_dir = $(expts_dir)/$(iwasaki_expt_dir_name)
iwasaki_proc_dir = $(iwasaki_expt_dir)/process
iwasaki_plot_dir = $(iwasaki_expt_dir)/plots
iwasaki_nn_dir = $(iwasaki_expt_dir)/lasagne_nn
iwasaki_lr_dir = $(iwasaki_expt_dir)/linreg
iwasaki_subdirs  = $(addprefix $(iwasaki_expt_dir)/,$(expt_subdirs))
iwasaki_subdirs_pattern = $(addprefix $(iwasaki_expt_dir)%,$(expt_subdirs))

iwasaki_gene_len_file = $(genome_dir)/gencode.v22.transcript.13cds10.lengths.txt
iwasaki_gene_seq_file = $(genome_dir)/gencode.v22.transcript.13cds10.fa

iwasaki_raw_fastq_files = \
	$(iwasaki_proc_dir)/SRR2075925.fastq \
	$(iwasaki_proc_dir)/SRR2075926.fastq 
iwasaki_trimmer_script = $(iwasaki_proc_dir)/trim_linker.pl
iwasaki_mapping_script = $(iwasaki_proc_dir)/iwasaki.sh

iwasaki_raw_sam_file = $(iwasaki_expt_dir)/process/iwasaki.transcript.sam
iwasaki_mapped_sam_file = $(iwasaki_expt_dir)/process/iwasaki.transcript.mapped.sam
iwasaki_sam_file = $(iwasaki_expt_dir)/process/iwasaki.transcript.mapped.wts.sam

iwasaki_plot_files = ${addprefix ${iwasaki_plot_dir}/,${plot_files}}

iwasaki_plot_files_pattern = \
	${addprefix ${iwasaki_plot_dir}/,${plot_files_pattern}}

iwasaki_27_30_cts_by_codon = $(iwasaki_proc_dir)/cts_by_codon.size.27.30.txt
iwasaki_27_30_outputs = $(iwasaki_proc_dir)/outputs.size.27.30.txt
iwasaki_27_30_te_bounds = $(iwasaki_proc_dir)/te_set_bounds.size.27.30.trunc.20.20.min_cts.200.min_cod.100.top.500.txt
iwasaki_27_30_te_data_table = $(iwasaki_proc_dir)/te_set_bounds.size.27.30.trunc.20.20.min_cts.200.min_cod.100.top.500.data_table.txt
iwasaki_27_30_tr_bounds = $(iwasaki_proc_dir)/tr_set_bounds.size.27.30.trunc.20.20.min_cts.200.min_cod.100.top.500.txt
iwasaki_27_30_tr_data_table = $(iwasaki_proc_dir)/tr_set_bounds.size.27.30.trunc.20.20.min_cts.200.min_cod.100.top.500.data_table.txt

iwasaki_27_30_proc_sam_files = \
        $(iwasaki_27_30_cts_by_codon) \
        $(iwasaki_27_30_outputs) \
        $(iwasaki_27_30_te_bounds) \
        $(iwasaki_27_30_te_data_table) \
        $(iwasaki_27_30_tr_bounds) \
        $(iwasaki_27_30_tr_data_table)

iwasaki_27_30_proc_sam_pattern = \
        $(subst .txt,%txt,$(iwasaki_27_30_cts_by_codon)) \
        $(subst .txt,%txt,$(iwasaki_27_30_outputs)) \
        $(subst .txt,%txt,$(iwasaki_27_30_te_bounds)) \
        $(subst .txt,%txt,$(iwasaki_27_30_te_data_table)) \
        $(subst .txt,%txt,$(iwasaki_27_30_tr_bounds)) \
        $(subst .txt,%txt,$(iwasaki_27_30_tr_data_table))

iwasaki_28_cts_by_codon = $(iwasaki_proc_dir)/cts_by_codon.size.28.28.txt
iwasaki_28_outputs = $(iwasaki_proc_dir)/outputs.size.28.28.txt
iwasaki_28_te_bounds = $(iwasaki_proc_dir)/te_set_bounds.size.28.28.trunc.20.20.min_cts.200.min_cod.100.top.300.txt
iwasaki_28_te_data_table = $(iwasaki_proc_dir)/te_set_bounds.size.28.28.trunc.20.20.min_cts.200.min_cod.100.top.300.data_table.txt
iwasaki_28_tr_bounds = $(iwasaki_proc_dir)/tr_set_bounds.size.28.28.trunc.20.20.min_cts.200.min_cod.100.top.300.txt
iwasaki_28_tr_data_table = $(iwasaki_proc_dir)/tr_set_bounds.size.28.28.trunc.20.20.min_cts.200.min_cod.100.top.300.data_table.txt

iwasaki_28_proc_sam_files = \
	$(iwasaki_28_cts_by_codon) \
	$(iwasaki_28_outputs) \
	$(iwasaki_28_te_bounds) \
	$(iwasaki_28_te_data_table) \
	$(iwasaki_28_tr_bounds) \
	$(iwasaki_28_tr_data_table)

iwasaki_28_proc_sam_pattern = \
	$(subst .txt,%txt,$(iwasaki_28_cts_by_codon)) \
	$(subst .txt,%txt,$(iwasaki_28_outputs)) \
	$(subst .txt,%txt,$(iwasaki_28_te_bounds)) \
	$(subst .txt,%txt,$(iwasaki_28_te_data_table)) \
	$(subst .txt,%txt,$(iwasaki_28_tr_bounds)) \
	$(subst .txt,%txt,$(iwasaki_28_tr_data_table))

leaveout_cods = $(shell seq -s ' ' -7 5)

iwasaki_full_model_files = $(foreach rep,$(10_reps),$(addprefix $(iwasaki_nn_dir)/full_cod_n7p5_nt_n21p17_rep$(rep)/,$(nn_10_epoch_files)))
iwasaki_full_model_pattern = $(addprefix $(iwasaki_nn_dir)/full_%/,$(nn_10_epoch_files))

iwasaki_leaveout_series_files = \
	$(foreach lo_cod,$(leaveout_cods),\
	$(foreach rep,$(10_reps), $(addprefix \
	$(iwasaki_nn_dir)/nocod$(lo_cod)_cod_n7p5_nt_n21p17_rep$(rep)/,\
	$(nn_10_epoch_files))))

iwasaki_leaveout_series_pattern = $(addprefix $(iwasaki_nn_dir)/nocod%/,\
	$(nn_10_epoch_files))

iwasaki_28mer_files = $(addprefix $(iwasaki_nn_dir)/s28_cod_n7p5_nt_n21p17/,$(nn_10_epoch_files))
iwasaki_28mer_pattern = $(addprefix $(iwasaki_nn_dir)/s28_%/,$(nn_10_epoch_files))

iwasaki_full_plot_files = $(addprefix $(iwasaki_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/plots/,$(nn_plot_files))
iwasaki_full_plot_pattern = $(addprefix $(iwasaki_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/plot%/,$(nn_plot_files))

iwasaki_28_plot_files = $(addprefix $(iwasaki_nn_dir)/s28_cod_n7p5_nt_n21p17/plots/,$(nn_plot_files))
iwasaki_28_plot_pattern = $(addprefix $(iwasaki_nn_dir)/s28_cod_n7p5_nt_n21p17/plot%/,$(nn_plot_files))

# Iwasaki results files variable definitions

iwasaki_results_dir = $(results_dir)/iwasaki
iwasaki_results_leaveout_dir = $(iwasaki_results_dir)/leaveout_series
iwasaki_results_full_dir = $(iwasaki_results_dir)/full_cod_n7p5_nt_n21p17_rep0
iwasaki_results_full_epoch_dir = $(iwasaki_results_full_dir)/epoch10
iwasaki_results_28_dir = $(iwasaki_results_dir)/s28_cod_n7p5_nt_n21p17
iwasaki_results_28_epoch_dir = $(iwasaki_results_28_dir)/epoch10

iwasaki_full_analysis_epoch_dir = $(iwasaki_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/epoch10
iwasaki_28_analysis_epoch_dir = $(iwasaki_nn_dir)/s28_cod_n7p5_nt_n21p17/epoch10

iwasaki_full_codon_scores_results_files = \
	$(iwasaki_results_full_epoch_dir)/codon_scores.tsv \
	$(iwasaki_results_full_epoch_dir)/codon_scores_colormap.pdf

iwasaki_full_codon_scores_results_files_pattern = \
	$(iwasaki_results_full_epoch_dir)/codon_scores%tsv \
	$(iwasaki_results_full_epoch_dir)/codon_scores_colormap%pdf

iwasaki_28_codon_scores_results_files = \
        $(iwasaki_results_28_epoch_dir)/codon_scores.tsv \
        $(iwasaki_results_28_epoch_dir)/codon_scores_colormap.pdf

iwasaki_28_codon_scores_results_files_pattern = \
        $(iwasaki_results_28_epoch_dir)/codon_scores%tsv \
        $(iwasaki_results_28_epoch_dir)/codon_scores_colormap%pdf

iwasaki_leaveout_mses_file = \
	$(iwasaki_results_leaveout_dir)/leaveout_mses.txt

### Green variable definitions

green_expt_dir_name = green
green_expt_dir = $(expts_dir)/$(green_expt_dir_name)
green_proc_dir = $(green_expt_dir)/process
green_plot_dir = $(green_expt_dir)/plots
green_nn_dir = $(green_expt_dir)/lasagne_nn
green_lr_dir = $(green_expt_dir)/linreg
green_subdirs  = $(addprefix $(green_expt_dir)/,$(expt_subdirs))
green_subdirs_pattern = $(addprefix $(green_expt_dir)%,$(expt_subdirs))

green_gene_len_file = $(genome_dir)/yeast_13cds10_lengths.txt
green_gene_seq_file = $(genome_dir)/yeast_13cds10.fa

green_raw_fastq_files = \
	$(green_proc_dir)/SRR5008134.fastq \
	$(green_proc_dir)/SRR5008135.fastq 
green_trimmer_script = $(green_proc_dir)/trim_linker_green.pl
green_mapping_script = $(green_proc_dir)/green.sh

green_raw_sam_file = $(green_expt_dir)/process/green.transcript.sam
green_mapped_sam_file = $(green_expt_dir)/process/green.transcript.mapped.sam
green_sam_file = $(green_expt_dir)/process/green.transcript.mapped.wts.sam

green_plot_files = ${addprefix ${green_plot_dir}/,${plot_files}}

green_plot_files_pattern = \
	${addprefix ${green_plot_dir}/,${plot_files_pattern}}

green_27_29_cts_by_codon = $(green_proc_dir)/cts_by_codon.size.27.29.txt
green_27_29_outputs = $(green_proc_dir)/outputs.size.27.29.txt
green_27_29_te_bounds = $(green_proc_dir)/te_set_bounds.size.27.29.trunc.20.20.min_cts.200.min_cod.100.top.500.txt
green_27_29_te_data_table = $(green_proc_dir)/te_set_bounds.size.27.29.trunc.20.20.min_cts.200.min_cod.100.top.500.data_table.txt
green_27_29_tr_bounds = $(green_proc_dir)/tr_set_bounds.size.27.29.trunc.20.20.min_cts.200.min_cod.100.top.500.txt
green_27_29_tr_data_table = $(green_proc_dir)/tr_set_bounds.size.27.29.trunc.20.20.min_cts.200.min_cod.100.top.500.data_table.txt

green_27_29_proc_sam_files = \
        $(green_27_29_cts_by_codon) \
        $(green_27_29_outputs) \
        $(green_27_29_te_bounds) \
        $(green_27_29_te_data_table) \
        $(green_27_29_tr_bounds) \
        $(green_27_29_tr_data_table)

green_27_29_proc_sam_pattern = \
        $(subst .txt,%txt,$(green_27_29_cts_by_codon)) \
        $(subst .txt,%txt,$(green_27_29_outputs)) \
        $(subst .txt,%txt,$(green_27_29_te_bounds)) \
        $(subst .txt,%txt,$(green_27_29_te_data_table)) \
        $(subst .txt,%txt,$(green_27_29_tr_bounds)) \
        $(subst .txt,%txt,$(green_27_29_tr_data_table))

green_28_cts_by_codon = $(green_proc_dir)/cts_by_codon.size.28.28.txt
green_28_outputs = $(green_proc_dir)/outputs.size.28.28.txt
green_28_te_bounds = $(green_proc_dir)/te_set_bounds.size.28.28.trunc.20.20.min_cts.200.min_cod.100.top.500.txt
green_28_te_data_table = $(green_proc_dir)/te_set_bounds.size.28.28.trunc.20.20.min_cts.200.min_cod.100.top.500.data_table.txt
green_28_tr_bounds = $(green_proc_dir)/tr_set_bounds.size.28.28.trunc.20.20.min_cts.200.min_cod.100.top.500.txt
green_28_tr_data_table = $(green_proc_dir)/tr_set_bounds.size.28.28.trunc.20.20.min_cts.200.min_cod.100.top.500.data_table.txt

green_28_proc_sam_files = \
	$(green_28_cts_by_codon) \
	$(green_28_outputs) \
	$(green_28_te_bounds) \
	$(green_28_te_data_table) \
	$(green_28_tr_bounds) \
	$(green_28_tr_data_table)

green_28_proc_sam_pattern = \
	$(subst .txt,%txt,$(green_28_cts_by_codon)) \
	$(subst .txt,%txt,$(green_28_outputs)) \
	$(subst .txt,%txt,$(green_28_te_bounds)) \
	$(subst .txt,%txt,$(green_28_te_data_table)) \
	$(subst .txt,%txt,$(green_28_tr_bounds)) \
	$(subst .txt,%txt,$(green_28_tr_data_table))

green_full_model_files = \
	$(foreach rep,$(10_reps),$(addprefix \
	$(green_nn_dir)/full_cod_n7p5_nt_n21p17_rep$(rep)/,\
	$(nn_20_epoch_files)))
green_full_model_pattern = $(addprefix $(green_nn_dir)/full_%/,\
	$(nn_20_epoch_files))

green_leaveout_series_files = \
	$(foreach lo_cod,$(leaveout_cods),\
	$(foreach rep,$(10_reps), $(addprefix \
	$(green_nn_dir)/nocod$(lo_cod)_cod_n7p5_nt_n21p17_rep$(rep)/,\
	$(nn_20_epoch_files))))
green_leaveout_series_pattern = $(addprefix $(green_nn_dir)/nocod%/,\
	$(nn_20_epoch_files))

green_28mer_files = $(addprefix $(green_nn_dir)/s28_cod_n7p5_nt_n21p17/,\
	$(nn_10_epoch_files))
green_28mer_pattern = $(addprefix $(green_nn_dir)/s28_%/,$(nn_10_epoch_files))

green_full_plot_files = $(addprefix $(green_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/plots/,$(nn_plot_files))
green_full_plot_pattern = $(addprefix $(green_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/plot%/,$(nn_plot_files))

green_28_plot_files = $(addprefix $(green_nn_dir)/s28_cod_n7p5_nt_n21p17/plots/,$(nn_plot_files))
green_28_plot_pattern = $(addprefix $(green_nn_dir)/s28_cod_n7p5_nt_n21p17/plot%/,$(nn_plot_files))

# Green results files variable definitions

green_results_dir = $(results_dir)/green
green_results_leaveout_dir = $(green_results_dir)/leaveout_series
green_results_full_dir = $(green_results_dir)/full_cod_n7p5_nt_n21p17_rep0
green_results_full_epoch_dir = $(green_results_full_dir)/epoch20
green_results_28_dir = $(green_results_dir)/s28_cod_n7p5_nt_n21p17
green_results_28_epoch_dir = $(green_results_28_dir)/epoch10

green_full_analysis_epoch_dir = $(green_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/epoch20
green_28_analysis_epoch_dir = $(green_nn_dir)/s28_cod_n7p5_nt_n21p17/epoch10

green_full_codon_scores_results_files = \
	$(green_results_full_epoch_dir)/codon_scores.tsv \
	$(green_results_full_epoch_dir)/codon_scores_colormap.pdf

green_full_codon_scores_results_files_pattern = \
	$(green_results_full_epoch_dir)/codon_scores%tsv \
	$(green_results_full_epoch_dir)/codon_scores_colormap%pdf

green_28_codon_scores_results_files = \
        $(green_results_28_epoch_dir)/codon_scores.tsv \
        $(green_results_28_epoch_dir)/codon_scores_colormap.pdf

green_28_codon_scores_results_files_pattern = \
        $(green_results_28_epoch_dir)/codon_scores%tsv \
        $(green_results_28_epoch_dir)/codon_scores_colormap%pdf

green_leaveout_mses_file = \
	$(green_results_leaveout_dir)/leaveout_mses.txt

### figure variable defintions

fig_files = \
	$(fig_dir)/figure_1B_scaledcts.pdf \
	$(fig_dir)/figure_1C_mse.pdf \
	$(fig_dir)/figure_1D_scatter.pdf \
	$(fig_dir)/figure_1E_indiv_gene.pdf \
	$(fig_dir)/figure_1F_binned_error.pdf \
	$(fig_dir)/figure_2A_codonmse.pdf \
	$(fig_dir)/figure_2B_heatmap.pdf \
	$(fig_dir)/figure_2C_tai.pdf \
	$(fig_dir)/figure_2D_wobble.pdf \
	$(fig_dir)/figure_2E_lareaucodonmse.pdf \
	$(fig_dir)/figure_2F_5prime.pdf \
	$(fig_dir)/figure_2G_asite.pdf \
	$(fig_dir)/figure_2H_cl2.pdf \
	$(fig_dir)/figure_2I_cl1.pdf \
	$(fig_dir)/figure_3B_citrine_dist.pdf \
	$(fig_dir)/figure_3C_facs.pdf \
	$(fig_dir)/figure_3D_te.pdf \
	$(fig_dir)/supp_table_codon_scores.csv \
	$(fig_dir)/supp_figure_mrna_qpcr.pdf \
	$(fig_dir)/supp_figure_facs.pdf \
	$(fig_dir)/supp_figure_greenmse.pdf \
	$(fig_dir)/supp_figure_iwasakimse.pdf \
	$(fig_dir)/supp_figure_cl1v2.pdf

### paper data files

paper_data_dir = $(results_dir)/paper_data
paper_data_fname = $(paper_data_dir)/paper_data.txt

.PHONY: install expts weinberg_expt iwasaki_expt iwasaki_expt weinberg_clean figures

################################
### all: Overall make command
################################

all: install expts paper_data

##########################################################
### Installation: Build basic top level file structures
##########################################################

# NOTE: This is where I test definitions in the makefile
test: 
	echo $(weinberg_plot_files_pattern)

# NOTE: Installation graveyard. Has been (mostly?) replaced by cloning git repo
#install: \
	$(lib_files) \
	$(genome_files) \
	$(repro_files) \
	$(struc_files) \
	$(wetlab_files) \
	$(results_dir) \
	$(fig_dir) \
	$(yeast_30_windows_seqs_file) \
	$(yeast_30_windows_str_scores_file)

# Make main directory
$(top_dir):
	mkdir $(top_dir)

#$(lib_dir): | $(top_dir)
#	if ! [ -a $(lib_dir) ] ; then mkdir $(lib_dir) ; fi;

#$(genome_dir): | $(top_dir)
#	mkdir $(genome_dir)

#$(struc_dir): | $(top_dir)
#	mkdir $(struc_dir)

#$(wetlab_dir): | $(top_dir)
#	mkdir $(wetlab_dir)

#$(repro_dir): | $(top_dir)
#	mkdir $(repro_dir)

#$(lib_files): | $(lib_dir)
#	cp /mnt/lareaulab/rtunney/Regression/rp_predict/$(subst $(lib_dir)/,,$@) $@

#$(genome_files): | $(genome_dir)
#	cp /mnt/lareaulab/rtunney/Regression/genome_data/$(subst $(genome_dir)/,,$@) $@

#$(struc_files): | $(struc_dir)
#	cp /mnt/lareaulab/rtunney/Regression/structure_data/$(subst $(struc_dir)/,,$@) $@

#$(wetlab_files): | $(wetlab_dir)
#	cp /mnt/lareaulab/rtunney/Regression/wetlab_data/$(subst $(wetlab_dir)/,,$@) $@

#$(repro_files): | $(repro_dir)
#	cp /mnt/lareaulab/rtunney/Regression/reproduce_scripts/$(subst $(repro_dir)/,,$@) $@

#$(expts_dir): | $(top_dir)
#	mkdir $(expts_dir)

$(results_dir): | $(top_dir)
	mkdir $(results_dir)

$(fig_dir): | $(results_dir)
	mkdir $(fig_dir)

$(yeast_30_windows_seqs_file): \
		$(yeast_30_windows_seqs_file_zip)
	gunzip -k $(yeast_30_windows_seqs_file_zip)

$(yeast_30_windows_str_scores_file): \
		$(yeast_30_windows_str_scores_file_zip)
	gunzip -k $(yeast_30_windows_str_scores_file_zip)

#########################################################
### Experiments: Reproduce analyses for RP experiments
#########################################################

expts: weinberg_expt iwasaki_expt lareau_expt green_expt

##########################
### Weinberg Experiment
##########################

weinberg_expt: \
		$(weinberg_sam_file) $(weinberg_plot_files) \
		$(weinberg_27_31_proc_sam_files) \
		$(weinberg_28_proc_sam_files) \
		$(weinberg_feat_nb_series_files) \
		$(weinberg_leaveout_series_files) \
		$(weinberg_28mer_files) \
		$(weinberg_fp_struc_files) \
		$(weinberg_max_struc_files) \
		$(weinberg_linreg_series_files) \
		$(weinberg_full_codon_scores_results_files) \
		$(weinberg_28_codon_scores_results_files) \
		$(weinberg_feat_nb_mses_file) \
		$(weinberg_leaveout_mses_file) \
		$(weinberg_feat_nb_linreg_mses_file) \
		$(weinberg_struc_mses_file) \
		$(weinberg_struc_plot_files) \
		$(weinberg_full_plot_files) \
		$(weinberg_28_plot_files) \
		$(weinberg_struc_gene_plot_dir) \
		$(weinberg_final_model_y_te) \
		$(weinberg_final_model_y_te_hat)\
		$(weinberg_results_final_model_y_te)\
		$(weinberg_results_final_model_y_te_hat)\
		$(weinberg_opt_series_file) \
		| $(weinberg_expt_dir) $(weinberg_subdirs)

$(weinberg_expt_dir): | $(expts_dir)
	mkdir $(weinberg_expt_dir)

$(weinberg_proc_dir): | $(weinberg_expt_dir)
	mkdir $(weinberg_proc_dir)

$(weinberg_plot_dir): | $(weinberg_expt_dir)
	mkdir $(weinberg_plot_dir)

$(weinberg_nn_dir): | $(weinberg_expt_dir)
	mkdir $(weinberg_nn_dir)

$(weinberg_lr_dir): | $(weinberg_expt_dir)
	mkdir $(weinberg_lr_dir)

$(weinberg_raw_fastq_file) : \
		| $(weinberg_proc_dir)
	fastq-dump SRR1049521 -O $(weinberg_proc_dir)

#NOTE: This is missing a bunch of prereqs
$(weinberg_raw_sam_file): \
		| $(weinberg_raw_fastq_file) \
		$(weinberg_trimmer_script) \
		$(weinberg_mapping_script) \
		$(weinberg_proc_dir)
	bash $(weinberg_mapping_script) $(genome_dir) $(weinberg_proc_dir)

$(weinberg_sam_file): \
		| $(weinberg_raw_sam_file) \
		$(weinberg_proc_dir) $(weinberg_expt_dir)
	python $(repro_dir)/process_data.py edit_sam_file weinberg \
		$(weinberg_expt_dir) $(weinberg_raw_sam_file)

$(weinberg_plot_files_pattern): \
		| $(weinberg_sam_file) $(yeast_gene_len_file) \
		${weinberg_plot_dir} $(weinberg_expt_dir)
	python $(repro_dir)/process_data.py size_and_frame_analysis weinberg \
		$(weinberg_expt_dir) $(weinberg_sam_file) \
		$(yeast_gene_len_file)

$(weinberg_27_31_proc_sam_pattern): \
		| $(weinberg_sam_file) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(weinberg_proc_dir) $(weinberg_expt_dir)
	python $(repro_dir)/process_data.py process_sam_file weinberg \
		$(weinberg_expt_dir) $(weinberg_sam_file) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(yeast_paralogs_file)

$(weinberg_28_proc_sam_pattern): \
		| $(weinberg_sam_file) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(weinberg_proc_dir) $(weinberg_expt_dir)
	python $(repro_dir)/process_data.py process_sam_file_28mers weinberg \
		$(weinberg_expt_dir) $(weinberg_sam_file) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(yeast_paralogs_file)

#Leaveout series pattern depends on "full_" at beginning of model name to 
#       determine when this rule should be applied (cf. leaveout series)

$(weinberg_feat_nb_series_pattern): \
		| $(weinberg_sam_file) \
		$(weinberg_27_31_proc_sam_files) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(weinberg_27_31_tr_bounds) $(weinberg_27_31_te_bounds) \
		$(weinberg_27_31_outputs) $(weinberg_nn_dir) \
		$(weinberg_expt_dir)
	$(eval model_name=$(firstword $(subst _rep, ,$*)))
	$(eval model_rep=$(lastword $(subst _rep, ,$*)))
	echo $*
	echo ...Making model $(model_name) rep $(model_rep)
	echo
	python $(repro_dir)/feat_neighborhood_nn_series.py \
		$(model_name) $(model_rep) \
		$(weinberg_expt_dir) $(weinberg_sam_file) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(weinberg_27_31_tr_bounds) $(weinberg_27_31_te_bounds) \
		$(weinberg_27_31_outputs) 30

#Leaveout series pattern depends on "nocod" at beginning of model name to 
#	determine when this rule should be applied (cf. feat nb series)
#Note: We have to stick "nocod" back in the beginning of the model name
#	when parsing the pattern, to feed the model name to python script

$(weinberg_leaveout_series_pattern): \
		| $(weinberg_sam_file) \
		$(weinberg_27_31_proc_sam_files) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(weinberg_27_31_tr_bounds) $(weinberg_27_31_te_bounds) \
		$(weinberg_27_31_outputs) $(weinberg_nn_dir) \
		$(weinberg_expt_dir)
	$(eval model_name=nocod$(firstword $(subst _rep, ,$*)))
	$(eval model_rep=$(lastword $(subst _rep, ,$*)))
	echo
	echo ...Making model $(model_name) rep $(model_rep)
	echo
	python $(repro_dir)/leaveout_series.py \
		$(model_name) $(model_rep) \
		$(weinberg_expt_dir) $(weinberg_sam_file) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(weinberg_27_31_tr_bounds) $(weinberg_27_31_te_bounds) \
		$(weinberg_27_31_outputs) 30

$(weinberg_28mer_pattern): \
		| $(weinberg_sam_file) \
		$(weinberg_28_proc_sam_files) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(weinberg_28_tr_bounds) $(weinberg_28_te_bounds) \
		$(weinberg_28_outputs) $(weinberg_nn_dir) \
		$(weinberg_expt_dir)
	$(eval model_name=s28_$*)
	echo
	echo ...Making model $(model_name)
	echo
	python $(repro_dir)/28mer_models.py \
		$(model_name) \
		$(weinberg_expt_dir) $(weinberg_sam_file) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(weinberg_28_tr_bounds) $(weinberg_28_te_bounds) \
		$(weinberg_28_outputs) 25

$(weinberg_fp_struc_pattern): \
		| $(weinberg_sam_file) \
		$(weinberg_27_31_proc_sam_files) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(yeast_30_windows_str_scores_file) \
		$(weinberg_27_31_tr_bounds) $(weinberg_27_31_te_bounds) \
		$(weinberg_27_31_outputs) $(weinberg_nn_dir) \
		$(weinberg_expt_dir)
	$(eval model_name=$(firstword $(subst _rep, ,$*)))
	$(eval model_rep=$(lastword $(subst _rep, ,$*)))
	echo $*
	echo ...Making model $(model_name) rep $(model_rep)
	echo
	python $(repro_dir)/struc_series.py \
		$(model_name) $(model_rep) \
		$(weinberg_expt_dir) $(weinberg_sam_file) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(yeast_30_windows_str_scores_file) \
		$(weinberg_27_31_tr_bounds) $(weinberg_27_31_te_bounds) \
		$(weinberg_27_31_outputs) 30

$(weinberg_max_struc_pattern): \
		| $(weinberg_sam_file) \
		$(weinberg_27_31_proc_sam_files) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(yeast_30_windows_str_scores_file) \
		$(weinberg_27_31_tr_bounds) $(weinberg_27_31_te_bounds) \
		$(weinberg_27_31_outputs) $(weinberg_nn_dir) \
		$(weinberg_expt_dir)
	$(eval model_name=$(firstword $(subst _rep, ,$*)))
	$(eval model_rep=$(lastword $(subst _rep, ,$*)))
	echo $*
	echo ...Making model $(model_name) rep $(model_rep)
	echo
	python $(repro_dir)/struc_series.py \
		$(model_name) $(model_rep) \
		$(weinberg_expt_dir) $(weinberg_sam_file) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(yeast_30_windows_str_scores_file) \
		$(weinberg_27_31_tr_bounds) $(weinberg_27_31_te_bounds) \
		$(weinberg_27_31_outputs) 35

$(weinberg_linreg_series_pattern): \
		| $(weinberg_sam_file) \
		$(weinberg_27_31_proc_sam_files) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(weinberg_27_31_tr_bounds) $(weinberg_27_31_te_bounds) \
		$(weinberg_27_31_outputs) $(weinberg_lr_dir) \
		$(weinberg_expt_dir)
	$(eval model_name=$(lastword $(subst lr_, ,$*)))
	echo $*
	echo ...Making linreg model $(model_name)
	echo
	python $(repro_dir)/feat_neighborhood_linreg_series.py \
		$(model_name) $(weinberg_expt_dir) $(yeast_gene_len_file) \
		$(yeast_gene_seq_file) $(weinberg_27_31_tr_bounds) \
		$(weinberg_27_31_te_bounds) $(weinberg_27_31_outputs)

# Make weinberg results directory and files
$(weinberg_results_dir): | $(results_dir)
	mkdir $(weinberg_results_dir)

$(weinberg_results_feat_neighborhood_dir): | $(weinberg_results_dir)
	mkdir $(weinberg_results_feat_neighborhood_dir)

$(weinberg_results_leaveout_dir): | $(weinberg_results_dir)
	mkdir $(weinberg_results_leaveout_dir)

$(weinberg_results_struc_series_dir): | $(weinberg_results_dir)
	mkdir $(weinberg_results_struc_series_dir)

$(weinberg_results_struc_dir): | $(weinberg_results_dir)
	mkdir $(weinberg_results_struc_dir)

$(weinberg_results_struc_epoch_dir): | $(weinberg_results_struc_dir)
	mkdir $(weinberg_results_struc_epoch_dir)

$(weinberg_results_full_dir): | $(weinberg_results_dir)
	mkdir $(weinberg_results_full_dir)

$(weinberg_results_full_epoch_dir): | $(weinberg_results_full_dir)
	mkdir $(weinberg_results_full_epoch_dir)

$(weinberg_results_28_dir): | $(weinberg_results_dir)
	mkdir $(weinberg_results_28_dir)

$(weinberg_results_28_epoch_dir): | $(weinberg_results_28_dir)
	mkdir $(weinberg_results_28_epoch_dir)

$(weinberg_results_opt_model_dir): | $(weinberg_results_dir)
	mkdir $(weinberg_results_opt_model_dir)

$(weinberg_results_opt_model_epoch_dir): | $(weinberg_results_opt_model_dir)
	mkdir $(weinberg_results_opt_model_epoch_dir)

#NOTE: We put an example .pkl file as a dependency so we can be sure this
# epoch is made. We can't put the epoch file as a dependency bc there's no
# rule to create it, and then make can't place the targets properly in its DAG.

$(weinberg_full_codon_scores_results_files_pattern): \
		| $(weinberg_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/epoch30/te_cost_by_epoch.pkl \
		$(weinberg_results_full_epoch_dir) 
	python $(repro_dir)/codon_scores.py \
		$(weinberg_nn_dir)/full_cod_n7p5_nt_n21p17_rep0 30
	cp $(weinberg_full_analysis_epoch_dir)/codon_scores.tsv \
		$(weinberg_results_full_epoch_dir)/codon_scores.tsv
	cp $(weinberg_full_analysis_epoch_dir)/codon_scores_colormap.pdf \
		$(weinberg_results_full_epoch_dir)/codon_scores_colormap.pdf

#NOTE: We put an example .pkl file as a dependency so we can be sure this
# epoch is made. We can't put the epoch file as a dependency bc there's no 
# rule to create it, and then make can't place the targets properly in its DAG.

$(weinberg_28_codon_scores_results_files_pattern): \
		| $(weinberg_nn_dir)/s28_cod_n7p5_nt_n21p17/epoch25/te_cost_by_epoch.pkl \
		$(weinberg_results_28_epoch_dir)
	python $(repro_dir)/codon_scores.py \
		$(weinberg_nn_dir)/s28_cod_n7p5_nt_n21p17 25
	cp $(weinberg_28_analysis_epoch_dir)/codon_scores.tsv \
		$(weinberg_results_28_epoch_dir)/codon_scores.tsv
	cp $(weinberg_28_analysis_epoch_dir)/codon_scores_colormap.pdf \
		$(weinberg_results_28_epoch_dir)/codon_scores_colormap.pdf

$(weinberg_feat_nb_mses_file): \
		| $(weinberg_feat_nb_series_files) \
		$(weinberg_results_feat_neighborhood_dir)
	python $(repro_dir)/aggregate_mses.py aggregate_w_diffs \
		$(weinberg_nn_dir) 30 10 \
		$(weinberg_feat_nb_mses_file) full_cod_n7p5_nt_n21p17 \
		$(addprefix full_,$(feat_neighborhoods))

$(weinberg_leaveout_mses_file): \
		| $(weinberg_leaveout_series_files) \
		$(weinberg_results_leaveout_dir)
	python $(repro_dir)/aggregate_mses.py aggregate_w_diffs \
		$(weinberg_nn_dir) 30 10 \
		$(weinberg_leaveout_mses_file) full_cod_n7p5_nt_n21p17 \
		$(addsuffix _cod_n7p5_nt_n21p17, \
			$(addprefix nocod,$(leaveout_cods)))

$(weinberg_struc_mses_file): \
		| $(weinberg_fp_struc_files) \
		$(weinberg_max_struc_files) \
		$(weinberg_results_struc_series_dir)
	python $(repro_dir)/aggregate_mses.py aggregate\
		$(weinberg_nn_dir) 30 10 \
		$(weinberg_struc_mses_file) str_n17n15_cod_n7p5_nt_n21p17 \
		max_str_p13p42_cod_n7p5_nt_n21p17 

$(weinberg_feat_nb_linreg_mses_file): \
		| $(weinberg_linreg_series_files) \
		$(weinberg_results_leaveout_dir)
	python $(repro_dir)/aggregate_linreg_mses.py $(weinberg_lr_dir) \
		$(weinberg_feat_nb_linreg_mses_file) \
		$(addprefix lr_,$(feat_neighborhoods))

$(weinberg_struc_plot_pattern): \
		| $(repro_dir)/plot_nn.py \
		$(weinberg_results_struc_epoch_dir) \
		$(weinberg_nn_dir)/str_n17n15_cod_n7p5_nt_n21p17_rep0/epoch30/te_cost_by_epoch.pkl 
	echo
	echo "Plotting nn str_n17n15_cod_n7p5_nt_n21p17_rep0"
	echo
	python $(repro_dir)/plot_nn.py \
		$(weinberg_expt_dir) str_n17n15_cod_n7p5_nt_n21p17_rep0 30
	cp -r $(weinberg_nn_dir)/str_n17n15_cod_n7p5_nt_n21p17_rep0/plots \
		$(weinberg_results_struc_dir)

$(weinberg_full_plot_pattern): \
		| $(repro_dir)/plot_nn.py \
		$(weinberg_results_full_epoch_dir) \
		$(weinberg_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/epoch30/te_cost_by_epoch.pkl 
	echo
	echo "Plotting nn full_cod_n7p5_nt_n21p17_rep0"
	echo
	python $(repro_dir)/plot_nn.py \
		$(weinberg_expt_dir) full_cod_n7p5_nt_n21p17_rep0 30
	cp -r $(weinberg_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/plots \
		$(weinberg_results_full_dir)

$(weinberg_28_plot_pattern): \
		| $(repro_dir)/plot_nn.py \
		$(weinberg_results_28_epoch_dir) \
		$(weinberg_nn_dir)/s28_cod_n7p5_nt_n21p17/epoch25/te_cost_by_epoch.pkl 
	echo
	echo "Plotting nn s28_n17n15_cod_n7p5_nt_n21p17"
	echo
	python $(repro_dir)/plot_nn.py \
		$(weinberg_expt_dir) s28_cod_n7p5_nt_n21p17 25
	cp -r $(weinberg_nn_dir)/s28_cod_n7p5_nt_n21p17/plots \
		$(weinberg_results_28_dir)

$(weinberg_struc_gene_plot_dir): \
		| $(weinberg_nn_dir)/str_n17n15_cod_n7p5_nt_n21p17_rep0/epoch30/te_cost_by_epoch.pkl \
		$(weinberg_results_struc_epoch_dir)
	echo
	echo "Plotting str_n17n15_cod_n7p5_nt_n21p17_rep0 genes"
	echo
	python $(repro_dir)/plot_genes.py \
		$(weinberg_nn_dir)/str_n17n15_cod_n7p5_nt_n21p17_rep0 30 \
		$(weinberg_expt_dir)/process/te_set_bounds.size.27.31.trunc.20.20.min_cts.200.min_cod.100.top.500.data_table.txt \
		$(yeast_gene_symbol_file) \
		$(weinberg_nn_dir)/str_n17n15_cod_n7p5_nt_n21p17_rep0/gene_plots

$(weinberg_final_model_y_te): \
		| $(weinberg_nn_dir)/str_n17n15_cod_n7p5_nt_n21p17_rep0/init_data/y_te.pkl \
		$(repro_dir)/pkl2txt.py
	python $(repro_dir)/pkl2txt.py \
		$(weinberg_nn_dir)/str_n17n15_cod_n7p5_nt_n21p17_rep0/init_data/y_te.pkl \
		$(weinberg_final_model_y_te)

$(weinberg_final_model_y_te_hat): \
		| $(weinberg_nn_dir)/str_n17n15_cod_n7p5_nt_n21p17_rep0/epoch30/y_te_hat.pkl \
		$(repro_dir)/pkl2txt.py
	python $(repro_dir)/pkl2txt.py \
		$(weinberg_nn_dir)/str_n17n15_cod_n7p5_nt_n21p17_rep0/epoch30/y_te_hat.pkl \
		$(weinberg_final_model_y_te_hat)

$(weinberg_results_final_model_y_te): \
		| $(weinberg_final_model_y_te) \
		$(weinberg_results_struc_epoch_dir) 
	cp $(weinberg_final_model_y_te) \
		$(weinberg_results_struc_epoch_dir)/y_te.txt

$(weinberg_results_final_model_y_te_hat): \
		| $(weinberg_final_model_y_te_hat) \
		$(weinberg_results_struc_epoch_dir) 
	cp $(weinberg_final_model_y_te_hat) \
		$(weinberg_results_struc_epoch_dir)/y_te_hat.txt

$(weinberg_opt_series_file): \
		| $(weinberg_nn_dir)/full_cod_n3p2_nt_n9p8_rep0/init_data/init_data.pkl \
		$(repro_dir)/optimize_cds.py \
		$(weinberg_results_opt_model_epoch_dir) 
	python $(repro_dir)/optimize_cds.py \
		$(weinberg_nn_dir)/full_cod_n3p2_nt_n9p8_rep0 30 True \
		$(citrine_aa_seq) $(weinberg_opt_series_file) 

weinberg_clean: 
	rm $(weinberg_mapped_sam_file)

##########################
### Lareau Experiment
##########################

lareau_expt: $(lareau_sam_file) ${lareau_plot_files} \
		$(lareau_27_29_proc_sam_files) \
		$(lareau_28_proc_sam_files) \
		$(lareau_full_model_files) \
		$(lareau_leaveout_series_files) \
		$(lareau_28mer_files) \
		$(lareau_full_codon_scores_results_files) \
		$(lareau_28_codon_scores_results_files) \
		$(lareau_leaveout_mses_file) \
		$(lareau_full_plot_files) \
		$(lareau_28_plot_files) \
		| $(lareau_expt_dir) $(lareau_subdirs)

$(lareau_expt_dir): | $(expts_dir)
	mkdir $(lareau_expt_dir)

$(lareau_proc_dir): | $(lareau_expt_dir)
	mkdir $(lareau_proc_dir)

$(lareau_plot_dir): | $(lareau_expt_dir)
	mkdir $(lareau_plot_dir)

$(lareau_nn_dir): | $(lareau_expt_dir)
	mkdir $(lareau_nn_dir)

$(lareau_lr_dir): | $(lareau_expt_dir)
	mkdir $(lareau_lr_dir)

$(lareau_raw_fastq_file_1): | $(lareau_proc_dir)
	wget storage.ingolia-lab.org/lareaulab/iXnos/LLMG004_S31_L007_R1_001.fastq.gz -O $(lareau_proc_dir)/LLMG004_S31_L007_R1_001.fastq.gz 
	gunzip -k $(lareau_proc_dir)/LLMG004_S31_L007_R1_001.fastq.gz 

$(lareau_raw_fastq_file_2): | $(lareau_proc_dir)
	wget storage.ingolia-lab.org/lareaulab/iXnos/LLMG005_S1_L001_R1_001.fastq.gz -O $(lareau_proc_dir)/LLMG005_S1_L001_R1_001.fastq.gz
	gunzip -k $(lareau_proc_dir)/LLMG005_S1_L001_R1_001.fastq.gz 

#NOTE: This is missing a bunch of prereqs
$(lareau_raw_sam_file): \
		| $(lareau_raw_fastq_files) \
		$(lareau_trimmer_script) \
		$(lareau_mapping_script) \
		$(lareau_proc_dir)
	bash $(lareau_mapping_script) $(genome_dir) $(lareau_proc_dir)

$(lareau_sam_file): | $(lareau_raw_sam_file) \
		$(lareau_proc_dir) $(lareau_expt_dir)
	python $(repro_dir)/process_data.py edit_sam_file lareau \
		$(lareau_expt_dir) $(lareau_raw_sam_file)

${lareau_plot_files_pattern}: \
		| $(lareau_sam_file) $(yeast_gene_len_file) \
		${lareau_plot_dir} $(lareau_expt_dir)
	python $(repro_dir)/process_data.py size_and_frame_analysis lareau \
		$(lareau_expt_dir) $(lareau_sam_file) \
		$(yeast_gene_len_file)

$(lareau_27_29_proc_sam_pattern): \
		| $(lareau_sam_file) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(lareau_proc_dir) $(lareau_expt_dir)
	python $(repro_dir)/process_data.py process_sam_file lareau \
		$(lareau_expt_dir) $(lareau_sam_file) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(yeast_paralogs_file)

$(lareau_28_proc_sam_pattern): \
                | $(lareau_sam_file) \
                $(yeast_gene_len_file) $(yeast_gene_seq_file) \
                $(lareau_proc_dir) $(lareau_expt_dir)
	python $(repro_dir)/process_data.py process_sam_file_28mers lareau \
		$(lareau_expt_dir) $(lareau_sam_file) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(yeast_paralogs_file)

$(lareau_full_model_pattern): \
		| $(lareau_sam_file) \
		$(lareau_27_29_proc_sam_files) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(lareau_27_29_tr_bounds) $(lareau_27_29_te_bounds) \
		$(lareau_27_29_outputs) $(lareau_nn_dir) \
		$(lareau_expt_dir)
	$(eval model_name=$(firstword $(subst _rep, ,$*)))
	$(eval model_rep=$(lastword $(subst _rep, ,$*)))
	echo
	echo ...Making model $(model_name) rep $(model_rep)
	echo
	python $(repro_dir)/feat_neighborhood_nn_series.py \
		$(model_name) $(model_rep) \
		$(lareau_expt_dir) $(lareau_sam_file) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(lareau_27_29_tr_bounds) $(lareau_27_29_te_bounds) \
		$(lareau_27_29_outputs) 15

$(lareau_leaveout_series_pattern): \
		| $(lareau_sam_file) \
		$(lareau_27_29_proc_sam_files) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(lareau_27_29_tr_bounds) $(lareau_27_29_te_bounds) \
		$(lareau_27_29_outputs) $(lareau_nn_dir) \
		$(lareau_expt_dir)
	$(eval model_name=nocod$(firstword $(subst _rep, ,$*)))
	$(eval model_rep=$(lastword $(subst _rep, ,$*)))
	echo
	echo ...Making model $(model_name) rep $(model_rep)
	echo
	python $(repro_dir)/leaveout_series.py \
		$(model_name) $(model_rep) \
		$(lareau_expt_dir) $(lareau_sam_file) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(lareau_27_29_tr_bounds) $(lareau_27_29_te_bounds) \
		$(lareau_27_29_outputs) 15

$(lareau_28mer_pattern): \
		| $(lareau_sam_file) \
		$(lareau_28_proc_sam_files) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(lareau_28_tr_bounds) $(lareau_28_te_bounds) \
		$(lareau_28_outputs) $(lareau_nn_dir) \
		$(lareau_expt_dir)
	$(eval model_name=s28_$*)
	echo
	echo ...Making model $(model_name)
	echo
	python $(repro_dir)/28mer_models.py \
		$(model_name) \
		$(lareau_expt_dir) $(lareau_sam_file) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(lareau_28_tr_bounds) $(lareau_28_te_bounds) \
		$(lareau_28_outputs) 10

# Make lareau results directory and files
$(lareau_results_dir): | $(results_dir)
	mkdir $(lareau_results_dir)

$(lareau_results_feat_neighborhood_dir): | $(lareau_results_dir)
	mkdir $(lareau_results_feat_neighborhood_dir)

$(lareau_results_leaveout_dir): | $(lareau_results_dir)
	mkdir $(lareau_results_leaveout_dir)

$(lareau_results_full_dir): | $(lareau_results_dir)
	mkdir $(lareau_results_full_dir)

$(lareau_results_full_epoch_dir): | $(lareau_results_full_dir)
	mkdir $(lareau_results_full_epoch_dir)

$(lareau_results_28_dir): | $(lareau_results_dir)
	mkdir $(lareau_results_28_dir)

$(lareau_results_28_epoch_dir): | $(lareau_results_28_dir)
	mkdir $(lareau_results_28_epoch_dir)

#NOTE: We put an example .pkl file as a dependency so we can be sure this
# epoch is made. We can't put the epoch file as a dependency bc there's no 
# rule to create it, and then make can't place the targets properly in its DAG. 

$(lareau_full_codon_scores_results_files_pattern): \
		| $(lareau_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/epoch15/te_cost_by_epoch.pkl \
		$(lareau_results_full_epoch_dir)
	python $(repro_dir)/codon_scores.py \
		$(lareau_nn_dir)/full_cod_n7p5_nt_n21p17_rep0 15
	cp $(lareau_full_analysis_epoch_dir)/codon_scores.tsv \
		$(lareau_results_full_epoch_dir)/codon_scores.tsv
	cp $(lareau_full_analysis_epoch_dir)/codon_scores_colormap.pdf \
		$(lareau_results_full_epoch_dir)/codon_scores_colormap.pdf

#NOTE: We put an example .pkl file as a dependency so we can be sure this
# epoch is made. We can't put the epoch file as a dependency bc there's no
# rule to create it, and then make can't place the targets properly in its DAG. 

$(lareau_28_codon_scores_results_files_pattern): \
		| $(lareau_nn_dir)/s28_cod_n7p5_nt_n21p17/epoch10/te_cost_by_epoch.pkl \
		$(lareau_results_28_epoch_dir)
	python $(repro_dir)/codon_scores.py \
		$(lareau_nn_dir)/s28_cod_n7p5_nt_n21p17 10
	cp $(lareau_28_analysis_epoch_dir)/codon_scores.tsv \
		$(lareau_results_28_epoch_dir)/codon_scores.tsv
	cp $(lareau_28_analysis_epoch_dir)/codon_scores_colormap.pdf \
		$(lareau_results_28_epoch_dir)/codon_scores_colormap.pdf

$(lareau_leaveout_mses_file) : \
		| $(lareau_leaveout_series_files) \
		$(lareau_results_leaveout_dir)
	python $(repro_dir)/aggregate_mses.py aggregate_w_diffs \
		$(lareau_nn_dir) 15 10 \
		$(lareau_leaveout_mses_file) full_cod_n7p5_nt_n21p17 \
		$(addsuffix _cod_n7p5_nt_n21p17, \
			$(addprefix nocod,$(leaveout_cods)))

$(lareau_full_plot_pattern): \
		| $(repro_dir)/plot_nn.py \
		$(lareau_results_full_epoch_dir) \
		$(lareau_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/epoch15/te_cost_by_epoch.pkl
	echo
	echo "Plotting nn full_cod_n7p5_n21p17_rep0"
	echo
	python $(repro_dir)/plot_nn.py \
		$(lareau_expt_dir) full_cod_n7p5_nt_n21p17_rep0 15
	cp -r $(lareau_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/plots \
		$(lareau_results_full_dir)

$(lareau_28_plot_pattern): \
		| $(repro_dir)/plot_nn.py \
		$(lareau_results_28_epoch_dir) \
		$(lareau_nn_dir)/s28_cod_n7p5_nt_n21p17/epoch10/te_cost_by_epoch.pkl
	echo
	echo "Plotting nn s28_cod_n7p5_n21p17_rep0"
	echo
	python $(repro_dir)/plot_nn.py \
		$(lareau_expt_dir) s28_cod_n7p5_nt_n21p17 10
	cp -r $(lareau_nn_dir)/s28_cod_n7p5_nt_n21p17/plots \
		$(lareau_results_28_dir)

lareau_clean: 
	rm $(lareau_mapped_sam_file)

#########################
### Iwasaki Experiment
#########################

iwasaki_expt: $(iwasaki_sam_file) ${iwasaki_plot_files} \
		$(iwasaki_27_30_proc_sam_files) \
		$(iwasaki_28_proc_sam_files) \
		$(iwasaki_full_model_files) \
		$(iwasaki_leaveout_series_files) \
		$(iwasaki_28mer_files) \
		$(iwasaki_full_codon_scores_results_files) \
		$(iwasaki_28_codon_scores_results_files) \
		$(iwasaki_leaveout_mses_file) \
		$(iwasaki_full_plot_files) \
		$(iwasaki_28_plot_files) \
		| $(iwasaki_expt_dir) $(iwasaki_subdirs)

$(iwasaki_expt_dir): | $(expts_dir)
	mkdir $(iwasaki_expt_dir)

$(iwasaki_proc_dir): | $(iwasaki_expt_dir)
	mkdir $(iwasaki_proc_dir)

$(iwasaki_plot_dir): | $(iwasaki_expt_dir)
	mkdir $(iwasaki_plot_dir)

$(iwasaki_nn_dir): | $(iwasaki_expt_dir)
	mkdir $(iwasaki_nn_dir)

$(iwasaki_lr_dir): | $(iwasaki_expt_dir)
	mkdir $(iwasaki_lr_dir)

$(iwasaki_raw_fastq_files): \
		| $(iwasaki_proc_dir)
	fastq-dump SRR2075925 -O $(iwasaki_proc_dir)
	fastq-dump SRR2075926 -O $(iwasaki_proc_dir)

#NOTE: This is missing a bunch of prereqs
$(iwasaki_raw_sam_file): \
		| $(iwasaki_raw_fastq_files) \
		$(iwasaki_trimmer_script) \
		$(iwasaki_mapping_script) \
		$(iwasaki_proc_dir)
	bash $(iwasaki_mapping_script) $(genome_dir) $(iwasaki_proc_dir)

$(iwasaki_sam_file): \
		| $(iwasaki_raw_sam_file) $(repro_dir)/process_data.py \
		$(iwasaki_proc_dir) $(iwasaki_expt_dir)
	python $(repro_dir)/process_data.py edit_sam_file iwasaki \
		$(iwasaki_expt_dir) $(iwasaki_raw_sam_file)

${iwasaki_plot_files_pattern}: \
		| $(iwasaki_sam_file) $(human_gene_len_file) \
		${iwasaki_plot_dir} $(iwasaki_expt_dir)
	python $(repro_dir)/process_data.py size_and_frame_analysis iwasaki \
		$(iwasaki_expt_dir) $(iwasaki_sam_file) \
		$(human_gene_len_file)

$(iwasaki_27_30_proc_sam_pattern): \
		| $(iwasaki_sam_file) \
		$(human_gene_len_file) $(human_gene_seq_file) \
		$(iwasaki_proc_dir) $(iwasaki_expt_dir)
	python $(repro_dir)/process_data.py process_sam_file iwasaki \
		$(iwasaki_expt_dir) $(iwasaki_sam_file) \
		$(human_gene_len_file) $(human_gene_seq_file)

$(iwasaki_28_proc_sam_pattern): \
                | $(iwasaki_sam_file) \
                $(human_gene_len_file) $(human_gene_seq_file) \
                $(iwasaki_proc_dir) $(iwasaki_expt_dir)
	python $(repro_dir)/process_data.py process_sam_file_28mers iwasaki \
		$(iwasaki_expt_dir) $(iwasaki_sam_file) \
		$(human_gene_len_file) $(human_gene_seq_file)

$(iwasaki_full_model_pattern): \
		| $(iwasaki_sam_file) \
		$(iwasaki_27_30_proc_sam_files) \
		$(human_gene_len_file) $(human_gene_seq_file) \
		$(iwasaki_27_30_tr_bounds) $(iwasaki_27_30_te_bounds) \
		$(iwasaki_27_30_outputs) $(iwasaki_nn_dir) \
		$(iwasaki_expt_dir)
	$(eval model_name=$(firstword $(subst _rep, ,$*)))
	$(eval model_rep=$(lastword $(subst _rep, ,$*)))
	echo
	echo ...Making model $(model_name) rep $(model_rep)
	echo
	python $(repro_dir)/feat_neighborhood_nn_series.py \
		$(model_name) $(model_rep) \
		$(iwasaki_expt_dir) $(iwasaki_sam_file) \
		$(human_gene_len_file) $(human_gene_seq_file) \
		$(iwasaki_27_30_tr_bounds) $(iwasaki_27_30_te_bounds) \
		$(iwasaki_27_30_outputs) 10

$(iwasaki_leaveout_series_pattern): \
		| $(iwasaki_sam_file) \
		$(iwasaki_27_30_proc_sam_files) \
		$(human_gene_len_file) $(human_gene_seq_file) \
		$(iwasaki_27_30_tr_bounds) $(iwasaki_27_30_te_bounds) \
		$(iwasaki_27_30_outputs) $(iwasaki_nn_dir) \
		$(iwasaki_expt_dir)
	$(eval model_name=nocod$(firstword $(subst _rep, ,$*)))
	$(eval model_rep=$(lastword $(subst _rep, ,$*)))
	echo
	echo ...Making model $(model_name) rep $(model_rep)
	echo
	python $(repro_dir)/leaveout_series.py \
		$(model_name) $(model_rep) \
		$(iwasaki_expt_dir) $(iwasaki_sam_file) \
		$(human_gene_len_file) $(human_gene_seq_file) \
		$(iwasaki_27_30_tr_bounds) $(iwasaki_27_30_te_bounds) \
		$(iwasaki_27_30_outputs) 10

$(iwasaki_28mer_pattern): \
		| $(iwasaki_sam_file) \
		$(iwasaki_28_proc_sam_files) \
		$(human_gene_len_file) $(human_gene_seq_file) \
		$(iwasaki_28_tr_bounds) $(iwasaki_28_te_bounds) \
		$(iwasaki_28_outputs) $(iwasaki_nn_dir) \
		$(iwasaki_expt_dir)
	$(eval model_name=s28_$*)
	echo
	echo ...Making model $(model_name)
	echo
	python $(repro_dir)/28mer_models.py \
		$(model_name) \
		$(iwasaki_expt_dir) $(iwasaki_sam_file) \
		$(human_gene_len_file) $(human_gene_seq_file) \
		$(iwasaki_28_tr_bounds) $(iwasaki_28_te_bounds) \
		$(iwasaki_28_outputs) 10

# Make iwasaki results directory and files
$(iwasaki_results_dir): | $(results_dir)
	mkdir $(iwasaki_results_dir)

$(iwasaki_results_feat_neighborhood_dir): | $(iwasaki_results_dir)
	mkdir $(iwasaki_results_feat_neighborhood_dir)

$(iwasaki_results_leaveout_dir): | $(iwasaki_results_dir)
	mkdir $(iwasaki_results_leaveout_dir)

$(iwasaki_results_full_dir): | $(iwasaki_results_dir)
	mkdir $(iwasaki_results_full_dir)

$(iwasaki_results_full_epoch_dir): | $(iwasaki_results_full_dir)
	mkdir $(iwasaki_results_full_epoch_dir)

$(iwasaki_results_28_dir): | $(iwasaki_results_dir)
	mkdir $(iwasaki_results_28_dir)

$(iwasaki_results_28_epoch_dir): | $(iwasaki_results_28_dir)
	mkdir $(iwasaki_results_28_epoch_dir)

#NOTE: We put an example .pkl file as a dependency so we can be sure this
# epoch is made. We can't put the epoch file as a dependency bc there's no rule
# to create it, and then make can't place the targets properly in its DAG. 

$(iwasaki_full_codon_scores_results_files_pattern): \
		| $(iwasaki_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/epoch10/te_cost_by_epoch.pkl \
		$(iwasaki_results_full_epoch_dir)
	python $(repro_dir)/codon_scores.py \
		$(iwasaki_nn_dir)/full_cod_n7p5_nt_n21p17_rep0 10
	cp $(iwasaki_full_analysis_epoch_dir)/codon_scores.tsv \
		$(iwasaki_results_full_epoch_dir)/codon_scores.tsv
	cp $(iwasaki_full_analysis_epoch_dir)/codon_scores_colormap.pdf \
		$(iwasaki_results_full_epoch_dir)/codon_scores_colormap.pdf

#NOTE: We put an example .pkl file as a dependency so we can be sure this
# epoch is made. We can't put the epoch file as a dependency bc there's no
# rule to create it, and then make can't place the targets properly in its DAG. 

$(iwasaki_28_codon_scores_results_files_pattern): \
		| $(iwasaki_nn_dir)/s28_cod_n7p5_nt_n21p17/epoch10/te_cost_by_epoch.pkl \
		$(iwasaki_results_28_epoch_dir)
	python $(repro_dir)/codon_scores.py \
		$(iwasaki_nn_dir)/s28_cod_n7p5_nt_n21p17 10
	cp $(iwasaki_28_analysis_epoch_dir)/codon_scores.tsv \
		$(iwasaki_results_28_epoch_dir)/codon_scores.tsv
	cp $(iwasaki_28_analysis_epoch_dir)/codon_scores_colormap.pdf \
		$(iwasaki_results_28_epoch_dir)/codon_scores_colormap.pdf

$(iwasaki_leaveout_mses_file) : \
		| $(iwasaki_leaveout_series_files) \
		$(iwasaki_results_leaveout_dir)
	python $(repro_dir)/aggregate_mses.py aggregate_w_diffs \
		$(iwasaki_nn_dir) 10 10 \
		$(iwasaki_leaveout_mses_file) full_cod_n7p5_nt_n21p17 \
		$(addsuffix _cod_n7p5_nt_n21p17, \
			$(addprefix nocod,$(leaveout_cods)))

$(iwasaki_full_plot_pattern): \
		| $(repro_dir)/plot_nn.py \
		$(iwasaki_results_full_epoch_dir) \
		$(iwasaki_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/epoch10/te_cost_by_epoch.pkl
	echo
	echo "Plotting nn full_cod_n7p5_n21p17_rep0"
	echo
	python $(repro_dir)/plot_nn.py \
		$(iwasaki_expt_dir) full_cod_n7p5_nt_n21p17_rep0 10
	cp -r $(iwasaki_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/plots \
		$(iwasaki_results_full_dir)

$(iwasaki_28_plot_pattern): \
		| $(repro_dir)/plot_nn.py \
		$(iwasaki_results_28_epoch_dir) \
		$(iwasaki_nn_dir)/s28_cod_n7p5_nt_n21p17/epoch10/te_cost_by_epoch.pkl 
	echo
	echo "Plotting nn s28_cod_n7p5_n21p17_rep0"
	echo
	python $(repro_dir)/plot_nn.py \
		$(iwasaki_expt_dir) s28_cod_n7p5_nt_n21p17 10
	cp -r $(iwasaki_nn_dir)/s28_cod_n7p5_nt_n21p17/plots \
		$(iwasaki_results_28_dir)

iwasaki_clean: 
	rm $(iwasaki_mapped_sam_file)

##########################
### Green Experiment
##########################

green_expt: $(green_sam_file) ${green_plot_files} \
		$(green_27_29_proc_sam_files) \
		$(green_28_proc_sam_files) \
		$(green_full_model_files) \
		$(green_leaveout_series_files) \
		$(green_28mer_files) \
		$(green_full_codon_scores_results_files) \
		$(green_28_codon_scores_results_files) \
		$(green_leaveout_mses_file) \
		$(green_full_plot_files) \
		$(green_28_plot_files) \
		| $(green_expt_dir) $(green_subdirs)

$(green_expt_dir): | $(expts_dir)
	mkdir $(green_expt_dir)

$(green_proc_dir): | $(green_expt_dir)
	mkdir $(green_proc_dir)

$(green_plot_dir): | $(green_expt_dir)
	mkdir $(green_plot_dir)

$(green_nn_dir): | $(green_expt_dir)
	mkdir $(green_nn_dir)

$(green_lr_dir): | $(green_expt_dir)
	mkdir $(green_lr_dir)

#$(green_raw_sam_file): | $(green_proc_dir)
#	samtools view /mnt/lareaulab/lareau/YeastFootprints/WuGreen2017/wu_green_wt/wu_green_wt.transcript.bam > $(green_raw_sam_file)

$(green_raw_fastq_files): \
		| $(green_proc_dir) 
	fastq-dump SRR5008134 -O $(green_proc_dir)
	fastq-dump SRR5008135 -O $(green_proc_dir)

#NOTE: This is missing a bunch of prereqs
$(green_raw_sam_file): \
		| $(green_raw_fastq_files) \
		$(green_trimmer_script) \
		$(green_mapping_script) \
		$(green_proc_dir)
	bash $(green_mapping_script) $(genome_dir) $(green_proc_dir)

$(green_sam_file): $(green_raw_sam_file) \
		| $(green_proc_dir) $(green_expt_dir)
	python $(repro_dir)/process_data.py edit_sam_file green \
		$(green_expt_dir) $(green_raw_sam_file)

${green_plot_files_pattern}: \
		| $(green_sam_file) $(yeast_gene_len_file) \
		${green_plot_dir} $(green_expt_dir)
	python $(repro_dir)/process_data.py size_and_frame_analysis green \
		$(green_expt_dir) $(green_sam_file) \
		$(yeast_gene_len_file)

$(green_27_29_proc_sam_pattern): \
		| $(green_sam_file) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(green_proc_dir) $(green_expt_dir)
	python $(repro_dir)/process_data.py process_sam_file green \
		$(green_expt_dir) $(green_sam_file) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(yeast_paralogs_file)

$(green_28_proc_sam_pattern): \
                | $(green_sam_file) \
                $(yeast_gene_len_file) $(yeast_gene_seq_file) \
                $(green_proc_dir) $(green_expt_dir)
	python $(repro_dir)/process_data.py process_sam_file_28mers green \
		$(green_expt_dir) $(green_sam_file) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(yeast_paralogs_file)

$(green_full_model_pattern): \
		| $(green_sam_file) \
		$(green_27_29_proc_sam_files) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(green_27_29_tr_bounds) $(green_27_29_te_bounds) \
		$(green_27_29_outputs) $(green_nn_dir) \
		$(green_expt_dir)
	$(eval model_name=$(firstword $(subst _rep, ,$*)))
	$(eval model_rep=$(lastword $(subst _rep, ,$*)))
	echo
	echo ...Making model $(model_name) rep $(model_rep)
	echo
	python $(repro_dir)/feat_neighborhood_nn_series.py \
		$(model_name) $(model_rep) \
		$(green_expt_dir) $(green_sam_file) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(green_27_29_tr_bounds) $(green_27_29_te_bounds) \
		$(green_27_29_outputs) 20

$(green_leaveout_series_pattern): \
		| $(green_sam_file) \
		$(green_27_29_proc_sam_files) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(green_27_29_tr_bounds) $(green_27_29_te_bounds) \
		$(green_27_29_outputs) $(green_nn_dir) \
		$(green_expt_dir)
	$(eval model_name=nocod$(firstword $(subst _rep, ,$*)))
	$(eval model_rep=$(lastword $(subst _rep, ,$*)))
	echo
	echo ...Making model $(model_name) rep $(model_rep)
	echo
	python $(repro_dir)/leaveout_series.py \
		$(model_name) $(model_rep) \
		$(green_expt_dir) $(green_sam_file) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(green_27_29_tr_bounds) $(green_27_29_te_bounds) \
		$(green_27_29_outputs) 20

$(green_28mer_pattern): \
		| $(green_sam_file) \
		$(green_28_proc_sam_files) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(green_28_tr_bounds) $(green_28_te_bounds) \
		$(green_28_outputs) $(green_nn_dir) \
		$(green_expt_dir)
	$(eval model_name=s28_$*)
	echo
	echo ...Making model $(model_name)
	echo
	python $(repro_dir)/28mer_models.py \
		$(model_name) \
		$(green_expt_dir) $(green_sam_file) \
		$(yeast_gene_len_file) $(yeast_gene_seq_file) \
		$(green_28_tr_bounds) $(green_28_te_bounds) \
		$(green_28_outputs) 10

# Make green results directory and files
$(green_results_dir): | $(results_dir)
	mkdir $(green_results_dir)

$(green_results_feat_neighborhood_dir): | $(green_results_dir)
	mkdir $(green_results_feat_neighborhood_dir)

$(green_results_leaveout_dir): | $(green_results_dir)
	mkdir $(green_results_leaveout_dir)

$(green_results_full_dir): | $(green_results_dir)
	mkdir $(green_results_full_dir)

$(green_results_full_epoch_dir): | $(green_results_full_dir)
	mkdir $(green_results_full_epoch_dir)

$(green_results_28_dir): | $(green_results_dir)
	mkdir $(green_results_28_dir)

$(green_results_28_epoch_dir): | $(green_results_28_dir)
	mkdir $(green_results_28_epoch_dir)

#NOTE: We put an example .pkl file as a dependency so we can be sure this
# epoch is made. We can't put the epoch file as a dependency bc there's no
# rule to create it, and then make can't place the targets properly in its DAG. 

$(green_full_codon_scores_results_files_pattern): \
		| $(green_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/epoch20/te_cost_by_epoch.pkl \
		$(green_results_full_epoch_dir)
	python $(repro_dir)/codon_scores.py \
		$(green_nn_dir)/full_cod_n7p5_nt_n21p17_rep0 20
	cp $(green_full_analysis_epoch_dir)/codon_scores.tsv \
		$(green_results_full_epoch_dir)/codon_scores.tsv
	cp $(green_full_analysis_epoch_dir)/codon_scores_colormap.pdf \
		$(green_results_full_epoch_dir)/codon_scores_colormap.pdf

#NOTE: We put an example .pkl file as a dependency so we can be sure this
# epoch is made. We can't put the epoch file as a dependency bc there's no
# rule to create it, and then make can't place the targets properly in its DAG. 

$(green_28_codon_scores_results_files_pattern): \
		| $(green_nn_dir)/s28_cod_n7p5_nt_n21p17/epoch10/te_cost_by_epoch.pkl \
		$(green_results_28_epoch_dir)
	python $(repro_dir)/codon_scores.py \
		$(green_nn_dir)/s28_cod_n7p5_nt_n21p17 10
	cp $(green_28_analysis_epoch_dir)/codon_scores.tsv \
		$(green_results_28_epoch_dir)/codon_scores.tsv
	cp $(green_28_analysis_epoch_dir)/codon_scores_colormap.pdf \
		$(green_results_28_epoch_dir)/codon_scores_colormap.pdf

$(green_leaveout_mses_file) : \
		| $(green_leaveout_series_files) \
		$(green_results_leaveout_dir)
	python $(repro_dir)/aggregate_mses.py aggregate_w_diffs \
		$(green_nn_dir) 20 10 \
		$(green_leaveout_mses_file) full_cod_n7p5_nt_n21p17 \
		$(addsuffix _cod_n7p5_nt_n21p17, \
			$(addprefix nocod,$(leaveout_cods)))

$(green_full_plot_pattern): \
		| $(repro_dir)/plot_nn.py \
		$(green_results_full_epoch_dir) \
		$(green_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/epoch20/te_cost_by_epoch.pkl
	echo
	echo "Plotting nn full_cod_n7p5_n21p17_rep0"
	echo
	python $(repro_dir)/plot_nn.py \
		$(green_expt_dir) full_cod_n7p5_nt_n21p17_rep0 20
	cp -r $(green_nn_dir)/full_cod_n7p5_nt_n21p17_rep0/plots \
		$(green_results_full_dir)

$(green_28_plot_pattern): \
		| $(repro_dir)/plot_nn.py \
		$(green_results_28_epoch_dir) \
		$(green_nn_dir)/s28_cod_n7p5_nt_n21p17/epoch10/te_cost_by_epoch.pkl
	echo
	echo "Plotting nn s28_cod_n7p5_n21p17_rep0"
	echo
	python $(repro_dir)/plot_nn.py \
		$(green_expt_dir) s28_cod_n7p5_nt_n21p17 10
	cp -r $(green_nn_dir)/s28_cod_n7p5_nt_n21p17/plots \
		$(green_results_28_dir)

green_clean: 
	rm $(green_mapped_sam_file)

##############
### Figures
##############

figures: $(fig_files) | $(fig_dir)

$(facs_data_file): | $(facs_data_zip_file)
	gunzip -c $(facs_data_zip_file) > $(facs_data_file)

#NOTE: Update first two parameter entries
$(fig_dir)/figure_1B_scaledcts.pdf: \
		| $(weinberg_27_31_te_data_table) \
		$(weinberg_results_final_model_y_te_hat) \
		$(repro_dir)/figure_1B_scaledcts.R $(fig_dir)
	Rscript $(repro_dir)/figure_1B_scaledcts.R \
		YLR044C 40 \
		$(weinberg_27_31_te_data_table) \
		$(fig_dir)/figure_1B_scaledcts.pdf

$(fig_dir)/figure_1C_mse.pdf: \
		| $(weinberg_feat_nb_mses_file) \
		$(weinberg_feat_nb_linreg_mses_file) \
		$(weinberg_struc_mses_file) \
		$(repro_dir)/figure_1C_mse.R $(fig_dir)
	Rscript $(repro_dir)/figure_1C_mse.R \
		$(weinberg_feat_nb_mses_file) \
		$(weinberg_feat_nb_linreg_mses_file) \
		$(weinberg_struc_mses_file) \
		$(fig_dir)/figure_1C_mse.pdf

#NOTE: Think about where to put the y_te/y_te_hat files
$(fig_dir)/figure_1D_scatter.pdf: \
		| $(weinberg_results_final_model_y_te) \
		$(weinberg_results_final_model_y_te_hat) \
		$(repro_dir)/figure_1D_scatter.R $(fig_dir)
	Rscript $(repro_dir)/figure_1D_scatter.R \
		$(weinberg_results_final_model_y_te) \
		$(weinberg_results_final_model_y_te_hat) \
		$(fig_dir)/figure_1D_scatter.pdf

#NOTE: Update first parameter entries
$(fig_dir)/figure_1E_indiv_gene.pdf: \
		| $(weinberg_27_31_te_data_table) \
		$(weinberg_results_final_model_y_te_hat) \
		$(yeast_gene_symbol_file) \
		$(repro_dir)/figure_1E_indiv_gene.R $(fig_dir)
	Rscript $(repro_dir)/figure_1E_indiv_gene.R \
		YLR044C \
		$(weinberg_27_31_te_data_table) \
		$(weinberg_results_final_model_y_te_hat) \
		$(yeast_gene_symbol_file) \
		$(fig_dir)/figure_1E_indiv_gene.pdf

$(fig_dir)/figure_1F_binned_error.pdf: \
		| $(weinberg_results_final_model_y_te) \
		$(weinberg_results_final_model_y_te_hat) \
		$(repro_dir)/figure_1F_binned_error.R $(fig_dir)
	Rscript $(repro_dir)/figure_1F_binned_error.R \
		$(weinberg_results_final_model_y_te) \
		$(weinberg_results_final_model_y_te_hat) \
		$(fig_dir)/figure_1F_binned_error.pdf

$(fig_dir)/figure_2A_codonmse.pdf: \
		| $(weinberg_leaveout_mses_file) \
		$(repro_dir)/figure_2A_codonmse.R $(fig_dir)
	Rscript $(repro_dir)/figure_2A_codonmse.R $(weinberg_leaveout_mses_file) \
		$(fig_dir)/figure_2A_codonmse.pdf

$(fig_dir)/figure_2B_heatmap.pdf: \
		| $(weinberg_results_full_epoch_dir)/codon_scores.tsv \
		$(repro_dir)/figure_2B_heatmap.R $(fig_dir)
	Rscript $(repro_dir)/figure_2B_heatmap.R \
		$(weinberg_results_full_epoch_dir)/codon_scores.tsv \
		$(fig_dir)/figure_2B_heatmap.pdf

#Check that liana wants this without struc features
$(fig_dir)/figure_2C_tai.pdf: \
		| $(weinberg_results_full_epoch_dir)/codon_scores.tsv \
		$(yeast_codon_props_file) \
		$(repro_dir)/figure_2C_tai.R $(fig_dir)
	Rscript $(repro_dir)/figure_2C_tai.R \
		$(weinberg_results_full_epoch_dir)/codon_scores.tsv \
		$(yeast_codon_props_file) \
		$(fig_dir)/figure_2C_tai.pdf

#Check that liana wants this without struc features
$(fig_dir)/figure_2D_wobble.pdf: \
		| $(weinberg_results_full_epoch_dir)/codon_scores.tsv \
		$(yeast_codon_anticodon_file) \
		$(repro_dir)/figure_2D_wobble.R $(fig_dir)
	Rscript $(repro_dir)/figure_2D_wobble.R \
		$(weinberg_results_full_epoch_dir)/codon_scores.tsv \
		$(yeast_codon_anticodon_file) \
		$(fig_dir)/figure_2D_wobble.pdf

#Check that liana wants this for full models
$(fig_dir)/figure_2E_lareaucodonmse.pdf: \
		| $(lareau_leaveout_mses_file) \
		$(repro_dir)/figure_2E_lareaucodonmse.R $(fig_dir)
	Rscript $(repro_dir)/figure_2E_lareaucodonmse.R $(lareau_leaveout_mses_file) \
		$(fig_dir)/figure_2E_lareaucodonmse.pdf

#Check that liana wants this for 28mer models
$(fig_dir)/figure_2F_5prime.pdf: \
		| $(lareau_results_28_epoch_dir)/codon_scores.tsv \
		$(iwasaki_results_28_epoch_dir)/codon_scores.tsv \
		$(repro_dir)/figure_2F_5prime.R $(fig_dir)
	Rscript $(repro_dir)/figure_2F_5prime.R \
		$(lareau_results_28_epoch_dir)/codon_scores.tsv \
		$(iwasaki_results_28_epoch_dir)/codon_scores.tsv \
		$(fig_dir)/figure_2F_5prime.pdf

#Check that liana wants this for 28mer models
$(fig_dir)/figure_2G_asite.pdf: \
		| $(lareau_results_28_epoch_dir)/codon_scores.tsv \
		$(iwasaki_results_28_epoch_dir)/codon_scores.tsv \
		$(repro_dir)/figure_2G_asite.R $(fig_dir)
	Rscript $(repro_dir)/figure_2G_asite.R \
		$(lareau_results_28_epoch_dir)/codon_scores.tsv \
		$(iwasaki_results_28_epoch_dir)/codon_scores.tsv \
		$(fig_dir)/figure_2G_asite.pdf

$(fig_dir)/figure_2H_cl2.pdf: \
		| $(lareau_results_28_epoch_dir)/codon_scores.tsv \
		$(green_results_28_epoch_dir)/codon_scores.tsv \
		$(wetlab_dir)/circligase_qpcr.csv \
		$(repro_dir)/figure_2H_cl2.R $(fig_dir)
	Rscript $(repro_dir)/figure_2H_cl2.R \
		$(lareau_results_28_epoch_dir)/codon_scores.tsv \
		$(green_results_28_epoch_dir)/codon_scores.tsv \
		$(wetlab_dir)/circligase_qpcr.csv \
		$(fig_dir)/figure_2H_cl2.pdf

$(fig_dir)/figure_2I_cl1.pdf: \
		| $(lareau_results_28_epoch_dir)/codon_scores.tsv \
		$(green_results_28_epoch_dir)/codon_scores.tsv \
		$(wetlab_dir)/circligase_qpcr.csv \
		$(repro_dir)/figure_2I_cl1.R $(fig_dir)
	Rscript $(repro_dir)/figure_2I_cl1.R \
		$(lareau_results_28_epoch_dir)/codon_scores.tsv \
		$(green_results_28_epoch_dir)/codon_scores.tsv \
		$(wetlab_dir)/circligase_qpcr.csv \
		$(fig_dir)/figure_2I_cl1.pdf

$(fig_dir)/figure_3B_citrine_dist.pdf: \
		| $(repro_dir)/figure_3B_citrine_dist.R \
		$(paper_data_dir)/random_citrine_score_dist.txt \
		$(paper_data_dir)/citrine_construct_scores.txt 
	Rscript $(repro_dir)/figure_3B_citrine_dist.R \
		$(paper_data_dir)/random_citrine_score_dist.txt \
		$(paper_data_dir)/citrine_construct_scores.txt \
		$(fig_dir)/figure_3B_citrine_dist.pdf

$(fig_dir)/figure_3C_facs.pdf: \
		| $(repro_dir)/figure_3C_facs.R \
		$(facs_data_file) \
		$(paper_data_dir)/citrine_construct_scores.txt $(fig_dir)
	Rscript $(repro_dir)/figure_3C_facs.R \
		$(facs_data_file) \
		$(paper_data_dir)/citrine_construct_scores.txt \
		$(fig_dir)/figure_3C_facs.pdf

$(fig_dir)/figure_3D_te.pdf: \
		| $(repro_dir)/figure_3D_te.R \
		$(paper_data_dir)/citrine_construct_scores.txt \
		$(mrna_data_file) \
		$(facs_data_file)
	Rscript $(repro_dir)/figure_3D_te.R \
		$(paper_data_dir)/citrine_construct_scores.txt \
		$(mrna_data_file) \
		$(facs_data_file) \
		$(fig_dir)/figure_3D_te.pdf

$(fig_dir)/supp_table_codon_scores.csv: \
		| $(repro_dir)/supp_table_codon_scores.py \
		$(weinberg_results_full_epoch_dir)/codon_scores.tsv \
		$(fig_dir)
	python $(repro_dir)/supp_table_codon_scores.py \
		$(weinberg_results_full_epoch_dir)/codon_scores.tsv \
		$(fig_dir)/supp_table_codon_scores.csv

$(fig_dir)/supp_figure_mrna_qpcr.pdf: \
		| $(repro_dir)/supp_figure_mrna_qpcr.R \
		$(paper_data_dir)/citrine_construct_scores.txt \
		$(mrna_data_file) 
	Rscript $(repro_dir)/supp_figure_mrna_qpcr.R \
		$(paper_data_dir)/citrine_construct_scores.txt \
		$(mrna_data_file) \
		$(fig_dir)/supp_figure_mrna_qpcr.pdf

$(fig_dir)/supp_figure_facs.pdf: \
		| $(repro_dir)/supp_figure_facs.R \
		$(facs_data_file) \
		$(paper_data_dir)/citrine_construct_scores.txt
	Rscript $(repro_dir)/supp_figure_facs.R \
		$(facs_data_file) \
		$(paper_data_dir)/citrine_construct_scores.txt \
		$(fig_dir)/supp_figure_facs.pdf

$(fig_dir)/supp_figure_greenmse.pdf: \
		| $(green_leaveout_mses_file) \
		$(repro_dir)/supp_figure_greenmse.R $(fig_dir)
	Rscript $(repro_dir)/supp_figure_greenmse.R \
		$(green_leaveout_mses_file) \
		$(fig_dir)/supp_figure_greenmse.pdf

$(fig_dir)/supp_figure_iwasakimse.pdf: \
		| $(iwasaki_leaveout_mses_file) \
		$(repro_dir)/supp_figure_iwasakimse.R $(fig_dir)
	Rscript $(repro_dir)/supp_figure_iwasakimse.R \
		$(iwasaki_leaveout_mses_file) \
		$(fig_dir)/supp_figure_iwasakimse.pdf

$(fig_dir)/supp_figure_cl1v2.pdf: \
		| $(lareau_results_28_epoch_dir)/codon_scores.tsv \
		$(green_results_28_epoch_dir)/codon_scores.tsv \
		$(repro_dir)/supp_figure_cl1v2.R $(fig_dir)
	Rscript $(repro_dir)/supp_figure_cl1v2.R \
		$(lareau_results_28_epoch_dir)/codon_scores.tsv \
		$(green_results_28_epoch_dir)/codon_scores.tsv \
		$(fig_dir)/supp_figure_cl1v2.pdf

#################
### Paper data
#################

paper_data: \
	$(paper_data_dir) \
	$(paper_data_fname) \

$(paper_data_dir): | $(results_dir)
	mkdir $(paper_data_dir)

$(paper_data_fname): \
		| $(weinberg_nn_dir)/str_n17n15_cod_n7p5_nt_n21p17_rep0/epoch30/te_cost_by_epoch.pkl \
		$(weinberg_feat_nb_mses_file) \
		$(weinberg_leaveout_mses_file) \
		$(lareau_results_28_epoch_dir)/codon_scores.tsv \
		$(iwasaki_results_28_epoch_dir)/codon_scores.tsv \
		$(repro_dir)/paper_data.py \
		$(paper_data_dir)
	python $(repro_dir)/paper_data.py \
		$(weinberg_nn_dir)/str_n17n15_cod_n7p5_nt_n21p17_rep0 30 \
		$(weinberg_feat_nb_mses_file) \
		$(weinberg_leaveout_mses_file) \
		$(weinberg_struc_mses_file) \
		$(weinberg_results_full_epoch_dir)/codon_scores.tsv \
		$(lareau_results_28_epoch_dir)/codon_scores.tsv \
		$(iwasaki_results_28_epoch_dir)/codon_scores.tsv \
		$(paper_data_fname)
