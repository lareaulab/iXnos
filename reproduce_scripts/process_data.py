import iXnos.interface as inter
import sys

if __name__ == "__main__":
    all_actions = [
        "edit_sam_file", "size_and_frame_analysis", 
        "process_sam_file", "process_sam_file_28mers"]
    action = sys.argv[1]
    expt = sys.argv[2]
    expt_dir = sys.argv[3]
    sam_fname = sys.argv[4]
    if action in all_actions[1:4]:
        gene_len_fname = sys.argv[5]
    if action in all_actions[2:4]:
        gene_seq_fname = sys.argv[6]
    if action in all_actions[2:4]:
        if len(sys.argv) >= 8:
            paralog_groups_fname = sys.argv[7]
        else: 
            paralog_groups_fname = False
            print "Warning: May want to input paralog groups file to select" +\
                  " for most abundant paralog"
 

    if expt == "weinberg":
        cod_trunc_5p = 20
        cod_trunc_3p = 20
        shift_dict = {
            27:{0:False, 1:14, 2:False}, 28:{0:15, 1:False, 2:False}, 
            29:{0:15, 1:False, 2:16}, 30:{0:15, 1:False, 2:16}, 
            31:{0:15, 1:False, 2:16}}
        min_fp_size = 27
        max_fp_size = 31
        num_tr_genes = 333
        num_te_genes = 167
        min_cts_per_gene = 200
        min_cod_w_data = 100

        if action == "edit_sam_file":
            inter.edit_sam_file(
                expt_dir, sam_fname, filter_unmapped=True, 
                sam_add_simple_map_wts=False, RSEM_add_map_wts=True)

        elif action == "size_and_frame_analysis":
            inter.do_frame_and_size_analysis(
                expt_dir, sam_fname, gene_len_fname, min_size=13, 
                max_size=36, verbose=False)

        elif action == "process_sam_file":
            inter.process_sam_file(
                expt_dir, sam_fname, gene_seq_fname, gene_len_fname, 
                shift_dict, cod_trunc_5p, cod_trunc_3p, min_fp_size, 
                max_fp_size, num_tr_genes, num_te_genes, min_cts_per_gene, 
                min_cod_w_data, paralog_groups_fname=paralog_groups_fname, 
                overwrite=False)

        elif action == "process_sam_file_28mers":
            min_fp_size = 28
            max_fp_size = 28
            inter.process_sam_file(
                expt_dir, sam_fname, gene_seq_fname, gene_len_fname, 
                shift_dict, cod_trunc_5p, cod_trunc_3p, min_fp_size, 
                max_fp_size, num_tr_genes, num_te_genes, min_cts_per_gene, 
                min_cod_w_data, paralog_groups_fname=paralog_groups_fname,
                overwrite=False)

    elif expt == "lareau":
        cod_trunc_5p = 20
        cod_trunc_3p = 20
        shift_dict = {27:{0:15, 1:14, 2:False}, 28:{0:15, 1:False, 2:False}, 
            29:{0:15, 1:False, 2:16}}
        min_fp_size = 27
        max_fp_size = 29
        num_tr_genes = 333
        num_te_genes = 167
        min_cts_per_gene = 200
        min_cod_w_data = 100

        if action == "edit_sam_file":
            inter.edit_sam_file(
                expt_dir, sam_fname, filter_unmapped=True,
                sam_add_simple_map_wts=False, RSEM_add_map_wts=True)

        elif action == "size_and_frame_analysis":
            inter.do_frame_and_size_analysis(
                expt_dir, sam_fname, gene_len_fname, min_size=13,
                max_size=36, verbose=False)

        elif action == "process_sam_file":
            inter.process_sam_file(
                expt_dir, sam_fname, gene_seq_fname, gene_len_fname,
                shift_dict, cod_trunc_5p, cod_trunc_3p, min_fp_size,
                max_fp_size, num_tr_genes, num_te_genes, min_cts_per_gene,
                min_cod_w_data, paralog_groups_fname=paralog_groups_fname, 
                overwrite=False)

        elif action == "process_sam_file_28mers":
            min_fp_size = 28
            max_fp_size = 28
            inter.process_sam_file(
                expt_dir, sam_fname, gene_seq_fname, gene_len_fname, 
                shift_dict, cod_trunc_5p, cod_trunc_3p, min_fp_size, 
                max_fp_size, num_tr_genes, num_te_genes, min_cts_per_gene, 
                min_cod_w_data, paralog_groups_fname=paralog_groups_fname, 
                overwrite=False)

    elif expt == "iwasaki":
        cod_trunc_5p = 20
        cod_trunc_3p = 20
        shift_dict = {27:{0:15, 1:14, 2:False}, 28:{0:15, 1:False, 2:False}, 
            29:{0:15, 1:False, 2:False}, 30:{0:15, 1:False, 2:16}}
        min_fp_size = 27
        max_fp_size = 30
        num_tr_genes = 333
        num_te_genes = 167
        min_cts_per_gene = 200
        min_cod_w_data = 100

        if action == "edit_sam_file":
            inter.edit_sam_file(
                expt_dir, sam_fname, filter_unmapped=True,
                sam_add_simple_map_wts=False, RSEM_add_map_wts=True)

        elif action == "size_and_frame_analysis":
            inter.do_frame_and_size_analysis(
                expt_dir, sam_fname, gene_len_fname, min_size=13,
                max_size=36, verbose=False)

        elif action == "process_sam_file":
            inter.process_sam_file(
                expt_dir, sam_fname, gene_seq_fname, gene_len_fname,
                shift_dict, cod_trunc_5p, cod_trunc_3p, min_fp_size,
                max_fp_size, num_tr_genes, num_te_genes, min_cts_per_gene,
                min_cod_w_data, paralog_groups_fname=paralog_groups_fname, 
                overwrite=False)

        elif action == "process_sam_file_28mers":
            min_fp_size = 28
            max_fp_size = 28
            num_tr_genes = 200
            num_te_genes = 100
            inter.process_sam_file(
                expt_dir, sam_fname, gene_seq_fname, gene_len_fname, 
                shift_dict, cod_trunc_5p, cod_trunc_3p, min_fp_size, 
                max_fp_size, num_tr_genes, num_te_genes, min_cts_per_gene, 
                min_cod_w_data, paralog_groups_fname=paralog_groups_fname,
                overwrite=False)
    
    elif expt == "green":
        cod_trunc_5p = 20
        cod_trunc_3p = 20
        shift_dict = {27:{0:15, 1:14, 2:False}, 28:{0:15, 1:False, 2:False},
            29:{0:15, 1:False, 2:16}}
        min_fp_size = 27
        max_fp_size = 29
        num_tr_genes = 333
        num_te_genes = 167
        min_cts_per_gene = 200
        min_cod_w_data = 100

        if action == "edit_sam_file":
            inter.edit_sam_file(
                expt_dir, sam_fname, filter_unmapped=True,
                sam_add_simple_map_wts=False, RSEM_add_map_wts=True)

        elif action == "size_and_frame_analysis":
            inter.do_frame_and_size_analysis(
                expt_dir, sam_fname, gene_len_fname, min_size=13,
                max_size=36, verbose=False)

        elif action == "process_sam_file":
            inter.process_sam_file(
                expt_dir, sam_fname, gene_seq_fname, gene_len_fname,
                shift_dict, cod_trunc_5p, cod_trunc_3p, min_fp_size,
                max_fp_size, num_tr_genes, num_te_genes, min_cts_per_gene,
                min_cod_w_data, paralog_groups_fname=paralog_groups_fname, 
                overwrite=False)

        elif action == "process_sam_file_28mers":
            min_fp_size = 28
            max_fp_size = 28
            inter.process_sam_file(
                expt_dir, sam_fname, gene_seq_fname, gene_len_fname, 
                shift_dict, cod_trunc_5p, cod_trunc_3p, min_fp_size, 
                max_fp_size, num_tr_genes, num_te_genes, min_cts_per_gene, 
                min_cod_w_data, paralog_groups_fname=paralog_groups_fname, 
                overwrite=False)

