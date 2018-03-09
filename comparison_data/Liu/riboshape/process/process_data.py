# Global parameters
cts_by_codon_fname = "cts_by_codon.size.27.31.txt"
gene_seq_fname = "scer.transcripts.13cds10.fa"
trunc_len_5p = 20
trunc_len_3p = 20
pad_len_5p = 13
pad_len_3p = 10
min_trunc_len = 50
min_cts_per_gene = 200
min_cod_w_data = 100

if __name__ == "__main__":
    # Open input files
    cts_by_codon_file = open("cts_by_codon.size.27.31.txt", "r")
    gene_seq_file = open("scer.transcripts.13cds10.fa", "r")

    # Open output files
    GeneName_file = open("GeneName.txt", "w")
    Asitecount_file = open("Asitecount.txt", "w")
    asite_density_file = open("asite_density.txt", "w")
    CDS_file = open("CDS.txt", "w")
    codon_file = open("codon.txt", "w")
    distribution_file = open("distribution.txt", "w")
    
    # Write codon file
    codons = [x + y + z for x in "ACGT" for y in "ACGT" for z in "ACGT"]
    cod2idx = {cod:i for cod, i in zip(codons, range(len(codons)))}
    for codon in codons: 
        codon_file.write(codon + "\n")

    # Data structure to remember canonical ordering of genes for Liu program
    gene2idx = {}
    idx_ctr = 0
    
    # Iterate over cts_by_codon file
    # Each line contains gene_name as first field, cts_by_codon as 
    #    successive fields
    for line in cts_by_codon_file:
        line = line.strip().split()
        # Process gene_name and cts by codon
        gene_name = line[0]
        gene_cts_by_codon = [float(elt) for elt in line[1:]]
        # Get length of gene in codons
        gene_num_cod = len(gene_cts_by_codon)
        # Check that length of gene is greater than truncations
        if gene_num_cod - trunc_len_5p - trunc_len_3p < min_trunc_len:
            continue
        # Check that gene has sufficient codons with data
        cod_w_data = sum([ct > 0 for ct in gene_cts_by_codon[trunc_len_5p:-trunc_len_3p]])
        if cod_w_data < min_cod_w_data:
            continue
        # Check that gene has sufficient counts after truncation
        gene_tot_cts = sum(gene_cts_by_codon[trunc_len_5p:-trunc_len_3p])
        if gene_tot_cts < min_cts_per_gene:
            continue
        # Assign gene index
        gene2idx[gene_name] = idx_ctr
        idx_ctr += 1
        # Truncate gene_cts_by_codon, compute new length
        trunc_cts_by_codon = \
            gene_cts_by_codon[trunc_len_5p:gene_num_cod-trunc_len_3p]
        trunc_num_cod = len(trunc_cts_by_codon)
        # Compute truncated gene density
        trunc_density = sum(trunc_cts_by_codon)/trunc_num_cod
        # Write to output files
        GeneName_file.write(gene_name + "\n")
        Asitecount_line = "\t".join([str(elt) for elt in trunc_cts_by_codon]) + "\n"
        Asitecount_file.write(Asitecount_line)
        asite_density_file.write(str(trunc_density) + "\n")

    # Get idx2gene dict
    idx2gene = {value:key for key, value in gene2idx.iteritems()}
        
    # Make dict for gene CDSs
    cds_dict = {}
    trunc_cds_dict = {}

    # Populate CDS dict
    first_gene = True
    for line in gene_seq_file:
        if line[0] == ">":
            if not first_gene:
                # NOTE: Be sure fasta gene seqs are longer than padding lens
                gene_cds = gene_seq[pad_len_5p:-pad_len_3p]
                trunc_cds = gene_cds[3*trunc_len_5p:-3*trunc_len_3p]
                cds_dict[gene_name] = gene_cds
                trunc_cds_dict[gene_name] = trunc_cds
            gene_name = line.strip().split()[0][1:]
            first_gene = False
            gene_seq = ""
        else:
            gene_seq += line.strip()
    gene_cds = gene_seq[pad_len_5p:-pad_len_3p]
    trunc_cds = gene_cds[3*trunc_len_5p:-3*trunc_len_3p]
    cds_dict[gene_name] = gene_cds
    trunc_cds_dict[gene_name] = trunc_cds
        
    # Get num genes long enough to truncate
    num_genes = len(gene2idx)

    for idx in range(num_genes):
        # Store gene_name and truncated CDS
        gene_name = idx2gene[idx]
        trunc_cds = trunc_cds_dict[gene_name]
        # Check length of truncated CDS
        if len(trunc_cds) % 3 != 0:
            print "Error! Truncated cds should have length multiple of 3"
            print gene_name, len(trunc_cds)
        # Compute CDS binary matrix representation
        num_cods = len(trunc_cds) / 3
        cds_matrix = [[0 for j in xrange(num_cods)] 
                         for i in xrange(len(codons))]
        for j in xrange(num_cods):
            cod = trunc_cds[3*j:3*(j+1)]
            cod_idx = cod2idx[cod]
            cds_matrix[cod_idx][j] = 1
        # Write to output files
        CDS_file.write(trunc_cds + "\n")
        for row in cds_matrix:
            line = " ".join([str(elt) for elt in row]) + "\n"
            distribution_file.write(line)
        distribution_file.write("---\n")
