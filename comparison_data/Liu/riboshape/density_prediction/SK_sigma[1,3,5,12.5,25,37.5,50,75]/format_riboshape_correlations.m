dir_path = '/mnt/lareaulab/rtunney/iXnos/comparison_data/Liu';
SK_path = strcat(dir_path, '/riboshape/density_prediction/SK_sigma[1            3            5         12.5           25         37.5           50           75]');

raw_data_fname = strcat(dir_path, '/riboshape/data/data.liu_chx_genes.mat');

fname1 = strcat(SK_path, '/asite_lengthmin100_lengthmax210.mat');
fname2 = strcat(SK_path, '/asite_lengthmin211_lengthmax460.mat');
fname3 = strcat(SK_path, '/asite_lengthmin461_lengthmax710.mat');
fname4 = strcat(SK_path, '/asite_lengthmin711_lengthmax960.mat');
fname5 = strcat(SK_path, '/asite_lengthmin961_lengthmax4871.mat');

out_fname = strcat(SK_path, '/riboshape_corrs.txt');

for filename = {fname1, fname2, fname3, fname4, fname5}
    filename = char(filename);
    load(filename);
    filename_split = strsplit(filename, '.mat');
    corrs_outfname = strcat(filename_split(1), '.subspace_corrs.txt');
    gene_outfname = strcat(filename_split(1), '.gene_name.txt');
    writetable(cell2table(GeneName), char(gene_outfname))
    dlmwrite(char(corrs_outfname), correlation, '\t')
end
