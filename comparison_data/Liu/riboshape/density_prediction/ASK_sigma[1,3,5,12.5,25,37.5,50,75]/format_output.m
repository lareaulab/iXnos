dir_path = '/mnt/lareaulab/rtunney/iXnos/comparison_data/Liu';
ASK_path = strcat(dir_path, '/riboshape/density_prediction/ASK_sigma[1            3            5         12.5           25         37.5           50           75]');

raw_data_fname = strcat(dir_path, '/riboshape/data/data.liu_chx_genes.mat');
fname1 = strcat(ASK_path, '/asite_lengthmin100_lengthmax210.mat');
fname2 = strcat(ASK_path, '/asite_lengthmin211_lengthmax460.mat');
fname3 = strcat(ASK_path, '/asite_lengthmin461_lengthmax710.mat');
fname4 = strcat(ASK_path, '/asite_lengthmin711_lengthmax960.mat');
fname5 = strcat(ASK_path, '/asite_lengthmin961_lengthmax4871.mat');

out_fname = strcat(ASK_path, '/corrs_by_gene.txt');

% Load raw data
RAW_DATA = load(raw_data_fname);

% Open out file
out_fID = fopen(out_fname, 'w');

% Load trained model predictions
%for filename = {fname1, fname2, fname3, fname4, fname5}
for filename = {fname1, fname2, fname3, fname4, fname5}
    % Load trained model for a bin of genes by length
    TRAINED_DATA = load(char(filename));
    % Get number of genes and subspaces for bin
    NumBinGenes = length(TRAINED_DATA.GeneName);
    NumSubspaces = size(TRAINED_DATA.waveform, 2);
    % For each gene in bin
    for TrainIdx=1:NumBinGenes
        % Make string to write to out_file
        out_string = '';
        % Get name of gene add to out_string
        gene = TRAINED_DATA.GeneName{TrainIdx};
        out_string = strcat(out_string, gene);
        % Get index of gene in raw data set
        RawIdx = find(strcmp(RAW_DATA.GeneName, gene));
        % Get A site count distribution of gene in raw data set
        CtsByCodon = RAW_DATA.Asitecount{RawIdx};
        CodonDist = CtsByCodon/sum(RAW_DATA.Asitecount{RawIdx});
        for SubspaceIdx=1:NumSubspaces
            % Get correlation between CtsByCodon and regression fit in subspace
            pearson_r = corr(CodonDist, TRAINED_DATA.waveform{TrainIdx,SubspaceIdx}(:,2));
            if isnan(pearson_r)
                pearson_r = 'NA';
            end
            out_string = strcat(out_string, '\t', string(pearson_r));
        end
        out_string = strcat(out_string, '\n');
        fprintf(out_fID, out_string);
    end
end
