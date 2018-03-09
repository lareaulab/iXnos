% Load Liu CHX genes
LiuChxGenes = load('../../riboshape/data/chxdata.mat', 'GeneName');
LiuChxGenes = LiuChxGenes.GeneName;

% Make codon data structure
% Open codon file
codon_fname = 'codon.txt';
codon_fID = fopen(codon_fname, 'r');
% Init codon cell and idx
codon = {};
i = 1;
% Load first line
line = fgetl(codon_fID);
% While line is a codon
while ischar(line)
    % Add new codon to data structure, increment index, load next codon
    codon{i} = line;
    i = i + 1;
    line = fgetl(codon_fID);
end
num_cods = i-1;
% Reshape codon cell to be a column
codon = reshape(codon, [num_cods,1]);
% Close codon file
fclose(codon_fID);

% Make asite_density data structure
asite_density_fname = 'asite_density.txt';
asite_density_fID = fopen(asite_density_fname, 'r');
formatSpec = '%f';
asite_density = fscanf(asite_density_fID, formatSpec);
fclose(asite_density_fID);

% Make GeneName data structure
% Open GeneName file
GeneName_fname = 'GeneName.txt';
GeneName_fID = fopen(GeneName_fname, 'r');
% Init GeneName cell and idx
GeneName = {};
i = 1;
% Load first line
line = fgetl(GeneName_fID);
% While line is a GeneName
while ischar(line)
    % Add new GeneName to data structure, increment index, load next GeneName
    GeneName{i} = line;
    i = i + 1;
    line = fgetl(GeneName_fID);
end
num_GeneName = i-1;
% Reshape GeneName cell to be a column
GeneName = reshape(GeneName, [num_GeneName,1]);
% Close codon file
fclose(GeneName_fID);

% Make CDS data structure
% Open CDS file
CDS_fname = 'CDS.txt';
CDS_fID = fopen(CDS_fname, 'r');
% Init CDS cell and idx
CDS = {};
i = 1;
% Load first line
line = fgetl(CDS_fID);
% While line is a CDS
while ischar(line)
    % Add new CDS to data structure, increment index, load next CDS
    CDS{i} = line;
    i = i + 1;
    line = fgetl(CDS_fID);
end
num_CDS = i-1;
% Reshape CDS cell to be a column
CDS = reshape(CDS, [num_CDS,1]);
% Close codon file
fclose(CDS_fID);

% Make distribution data structure
% Open distribution file
distribution_fname = 'distribution.txt';
distribution_fID = fopen(distribution_fname, 'r');
% Init distribution cell and idx
distribution = {};
i = 1;
% Load first CDS distribution
bin_vec = fscanf(distribution_fID, '%f');
sep_line = fgetl(distribution_fID);
% while sep_line is ---
while ischar(sep_line)
    num_cod = length(bin_vec)/64;
    bin_mat = transpose(reshape(bin_vec, [num_cod, 64]));
    distribution{i} = bin_mat;
    i = i + 1;
    bin_vec = fscanf(distribution_fID, '%f');
    sep_line = fgetl(distribution_fID);
end
num_genes = i-1;
% Reshape distribution cell to be a column
distribution = reshape(distribution, [num_genes,1]);
% Close codon file
fclose(distribution_fID);

% Make Asitecount data structure
% Open Asitecount file
Asitecount_fname = 'Asitecount.txt';
Asitecount_fID = fopen(Asitecount_fname, 'r');
% Init distribution cell and idx
Asitecount = {};
i = 1;
% Load first Asitecount line
line = fgetl(Asitecount_fID);
while ischar(line)
    cts_by_codon = str2num(line);
    num_cod = length(cts_by_codon);
    cts_by_codon = reshape(cts_by_codon, [num_cod,1]);
    Asitecount{i} = cts_by_codon;
    i = i + 1;
    line = fgetl(Asitecount_fID);
end
num_genes = i-1;
% Reshape Asitecount cell to be a column
Asitecount = reshape(Asitecount, [num_genes,1]);
% Close codon file
fclose(Asitecount_fID);

disp(min(cellfun('length', Asitecount)))
disp(max(cellfun('length', Asitecount)))

% Get gene indexes to save
keep_gene_idxs = [];

j = 1 
for i=1:length(LiuChxGenes)
    keep_gene_idx = find(cellfun(@(x) isequal(x, LiuChxGenes{i}), GeneName));
    %disp(keep_gene_idx);
    if length(keep_gene_idx) > 0
        keep_gene_idxs(j) = keep_gene_idx;
        j = j + 1;
    else
        disp(keep_gene_idx);
    end
end

disp(length(keep_gene_idxs));

% Keep only Liu CHX genes
GeneName = GeneName(keep_gene_idxs);
Asitecount = Asitecount(keep_gene_idxs);
asite_density = asite_density(keep_gene_idxs);
CDS = CDS(keep_gene_idxs);
distribution = distribution(keep_gene_idxs);

% Save variables to file
save('../../riboshape/data/data.liu_chx_genes.mat', 'GeneName', 'Asitecount', 'asite_density', 'CDS', 'codon', 'distribution');
