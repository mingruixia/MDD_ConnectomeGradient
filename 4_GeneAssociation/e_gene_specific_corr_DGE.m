%% process pipeline
%% define work dir
data_dir = 'D:\Data\DIDA-MDD\gradient_analysis\analysis2\';
script_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\3.Scripts\';
figure_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\5.Figures\ForRevision1\';

% load DGE for mdd, asd, scz,and bd
% DGE data was download from https://science.sciencemag.org/highwire/filestream/705756/field_highwire_adjunct_files/1/aad6469_Gandal_SM_Data-Table-S1.xlsx
dge_data_all = importdata([script_dir,'aad6469_Gandal_SM_Data-Table-S1_DGE.xlsx']);
gene_name_dge = dge_data_all.textdata(2:end,3);
dge_disorder = cell(6,3);
dge_disorder{1,1} = 'asd';
dge_disorder{1,2} = dge_data_all.data(:,5);
dge_disorder{2,1} = 'scz';
dge_disorder{2,2} = dge_data_all.data(:,8);
dge_disorder{3,1} = 'bd';
dge_disorder{3,2} = dge_data_all.data(:,11);
dge_disorder{4,1} = 'mdd';
dge_disorder{4,2} = dge_data_all.data(:,14);
dge_disorder{5,1} = 'add';
dge_disorder{5,2} = dge_data_all.data(:,17);
dge_disorder{6,1} = 'ibd';
dge_disorder{6,2} = dge_data_all.data(:,20);
for i = 1:6
    dge_disorder{i,3} = dge_data_all.data(:,i * 3 + 4);
end

% load pls weight
pls_weight_data = importdata([data_dir,'GeneAssociation_Main\PLS1_geneWeights.csv']);
gene_name_pls = pls_weight_data.textdata;
pls_weight = pls_weight_data.data(:,2);

% match gene
gene_match_id = zeros(length(gene_name_pls),1);
for i = 1:length(gene_name_pls)
    id_tmp = find(strcmp(gene_name_dge,gene_name_pls{i}));
    if ~isempty(id_tmp)
        gene_match_id(i) = id_tmp(1);
    end
end

gene_name_pls_match = gene_name_pls;
gene_name_pls_match(gene_match_id==0) = [];
pls_weight_match = pls_weight;
pls_weight_match(gene_match_id==0) = [];
gene_match_id(gene_match_id==0) = [];
gene_name_dge_match = gene_name_dge(gene_match_id);
dge_disorder_match = cell(6,3);
for i = 1:6
    dge_disorder_match{i,1} = dge_disorder{i,1};
    dge_disorder_match{i,2} = dge_disorder{i,2}(gene_match_id);
    dge_disorder_match{i,3} = dge_disorder{i,3}(gene_match_id);
end

% correlation
r = zeros(6,1);
p = zeros(6,1);

for i = 1:6
    ex_id = ~isnan(dge_disorder_match{i,2})&dge_disorder_match{i,3}<0.05;
%     ex_id = ~isnan(dge_disorder_match{i,2})&pls_weight_match<-3;
%     id_disorder = 1:6;
%     id_disorder(id_disorder == i) = [];
%     for j = 1:5
%         ex_id(dge_disorder_match{j,3}<0.05) =0;
%     end
        
    x = pls_weight_match(ex_id);
    y = dge_disorder_match{i,2}(ex_id);    
    r(i) = corr(x,y); 
    subplot(2,3,i);
    subtitle(dge_disorder_match{i,1})
    plot(y,x,'.');
    lsline;
    title(dge_disorder{i,1})
    r_surr = zeros(10000,1);
    for j = 1:10000
        perm_id = randperm(length(x));
        r_surr(j) = corr(x,y(perm_id));
    end
    if r(i) > 0
        p(i) = length(find(r_surr>(r(i))))/10000;
    else
        p(i) = length(find(r_surr<(r(i))))/10000;
    end
end





