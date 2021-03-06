%% Downsample preprocessed images
cd D:\Data\DIDA-MDD\gradient_analysis\analysis\FunImgARWSDCFB
list = dir;
list(1:2) = [];
parfor i = 1:length(list)    
    cd(['D:\Data\DIDA-MDD\gradient_analysis\analysis\FunImgARWSDCFB\',list(i).name]);
    filename = dir('*.nii');    
    x_reslice('D:\Data\DIDA-MDD\gradient_analysis\analysis\Reslice_group_mask.nii',filename.name,4);
	mkdir(['D:\Data\DIDA-MDD\gradient_analysis\analysis\Resliced3\',list(i).name]);
    movefile('r*.nii',['D:\Data\DIDA-MDD\gradient_analysis\analysis\Resliced3\',list(i).name]);   
end


%% Generate connectivity matrix and calculate gradients
emb_all = cell(length(list),1);
res_all = cell(length(list),1);
for i = 1:length(list)
    cd(['D:\Data\DIDA-MDD\gradient_analysis\analysis\Resliced\',list(i).name]);
    filename = dir('*.nii');
    M = x_gen_matrix_voxel('D:\Data\DIDA-MDD\gradient_analysis\analysis\Reslice_group_mask.nii',filename.name);
    n = length(M);
    M_spar = M;
    tmp = sort(M);
    tmp = M - repmat(tmp(round(n*0.9),:),n,1);
    M_spar(tmp<0) = 0;
    M_spar = M_spar';    
    M_cos = 1 - squareform(pdist(M_spar,'cosine'));
    M_normalized = 1 - acos(M_cos)/pi;
    [emb,res] = x_compute_diffusion_map(M_normalized,0.5);
	filename = ['D:\Data\DIDA-MDD\gradient_analysis\analysis\Results\',list(i).name,'.mat'];
    x_savefile(filename,emb,res);
	emb_all{i} = emb;
	res_all{i} = res;
end


%% Align gradient maps across subjects
[realigned, xfms] = mica_iterativeAlignment(emb_all,100);
realigned = real(realigned);
xfms = cellfun(@real,xfms,'UniformOutput',false);
gradient_emb = cell(3,1);
for i = 1:3
    gradient_emb{i} = squeeze(realigned(:,i,:));
end


%% Determine adjusted gradient sequence
seq = zeros(n_sub,30);
for i = 1:n_sub
    tmp = abs(xfms{i});
    [~,I] = sort(tmp,'descend');
    seq_tmp = zeros(5,1);
    seq_tmp(1) = I(1);
    for j = 2:30
        for k = 1:30
            if isempty(find(seq_tmp == I(k,j), 1))
                seq_tmp(j) = I(k,j);
                break;
            end
        end
    end
    seq(i,:) = seq_tmp;
end


%% Calculate explained variance
exprate = zeros(n_sub,3);
for i = 1:n_sub
    tmp = res_all{i}.lambdas./sum(res_all{i}.lambdas);
    exprate(i,:) = tmp(seq(i,:));
end


%% Reorder gradient according to explained variance
exprate_mean = mean(exprate);
[~,I] = sort(exprate_mean,'descend');
gradient_emb_reordered = cell(3,1);
for i = 1:3
    gradient_emb_reordered{i} = gradient_emb{I(i)};
end
exprate_reordered = exprate(:,I);


%% Correct for center effect
variables = load('SiteGroupAgeSex.txt');
gradient_emb_correct = cell(3,1);
for i = 1:3
    gradient_emb_correct{i} = combat(gradient_emb_reordered{i},variables(:,1)',variables(:,2:4));
end


%% Calculate gradient range
n_sub = size(gradient_emb_correct{1},2);
emb_range = zeros(n_sub,3);
for i = 1:3
    for j = 1:n_sub
        emb_range(j,i) = max(gradient_emb_correct{i}(:,j)) - min(gradient_emb_correct{i}(:,j));
    end
end


%% Calculate gradient spatial variance
emb_std = zeros(n_sub,3);
for i = 1:3
    for j = 1:n_sub
        emb_std(j,i) = std(gradient_emb_correct{i}(:,j));
    end
end


%% Generate nifti file for gradient maps
hdr_mask = spm_vol('D:\Data\DIDA-MDD\gradient_analysis\code_test\Reslice_group_mask.nii');
vol_mask = spm_read_vols(hdr_mask);
hdr = hdr_mask;
hdr.dt(1) = 64;
ind = find(vol_mask);
for i = 1:3    
    mkdir(['D:\Data\DIDA-MDD\gradient_analysis\analysis\Corrected_emb\grad0',num2str(i)]);
    cd(['D:\Data\DIDA-MDD\gradient_analysis\analysis\Corrected_emb\grad0',num2str(i)]);
    for j = 1:length(list)        
        vol = zeros(hdr.dim);
        vol(ind) = gradient_emb_correct{i}(:,j);
        hdr.fname = ['g',num2str(i),'_',list(j).name,'.nii'];
        spm_write_vol(hdr,vol);
    end
end

%% Generate group mean gradient maps
for i = 1:3
    x_mean_image('['D:\Data\DIDA-MDD\gradient_analysis\analysis\Corrected_emb\grad0',num2str(i),'\HC'],['mean_hc_g',num2str(i),'.nii']);
    x_mean_image('['D:\Data\DIDA-MDD\gradient_analysis\analysis\Corrected_emb\grad0',num2str(i),'\MDD'],['mean_mdd_g',num2str(i),'.nii']);
end

%% Examine spatial correlation between the group-averaged maps of hc and mdd
% run surrogate_spatial.py to generate surrogate maps for hc
hdr_mask = spm_vol('D:\Data\DIDA-MDD\gradient_analysis\analysis\Reslice_group_mask.nii');
vol_mask = spm_read_vols(hdr_mask);
ind = find(vol_mask);
real_phi = zeros(3,1);
p = zeros(3,1);
for i = 1:3
    hdr_hc = spm_vol(['mean_hc_g',num2str(i),'.nii']);
    vol_hc = spm_read_vols(hdr_hc);
    hdr_mdd = spm_vol(['mean_mdd_g',num2str(i),'.nii']);
    vol_mdd = spm_read_vols(hdr_mdd);
    g_hc = vol_hc(ind);
    g_mdd = vol_mdd(ind);
    real_phi(i) = corr(g_hc,g_mdd,'type','spearman');
    load(['surrogate_maps_hc_g',num2str(i),'.mat']);
    surr_phi = corr(surrogate_maps',g_mdd,'type','spearman');
    p(i) = length(find(surr_phi>real_phi(i)))/10000;
end



