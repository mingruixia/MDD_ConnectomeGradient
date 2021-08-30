%% process pipeline
%% define work dir
data_dir = 'D:\Data\DIDA-MDD\gradient_analysis\analysis2\';
script_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\3.Scripts\';
figure_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\5.Figures\ForRevision1\';


%% Generate connectivity matrix and calculate gradients
cd([data_dir,'Resliced']);
list = dir;
list(1:2) = [];

for i = 1:length(list)
% for i = [1,654,678]
    tic
    cd([data_dir,'Resliced\',list(i).name]);
    filename = dir('*.nii');
    M = x_gen_matrix_voxel([script_dir,'Reslice_group_mask.nii'],filename.name);
    M = M - diag(diag(M));
    M = 0.5*log((1+M)./(1-M));
    gm_k = GradientMaps('kernel','na','approach','dm');
    gm_k = gm_k.fit(M);
    ind_dir = [data_dir,'Gradient_zMat\',list(i).name];
    mkdir(ind_dir);
    save([ind_dir,'\gm_k.mat'],'gm_k');
    toc
end

%% arrange emb
cd([data_dir,'Gradient_zMat']);
list = dir;
list(1:2) = [];
n_sub = length(list);

emb_all = cell(n_sub,1);
lambda_all = cell(n_sub,1);
for i = 1:n_sub
    cd([data_dir,'Gradient_zMat\',list(i).name]);
    load('gm_k.mat');
    emb_all{i} = gm_k.gradients{1};
    lambda_all{i} = gm_k.lambda{1};    
end

% aligned across subjects
[realigned, xfms] = mica_iterativeAlignment(emb_all,100);
realigned = real(realigned);
xfms = cellfun(@real,xfms,'UniformOutput',false);

gradient_emb = cell(10,1);
for i = 1:10
    gradient_emb{i} = squeeze(realigned(:,i,:));
end

%calculate order sequence
seq = zeros(n_sub,10);
for i = 1:n_sub
    tmp = abs(xfms{i});
    [~,I] = sort(tmp,'descend');
    seq_tmp = zeros(5,1);
    seq_tmp(1) = I(1);
    for j = 2:10
        for k = 1:10
            if isempty(find(seq_tmp == I(k,j), 1))
                seq_tmp(j) = I(k,j);
                break;
            end
        end
    end
    seq(i,:) = seq_tmp;
end

%calculate explaination rate
exprate = zeros(n_sub,10);
for i = 1:n_sub
    tmp = lambda_all{i}./sum(lambda_all{i});
    exprate(i,:) = tmp(seq(i,:));
end


% reorder gradient according to explaination ratio
exprate_mean = mean(exprate);
[~,I] = sort(exprate_mean,'descend');
gradient_emb_reordered = cell(10,1);
for i = 1:10
    gradient_emb_reordered{i} = gradient_emb{I(i)};
end
exprate_reordered = exprate(:,I);

% harmonization
sub_info = xlsread([script_dir,'All_sub_info.xlsx']);
gradient_emb_correct = cell(3,1);
for i = 1:3
    gradient_emb_correct{i} = combat(gradient_emb_reordered{i},sub_info(:,2)',sub_info(:,[1 3 4]));
end

% visual check mean gradient
hdr_mask = spm_vol([script_dir,'Reslice_group_mask.nii']);
vol_mask = spm_read_vols(hdr_mask);
hdr = hdr_mask;
hdr.dt(1) = 64;
ind = find(vol_mask);
cd(data_dir);
for i = 1:3
    vol = zeros(hdr.dim);
    vol(ind) = mean(gradient_emb_correct{i},2);
    hdr.fname = ['mean_g_',num2str(i),'.nii'];    
    spm_write_vol(hdr,vol);
end

for i = 1:2
    gradient_emb_correct{i} = -gradient_emb_correct{i};
end
save([data_dir,'\Validation_zMat\gradient_emb_correct(1-3).mat'],'gradient_emb_correct');

% correct explaination
exprate_corrected = combat(exprate_reordered',sub_info(:,2)',sub_info(:,[1 3 4]));
exprate_corrected = exprate_corrected';
save([data_dir,'Validation_zMat\exprate_corrected.mat'], 'exprate_corrected');

% calculate gradient range
emb_range = zeros(n_sub,3);
for i = 1:3
    for j = 1:n_sub
        emb_range(j,i) = max(gradient_emb_correct{i}(:,j)) - min(gradient_emb_correct{i}(:,j));
    end
end
save([data_dir,'Validation_zMat\emb_range.mat'],'emb_range');

% calculate variation
emb_std = zeros(n_sub,3);
for i = 1:3
    for j = 1:n_sub
        emb_std(j,i) = std(gradient_emb_correct{i}(:,j));
    end
end
save([data_dir,'Validation_zMat\emb_std.mat'],'emb_std');

%% Between-group difference for global topology

% select subject
des = [sub_info(:,1),sub_info(:,3),sub_info(:,4)];
id_hc = find(sub_info(:,1)==1);
id_mdd = find(sub_info(:,1)==2);
n_mdd = length(id_mdd);
n_hc = length(id_hc);
disp([n_mdd,n_hc]);
% explained ratio
stat_exp = zeros(3,7);
for i = 1:3
    stat_exp(i,1) = mean(exprate_corrected(id_mdd,i));
    stat_exp(i,2) = std(exprate_corrected(id_mdd,i));
    stat_exp(i,3) = mean(exprate_corrected(id_hc,i));
    stat_exp(i,4) = std(exprate_corrected(id_hc,i));
    stat_result = regstats(exprate_corrected(:,i),des,'linear',{'tstat','r'});
    stat_exp(i,5) = stat_result.tstat.t(2);
    stat_exp(i,6) = stat_result.tstat.t(2) * sqrt(1/n_mdd + 1/n_hc);
    stat_exp(i,7) = stat_result.tstat.pval(2);
end
disp(stat_exp)

% gradient range
stat_range = zeros(3,7);
for i = 1:3
    stat_range(i,1) = mean(emb_range(id_mdd,i));
    stat_range(i,2) = std(emb_range(id_mdd,i));
    stat_range(i,3) = mean(emb_range(id_hc,i));
    stat_range(i,4) = std(emb_range(id_hc,i));
    stat_result = regstats(emb_range(:,i),des,'linear',{'tstat','r'});
    stat_range(i,5) = stat_result.tstat.t(2);
    stat_range(i,6) = stat_result.tstat.t(2) * sqrt(1/n_mdd + 1/n_hc);
    stat_range(i,7) = stat_result.tstat.pval(2);
end
disp(stat_range)

% gradient variance
stat_std = zeros(3,7);
for i = 1:3
    stat_std(i,1) = mean(emb_std(id_mdd,i));
    stat_std(i,2) = std(emb_std(id_mdd,i));
    stat_std(i,3) = mean(emb_std(id_hc,i));
    stat_std(i,4) = std(emb_std(id_hc,i));
    stat_result = regstats(emb_std(:,i),des,'linear',{'tstat','r'});
    stat_std(i,5) = stat_result.tstat.t(2);
    stat_std(i,6) = stat_result.tstat.t(2) * sqrt(1/n_mdd + 1/n_hc);
    stat_std(i,7) = stat_result.tstat.pval(2);
end
disp(stat_std)

save([data_dir,'Validation_zMat\Topography_stat.mat'],'stat_exp','stat_range','stat_exp');




%% Reginal difference
% estimate group difference


[n_voxel,n_sub] = size(gradient_emb_correct{1});
T_value = zeros(1,n_voxel);
res_map = cell(3,1);

for i = 1:3
    res_map{i} = zeros(n_voxel,n_sub);
    for j = 1:n_voxel
        stat_result = regstats(gradient_emb_correct{i}(j,:)',des,'linear',{'tstat','r'});
        T_value(i,j) = stat_result.tstat.t(2);
        res_map{i}(j,:) = stat_result.r;
    end
end
dfe = stat_result.tstat.dfe;
Z_value = spm_t2z(T_value,dfe);

% estimate GRF correction, save cluster nii, and draw pic
hdr_mask = spm_vol([script_dir,'Reslice_group_mask.nii']);
vol_mask = spm_read_vols(hdr_mask);
ind = find(vol_mask);
vox = [4 4 4];
voxel_p = 0.001;
cluster_p = 0.05;
tail = 2;
zthrd = norminv(1 - voxel_p/tail);
surfacefile = 'D:\Dropbox\Projects\BrainNet\otherfile\Surface\FSaverage_midthickness_32K.nv';
cd([data_dir,'Validation_zMat']);

for i = 1:3
    R_volume = zeros([hdr_mask.dim,n_sub]);
    for j = 1:n_sub
        tmp_vol = zeros(hdr_mask.dim);
        tmp_vol(ind) = res_map{i}(:,j);
        R_volume(:,:,:,j) = tmp_vol;
    end
    [cluster_size,dlh,fwhm] = x_GRF(R_volume,dfe,vol_mask,vox,voxel_p,cluster_p,tail);
    Z_vol = zeros(hdr_mask.dim);
    Z_vol(ind) = Z_value(i,:);
    
    
    hdr_z = hdr_mask;
    hdr_z.fname = ['g',num2str(i),'_T2_z.nii'];
    hdr_z.descrip = sprintf('SPM{Z_[1]}{dLh_%f}{FWHMx_%fFWHMy_%fFWHMz_%fmm}',...
        dlh, fwhm(1), fwhm(2), fwhm(3));
    hdr_z.dt(1) = 64;
    spm_write_vol(hdr_z,Z_vol);
    
    Z_vol(Z_vol < zthrd & Z_vol > -zthrd) = 0;
    [L, num] = bwlabeln(Z_vol,26);
    n = 0;
    for x = 1:num
        theCurrentCluster = L == x;
        if length(find(theCurrentCluster)) <= cluster_size
            n = n + 1;
            Z_vol(logical(theCurrentCluster)) = 0;
        end
    end
    hdr_z.fname = ['g',num2str(i),'_T2_z_cluster.nii'];
    spm_write_vol(hdr_z,Z_vol);
    niifile = hdr_z.fname;
    txtfile = ['g',num2str(i),'_T2_z_cluster.txt'];
    
    hdr_z.fname = ['g',num2str(i),'_T2_z_cluster_mask.nii'];
    Z_vol(Z_vol~=0) = 1;
    spm_write_vol(hdr_z,Z_vol);
    
    BrainNet_nii2txt(surfacefile,niifile,txtfile);
    BrainNet_MapCfg('D:\Dropbox\Projects\BrainNet\otherfile\Surface\FSaverage_inflated_32K.nv',...
        txtfile,[script_dir,'BNVCfg_T2_cluster.mat'],['g',num2str(i),'_T2_z_cluster.tif']);
end

for i = 1:3
    txt = Z_value(i,:);
    save(['g',num2str(i),'_T2_z.txt'],'txt','-ascii');
end

