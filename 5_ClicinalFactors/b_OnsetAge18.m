%% process pipeline
%% define work dir
data_dir = 'D:\Data\DIDA-MDD\gradient_analysis\analysis2\';
script_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\3.Scripts\';
figure_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\5.Figures\ForRevision1\';



%% Between-group difference for global topology
% select subject
id_OAa18 = find(sub_info(:,1)==2&sub_info(:,7)>18);
id_OAb18 = find(sub_info(:,1)==2&((sub_info(:,7)<=18&sub_info(:,7)~=-1)|sub_info(:,3)<=18));

sub_id = [id_OAa18;id_OAb18];

n_a18 = length(id_OAa18);
n_b18 = length(id_OAb18);

des = [[ones(n_a18,1);zeros(n_b18,1)],sub_info(sub_id,3),sub_info(sub_id,4)];

% explained ratio
stat_exp = zeros(3,7);
for i = 1:3
    stat_exp(i,1) = mean(exprate_corrected(id_OAa18,i));
    stat_exp(i,2) = std(exprate_corrected(id_OAa18,i));
    stat_exp(i,3) = mean(exprate_corrected(id_OAb18,i));
    stat_exp(i,4) = std(exprate_corrected(id_OAb18,i));    
    stat_result = regstats(exprate_corrected(sub_id,i),des,'linear',{'tstat','r'});
    stat_exp(i,5) = stat_result.tstat.t(2);
    stat_exp(i,6) = stat_result.tstat.t(2) * sqrt(1/n_a18 + 1/n_b18);
    stat_exp(i,7) = stat_result.tstat.pval(2);
end
disp(stat_exp)

% gradient range
stat_range = zeros(3,7);
for i = 1:3
    stat_range(i,1) = mean(emb_range(id_OAa18,i));
    stat_range(i,2) = std(emb_range(id_OAa18,i));
    stat_range(i,3) = mean(emb_range(id_OAb18,i));
    stat_range(i,4) = std(emb_range(id_OAb18,i));    
    stat_result = regstats(emb_range(sub_id,i),des,'linear',{'tstat','r'});
    stat_range(i,5) = stat_result.tstat.t(2);
    stat_range(i,6) = stat_result.tstat.t(2) * sqrt(1/n_a18 + 1/n_b18);
    stat_range(i,7) = stat_result.tstat.pval(2);
end
disp(stat_range)

% gradient variance
stat_std = zeros(3,7);
for i = 1:3
    stat_std(i,1) = mean(emb_std(id_OAa18,i));
    stat_std(i,2) = std(emb_std(id_OAa18,i));
    stat_std(i,3) = mean(emb_std(id_OAb18,i));
    stat_std(i,4) = std(emb_std(id_OAb18,i));    
    stat_result = regstats(emb_std(sub_id,i),des,'linear',{'tstat','r'});
    stat_std(i,5) = stat_result.tstat.t(2);
    stat_std(i,6) = stat_result.tstat.t(2) * sqrt(1/n_a18 + 1/n_b18);
    stat_std(i,7) = stat_result.tstat.pval(2);
end
disp(stat_std)

save([data_dir,'BetweenGroupDiff_OA18\Topography_stat.mat'],'stat_exp','stat_range','stat_exp');

%% Reginal difference
% estimate group difference
n_sub = length(sub_id);
[n_voxel,~] = size(gradient_emb_correct{1});
T_value = zeros(1,n_voxel);
res_map = cell(3,1);

for i = 1:3
    res_map{i} = zeros(n_voxel,n_sub);
    for j = 1:n_voxel
        stat_result = regstats(gradient_emb_correct{i}(j,sub_id)',des,'linear',{'tstat','r'});
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
cd([data_dir,'BetweenGroupDiff_OA18']);

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
