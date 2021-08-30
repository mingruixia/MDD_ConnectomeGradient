%% NeuroSynth
%% define work dir
data_dir = 'D:\Data\DIDA-MDD\gradient_analysis\analysis2\';
script_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\3.Scripts\';
figure_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\5.Figures\ForRevision1\';

%% generate between-group mask
cd([data_dir,'NeuroSynth']);
for i = 1:3
    hdr = spm_vol([data_dir,'BetweenGroupDiff\g',num2str(i),'_T2_z_cluster.nii']);
    vol = spm_read_vols(hdr);
    vol_pos = vol;
    vol_pos(vol<0) = 0;
    hdr_pos = hdr;
    hdr_pos.fname = ['g',num2str(i),'_t2_z_cluster_pos.nii'];
    spm_write_vol(hdr_pos,vol_pos);
    vol_neg = abs(vol);
    vol_neg(vol>0) = 0;
    hdr_neg = hdr;
    hdr_neg.fname = ['g',num2str(i),'_t2_z_cluster_neg.nii'];
    spm_write_vol(hdr_neg,vol_neg);
    x_reslice([script_dir,'temlate_for_neurosynth.nii'],hdr_pos.fname,1);
    x_reslice([script_dir,'temlate_for_neurosynth.nii'],hdr_neg.fname,1);
end

%% using neurosynth (python) to identify cognitive terms

%% generate surrogate Z-maps for corgnitive terms
hdr_mask = spm_vol([script_dir,'Reslice_group_mask.nii']);
vol_mask = spm_read_vols(hdr_mask);
ind = find(vol_mask);
nVoxels = length(ind);
real_g1_z_hdr = spm_vol([data_dir,'BetweenGroupDiff\g1_T2_z.nii']);

tok=regexp(real_g1_z_hdr.descrip, '\{dLh_(.*?)\}\{FWHMx_(.*?)FWHMy_(.*?)FWHMz_(.*?)mm\}',...
    'tokens');
if isempty(tok) || numel(tok)~=1
    return;
end
dLh=str2double(tok{1}{1});

VoxelPThreshold = 0.001;
ClusterPThreshold = 0.05/2;

zThrd=norminv(1 - VoxelPThreshold/2);
D=3;
Em = nVoxels * (2*pi)^(-(D+1)/2) * dLh * (zThrd*zThrd-1)^((D-1)/2) * exp(-zThrd*zThrd/2);
EN = nVoxels * (1-normcdf(zThrd)); 
Beta = ((gamma(D/2+1)*Em)/(EN)) ^ (2/D); 
pTemp=1;
ClusterSize=0;
while pTemp >= ClusterPThreshold
    ClusterSize=ClusterSize+1;
    pTemp = 1 - exp(-Em * exp(-Beta * ClusterSize^(2/D))); 
end

load([data_dir,'surrogate_maps_g1_z\surrogate_maps_g1_z_resample.mat']);

parfor j = 1:10000
    id = ['0000',num2str(j)];
    pos_flag = 0;
    neg_flag = 0;
    cd([data_dir,'NeuroSynth\surr_g1_pos']);
    if exist(['rsurr_g1_z_pos_',id(end-4:end),'.nii'],'file')
        pos_flag = 1;
    end
    cd([data_dir,'NeuroSynth\surr_g1_neg']);
    if exist(['rsurr_g1_z_neg_',id(end-4:end),'.nii'],'file')
        neg_flag = 1;
    end
    
    if pos_flag*neg_flag == 0
        vol_surr = zeros(hdr_mask.dim);
        vol_surr(ind) = surrogate_maps(j,:);
        vol_surr(vol_surr < zThrd & vol_surr > -zThrd) = 0;
        [L, num] = bwlabeln(vol_surr,26);
        n = 0;
        for x = 1:num
            theCurrentCluster = L == x;
            if length(find(theCurrentCluster)) <= ClusterSize
                n = n + 1;
                vol_surr(logical(theCurrentCluster)) = 0;
            end
        end
        if pos_flag == 0
            vol_pos = vol_surr;
            vol_pos(vol_surr<0) = 0;
            hdr_pos = real_g1_z_hdr;
            hdr_pos.fname = ['surr_g1_z_pos_',id(end-4:end),'.nii'];
            cd([data_dir,'NeuroSynth\surr_g1_pos']);
            spm_write_vol(hdr_pos,vol_pos);
            x_reslice([script_dir,'temlate_for_neurosynth.nii'],['surr_g1_z_pos_',id(end-4:end),'.nii'],1);
            delete(['surr_g1_z_pos_',id(end-4:end),'.nii']);
        end
        
        if neg_flag == 0
            vol_neg = -vol_surr;
            vol_neg(vol_surr>0) = 0;
            hdr_neg = real_g1_z_hdr;
            hdr_neg.fname = ['surr_g1_z_neg_',id(end-4:end),'.nii'];
            cd([data_dir,'NeuroSynth\surr_g1_neg']);
            spm_write_vol(hdr_neg,vol_neg);
            x_reslice([script_dir,'temlate_for_neurosynth.nii'],['surr_g1_z_neg_',id(end-4:end),'.nii'],1);
            delete(['surr_g1_z_neg_',id(end-4:end),'.nii']);
        end
    end
end

%g2
real_g2_z_hdr = spm_vol([data_dir,'BetweenGroupDiff\g2_T2_z.nii']);
tok=regexp(real_g2_z_hdr.descrip, '\{dLh_(.*?)\}\{FWHMx_(.*?)FWHMy_(.*?)FWHMz_(.*?)mm\}',...
    'tokens');
if isempty(tok) || numel(tok)~=1
    return;
end
dLh=str2double(tok{1}{1});
VoxelPThreshold = 0.001;
ClusterPThreshold = 0.05/2;
zThrd=norminv(1 - VoxelPThreshold/2);
D=3;
Em = nVoxels * (2*pi)^(-(D+1)/2) * dLh * (zThrd*zThrd-1)^((D-1)/2) * exp(-zThrd*zThrd/2);
EN = nVoxels * (1-normcdf(zThrd)); 
Beta = ((gamma(D/2+1)*Em)/(EN)) ^ (2/D); 
pTemp=1;
ClusterSize=0;
while pTemp >= ClusterPThreshold
    ClusterSize=ClusterSize+1;
    pTemp = 1 - exp(-Em * exp(-Beta * ClusterSize^(2/D))); 
end
load([data_dir,'surrogate_maps_g2_z\surrogate_maps_g2_z_resample.mat']);
parfor j = 1:10000
    id = ['0000',num2str(j)];
    pos_flag = 0;
    neg_flag = 0;
    cd([data_dir,'NeuroSynth\surr_g2_pos']);
    if exist(['rsurr_g2_z_pos_',id(end-4:end),'.nii'],'file')
        pos_flag = 1;
    end
    cd([data_dir,'NeuroSynth\surr_g2_neg']);
    if exist(['rsurr_g2_z_neg_',id(end-4:end),'.nii'],'file')
        neg_flag = 1;
    end
    if pos_flag*neg_flag == 0
        vol_surr = zeros(hdr_mask.dim);
        vol_surr(ind) = surrogate_maps(j,:);
        vol_surr(vol_surr < zThrd & vol_surr > -zThrd) = 0;
        [L, num] = bwlabeln(vol_surr,26);
        n = 0;
        for x = 1:num
            theCurrentCluster = L == x;
            if length(find(theCurrentCluster)) <= ClusterSize
                n = n + 1;
                vol_surr(logical(theCurrentCluster)) = 0;
            end
        end
        if pos_flag == 0
            vol_pos = vol_surr;
            vol_pos(vol_surr<0) = 0;
            hdr_pos = real_g2_z_hdr;
            hdr_pos.fname = ['surr_g2_z_pos_',id(end-4:end),'.nii'];
            cd([data_dir,'NeuroSynth\surr_g2_pos']);
            spm_write_vol(hdr_pos,vol_pos);
            x_reslice([script_dir,'temlate_for_neurosynth.nii'],['surr_g2_z_pos_',id(end-4:end),'.nii'],1);
            delete(['surr_g2_z_pos_',id(end-4:end),'.nii']);
        end
        if neg_flag == 0
            vol_neg = -vol_surr;
            vol_neg(vol_surr>0) = 0;
            hdr_neg = real_g2_z_hdr;
            hdr_neg.fname = ['surr_g2_z_neg_',id(end-4:end),'.nii'];
            cd([data_dir,'NeuroSynth\surr_g2_neg']);
            spm_write_vol(hdr_neg,vol_neg);
            x_reslice([script_dir,'temlate_for_neurosynth.nii'],['surr_g2_z_neg_',id(end-4:end),'.nii'],1);
            delete(['surr_g2_z_neg_',id(end-4:end),'.nii']);
        end
    end
end

%g3
real_g3_z_hdr = spm_vol([data_dir,'BetweenGroupDiff\g3_T2_z.nii']);
tok=regexp(real_g3_z_hdr.descrip, '\{dLh_(.*?)\}\{FWHMx_(.*?)FWHMy_(.*?)FWHMz_(.*?)mm\}',...
    'tokens');
if isempty(tok) || numel(tok)~=1
    return;
end
dLh=str2double(tok{1}{1});
VoxelPThreshold = 0.001;
ClusterPThreshold = 0.05/2;
zThrd=norminv(1 - VoxelPThreshold/2);
D=3;
Em = nVoxels * (2*pi)^(-(D+1)/2) * dLh * (zThrd*zThrd-1)^((D-1)/2) * exp(-zThrd*zThrd/2);
EN = nVoxels * (1-normcdf(zThrd)); 
Beta = ((gamma(D/2+1)*Em)/(EN)) ^ (2/D);
pTemp=1;
ClusterSize=0;
while pTemp >= ClusterPThreshold
    ClusterSize=ClusterSize+1;
    pTemp = 1 - exp(-Em * exp(-Beta * ClusterSize^(2/D)));
end
load([data_dir,'surrogate_maps_g3_z\surrogate_maps_g3_z_resample.mat']);
parfor j = 1:10000
    id = ['0000',num2str(j)];
    pos_flag = 0;
    neg_flag = 0;
    cd([data_dir,'NeuroSynth\surr_g3_pos']);
    if exist(['rsurr_g3_z_pos_',id(end-4:end),'.nii'],'file')
        pos_flag = 1;
    end
    cd([data_dir,'NeuroSynth\surr_g3_neg']);
    if exist(['rsurr_g3_z_neg_',id(end-4:end),'.nii'],'file')
        neg_flag = 1;
    end
    if pos_flag*neg_flag == 0
        vol_surr = zeros(hdr_mask.dim);
        vol_surr(ind) = surrogate_maps(j,:);
        vol_surr(vol_surr < zThrd & vol_surr > -zThrd) = 0;
        [L, num] = bwlabeln(vol_surr,26);
        n = 0;
        for x = 1:num
            theCurrentCluster = L == x;
            if length(find(theCurrentCluster)) <= ClusterSize
                n = n + 1;
                vol_surr(logical(theCurrentCluster)) = 0;
            end
        end
        if pos_flag == 0
            vol_pos = vol_surr;
            vol_pos(vol_surr<0) = 0;
            hdr_pos = real_g3_z_hdr;
            hdr_pos.fname = ['surr_g3_z_pos_',id(end-4:end),'.nii'];
            cd([data_dir,'NeuroSynth\surr_g3_pos']);
            spm_write_vol(hdr_pos,vol_pos);
            x_reslice([script_dir,'temlate_for_neurosynth.nii'],['surr_g3_z_pos_',id(end-4:end),'.nii'],1);
            delete(['surr_g3_z_pos_',id(end-4:end),'.nii']);
        end
        if neg_flag == 0
            vol_neg = -vol_surr;
            vol_neg(vol_surr>0) = 0;
            hdr_neg = real_g3_z_hdr;
            hdr_neg.fname = ['surr_g3_z_neg_',id(end-4:end),'.nii'];
            cd([data_dir,'NeuroSynth\surr_g3_neg']);
            spm_write_vol(hdr_neg,vol_neg);
            x_reslice([script_dir,'temlate_for_neurosynth.nii'],['surr_g3_z_neg_',id(end-4:end),'.nii'],1);
            delete(['surr_g3_z_neg_',id(end-4:end),'.nii']);
        end
    end
end

% determine p
for j = 1:3
    fid = fopen([data_dir,'NeuroSynth\real_r_g',num2str(j),'_pos.txt']);
    c = textscan(fid, '%s %s','delimiter',',');
    fclose(fid);
    real_r = str2double(c{1,2}(2:end));
    surr_r = zeros(30,10000);
    parfor i = 1:10000
        id = ['0000',num2str(i)];
        fid = fopen([data_dir,'NeuroSynth\surr_g',num2str(j),'_pos\rsurr_g',num2str(j),'_z_pos_',id(end-4:end),'.txt']);
        c = textscan(fid, '%s %s','delimiter',',');
        fclose(fid);
        surr_r(:,i) = str2double(c{1,2}(2:end));
    end
    p_cogn_term = sum(gt(surr_r,real_r),2)/10000;
    [~,~,rnk] = unique(p_cogn_term);
    q_fdr = (p_cogn_term * length(p_cogn_term)) ./ rnk;
    save([data_dir,'NeuroSynth\cog_term_g',num2str(j),'_pos.mat'],'real_r','surr_r','p_cogn_term','q_fdr');
    
    fid = fopen([data_dir,'NeuroSynth\real_r_g',num2str(j),'_neg.txt']);
    c = textscan(fid, '%s %s','delimiter',',');
    fclose(fid);
    real_r = str2double(c{1,2}(2:end));
    surr_r = zeros(30,10000);
    parfor i = 1:10000
        id = ['0000',num2str(i)];
        fid = fopen([data_dir,'NeuroSynth\surr_g',num2str(j),'_neg\rsurr_g',num2str(j),'_z_neg_',id(end-4:end),'.txt']);
        c = textscan(fid, '%s %s','delimiter',',');
        fclose(fid);
        surr_r(:,i) = str2double(c{1,2}(2:end));
    end
    p_cogn_term = sum(gt(surr_r,real_r),2)/10000;
    [~,~,rnk] = unique(p_cogn_term);
    q_fdr = (p_cogn_term * length(p_cogn_term)) ./ rnk;
    save([data_dir,'NeuroSynth\cog_term_g',num2str(j),'_neg.mat'],'real_r','surr_r','p_cogn_term','q_fdr');
end