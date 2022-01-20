%% process pipeline
%% define work dir
data_dir = 'D:\Data\DIDA-MDD\gradient_analysis\analysis2\';
script_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\3.Scripts\';
figure_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\5.Figures\ForRevision1\';

%% Calculate voxel-wise JaccardInd

main_finding_clusters = cell(3,1);

main_finding_clusters{1} = spm_read_vols(spm_vol([data_dir,'BetweenGroupDiff\g1_T2_z_cluster_mask.nii']));
main_finding_clusters{2} = spm_read_vols(spm_vol([data_dir,'BetweenGroupDiff\g2_T2_z_cluster_mask.nii']));
main_finding_clusters{3} = spm_read_vols(spm_vol([data_dir,'BetweenGroupDiff\g3_T2_z_cluster_mask.nii']));

jarccard_all = zeros(15,1);

% site
for i = 1:10
    num = ['0',num2str(i)];
    validation_dir = ['Validation_LeaveOneSiteOut\site',num(end-1:end),'\'];
    validation_finding_clusters = cell(3,1);
    for j = 1:3
        validation_finding_clusters{j} = spm_read_vols(spm_vol([data_dir,validation_dir,'g',num2str(j),'_T2_z_cluster_mask.nii']));
    end
    pooled_intersection = 0;
    pooled_union = 0;
    
    for j = 1:3
        vol_intersection = main_finding_clusters{j} .* validation_finding_clusters{j};
        vol_union = sign(main_finding_clusters{j} + validation_finding_clusters{j});
        
        
        pooled_intersection = pooled_intersection + sum(vol_intersection(:));
        pooled_union = pooled_union + sum(vol_union(:));
    end
    
    jarccard_all(i) = pooled_intersection / pooled_union;
end

% other
vali_name = {'adult','mFDcov','jointemb','zMat','mFD025'};
for i = 1:5
    validation_dir = ['Validation_',vali_name{i},'\'];
    validation_finding_clusters = cell(3,1);
    for j = 1:3
        validation_finding_clusters{j} = spm_read_vols(spm_vol([data_dir,validation_dir,'g',num2str(j),'_T2_z_cluster_mask.nii']));
    end
    pooled_intersection = 0;
    pooled_union = 0;    
    for j = 1:3
        vol_intersection = main_finding_clusters{j} .* validation_finding_clusters{j};
        vol_union = sign(main_finding_clusters{j} + validation_finding_clusters{j});
        pooled_intersection = pooled_intersection + sum(vol_intersection(:));
        pooled_union = pooled_union + sum(vol_union(:));
    end    
    jarccard_all(i+10) = pooled_intersection / pooled_union;
end
