% load displacement in gradient space as features
cd D:\Data\DIDA-MDD\gradient_analysis\analysis\
mask_hdr = spm_vol('conj_max.nii');
mask_vol = spm_read_vols(mask_hdr);
ind = find(mask_vol);
cd D:\Data\DIDA-MDD\gradient_analysis\analysis\dis_map\all
list = textread('PatientList.txt');
dis_data = zeros(length(list),length(ind));
for i = 1:length(list)
    vol = spm_read_vols(spm_vol(list{i}));
    dis_data(i,:) = vol(ind);
end

% perform SVR analysis
cd D:\Data\DIDA-MDD\gradient_analysis\analysis\
Subjects_Data = dis_data;
Subjects_Scores = load('hdrscbt_mdd.txt')';
Covariates = load('cov_hdrs.txt');
FoldQuantity = 10;
Pre_Method = 'Normalize';
C_Range = power(2, -5:10);
Weight_Flag = 1;
Permutation_Flag = 0;
ResultantFolder = 'D:\Data\DIDA-MDD\gradient_analysis\analysis\SVR_hdrs_dis';

Prediction = SVR_NFolds_Sort_CSelect(Subjects_Data, Subjects_Scores, Covariates, FoldQuantity, Pre_Method, C_Range, Weight_Flag, Permutation_Flag, ResultantFolder);
