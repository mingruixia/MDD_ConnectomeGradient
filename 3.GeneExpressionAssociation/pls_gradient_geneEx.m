clear,clc

V=spm_vol('4mm_Glasser360.nii');
[Y,XYZmm]=spm_read_vols(V);

V1=spm_vol('g1_Z_T2.nii');
[Y1,XYZmm]=spm_read_vols(V1);

V2=spm_vol('g2_Z_T2.nii');
[Y2,XYZmm]=spm_read_vols(V2);

V3=spm_vol('g3_Z_T2.nii');
[Y3,XYZmm]=spm_read_vols(V3);

for i=1:360
    gradients(i,1)=mean(Y1(Y==i));
    gradients(i,2)=-mean(Y2(Y==i));
    gradients(i,3)=mean(Y3(Y==i));
    overlap_rate(i,1)=length(find(Y1(Y==i)))/length(Y1(Y==i));
end
load('100DS360scaledRobustSigmoidNSGRNAseqQC1LRcortex_ROI_NOdistCorrEuclidean.mat','parcelExpression')
load('100DS360scaledRobustSigmoidNSGRNAseqQC1LRcortex_ROI_NOdistCorrEuclidean.mat','probeInformation')

temp1=find(overlap_rate<0.5);
temp2=find(isnan(parcelExpression(:,2)));
missingdata_regions=union(temp1,temp2);
region_ind=setdiff(parcelExpression(:,1),missingdata_regions);
clear Y Y1 Y2 Y3 V V1 V2 V3 XYZmm temp1 temp2

group_express=parcelExpression(region_ind,2:end);
gene_name = probeInformation.GeneSymbol;

GENEdata=group_express;
MRIdata=gradients(region_ind,:);
PLS_calculate_stats_Jin(MRIdata, GENEdata, 'D:\MX\SynologyDrive\Projects\2019_MDD_gradient\4.Results\3.GeneExpression\PLS_new');
[PLS1_score] = PLS_bootstrap_Jin(MRIdata, GENEdata, gene_name, region_ind,'D:\MX\SynologyDrive\Projects\2019_MDD_gradient\4.Results\3.GeneExpression\PLS_new');

response_var_file= gradients(region_ind,:);
predictor_var_file = group_express;