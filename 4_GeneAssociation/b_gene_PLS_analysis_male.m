%% Gene Association Analysis in male subjects
%% define work dir
data_dir = 'D:\Data\DIDA-MDD\gradient_analysis\analysis2\';
script_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\3.Scripts\';
figure_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\5.Figures\ForRevision1\';


%% load parcellation,Z-maps,gene
hdr_par = spm_vol([script_dir,'4mm_Glasser360.nii']);
vol_par =spm_read_vols(hdr_par);

V1=spm_vol([data_dir,'BetweenGroupDiff_male\g1_T2_z.nii']);
Y1=spm_read_vols(V1);

gradients = zeros(360,1);
overlap_rate = zeros(360,1);
for i=1:360
    gradients(i,1)=mean(Y1(vol_par==i));
    overlap_rate(i,1)=length(find(Y1(vol_par==i)))/length(Y1(vol_par==i));
end
load([script_dir,'100DS360scaledRobustSigmoidNSGRNAseqQC1LRcortex_ROI_NOdistCorrEuclidean_male_donors.mat'],'parcelExpression')
load([script_dir,'100DS360scaledRobustSigmoidNSGRNAseqQC1LRcortex_ROI_NOdistCorrEuclidean_male_donors.mat'],'probeInformation')

%% remove missing roi
temp1=find(overlap_rate<0.5);
temp2=find(isnan(parcelExpression(:,2)));
missingdata_regions=union(temp1,temp2);
region_ind=setdiff(parcelExpression(:,1),missingdata_regions);

group_express=parcelExpression(region_ind,2:end);
gene_name = probeInformation.GeneSymbol;

GENEdata=group_express;
MRIdata=gradients(region_ind);

%% PLS_calculation
Y = zscore(MRIdata);
dim = 10;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(GENEdata,Y,dim,'CV',dim);
temp=cumsum(100*PCTVAR(2,1:dim));
Rsquared = temp(dim);

%align PLS components with desired direction%
R1 = corr([XS(:,1),XS(:,2),XS(:,3)],MRIdata);
if R1(1,1)<0
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0
    XS(:,2)=-1*XS(:,2);
end
if R1(3,1)<0
    XS(:,3)=-1*XS(:,3);
end

% permutation test
load([data_dir,'surrogate_maps_g1_z_male\surrogate_maps_g1_z_resample_male.mat']); % using gen_surrogate_map_for g1z.py to genderate surrogate maps
hdr_mask = spm_vol([script_dir,'Reslice_group_mask.nii']);
vol_mask = spm_read_vols(hdr_mask);
ind = find(vol_mask);
PCTVARrand = zeros(10000,10);
Rsq = zeros(10000,1);
parfor j = 1:10000
    disp(j);
    gradients_sur = zeros(360,1);
    Y1 = zeros(hdr_par.dim);    
    Y1(ind) = surrogate_maps(j,:);
    for i=1:360
        gradients_sur(i,1)=mean(Y1(vol_par==i));
    end
    MRIdata_s = gradients_sur(region_ind);
    [XLr,YLr,XSr,YSr,BETAr,PCTVARr,MSEr,statsr]=plsregress(GENEdata,MRIdata_s,dim);
    PCTVARrand(j,:)=PCTVARr(2,:);
    temp=cumsum(100*PCTVARr(2,1:dim));
    Rsq(j) = temp(dim);    
end
p_single = zeros(1,10);
for l=1:dim
    p_single(l)=length(find(PCTVARrand(:,l)>=PCTVAR(2,l)))/10000;
end
p_cum = length(find(Rsq>=Rsquared))/10000;
myStats=[PCTVAR; p_single];
csvwrite([data_dir,'GeneAssociation_male_donor\PLS_stats.csv'],myStats);


%% calculate corrected weight
gene_name = probeInformation.GeneSymbol;
geneindex=1:size(GENEdata,2);
genes = gene_name;
bootnum=10000;
X=GENEdata;
Y=zscore(MRIdata);
dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

[R1,p1]=corr(XS(:,1),MRIdata);
if R1(1,1)<0
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
[PLS1w,x1] = sort(stats.W(:,1),'descend');
PLS1ids=genes(x1);
geneindex1=geneindex(x1);
PLS1_ROIscores_280=XS(:,1);
save([data_dir,'GeneAssociation_male_donor\PLS1_ROIscore.mat'],'PLS1_ROIscores_280');
csvwrite([data_dir,'GeneAssociation_male_donor\PLS1_ROIscores.csv'],XS(:,1));
PLS1_score=XS(:,1);

[R2,p2]=corr(XS(:,2),MRIdata);
if R2(1,1)<0
    stats.W(:,2)=-1*stats.W(:,2);
    XS(:,2)=-1*XS(:,2);
end
[PLS2w,x2] = sort(stats.W(:,2),'descend');
PLS2ids=genes(x2);
geneindex2=geneindex(x2);
PLS2_ROIscores_280=XS(:,2);
save([data_dir,'GeneAssociation_male_donor\PLS2_ROIscore.mat'],'PLS2_ROIscores_280');
csvwrite([data_dir,'GeneAssociation_male_donor\PLS2_ROIscores.csv'],XS(:,2));
PLS2_score=XS(:,2);

PLS1weights = zeros(length(genes),10000);
PLS2weights = zeros(length(genes),10000);

res = zeros(10000,length(region_ind));

parfor i=1:bootnum
    myresample = randsample(size(X,1),size(X,1),1);
    res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled regions
    Yr=Y(myresample,:); % define Y for resampled regions
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data

    temp=stats.W(:,1);%extract PLS1 weights
    newW=temp(x1); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS1weights(:,i) = newW;%store (ordered) weights from this bootstrap run
    
    temp=stats.W(:,2);%extract PLS2 weights
    newW=temp(x2); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS2w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS2weights(:,i) = newW; %store (ordered) weights from this bootstrap run    
end

PLS1sw = std(PLS1weights');
temp1=PLS1w./PLS1sw';
[Z1,ind1]=sort(temp1,'descend');
PLS1=PLS1ids(ind1);
geneindex1=geneindex1(ind1);
fid1 = fopen([data_dir,'GeneAssociation_male_donor\PLS1_geneWeights.csv'],'w');
for i=1:length(genes)
  fprintf(fid1,'%s, %d, %f\n', PLS1{i},geneindex1(i), Z1(i));
end
fclose(fid1);

PLS2sw = std(PLS2weights');
temp2=PLS2w./PLS2sw';
[Z2,ind2]=sort(temp2,'descend');
PLS2=PLS2ids(ind2);
geneindex2=geneindex2(ind2);
fid2 = fopen([data_dir,'GeneAssociation_male_donor\PLS2_geneWeights.csv'],'w');
for i=1:length(genes)
  fprintf(fid2,'%s, %d, %f\n', PLS2{i},geneindex2(i), Z2(i));
end
fclose(fid2);


