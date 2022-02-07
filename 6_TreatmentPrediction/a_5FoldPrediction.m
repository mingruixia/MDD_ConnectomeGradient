data_dir = 'D:\Data\DIDA-MDD\gradient_analysis\analysis2\';
script_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\3.Scripts\';
figure_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\5.Figures\ForRevision1\';

%% SVR using g1 map
clinicinfo = load([data_dir, 'Treatment_Prediction\clinicinfo.txt']);
maskfile = [script_dir,'Reslice_group_mask.nii'];
hdr = spm_vol(maskfile);
vol = spm_read_vols(hdr);
ind = find(vol);

cd([data_dir,'Treatment_Prediction\gradient_nii\g1'])
list = dir('*.nii');
Subjects_Data = zeros(length(ind),20);
for i = 1:length(list)
    hdr_g1 = spm_vol(list(i).name);
    vol_g1 = spm_read_vols(hdr_g1);
    Subjects_Data(:,i) = vol_g1(ind);
end

Subjects_Scores = clinicinfo(:,6);
Covariates = [];
FoldQuantity = 5;
Pre_Method = 'Normalize';
C_Range = power(2, -5:10);
Weight_Flag = 1;
Permutation_Flag = 0;
ResultantFolder = [data_dir,'gradient_nii\g1\SVR'];
Prediction = SVR_NFolds_Sort_CSelect(Subjects_Data', Subjects_Scores', Covariates, FoldQuantity, Pre_Method, C_Range, Weight_Flag, Permutation_Flag, ResultantFolder);
score_predicted = cell2mat(Prediction.Score');
id = cell2mat(Prediction.Origin_ID);
corr_res = corr(score_predicted,clinicinfo(id,6));

%% permutation test
per_corr = zeros(1000,1);
Covariates = [];
FoldQuantity = 5;
Pre_Method = 'Normalize';
C_Range = power(2, -5:10);
Weight_Flag = 1;
Permutation_Flag = 1;
parfor i = 1:1000
    Prediction = SVR_NFolds_Sort_CSelect(Subjects_Data', Subjects_Scores', Covariates, FoldQuantity, Pre_Method, C_Range, Weight_Flag, Permutation_Flag);
    a = cell2mat(Prediction.Score');
    id = cell2mat(Prediction.Origin_ID);
    per_corr(i) = corr(a,clinicinfo(id,6));
end
p = length(find(per_corr>corr_res))/1000;

%% Draw figures for correlations
close all
dotcolor = [115 130 184]/255;
linecolor = [0 0 0];
ylable1 = 'Predicted HDRS change';
xlable1 = 'Observed HDRS change';
x = clinicinfo(id,6);
y = a;
[xData, yData] = prepareCurveData(x,y);
ft = fittype( 'poly1' );
opts = fitoptions( ft );
opts.Lower = [-Inf -Inf];
opts.Upper = [Inf Inf];
[fitresult, gof] = fit( xData, yData, ft, opts );
h=plot( fitresult, xData, yData);
set(h(1),'Marker','.','MarkerSize',6,'Color',dotcolor)
set(h(2),'LineWidth',0.5,'Color',linecolor)
hold on
xFit = linspace(min(xData),max(xData),100);
yPredict = predint(fitresult,xFit,0.95,'functional','off');
fy = cat(2,yPredict(:,2)',flip(yPredict(:,1),1)')';
fx = cat(2,xFit,flip(xFit',1)')';
fill(fx,fy,[0.5 0.5 0.5],'EdgeAlpha',0,'FaceAlpha',0.3);
hold off
legend off
ylabel(ylable1);
xlabel(xlable1);
set(gca,'LineWidth',0.5);
set(gca,'FontName','Arial','FontSize',6);
set(gca,'YLim',[16,28]);
set(gca,'YTick',16:3:28);
set(gca,'XLim',[12,32]);
set(gca,'XTick',12:5:32);
t1 = text(18,18,'{\itr} = 0.652, {\itP} < 0.001','FontName','Arial','FontSize',6);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[0 0 3.8 3.4]);
grid off
box off
print(gcf,[figure_dir,'corr_SVR.tif'],'-dtiff','-r1000')

%% draw weight figure
hdr = spm_vol(maskfile);
vol = spm_read_vols(hdr);
ind = find(vol);
vol_weight = zeros(hdr.dim);
load([data_dir,'gradient_nii\g1\SVR\w_Brain.mat']);
vol_weight(ind) = abs(w_Brain);
hdr_weight = hdr;
hdr_weight.fname = 'abs_brain_weight.nii';
hdr_weight.dt(1) = 64;
spm_write_vol(hdr_weight,vol_weight);
surfacefile = 'D:\Dropbox\Projects\BrainNet\otherfile\Surface\FSaverage_midthickness_32K.nv';
txtfile = 'abs_brain_weight.txt';
BrainNet_nii2txt(surfacefile,'abs_brain_weight.nii',txtfile);
BrainNet_MapCfg('D:\Dropbox\Projects\BrainNet\otherfile\Surface\FSaverage_inflated_32K.nv',...
    txtfile,[script_dir,'BNVCfg_SVR_weight.mat'],[figure_dir,'SVR_brain_weight.tif']);

%% weight distribution
vol_mod = spm_read_vols(spm_vol([script_dir,'mod_mask.nii']));
ind_mod = vol_mod(find(vol_mod));
weight_mod = zeros(8,1);
abs_w_Brain = abs(w_Brain);
for i = 1:8
    weight_mod(i) = sum(abs_w_Brain(ind_mod==i));
end
ratio_mod = weight_mod / sum(weight_mod);

%% SVR using FCS map
clinicinfo = load([data_dir, 'Treatment_Prediction\clinicinfo.txt']);
maskfile = [script_dir,'Reslice_group_mask.nii'];
hdr = spm_vol(maskfile);
vol = spm_read_vols(hdr);
ind = find(vol);

cd([data_dir,'Treatment_Prediction\CpLpFCS\weighted'])
list = dir('*.nii');
Subjects_Data = zeros(length(ind),20);
for i = 1:length(list)
    hdr_fcs = spm_vol(list(i).name);
    vol_fcs = spm_read_vols(hdr_fcs);
    Subjects_Data(:,i) = vol_fcs(ind);
end

Subjects_Scores = clinicinfo(:,6);
Covariates = [];
FoldQuantity = 5;
Pre_Method = 'Normalize';
C_Range = power(2, -5:10);
Weight_Flag = 1;
Permutation_Flag = 0;
ResultantFolder = [data_dir,'Treatment_Prediction\CpLpFCS\SVR'];
Prediction = SVR_NFolds_Sort_CSelect(Subjects_Data', Subjects_Scores', Covariates, FoldQuantity, Pre_Method, C_Range, Weight_Flag, Permutation_Flag, ResultantFolder);
score_predicted = cell2mat(Prediction.Score');
id = cell2mat(Prediction.Origin_ID);
corr_res = corr(score_predicted,clinicinfo(id,6));

%% permutation test
per_corr = zeros(1000,1);
Covariates = [];
FoldQuantity = 5;
Pre_Method = 'Normalize';
C_Range = power(2, -5:10);
Weight_Flag = 1;
Permutation_Flag = 1;
parfor i = 1:1000
    Prediction = SVR_NFolds_Sort_CSelect(Subjects_Data', Subjects_Scores', Covariates, FoldQuantity, Pre_Method, C_Range, Weight_Flag, Permutation_Flag);
    a = cell2mat(Prediction.Score');
    id = cell2mat(Prediction.Origin_ID);
    per_corr(i) = corr(a,clinicinfo(id,6));
end
p = length(find(per_corr>corr_res))/1000;

%% correlation between Cp/Lp and treatment outcomes
Cp_all = zeros(20,1);
cd([data_dir,'Treatment_Prediction\CpLpFCS\weighted'])
list = dir('*cp.txt');
for i = 1:length(list)
    Cp_all(i) = load(list(i).name);
end
[r,p] = corr(Cp_all,clinicinfo(:,6))

Lp_all = zeros(20,1);
cd([data_dir,'Treatment_Prediction\CpLpFCS\weighted'])
list = dir('*lp.txt');
for i = 1:length(list)
    Lp_all(i) = load(list(i).name);
end
[r,p] = corr(Lp_all,clinicinfo(:,6))

