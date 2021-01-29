%% Gene Association Analysis
% load parcellation,Z-maps,gene
hdr_par = spm_vol('4mm_Glasser360.nii');
vol_par =spm_read_vols(hdr_par);

V1=spm_vol('g1_Z_T2.nii');
Y1=spm_read_vols(V1);

gradients = zeros(360,1);
overlap_rate = zeros(360,1);
for i=1:360
    gradients(i,1)=mean(Y1(vol_par==i));
    overlap_rate(i,1)=length(find(Y1(vol_par==i)))/length(Y1(vol_par==i));
end
load('100DS360scaledRobustSigmoidNSGRNAseqQC1LRcortex_ROI_NOdistCorrEuclidean.mat','parcelExpression')
load('100DS360scaledRobustSigmoidNSGRNAseqQC1LRcortex_ROI_NOdistCorrEuclidean.mat','probeInformation')

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
load surrogate_g1.mat; % using gen_surrogate_map_for g1z.py to genderate surrogate maps

hdr_mask = spm_vol('Reslice_group_mask.nii');
vol_mask = spm_read_vols(hdr_mask);
ind = find(vol_mask);
PCTVARrand = zeros(10000,10);
Rsq = zeros(10000,1);
parfor j = 1:10000
    disp(j);
    gradients_sur = zeros(360,1);
    Y1 = zeros(hdr_par.dim);
    
    Y1(ind) = surrogate_g1(j,:);
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
csvwrite('PLS_stats.csv',myStats);

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
save('PLS1_ROIscore.mat','PLS1_ROIscores_280');
csvwrite('PLS1_ROIscores.csv',XS(:,1));
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
save('PLS2_ROIscore.mat','PLS2_ROIscores_280');
csvwrite('PLS2_ROIscores.csv',XS(:,2));
PLS2_score=XS(:,2);

PLS1weights = zeros(10027,10000);
PLS2weights = zeros(10027,10000);

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
fid1 = fopen('PLS1_geneWeights.csv','w');
for i=1:length(genes)
  fprintf(fid1,'%s, %d, %f\n', PLS1{i},geneindex1(i), Z1(i));
end
fclose(fid1);

PLS2sw = std(PLS2weights');
temp2=PLS2w./PLS2sw';
[Z2,ind2]=sort(temp2,'descend');
PLS2=PLS2ids(ind2);
geneindex2=geneindex2(ind2);
fid2 = fopen('PLS2_geneWeights.csv','w');
for i=1:length(genes)
  fprintf(fid2,'%s, %d, %f\n', PLS2{i},geneindex2(i), Z2(i));
end
fclose(fid2);

%% Draw variance explanation
load('PLS_stats.csv');

py = plot(sort(PLS_stats(2,:)','descend'),'.-');
py.Color = [115 130 184]/255;
py.MarkerSize = 7;
hold on
plot(0:11,0.1*ones(1,12),'--','Color',[226,115,134]./255,'LineWidth',0.5);
hold off

xlabel('Component');
ylabel('Variance explained');
set(gca,'XLim',[0,11]);
set(gca,'YLim',[0,0.4],'YTick',0:0.1:0.4);
t1 = text(1,0.36,'***','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',7);
t2 = text(2,0.36,'**','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',7);
set(gca,'LineWidth',0.5);
set(gca,'FontName','Arial','FontSize',7);
box off

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 9 5.6]);
print(gcf,'Variance.tif','-dtiff','-r1000')

%% generate 360weighted image for BrainNet Viewer
gii1 = gifti('Glasser180_210P_L.label.gii');
gii2 = gifti('Glasser180_210P_R.label.gii');
ParcelLabel = double([gii1.cdata;gii2.cdata]);

load PLS1_ROIscore.mat
load region_ind.mat
PLS1_weight = zeros(length(ParcelLabel),1);
for i = 1:280
    PLS1_weight(ParcelLabel==region_ind(i)) = PLS1_ROIscores_280(i);
end
save('PLS1_weight.txt','PLS1_weight','-ascii');

load PLS2_ROIscore.mat
PLS2_weight = zeros(length(ParcelLabel),1);
for i = 1:280
    PLS2_weight(ParcelLabel==region_ind(i)) = PLS2_ROIscores_280(i);
end
save('PLS2_weight.txt','PLS2_weight','-ascii');

%% correlation between PLS score and G1 difference
load PLS1_ROIscore.mat
load PLS2_ROIscore.mat
load region_ind.mat
V=spm_vol('4mm_Glasser360.nii');
Y = spm_read_vols(V);

V1 = spm_vol('g1_Z_T2.nii');
Y1 = spm_read_vols(V1);

gradients = zeros(360,1);

for i=1:360
    gradients(i,1)=mean(Y1(Y==i));
end
corr_real = zeros(2,1);
corr_real(1) = corr(gradients(region_ind),PLS1_ROIscores_280);
corr_real(2) = corr(gradients(region_ind),PLS2_ROIscores_280);


load surrogate_g1.mat
hdr_mask = spm_vol('Reslice_group_mask.nii');
vol_mask = spm_read_vols(hdr_mask);
ind = find(vol_mask);
corr_surr = zeros(2,10000);
gradients_surr = zeros(360,1);
for j = 1:10000
    disp(j);
    Y_s = zeros(V.dim);
    Y_s(ind) = surrogate_g1(j,:);
    for i = 1:360
        gradients_surr(i) = mean(Y_s(Y==i));
    end
    corr_surr(1,j) = corr(gradients_surr(region_ind),PLS1_ROIscores_280);
    corr_surr(2,j) = corr(gradients_surr(region_ind),PLS2_ROIscores_280);
end
p = zeros(2,1);
p(1) = length(find(corr_surr(1,:)>corr_real(1)))/10000;
p(2) = length(find(corr_surr(2,:)>corr_real(1)))/10000;
save corr_surr.mat corr_surr corr_real p

%% Draw figures for correlations
dotcolor = [115 130 184]/255;
linecolor = [0 0 0];
ylable1 = 'PLS1 scores';
xlable1 = '{\itZ}-statistic of gradient 1';
[xData, yData] = prepareCurveData(gradients(region_ind,1), PLS1_ROIscores_360(region_ind));
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
set(gca,'XTick',-5:2.5:5);
t1 = text(-3,-0.17,'{\itr} = 0.551, {\itP} < 0.0001','FontName','Arial','FontSize',6);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[0 0 3.65 3.35]);
grid off
box off
print(gcf,'corr_PLS1_g1.tif','-dtiff','-r1000')

dotcolor = [115 130 184]/255;
linecolor = [0 0 0];
ylable1 = 'PLS1 scores';
xlable1 = '{\itZ}-statistic of gradient 2';
[xData, yData] = prepareCurveData(gradients(region_ind,2), PLS1_ROIscores_360(region_ind));
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
t1 = text(-3,-0.17,'{\itr} = 0.264, {\itP} = 0.016','FontName','Arial','FontSize',6);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[0 0 3.65 3.35]);
grid off
box off
print(gcf,'corr_PLS1_g2.tif','-dtiff','-r1000')