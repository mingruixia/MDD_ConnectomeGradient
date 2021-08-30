%% process pipeline
%% define work dir
data_dir = 'D:\Data\DIDA-MDD\gradient_analysis\analysis2\';
script_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\3.Scripts\';
figure_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\5.Figures\ForRevision1\';

%% check if specific mdd gene had lower z
clear ish_data_all
% load PLS weight
pls_weight_data = importdata([data_dir,'GeneAssociation_Main\PLS1_geneWeights.csv']);
% abs_pls_weight = abs(pls_weight_data.data(:,2));
pls_weight = pls_weight_data.data(:,2);

% load gene list
ish_data_all(1).name = 'MDD';
ish_data_all(2).name = 'AD';
ish_data_all(3).name = 'ASD';
ish_data_all(4).name = 'EP';
ish_data_all(5).name = 'ID';
ish_data_all(6).name = 'PD';
ish_data_all(7).name = 'SCZ';

% find gene id
for i = 1:length(ish_data_all)
    fid = fopen([script_dir,'HBA_ISH_genelist_',ish_data_all(i).name,'.txt']);
    ish_data_all(i).genelist = textscan(fid,'%s');
    fclose(fid);
    ish_data_all(i).gene_id = zeros(length(ish_data_all(i).genelist{1,1}),1);
    for j = 1:length(ish_data_all(i).gene_id)
        id_tmp = find(strcmp(ish_data_all(i).genelist{1,1}{j},pls_weight_data.textdata));
        if ~isempty(id_tmp)
            ish_data_all(i).gene_id(j) = id_tmp(1);
        end
    end
    ish_data_all(i).gene_id_plsfilter = ish_data_all(i).gene_id;
    ish_data_all(i).gene_id_plsfilter(ish_data_all(i).gene_id_plsfilter==0) = [];
end

% calculate z value
for i = 1:length(ish_data_all)
    ish_data_all(i).pls_z = pls_weight(ish_data_all(i).gene_id_plsfilter);
    surr_z = zeros(10000,1);
    for j = 1:10000
         surr_weight = pls_weight(randperm(length(pls_weight)));
         surr_z(j) = mean(surr_weight(ish_data_all(i).gene_id_plsfilter));
    end
    ish_data_all(i).pls_z_p = length(find(surr_z<mean(ish_data_all(i).pls_z)))/10000;
end

% between-disorder difference in z value
for i = 1:length(ish_data_all)
    real_diff = mean(pls_weight(ish_data_all(1).gene_id_plsfilter)) - mean(pls_weight(ish_data_all(i).gene_id_plsfilter));
    n_gene_mdd = length(ish_data_all(1).gene_id_plsfilter);
    n_gene_other = length(ish_data_all(i).gene_id_plsfilter);
    value_all = pls_weight([ish_data_all(1).gene_id_plsfilter;ish_data_all(i).gene_id_plsfilter]);
    id_all = [ones(n_gene_mdd,1);ones(n_gene_other,1)*2];
    surr_diff = zeros(10000,1);
    for j = 1:10000
        surr_id = id_all(randperm(length(id_all)));
        surr_diff(j) = mean(value_all(surr_id==1)) - mean(value_all(surr_id==2));
    end
    ish_data_all(i).pls_z_compmdd_p = length(find(surr_diff<real_diff))/10000;
end

%% compare correlation r
% load parcellation,Z-maps,gene
hdr_par = spm_vol([script_dir,'4mm_Glasser360.nii']);
vol_par =spm_read_vols(hdr_par);

V1=spm_vol([data_dir,'BetweenGroupDiff\g1_T2_z.nii']);
Y1=spm_read_vols(V1);

gradients = zeros(360,1);
overlap_rate = zeros(360,1);
for i=1:360
    gradients(i,1)=mean(Y1(vol_par==i));
    overlap_rate(i,1)=length(find(Y1(vol_par==i)))/length(Y1(vol_par==i));
end
load([script_dir,'100DS360scaledRobustSigmoidNSGRNAseqQC1LRcortex_ROI_NOdistCorrEuclidean.mat'],'parcelExpression')
load([script_dir,'100DS360scaledRobustSigmoidNSGRNAseqQC1LRcortex_ROI_NOdistCorrEuclidean.mat'],'probeInformation')

% remove missing roi
temp1=find(overlap_rate<0.5);
temp2=find(isnan(parcelExpression(:,2)));
missingdata_regions=union(temp1,temp2);
region_ind=setdiff(parcelExpression(:,1),missingdata_regions);

group_express=parcelExpression(region_ind,2:end);
gene_name = probeInformation.GeneSymbol;

GENEdata=group_express;
MRIdata=gradients(region_ind);

% corr
load([data_dir,'surrogate_maps_g1_z\surrogate_maps_g1_z_resample.mat']); % using gen_surrogate_map_for g1z.py to genderate surrogate maps
hdr_mask = spm_vol([script_dir,'Reslice_group_mask.nii']);
vol_mask = spm_read_vols(hdr_mask);
ind = find(vol_mask);
for i = 1:length(ish_data_all)
    ish_data_all(i).g1t2z_corr_r = corr(GENEdata(:,ish_data_all(i).gene_id_plsfilter),MRIdata);
    corr_surr = zeros(length(ish_data_all(i).gene_id_plsfilter),10000);
    for j = 1:10000
        gradients_surr = zeros(360,1);
        Y1 = zeros(hdr_par.dim);
        Y1(ind) = surrogate_maps(j,:);
        for k=1:360
            gradients_surr(k,1)=mean(Y1(vol_par==k));
        end
        MRIdata_s = gradients_surr(region_ind);
        corr_surr(:,j) = corr(GENEdata(:,ish_data_all(i).gene_id_plsfilter),gradients_surr(region_ind));
    end
    ish_data_all(i).g1t2z_corr_p = zeros(length(ish_data_all(i).g1t2z_corr_r),1);
    for j = 1:length(ish_data_all(i).g1t2z_corr_r)
        if ish_data_all(i).g1t2z_corr_r(j)>0
            ish_data_all(i).g1t2z_corr_p(j) = length(find(corr_surr(j,:)>ish_data_all(i).g1t2z_corr_r(j)))/10000;
        else
            ish_data_all(i).g1t2z_corr_p(j) = length(find(corr_surr(j,:)<ish_data_all(i).g1t2z_corr_r(j)))/10000;
        end
    end
    [~,~,rnk] = unique(ish_data_all(i).g1t2z_corr_p);
    ish_data_all(i).g1t2z_corr_q_fdr = (ish_data_all(i).g1t2z_corr_p * length(ish_data_all(i).g1t2z_corr_p)) ./ rnk;
end

% between-disorder difference in r value
for i = 1:length(ish_data_all)
    real_diff = mean(abs(ish_data_all(1).g1t2z_corr_r)) - mean(abs(ish_data_all(i).g1t2z_corr_r));
    n_gene_mdd = length(ish_data_all(1).g1t2z_corr_r);
    n_gene_other = length(ish_data_all(i).g1t2z_corr_r);
    value_all = abs([ish_data_all(1).g1t2z_corr_r;ish_data_all(i).g1t2z_corr_r]);
    id_all = [ones(n_gene_mdd,1);ones(n_gene_other,1)*2];
    surr_diff = zeros(10000,1);
    for j = 1:10000
        surr_id = id_all(randperm(length(id_all)));
        surr_diff(j) = mean(value_all(surr_id==1)) - mean(value_all(surr_id==2));
    end
    ish_data_all(i).g1t2z_corr_r_compmdd_p = length(find(surr_diff>real_diff))/10000;
end



%% percentage of significant corr
for i = 1:length(ish_data_all)
    ish_data_all(i).g1t2z_corr_signperc = length(find(ish_data_all(i).g1t2z_corr_q_fdr<0.05))/length(ish_data_all(i).g1t2z_corr_q_fdr);
end

% between-disorder difference in percentage
for j = 1:length(ish_data_all)
    p_all = [ish_data_all(1).g1t2z_corr_p;ish_data_all(j).g1t2z_corr_p];
    n_gene_mdd = length(ish_data_all(1).g1t2z_corr_r);
    n_gene_other = length(ish_data_all(j).g1t2z_corr_r);
    id_all = [ones(n_gene_mdd,1);ones(n_gene_other,1)*2];
    diff_per = ish_data_all(1).g1t2z_corr_signperc - ish_data_all(j).g1t2z_corr_signperc;
    surr_diff_per = zeros(10000,1);
    for i = 10000
        surr_id = id_all(randperm(length(id_all)));
        surr_p_mdd = p_all(surr_id==1);
        [~,~,rnk] = unique(surr_p_mdd);
        surr_q_mdd = (surr_p_mdd * length(surr_p_mdd)) ./ rnk;
        surr_perc_mdd = length(find(surr_q_mdd<0.05))/length(surr_q_mdd);
        
        surr_p_comp = p_all(surr_id==2);
        [~,~,rnk] = unique(surr_p_comp);
        surr_q_comp = (surr_p_comp * length(surr_p_comp)) ./ rnk;
        surr_perc_comp= length(find(surr_q_comp<0.05))/52;
        
        surr_diff_per(i) = surr_perc_mdd - surr_perc_comp;
    end
    ish_data_all(j).g1t2z_corr_signperc_compmdd_p = length(find(surr_diff_per>diff_per_mddvsscz))/10000;
end

save([data_dir,'ish_data.mat'],'ish_data_all');

%% Draw figure for z
mean_z = zeros(7,1);
std_z = zeros(7,1);
se_z = zeros(7,1);
for i = 1:7
    mean_z(i) = mean(ish_data_all(i).pls_z);
    std_z(i) = std(ish_data_all(i).pls_z);
    se_z(i) = std_z(i)/length(ish_data_all(i).pls_z);
end
close all
bar(mean_z,'FaceColor',[226,115,134]./255);
hold on
errorbar(mean_z,se_z,'_k');
hold off
set(gca,'YLim',[-3.5,1.5],'YTick',[-3.5,0,1.5]);

text(2,1.2,'**','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',8);
for i = 4:7
    text(i,1.2,'*','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',8);
end

set(gca,'XTickLabel',{ish_data_all(1).name,ish_data_all(2).name,...
    ish_data_all(3).name,ish_data_all(4).name,ish_data_all(5).name,...
    ish_data_all(6).name,ish_data_all(7).name});
ylabel('PLS1 weight');
box off
set(gca,'XTickLabelRotation',60);
set(gca,'LineWidth',0.5);
set(gca,'FontName','Arial','FontSize',6);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 5 2.5]);
print(gcf,[figure_dir,'ish_z.tif'],'-dtiff','-r1000')

%% Draw pie figure
sign_number = zeros(7,1);
nonsign_number = zeros(7,1);
for i = 1:7
    sign_number(i) = length(find(ish_data_all(i).g1t2z_corr_q_fdr<0.05));
    total = length(ish_data_all(i).g1t2z_corr_q_fdr);
    nonsign_number(i) = total - sign_number(i);
end

for i = 1:7
    close all
    p = pie([sign_number(i),nonsign_number(i)]);
    cm = [[226,115,134]./255;0.5,0.5,0.5];
    colormap(cm);
    %     title(ish_data_all(i).name,'FontSize',7)
    p(2).Visible = 'off';
    p(4).Visible = 'off';
    set(gca,'LineWidth',0.5);
    set(gca,'FontName','Arial','FontSize',7);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    if i == 1
        %         legend('Location','eastoutside','String',{'FDR {\itq} < 0.05','n.s.'});
        set(gcf,'Paperposition',[1 1 4 4]);
    else
        set(gcf,'Paperposition',[1 1 2 2]);
    end
    print(gcf,[figure_dir,'pie_',ish_data_all(i).name,'_.tif'],'-dtiff','-r1000')
end
