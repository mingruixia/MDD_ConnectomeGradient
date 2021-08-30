%% process pipeline
%% define work dir
data_dir = 'D:\Data\DIDA-MDD\gradient_analysis\analysis2\';
script_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\3.Scripts\';
figure_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\5.Figures\ForRevision1\';

%% sub-network differences in mean maps
load([data_dir,'gradient_emb_correct(1-3).mat']);
id_hc = sub_info(:,1)==1;
id_mdd = sub_info(:,1)==2;
mean_gradient_mdd = mean(gradient_emb_correct{1}(:,id_mdd),2)';
mean_gradient_hc = mean(gradient_emb_correct{1}(:,id_hc),2)';
vol_mod = spm_read_vols(spm_vol([script_dir,'mod_mask.nii']));
ind_mod = vol_mod(find(vol_mod));
gradient_score_subnetwork = cell(8,2);
for i = 1:8
    gradient_score_subnetwork{i,1} = mean_gradient_mdd(ind_mod==i);
    gradient_score_subnetwork{i,2} = mean_gradient_hc(ind_mod==i);
end

% test difference
t = zeros(8,1);
p = zeros(8,1);
df = zeros(8,1);
cohen_d = zeros(8,1);
for i = 1:8
    [h,p(i),ci,stats] = ttest(gradient_score_subnetwork{i,1},gradient_score_subnetwork{i,2});
    t(i) = stats.tstat;
    df(i) = stats.df;
    cohen_d(i) = mean(gradient_score_subnetwork{i,1}-gradient_score_subnetwork{i,2})/stats.sd;
end
[~,~,rnk] = unique(p);
q_fdr = (p * 8) ./ rnk;
mod_name = {'VIS','SMN','DAN','VAN','LIM','FPN','DMN','SUB'};
save([data_dir,'\SubNetworkDiffinMeanMap.mat'],'t','p','df','mod_name','q_fdr','cohen_d');

%% check skewness and kurtosis
id_hc = sub_info(:,1)==1;
id_mdd = sub_info(:,1)==2;
skew_stat = zeros(3,3,2);
for i = 1:3
skew_stat(i,1,1) = skewness(exprate_corrected(id_mdd,i));
skew_stat(i,1,2) = skewness(exprate_corrected(id_hc,i));
skew_stat(i,2,1) = skewness(emb_range(id_mdd,i));
skew_stat(i,2,2) = skewness(emb_range(id_hc,i));
skew_stat(i,3,1) = skewness(emb_std(id_mdd,i));
skew_stat(i,3,2) = skewness(emb_std(id_hc,i));
end

skew_voxel = zeros(3,18933,2);
for i = 1:3
    for j = 1:18933
        skew_voxel(i,j,1) = skewness(gradient_emb_correct{i}(j,id_mdd));
        skew_voxel(i,j,2) = skewness(gradient_emb_correct{i}(j,id_hc));
    end
end

id_hc = sub_info(:,1)==1;
id_mdd = sub_info(:,1)==2;
kurt_stat = zeros(3,3,2);
for i = 1:3
kurt_stat(i,1,1) = kurtosis(exprate_corrected(id_mdd,i));
kurt_stat(i,1,2) = kurtosis(exprate_corrected(id_hc,i));
kurt_stat(i,2,1) = kurtosis(emb_range(id_mdd,i));
kurt_stat(i,2,2) = kurtosis(emb_range(id_hc,i));
kurt_stat(i,3,1) = kurtosis(emb_std(id_mdd,i));
kurt_stat(i,3,2) = kurtosis(emb_std(id_hc,i));
end

kurt_voxel = zeros(3,18933,2);
for i = 1:3
    for j = 1:18933
        kurt_voxel(i,j,1) = kurtosis(gradient_emb_correct{i}(j,id_mdd));
        kurt_voxel(i,j,2) = kurtosis(gradient_emb_correct{i}(j,id_hc));
    end
end

%% check variance
var_stat = zeros(3,3);
var_p = zeros(3,3);
for i = 1:3
[~,var_p(i,1),~,stat] = vartest2(exprate_corrected(id_mdd,i),exprate_corrected(id_hc,i));
var_stat(i,1) = stat.fstat;
[~,var_p(i,2),~,stat] = vartest2(emb_range(id_mdd,i),emb_range(id_hc,i));
var_stat(i,2) = stat.fstat;
[~,var_p(i,3),~,stat] = vartest2(emb_std(id_mdd,i),emb_std(id_hc,i));
var_stat(i,3) = stat.fstat;
end

%% Between-group difference for global topology
% explained ratio
stat = zeros(3,7);
id_hc = sub_info(:,1)==1;
id_mdd = sub_info(:,1)==2;
des = [sub_info(:,1),sub_info(:,3:4)];
n_hc = length(find(id_hc));
n_mdd = length(find(id_mdd));
for i = 1:3
    stat(i,1) = mean(exprate_corrected(id_mdd,i));
    stat(i,2) = std(exprate_corrected(id_mdd,i));
    stat(i,3) = mean(exprate_corrected(id_hc,i));
    stat(i,4) = std(exprate_corrected(id_hc,i));    
    stat_result = regstats(exprate_corrected(:,i),des,'linear',{'tstat','r'});
    stat(i,5) = stat_result.tstat.t(2);
    stat(i,6) = stat_result.tstat.t(2) * sqrt(1/n_hc + 1/n_mdd);
    stat(i,7) = stat_result.tstat.pval(2);
end
disp(stat)

% gradient range
load([data_dir,'emb_range.mat']);
stat = zeros(3,7);
id_hc = sub_info(:,1)==1;
id_mdd = sub_info(:,1)==2;
des = [sub_info(:,1),sub_info(:,3:4)];
n_hc = length(find(id_hc));
n_mdd = length(find(id_mdd));
for i = 1:3
    stat(i,1) = mean(emb_range(id_mdd,i));
    stat(i,2) = std(emb_range(id_mdd,i));
    stat(i,3) = mean(emb_range(id_hc,i));
    stat(i,4) = std(emb_range(id_hc,i));    
    stat_result = regstats(emb_range(:,i),des,'linear',{'tstat','r'});
    stat(i,5) = stat_result.tstat.t(2);
    stat(i,6) = stat_result.tstat.t(2) * sqrt(1/n_hc + 1/n_mdd);
    stat(i,7) = stat_result.tstat.pval(2);
end
disp(stat)

% gradient variance
load([data_dir,'emb_std.mat']);
stat = zeros(3,7);
id_hc = sub_info(:,1)==1;
id_mdd = sub_info(:,1)==2;
des = [sub_info(:,1),sub_info(:,3:4)];
n_hc = length(find(id_hc));
n_mdd = length(find(id_mdd));
for i = 1:3
    stat(i,1) = mean(emb_std(id_mdd,i));
    stat(i,2) = std(emb_std(id_mdd,i));
    stat(i,3) = mean(emb_std(id_hc,i));
    stat(i,4) = std(emb_std(id_hc,i));    
    stat_result = regstats(emb_std(:,i),des,'linear',{'tstat','r'});
    stat(i,5) = stat_result.tstat.t(2);
    stat(i,6) = stat_result.tstat.t(2) * sqrt(1/n_hc + 1/n_mdd);
    stat(i,7) = stat_result.tstat.pval(2);
end
disp(stat)

%% draw figure for exp
% g1
close all
HC_exp_G1 = exprate_corrected(id_hc);
MDD_exp_G1 = exprate_corrected(id_mdd);
Data = cell(1,2);
Data{1} = MDD_exp_G1;
Data{2} = HC_exp_G1;
color_mdd = [226,115,134]./255;
color_hc = [0.7,0.7,0.7];
gretna_plot_violin(Data, {'MDD','HC'}, {'Explain ratio'}, 'meanstdfill');
h = gca;
h.Children(4).FaceColor = color_hc;
h.Children(5).FaceColor = color_mdd;
hold on
plot([1,2],[0.3,0.3],'-k','LineWidth',0.3);
set(gca,'YLim',[0,0.32],'YTick',[0,0.1 0.2 0.3]);
hold off
t1 = text(1.5,0.305,'***','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',8);
set(gca,'XTick',[1,2]);
set(gca,'XTickLabel',{});
ylabel('Gradient explained ratio');
legend off
set(gca,'LineWidth',0.5);
set(gca,'FontName','Arial','FontSize',6);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 3 3.5]);
print(gcf,[figure_dir,'ExpVio1.tif'],'-dtiff','-r1000')

%g2
close all
HC_exp_G2 = exprate_corrected(id_hc,2);
MDD_exp_G2 = exprate_corrected(id_mdd,2);
Data = cell(1,2);
Data{1} = MDD_exp_G2;
Data{2} = HC_exp_G2;
color_mdd = [226,115,134]./255;
color_hc = [0.7,0.7,0.7];
gretna_plot_violin(Data, {'MDD','HC'}, {'Explain ratio'}, 'meanstdfill');
h = gca;
h.Children(4).FaceColor = color_hc;
h.Children(5).FaceColor = color_mdd;
set(gca,'YLim',[0,0.32],'YTick',[0,0.1 0.2 0.3]);
t1 = text(1.5,0.305,'n.s.','FontName','Arial','HorizontalAlignment','Center','FontSize',5);
set(gca,'XTick',[1,2]);
set(gca,'XTickLabel',{});
ylabel('Explained ratio');
legend off
set(gca,'LineWidth',0.5);
set(gca,'FontName','Arial','FontSize',5);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 2 1.8]);
print(gcf,[figure_dir,'ExpVio2.tif'],'-dtiff','-r1000')

%g3
close all
HC_exp_G3 = exprate_corrected(id_hc,3);
MDD_exp_G3 = exprate_corrected(id_mdd,3);
Data = cell(1,2);
Data{1} = MDD_exp_G3;
Data{2} =  HC_exp_G3;
color_mdd = [226,115,134]./255;
color_hc = [0.7,0.7,0.7];
gretna_plot_violin(Data, {'MDD','HC'}, {'Explain ratio'}, 'meanstdfill');
h = gca;
h.Children(4).FaceColor = color_hc;
h.Children(5).FaceColor = color_mdd;
hold on
plot([1,2],[0.21,0.21],'-k','LineWidth',0.3);
set(gca,'YLim',[0,0.23],'YTick',[0,0.1 0.2 ]);
hold off
t1 = text(1.5,0.215,'*','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',8);
set(gca,'XTick',[1,2]);
set(gca,'XTickLabel',{});
ylabel('Explained ratio');
legend off
set(gca,'LineWidth',0.5);
set(gca,'FontName','Arial','FontSize',5);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 2 1.8]);
print(gcf,[figure_dir,'ExpVio3.tif'],'-dtiff','-r1000')

%% draw figure for range
%g1
close all
HC_rang_G1 = emb_range(id_hc,1);
MDD_rang_G1 = emb_range(id_mdd,1);
Data = cell(1,2);
Data{1} = MDD_rang_G1;
Data{2} = HC_rang_G1;
color_mdd = [226,115,134]./255;
color_hc = [0.7,0.7,0.7];
gretna_plot_violin(Data, {'MDD','HC'}, {'G1 rang'}, 'meanstdfill');
h = gca;
h.Children(4).FaceColor = color_hc;
h.Children(5).FaceColor = color_mdd;
hold on
plot([1,2],[0.38,0.38],'-k','LineWidth',0.3);
set(gca,'YLim',[0.08,0.39],'YTick',[0.1 0.2 0.3]);
hold off
t1 = text(1.5,0.385,'***','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',8);
set(gca,'XTick',[1,2]);
set(gca,'XTickLabel',{});
ylabel('Gradient range');
legend off
set(gca,'FontName','Arial','FontSize',6);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 3 3.5]);
print(gcf,[figure_dir,'RangVio1.tif'],'-dtiff','-r1000')

%g2
close all
HC_rang_G1 = emb_range(id_hc,2);
MDD_rang_G1 = emb_range(id_mdd,2);
Data = cell(1,2);
Data{1} = MDD_rang_G1;
Data{2} = HC_rang_G1;
color_mdd = [226,115,134]./255;
color_hc = [0.7,0.7,0.7];
gretna_plot_violin(Data, {'MDD','HC'}, {'G1 rang'}, 'meanstdfill');
h = gca;
h.Children(4).FaceColor = color_hc;
h.Children(5).FaceColor = color_mdd;
hold on
plot([1,2],[0.32,0.32],'-k','LineWidth',0.3);
set(gca,'YLim',[0.08,0.33],'YTick',[0.1 0.2 0.3 ]);
hold off
t1 = text(1.5,0.325,'***','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',8);
set(gca,'XTick',[1,2]);
set(gca,'XTickLabel',{});
ylabel('Gradient range');
legend off
set(gca,'FontName','Arial','FontSize',5);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 2 1.8]);
print(gcf,[figure_dir,'RangVio2.tif'],'-dtiff','-r1000')

%g3
close all
HC_rang_G1 = emb_range(id_hc,3);
MDD_rang_G1 = emb_range(id_mdd,3);
Data = cell(1,2);
Data{1} = MDD_rang_G1;
Data{2} = HC_rang_G1;
color_mdd = [226,115,134]./255;
color_hc = [0.7,0.7,0.7];
gretna_plot_violin(Data, {'MDD','HC'}, {'G3 rang'}, 'meanstdfill');
h = gca;
h.Children(4).FaceColor = color_hc;
h.Children(5).FaceColor = color_mdd;
hold on
plot([1,2],[0.28,0.28],'-k','LineWidth',0.3);
set(gca,'YLim',[0.05,0.3],'YTick',[0.1 0.2 ]);
hold off
t1 = text(1.5,0.285,'***','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',8);
set(gca,'XTick',[1,2]);
set(gca,'XTickLabel',{});
ylabel('Gradient range');
legend off
set(gca,'FontName','Arial','FontSize',5);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 2 1.8 ]);
print(gcf,[figure_dir,'RangVio3.tif'],'-dtiff','-r1000')

%% draw figure for variation
%g1
close all
HC_rang_G1 = emb_std(id_hc,1);
MDD_rang_G1 = emb_std(id_mdd,1);
Data = cell(1,2);
Data{1} = MDD_rang_G1;
Data{2} = HC_rang_G1;
color_mdd = [226,115,134]./255;
color_hc = [0.7,0.7,0.7];
gretna_plot_violin(Data, {'MDD','HC'}, {'G3 rang'}, 'meanstdfill');
h = gca;
h.Children(4).FaceColor = color_hc;
h.Children(5).FaceColor = color_mdd;
hold on
plot([1,2],[0.11,0.11],'-k','LineWidth',0.3);
set(gca,'YLim',[0,0.12],'YTick',[0,0.1 ]);
hold off
t1 = text(1.5,0.115,'***','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',8);
set(gca,'XTick',[1,2]);
set(gca,'XTickLabel',{});
ylabel('Gradient variance');
legend off
set(gca,'FontName','Arial','FontSize',7);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 3 3.5]);
print(gcf,[figure_dir,'StdVio1.tif'],'-dtiff','-r1000')

%g2
close all
HC_std_G2 = emb_std(id_hc,2);
MDD_std_G2 = emb_std(id_mdd,2);
Data = cell(1,2);
Data{1} = MDD_std_G2;
Data{2} = HC_std_G2;
color_hc = [0.7,0.7,0.7];
color_mdd = [226,115,134]./255;
gretna_plot_violin(Data, {'MDD','HC'}, {'G3 rang'}, 'meanstdfill');
h = gca;
h.Children(4).FaceColor = color_hc;
h.Children(5).FaceColor = color_mdd;
set(gca,'YLim',[0,0.12],'YTick',[0,0.1 ]);
t1 = text(1.5,0.105,'n.s.','FontName','Arial','HorizontalAlignment','Center','FontSize',5);
set(gca,'XTick',[1,2]);
set(gca,'XTickLabel',{});
ylabel('Gradient variation');
legend off
set(gca,'FontName','Arial','FontSize',5);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 2 1.8]);
print(gcf,[figure_dir,'StdVio2.tif'],'-dtiff','-r1000')

%g3
close all
HC_rang_G1 = emb_std(id_hc,3);
MDD_rang_G1 = emb_std(id_mdd,3);
Data = cell(1,2);
Data{1} = MDD_rang_G1;
Data{2} = HC_rang_G1;
color_mdd = [226,115,134]./255;
color_hc = [0.7,0.7,0.7];
gretna_plot_violin(Data, {'MDD','HC'}, {'G3 rang'}, 'meanstdfill');
h = gca;
h.Children(4).FaceColor = color_hc;
h.Children(5).FaceColor = color_mdd;
hold on
plot([1,2],[0.09,0.09],'-k','LineWidth',0.3);
set(gca,'YLim',[0,0.11],'YTick',[0,0.1]);
hold off
t1 = text(1.5,0.095,'***','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',8);
set(gca,'XTick',[1,2]);
set(gca,'XTickLabel',{});
ylabel('Gradient variance');
legend off
set(gca,'FontName','Arial','FontSize',5);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 2 1.8]);
print(gcf,[figure_dir,'StdVio3.tif'],'-dtiff','-r1000');

%% Reginal difference
% estimate group difference
[n_voxel,n_sub] = size(gradient_emb_correct{1});
T_value = zeros(3,n_voxel);
res_map = cell(3,1);
des = [sub_info(:,1),sub_info(:,3:4)];
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
mkdir([data_dir,'BetweenGroupDiff']);
cd([data_dir,'BetweenGroupDiff']);

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

%% run python code of Brainsmash to generate surrogate maps for between-group z maps



%% between-group difference among moduele
pos_perc_real = zeros(8,1);
neg_perc_real = zeros(8,1);

z_real = load([data_dir,'\BetweenGroupDiff\g1_T2_Z.txt']);
hdr_mask = spm_vol([script_dir,'Reslice_group_mask.nii']);
vol_mask = spm_read_vols(hdr_mask);
ind = find(vol_mask);

vol_real_z = zeros(hdr_mask.dim);
vol_real_z(ind) = z_real;

zThrd = norminv(1 - 0.001/2);
vol_real_z(vol_real_z < zThrd & vol_real_z > -zThrd) = 0;
[L, num] = bwlabeln(vol_real_z,26);
n = 0;
for x = 1:num
    theCurrentCluster = L == x;
    if length(find(theCurrentCluster)) <= 26
        n = n + 1;
        vol_real_z(logical(theCurrentCluster)) = 0;
    end
end

hdr_mod = spm_vol([script_dir,'mod_mask.nii']);
vol_mod = spm_read_vols(hdr_mod);

n_total_pos = length(find(vol_real_z>0));
n_total_neg = length(find(vol_real_z<0));

for i = 1:8
    mod_voxel = vol_real_z(vol_mod == i);
    pos_perc_real(i) = length(find(mod_voxel>0)) / n_total_pos;
    neg_perc_real(i) = length(find(mod_voxel<0)) / n_total_neg;
end

pos_perc_surr = zeros(8,10000);
neg_perc_surr = zeros(8,10000);
load([data_dir,'surrogate_maps_g1_z\surrogate_maps_g1_z_resample.mat']);
parfor j = 1:10000
    vol_surr = zeros(hdr_mask.dim);
    vol_surr(ind) = surrogate_maps(j,:);
    vol_surr(vol_surr < zThrd & vol_surr > -zThrd) = 0;
    [L, num] = bwlabeln(vol_surr,26);
    n = 0;
    for x = 1:num
        theCurrentCluster = L == x;
        if length(find(theCurrentCluster)) <= 26
            n = n + 1;
            vol_surr(logical(theCurrentCluster)) = 0;
        end
    end
    
    n_total_pos = length(find(vol_real_z>0));
    n_total_neg = length(find(vol_real_z<0));
    for i = 1:8
        mod_voxel = vol_surr(vol_mod == i);
        pos_perc_surr(i,j) = length(find(mod_voxel>0)) / n_total_pos;
        neg_perc_surr(i,j) = length(find(mod_voxel<0)) / n_total_neg;
    end
end

p_pos_perc = zeros(8,1);
p_neg_perc = zeros(8,1);
for i = 1:8
    p_pos_perc(i) = length(find(pos_perc_surr(i,:)>pos_perc_real(i)))/10000;
    p_neg_perc(i) = length(find(neg_perc_surr(i,:)>neg_perc_real(i)))/10000;
end

p_all = [p_neg_perc;p_pos_perc];
[~,~,rnk] = unique(p_all);
q_fdr = (p_all * 16) ./ rnk;

save([data_dir,'Z_module_distribtuion.mat'],'pos_perc_real','neg_perc_real',...
    'pos_perc_surr','neg_perc_surr','p_pos_perc','p_neg_perc','q_fdr');

% draw figure for modular distribution of z
mod_color = [162 81 172;120 154 192;64 152 50;224 101 254;221.4000  227.7000  180.9000;239 185 67;217 113 125;118 113 113]/255;
mod_label = {'VIS','SMN','DAN','VAN','LIB','FPN','DMN','SUB'};
for i = 1:8
    close all
    h_pos = histfit(pos_perc_surr(i,:),[],'kernel');
    h_pos(1).Visible = 'off';
    h_pos(2).Color = mod_color(i,:);
    h_pos(2).LineWidth = 0.5;
    hold on
    yFill = [h_pos(2).YData,zeros(1,length(h_pos(2).YData))];
    xFill = [h_pos(2).XData,fliplr(h_pos(2).XData)];
    fill(xFill,yFill,mod_color(i,:),'FaceAlpha',0.5,'EdgeColor',mod_color(i,:));
    p_pos = plot(pos_perc_real(i),0,'o');
    p_pos.MarkerFaceColor = mod_color(i,:);
    p_pos.MarkerEdgeColor = [0 0 0];
    p_pos.MarkerSize = 3;
    p_pos.LineWidth = 0.5;
    hold off
    set(gca,'XLim',[0,0.8]);
    set(gca,'XTick',[0,0.4,0.8]);
    set(gca,'XTickLabel',{'0%','40%','80%'});
    if i == 8
        set(gca,'YLim',[0,950]);
        set(gca,'YTick',[0,950]);
    else
        set(gca,'YLim',[0,400]);
        set(gca,'YTick',[0,400]);
    end
    view(-90, 90) %# Swap the axes
    set(gca, 'ydir', 'reverse');
    set(gca,'FontName','Arial','FontSize',7);
    set(gca,'box','off');
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'Paperposition',[0 0 1.3 2]);
    grid off
    print(gcf,[figure_dir,'pos_t2_dist_',mod_label{i},'.tif'],'-dtiff','-r1000')    
    
    close all
    h_neg = histfit(neg_perc_surr(i,:),[],'kernel');
    h_neg(1).Visible = 'off';
    h_neg(2).Color = mod_color(i,:);
    h_neg(2).LineWidth = 0.5;
    hold on
    yFill = [h_neg(2).YData,zeros(1,length(h_neg(2).YData))];
    xFill = [h_neg(2).XData,fliplr(h_neg(2).XData)];
    fill(xFill,yFill,mod_color(i,:),'FaceAlpha',0.5,'EdgeColor',mod_color(i,:));
    p_neg = plot(neg_perc_real(i),0,'o');
    p_neg.MarkerFaceColor = mod_color(i,:);
    p_neg.MarkerEdgeColor = [0 0 0];
    p_neg.MarkerSize = 3;
    p_neg.LineWidth = 0.5;
    hold off
    set(gca,'XLim',[0,0.8]);
    set(gca,'XTick',[0,0.4,0.8]);
    set(gca,'XTickLabel',{'0%','40%','80%'});
    if i == 8
        set(gca,'YLim',[0,950]);
        set(gca,'YTick',[0,950]);
    else
        set(gca,'YLim',[0,400]);
        set(gca,'YTick',[0,400]);
    end
    view(-90, 90) %# Swap the axes
    set(gca, 'ydir', 'reverse');
    set(gca,'FontName','Arial','FontSize',7);
    set(gca,'box','off');
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'Paperposition',[0 0 1.3 2]);
    grid off
    print(gcf,[figure_dir,'neg_t2_dist_',mod_label{i},'.tif'],'-dtiff','-r1000')
end

%% Check group-to-sex interaction
% explained ratio
stat = zeros(3,2);
des = [sub_info(:,1),sub_info(:,3:4)];
des(:,4) = des(:,1).*des(:,3);
for i = 1:3  
    stat_result = regstats(exprate_corrected(:,i),des,'linear',{'tstat','r'});
    stat(i,1) = stat_result.tstat.t(5);
    stat(i,2) = stat_result.tstat.pval(5);
end
disp(stat)

% gradient range
for i = 1:3
    stat_result = regstats(emb_range(:,i),des,'linear',{'tstat','r'});
    stat(i,1) = stat_result.tstat.t(5);
    stat(i,2) = stat_result.tstat.pval(5);    
end
disp(stat)

% gradient variance
for i = 1:3
    stat_result = regstats(emb_std(:,i),des,'linear',{'tstat','r'});
    stat(i,1) = stat_result.tstat.t(5);
    stat(i,2) = stat_result.tstat.pval(5);    
end
disp(stat)


% estimate group difference
[n_voxel,n_sub] = size(gradient_emb_correct{1});
T_value = zeros(3,n_voxel);
res_map = cell(3,1);

for i = 1:3
    res_map{i} = zeros(n_voxel,n_sub);
    for j = 1:n_voxel
        stat_result = regstats(gradient_emb_correct{i}(j,:)',des,'linear',{'tstat','r'});
        T_value(i,j) = stat_result.tstat.t(5);
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
mkdir([data_dir,'GroupSexInteraction']);
cd([data_dir,'GroupSexInteraction']);

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
    hdr_z.fname = ['g',num2str(i),'_inter_z.nii'];
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
    hdr_z.fname = ['g',num2str(i),'_inter_z_cluster.nii'];
    spm_write_vol(hdr_z,Z_vol);
    niifile = hdr_z.fname;
    txtfile = ['g',num2str(i),'_inter_z_cluster.txt'];
    
    hdr_z.fname = ['g',num2str(i),'_inter_z_cluster_mask.nii'];
    Z_vol(Z_vol~=0) = 1;
    spm_write_vol(hdr_z,Z_vol);
    
    BrainNet_nii2txt(surfacefile,niifile,txtfile);
    BrainNet_MapCfg('D:\Dropbox\Projects\BrainNet\otherfile\Surface\FSaverage_inflated_32K.nv',...
        txtfile,[script_dir,'BNVCfg_T2_cluster.mat'],['g',num2str(i),'_inter_z_cluster.tif']);
end

