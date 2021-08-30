%% process pipeline
%% define work dir
data_dir = 'D:\Data\DIDA-MDD\gradient_analysis\analysis2\';
script_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\3.Scripts\';
figure_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\5.Figures\ForRevision1\';



%% Between-group difference for global topology
% select subject
id_OAa21 = find(sub_info(:,1)==2&sub_info(:,7)>21);
id_OAb21 = find(sub_info(:,1)==2&((sub_info(:,7)<=21&sub_info(:,7)~=-1)|sub_info(:,3)<=21));

sub_id = [id_OAa21;id_OAb21];

n_a21 = length(id_OAa21);
n_b21 = length(id_OAb21);

des = [[ones(n_a21,1);zeros(n_b21,1)],sub_info(sub_id,3),sub_info(sub_id,4)];

% explained ratio
stat_exp = zeros(3,7);
for i = 1:3
    stat_exp(i,1) = mean(exprate_corrected(id_OAa21,i));
    stat_exp(i,2) = std(exprate_corrected(id_OAa21,i));
    stat_exp(i,3) = mean(exprate_corrected(id_OAb21,i));
    stat_exp(i,4) = std(exprate_corrected(id_OAb21,i));    
    stat_result = regstats(exprate_corrected(sub_id,i),des,'linear',{'tstat','r'});
    stat_exp(i,5) = stat_result.tstat.t(2);
    stat_exp(i,6) = stat_result.tstat.t(2) * sqrt(1/n_a21 + 1/n_b21);
    stat_exp(i,7) = stat_result.tstat.pval(2);
end
disp(stat_exp)

% gradient range
stat_range = zeros(3,7);
for i = 1:3
    stat_range(i,1) = mean(emb_range(id_OAa21,i));
    stat_range(i,2) = std(emb_range(id_OAa21,i));
    stat_range(i,3) = mean(emb_range(id_OAb21,i));
    stat_range(i,4) = std(emb_range(id_OAb21,i));    
    stat_result = regstats(emb_range(sub_id,i),des,'linear',{'tstat','r'});
    stat_range(i,5) = stat_result.tstat.t(2);
    stat_range(i,6) = stat_result.tstat.t(2) * sqrt(1/n_a21 + 1/n_b21);
    stat_range(i,7) = stat_result.tstat.pval(2);
end
disp(stat_range)

% gradient variance
stat_std = zeros(3,7);
for i = 1:3
    stat_std(i,1) = mean(emb_std(id_OAa21,i));
    stat_std(i,2) = std(emb_std(id_OAa21,i));
    stat_std(i,3) = mean(emb_std(id_OAb21,i));
    stat_std(i,4) = std(emb_std(id_OAb21,i));    
    stat_result = regstats(emb_std(sub_id,i),des,'linear',{'tstat','r'});
    stat_std(i,5) = stat_result.tstat.t(2);
    stat_std(i,6) = stat_result.tstat.t(2) * sqrt(1/n_a21 + 1/n_b21);
    stat_std(i,7) = stat_result.tstat.pval(2);
end
disp(stat_std)

save([data_dir,'BetweenGroupDiff_OA21\Topography_stat.mat'],'stat_exp','stat_range','stat_exp');

%% draw figure for range
% g1
close all
OAa_rang_G1 = emb_range(id_OAa21,1);
OAb_rang_G1 = emb_range(id_OAb21,1);
Data = cell(1,2);
Data{1} = OAa_rang_G1;
Data{2} = OAb_rang_G1;
color_OAa = [203,142,133]./255;
color_OAb = [242,222,189]./255;
gretna_plot_violin(Data, {'',''}, {'G1 rang'}, 'meanstdfill');
h = gca;
h.Children(4).FaceColor = color_OAb;
h.Children(5).FaceColor = color_OAa;
hold on
plot([1,2],[0.38,0.38],'-k','LineWidth',0.3);
set(gca,'YLim',[0.08,0.39],'YTick',[0.1 0.2 0.3]);
hold off
t1 = text(1.5,0.385,'***','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',8);
set(gca,'XTick',[1,2]);
set(gca,'XTickLabel',{});
ylabel('Gradient range');
legend off
set(gca,'Linewidth',0.5);
set(gca,'FontName','Arial','FontSize',7);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 4 4.5]);
print(gcf,[figure_dir,'OA21_G1_range.tif'],'-dtiff','-r1000');

%g2
close all
OAa_rang_G2 = emb_range(id_OAa21,2);
OAb_rang_G2 = emb_range(id_OAb21,2);
Data = cell(1,2);
Data{1} = OAa_rang_G2;
Data{2} = OAb_rang_G2;
color_OAa = [203,142,133]./255;
color_OAb = [242,222,189]./255;
gretna_plot_violin(Data, {'',''}, {'G2 rang'}, 'meanstdfill');
h = gca;
h.Children(4).FaceColor = color_OAb;
h.Children(5).FaceColor = color_OAa;
hold on
plot([1,2],[0.32,0.32],'-k','LineWidth',0.3);
set(gca,'YLim',[0.08,0.33],'YTick',[0.1 0.2 0.3 ]);
hold off
t1 = text(1.5,0.325,'**','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',8);
set(gca,'XTick',[1,2]);
set(gca,'XTickLabel',{'OA>21','OA≤21'});
ylabel('Gradient range');
legend off
set(gca,'Linewidth',0.5);
set(gca,'FontName','Arial','FontSize',6);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 4 4.5]);
print(gcf,[figure_dir,'OA21_G2_range.tif'],'-dtiff','-r1000');

%% draw figure for variation
%g1
close all
OAa_var_G1 = emb_std(id_OAa21,1);
OAb_var_G1 = emb_std(id_OAb21,1);
Data = cell(1,2);
Data{1} = OAa_var_G1;
Data{2} = OAb_var_G1;
color_OAa = [203,142,133]./255;
color_OAb = [242,222,189]./255;
gretna_plot_violin(Data, {' ',' '}, {'G1 variance'}, 'meanstdfill');
h = gca;
h.Children(4).FaceColor = color_OAb;
h.Children(5).FaceColor = color_OAa;
hold on
plot([1,2],[0.11,0.11],'-k','LineWidth',0.3);
set(gca,'YLim',[0,0.12],'YTick',[0,0.1 ]);
hold off
t1 = text(1.5,0.115,'*','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',8);
set(gca,'XTick',[1,2]);
set(gca,'XTickLabel',{});
ylabel('Gradient variance');
legend off
set(gca,'Linewidth',0.5);
set(gca,'FontName','Arial','FontSize',6);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 4 4.5]);
print(gcf,[figure_dir,'OA21_G1_var.tif'],'-dtiff','-r1000')

%g2
close all
OAa_var_G1 = emb_std(id_OAa21,2);
OAb_var_G1 = emb_std(id_OAb21,2);
Data = cell(1,2);
Data{1} = OAa_var_G1;
Data{2} = OAb_var_G1;
color_OAa = [203,142,133]./255;
color_OAb = [242,222,189]./255;
gretna_plot_violin(Data, {' ',' '}, {'G2 variance'}, 'meanstdfill');
h = gca;
h.Children(4).FaceColor = color_OAb;
h.Children(5).FaceColor = color_OAa;
hold on
plot([1,2],[0.11,0.11],'-k','LineWidth',0.3);
set(gca,'YLim',[0,0.12],'YTick',[0,0.1 ]);
hold off
t1 = text(1.5,0.115,'*','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',8);
set(gca,'XTick',[1,2]);
set(gca,'XTickLabel',{'OA>21','OA≤21'});
ylabel('Gradient variance');
legend off
set(gca,'Linewidth',0.5);
set(gca,'FontName','Arial','FontSize',6);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 4 4.5]);
print(gcf,[figure_dir,'OA21_G2_var.tif'],'-dtiff','-r1000')


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
cd([data_dir,'BetweenGroupDiff_OA21']);

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
