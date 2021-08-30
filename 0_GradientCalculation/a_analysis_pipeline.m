%% process pipeline
%% define work dir
data_dir = 'D:\Data\DIDA-MDD\gradient_analysis\analysis2\';
script_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\3.Scripts\';
figure_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\5.Figures\ForRevision1\';

%% Downsample preprocessed images
cd([data_dir,'FunImgARWSDCFB']);
list = dir;
list(1:2) = [];

for i = 1:length(list)    
    cd([data_dir,'FunImgARWSDCFB\',list(i).name]);
    filename = dir('*.nii');
    x_reslice([script_dir,'Reslice_group_mask.nii'],filename.name,4);
    mkdir([data_dir,'Resliced\',list(i).name]);
    movefile('r*.nii',[data_dir,'\Resliced\',list(i).name]);
end

%% Generate connectivity matrix and calculate gradients
cd([data_dir,'Resliced']);
list = dir;
list(1:2) = [];

for i = 1:length(list)
% for i = [654,678]
    tic
    cd([data_dir,'Resliced\',list(i).name]);
    filename = dir('*.nii');
    M = x_gen_matrix_voxel([script_dir,'Reslice_group_mask.nii'],filename.name);
    n = length(M);
    M_spar = M;
    tmp = sort(M);
    tmp = M - repmat(tmp(round(n*0.9),:),n,1);
    M_spar(tmp<0) = 0;
    M_spar = M_spar';
    
    M_cos = 1 - squareform(pdist(M_spar,'cosine'));
    M_normalized = 1 - acos(M_cos)/pi;
    [embedding,result] = x_compute_diffusion_map(M_normalized,0.5,30);
    
    ind_dir = [data_dir,'Gradient_SameLength\',list(i).name];
    mkdir(ind_dir);
        filename = [ind_dir,'\gradient.mat'];
        save(filename,'embedding','result');
    toc
end

%% arrange emb
cd([data_dir,'Gradient_SameLength']);
list = dir;
list(1:2) = [];
n_sub = length(list);

emb_all = cell(n_sub,1);
res_all = cell(n_sub,1);
for i = 1:n_sub
    cd([data_dir,'Gradient_SameLength\',list(i).name]);
    load('gradient.mat');
    emb_all{i} = embedding;
    res_all{i} = result;    
end

% aligned across subjects
[realigned, xfms] = mica_iterativeAlignment(emb_all,100);
realigned = real(realigned);
xfms = cellfun(@real,xfms,'UniformOutput',false);

gradient_emb = cell(30,1);
for i = 1:30
    gradient_emb{i} = squeeze(realigned(:,i,:));
end

%calculate order sequence
seq = zeros(n_sub,30);
for i = 1:n_sub
    tmp = abs(xfms{i});
    [~,I] = sort(tmp,'descend');
    seq_tmp = zeros(5,1);
    seq_tmp(1) = I(1);
    for j = 2:30
        for k = 1:30
            if isempty(find(seq_tmp == I(k,j), 1))
                seq_tmp(j) = I(k,j);
                break;
            end
        end
    end
    seq(i,:) = seq_tmp;
end

%calculate explaination rate
exprate = zeros(n_sub,30);
for i = 1:n_sub
    tmp = res_all{i}.lambdas./sum(res_all{i}.lambdas);
    exprate(i,:) = tmp(seq(i,:));
end


% reorder gradient according to explaination ratio
exprate_mean = mean(exprate);
[~,I] = sort(exprate_mean,'descend');
gradient_emb_reordered = cell(30,1);
for i = 1:30
    gradient_emb_reordered{i} = gradient_emb{I(i)};
end

exprate_reordered = exprate(:,I);

% harmonization
sub_info = xlsread([script_dir,'All_sub_info.xlsx']);

gradient_emb_correct = cell(3,1);
for i = 1:3
    gradient_emb_correct{i} = combat(gradient_emb_reordered{i},sub_info(:,2)',sub_info(:,[1 3 4]));
end

% visual check mean gradient
hdr_mask = spm_vol([script_dir,'Reslice_group_mask.nii']);
vol_mask = spm_read_vols(hdr_mask);
hdr = hdr_mask;
hdr.dt(1) = 64;
ind = find(vol_mask);
cd(data_dir);
for i = 1:3
    vol = zeros(hdr.dim);
    vol(ind) = mean(gradient_emb_correct{i},2);
    hdr.fname = ['mean_g_',num2str(i),'.nii'];    
    spm_write_vol(hdr,vol);
end

for i = 1:3
    gradient_emb_correct{i} = -gradient_emb_correct{i};
end
save([data_dir,'\gradient_emb_correct(1-3).mat'],'gradient_emb_correct');

% generate mean map for both groups
hdr_mask = spm_vol([script_dir,'Reslice_group_mask.nii']);
vol_mask = spm_read_vols(hdr_mask);
hdr = hdr_mask;
hdr.dt(1) = 64;
ind = find(vol_mask);
cd(data_dir);
for i = 1:3
    vol = zeros(hdr.dim);
    vol(ind) = mean(gradient_emb_correct{i}(:,sub_info(:,1)==2),2);
    hdr.fname = ['mean_g_',num2str(i),'_mdd.nii'];    
    spm_write_vol(hdr,vol);
    vol(ind) = mean(gradient_emb_correct{i}(:,sub_info(:,1)==1),2);
    hdr.fname = ['mean_g_',num2str(i),'_hc.nii'];   
    spm_write_vol(hdr,vol);
end


% correct explaination
exprate_corrected = combat(exprate_reordered',sub_info(:,2)',sub_info(:,[1 3 4]));
exprate_corrected = exprate_corrected';
save([data_dir,'exprate_corrected.mat'], 'exprate_corrected');
mean(exprate_corrected(:,1))
std(exprate_corrected(:,1))
mean(exprate_corrected(sub_info(:,1)==2,1))
std(exprate_corrected(sub_info(:,1)==2,1))
mean(exprate_corrected(sub_info(:,1)==1,1))
std(exprate_corrected(sub_info(:,1)==1,1))


% draw explaination ratio
close all
mean_exprate_hc = mean(exprate_corrected(sub_info(:,1)==1,:));
mean_exprate_mdd = mean(exprate_corrected(sub_info(:,1)==2,:));
color_hc = [0.7,0.7,0.7];
color_mdd = [226,115,134]./255;
ph = plot(1.2:1:30.2,mean_exprate_hc,'.-','color',color_hc,'linewidth',0.5);
hold on
pm = plot(0.8:1:29.8,mean_exprate_mdd,'.-','color',color_mdd,'linewidth',0.5);
hold off
ph.MarkerSize = 10;
pm.MarkerSize = 10;
l = legend;
l.String = {'Controls','MDD'};
xlabel('Component of diffusion embedding');
ylabel('Explained ratio');
set(gca,'LineWidth',0.5);
set(gca,'FontName','Arial','FontSize',6);
set(gca,'YLim',[0,0.15]);
set(gca,'XLim',[0,31]);
box off
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 8.9 8]);
print(gcf,[figure_dir,'ComponentExpratio.tif'],'-dtiff','-r1000');

close all
cum_exprate_hc = cumsum(mean_exprate_hc);
cum_exprate_mdd = cumsum(mean_exprate_mdd);
ph = plot(1.2:1:30.2,cum_exprate_hc,'.-','color',color_hc,'linewidth',0.5);
hold on
pm = plot(0.8:1:29.8,cum_exprate_mdd,'.-','color',color_mdd,'linewidth',0.5);
hold off
ph.MarkerSize = 10;
pm.MarkerSize = 10;
l = legend;
l.String = {'Controls','MDD'};
l.Location = 'northeast';
xlabel('Component of diffusion embedding');
ylabel('Cumulative explained ratio');
set(gca,'LineWidth',0.5);
set(gca,'FontName','Arial','FontSize',6);
set(gca,'YLim',[0,1],'YTick',0:0.25:1);
set(gca,'XLim',[0,31]);
box off
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 8.9 8]);
print(gcf,[figure_dir,'CompCum.tif'],'-dtiff','-r1000');

% calculate gradient range
emb_range = zeros(n_sub,3);
for i = 1:3
    for j = 1:n_sub
        emb_range(j,i) = max(gradient_emb_correct{i}(:,j)) - min(gradient_emb_correct{i}(:,j));
    end
end
save([data_dir,'emb_range.mat'],'emb_range');

% calculate variation
emb_std = zeros(n_sub,3);
for i = 1:3
    for j = 1:n_sub
        emb_std(j,i) = std(gradient_emb_correct{i}(:,j));
    end
end
save([data_dir,'emb_std.mat'],'emb_std');
