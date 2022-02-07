%% process pipeline
%% define work dir
data_dir = 'D:\Data\DIDA-MDD\gradient_analysis\analysis2\';
script_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\3.Scripts\';
figure_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\5.Figures\ForRevision3\';



%% Between-group difference for global topology of g1 (range, std)
%% control for mean FD
sub_info = xlsread([script_dir,'All_sub_info.xlsx']);
% select subject
id_OAa21 = find(sub_info(:,1)==2&sub_info(:,7)>21);
id_OAb21 = find(sub_info(:,1)==2&((sub_info(:,7)<=21&sub_info(:,7)~=-1)|sub_info(:,3)<=21));

sub_id = [id_OAa21;id_OAb21];

n_a21 = length(id_OAa21);
n_b21 = length(id_OAb21);

des = [[ones(n_a21,1);zeros(n_b21,1)],sub_info(sub_id,3),sub_info(sub_id,4),sub_info(sub_id,12)];
% gradient range
load([data_dir,'emb_range.mat']);
stat_range = zeros(1,7);
    stat_range(1) = mean(emb_range(id_OAa21,1));
    stat_range(2) = std(emb_range(id_OAa21,1));
    stat_range(3) = mean(emb_range(id_OAb21,1));
    stat_range(4) = std(emb_range(id_OAb21,1));    
    stat_result = regstats(emb_range(sub_id,1),des,'linear',{'tstat','r'});
    stat_range(5) = stat_result.tstat.t(2);
    stat_range(6) = stat_result.tstat.t(2) * sqrt(1/n_a21 + 1/n_b21);
    stat_range(7) = stat_result.tstat.pval(2);
disp(stat_range)

% gradient variance
load([data_dir,'emb_std.mat']);
stat_std = zeros(1,7);
    stat_std(1,1) = mean(emb_std(id_OAa21,1));
    stat_std(1,2) = std(emb_std(id_OAa21,1));
    stat_std(1,3) = mean(emb_std(id_OAb21,1));
    stat_std(1,4) = std(emb_std(id_OAb21,1));    
    stat_result = regstats(emb_std(sub_id,1),des,'linear',{'tstat','r'});
    stat_std(1,5) = stat_result.tstat.t(2);
    stat_std(1,6) = stat_result.tstat.t(2) * sqrt(1/n_a21 + 1/n_b21);
    stat_std(1,7) = stat_result.tstat.pval(2);
disp(stat_std)

%% Joint embedding
% gradient range
des = [[ones(n_a21,1);zeros(n_b21,1)],sub_info(sub_id,3),sub_info(sub_id,4)];
load([data_dir,'Validation_jointemb\emb_range.mat']);
stat_range = zeros(1,7);
    stat_range(1) = mean(emb_range(id_OAa21,1));
    stat_range(2) = std(emb_range(id_OAa21,1));
    stat_range(3) = mean(emb_range(id_OAb21,1));
    stat_range(4) = std(emb_range(id_OAb21,1));    
    stat_result = regstats(emb_range(sub_id,1),des,'linear',{'tstat','r'});
    stat_range(5) = stat_result.tstat.t(2);
    stat_range(6) = stat_result.tstat.t(2) * sqrt(1/n_a21 + 1/n_b21);
    stat_range(7) = stat_result.tstat.pval(2);
disp(stat_range)

% gradient variance
load([data_dir,'Validation_jointemb\emb_std.mat']);
stat_std = zeros(1,7);
    stat_std(1,1) = mean(emb_std(id_OAa21,1));
    stat_std(1,2) = std(emb_std(id_OAa21,1));
    stat_std(1,3) = mean(emb_std(id_OAb21,1));
    stat_std(1,4) = std(emb_std(id_OAb21,1));    
    stat_result = regstats(emb_std(sub_id,1),des,'linear',{'tstat','r'});
    stat_std(1,5) = stat_result.tstat.t(2);
    stat_std(1,6) = stat_result.tstat.t(2) * sqrt(1/n_a21 + 1/n_b21);
    stat_std(1,7) = stat_result.tstat.pval(2);
disp(stat_std)

%% Fishers Z Matrix
% gradient range
load([data_dir,'Validation_zMat\emb_range.mat']);
stat_range = zeros(1,7);
    stat_range(1) = mean(emb_range(id_OAa21,1));
    stat_range(2) = std(emb_range(id_OAa21,1));
    stat_range(3) = mean(emb_range(id_OAb21,1));
    stat_range(4) = std(emb_range(id_OAb21,1));    
    stat_result = regstats(emb_range(sub_id,1),des,'linear',{'tstat','r'});
    stat_range(5) = stat_result.tstat.t(2);
    stat_range(6) = stat_result.tstat.t(2) * sqrt(1/n_a21 + 1/n_b21);
    stat_range(7) = stat_result.tstat.pval(2);
disp(stat_range)

% gradient variance
load([data_dir,'Validation_zMat\emb_std.mat']);
stat_std = zeros(1,7);
    stat_std(1,1) = mean(emb_std(id_OAa21,1));
    stat_std(1,2) = std(emb_std(id_OAa21,1));
    stat_std(1,3) = mean(emb_std(id_OAb21,1));
    stat_std(1,4) = std(emb_std(id_OAb21,1));    
    stat_result = regstats(emb_std(sub_id,1),des,'linear',{'tstat','r'});
    stat_std(1,5) = stat_result.tstat.t(2);
    stat_std(1,6) = stat_result.tstat.t(2) * sqrt(1/n_a21 + 1/n_b21);
    stat_std(1,7) = stat_result.tstat.pval(2);
disp(stat_std)

%% in participants with mFD<0.25
sub_info = xlsread([script_dir,'All_sub_info.xlsx']);
% select subject
id_OAa21 = find(sub_info(:,1)==2&sub_info(:,7)>21&sub_info(:,12)<0.25);
id_OAb21 = find(sub_info(:,1)==2&sub_info(:,12)<0.25&((sub_info(:,7)<=21&sub_info(:,7)~=-1)|sub_info(:,3)<=21));

sub_id = [id_OAa21;id_OAb21];

n_a21 = length(id_OAa21);
n_b21 = length(id_OAb21);

des = [[ones(n_a21,1);zeros(n_b21,1)],sub_info(sub_id,3),sub_info(sub_id,4)];
% gradient range
load([data_dir,'emb_range.mat']);
stat_range = zeros(1,7);
    stat_range(1) = mean(emb_range(id_OAa21,1));
    stat_range(2) = std(emb_range(id_OAa21,1));
    stat_range(3) = mean(emb_range(id_OAb21,1));
    stat_range(4) = std(emb_range(id_OAb21,1));    
    stat_result = regstats(emb_range(sub_id,1),des,'linear',{'tstat','r'});
    stat_range(5) = stat_result.tstat.t(2);
    stat_range(6) = stat_result.tstat.t(2) * sqrt(1/n_a21 + 1/n_b21);
    stat_range(7) = stat_result.tstat.pval(2);
disp(stat_range)

% gradient variance
load([data_dir,'emb_std.mat']);
stat_std = zeros(1,7);
    stat_std(1,1) = mean(emb_std(id_OAa21,1));
    stat_std(1,2) = std(emb_std(id_OAa21,1));
    stat_std(1,3) = mean(emb_std(id_OAb21,1));
    stat_std(1,4) = std(emb_std(id_OAb21,1));    
    stat_result = regstats(emb_std(sub_id,1),des,'linear',{'tstat','r'});
    stat_std(1,5) = stat_result.tstat.t(2);
    stat_std(1,6) = stat_result.tstat.t(2) * sqrt(1/n_a21 + 1/n_b21);
    stat_std(1,7) = stat_result.tstat.pval(2);
disp(stat_std)

%% Permutation test
% select subject
id_OAa21 = find(sub_info(:,1)==2&sub_info(:,7)>21);
id_OAb21 = find(sub_info(:,1)==2&((sub_info(:,7)<=21&sub_info(:,7)~=-1)|sub_info(:,3)<=21));

sub_id = [id_OAa21;id_OAb21];

n_a21 = length(id_OAa21);
n_b21 = length(id_OAb21);

des = [[ones(n_a21,1);zeros(n_b21,1)],sub_info(sub_id,3),sub_info(sub_id,4)];
% gradient range
load([data_dir,'emb_range.mat']);
stat_range = zeros(2,1);
    stat_result = regstats(emb_range(sub_id,1),des,'linear',{'tstat','r'});
    stat_range(1) = stat_result.tstat.t(2);
    t_rand = zeros(10000,1);
    for j = 1:10000
        des_rand = des;
        des_rand(:,1) = des(randperm(size(des,1)),1);
        stat_result = regstats(emb_range(sub_id,1),des_rand,'linear',{'tstat','r'});
        t_rand(j) = stat_result.tstat.t(2);
    end
    if stat_range(1)>0
        stat_range(2) = length(find(t_rand>stat_range(1)))/10000;
    else
        stat_range(2) = length(find(t_rand<stat_range(1)))/10000;
    end
disp(stat_range)

% gradient variance
load([data_dir,'emb_std.mat']);
stat_std = zeros(2,1);
    stat_result = regstats(emb_std(sub_id,1),des,'linear',{'tstat','r'});
    stat_std(1) = stat_result.tstat.t(2);
    t_rand = zeros(10000,1);
    for j = 1:10000
        des_rand = des;
        des_rand(:,1) = des(randperm(size(des,1)),1);
        stat_result = regstats(emb_std(sub_id,1),des_rand,'linear',{'tstat','r'});
        t_rand(j) = stat_result.tstat.t(2);
    end
    if stat_std(1)>0
        stat_std(2) = length(find(t_rand>stat_std(1)))/10000;
    else
        stat_std(2) = length(find(t_rand<stat_std(1)))/10000;
    end
disp(stat_std)