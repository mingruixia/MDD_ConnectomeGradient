%% load displacement in gradient space as features
load('SVR_Subjects_Data.mat');
load('SVR_Subjects_Scores.mat');
load('SVR_Covariates.mat')

%% perform SVR analysis
FoldQuantity = 10;
Pre_Method = 'Normalize';
C_Range = power(2, -5:10);
Weight_Flag = 1;
Permutation_Flag = 0;
ResultantFolder = 'D:\Data\DIDA-MDD\gradient_analysis\analysis\SVR_hdrs_dis';

Prediction = SVR_NFolds_Sort_CSelect(Subjects_Data, Subjects_Scores, Covariates, FoldQuantity, Pre_Method, C_Range, Weight_Flag, Permutation_Flag, ResultantFolder);

%% permutation, this was run on a cluster using PSOM
addpath('/home/mxia/matlab/libsvm-3.24/linux/');
addpath('/home/mxia/matlab/surfstat');
addpath('/home/mxia/matlab/Pattern_Regression-master/SVR');

load('/HeLabData2/mxia/DIDA-MDD_grediant/analysis/SVR_Subjects_Data.mat');
load('/HeLabData2/mxia/DIDA-MDD_grediant/analysis/SVR_Subjects_Scores.mat');
load('/HeLabData2/mxia/DIDA-MDD_grediant/analysis/SVR_Covariates.mat')

Opt.mode = 'qsub';
Opt.max_queued = 300;
Opt.flag_pause = false;
Opt.path_logs = '/HeLabData2/mxia/DIDA-MDD_grediant/analysis/psom_log';

for n = 1:10000
    Job_name = ['g',num2str(n)];
    pipeline.(Job_name).command = 'Prediction = SVR_NFolds_Sort_CSelect(opt.Subjects_Data, opt.Subjects_Scores, opt.Covariates, opt.FoldQuantity, opt.Pre_Method, opt.C_Range, opt.Weight_Flag, opt.Permutation_Flag, opt.ResultantFolder);';
    num = ['0000',num2str(n)];
    pipeline.(Job_name).opt.Subjects_Data = Subjects_Data;
    pipeline.(Job_name).opt.Subjects_Scores = Subjects_Scores;
    pipeline.(Job_name).opt.Covariates = Covariates;
    pipeline.(Job_name).opt.FoldQuantity = 10;
    pipeline.(Job_name).opt.Pre_Method = 'Normalize';
    pipeline.(Job_name).opt.C_Range = power(2, -5:10);
    pipeline.(Job_name).opt.Weight_Flag = 1;
    pipeline.(Job_name).opt.Permutation_Flag = 1;
    pipeline.(Job_name).opt.ResultantFolder = ['/HeLabData2/mxia/DIDA-MDD_grediant/analysis/SVR_rand/',num(end-4:end)];
end

psom_run_pipeline(pipeline,Opt);
