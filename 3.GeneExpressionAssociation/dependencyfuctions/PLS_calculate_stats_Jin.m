function PLS_calculate_stats_Jin(response_var_file, predictor_var_file, output_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the PLS calculate stats function with the following arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% response_var_file ------ full path to the PLS_MRI_response_vars.csv file
%%%                           that is created by the NSPN_CorticalMyelination 
%%%                           wrapper script
%%% predictor_var_file ----- full path to the PLS_gene_predictor_vars.csv file
%%%                           that is provided as raw data
%%% output_dir ------------- where to save the PLS_stats file (for PLS1 and PLS2 together)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Re-run PLS to get explained variance and associated stats')

%import response variables
% importdata(response_var_file);

%and store the response variables in matrix Y
MRIdata=response_var_file;

%import predictor variables
% importdata(predictor_var_file);
GENEdata=predictor_var_file;
geneindex=1:size(GENEdata,2);
load unique_gene.mat

%DO PLS in 2 dimensions (with 2 components) 
Y=zscore(MRIdata);
dim=10;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(GENEdata,Y,dim,'CV',dim);
 temp=cumsum(100*PCTVAR(2,1:dim));
 Rsquared = temp(dim);

%align PLS components with desired direction%
[R1,p1]=corr([XS(:,1),XS(:,2),XS(:,3)],MRIdata(:,1:3));
if R1(1,1)<0
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0
    XS(:,2)=-1*XS(:,2);
end
if R1(3,1)<0
    XS(:,3)=-1*XS(:,3);
end
%calculate correlations of PLS components with MRI variables
[R1,p1]=corr(XS(:,1),MRIdata);
[R2,p2]=corr(XS(:,2),MRIdata);
[R3,p3]=corr(XS(:,3),MRIdata);
a=[R1',p1',R2',p2'];

%assess significance of PLS result
for j=1:1000
    j
    order=randperm(size(Y,1));
    Yp=Y(order,:);
    [XLr,YLr,XSr,YSr,BETAr,PCTVARr,MSEr,statsr]=plsregress(GENEdata,Yp,dim);
     PCTVARrand(j,:)=PCTVARr(2,:);
     temp=cumsum(100*PCTVARr(2,1:dim));
     Rsq(j) = temp(dim);
end
for l=1:dim
p_single(l)=length(find(PCTVARrand(:,l)>=PCTVAR(2,l)))/j
end
p_cum=length(find(Rsq>=Rsquared))/j;

% plot histogram
% hist(Rsq,10)
% hold on
% plot(Rsquared,20,'.r','MarkerSize',15)
% set(gca,'Fontsize',14)
% xlabel('R squared','FontSize',14);
% ylabel('Permuted runs','FontSize',14);
% title('p<0.0001')

%save stats
myStats=[PCTVAR; p_single];
csvwrite(fullfile(output_dir,'PLS_stats.csv'),myStats);


