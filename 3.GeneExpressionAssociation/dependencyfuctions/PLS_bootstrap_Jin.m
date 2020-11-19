function [PLS1_score] = PLS_bootstrap_Jin(response_var_file, predictor_var_file,gene_name, region_ind,output_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the PLS bootstrap function with the following arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% response_var_file ------ full path to the PLS_MRI_response_vars.csv file
%%%                           that is created by the NSPN_CorticalMyelination 
%%%                           wrapper script
%%% predictor_var_file ----- full path to the PLS_gene_predictor_vars.csv file
%%%                           that is provided as raw data
%%% output_dir ------------- where to save the PLS_geneWeights and PLS_ROIscores 
%%%                           files (for PLS1 and PLS2 separately)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Running PLS')

%import response variables
% importdata(response_var_file);

%unwrap and tidy MRI response variable names and store the response variables in matrix Y
MRIdata=response_var_file;

%import predictor variables
GENEdata=predictor_var_file;
geneindex=1:size(GENEdata,2);
% load /Users/jinjin/Documents/MDD_gradients/unique_gene.mat
genes = gene_name;

%number of bootstrap iterations
bootnum=10000;

%DO PLS in 2 dimensions (with 2 components) 
X=GENEdata;
Y=zscore(MRIdata);
dim=1;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

%store regions IDs and weights in descending order of weight for both
%components
[R1,p1]=corr([XS(:,1)],MRIdata);
% [R2,p2]=corr([XS(:,2)],MRIdata);
% [R3,p3]=corr([XS(:,3)],MRIdata);
% [R4,p4]=corr([XS(:,4)],MRIdata);
% [R5,p5]=corr([XS(:,5)],MRIdata);
%align PLS components with desired direction%
if R1(1,1)<0
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
% if R2(1,1)<0
%     stats.W(:,2)=-1*stats.W(:,2);
%     XS(:,2)=-1*XS(:,2);
% end
% if R3(1,1)<0
%     stats.W(:,3)=-1*stats.W(:,3);
%     XS(:,3)=-1*XS(:,3);
% end
% if R4(1,1)<0
%     stats.W(:,4)=-1*stats.W(:,4);
%     XS(:,4)=-1*XS(:,4);
% end
% if R5(1,1)<0
%     stats.W(:,5)=-1*stats.W(:,5);
%     XS(:,5)=-1*XS(:,5);
% end
[PLS1w,x1] = sort(stats.W(:,1),'descend');
PLS1ids=genes(x1);
geneindex1=geneindex(x1);
% [PLS2w,x2] = sort(stats.W(:,2),'descend');
% PLS2ids=genes(x2);
% geneindex2=geneindex(x2);
% [PLS3w,x3] = sort(stats.W(:,3),'descend');
% PLS3ids=genes(x3);
% geneindex3=geneindex(x3);
% [PLS4w,x4] = sort(stats.W(:,4),'descend');
% PLS4ids=genes(x4);
% geneindex4=geneindex(x4);
% [PLS5w,x5] = sort(stats.W(:,5),'descend');
% PLS5ids=genes(x5);
% geneindex5=geneindex(x5);

PLS1_ROIscores_360=zeros(360,1);
% PLS2_ROIscores_360=zeros(360,1);
% PLS3_ROIscores_360=zeros(360,1);
% PLS4_ROIscores_360=zeros(360,1);
% PLS5_ROIscores_360=zeros(360,1);
% load M:\dynamic_gene_final\geneExpression\sample_matching\region_ind.mat
PLS1_ROIscores_360(region_ind)=XS(:,1);
% PLS2_ROIscores_360(region_ind)=XS(:,2);
% PLS3_ROIscores_360(region_ind)=XS(:,3);
% PLS4_ROIscores_360(region_ind)=XS(:,4);
% PLS5_ROIscores_360(region_ind)=XS(:,5);
save([output_dir 'PLS1_ROIscore.mat'],'PLS1_ROIscores_360')
% save([output_dir 'PLS2_ROIscore.mat'],'PLS2_ROIscores_360')
% save([output_dir 'PLS3_ROIscore.mat'],'PLS3_ROIscores_360')
% save([output_dir 'PLS4_ROIscore.mat'],'PLS4_ROIscores_360')
% save([output_dir 'PLS5_ROIscore.mat'],'PLS5_ROIscores_360')

%print out results
csvwrite(fullfile(output_dir,'PLS1_ROIscores.csv'),XS(:,1));
% csvwrite(fullfile(output_dir,'PLS2_ROIscores.csv'),XS(:,2));
% csvwrite(fullfile(output_dir,'PLS3_ROIscores.csv'),XS(:,3));
% csvwrite(fullfile(output_dir,'PLS4_ROIscores.csv'),XS(:,4));
% csvwrite(fullfile(output_dir,'PLS5_ROIscores.csv'),XS(:,5));

PLS1_score=XS(:,1);
% PLS2_score=XS(:,2);
% PLS3_score=XS(:,3);
% PLS4_score=XS(:,4);
% PLS5_score=XS(:,5);

%define variables for storing the (ordered) weights from all bootstrap runs
PLS1weights=[];
% PLS2weights=[];
% PLS3weights=[];
% PLS4weights=[];
% PLS5weights=[];

%start bootstrap
disp('  Bootstrapping - could take a while')
for i=1:bootnum
    i
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
    PLS1weights=[PLS1weights,newW];%store (ordered) weights from this bootstrap run
    
%     temp=stats.W(:,2);%extract PLS2 weights
%     newW=temp(x2); %order the newly obtained weights the same way as initial PLS 
%     if corr(PLS2w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
%         newW=-1*newW;
%     end
%     PLS2weights=[PLS2weights,newW]; %store (ordered) weights from this bootstrap run    
%     
%     temp=stats.W(:,3);%extract PLS2 weights
%     newW=temp(x3); %order the newly obtained weights the same way as initial PLS 
%     if corr(PLS3w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
%         newW=-1*newW;
%     end
%     PLS3weights=[PLS3weights,newW]; %store (ordered) weights from this bootstrap run    
    
%     temp=stats.W(:,4);%extract PLS2 weights
%     newW=temp(x4); %order the newly obtained weights the same way as initial PLS 
%     if corr(PLS4w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
%         newW=-1*newW;
%     end
%     PLS4weights=[PLS4weights,newW]; %store (ordered) weights from this bootstrap run    
%     
%     temp=stats.W(:,5);%extract PLS2 weights
%     newW=temp(x5); %order the newly obtained weights the same way as initial PLS 
%     if corr(PLS5w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
%         newW=-1*newW;
%     end
%     PLS5weights=[PLS5weights,newW]; %store (ordered) weights from this bootstrap run    
    
end

%get standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');
% PLS2sw=std(PLS2weights');
% PLS3sw=std(PLS3weights');
% PLS4sw=std(PLS4weights');
% PLS5sw=std(PLS5weights');

%get bootstrap weights
temp1=PLS1w./PLS1sw';
% temp2=PLS2w./PLS2sw';
% temp3=PLS3w./PLS3sw';
% temp4=PLS4w./PLS4sw';
% temp5=PLS5w./PLS5sw';

%order bootstrap weights (Z) and names of regions
[Z1 ind1]=sort(temp1,'ascend');
PLS1=PLS1ids(ind1);
geneindex1=geneindex1(ind1);

% [Z2 ind2]=sort(temp2,'descend');
% PLS2=PLS2ids(ind2);
% geneindex2=geneindex2(ind2);
% 
% [Z3 ind3]=sort(temp3,'descend');
% PLS3=PLS3ids(ind3);
% geneindex3=geneindex3(ind3);

% [Z4 ind4]=sort(temp4,'descend');
% PLS4=PLS4ids(ind4);
% geneindex4=geneindex4(ind4);
% 
% [Z5 ind5]=sort(temp5,'descend');
% PLS5=PLS5ids(ind5);
% geneindex5=geneindex5(ind5);

%print out results
fid1 = fopen(fullfile(output_dir,'PLS1_geneWeights_inv.csv'),'w');
for i=1:length(genes)
  fprintf(fid1,'%s, %d, %f\n', PLS1{i}, geneindex1(i), Z1(i));
end
fclose(fid1);

% fid2 = fopen(fullfile(output_dir,'PLS2_geneWeights.csv'),'w');
% for i=1:length(genes)
%   fprintf(fid2,'%s, %d, %f\n', PLS2{i},geneindex2(i), Z2(i));
% end
% fclose(fid2);
% 
% fid3 = fopen(fullfile(output_dir,'PLS3_geneWeights.csv'),'w');
% for i=1:length(genes)
%   fprintf(fid3,'%s, %d, %f\n', PLS3{i},geneindex3(i), Z3(i));
% end
% fclose(fid3);

% fid4 = fopen(fullfile(output_dir,'PLS4_geneWeights.csv'),'w');
% for i=1:length(genes)
%   fprintf(fid4,'%s, %d, %f\n', PLS4{i},geneindex4(i), Z4(i));
% end
% fclose(fid4);
% 
% fid5 = fopen(fullfile(output_dir,'PLS5_geneWeights.csv'),'w');
% for i=1:length(genes)
%   fprintf(fid5,'%s, %d, %f\n', PLS5{i},geneindex5(i), Z5(i));
% end
% fclose(fid5);
