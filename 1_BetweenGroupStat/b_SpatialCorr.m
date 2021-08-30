%% spatial correlation 
%% define work dir
data_dir = 'D:\Data\DIDA-MDD\gradient_analysis\analysis2\';
script_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\3.Scripts\';

%% real correlation
load([data_dir,'gradient_emb_correct(1-3).mat']);
mean_gradient_mdd = zeros(3,18933);
mean_gradient_hc = zeros(3,18933);
for i = 1:3
    mean_gradient_mdd(i,:) = mean(gradient_emb_correct{i}(:,sub_info(:,1)==2),2);
    mean_gradient_hc(i,:) = mean(gradient_emb_correct{i}(:,sub_info(:,1)==1),2);
end
real_spatial_r = zeros(3,1);
for i = 1:3
    real_spatial_r(i) = corr(mean_gradient_hc(i,:)',mean_gradient_mdd(i,:)','type','spearman');
    dotcolor = [115 130 184]/255;
    linecolor = [0 0 0];
    ylable1 = 'Gradient scores of MDD';
    xlable1 = 'Gradient scores of controls';
    [xData, yData] = prepareCurveData(mean_gradient_hc(i,:), mean_gradient_mdd(i,:));
    ft = fittype( 'poly1' );
    opts = fitoptions( ft );
    opts.Lower = [-Inf -Inf];
    opts.Upper = [Inf Inf];
    [fitresult, gof] = fit( xData, yData, ft, opts );
    close all
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
    t1 = text(-0.06,-0.1,'{\it\rho} = 0.999, {\itP} < 0.00001','FontName','Arial','FontSize',6);
    set(gca,'XLim',[-0.12 0.08],'XTick',-0.12:0.05:0.08);
    set(gca,'YLim',[-0.12 0.08],'YTick',-0.12:0.05:0.08);
    set(gca,'LineWidth',0.5);
    set(gca,'FontName','Arial','FontSize',7);
    set(gca,'box','off');
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'Paperposition',[0 0 5.93 5.5]);
    grid off
    print(gcf,[figure_dir,'SpatialCorrG',num2str(i),'.tif'],'-dtiff','-r1000')
end

mean_g1_hc = mean_gradient_hc(1,:);
mkdir([data_dir,'\surrogate_maps_g1_hc']);
save([data_dir,'\surrogate_maps_g1_hc\mean_g1_hc.txt'],'mean_g1_hc','-ascii');

mean_g2_hc = mean_gradient_hc(2,:);
mkdir([data_dir,'\surrogate_maps_g2_hc']);
save([data_dir,'\surrogate_maps_g2_hc\mean_g2_hc.txt'],'mean_g2_hc','-ascii');

mean_g3_hc = mean_gradient_hc(3,:);
mkdir([data_dir,'\surrogate_maps_g3_hc']);
save([data_dir,'\surrogate_maps_g3_hc\mean_g3_hc.txt'],'mean_g3_hc','-ascii');

%% run python code of Brainsmash to generate surrogate maps for HC mean gradient

%% spatial corr p permutation
surrogate_r = zeros(3,10000);
p = ones(3,1);
for i = 1:3
    load([data_dir,'surrogate_maps_g',num2str(i),'_hc\surrogate_maps_g',num2str(i),'_hc_resample.mat']);
    surrogate_r(i,:) = corr(surrogate_maps',mean_gradient_mdd(i,:)','type','spearman');
    p(i) = length(find(surrogate_r(i,:)>real_spatial_r(i)))/10000;
end
    
    
    
    



