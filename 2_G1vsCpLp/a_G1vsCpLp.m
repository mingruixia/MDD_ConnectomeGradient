%% Association with CpLp
%% define work dir
data_dir = 'D:\Data\DIDA-MDD\gradient_analysis\analysis2\';
script_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\3.Scripts\';
figure_dir = 'D:\SynologyDrive\ds1618+\SynologyDrive\Projects\2019_MDD_gradient\5.Figures\ForRevision1\';

%% extract Cp Lp
% load original data, calculated by using PAGANI
cd([data_dir,'CpLp\weighted']);
filelist_cp = dir('*cp.txt');
n_sub = length(filelist_cp);
cp_all = zeros(n_sub,1);
for i = 1:n_sub
    cp_all(i) = load(filelist_cp(i).name);
end

filelist_lp = dir('*lp.txt');
lp_all = zeros(n_sub,1);
for i = 1:n_sub
    lp_all(i) = load(filelist_lp(i).name);
end


%% Between-Group diff in Cp Lp
sub_info = xlsread([script_dir,'All_sub_info.xlsx']);
stat = zeros(2,7);
id_hc = sub_info(:,1)==1;
id_mdd = sub_info(:,1)==2;
des = [sub_info(:,1),sub_info(:,3:4)];
n_hc = length(find(id_hc));
n_mdd = length(find(id_mdd));
for i = 1:2
    stat(i,1) = mean(cp_lp_all(id_mdd,i));
    stat(i,2) = std(cp_lp_all(id_mdd,i));
    stat(i,3) = mean(cp_lp_all(id_hc,i));
    stat(i,4) = std(cp_lp_all(id_hc,i));    
    stat_result = regstats(cp_lp_all(:,i),des,'linear',{'tstat','r'});
    stat(i,5) = stat_result.tstat.t(2);
    stat(i,6) = stat_result.tstat.t(2) * sqrt(1/n_hc + 1/n_mdd);
    stat(i,7) = stat_result.tstat.pval(2);
end
disp(stat)

%% draw between-group figure for Cp
close all
sub_info = xlsread([script_dir,'All_sub_info.xlsx']);
id_hc = sub_info(:,1)==1;
id_mdd = sub_info(:,1)==2;
HC_cp = cp_lp_all(id_hc,1);
MDD_cp = cp_lp_all(id_mdd,1);
Data = cell(1,2);
Data{1} = MDD_cp;
Data{2} = HC_cp;
color_mdd = [226,115,134]./255;
color_hc = [0.7,0.7,0.7];
gretna_plot_violin(Data, {'MDD','HC'}, {'cp'}, 'meanstdfill');
h = gca;
h.Children(4).FaceColor = color_hc;
h.Children(5).FaceColor = color_mdd;
hold on
plot([1,2],[0.65,0.65],'-k','LineWidth',0.3);
set(gca,'YLim',[0.1,0.7],'YTick',[0.1,0.4 0.7]);
hold off
t1 = text(1.5,0.67,'***','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',8);
set(gca,'XTick',[1,2]);
set(gca,'XTickLabel',{'MDD','Controls'});
ylabel('Clustering coefficient');
legend off
set(gca,'LineWidth',0.5);
set(gca,'FontName','Arial','FontSize',6);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 3 3.5]);
print(gcf,[figure_dir,'Cp.tif'],'-dtiff','-r1000')

%% Corr G1 and Cp
stat = zeros(3,2);
[stat(1,1),stat(1,2)] = partialcorr(exprate_corrected(id_mdd,1),cp_lp_all(id_mdd,1),sub_info(id_mdd,3:4));
[stat(2,1),stat(2,2)] = partialcorr(emb_range(id_mdd,1),cp_lp_all(id_mdd,1),sub_info(id_mdd,3:4));
[stat(3,1),stat(3,2)] = partialcorr(emb_std(id_mdd,1),cp_lp_all(id_mdd,1),sub_info(id_mdd,3:4));
disp(stat);

% G1_exp vs Cp
close all
dotcolor = [115 130 184]/255;
linecolor = [0 0 0];
ylable1 = 'Clustering coefficient';
xlable1 = 'Gradient explained ratio';
x = gretna_glm(exprate_corrected(id_mdd,1),sub_info(id_mdd,3:4),'r');
y = gretna_glm(cp_lp_all(id_mdd,1),sub_info(id_mdd,3:4),'r');
[xData, yData] = prepareCurveData(x.r, y.r);
ft = fittype( 'poly1' );
opts = fitoptions( ft );
opts.Lower = [-Inf -Inf];
opts.Upper = [Inf Inf];
[fitresult, gof] = fit( xData, yData, ft, opts );
h=plot( fitresult, xData, yData);
set(h(1),'Marker','.','MarkerSize',3,'Color',dotcolor)
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
t1 = text(0.13,0.15,'{\itr} = 0.336, {\itP} < 0.001','FontName','Arial','FontSize',4);
set(gca,'YLim',[0.1,0.7]);
set(gca,'YTick',0.1:0.3:0.7);
set(gca,'XLim',[0,0.28]);
set(gca,'XTick',0:0.14:0.28);
set(gca,'LineWidth',0.5);
set(gca,'FontName','Arial','FontSize',6);
set(gca,'box','off');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[0 0 3.5 3.5]);
grid off
print(gcf,[figure_dir,'G1_exp_vs_Cp.tif'],'-dtiff','-r1000')

% G1_range vs Cp
close all
dotcolor = [115 130 184]/255;
linecolor = [0 0 0];
ylable1 = 'Clustering coefficient';
xlable1 = 'Gradient range';
x = gretna_glm(emb_range(id_mdd,1),sub_info(id_mdd,3:4),'r');
y = gretna_glm(cp_lp_all(id_mdd,1),sub_info(id_mdd,3:4),'r');
[xData, yData] = prepareCurveData(x.r, y.r);
ft = fittype( 'poly1' );
opts = fitoptions( ft );
opts.Lower = [-Inf -Inf];
opts.Upper = [Inf Inf];
[fitresult, gof] = fit( xData, yData, ft, opts );
h=plot( fitresult, xData, yData);
set(h(1),'Marker','.','MarkerSize',3,'Color',dotcolor)
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
t1 = text(0.20,0.15,'{\itr} = 0.257, {\itP} < 0.001','FontName','Arial','FontSize',4);
set(gca,'YLim',[0.1,0.7]);
set(gca,'YTick',0.1:0.3:0.7);
set(gca,'XLim',[0.1,0.36]);
set(gca,'XTick',0.1:0.13:0.36);
set(gca,'LineWidth',0.5);
set(gca,'FontName','Arial','FontSize',6);
set(gca,'box','off');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[0 0 3.5 3.5]);
grid off
print(gcf,[figure_dir,'G1_range_vs_Cp.tif'],'-dtiff','-r1000')

% G1_std vs Cp
close all
dotcolor = [115 130 184]/255;
linecolor = [0 0 0];
ylable1 = 'Clustering coefficient';
xlable1 = 'Gradient variance';
x = gretna_glm(emb_std(id_mdd,1),sub_info(id_mdd,3:4),'r');
y = gretna_glm(cp_lp_all(id_mdd,1),sub_info(id_mdd,3:4),'r');
[xData, yData] = prepareCurveData(x.r, y.r);
ft = fittype( 'poly1' );
opts = fitoptions( ft );
opts.Lower = [-Inf -Inf];
opts.Upper = [Inf Inf];
[fitresult, gof] = fit( xData, yData, ft, opts );
h=plot( fitresult, xData, yData);
set(h(1),'Marker','.','MarkerSize',3,'Color',dotcolor)
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
t1 = text(0.05,0.15,'{\itr} = 0.356, {\itP} < 0.001','FontName','Arial','FontSize',4);
set(gca,'YLim',[0.1,0.7]);
set(gca,'YTick',0.1:0.3:0.7);
set(gca,'XLim',[0.02,0.1]);
set(gca,'XTick',0.02:0.04:0.1);
set(gca,'LineWidth',0.5);
set(gca,'FontName','Arial','FontSize',6);
set(gca,'box','off');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[0 0 3.5 3.5]);
grid off
print(gcf,[figure_dir,'G1_std_vs_Cp.tif'],'-dtiff','-r1000')