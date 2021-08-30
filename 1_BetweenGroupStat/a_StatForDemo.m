%% Demographic Statistic
% number
for i = 1:10
    disp(length(find(sub_info(sub_info(:,2)==i,1)==2)));
    disp(length(find(sub_info(sub_info(:,2)==i,1)==1)));
end
disp(length(find(sub_info(:,1)==2)));
disp(length(find(sub_info(:,1)==1)));

% age
stat = zeros(11,6);
for i = 1:11
    if i ~= 11
        stat(i,1) = mean(sub_info(sub_info(:,2)==i&sub_info(:,1)==2,3));
        stat(i,2) = std(sub_info(sub_info(:,2)==i&sub_info(:,1)==2,3));
        stat(i,3) = mean(sub_info(sub_info(:,2)==i&sub_info(:,1)==1,3));
        stat(i,4) = std(sub_info(sub_info(:,2)==i&sub_info(:,1)==1,3));
        [H,stat(i,6),CI,STATS] = ttest2(sub_info(sub_info(:,2)==i&sub_info(:,1)==2,3),sub_info(sub_info(:,2)==i&sub_info(:,1)==1,3));
        stat(i,5) = STATS.tstat;
    else
        stat(i,1) = mean(sub_info(sub_info(:,1)==2,3));
        stat(i,2) = std(sub_info(sub_info(:,1)==2,3));
        stat(i,3) = mean(sub_info(sub_info(:,1)==1,3));
        stat(i,4) = std(sub_info(sub_info(:,1)==1,3));
        [H,stat(i,6),CI,STATS] = ttest2(sub_info(sub_info(:,1)==2,3),sub_info(sub_info(:,1)==1,3));
        stat(i,5) = STATS.tstat;
    end
end
disp(stat);

% sex
stat = zeros(11,6);
for i = 1:11
    if i ~= 11
        stat(i,1) = length(find(sub_info(sub_info(:,2)==i&sub_info(:,1)==2,4)==1));
        stat(i,2) = length(find(sub_info(sub_info(:,2)==i&sub_info(:,1)==2,4)==-1));
        stat(i,3) = length(find(sub_info(sub_info(:,2)==i&sub_info(:,1)==1,4)==1));
        stat(i,4) = length(find(sub_info(sub_info(:,2)==i&sub_info(:,1)==1,4)==-1));        
        [stat(i,6),stat(i,5)] = chi2test([stat(i,1),stat(i,2);stat(i,3),stat(i,4)]);
    else
        stat(i,1) = length(find(sub_info(sub_info(:,1)==2,4)==1));
        stat(i,2) = length(find(sub_info(sub_info(:,1)==2,4)==-1));
        stat(i,3) = length(find(sub_info(sub_info(:,1)==1,4)==1));
        stat(i,4) = length(find(sub_info(sub_info(:,1)==1,4)==-1));        
        [stat(i,6),stat(i,5)] = chi2test([stat(i,1),stat(i,2);stat(i,3),stat(i,4)]);
    end
end
disp(stat);

% edu
stat = zeros(11,6);
for i = 1:11
    if i ~= 11
        id_mdd = sub_info(:,2)==i&sub_info(:,1)==2&sub_info(:,5)~=-1;
        id_hc = sub_info(:,2)==i&sub_info(:,1)==1&sub_info(:,5)~=-1;
    else
        id_mdd = sub_info(:,1)==2&sub_info(:,5)~=-1;
        id_hc = sub_info(:,1)==1&sub_info(:,5)~=-1;
    end
        stat(i,1) = mean(sub_info(id_mdd,5));
        stat(i,2) = std(sub_info(id_mdd,5));
        stat(i,3) = mean(sub_info(id_hc,5));
        stat(i,4) = std(sub_info(id_hc,5));
        [H,stat(i,6),CI,STATS] = ttest2(sub_info(id_mdd,5),sub_info(id_hc,5));
        stat(i,5) = STATS.tstat;    
end
disp(stat);

% duration
stat = zeros(11,6);
feature_id = 6;
for i = 1:11
    if i ~= 11
        id_mdd = sub_info(:,2)==i&sub_info(:,1)==2&sub_info(:,feature_id)~=-1;
        id_hc = sub_info(:,2)==i&sub_info(:,1)==1&sub_info(:,feature_id)~=-1;
    else
        id_mdd = sub_info(:,1)==2&sub_info(:,feature_id)~=-1;
        id_hc = sub_info(:,1)==1&sub_info(:,feature_id)~=-1;
    end
        stat(i,1) = mean(sub_info(id_mdd,feature_id));
        stat(i,2) = std(sub_info(id_mdd,feature_id));
        stat(i,3) = mean(sub_info(id_hc,feature_id));
        stat(i,4) = std(sub_info(id_hc,feature_id));
        [H,stat(i,6),CI,STATS] = ttest2(sub_info(id_mdd,feature_id),sub_info(id_hc,feature_id));
        stat(i,5) = STATS.tstat;    
end
disp(stat);
    
% first episode
stat = zeros(11,6);
feature_id = 8;
for i = 1:11
    if i ~= 11
        id_mdd = sub_info(:,2)==i&sub_info(:,1)==2&sub_info(:,feature_id)~=-1;
    else
        id_mdd = sub_info(:,1)==2&sub_info(:,feature_id)~=-1;
    end
        stat(i,1) = length(find(sub_info(id_mdd,feature_id)==1));
        stat(i,2) = length(find(sub_info(id_mdd,feature_id)==0));

end
disp(stat);


%medicated
stat = zeros(11,6);
feature_id = 10;
for i = 1:11
    if i ~= 11
        id_mdd = sub_info(:,2)==i&sub_info(:,1)==2&sub_info(:,feature_id)~=-1;
    else
        id_mdd = sub_info(:,1)==2&sub_info(:,feature_id)~=-1;
    end
        stat(i,1) = length(find(sub_info(id_mdd,feature_id)==1));
        stat(i,2) = length(find(sub_info(id_mdd,feature_id)==0));    
end
disp(stat);

% hdrs
stat = zeros(11,6);
feature_id = 11;
for i = 1:11
    if i ~= 11
        id_mdd = sub_info(:,2)==i&sub_info(:,1)==2&sub_info(:,feature_id)~=-1;
    else
        id_mdd = sub_info(:,1)==2&sub_info(:,feature_id)~=-1;
    end
        stat(i,1) = mean(sub_info(id_mdd,feature_id));
        stat(i,2) = std(sub_info(id_mdd,feature_id));

end
disp(stat);


% mFD
stat = zeros(11,6);
feature_id = 12;
for i = 1:11
    if i ~= 11
        id_mdd = sub_info(:,2)==i&sub_info(:,1)==2&sub_info(:,feature_id)~=-1;
        id_hc = sub_info(:,2)==i&sub_info(:,1)==1&sub_info(:,feature_id)~=-1;
    else
        id_mdd = sub_info(:,1)==2&sub_info(:,feature_id)~=-1;
        id_hc = sub_info(:,1)==1&sub_info(:,feature_id)~=-1;
    end
        stat(i,1) = mean(sub_info(id_mdd,feature_id));
        stat(i,2) = std(sub_info(id_mdd,feature_id));
        stat(i,3) = mean(sub_info(id_hc,feature_id));
        stat(i,4) = std(sub_info(id_hc,feature_id));
        [H,stat(i,6),CI,STATS] = ttest2(sub_info(id_mdd,feature_id),sub_info(id_hc,feature_id));
        stat(i,5) = STATS.tstat;    
end
disp(stat);

% ANOVA for HDRS
feature_id = 11;
id_mdd = sub_info(:,1)==2&sub_info(:,feature_id)~=-1;
[p,anovatab,stats] = anova1(sub_info(id_mdd,feature_id),sub_info(id_mdd,2));
eta_squared = anovatab{2,2} / anovatab{4,2};
c = multcompare(stats,'CType','bonferroni');