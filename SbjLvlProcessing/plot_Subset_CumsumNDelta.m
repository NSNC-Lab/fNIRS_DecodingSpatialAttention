% to be run after preprocessFNIRS06_FInteraction_MultiOnly.m

function plot_Subset_CumsumNDelta(s,numClasses,fefOnlyOp,saveOp)

fs = 50;
timeLen = [0.06 0.1 0.2 0.5 1 1.5 2 3 4 5].*fs;
deltaLen = [5 11 25 49];

nrept = 10;
kFold = 5;

sbjNum = s.name;

processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];

if numClasses == 2
    if fefOnlyOp
        fileName = [processedDataDir filesep 'performance_LOFO_Sbst_FEF_LR.mat'];
    else
        fileName = [processedDataDir filesep 'performance_LOFO_Sbst_LR.mat'];
    end
end

load(fileName,'performanceLDALedoitHbT_All','performanceLDALedoitHbT_Final',...
    'chnIdxSelect','performanceLDALedoitHbT_Delta1','performanceLDALedoitHbT_Delta');

performanceLDALedoitHbT_All_Mean = mean(performanceLDALedoitHbT_All,1);
performanceLDALedoitHbT_Subset_Mean = mean(performanceLDALedoitHbT_Final,1);
performanceLDALedoitHbT_Delta_Mean = mean(performanceLDALedoitHbT_Delta,1);
performanceLDALedoitHbT_Delta1_Mean = mean(performanceLDALedoitHbT_Delta1,1);

% diff_FEF_HbT = performanceLDALedoitHbTAll_Mean - performanceLDALedoitHbT_Delta_Mean;
% diff_IPS_HbT = performanceLDALedoitHbTIPS_Mean - performanceLDALedoitHbTIPS_3FLO_Mean;

performanceLDALedoitHbT_All_CI = zeros(length(timeLen),2);
performanceLDALedoitHbT_Subset_CI = zeros(length(timeLen),2);
performanceLDALedoitHbT_Delta_CI = zeros(length(deltaLen),2);
performanceLDALedoitHbT_Delta1_CI = zeros(length(deltaLen),2);

conf = 0.95;



for i = 1:length(timeLen)

    temp = performanceLDALedoitHbT_All(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDALedoitHbT_All_CI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceLDALedoitHbT_Final(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDALedoitHbT_Subset_CI(i,:) = calcCITDist(nrept*kFold,se,conf);

end

for i = 1:length(deltaLen)
    
    temp = performanceLDALedoitHbT_Delta(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDALedoitHbT_Delta_CI(i,:) = calcCITDist(nrept*kFold,se,conf);

    temp = performanceLDALedoitHbT_Delta1(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDALedoitHbT_Delta1_CI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
end

cmap = loadDefaultColors(1);

figure('units','normalized','outerposition',[0 0 1 1]); hold on;
%subplot(1,2,1);hold on;
errorbar(timeLen./fs,performanceLDALedoitHbT_All_Mean,performanceLDALedoitHbT_All_CI(:,2),'Color',cmap(1,:));hold on;
errorbar(timeLen./fs,performanceLDALedoitHbT_Subset_Mean,performanceLDALedoitHbT_Subset_CI(:,2),'Color',cmap(2,:));hold on;
errorbar(deltaLen./fs,performanceLDALedoitHbT_Delta_Mean,performanceLDALedoitHbT_Delta_CI(:,2),'Color',cmap(3,:));hold on;
errorbar(deltaLen./fs,performanceLDALedoitHbT_Delta1_Mean,performanceLDALedoitHbT_Delta1_CI(:,2),'Color',cmap(4,:));hold on;
legend({'All','Subset Sum','Subset Delta 0','Subset Delta 1'},'Location','Southeast');
%legend({'Subset Sum','Subset Delta 0','Subset Delta 1'},'Location','Southeast');
title(sbjNum);

figSaveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' ...
        num2str(sbjNum) '\Figures\Classification'];
fn = sprintf('Performance_SubsetNDelta');

if saveOp
    
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
    
end

end