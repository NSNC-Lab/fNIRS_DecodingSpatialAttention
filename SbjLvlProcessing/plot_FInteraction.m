% to be run after preprocessFNIRS06_FInteraction_MultiOnly.m

function plot_FInteraction(s,numClasses,fefOnlyOp,saveOp)

fs = 50;
timeLen = [0.1 0.2 0.5 1 1.5 2 3 4 5].*fs;

nrept = 10;
kFold = 5;

sbjNum = s.name;

processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];

if numClasses == 2
    if fefOnlyOp
        fileName = [processedDataDir filesep 'performance_3FLO_FEF_LR.mat'];
    else
        fileName = [processedDataDir filesep 'performance_3FLO_FEF_LR.mat'];
    end
else
    if fefOnlyOp
        fileName = [processedDataDir filesep 'performance_3FLO_FEF.mat'];
    else
        fileName = [processedDataDir filesep 'performance_3FLO_FEF.mat'];
    end
end

load(fileName,'performanceLDALedoitHbT_All','performanceLDALedoitHbT_3FLO',...
    'performanceLDALedoitHbT_IPS','performanceLDALedoitHbT_IPS_3FLO');

performanceLDALedoitHbTAll_Mean = mean(performanceLDALedoitHbT_All,1);
performanceLDALedoitHbT_3FLO_Mean = mean(performanceLDALedoitHbT_3FLO,1);
performanceLDALedoitHbTIPS_Mean = mean(performanceLDALedoitHbT_IPS,1);
performanceLDALedoitHbTIPS_3FLO_Mean = mean(performanceLDALedoitHbT_IPS_3FLO,1);

diff_FEF_HbT = performanceLDALedoitHbTAll_Mean - performanceLDALedoitHbT_3FLO_Mean;
diff_IPS_HbT = performanceLDALedoitHbTIPS_Mean - performanceLDALedoitHbTIPS_3FLO_Mean;

performanceLDALedoitHbTAll_CI = zeros(length(timeLen),2);
performanceLDALedoitHbTAll_3FLO_CI = zeros(length(timeLen),2);
performanceLDALedoitHbTIPS_CI = zeros(length(timeLen),2);
performanceLDALedoitHbTIPS_3FLO_CI = zeros(length(timeLen),2);

conf = 0.95;

for i = 1:length(timeLen)

    temp = performanceLDALedoitHbT_All(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDALedoitHbTAll_CI(i,:) = calcCITDist(nrept*kFold,se,conf);

    temp = performanceLDALedoitHbT_3FLO(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDALedoitHbTAll_3FLO_CI(i,:) = calcCITDist(nrept*kFold,se,conf);

    temp = performanceLDALedoitHbT_IPS(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDALedoitHbTIPS_CI(i,:) = calcCITDist(nrept*kFold,se,conf);

    temp = performanceLDALedoitHbT_IPS_3FLO(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDALedoitHbTIPS_3FLO_CI(i,:) = calcCITDist(nrept*kFold,se,conf);

end

cmap = loadDefaultColors(1);

figure('units','normalized','outerposition',[0 0 1 1]); hold on;
subplot(1,2,1);hold on;
errorbar(timeLen./fs,performanceLDALedoitHbTAll_Mean,performanceLDALedoitHbTAll_CI(:,2),'Color',cmap(1,:));hold on;
errorbar(timeLen./fs,performanceLDALedoitHbT_3FLO_Mean,performanceLDALedoitHbTAll_3FLO_CI(:,2),'Color',cmap(2,:));hold on;
errorbar(timeLen./fs,performanceLDALedoitHbTIPS_Mean,performanceLDALedoitHbTIPS_CI(:,2),'Color',cmap(3,:));hold on;
errorbar(timeLen./fs,performanceLDALedoitHbTIPS_3FLO_Mean,performanceLDALedoitHbTIPS_3FLO_CI(:,2),'Color',cmap(4,:));hold on;
legend({'All FEF','FEF 3FLO','All IPS','IPS 3FLO'},'Location','Southeast');
title(sbjNum);

subplot(1,2,2);hold on;
bar([1:1:9],[diff_FEF_HbT; diff_IPS_HbT]);
set(gca,'XTick',(1:1:9));
set(gca, 'XTickLabel', {'0.1s','0.2s','0.5s','1s','1.5s','2s','3s','4s','5s'});
title('All - 3 Features Left Out');
legend({'FEF','IPS'});

figSaveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' ...
        num2str(sbjNum) '\Figures\Classification'];
fn = sprintf('Performance_FInteraction_3FLO');

if saveOp
    
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
    
end

end