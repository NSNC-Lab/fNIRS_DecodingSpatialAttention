function plot_MultiOnly_DiffStartT_AllChn_CompClassifiers_Final(sbjNum,numClasses,saveOp)

processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];

if numClasses == 2
    savePerfFN = 'performanceLinearDiscriminantUpdated_LR.mat';
else
    savePerfFN = 'performanceLinearDiscriminantUpdated.mat';
end
load([processedDataDir filesep savePerfFN]);

conf = 0.95;

performanceLDALedoitHbOCI = zeros(length(timePt),2);
performanceLDALedoitHbRCI = zeros(length(timePt),2);
performanceLDALedoitHbTCI = zeros(length(timePt),2);

performanceLDACERNNHbOCI = zeros(length(timePt),2);
performanceLDACERNNHbRCI = zeros(length(timePt),2);
performanceLDACERNNHbTCI = zeros(length(timePt),2);

performanceLogRegHbOCI = zeros(length(timePt),2);
performanceLogRegHbRCI = zeros(length(timePt),2);
performanceLogRegHbTCI = zeros(length(timePt),2);

performanceSVMHbOCI = zeros(length(timePt),2);
performanceSVMHbRCI = zeros(length(timePt),2);
performanceSVMHbTCI = zeros(length(timePt),2);

performanceBaggingHbOCI = zeros(length(timePt),2);
performanceBaggingHbRCI = zeros(length(timePt),2);
performanceBaggingHbTCI = zeros(length(timePt),2);

performanceBoostedHbOCI = zeros(length(timePt),2);
performanceBoostedHbRCI = zeros(length(timePt),2);
performanceBoostedHbTCI = zeros(length(timePt),2);

for i = 2:length(timePt)
    
    temp = performanceLDALedoitHbOFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDALedoitHbOMean(i) = mean(temp(:));
    performanceLDALedoitHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceLDALedoitHbRFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDALedoitHbRMean(i) = mean(temp(:));
    performanceLDALedoitHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceLDALedoitHbTFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDALedoitHbTMean(i) = mean(temp(:));
    performanceLDALedoitHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % CERNN
    temp = performanceLDACERNNHbOFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDACERNNHbOMean(i) = mean(temp(:));
    performanceLDACERNNHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceLDACERNNHbRFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDACERNNHbRMean(i) = mean(temp(:));
    performanceLDACERNNHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceLDACERNNHbTFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDACERNNHbTMean(i) = mean(temp(:));
    performanceLDACERNNHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % Log
    temp = performanceLogRegHbOFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLogRegHbOMean(i) = mean(temp(:));
    performanceLogRegHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceLogRegHbRFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLogRegHbRMean(i) = mean(temp(:));
    performanceLogRegHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceLogRegHbTFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLogRegHbTMean(i) = mean(temp(:));
    performanceLogRegHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % SVM
    temp = performanceSVMHbOFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceSVMHbOMean(i) = mean(temp(:));
    performanceSVMHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceSVMHbRFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceSVMHbRMean(i) = mean(temp(:));
    performanceSVMHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceSVMHbTFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceSVMHbTMean(i) = mean(temp(:));
    performanceSVMHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % Bagging
    temp = performanceBaggingHbOFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceBaggingHbOMean(i) = mean(temp(:));
    performanceBaggingHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceBaggingHbRFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceBaggingHbRMean(i) = mean(temp(:));
    performanceBaggingHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceBaggingHbTFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceBaggingHbTMean(i) = mean(temp(:));
    performanceBaggingHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % Boosting
    temp = performanceBoostedHbOFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceBoostedHbOMean(i) = mean(temp(:));
    performanceBoostedHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceBoostedHbRFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceBoostedHbRMean(i) = mean(temp(:));
    performanceBoostedHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceBoostedHbTFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceBoostedHbTMean(i) = mean(temp(:));
    performanceBoostedHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
end

numClassifiers = 6;
cmap = jet(numClassifiers);

subplot(1,3,1);hold on;
% for j = 1:size(performanceArrMultiHbO,1)
%     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
%     plot(timePt./fs-2,performanceArrMultiHbO(j,:),'Color',cmap(j,:));
% end
% plot(timePt./fs-2,performanceArrMultiHbO_95,'Color',cmap(1,:));
% plot(timePt./fs-2,performanceArrMultiHbO_85,'Color',cmap(2,:));
% plot(timePt./fs-2,performanceArrMultiHbO_75,'Color',cmap(3,:));
% plot(timePt./fs-2,performanceArrMultiHbO_65,'Color',cmap(4,:));
% plot(timePt./fs-2,performanceArrMultiHbO_55,'Color',cmap(5,:));
% plot(timePt./fs-2,performanceLDALedoitHbO,'Color',cmap(6,:));
% plot(timePt./fs-2,performanceLDACERNNHbO,'Color',cmap(7,:));
% % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
% % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
% plot(timePt./fs-2,performanceLogRegHbO,'Color',cmap(8,:));
% plot(timePt./fs-2,performanceSVMHbO,'Color',cmap(9,:));

errorbar(timePt./fs-2,performanceLDALedoitHbOMean,performanceLDALedoitHbOCI(:,2),'Color',cmap(1,:));
errorbar(timePt./fs-2,performanceLDACERNNHbOMean,performanceLDACERNNHbOCI(:,2),'Color',cmap(2,:));
% plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
% plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
errorbar(timePt./fs-2,performanceLogRegHbOMean,performanceLogRegHbOCI(:,2),'Color',cmap(3,:));
errorbar(timePt./fs-2,performanceSVMHbOMean,performanceSVMHbOCI(:,2),'Color',cmap(4,:));
errorbar(timePt./fs-2,performanceBaggingHbOMean,performanceBaggingHbOCI(:,2),'Color',cmap(5,:));
errorbar(timePt./fs-2,performanceBoostedHbOMean,performanceBoostedHbOCI(:,2),'Color',cmap(6,:));

ylim(yLimAxis);
xlim([0 7]);
title(sprintf('Sbj %s: Δ[HbO] Multi',num2str(sbjNum)));
% legend({'PCA 95%','PCA 85%','PCA 75%','PCA 65%','PCA 55%','LDAShrink',...
%     'QDAShrink','LDAGD','LogReg','SVM'});
% legend({'PCA 95%','PCA 85%','PCA 75%','PCA 65%','PCA 55%','LDAShrink',...
%     'LDA CERNN','LogReg','SVM'});
legend({'LDAShrink','LDA CERNN','LogReg','SVM','Bagging','Boosted'});
xlabel('Time [s]');ylabel('Accuracy');
hold off;

subplot(1,3,2);hold on;
% for j = 1:size(performanceArrMultiHbR,1)
%     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
%     plot(timePt./fs-2,performanceArrMultiHbR(j,:),'Color',cmap(j,:));
% end
% plot(timePt./fs-2,performanceArrMultiHbR_95,'Color',cmap(1,:));
% plot(timePt./fs-2,performanceArrMultiHbR_85,'Color',cmap(2,:));
% plot(timePt./fs-2,performanceArrMultiHbR_75,'Color',cmap(3,:));
% plot(timePt./fs-2,performanceArrMultiHbR_65,'Color',cmap(4,:));
% plot(timePt./fs-2,performanceArrMultiHbR_55,'Color',cmap(5,:));
% plot(timePt./fs-2,performanceLDALedoitHbR,'Color',cmap(6,:));
% plot(timePt./fs-2,performanceLDACERNNHbR,'Color',cmap(7,:));
% % plot(timePt./fs-2,performanceQDALedoitHbR,'Color',cmap(7,:));
% % plot(timePt./fs-2,performanceLDAGDHbR,'Color',cmap(8,:));
% plot(timePt./fs-2,performanceLogRegHbR,'Color',cmap(8,:));
% plot(timePt./fs-2,performanceSVMHbR,'Color',cmap(9,:));

errorbar(timePt./fs-2,performanceLDALedoitHbRMean,performanceLDALedoitHbRCI(:,2),'Color',cmap(1,:));
errorbar(timePt./fs-2,performanceLDACERNNHbRMean,performanceLDACERNNHbRCI(:,2),'Color',cmap(2,:));
% plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
% plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
errorbar(timePt./fs-2,performanceLogRegHbRMean,performanceLogRegHbRCI(:,2),'Color',cmap(3,:));
errorbar(timePt./fs-2,performanceSVMHbRMean,performanceSVMHbRCI(:,2),'Color',cmap(4,:));
errorbar(timePt./fs-2,performanceBaggingHbRMean,performanceBaggingHbRCI(:,2),'Color',cmap(5,:));
errorbar(timePt./fs-2,performanceBoostedHbRMean,performanceBoostedHbRCI(:,2),'Color',cmap(6,:));

ylim(yLimAxis);
xlim([0 7]);
title(sprintf('Sbj %s: Δ[HbR] Multi',num2str(sbjNum)));
%legend('S1D1','S1D2','S1D3','S1D4');
xlabel('Time [s]');ylabel('Accuracy');
hold off;

subplot(1,3,3);hold on;
% for j = 1:size(performanceArrMultiHbT,1)
%     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
%     plot(timePt./fs-2,performanceArrMultiHbT(j,:),'Color',cmap(j,:));
% end
% plot(timePt./fs-2,performanceArrMultiHbT_95,'Color',cmap(1,:));
% plot(timePt./fs-2,performanceArrMultiHbT_85,'Color',cmap(2,:));
% plot(timePt./fs-2,performanceArrMultiHbT_75,'Color',cmap(3,:));
% plot(timePt./fs-2,performanceArrMultiHbT_65,'Color',cmap(4,:));
% plot(timePt./fs-2,performanceArrMultiHbT_55,'Color',cmap(5,:));
% plot(timePt./fs-2,performanceLDALedoitHbT,'Color',cmap(6,:));
% plot(timePt./fs-2,performanceLDACERNNHbT,'Color',cmap(7,:));
% % plot(timePt./fs-2,performanceQDALedoitHbT,'Color',cmap(7,:));
% % plot(timePt./fs-2,performanceLDAGDHbT,'Color',cmap(8,:));
% plot(timePt./fs-2,performanceLogRegHbT,'Color',cmap(8,:));
% plot(timePt./fs-2,performanceSVMHbT,'Color',cmap(9,:));

errorbar(timePt./fs-2,performanceLDALedoitHbTMean,performanceLDALedoitHbTCI(:,2),'Color',cmap(1,:));
errorbar(timePt./fs-2,performanceLDACERNNHbTMean,performanceLDACERNNHbTCI(:,2),'Color',cmap(2,:));
% plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
% plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
errorbar(timePt./fs-2,performanceLogRegHbTMean,performanceLogRegHbTCI(:,2),'Color',cmap(3,:));
errorbar(timePt./fs-2,performanceSVMHbTMean,performanceSVMHbTCI(:,2),'Color',cmap(4,:));
errorbar(timePt./fs-2,performanceBaggingHbTMean,performanceBaggingHbTCI(:,2),'Color',cmap(5,:));
errorbar(timePt./fs-2,performanceBoostedHbTMean,performanceBoostedHbTCI(:,2),'Color',cmap(6,:));

ylim(yLimAxis);
xlim([0 7]);
title(sprintf('Sbj %s: Δ[HbT] Multi',num2str(sbjNum)));
%legend('S1D1','S1D2','S1D3','S1D4');
xlabel('Time [s]');ylabel('Accuracy');
hold off;

% subplot(1,4,4);hold on;
% for j = 1:size(performanceArrMultiHbT,1)
%     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
%     plot(timePt./fs-2,zeros(1,length(timePt./fs-2)),'Color',cmap(j,:));
% end
% %legend(chnName(1:30));
% annotation('textbox', [0.8, 0.05, 0.1, 0.1], 'String', ['Score: ' num2str(behScore) '%'],'FontSize',8);
% annotation('textbox', [0.71, 0.2, 0.1, 0.1],'String','calcCumSum\_MultiOnly\_DiffStartT\_AllChn\_AllBasis.m','FontSize',8);
% title(sprintf('All Channels. Diff Start T: %s-Class',num2str(numClasses)));
% hold off;

if numClasses == 2
    fn = sprintf('PerformanceCumsumVsTime_DiffStartT_LR_AllChns_DiffClassifiers');
else
    fn = sprintf('PerformanceCumsumVsTime_DiffStartT_AllChns_DiffClassifiers');
end

if saveOp == 1
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
end

end