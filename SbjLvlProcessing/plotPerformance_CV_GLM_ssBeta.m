function plotPerformance_CV_GLM_ssBeta(sbjNum,numClasses,saveOp,rejTrOp)

processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
figSaveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' ...
    num2str(sbjNum) '\Figures\Classification'];

if ~exist(figSaveDir,'dir')
    mkdir(figSaveDir);
end

fs = 50;
%timePt = (0:0.25*fs:5*fs)+2*fs;
timePt = (-2*fs:0.25*fs:5*fs);
yLimAxis = [0 1];

if numClasses == 2
    if rejTrOp
        fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_LR_RejTr_SNR_1.5_OrigHomer3.mat'];
        %fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_LR_RejTr_SNR_1.5.mat'];
        %fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_LR_RejTr_SNR_1.5_OrigHomer3.mat'];
    else
        fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_LR.mat'];
    end
else
    if rejTrOp
        fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_RejTr_SNR_1.5.mat'];
    else
        fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta.mat'];
    end
end
% folds x time
load(fileName,'performanceLDALedoitHbO','performanceLDACERNNHbO',...
    'performanceLogRegHbO','performanceSVMHbO','performanceBaggingHbO',...
    'performanceBoostedHbO','performanceLDALedoitHbR','performanceLDACERNNHbR',...
    'performanceLogRegHbR','performanceSVMHbR','performanceBaggingHbR',...
    'performanceBoostedHbR','performanceLDALedoitHbT','performanceLDACERNNHbT',...
    'performanceLogRegHbT','performanceSVMHbT','performanceBaggingHbT',...
    'performanceBoostedHbT');

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

conf = 0.95;
nrept = 10;
kFold = 5;

for i = 1:length(timePt)
    
    temp = performanceLDALedoitHbO(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceLDALedoitHbOMean(i) = mean(temp(:));
    
    performanceLDALedoitHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceLDALedoitHbR(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceLDALedoitHbRMean(i) = mean(temp(:));
    
    performanceLDALedoitHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceLDALedoitHbT(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceLDALedoitHbTMean(i) = mean(temp(:));
    
    performanceLDALedoitHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % CERNN
    temp = performanceLDACERNNHbO(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceLDACERNNHbOMean(i) = mean(temp(:));
    
    performanceLDACERNNHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceLDACERNNHbR(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceLDACERNNHbRMean(i) = mean(temp(:));
    
    performanceLDACERNNHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceLDACERNNHbT(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceLDACERNNHbTMean(i) = mean(temp(:));
    
    performanceLDACERNNHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % Log
    temp = performanceLogRegHbO(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceLogRegHbOMean(i) = mean(temp(:));
    
    performanceLogRegHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceLogRegHbR(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceLogRegHbRMean(i) = mean(temp(:));
    
    performanceLogRegHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceLogRegHbT(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceLogRegHbTMean(i) = mean(temp(:));
    
    performanceLogRegHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % SVM
    temp = performanceSVMHbO(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceSVMHbOMean(i) = mean(temp(:));
    
    performanceSVMHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceSVMHbR(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceSVMHbRMean(i) = mean(temp(:));
    
    performanceSVMHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceSVMHbT(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceSVMHbTMean(i) = mean(temp(:));
    
    performanceSVMHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % Bagging
    temp = performanceBaggingHbO(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBaggingHbOMean(i) = mean(temp(:));
    
    performanceBaggingHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceBaggingHbR(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBaggingHbRMean(i) = mean(temp(:));
    
    performanceBaggingHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceBaggingHbT(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBaggingHbTMean(i) = mean(temp(:));
    
    performanceBaggingHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % Boosting
    temp = performanceBoostedHbO(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbOMean(i) = mean(temp(:));
    
    performanceBoostedHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceBoostedHbR(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbRMean(i) = mean(temp(:));
    
    performanceBoostedHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceBoostedHbT(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbTMean(i) = mean(temp(:));
    performanceBoostedHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
end

performanceLDALedoitHbOMean = mean(performanceLDALedoitHbO,1);
performanceLDALedoitHbRMean = mean(performanceLDALedoitHbR,1);
performanceLDALedoitHbTMean = mean(performanceLDALedoitHbT,1);

performanceLDACERNNHbOMean = mean(performanceLDACERNNHbO,1);
performanceLDACERNNHbRMean = mean(performanceLDACERNNHbR,1);
performanceLDACERNNHbTMean = mean(performanceLDACERNNHbT,1);

performanceLogRegHbOMean = mean(performanceLogRegHbO,1);
performanceLogRegHbRMean = mean(performanceLogRegHbR,1);
performanceLogRegHbTMean = mean(performanceLogRegHbT,1);

performanceSVMHbOMean = mean(performanceSVMHbO,1);
performanceSVMHbRMean = mean(performanceSVMHbR,1);
performanceSVMHbTMean = mean(performanceSVMHbT,1);

performanceBaggingHbOMean = mean(performanceBaggingHbO,1);
performanceBaggingHbRMean = mean(performanceBaggingHbR,1);
performanceBaggingHbTMean = mean(performanceBaggingHbT,1);

performanceBoostedHbOMean = mean(performanceBoostedHbO,1);
performanceBoostedHbRMean = mean(performanceBoostedHbR,1);
performanceBoostedHbTMean = mean(performanceBoostedHbT,1);

numClassifiers = 6;
cmap = jet(numClassifiers);

figure('units','normalized','outerposition',[0 0 1 1]); hold on;
subplot(1,3,1);hold on;

errorbar(timePt./fs,performanceLDALedoitHbOMean,performanceLDALedoitHbOCI(:,2),'Color',cmap(1,:));
errorbar(timePt./fs,performanceLDACERNNHbOMean,performanceLDACERNNHbOCI(:,2),'Color',cmap(2,:));
% plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
% plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
errorbar(timePt./fs,performanceLogRegHbOMean,performanceLogRegHbOCI(:,2),'Color',cmap(3,:));
errorbar(timePt./fs,performanceSVMHbOMean,performanceSVMHbOCI(:,2),'Color',cmap(4,:));
errorbar(timePt./fs,performanceBaggingHbOMean,performanceBaggingHbOCI(:,2),'Color',cmap(5,:));
errorbar(timePt./fs,performanceBoostedHbOMean,performanceBoostedHbOCI(:,2),'Color',cmap(6,:));

ylim(yLimAxis);
xlim([-2 7]);
title(sprintf('Sbj %s: Δ[HbO] Multi',num2str(sbjNum)));
% legend({'PCA 95%','PCA 85%','PCA 75%','PCA 65%','PCA 55%','LDAShrink',...
%     'QDAShrink','LDAGD','LogReg','SVM'});
% legend({'PCA 95%','PCA 85%','PCA 75%','PCA 65%','PCA 55%','LDAShrink',...
%     'LDA CERNN','LogReg','SVM'});
legend({'LDAShrink','LDA CERNN','LogReg','SVM','Bagging','Boosted'},'Location','southeast');
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

errorbar(timePt./fs,performanceLDALedoitHbRMean,performanceLDALedoitHbRCI(:,2),'Color',cmap(1,:));
errorbar(timePt./fs,performanceLDACERNNHbRMean,performanceLDACERNNHbRCI(:,2),'Color',cmap(2,:));
% plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
% plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
errorbar(timePt./fs,performanceLogRegHbRMean,performanceLogRegHbRCI(:,2),'Color',cmap(3,:));
errorbar(timePt./fs,performanceSVMHbRMean,performanceSVMHbRCI(:,2),'Color',cmap(4,:));
errorbar(timePt./fs,performanceBaggingHbRMean,performanceBaggingHbRCI(:,2),'Color',cmap(5,:));
errorbar(timePt./fs,performanceBoostedHbRMean,performanceBoostedHbRCI(:,2),'Color',cmap(6,:));

ylim(yLimAxis);
xlim([-2 7]);
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

errorbar(timePt./fs,performanceLDALedoitHbTMean,performanceLDALedoitHbTCI(:,2),'Color',cmap(1,:));
errorbar(timePt./fs,performanceLDACERNNHbTMean,performanceLDACERNNHbTCI(:,2),'Color',cmap(2,:));
% plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
% plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
errorbar(timePt./fs,performanceLogRegHbTMean,performanceLogRegHbTCI(:,2),'Color',cmap(3,:));
errorbar(timePt./fs,performanceSVMHbTMean,performanceSVMHbTCI(:,2),'Color',cmap(4,:));
errorbar(timePt./fs,performanceBaggingHbTMean,performanceBaggingHbTCI(:,2),'Color',cmap(5,:));
errorbar(timePt./fs,performanceBoostedHbTMean,performanceBoostedHbTCI(:,2),'Color',cmap(6,:));

ylim(yLimAxis);
xlim([-2 7]);
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
    if rejTrOp
        fn = sprintf('PerformanceCumsumVsTime_DiffStartT_LR_AllChns_SSBetaGLMCV_DiffClassifiers_RejTr');
    else
        fn = sprintf('PerformanceCumsumVsTime_DiffStartT_LR_AllChns_SSBetaGLMCV_DiffClassifiers');
    end
    %fn = sprintf('PerformanceCumsumVsTime_DiffStartT_SSBeta_LR_AllChns_DiffClassifiers');
else
    if rejTrOp
        fn = sprintf('PerformanceCumsumVsTime_DiffStartT_AllChns_SSBetaGLMCV_DiffClassifiers_RejTr');
    else
        fn = sprintf('PerformanceCumsumVsTime_DiffStartT_AllChns_SSBetaGLMCV_DiffClassifiers');
    end
    %fn = sprintf('PerformanceCumsumVsTime_DiffStartT_SSBeta_AllChns_DiffClassifiers');
end

if saveOp == 1
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
end

end