% For figure 7b in Draft 3
% Use sbj 12 and 19 as example
function grpPlotPerf_Fig7b(sbjList,numClasses,saveOp,rejTrOp)

fs = 50;
%timePt = (0:0.25*fs:5*fs)+2*fs;
timePt = (-2*fs:0.25*fs:5*fs);
performanceLDALedoitHbTMean = zeros(length(sbjList),length(timePt));
performanceLDALedoitHbTCI = zeros(length(sbjList),length(timePt),2);

figSaveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll\Figures\Classification'];

for i1 = 1:length(sbjList)
    sbjNum = sbjList{i1};

    processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
    

    if ~exist(figSaveDir,'dir')
        mkdir(figSaveDir);
    end

    
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
    load(fileName,'performanceLDALedoitHbT');

    %performanceLDALedoitHbOCI = zeros(length(timePt),2);
    %performanceLDALedoitHbRCI = zeros(length(timePt),2);
    %performanceLDALedoitHbTCI = zeros(length(timePt),2);

%     performanceLDACERNNHbOCI = zeros(length(timePt),2);
%     performanceLDACERNNHbRCI = zeros(length(timePt),2);
%     performanceLDACERNNHbTCI = zeros(length(timePt),2);
% 
%     performanceLogRegHbOCI = zeros(length(timePt),2);
%     performanceLogRegHbRCI = zeros(length(timePt),2);
%     performanceLogRegHbTCI = zeros(length(timePt),2);
% 
%     performanceSVMHbOCI = zeros(length(timePt),2);
%     performanceSVMHbRCI = zeros(length(timePt),2);
%     performanceSVMHbTCI = zeros(length(timePt),2);
% 
%     performanceBaggingHbOCI = zeros(length(timePt),2);
%     performanceBaggingHbRCI = zeros(length(timePt),2);
%     performanceBaggingHbTCI = zeros(length(timePt),2);
% 
%     performanceBoostedHbOCI = zeros(length(timePt),2);
%     performanceBoostedHbRCI = zeros(length(timePt),2);
%     performanceBoostedHbTCI = zeros(length(timePt),2);

    conf = 0.95;
    nrept = 10;
    kFold = 5;

    for i = 1:length(timePt)

%         temp = performanceLDALedoitHbO(:,i);
%         se = std(temp(:))/sqrt(nrept*kFold);
%         %performanceLDALedoitHbOMean(i) = mean(temp(:));
% 
%         performanceLDALedoitHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
% 
%         temp = performanceLDALedoitHbR(:,i);
%         se = std(temp(:))/sqrt(nrept*kFold);
%         %performanceLDALedoitHbRMean(i) = mean(temp(:));
% 
%         performanceLDALedoitHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        temp = performanceLDALedoitHbT(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        %performanceLDALedoitHbTMean(i) = mean(temp(:));

        performanceLDALedoitHbTCI(i1,i,:) = calcCITDist(nrept*kFold,se,conf);

%         % CERNN
%         temp = performanceLDACERNNHbO(:,i);
%         se = std(temp(:))/sqrt(nrept*kFold);
%         %performanceLDACERNNHbOMean(i) = mean(temp(:));
% 
%         performanceLDACERNNHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
% 
%         temp = performanceLDACERNNHbR(:,i);
%         se = std(temp(:))/sqrt(nrept*kFold);
%         %performanceLDACERNNHbRMean(i) = mean(temp(:));
% 
%         performanceLDACERNNHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
% 
%         temp = performanceLDACERNNHbT(:,i);
%         se = std(temp(:))/sqrt(nrept*kFold);
%         %performanceLDACERNNHbTMean(i) = mean(temp(:));
% 
%         performanceLDACERNNHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
% 
%         % Log
%         temp = performanceLogRegHbO(:,i);
%         se = std(temp(:))/sqrt(nrept*kFold);
%         %performanceLogRegHbOMean(i) = mean(temp(:));
% 
%         performanceLogRegHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
% 
%         temp = performanceLogRegHbR(:,i);
%         se = std(temp(:))/sqrt(nrept*kFold);
%         %performanceLogRegHbRMean(i) = mean(temp(:));
% 
%         performanceLogRegHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
% 
%         temp = performanceLogRegHbT(:,i);
%         se = std(temp(:))/sqrt(nrept*kFold);
%         %performanceLogRegHbTMean(i) = mean(temp(:));
% 
%         performanceLogRegHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
% 
%         % SVM
%         temp = performanceSVMHbO(:,i);
%         se = std(temp(:))/sqrt(nrept*kFold);
%         %performanceSVMHbOMean(i) = mean(temp(:));
% 
%         performanceSVMHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
% 
%         temp = performanceSVMHbR(:,i);
%         se = std(temp(:))/sqrt(nrept*kFold);
%         %performanceSVMHbRMean(i) = mean(temp(:));
% 
%         performanceSVMHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
% 
%         temp = performanceSVMHbT(:,i);
%         se = std(temp(:))/sqrt(nrept*kFold);
%         %performanceSVMHbTMean(i) = mean(temp(:));
% 
%         performanceSVMHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
% 
%         % Bagging
%         temp = performanceBaggingHbO(:,i);
%         se = std(temp(:))/sqrt(nrept*kFold);
%         %performanceBaggingHbOMean(i) = mean(temp(:));
% 
%         performanceBaggingHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
% 
%         temp = performanceBaggingHbR(:,i);
%         se = std(temp(:))/sqrt(nrept*kFold);
%         %performanceBaggingHbRMean(i) = mean(temp(:));
% 
%         performanceBaggingHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
% 
%         temp = performanceBaggingHbT(:,i);
%         se = std(temp(:))/sqrt(nrept*kFold);
%         %performanceBaggingHbTMean(i) = mean(temp(:));
% 
%         performanceBaggingHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
% 
%         % Boosting
%         temp = performanceBoostedHbO(:,i);
%         se = std(temp(:))/sqrt(nrept*kFold);
%         %performanceBoostedHbOMean(i) = mean(temp(:));
% 
%         performanceBoostedHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
% 
%         temp = performanceBoostedHbR(:,i);
%         se = std(temp(:))/sqrt(nrept*kFold);
%         %performanceBoostedHbRMean(i) = mean(temp(:));
% 
%         performanceBoostedHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
% 
%         temp = performanceBoostedHbT(:,i);
%         se = std(temp(:))/sqrt(nrept*kFold);
%         %performanceBoostedHbTMean(i) = mean(temp(:));
%         performanceBoostedHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);

    end

    %performanceLDALedoitHbOMean = mean(performanceLDALedoitHbO,1);
    %performanceLDALedoitHbRMean = mean(performanceLDALedoitHbR,1);
    performanceLDALedoitHbTMean(i1,:) = mean(performanceLDALedoitHbT,1);

%     performanceLDACERNNHbOMean = mean(performanceLDACERNNHbO,1);
%     performanceLDACERNNHbRMean = mean(performanceLDACERNNHbR,1);
%     performanceLDACERNNHbTMean = mean(performanceLDACERNNHbT,1);
% 
%     performanceLogRegHbOMean = mean(performanceLogRegHbO,1);
%     performanceLogRegHbRMean = mean(performanceLogRegHbR,1);
%     performanceLogRegHbTMean = mean(performanceLogRegHbT,1);
% 
%     performanceSVMHbOMean = mean(performanceSVMHbO,1);
%     performanceSVMHbRMean = mean(performanceSVMHbR,1);
%     performanceSVMHbTMean = mean(performanceSVMHbT,1);
% 
%     performanceBaggingHbOMean = mean(performanceBaggingHbO,1);
%     performanceBaggingHbRMean = mean(performanceBaggingHbR,1);
%     performanceBaggingHbTMean = mean(performanceBaggingHbT,1);
% 
%     performanceBoostedHbOMean = mean(performanceBoostedHbO,1);
%     performanceBoostedHbRMean = mean(performanceBoostedHbR,1);
%     performanceBoostedHbTMean = mean(performanceBoostedHbT,1);
end
numSbj = length(sbjList);
cmap = loadDefaultColors(1);

%figure('units','normalized','outerposition',[0 0 1 1]); hold on;
figure();
%subplot(1,3,1);hold on;

for i = 1:length(sbjList)
    errorbar(timePt./fs,performanceLDALedoitHbTMean(i,:),performanceLDALedoitHbTCI(i,:,1),performanceLDALedoitHbTCI(i,:,2),'Color',cmap(i,:));hold on;
end
%     errorbar(timePt./fs,performanceLDACERNNHbOMean,performanceLDACERNNHbOCI(:,2),'Color',cmap(2,:));
%     % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
%     % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
%     errorbar(timePt./fs,performanceLogRegHbOMean,performanceLogRegHbOCI(:,2),'Color',cmap(3,:));
%     errorbar(timePt./fs,performanceSVMHbOMean,performanceSVMHbOCI(:,2),'Color',cmap(4,:));
%     errorbar(timePt./fs,performanceBaggingHbOMean,performanceBaggingHbOCI(:,2),'Color',cmap(5,:));
%     errorbar(timePt./fs,performanceBoostedHbOMean,performanceBoostedHbOCI(:,2),'Color',cmap(6,:));

ylim(yLimAxis);
xlim([-2 7]);
%title(sprintf('Δ[HbT] for Selected Subjects'));
% legend({'PCA 95%','PCA 85%','PCA 75%','PCA 65%','PCA 55%','LDAShrink',...
%     'QDAShrink','LDAGD','LogReg','SVM'});
% legend({'PCA 95%','PCA 85%','PCA 75%','PCA 65%','PCA 55%','LDAShrink',...
%     'LDA CERNN','LogReg','SVM'});
%legend({'LDAShrink','LDA CERNN','LogReg','SVM','Bagging','Boosted'},'Location','southeast');
legend({'Sbj 08','Sbj 12','Sbj 13','Sbj 16','Sbj 19','Sbj 21','Sbj 24'},'Location','southeast');
xlabel('Time [s]');ylabel('Accuracy');
hold off;


fn = sprintf('Fig7b_%sClasses',num2str(numClasses));

if saveOp == 1
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
end

end