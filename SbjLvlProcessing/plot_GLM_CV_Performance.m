function plot_GLM_CV_Performance(sbjNum,respData,numClasses)

fs = 50;
timePt = (0:0.5*fs:6*fs)+2*fs;
%trialTime = [-2 15];

saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
rawDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
behFN = [rawDataDir filesep 'responses_' sbjNum];
load([saveDir filesep respData '.mat'],'fixedMaskerList','indexMoviesTest','maskerMovies','numTrials','uniqueMovies');

if strcmp(sbjNum,'15')
    indexMoviesTest = indexMoviesTest(2:end,:);
end

load(behFN,'responsesV','responsesA','correctRespV','correctRespA');
idx = (responsesA==correctRespA).*(responsesV==correctRespV);
indexMoviesTest = indexMoviesTest(logical(idx),:);

if numClasses==2
    idx = (indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1) | (indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1);
    %indexMoviesTestCond = indexMoviesTest(idx,:);
    numTr = sum(logical(idx));
    %indexMoviesTestCond(:,6)=1:size(indexMoviesTestCond,1);
else
    idx = (indexMoviesTest(:,5)==1);
    %indexMoviesTestCond = indexMoviesTest(indexMoviesTest(:,5)==indexMoviesTest(i1,5),:);
    numTr = sum(logical(idx));
    %indexMoviesTestCond(:,6)=1:size(indexMoviesTestCond,1);
end

fileName = [processedDataDir filesep 'GLM_CV_Performances.mat'];
load(fileName,'performanceArr_LDA_LW_HbO','performanceArr_LDA_CERNN_HbO',...
    'performanceArr_Logistics_HbO','performanceArr_SVM_Lin_HbO',...
    'performanceArr_Bagging_HbO','performanceArr_Boosting_HbO',...
    'performanceArr_LDA_LW_HbR','performanceArr_LDA_CERNN_HbR',...
    'performanceArr_Logistics_HbR','performanceArr_SVM_Lin_HbR',...
    'performanceArr_Bagging_HbR','performanceArr_Boosting_HbR',...
    'performanceArr_LDA_LW_HbT','performanceArr_LDA_CERNN_HbT',...
    'performanceArr_Logistics_HbT','performanceArr_SVM_Lin_HbT',...
    'performanceArr_Bagging_HbT','performanceArr_Boosting_HbT');

performanceArr_LDA_LW_HbO = performanceArr_LDA_LW_HbO(1:numTr,:);
performanceArr_LDA_CERNN_HbO = performanceArr_LDA_CERNN_HbO(1:numTr,:);
performanceArr_Logistics_HbO = performanceArr_Logistics_HbO(1:numTr,:);
performanceArr_SVM_Lin_HbO = performanceArr_SVM_Lin_HbO(1:numTr,:);
performanceArr_Bagging_HbO = performanceArr_Bagging_HbO(1:numTr,:);
performanceArr_Boosting_HbO = performanceArr_Boosting_HbO(1:numTr,:);

performanceArr_LDA_LW_HbR = performanceArr_LDA_LW_HbR(1:numTr,:);
performanceArr_LDA_CERNN_HbR = performanceArr_LDA_CERNN_HbR(1:numTr,:);
performanceArr_Logistics_HbR = performanceArr_Logistics_HbR(1:numTr,:);
performanceArr_SVM_Lin_HbR = performanceArr_SVM_Lin_HbR(1:numTr,:);
performanceArr_Bagging_HbR = performanceArr_Bagging_HbR(1:numTr,:);
performanceArr_Boosting_HbR = performanceArr_Boosting_HbR(1:numTr,:);

performanceArr_LDA_LW_HbT = performanceArr_LDA_LW_HbT(1:numTr,:);
performanceArr_LDA_CERNN_HbT = performanceArr_LDA_CERNN_HbT(1:numTr,:);
performanceArr_Logistics_HbT = performanceArr_Logistics_HbT(1:numTr,:);
performanceArr_SVM_Lin_HbT = performanceArr_SVM_Lin_HbT(1:numTr,:);
performanceArr_Bagging_HbT = performanceArr_Bagging_HbT(1:numTr,:);
performanceArr_Boosting_HbT = performanceArr_Boosting_HbT(1:numTr,:);

avgPerformance_LDA_LW_HbO = squeeze(mean(performanceArr_LDA_LW_HbO,1));
avgPerformance_LDA_CERNN_HbO = squeeze(mean(performanceArr_LDA_CERNN_HbO,1));
avgPerformance_Logistics_HbO = squeeze(mean(performanceArr_Logistics_HbO,1));
avgPerformance_SVM_Lin_HbO = squeeze(mean(performanceArr_SVM_Lin_HbO,1));
avgPerformance_Bagging_HbO = squeeze(mean(performanceArr_Bagging_HbO,1));
avgPerformance_Boosting_HbO = squeeze(mean(performanceArr_Boosting_HbO,1));

avgPerformance_LDA_LW_HbR = squeeze(mean(performanceArr_LDA_LW_HbR,1));
avgPerformance_LDA_CERNN_HbR = squeeze(mean(performanceArr_LDA_CERNN_HbR,1));
avgPerformance_Logistics_HbR = squeeze(mean(performanceArr_Logistics_HbR,1));
avgPerformance_SVM_Lin_HbR = squeeze(mean(performanceArr_SVM_Lin_HbR,1));
avgPerformance_Bagging_HbR = squeeze(mean(performanceArr_Bagging_HbR,1));
avgPerformance_Boosting_HbR = squeeze(mean(performanceArr_Boosting_HbR,1));

avgPerformance_LDA_LW_HbT = squeeze(mean(performanceArr_LDA_LW_HbT,1));
avgPerformance_LDA_CERNN_HbT = squeeze(mean(performanceArr_LDA_CERNN_HbT,1));
avgPerformance_Logistics_HbT = squeeze(mean(performanceArr_Logistics_HbT,1));
avgPerformance_SVM_Lin_HbT = squeeze(mean(performanceArr_SVM_Lin_HbT,1));
avgPerformance_Bagging_HbT = squeeze(mean(performanceArr_Bagging_HbT,1));
avgPerformance_Boosting_HbT = squeeze(mean(performanceArr_Boosting_HbT,1));

numClassifiers = 6;
cmap = jet(numClassifiers);

figure('units','normalized','outerposition',[0 0 1 1]);hold on;
subplot(1,3,1);hold on;
plot(timePt./fs-2,avgPerformance_LDA_LW_HbO,'Color',cmap(1,:));
plot(timePt./fs-2,avgPerformance_LDA_CERNN_HbO,'Color',cmap(2,:));
plot(timePt./fs-2,avgPerformance_Logistics_HbO,'Color',cmap(3,:));
plot(timePt./fs-2,avgPerformance_SVM_Lin_HbO,'Color',cmap(4,:));
plot(timePt./fs-2,avgPerformance_Bagging_HbO,'Color',cmap(5,:));
plot(timePt./fs-2,avgPerformance_Boosting_HbO,'Color',cmap(6,:));

ylim([0 1]);
title(sprintf('Sbj %s: Δ[HbO] Multi',num2str(sbjNum)));
legend({'LDAShrink','LDA CERNN','LogReg','SVM','Bagging','Boosted'});
xlabel('Time [s]');ylabel('Accuracy');
hold off;

subplot(1,3,2);hold on;
plot(timePt./fs-2,avgPerformance_LDA_LW_HbR,'Color',cmap(1,:));
plot(timePt./fs-2,avgPerformance_LDA_CERNN_HbR,'Color',cmap(2,:));
plot(timePt./fs-2,avgPerformance_Logistics_HbR,'Color',cmap(3,:));
plot(timePt./fs-2,avgPerformance_SVM_Lin_HbR,'Color',cmap(4,:));
plot(timePt./fs-2,avgPerformance_Bagging_HbR,'Color',cmap(5,:));
plot(timePt./fs-2,avgPerformance_Boosting_HbR,'Color',cmap(6,:));

ylim([0 1]);
title(sprintf('Sbj %s: Δ[HbR] Multi',num2str(sbjNum)));
legend({'LDAShrink','LDA CERNN','LogReg','SVM','Bagging','Boosted'});
xlabel('Time [s]');ylabel('Accuracy');
hold off;

subplot(1,3,3);hold on;
plot(timePt./fs-2,avgPerformance_LDA_LW_HbT,'Color',cmap(1,:));
plot(timePt./fs-2,avgPerformance_LDA_CERNN_HbT,'Color',cmap(2,:));
plot(timePt./fs-2,avgPerformance_Logistics_HbT,'Color',cmap(3,:));
plot(timePt./fs-2,avgPerformance_SVM_Lin_HbT,'Color',cmap(4,:));
plot(timePt./fs-2,avgPerformance_Bagging_HbT,'Color',cmap(5,:));
plot(timePt./fs-2,avgPerformance_Boosting_HbT,'Color',cmap(6,:));

ylim([0 1]);
title(sprintf('Sbj %s: Δ[HbT] Multi',num2str(sbjNum)));
legend({'LDAShrink','LDA CERNN','LogReg','SVM','Bagging','Boosted'});
xlabel('Time [s]');ylabel('Accuracy');
hold off;

if numClasses == 2
    fn = sprintf('PerformanceCumsumVsTime_GLM_CV_LR_AllChns_DiffClassifiers');
else
    fn = sprintf('PerformanceCumsumVsTime_GLM_CV_AllChns_DiffClassifiers');
end

end