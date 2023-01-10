% Bar chart for more classifier option, fixed at 2s after cue onset.
% Not for publication. For validation against strange logistic
% regression/linear SVM results

function plotPerformance_CV_GLM_ssBeta_AllClassifiers(sbjNum,numClasses,saveOp)

processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
figSaveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' ...
    num2str(sbjNum) '\Figures\Classification'];

fs = 50;
%timePt = (0:0.25*fs:5*fs)+2*fs;
timePt = (2+2)*fs;
yLimAxis = [0 1];

if numClasses == 2
    fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_LR_AllClassifiers.mat'];
else
    fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_AllClassifiers.mat'];
end
% folds x time
load(fileName,'performanceStrHbO','performanceStrHbR',...
    'performanceStrHbT');

performanceLDALedoitHbO = performanceStrHbO.performanceLDALedoit;
performanceLDACERNNHbO = performanceStrHbO.performanceLDACERNN;
performanceLogRegHbO = performanceStrHbO.performanceLogReg;
performanceSVMHbO = performanceStrHbO.performanceSVM;
performanceBaggingHbO = performanceStrHbO.performanceBagging;
performanceBoostedHbO = performanceStrHbO.performanceBoosted;

performanceCoarseKNNHbO = performanceStrHbO.performanceCoarseKNN;
performanceCoarseTreeHbO = performanceStrHbO.performanceCoarseTree;
performanceCosineKNNHbO = performanceStrHbO.performanceCosineKNN;
performanceCubicKNNHbO = performanceStrHbO.performanceCubicKNN;
performanceCubicSVMHbO = performanceStrHbO.performanceCubicSVM;
performanceFineGaussianSVMHbO = performanceStrHbO.performanceFineGaussianSVM;
performanceFineKNNHbO = performanceStrHbO.performanceFineKNN;
performanceFineTreeHbO = performanceStrHbO.performanceFineTree;
performanceGaussianNaiveBayesHbO = performanceStrHbO.performanceGaussianNaiveBayes;
performanceKernelNaiveBayesHbO = performanceStrHbO.performanceKernelNaiveBayes;
performanceWeightedKNNHbO = performanceStrHbO.performanceWeightedKNN;
performanceRandomHbO = performanceStrHbO.performanceRandom;

% HbR
performanceLDALedoitHbR = performanceStrHbR.performanceLDALedoit;
performanceLDACERNNHbR = performanceStrHbR.performanceLDACERNN;
performanceLogRegHbR = performanceStrHbR.performanceLogReg;
performanceSVMHbR = performanceStrHbR.performanceSVM;
performanceBaggingHbR = performanceStrHbR.performanceBagging;
performanceBoostedHbR = performanceStrHbR.performanceBoosted;

performanceCoarseKNNHbR = performanceStrHbR.performanceCoarseKNN;
performanceCoarseTreeHbR = performanceStrHbR.performanceCoarseTree;
performanceCosineKNNHbR = performanceStrHbR.performanceCosineKNN;
performanceCubicKNNHbR = performanceStrHbR.performanceCubicKNN;
performanceCubicSVMHbR = performanceStrHbR.performanceCubicSVM;
performanceFineGaussianSVMHbR = performanceStrHbR.performanceFineGaussianSVM;
performanceFineKNNHbR = performanceStrHbR.performanceFineKNN;
performanceFineTreeHbR = performanceStrHbR.performanceFineTree;
performanceGaussianNaiveBayesHbR = performanceStrHbR.performanceGaussianNaiveBayes;
performanceKernelNaiveBayesHbR = performanceStrHbR.performanceKernelNaiveBayes;
performanceWeightedKNNHbR = performanceStrHbR.performanceWeightedKNN;
performanceRandomHbR = performanceStrHbR.performanceRandom;

% HbT
performanceLDALedoitHbT = performanceStrHbR.performanceLDALedoit;
performanceLDACERNNHbT = performanceStrHbR.performanceLDACERNN;
performanceLogRegHbT = performanceStrHbR.performanceLogReg;
performanceSVMHbT = performanceStrHbR.performanceSVM;
performanceBaggingHbT = performanceStrHbR.performanceBagging;
performanceBoostedHbT = performanceStrHbR.performanceBoosted;

performanceCoarseKNNHbT = performanceStrHbT.performanceCoarseKNN;
performanceCoarseTreeHbT = performanceStrHbT.performanceCoarseTree;
performanceCosineKNNHbT = performanceStrHbT.performanceCosineKNN;
performanceCubicKNNHbT = performanceStrHbT.performanceCubicKNN;
performanceCubicSVMHbT = performanceStrHbT.performanceCubicSVM;
performanceFineGaussianSVMHbT = performanceStrHbT.performanceFineGaussianSVM;
performanceFineKNNHbT = performanceStrHbT.performanceFineKNN;
performanceFineTreeHbT = performanceStrHbT.performanceFineTree;
performanceGaussianNaiveBayesHbT = performanceStrHbT.performanceGaussianNaiveBayes;
performanceKernelNaiveBayesHbT = performanceStrHbT.performanceKernelNaiveBayes;
performanceWeightedKNNHbT = performanceStrHbT.performanceWeightedKNN;
performanceRandomHbT = performanceStrHbT.performanceRandom;

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


performanceCoarseKNNHbOCI = zeros(length(timePt),2);
performanceCoarseKNNHbRCI = zeros(length(timePt),2);
performanceCoarseKNNHbTCI = zeros(length(timePt),2);

performanceCoarseTreeHbOCI = zeros(length(timePt),2);
performanceCoarseTreeHbRCI = zeros(length(timePt),2);
performanceCoarseTreeHbTCI = zeros(length(timePt),2);

performanceCosineKNNHbOCI = zeros(length(timePt),2);
performanceCosineKNNHbRCI = zeros(length(timePt),2);
performanceCosineKNNHbTCI = zeros(length(timePt),2);

performanceCubicKNNHbOCI = zeros(length(timePt),2);
performanceCubicKNNHbRCI = zeros(length(timePt),2);
performanceCubicKNNHbTCI = zeros(length(timePt),2);

performanceCubicSVMHbOCI = zeros(length(timePt),2);
performanceCubicSVMHbRCI = zeros(length(timePt),2);
performanceCubicSVMHbTCI = zeros(length(timePt),2);

performanceFineGaussianSVMHbOCI = zeros(length(timePt),2);
performanceFineGaussianSVMHbRCI = zeros(length(timePt),2);
performanceFineGaussianSVMHbTCI = zeros(length(timePt),2);

performanceFineKNNHbOCI = zeros(length(timePt),2);
performanceFineKNNHbRCI = zeros(length(timePt),2);
performanceFineKNNHbTCI = zeros(length(timePt),2);

performanceFineTreeHbOCI = zeros(length(timePt),2);
performanceFineTreeHbRCI = zeros(length(timePt),2);
performanceFineTreeHbTCI = zeros(length(timePt),2);

performanceGaussianNaiveBayesHbOCI = zeros(length(timePt),2);
performanceGaussianNaiveBayesHbRCI = zeros(length(timePt),2);
performanceGaussianNaiveBayesHbTCI = zeros(length(timePt),2);

performanceKernelNaiveBayesHbOCI = zeros(length(timePt),2);
performanceKernelNaiveBayesHbRCI = zeros(length(timePt),2);
performanceKernelNaiveBayesHbTCI = zeros(length(timePt),2);

performanceWeightedKNNHbOCI = zeros(length(timePt),2);
performanceWeightedKNNHbRCI = zeros(length(timePt),2);
performanceWeightedKNNHbTCI = zeros(length(timePt),2);

performanceRandomHbOCI = zeros(length(timePt),2);
performanceRandomHbRCI = zeros(length(timePt),2);
performanceRandomHbTCI = zeros(length(timePt),2);

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
    
    %Coarse KNN
        temp = performanceCoarseKNNHbO(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbOMean(i) = mean(temp(:));
    
    performanceCoarseKNNHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceCoarseKNNHbR(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbRMean(i) = mean(temp(:));
    
    performanceCoarseKNNHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceCoarseKNNHbT(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbTMean(i) = mean(temp(:));
    performanceCoarseKNNHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % Coarse Tree
        temp = performanceCoarseTreeHbO(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbOMean(i) = mean(temp(:));
    
    performanceCoarseTreeHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceCoarseTreeHbR(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbRMean(i) = mean(temp(:));
    
    performanceCoarseTreeHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceCoarseTreeHbT(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbTMean(i) = mean(temp(:));
    performanceCoarseTreeHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % Cosine KNN
        temp = performanceCosineKNNHbO(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbOMean(i) = mean(temp(:));
    
    performanceCosineKNNHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceCosineKNNHbR(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbRMean(i) = mean(temp(:));
    
    performanceCosineKNNHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceCosineKNNHbT(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbTMean(i) = mean(temp(:));
    performanceCosineKNNHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % Cubic KNN
        temp = performanceCubicKNNHbO(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbOMean(i) = mean(temp(:));
    
    performanceCubicKNNHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceCubicKNNHbR(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbRMean(i) = mean(temp(:));
    
    performanceCubicKNNHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceCubicKNNHbT(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbTMean(i) = mean(temp(:));
    performanceCubicKNNHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    %Cubic SVM
        temp = performanceCubicSVMHbO(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbOMean(i) = mean(temp(:));
    
    performanceCubicSVMHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceCubicSVMHbR(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbRMean(i) = mean(temp(:));
    
    performanceCubicSVMHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceCubicSVMHbT(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbTMean(i) = mean(temp(:));
    performanceCubicSVMHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % Fine Gaussian SVM
        temp = performanceFineGaussianSVMHbO(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbOMean(i) = mean(temp(:));
    
    performanceFineGaussianSVMHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceFineGaussianSVMHbR(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbRMean(i) = mean(temp(:));
    
    performanceFineGaussianSVMHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceFineGaussianSVMHbT(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbTMean(i) = mean(temp(:));
    performanceFineGaussianSVMHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % Fine KNN
        temp = performanceFineKNNHbO(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbOMean(i) = mean(temp(:));
    
    performanceFineKNNHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceFineKNNHbR(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbRMean(i) = mean(temp(:));
    
    performanceFineKNNHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceFineKNNHbT(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbTMean(i) = mean(temp(:));
    performanceFineKNNHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % Fine Tree
        temp = performanceFineTreeHbO(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbOMean(i) = mean(temp(:));
    
    performanceFineTreeHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceFineTreeHbR(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbRMean(i) = mean(temp(:));
    
    performanceFineTreeHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceFineTreeHbT(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbTMean(i) = mean(temp(:));
    performanceFineTreeHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % Gaussian Naive Bayes
        temp = performanceGaussianNaiveBayesHbO(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbOMean(i) = mean(temp(:));
    
    performanceGaussianNaiveBayesHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceGaussianNaiveBayesHbR(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbRMean(i) = mean(temp(:));
    
    performanceGaussianNaiveBayesHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceGaussianNaiveBayesHbT(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbTMean(i) = mean(temp(:));
    performanceGaussianNaiveBayesHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % Kernel Naive Bayes
        temp = performanceKernelNaiveBayesHbO(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbOMean(i) = mean(temp(:));
    
    performanceKernelNaiveBayesHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceKernelNaiveBayesHbR(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbRMean(i) = mean(temp(:));
    
    performanceKernelNaiveBayesHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceKernelNaiveBayesHbT(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbTMean(i) = mean(temp(:));
    performanceKernelNaiveBayesHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % Weighted KNN
        temp = performanceWeightedKNNHbO(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbOMean(i) = mean(temp(:));
    
    performanceWeightedKNNHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceWeightedKNNHbR(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbRMean(i) = mean(temp(:));
    
    performanceWeightedKNNHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceWeightedKNNHbT(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbTMean(i) = mean(temp(:));
    performanceWeightedKNNHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % Random
        temp = performanceRandomHbO(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbOMean(i) = mean(temp(:));
    
    performanceRandomHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceRandomHbR(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbRMean(i) = mean(temp(:));
    
    performanceRandomHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceRandomHbT(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    %performanceBoostedHbTMean(i) = mean(temp(:));
    performanceRandomHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
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

performanceCoarseKNNHbOMean = mean(performanceCoarseKNNHbO,1);
performanceCoarseKNNHbRMean = mean(performanceCoarseKNNHbR,1);
performanceCoarseKNNHbTMean = mean(performanceCoarseKNNHbT,1);

performanceCoarseTreeHbOMean = mean(performanceCoarseTreeHbO,1);
performanceCoarseTreeHbRMean = mean(performanceCoarseTreeHbR,1);
performanceCoarseTreeHbTMean = mean(performanceCoarseTreeHbT,1);

performanceCosineKNNHbOMean = mean(performanceCosineKNNHbO,1);
performanceCosineKNNHbRMean = mean(performanceCosineKNNHbR,1);
performanceCosineKNNHbTMean = mean(performanceCosineKNNHbT,1);

performanceCubicKNNHbOMean = mean(performanceCubicKNNHbO,1);
performanceCubicKNNHbRMean = mean(performanceCubicKNNHbR,1);
performanceCubicKNNHbTMean = mean(performanceCubicKNNHbT,1);

performanceCubicSVMHbOMean = mean(performanceCubicSVMHbO,1);
performanceCubicSVMHbRMean = mean(performanceCubicSVMHbR,1);
performanceCubicSVMHbTMean = mean(performanceCubicSVMHbT,1);

performanceFineGaussianSVMHbOMean = mean(performanceFineGaussianSVMHbO,1);
performanceFineGaussianSVMHbRMean = mean(performanceFineGaussianSVMHbR,1);
performanceFineGaussianSVMHbTMean = mean(performanceFineGaussianSVMHbT,1);

performanceFineKNNHbOMean = mean(performanceFineKNNHbO,1);
performanceFineKNNHbRMean = mean(performanceFineKNNHbR,1);
performanceFineKNNHbTMean = mean(performanceFineKNNHbT,1);

performanceFineTreeHbOMean = mean(performanceFineTreeHbO,1);
performanceFineTreeHbRMean = mean(performanceFineTreeHbR,1);
performanceFineTreeHbTMean = mean(performanceFineTreeHbT,1);

performanceGaussianNaiveBayesHbOMean = mean(performanceGaussianNaiveBayesHbO,1);
performanceGaussianNaiveBayesHbRMean = mean(performanceGaussianNaiveBayesHbR,1);
performanceGaussianNaiveBayesHbTMean = mean(performanceGaussianNaiveBayesHbT,1);

performanceKernelNaiveBayesHbOMean = mean(performanceKernelNaiveBayesHbO,1);
performanceKernelNaiveBayesHbRMean = mean(performanceKernelNaiveBayesHbR,1);
performanceKernelNaiveBayesHbTMean = mean(performanceKernelNaiveBayesHbT,1);

performanceWeightedKNNHbOMean = mean(performanceWeightedKNNHbO,1);
performanceWeightedKNNHbRMean = mean(performanceWeightedKNNHbR,1);
performanceWeightedKNNHbTMean = mean(performanceWeightedKNNHbT,1);

performanceRandomHbOMean = mean(performanceRandomHbO,1);
performanceRandomHbRMean = mean(performanceRandomHbR,1);
performanceRandomHbTMean = mean(performanceRandomHbT,1);

res_HbO = [performanceLDALedoitHbOMean performanceLDACERNNHbOMean ...
    performanceLogRegHbOMean performanceSVMHbOMean...
    performanceBaggingHbOMean performanceBoostedHbOMean...
    performanceCoarseKNNHbOMean performanceCoarseTreeHbOMean...
    performanceCosineKNNHbOMean performanceCubicKNNHbOMean...
    performanceCubicSVMHbOMean performanceFineGaussianSVMHbOMean...
    performanceFineKNNHbOMean performanceFineTreeHbOMean...
    performanceGaussianNaiveBayesHbOMean performanceKernelNaiveBayesHbOMean...
    performanceWeightedKNNHbOMean performanceRandomHbOMean];

res_HbR = [performanceLDALedoitHbRMean performanceLDACERNNHbRMean ...
    performanceLogRegHbRMean performanceSVMHbRMean...
    performanceBaggingHbRMean performanceBoostedHbRMean...
    performanceCoarseKNNHbRMean performanceCoarseTreeHbRMean...
    performanceCosineKNNHbRMean performanceCubicKNNHbRMean...
    performanceCubicSVMHbRMean performanceFineGaussianSVMHbRMean...
    performanceFineKNNHbRMean performanceFineTreeHbRMean...
    performanceGaussianNaiveBayesHbRMean performanceKernelNaiveBayesHbRMean...
    performanceWeightedKNNHbRMean performanceRandomHbRMean];

res_HbT = [performanceLDALedoitHbTMean performanceLDACERNNHbTMean ...
    performanceLogRegHbTMean performanceSVMHbTMean...
    performanceBaggingHbTMean performanceBoostedHbTMean...
    performanceCoarseKNNHbTMean performanceCoarseTreeHbTMean...
    performanceCosineKNNHbTMean performanceCubicKNNHbTMean...
    performanceCubicSVMHbTMean performanceFineGaussianSVMHbTMean...
    performanceFineKNNHbTMean performanceFineTreeHbTMean...
    performanceGaussianNaiveBayesHbTMean performanceKernelNaiveBayesHbTMean...
    performanceWeightedKNNHbTMean performanceRandomHbTMean];
    
[~, origInd_HbO] = sort(res_HbO,'descend');
[res_HbR_Sorted, origInd_HbR] = sort(res_HbR,'descend');
[res_HbT_Sorted, origInd_HbT] = sort(res_HbT,'descend');

numClassifiers = length(res_HbO);
%cmap = jet(numClassifiers);

ClassifierLabels = ({'LDALedoit',...
    'LDACERNN','LogReg'...
    'LinearSVM','Bagging','Boosting'...
    'CoarseKNN','CoarseTree'...
    'CosineKNN','CubicKNN',...
    'CubicSVM','FineGaussianSVM',...
    'FineKNN','FineTree',...
    'GaussianNaiveBayes',...
    'KernelNaiveBayes',...
    'WeightedKNN','Random'});

% Bar plot
ClassifierLabels = categorical(ClassifierLabels);
ClassifierLabels = reordercats(ClassifierLabels,string(ClassifierLabels));

%bar(XSorted,perfVar);
hold on;
%b = bar(1:numClassifiers,res_HbO(origInd_HbO),'FaceColor','flat');
bar(1:numClassifiers,res_HbO(origInd_HbO),'FaceColor','flat');

% for i = 1:length(LSChnList)
%     % acceptedChns has 36 chns, fixed it down to 30 chns
%     if ~acceptedChns(i)
%         b.CData(origInd(i),:) = [0 0.8 0.8];
%     end
% end

set(gca,'XTick',[1:1:numClassifiers]);
set(gca, 'XTickLabel', ClassifierLabels(origInd_HbO));
xtickangle(45);
ylabel('Accuracy [Decimal]');
title('HbO 2s');
ylim([0 1]);









% % Line plot
% subplot(1,3,1);hold on;
% 
% errorbar(timePt./fs-2,performanceLDALedoitHbOMean,performanceLDALedoitHbOCI(:,2),'Color',cmap(1,:));
% errorbar(timePt./fs-2,performanceLDACERNNHbOMean,performanceLDACERNNHbOCI(:,2),'Color',cmap(2,:));
% % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
% % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
% errorbar(timePt./fs-2,performanceLogRegHbOMean,performanceLogRegHbOCI(:,2),'Color',cmap(3,:));
% errorbar(timePt./fs-2,performanceSVMHbOMean,performanceSVMHbOCI(:,2),'Color',cmap(4,:));
% errorbar(timePt./fs-2,performanceBaggingHbOMean,performanceBaggingHbOCI(:,2),'Color',cmap(5,:));
% errorbar(timePt./fs-2,performanceBoostedHbOMean,performanceBoostedHbOCI(:,2),'Color',cmap(6,:));
% 
% ylim(yLimAxis);
% xlim([0 7]);
% title(sprintf('Sbj %s: Δ[HbO] Multi',num2str(sbjNum)));
% % legend({'PCA 95%','PCA 85%','PCA 75%','PCA 65%','PCA 55%','LDAShrink',...
% %     'QDAShrink','LDAGD','LogReg','SVM'});
% % legend({'PCA 95%','PCA 85%','PCA 75%','PCA 65%','PCA 55%','LDAShrink',...
% %     'LDA CERNN','LogReg','SVM'});
% legend({'LDAShrink','LDA CERNN','LogReg','SVM','Bagging','Boosted'});
% xlabel('Time [s]');ylabel('Accuracy');
% hold off;
% 
% subplot(1,3,2);hold on;
% % for j = 1:size(performanceArrMultiHbR,1)
% %     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
% %     plot(timePt./fs-2,performanceArrMultiHbR(j,:),'Color',cmap(j,:));
% % end
% % plot(timePt./fs-2,performanceArrMultiHbR_95,'Color',cmap(1,:));
% % plot(timePt./fs-2,performanceArrMultiHbR_85,'Color',cmap(2,:));
% % plot(timePt./fs-2,performanceArrMultiHbR_75,'Color',cmap(3,:));
% % plot(timePt./fs-2,performanceArrMultiHbR_65,'Color',cmap(4,:));
% % plot(timePt./fs-2,performanceArrMultiHbR_55,'Color',cmap(5,:));
% % plot(timePt./fs-2,performanceLDALedoitHbR,'Color',cmap(6,:));
% % plot(timePt./fs-2,performanceLDACERNNHbR,'Color',cmap(7,:));
% % % plot(timePt./fs-2,performanceQDALedoitHbR,'Color',cmap(7,:));
% % % plot(timePt./fs-2,performanceLDAGDHbR,'Color',cmap(8,:));
% % plot(timePt./fs-2,performanceLogRegHbR,'Color',cmap(8,:));
% % plot(timePt./fs-2,performanceSVMHbR,'Color',cmap(9,:));
% 
% errorbar(timePt./fs-2,performanceLDALedoitHbRMean,performanceLDALedoitHbRCI(:,2),'Color',cmap(1,:));
% errorbar(timePt./fs-2,performanceLDACERNNHbRMean,performanceLDACERNNHbRCI(:,2),'Color',cmap(2,:));
% % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
% % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
% errorbar(timePt./fs-2,performanceLogRegHbRMean,performanceLogRegHbRCI(:,2),'Color',cmap(3,:));
% errorbar(timePt./fs-2,performanceSVMHbRMean,performanceSVMHbRCI(:,2),'Color',cmap(4,:));
% errorbar(timePt./fs-2,performanceBaggingHbRMean,performanceBaggingHbRCI(:,2),'Color',cmap(5,:));
% errorbar(timePt./fs-2,performanceBoostedHbRMean,performanceBoostedHbRCI(:,2),'Color',cmap(6,:));
% 
% ylim(yLimAxis);
% xlim([0 7]);
% title(sprintf('Sbj %s: Δ[HbR] Multi',num2str(sbjNum)));
% %legend('S1D1','S1D2','S1D3','S1D4');
% xlabel('Time [s]');ylabel('Accuracy');
% hold off;
% 
% subplot(1,3,3);hold on;
% % for j = 1:size(performanceArrMultiHbT,1)
% %     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
% %     plot(timePt./fs-2,performanceArrMultiHbT(j,:),'Color',cmap(j,:));
% % end
% % plot(timePt./fs-2,performanceArrMultiHbT_95,'Color',cmap(1,:));
% % plot(timePt./fs-2,performanceArrMultiHbT_85,'Color',cmap(2,:));
% % plot(timePt./fs-2,performanceArrMultiHbT_75,'Color',cmap(3,:));
% % plot(timePt./fs-2,performanceArrMultiHbT_65,'Color',cmap(4,:));
% % plot(timePt./fs-2,performanceArrMultiHbT_55,'Color',cmap(5,:));
% % plot(timePt./fs-2,performanceLDALedoitHbT,'Color',cmap(6,:));
% % plot(timePt./fs-2,performanceLDACERNNHbT,'Color',cmap(7,:));
% % % plot(timePt./fs-2,performanceQDALedoitHbT,'Color',cmap(7,:));
% % % plot(timePt./fs-2,performanceLDAGDHbT,'Color',cmap(8,:));
% % plot(timePt./fs-2,performanceLogRegHbT,'Color',cmap(8,:));
% % plot(timePt./fs-2,performanceSVMHbT,'Color',cmap(9,:));
% 
% errorbar(timePt./fs-2,performanceLDALedoitHbTMean,performanceLDALedoitHbTCI(:,2),'Color',cmap(1,:));
% errorbar(timePt./fs-2,performanceLDACERNNHbTMean,performanceLDACERNNHbTCI(:,2),'Color',cmap(2,:));
% % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
% % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
% errorbar(timePt./fs-2,performanceLogRegHbTMean,performanceLogRegHbTCI(:,2),'Color',cmap(3,:));
% errorbar(timePt./fs-2,performanceSVMHbTMean,performanceSVMHbTCI(:,2),'Color',cmap(4,:));
% errorbar(timePt./fs-2,performanceBaggingHbTMean,performanceBaggingHbTCI(:,2),'Color',cmap(5,:));
% errorbar(timePt./fs-2,performanceBoostedHbTMean,performanceBoostedHbTCI(:,2),'Color',cmap(6,:));
% 
% ylim(yLimAxis);
% xlim([0 7]);
% title(sprintf('Sbj %s: Δ[HbT] Multi',num2str(sbjNum)));
% %legend('S1D1','S1D2','S1D3','S1D4');
% xlabel('Time [s]');ylabel('Accuracy');
% hold off;
% 
% % subplot(1,4,4);hold on;
% % for j = 1:size(performanceArrMultiHbT,1)
% %     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
% %     plot(timePt./fs-2,zeros(1,length(timePt./fs-2)),'Color',cmap(j,:));
% % end
% % %legend(chnName(1:30));
% % annotation('textbox', [0.8, 0.05, 0.1, 0.1], 'String', ['Score: ' num2str(behScore) '%'],'FontSize',8);
% % annotation('textbox', [0.71, 0.2, 0.1, 0.1],'String','calcCumSum\_MultiOnly\_DiffStartT\_AllChn\_AllBasis.m','FontSize',8);
% % title(sprintf('All Channels. Diff Start T: %s-Class',num2str(numClasses)));
% % hold off;

if numClasses == 2
    fn = sprintf('PerformanceCumsumVsTime_DiffStartT_LR_AllChns_SSBetaGLMCV_DiffClassifiers');
    %fn = sprintf('PerformanceCumsumVsTime_DiffStartT_SSBeta_LR_AllChns_DiffClassifiers');
else
    fn = sprintf('PerformanceCumsumVsTime_DiffStartT_AllChns_SSBetaGLMCV_DiffClassifiers');
    %fn = sprintf('PerformanceCumsumVsTime_DiffStartT_SSBeta_AllChns_DiffClassifiers');
end

if saveOp == 1
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
end

end