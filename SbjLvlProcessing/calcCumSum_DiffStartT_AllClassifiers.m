% This is used in
% preprocessFNIRS06_CV_GLM_ssBeta_MultiOnly_AllClassifiers.m
%
% Add more classifiers to validate logistic regression and svm unusual
% results. Only test one time point to expedite process. Not for
% publication
%
% Each array is channels x time x trial
% This use cross-validation where GLM is computed for each training fold
% and SS beta coefficients from training folds are used for test fold.

function performanceStr = ...
    calcCumSum_DiffStartT_AllClassifiers(sbjNum,trialsTr,trialsTst,movieListTrain,movieListTest,numClasses,mlActAuto)

%rawDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];

fs = 50;
%timePt = (0:0.25*fs:5*fs)+2*fs;
timePt = (2+2)*fs;
zeroT = 2*fs;
numChn = 30;

if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
    isSDNew = 0;
else
    isSDNew = 1;
end

%figure(1);hold on;
%figure('units','normalized','outerposition',[0 0 1 1]);hold on;

% if strcmp(sbjNum,'15')
%     movieListTrain = movieListTrain(2:end,:);
% end

% Offset
trialsTr = offsetTrials(trialsTr,zeroT);
trialsTst = offsetTrials(trialsTst,zeroT);

% after convert2SD2, cut down from 42 chns to 36 chns. Remove old chns
if ~isSDNew
    [trialsTr,mlList] = convert2SD2(trialsTr,mlActAuto{1});
    [trialsTst] = convert2SD2(trialsTst,mlActAuto{1});
else
    mlList = mlActAuto{1};
end

% after selectRS, cut down from 36 to 30 chns. Remove SS. Both steps can be
% implemented simultaneously using bit logic.
[trialsTr] = selectRS(trialsTr,1,mlList);
[trialsTst,mlList] = selectRS(trialsTst,1,mlList);

% behFN = [rawDataDir filesep 'responses_' sbjNum];

% [trialsTr,movieIdx] = keepCorrectTrials_Special(sbjNum,trialsTr,behFN,numClasses,movieListTrain);
% [trialsTst,movieListTest] = keepCorrectTrials_Special(sbjNum,trialsTst,behFN,numClasses,movieListTest);

performanceLDALedoit = zeros(1,length(timePt));
performanceLDACERNN = zeros(1,length(timePt));
performanceLogReg = zeros(1,length(timePt));
performanceSVM = zeros(1,length(timePt));
performanceBagging = zeros(1,length(timePt));
performanceBoosted = zeros(1,length(timePt));

% Extra classifiers
performanceCoarseKNN = zeros(1,length(timePt));
performanceCoarseTree = zeros(1,length(timePt));
performanceCosineKNN = zeros(1,length(timePt));
performanceCubicKNN = zeros(1,length(timePt));
performanceCubicSVM = zeros(1,length(timePt));
performanceFineGaussianSVM = zeros(1,length(timePt));
performanceFineKNN = zeros(1,length(timePt));
performanceFineTree = zeros(1,length(timePt));
performanceGaussianNaiveBayes = zeros(1,length(timePt));
performanceKernelNaiveBayes = zeros(1,length(timePt));
performanceWeightedKNN = zeros(1,length(timePt));
performanceRandom = zeros(1,length(timePt));

% For fun. This is useful for covariance matrices, which are symmetric positive
% definite matrices
% However, less effective than shrinked LDA.
% https://hal.archives-ouvertes.fr/hal-00602700/document
%performanceRiemannGeometry = zeros(1,length(timePt));


lenT = 1*fs;

for i2 = 1:length(timePt)
    
    temp = cumsum(trialsTr(:,timePt(i2):timePt(i2)+lenT,:),2);
    cumsumTr = squeeze(temp(:,end,:));
    
    temp = cumsum(trialsTst(:,timePt(i2):timePt(i2)+lenT,:),2);
    cumsumTst = squeeze(temp(:,end,:));
    
    %[performanceLDALedoitHbO(1,i2)] = ...
    %    train_RLDA_Ledoit_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieIdx,movieListTest);
    
    [performanceLDALedoit(1,i2)] = train_RLDA_Ledoit_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    [performanceLDACERNN(1,i2)] = train_RLDA_CERNN_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    [performanceLogReg(1,i2)] = trainClassifierLogisticRegressionHRF_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    [performanceSVM(1,i2)] = trainClassifierSVMHRF_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    [performanceBagging(1,i2)] = trainClassifierBaggingTreesSingleHRF_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    [performanceBoosted(1,i2)] = trainClassifierBoostedTreesSingleHRF_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    
    % extra classifiers
    performanceCoarseKNN(1,i2) = trainClassifierCoarseKNN_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    performanceCoarseTree(1,i2) = trainClassifierCoarseTree_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    performanceCosineKNN(1,i2) = trainClassifierCosineKNN_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    performanceCubicKNN(1,i2) = trainCubicKNN_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    performanceCubicSVM(1,i2) = trainCubicSVM_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    performanceFineGaussianSVM(1,i2) = trainFineGaussianSVM_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    performanceFineKNN(1,i2) = trainFineKNN_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    performanceFineTree(1,i2) = trainFineTree_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    performanceGaussianNaiveBayes(1,i2) = trainGaussianNaiveBayes_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    performanceKernelNaiveBayes(1,i2) = trainKernelNaiveBayes_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    performanceWeightedKNN(1,i2) = trainWeightedKNN_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    performanceRandom(1,i2) = trainRandom(cumsumTst,numClasses,movieListTest);
end

performanceStr.performanceLDALedoit = performanceLDALedoit;
performanceStr.performanceLDACERNN = performanceLDACERNN;
performanceStr.performanceLogReg = performanceLogReg;
performanceStr.performanceSVM = performanceSVM;
performanceStr.performanceBagging = performanceBagging;
performanceStr.performanceBoosted = performanceBoosted;

performanceStr.performanceCoarseKNN = performanceCoarseKNN;
performanceStr.performanceCoarseTree = performanceCoarseTree;
performanceStr.performanceCosineKNN = performanceCosineKNN;
performanceStr.performanceCubicKNN = performanceCubicKNN;
performanceStr.performanceCubicSVM = performanceCubicSVM;
performanceStr.performanceFineGaussianSVM = performanceFineGaussianSVM;
performanceStr.performanceFineKNN = performanceFineKNN;
performanceStr.performanceFineTree = performanceFineTree;
performanceStr.performanceGaussianNaiveBayes = performanceGaussianNaiveBayes;
performanceStr.performanceKernelNaiveBayes = performanceKernelNaiveBayes;
performanceStr.performanceWeightedKNN = performanceWeightedKNN;
performanceStr.performanceRandom = performanceRandom;

end