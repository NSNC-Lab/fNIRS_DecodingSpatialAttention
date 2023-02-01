% STATUS: Active
% 
% SYNTAX:
% [performanceLDALedoit, performanceLDACERNN,...
%   performanceCosineKNN, performanceSVM,...
%   performanceBagging, performanceBoosted] = ...
%   calcCumSum_DiffStartT_AllChn_CompClassifiers_Final_NoSave(sbjNum,...
%       trialsTr,trialsTst,movieListTrain,movieListTest,numClasses,...
%       mlActAuto)
% 
% DESCRIPTION:
% Perform cross validation where SS beta coefficients from training fold 
%   is passed to test fold. Perform all-channel classification using 6 
%   different classifiers:
%       1) linear discriminant analysis using linear estimator of
%       covariance matrix.
%       2) linear discriminant analysis using non-linear estimator of
%       covariance matrix.
%       3) KNN using 1 minus cosine similarity as distance metric
%       4) SVM using cubic polynomial kernel
%       5) random forest ensemble
%       6) boosting ensemble
% 
% RESTRICTION:
% None.
% 
% INPUTS:
% sbjNum - string: subject ID
% trialsTr - training dataset. 3D double array: channel x time x trial
% trialsTst - test dataset. 3D double array: channel x time x trial
% movieListTrain - trials info for training dataset. numTrials x 5 double array:
%       col 1: index of target movies in uniqueMovies
%       col 2: index of spatial location
%       col 3: boolean: masker is fixed or random
%       col 4: index of masker movies in fixedMaskerList(?)
%       col 5: boolean: condition is target-alone or target+maskers
% movieListTest - trials info for test dataset. same structure as movieListTrain
% numClasses - int: number of classes for classification.
% mlActAuto - 1x1 cell array containing 1D int array of channel list of 2
%   different wavelengths
%
% RETURNED VARIABLES:
% performanceLDALedoit - decoding performance of rLDA classifier for 
%   each decision window. 1D double array: 1 x time windows
% performanceLDACERNN - decoding performance of CERNN LDA classifier for 
%   each decision window. 1D double array: 1 x time windows
% performanceCosineKNN - decoding performance of KNN classifier for 
%   each decision window. 1D double array: 1 x time windows
% performanceSVM - decoding performance of SVM classifier for 
%   each decision window. 1D double array: 1 x time windows
% performanceBagging - decoding performance of random forest ensemble 
%   classifier for each decision window. 1D double array: 1 x time windows
% performanceBoosted - decoding performance of boosting ensemble
%   classifier for each decision window. 1D double array: 1 x time windows
% 
% FILES SAVED:
% None.
% 
% PLOTTING:
% None.

function [performanceLDALedoit, performanceLDACERNN,...
    performanceCosineKNN, performanceSVM,...
    performanceBagging, performanceBoosted] = ...
    calcCumSum_DiffStartT_AllChn_CompClassifiers_Final_NoSave(sbjNum,trialsTr,trialsTst,movieListTrain,movieListTest,numClasses,mlActAuto)

rawDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];

fs = 50;
%timePt = (0:0.25*fs:5*fs)+2*fs;
timePt = (0:0.25*fs:7*fs);
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
performanceCosineKNN = zeros(1,length(timePt));
%performanceLogReg = zeros(1,length(timePt));
performanceSVM = zeros(1,length(timePt));
performanceBagging = zeros(1,length(timePt));
performanceBoosted = zeros(1,length(timePt));

lenT = 0.5*fs;

for i2 = 2:length(timePt)
    
    temp = cumsum(trialsTr(:,timePt(i2):timePt(i2)+lenT,:),2);
    cumsumTr = squeeze(temp(:,end,:));
    
    temp = cumsum(trialsTst(:,timePt(i2):timePt(i2)+lenT,:),2);
    cumsumTst = squeeze(temp(:,end,:));
    
    %[performanceLDALedoit(1,i2)] = ...
    %    train_RLDA_Ledoit_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieIdx,movieListTest);
    
    [performanceLDALedoit(1,i2)] = train_RLDA_Ledoit_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    [performanceLDACERNN(1,i2)] = train_RLDA_CERNN_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    % replace this one with Cosine KNN
    %[performanceLogReg(1,i2)] = trainClassifierLogisticRegressionHRF_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    [performanceCosineKNN(1,i2)] = trainClassifierCosineKNN_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    % replace this one with cubic KNN
    %[performanceSVM(1,i2)] = trainClassifierSVMHRF_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    [performanceSVM(1,i2)] = trainCubicSVM_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    [performanceBagging(1,i2)] = trainClassifierBaggingTreesSingleHRF_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    [performanceBoosted(1,i2)] = trainClassifierBoostedTreesSingleHRF_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    
end

end