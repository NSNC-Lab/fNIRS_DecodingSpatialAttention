% STATUS: Inactive(?)
% 
% SYNTAX:
% [performanceLDALedoitHbO] = ...
%   calcCumSum_DiffTLen_AllChn_CompClassifiers_Final_NoSave(sbjNum,...
%       trialsTr,trialsTst,movieListTrain,movieListTest,numClasses,...
%       mlActAuto,timeLen)
% 
% DESCRIPTION:
% Perform cross validation where SS beta coefficients from training fold 
%   is passed to test fold. Perform all-channel classification using rLDA
%   classifier. Test different decision window lengths.
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
% timeLen - 1D double array of window lengths to test, in sec.
%
% RETURNED VARIABLES:
% performanceLDALedoitHbO - decoding performance of rLDA classifier for 
%   each decision window. 1D double array: 1 x time windows
% performanceLDACERNNHbO - decoding performance of CERNN LDA classifier for 
%   each decision window. 1D double array: 1 x time windows
% performanceCosineKNN - decoding performance of KNN classifier for 
%   each decision window. 1D double array: 1 x time windows
% performanceSVMHbO - decoding performance of SVM classifier for 
%   each decision window. 1D double array: 1 x time windows
% performanceBaggingHbO - decoding performance of random forest ensemble 
%   classifier for each decision window. 1D double array: 1 x time windows
% performanceBoostedHbO - decoding performance of boosting ensemble
%   classifier for each decision window. 1D double array: 1 x time windows
% 
% FILES SAVED:
% None.
% 
% PLOTTING:
% None.

function [performanceLDALedoitHbO] = ...
    calcCumSum_DiffTLen_AllChn_CompClassifiers_Final_NoSave(sbjNum,trialsTr,trialsTst,movieListTrain,movieListTest,numClasses,mlActAuto,timeLen)

%rawDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];

fs = 50;
if nargin == 7
    timeLen = [0.1 0.2 0.5 1 1.5 2 3 4 5].*fs;
end
%startT = 4*fs;
startT = 2*fs;
zeroT = 2*fs;


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

numChn = size(trialsTr,1);

performanceLDALedoitHbO = zeros(1,length(timeLen));
% performanceLDACERNNHbO = zeros(1,length(timeLen));
% performanceLogRegHbO = zeros(1,length(timeLen));
% performanceSVMHbO = zeros(1,length(timeLen));
% performanceBaggingHbO = zeros(1,length(timeLen));
% performanceBoostedHbO = zeros(1,length(timeLen));

for i2 = 1:length(timeLen)
    
    temp = cumsum(trialsTr(:,startT:startT+timeLen(i2),:),2);
    cumsumTr = squeeze(temp(:,end,:));
    
    temp = cumsum(trialsTst(:,startT:startT+timeLen(i2),:),2);
    cumsumTst = squeeze(temp(:,end,:));
    
    %[performanceLDALedoitHbO(1,i2)] = ...
    %    train_RLDA_Ledoit_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieIdx,movieListTest);
    
    [performanceLDALedoitHbO(1,i2)] = train_RLDA_Ledoit_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
%     [performanceLDACERNNHbO(1,i2)] = train_RLDA_CERNN_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
%     [performanceLogRegHbO(1,i2)] = trainClassifierLogisticRegressionHRF_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
%     [performanceSVMHbO(1,i2)] = trainClassifierSVMHRF_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
%     [performanceBaggingHbO(1,i2)] = trainClassifierBaggingTreesSingleHRF_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
%     [performanceBoostedHbO(1,i2)] = trainClassifierBoostedTreesSingleHRF_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieListTrain,movieListTest);
    
end

end