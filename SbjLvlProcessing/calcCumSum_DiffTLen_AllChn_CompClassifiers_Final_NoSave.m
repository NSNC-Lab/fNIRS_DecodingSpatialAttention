% Add repts and CI. Used for publication
%
% Each array is channels x time x trial
% This use cross-validation where GLM is computed for each training fold
% and SS beta coefficients from training folds are used for test fold.

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