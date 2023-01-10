% Add repts and CI. Used for publication
%
% Each array is channels x time x trial
% This use cross-validation where GLM is computed for each training fold
% and SS beta coefficients from training folds are used for test fold.
%
% Used for subset of channels
% DiffTLen!!!

function [performanceLDALedoitHbO] = ...
    calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trialsTr,trialsTst,movieListTrain,movieListTest,numClasses,mlActAuto)

rawDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];

fs = 50;
timeLen = [0.1 0.2 0.5 1 1.5 2 3 4 5].*fs;
%startT = 4*fs;
startT = 2*fs;
zeroT = 2*fs;
numChn = length(mlActAuto)/2;

% if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
%     isSDNew = 0;
% else
%     isSDNew = 1;
% end

%figure(1);hold on;
%figure('units','normalized','outerposition',[0 0 1 1]);hold on;

% if strcmp(sbjNum,'15')
%     movieListTrain = movieListTrain(2:end,:);
% end

% Offset
trialsTr = offsetTrials(trialsTr,zeroT);
trialsTst = offsetTrials(trialsTst,zeroT);

% after convert2SD2, cut down from 42 chns to 36 chns. Remove old chns
% if ~isSDNew
%     [trialsTr,mlList] = convert2SD2(trialsTr,mlActAuto{1});
%     [trialsTst] = convert2SD2(trialsTst,mlActAuto{1});
% else
%     mlList = mlActAuto{1};
% end

%mlList = mlActAuto;

% after selectRS, cut down from 36 to 30 chns. Remove SS. Both steps can be
% implemented simultaneously using bit logic.
% [trialsTr] = selectRS(trialsTr,1,mlList);
% [trialsTst,mlList] = selectRS(trialsTst,1,mlList);

% behFN = [rawDataDir filesep 'responses_' sbjNum];

% [trialsTr,movieIdx] = keepCorrectTrials_Special(sbjNum,trialsTr,behFN,numClasses,movieListTrain);
% [trialsTst,movieListTest] = keepCorrectTrials_Special(sbjNum,trialsTst,behFN,numClasses,movieListTest);

performanceLDALedoitHbO = zeros(1,length(timeLen));
% performanceLDACERNNHbO = zeros(1,length(timePt));
% performanceLogRegHbO = zeros(1,length(timePt));
% performanceSVMHbO = zeros(1,length(timePt));
% performanceBaggingHbO = zeros(1,length(timePt));
% performanceBoostedHbO = zeros(1,length(timePt));

lenT = 1*fs;

for i2 = 1:length(timeLen)
    
    temp = cumsum(trialsTr(:,startT:startT+timeLen(i2),:),2);
    if size(temp,1)==1
        cumsumTr = squeeze(temp(:,end,:))';
    else
        cumsumTr = squeeze(temp(:,end,:));
    end
    
    temp = cumsum(trialsTst(:,startT:startT+timeLen(i2),:),2);
    if size(temp,1)==1
        cumsumTst = squeeze(temp(:,end,:))';
    else
        cumsumTst = squeeze(temp(:,end,:));
    end
    
    %[performanceLDALedoitHbO(1,i2)] = ...
    %    train_RLDA_Ledoit_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieIdx,movieListTest);
    
    [performanceLDALedoitHbO(1,i2)] = train_RLDA_Ledoit_TrTst(cumsumTr,cumsumTst,mlActAuto,numChn,numClasses,movieListTrain,movieListTest);
%     [performanceLDACERNNHbO(1,i2)] = train_RLDA_CERNN_TrTst(cumsumTr,cumsumTst,mlActAuto,numChn,numClasses,movieListTrain,movieListTest);
%     [performanceLogRegHbO(1,i2)] = trainClassifierLogisticRegressionHRF_TrTst(cumsumTr,cumsumTst,mlActAuto,numChn,numClasses,movieListTrain,movieListTest);
%     [performanceSVMHbO(1,i2)] = trainClassifierSVMHRF_TrTst(cumsumTr,cumsumTst,mlActAuto,numChn,numClasses,movieListTrain,movieListTest);
%     [performanceBaggingHbO(1,i2)] = trainClassifierBaggingTreesSingleHRF_TrTst(cumsumTr,cumsumTst,mlActAuto,numChn,numClasses,movieListTrain,movieListTest);
%     [performanceBoostedHbO(1,i2)] = trainClassifierBoostedTreesSingleHRF_TrTst(cumsumTr,cumsumTst,mlActAuto,numChn,numClasses,movieListTrain,movieListTest);
    
end

end