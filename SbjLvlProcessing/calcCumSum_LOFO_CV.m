% Add repts and CI. Used for publication
%
% Each array is channels x time x trial
% This use cross-validation where GLM is computed for each training fold
% and SS beta coefficients from training folds are used for test fold.
%
% Used for subset of channels
% DiffTLen!!!
% baselinePerf is classification accuracy for comparison. Here it is
% all-channel classification for leave-one-feature-out
% This is regular cross-validation.

function [perf_Subset,predLabels_Subset] = ...
    calcCumSum_LOFO_CV(sbjNum,trialsTr,trialsTst,movieListTrain,movieListTest,numClasses,idxChn,mlActAuto,timePt)

rawDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];

fs = 50;
lenT = 0.5*fs;
startT = 2*fs;
zeroT = 2*fs;
thisMLAll = mlActAuto([idxChn idxChn+(length(mlActAuto)/2)]);
numChn = length(thisMLAll)/2;
kFold = 5;

cvp = cvpartition(size(trialsTr,3), 'KFold',kFold);

% Offset
trialsTr = offsetTrials(trialsTr,zeroT);
trialsTst = offsetTrials(trialsTst,zeroT);

predictedLabels = zeros(length(timePt),size(trialsTr,3));
predictedLabelsLOFO = zeros(length(timePt),size(idxChn,2),size(trialsTr,3));
trPerf_LOFO = zeros(length(timePt),size(idxChn,2));
perf_Subset = zeros(1,length(timePt));
predLabels_Subset = zeros(length(timePt),size(trialsTst,3));

for i1 = 1:kFold
    
    trialsInnerTrain = trialsTr(:,:,cvp.training(i1));
    trialsInnerTest = trialsTr(:,:,cvp.test(i1));
    
    for iT = 1:length(timePt)
        startT = timePt(iT);
        if startT == 0
            startT = 1;
        end
    
        temp = cumsum(trialsInnerTrain(idxChn,startT:startT+lenT,:),2);
        if size(temp,1)==1
            cumsumTr = squeeze(temp(:,end,:))';
        else
            cumsumTr = squeeze(temp(:,end,:));
        end

        temp = cumsum(trialsInnerTest(idxChn,startT:startT+lenT,:),2);
        if size(temp,1)==1
            cumsumTst = squeeze(temp(:,end,:))';
        else
            cumsumTst = squeeze(temp(:,end,:));
        end

        % baseline classification for all-channels
        [~,predictedLabels(iT,cvp.test(i1))] = train_RLDA_Ledoit_TrTst(cumsumTr,cumsumTst,thisMLAll,numChn,numClasses,movieListTrain(cvp.training(i1),:),movieListTrain(cvp.test(i1),:));
    
    end
        
    for i2 = 1:size(idxChn,2)

        idxLOFO = idxChn;
        idxLOFO(i2) = [];

        thisTr = trialsInnerTrain(idxLOFO,:,:);
        thisTst = trialsInnerTest(idxLOFO,:,:);
        
        for iT = 1:length(timePt)
            startT = timePt(iT);
            if startT == 0
                startT = 1;
            end

            temp = cumsum(thisTr(:,startT:startT+lenT,:),2);
            if size(temp,1)==1
                cumsumTr = squeeze(temp(:,end,:))';
            else
                cumsumTr = squeeze(temp(:,end,:));
            end

            temp = cumsum(thisTst(:,startT:startT+lenT,:),2);
            if size(temp,1)==1
                cumsumTst = squeeze(temp(:,end,:))';
            else
                cumsumTst = squeeze(temp(:,end,:));
            end

            %[performanceLDALedoitHbO(1,i2)] = ...
            %    train_RLDA_Ledoit_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieIdx,movieListTest);

            thisMLTemp = thisMLAll;
            thisMLTemp([i2 i2+numChn]) = [];

            [~,predictedLabelsLOFO(iT,i2,cvp.test(i1))] = train_RLDA_Ledoit_TrTst(cumsumTr,cumsumTst,thisMLTemp,numChn-1,numClasses,movieListTrain(cvp.training(i1),:),movieListTrain(cvp.test(i1),:));

        end
            
    end
        
end

trLabels_All = squeeze(predictedLabels==movieListTrain(:,2)');
trPerf = sum(trLabels_All,2)/size(trLabels_All,2);
for iT = 1:length(timePt)

    startT = timePt(iT);
    if startT == 0
        startT = 1;
    end    
    
    trLabels_LOFO = squeeze(predictedLabelsLOFO(iT,:,:))==movieListTrain(:,2)';
    trPerf_LOFO(iT,:) = sum(trLabels_LOFO,2)/size(trLabels_LOFO,2);

    idxFeatFilt = idxChn(trPerf(iT) - trPerf_LOFO(iT,:) >= 0);
    
    trialsTrSubset = trialsTr(idxFeatFilt,:,:);
    trialsTstSubset = trialsTst(idxFeatFilt,:,:);

    thisML = mlActAuto;
    thisMLSubset = thisML([idxFeatFilt idxFeatFilt+numChn]);
    numChnSubset = length(thisMLSubset)/2;

    temp = cumsum(trialsTrSubset(:,startT:startT+lenT,:),2);
    if size(temp,1)==1
        cumsumTr = squeeze(temp(:,end,:))';
    else
        cumsumTr = squeeze(temp(:,end,:));
    end

    temp = cumsum(trialsTstSubset(:,startT:startT+lenT,:),2);
    if size(temp,1)==1
        cumsumTst = squeeze(temp(:,end,:))';
    else
        cumsumTst = squeeze(temp(:,end,:));
    end

    [perf_Subset(iT),predLabels_Subset(iT,:)] = train_RLDA_Ledoit_TrTst(cumsumTr,cumsumTst,thisMLSubset,numChnSubset,numClasses,movieListTrain,movieListTest);

end

end