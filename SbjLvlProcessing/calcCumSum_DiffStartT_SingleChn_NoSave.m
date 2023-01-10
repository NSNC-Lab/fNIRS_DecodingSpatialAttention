% Add repts and CI. Used for publication
%
% Each array is channels x time x trial
% This use cross-validation where GLM is computed for each training fold
% and SS beta coefficients from training folds are used for test fold.

function performanceLDA = calcCumSum_DiffStartT_SingleChn_NoSave(sbjNum,...
    trialsTr,trialsTst,movieListTrain,movieListTest,numClasses,mlActAuto)

fs = 50;
timePt = (-2*fs:0.25*fs:5*fs);
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

performanceLDA = zeros(size(trialsTr,1),length(timePt));

lenT = 1*fs;

%for i2 = 2:length(timePt)
t_Sample = [2*fs, 3*fs, 4*fs];

for i2 = 1:length(t_Sample)
    thisT = t_Sample(i2);
    %t_Sample = 4*fs;
    tidx = abs(timePt-thisT)<10*eps("double");
    
    temp = cumsum(trialsTr(:,timePt(tidx):timePt(tidx)+lenT,:),2);
    cumsumTr = squeeze(temp(:,end,:));
    
    temp = cumsum(trialsTst(:,timePt(tidx):timePt(tidx)+lenT,:),2);
    cumsumTst = squeeze(temp(:,end,:));
    
    %[performanceLDALedoitHbO(1,i2)] = ...
    %    train_RLDA_Ledoit_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieIdx,movieListTest);
    
    [performanceLDA(:,tidx)] = trainLDA_SingleChannel_TrTst(cumsumTr,cumsumTst,...
        mlList,numChn,numClasses,movieListTrain,movieListTest);
%end

end