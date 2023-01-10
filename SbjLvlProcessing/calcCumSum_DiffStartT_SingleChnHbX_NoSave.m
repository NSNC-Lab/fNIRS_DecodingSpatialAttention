% Add repts and CI. Used for publication
%
% Each array is channels x time x trial
% This use cross-validation where GLM is computed for each training fold
% and SS beta coefficients from training folds are used for test fold.

function performanceLDA = calcCumSum_DiffStartT_SingleChnHbX_NoSave(sbjNum,...
    trialsTrHbO,trialsTstHbO,trialsTrHbR,trialsTstHbR,trialsTrHbT,trialsTstHbT,...
    movieListTrain,movieListTest,numClasses,mlActAuto)

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
trialsTrHbO = offsetTrials(trialsTrHbO,zeroT);
trialsTstHbO = offsetTrials(trialsTstHbO,zeroT);

trialsTrHbR = offsetTrials(trialsTrHbR,zeroT);
trialsTstHbR = offsetTrials(trialsTstHbR,zeroT);

trialsTrHbT = offsetTrials(trialsTrHbT,zeroT);
trialsTstHbT = offsetTrials(trialsTstHbT,zeroT);

% after convert2SD2, cut down from 42 chns to 36 chns. Remove old chns
if ~isSDNew
    [trialsTrHbO,mlList] = convert2SD2(trialsTrHbO,mlActAuto{1});
    [trialsTstHbO] = convert2SD2(trialsTstHbO,mlActAuto{1});
    
    [trialsTrHbR] = convert2SD2(trialsTrHbR,mlActAuto{1});
    [trialsTstHbR] = convert2SD2(trialsTstHbR,mlActAuto{1});
    
    [trialsTrHbT] = convert2SD2(trialsTrHbT,mlActAuto{1});
    [trialsTstHbT] = convert2SD2(trialsTstHbT,mlActAuto{1});
else
    mlList = mlActAuto{1};
end

% after selectRS, cut down from 36 to 30 chns. Remove SS. Both steps can be
% implemented simultaneously using bit logic.
[trialsTrHbO] = selectRS(trialsTrHbO,1,mlList);
[trialsTstHbO] = selectRS(trialsTstHbO,1,mlList);

[trialsTrHbR] = selectRS(trialsTrHbR,1,mlList);
[trialsTstHbR] = selectRS(trialsTstHbR,1,mlList);

[trialsTrHbT] = selectRS(trialsTrHbT,1,mlList);
[trialsTstHbT,mlList] = selectRS(trialsTstHbT,1,mlList);

% behFN = [rawDataDir filesep 'responses_' sbjNum];

% [trialsTr,movieIdx] = keepCorrectTrials_Special(sbjNum,trialsTr,behFN,numClasses,movieListTrain);
% [trialsTst,movieListTest] = keepCorrectTrials_Special(sbjNum,trialsTst,behFN,numClasses,movieListTest);

performanceLDA = zeros(size(trialsTrHbO,1),length(timePt));

lenT = 1*fs;

%for i2 = 2:length(timePt)
t_Sample = [2*fs, 3*fs, 4*fs];

for i2 = 1:length(t_Sample)
    thisT = t_Sample(i2);
    %t_Sample = 4*fs;
    tidx = abs(timePt-thisT)<10*eps("double");
    
    temp = cumsum(trialsTrHbO(:,timePt(tidx):timePt(tidx)+lenT,:),2);
    % % trials are channels x trial (p x n)
    cumsumHbOTr = squeeze(temp(:,end,:));
    
    temp = cumsum(trialsTstHbO(:,timePt(tidx):timePt(tidx)+lenT,:),2);
    cumsumHbOTst = squeeze(temp(:,end,:));
    
    temp = cumsum(trialsTrHbR(:,timePt(tidx):timePt(tidx)+lenT,:),2);
    cumsumHbRTr = squeeze(temp(:,end,:));
    
    temp = cumsum(trialsTstHbR(:,timePt(tidx):timePt(tidx)+lenT,:),2);
    cumsumHbRTst = squeeze(temp(:,end,:));
    
    temp = cumsum(trialsTrHbT(:,timePt(tidx):timePt(tidx)+lenT,:),2);
    cumsumHbTTr = squeeze(temp(:,end,:));
    
    temp = cumsum(trialsTstHbT(:,timePt(tidx):timePt(tidx)+lenT,:),2);
    cumsumHbTTst = squeeze(temp(:,end,:));
    
    
    
    %[performanceLDALedoitHbO(1,i2)] = ...
    %    train_RLDA_Ledoit_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,movieIdx,movieListTest);
    
    for i3 = 1:numChn
        
        % check dimension
        cumsumTr = [cumsumHbOTr(i3,:); cumsumHbRTr(i3,:); cumsumHbTTr(i3,:)];
        cumsumTst = [cumsumHbOTst(i3,:); cumsumHbRTst(i3,:); cumsumHbTTst(i3,:)];
        
        [performanceLDA(i3,tidx)] = trainLDA_SingleChannelHbX_TrTst(cumsumTr,cumsumTst,...
            numClasses,movieListTrain,movieListTest);
    end
%end

end