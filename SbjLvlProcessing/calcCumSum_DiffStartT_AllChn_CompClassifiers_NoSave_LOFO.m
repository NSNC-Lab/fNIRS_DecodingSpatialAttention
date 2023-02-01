% STATUS: Active(?)
% 
% SYNTAX:
% [performanceLDALedoit] = ...
%   calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_LOFO(sbjNum,...
%   trialsTr,trialsTst,movieListTrain,movieListTest,numClasses,idxChn,...
%   mlActAuto)
% 
% DESCRIPTION:
% Format data and perform leave-one-channel-out classification using cross 
%   validation where SS beta coefficients from training fold is passed to
%   test fold. Use rLDA classifier. Test different decision window lengths.
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
% idxChn - 1D double array: index of channel in 1st dimension of 
%   trialsTr/trialsTst to test
% mlActAuto - 1x1 cell array containing 1D int array of channel list of 2
%   different wavelengths
%
% RETURNED VARIABLES:
% performanceLDALedoit - decoding performance of rLDA classifier for 
%   each decision window. 1D double array: 1 x time windows
% 
% FILES SAVED:
% None.
% 
% PLOTTING:
% None.

function [performanceLDALedoit] = ...
    calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_LOFO(sbjNum,trialsTr,trialsTst,movieListTrain,movieListTest,numClasses,idxChn,mlActAuto)

rawDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];

fs = 50;
lenT = 1.5*fs;
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

performanceLDALedoit = zeros(1,size(idxChn,2));
% performanceLDACERNNHbO = zeros(1,length(timePt));
% performanceLogRegHbO = zeros(1,length(timePt));
% performanceSVMHbO = zeros(1,length(timePt));
% performanceBaggingHbO = zeros(1,length(timePt));
% performanceBoostedHbO = zeros(1,length(timePt));

for i2 = 1:size(idxChn,2)
    
    idxLOFO = idxChn;
    idxLOFO(i2) = [];
    
    thisTr = trialsTr(idxLOFO,:,:);
    thisTst = trialsTst(idxLOFO,:,:);
    
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
    
    thisML = mlActAuto;
    thisML([i2 i2+numChn]) = [];
    
    [performanceLDALedoit(1,i2)] = train_RLDA_Ledoit_TrTst(cumsumTr,cumsumTst,thisML,numChn-1,numClasses,movieListTrain,movieListTest);

end

end