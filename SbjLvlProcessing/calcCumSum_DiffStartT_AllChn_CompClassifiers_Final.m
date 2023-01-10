% Add repts and CI. Used for publication
%
% For sbj 08 and 10.
%
% All-channel classification

function calcCumSum_DiffStartT_AllChn_CompClassifiers_Final(sbjNum,numClasses,saveOp)
processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
rawDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
figSaveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' ...
    num2str(sbjNum) '\Figures\Classification'];

fs = 50;
timePt = (0:0.25*fs:5*fs)+2*fs;
zeroT = 2*fs;
numChn = 30;
nrept = 10;
kFold = 5;
cmap = jet(numChn);
if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
    isSDNew = 0;
else
    isSDNew = 1;
end

chnName = getChnName(numChn);

%figure(1);hold on;
figure('units','normalized','outerposition',[0 0 1 1]);hold on;

% Too lazy to modularize this
% Basis 1
%load([processedDataDir filesep 'intermediateOutputsCombined_Basis1.mat']);
load([processedDataDir filesep 'intermediateOutputsCombined_Basis1_SSBetaPrior.mat']);
% Each array is channels x time x trial
if numClasses == 2
    load([processedDataDir filesep 'singleTrialsUpdated_SSBeta_LR.mat'],'singleTrialHRFHbOM',...
        'singleTrialHRFHbRM',...
        'singleTrialHRFHbTM','indexMoviesTest');
    
    yLimAxis = [0 1];
else
    load([processedDataDir filesep 'singleTrialsUpdated_SSBeta.mat'],'singleTrialHRFHbOM',...
        'singleTrialHRFHbRM',...
        'singleTrialHRFHbTM','indexMoviesTest');
    yLimAxis = [0 1];
end

if strcmp(sbjNum,'15')
    indexMoviesTest = indexMoviesTest(2:end,:);
end

% Offset
singleTrialHRFHbOM = offsetTrials(singleTrialHRFHbOM,zeroT);
singleTrialHRFHbRM = offsetTrials(singleTrialHRFHbRM,zeroT);
singleTrialHRFHbTM = offsetTrials(singleTrialHRFHbTM,zeroT);

% after convert2SD2, cut down from 42 chns to 36 chns. Remove old chns
if ~isSDNew
    [singleTrialHRFHbOM,mlList] = convert2SD2(singleTrialHRFHbOM,mlActAuto{1});
    [singleTrialHRFHbRM,~] = convert2SD2(singleTrialHRFHbRM);
    [singleTrialHRFHbTM,~] = convert2SD2(singleTrialHRFHbTM);
else
    mlList = mlActAuto{1};
end

% after selectRS, cut down from 36 to 30 chns. Remove SS. Both steps can be
% implemented simultaneously using bit logic.
[singleTrialHRFHbOM,mlList] = selectRS(singleTrialHRFHbOM,1,mlList);
[singleTrialHRFHbRM,~] = selectRS(singleTrialHRFHbRM,1);
[singleTrialHRFHbTM,~] = selectRS(singleTrialHRFHbTM,1);

behFN = [rawDataDir filesep 'responses_' sbjNum];

[singleTrialHRFHbOM,movieIdx,behScore] = keepCorrectTrials_Special(sbjNum,singleTrialHRFHbOM,behFN,numClasses,indexMoviesTest);
[singleTrialHRFHbRM] = keepCorrectTrials_Special(sbjNum,singleTrialHRFHbRM,behFN,numClasses,indexMoviesTest);
[singleTrialHRFHbTM] = keepCorrectTrials_Special(sbjNum,singleTrialHRFHbTM,behFN,numClasses,indexMoviesTest);

performanceLDALedoitHbO = zeros(length(timePt),nrept);
performanceLDALedoitHbR = zeros(length(timePt),nrept);
performanceLDALedoitHbT = zeros(length(timePt),nrept);

performanceLDACERNNHbO = zeros(length(timePt),nrept);
performanceLDACERNNHbR = zeros(length(timePt),nrept);
performanceLDACERNNHbT = zeros(length(timePt),nrept);

% performanceQDALedoitHbO = zeros(length(timePt),nrept);
% performanceQDALedoitHbR = zeros(length(timePt),nrept);
% performanceQDALedoitHbT = zeros(length(timePt),nrept);
% 
% performanceLDAGDHbO = zeros(length(timePt),nrept);
% performanceLDAGDHbR = zeros(length(timePt),nrept);
% performanceLDAGDHbT = zeros(length(timePt),nrept);

performanceLogRegHbO = zeros(length(timePt),nrept);
performanceLogRegHbR = zeros(length(timePt),nrept);
performanceLogRegHbT = zeros(length(timePt),nrept);

performanceSVMHbO = zeros(length(timePt),nrept);
performanceSVMHbR = zeros(length(timePt),nrept);
performanceSVMHbT = zeros(length(timePt),nrept);

performanceBaggingHbO = zeros(length(timePt),nrept);
performanceBaggingHbR = zeros(length(timePt),nrept);
performanceBaggingHbT = zeros(length(timePt),nrept);

performanceBoostedHbO = zeros(length(timePt),nrept);
performanceBoostedHbR = zeros(length(timePt),nrept);
performanceBoostedHbT = zeros(length(timePt),nrept);

performanceLDALedoitHbOFold = zeros(length(timePt),nrept,kFold);
performanceLDALedoitHbRFold = zeros(length(timePt),nrept,kFold);
performanceLDALedoitHbTFold = zeros(length(timePt),nrept,kFold);

performanceLDACERNNHbOFold = zeros(length(timePt),nrept,kFold);
performanceLDACERNNHbRFold = zeros(length(timePt),nrept,kFold);
performanceLDACERNNHbTFold = zeros(length(timePt),nrept,kFold);

% performanceQDALedoitHbOFold = zeros(length(timePt),nrept,kFold);
% performanceQDALedoitHbRFold = zeros(length(timePt),nrept,kFold);
% performanceQDALedoitHbTFold = zeros(length(timePt),nrept,kFold);
% 
% performanceLDAGDHbOFold = zeros(length(timePt),nrept,kFold);
% performanceLDAGDHbRFold = zeros(length(timePt),nrept,kFold);
% performanceLDAGDHbTFold = zeros(length(timePt),nrept,kFold);

performanceLogRegHbOFold = zeros(length(timePt),nrept,kFold);
performanceLogRegHbRFold = zeros(length(timePt),nrept,kFold);
performanceLogRegHbTFold = zeros(length(timePt),nrept,kFold);

performanceSVMHbOFold = zeros(length(timePt),nrept,kFold);
performanceSVMHbRFold = zeros(length(timePt),nrept,kFold);
performanceSVMHbTFold = zeros(length(timePt),nrept,kFold);

performanceBaggingHbOFold = zeros(length(timePt),nrept,kFold);
performanceBaggingHbRFold = zeros(length(timePt),nrept,kFold);
performanceBaggingHbTFold = zeros(length(timePt),nrept,kFold);

performanceBoostedHbOFold = zeros(length(timePt),nrept,kFold);
performanceBoostedHbRFold = zeros(length(timePt),nrept,kFold);
performanceBoostedHbTFold = zeros(length(timePt),nrept,kFold);

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


% mean
performanceLDALedoitHbOMean = zeros(length(timePt),1);
performanceLDALedoitHbRMean = zeros(length(timePt),1);
performanceLDALedoitHbTMean = zeros(length(timePt),1);

performanceLDACERNNHbOMean = zeros(length(timePt),1);
performanceLDACERNNHbRMean = zeros(length(timePt),1);
performanceLDACERNNHbTMean = zeros(length(timePt),1);

performanceLogRegHbOMean = zeros(length(timePt),1);
performanceLogRegHbRMean = zeros(length(timePt),1);
performanceLogRegHbTMean = zeros(length(timePt),1);

performanceSVMHbOMean = zeros(length(timePt),1);
performanceSVMHbRMean = zeros(length(timePt),1);
performanceSVMHbTMean = zeros(length(timePt),1);

performanceBaggingHbOMean = zeros(length(timePt),1);
performanceBaggingHbRMean = zeros(length(timePt),1);
performanceBaggingHbTMean = zeros(length(timePt),1);

performanceBoostedHbOMean = zeros(length(timePt),1);
performanceBoostedHbRMean = zeros(length(timePt),1);
performanceBoostedHbTMean = zeros(length(timePt),1);


lenT = 1*fs;

for i2 = 2:length(timePt)
    
    temp = cumsum(singleTrialHRFHbOM(:,timePt(i2):timePt(i2)+lenT,:),2);
    cumsumHbOM = squeeze(temp(:,end,:));
    
%     temp = cumsum(singleTrialHRFHbRS(:,timePt(i):timePt(i)+lenT,:),2);
%     cumsumHbRS = squeeze(temp(:,end,:));
    
    temp = cumsum(singleTrialHRFHbRM(:,timePt(i2):timePt(i2)+lenT,:),2);
    cumsumHbRM = squeeze(temp(:,end,:));
    
%     temp = cumsum(singleTrialHRFHbTS(:,startT:timePt(i),:),2);
%     cumsumHbTS = squeeze(temp(:,end,:));
    
    temp = cumsum(singleTrialHRFHbTM(:,timePt(i2):timePt(i2)+lenT,:),2);
    cumsumHbTM = squeeze(temp(:,end,:));

%     [performanceArrMultiHbO_95(i),predRespMHbO] = trainClassifierLinearDiscriminantSingleHRFPCA(log(cumsumHbOM),0.95,mlList,numChn,numClasses,movieIdx);
%     [performanceArrMultiHbR_95(i),predRespMHbR] = trainClassifierLinearDiscriminantSingleHRFPCA(log(cumsumHbRM),0.95,mlList,numChn,numClasses,movieIdx);
%     [performanceArrMultiHbT_95(i),predRespMHbT] = trainClassifierLinearDiscriminantSingleHRFPCA(log(cumsumHbTM),0.95,mlList,numChn,numClasses,movieIdx);
%        
%     [performanceArrMultiHbO_85(i),predRespMHbO] = trainClassifierLinearDiscriminantSingleHRFPCA(log(cumsumHbOM),0.85,mlList,numChn,numClasses,movieIdx);
%     [performanceArrMultiHbR_85(i),predRespMHbR] = trainClassifierLinearDiscriminantSingleHRFPCA(log(cumsumHbRM),0.85,mlList,numChn,numClasses,movieIdx);
%     [performanceArrMultiHbT_85(i),predRespMHbT] = trainClassifierLinearDiscriminantSingleHRFPCA(log(cumsumHbTM),0.85,mlList,numChn,numClasses,movieIdx);
%        
%     [performanceArrMultiHbO_75(i),predRespMHbO] = trainClassifierLinearDiscriminantSingleHRFPCA(log(cumsumHbOM),0.75,mlList,numChn,numClasses,movieIdx);
%     [performanceArrMultiHbR_75(i),predRespMHbR] = trainClassifierLinearDiscriminantSingleHRFPCA(log(cumsumHbRM),0.75,mlList,numChn,numClasses,movieIdx);
%     [performanceArrMultiHbT_75(i),predRespMHbT] = trainClassifierLinearDiscriminantSingleHRFPCA(log(cumsumHbTM),0.75,mlList,numChn,numClasses,movieIdx);
%     
%     [performanceArrMultiHbO_65(i),predRespMHbO] = trainClassifierLinearDiscriminantSingleHRFPCA(log(cumsumHbOM),0.65,mlList,numChn,numClasses,movieIdx);
%     [performanceArrMultiHbR_65(i),predRespMHbR] = trainClassifierLinearDiscriminantSingleHRFPCA(log(cumsumHbRM),0.65,mlList,numChn,numClasses,movieIdx);
%     [performanceArrMultiHbT_65(i),predRespMHbT] = trainClassifierLinearDiscriminantSingleHRFPCA(log(cumsumHbTM),0.65,mlList,numChn,numClasses,movieIdx);
%        
%     [performanceArrMultiHbO_55(i),predRespMHbO] = trainClassifierLinearDiscriminantSingleHRFPCA(log(cumsumHbOM),0.55,mlList,numChn,numClasses,movieIdx);
%     [performanceArrMultiHbR_55(i),predRespMHbR] = trainClassifierLinearDiscriminantSingleHRFPCA(log(cumsumHbRM),0.55,mlList,numChn,numClasses,movieIdx);
%     [performanceArrMultiHbT_55(i),predRespMHbT] = trainClassifierLinearDiscriminantSingleHRFPCA(log(cumsumHbTM),0.55,mlList,numChn,numClasses,movieIdx);
    
    for i3 = 1:nrept

        [performanceLDALedoitHbO(i2,i3), performanceLDALedoitHbOFold(i2,i3,:)] = train_RLDA_Ledoit(cumsumHbOM,mlList,numChn,numClasses,movieIdx);
        [performanceLDALedoitHbR(i2,i3), performanceLDALedoitHbRFold(i2,i3,:)] = train_RLDA_Ledoit(cumsumHbRM,mlList,numChn,numClasses,movieIdx);
        [performanceLDALedoitHbT(i2,i3), performanceLDALedoitHbTFold(i2,i3,:)] = train_RLDA_Ledoit(cumsumHbTM,mlList,numChn,numClasses,movieIdx);

        [performanceLDACERNNHbO(i2,i3),performanceLDACERNNHbOFold(i2,i3,:)] = train_RLDA_CERNN(cumsumHbOM,mlList,numChn,numClasses,movieIdx);
        [performanceLDACERNNHbR(i2,i3),performanceLDACERNNHbRFold(i2,i3,:)] = train_RLDA_CERNN(cumsumHbRM,mlList,numChn,numClasses,movieIdx);
        [performanceLDACERNNHbT(i2,i3),performanceLDACERNNHbTFold(i2,i3,:)] = train_RLDA_CERNN(cumsumHbTM,mlList,numChn,numClasses,movieIdx);

    %     performanceQDALedoitHbO(i) = train_RQDA_Ledoit(cumsumHbOM,mlList,numChn,numClasses,movieIdx);
    %     performanceQDALedoitHbR(i) = train_RQDA_Ledoit(cumsumHbRM,mlList,numChn,numClasses,movieIdx);
    %     performanceQDALedoitHbT(i) = train_RQDA_Ledoit(cumsumHbTM,mlList,numChn,numClasses,movieIdx);

    %     performanceLDAGDHbO(i) = train_RLDA_Matlab(cumsumHbOM,mlList,numChn,numClasses,movieIdx);
    %     performanceLDAGDHbR(i) = train_RLDA_Matlab(cumsumHbRM,mlList,numChn,numClasses,movieIdx);
    %     performanceLDAGDHbT(i) = train_RLDA_Matlab(cumsumHbTM,mlList,numChn,numClasses,movieIdx);

        [performanceLogRegHbO(i2,i3),performanceLogRegHbOFold(i2,i3,:)] = trainClassifierLogisticRegressionHRF(cumsumHbOM,mlList,numChn,numClasses,movieIdx);
        [performanceLogRegHbR(i2,i3),performanceLogRegHbRFold(i2,i3,:)] = trainClassifierLogisticRegressionHRF(cumsumHbRM,mlList,numChn,numClasses,movieIdx);
        [performanceLogRegHbT(i2,i3),performanceLogRegHbTFold(i2,i3,:)] = trainClassifierLogisticRegressionHRF(cumsumHbTM,mlList,numChn,numClasses,movieIdx);


    %     [performanceArrMultiHbO(i),predRespMHbO] = trainClassifierLinearDiscriminantSingleHRFAllChns(cumsumHbOM,mlActAuto{1},numClasses,movieIdx);
    %     [performanceArrMultiHbR(i),predRespMHbR] = trainClassifierLinearDiscriminantSingleHRFAllChns(cumsumHbRM,mlActAuto{1},numClasses,movieIdx);
    %     [performanceArrMultiHbT(i),predRespMHbT] = trainClassifierLinearDiscriminantSingleHRFAllChns(cumsumHbTM,mlActAuto{1},numClasses,movieIdx);

        [performanceSVMHbO(i2,i3),performanceSVMHbOFold(i2,i3,:)] = trainClassifierSVMHRF(cumsumHbOM,mlList,numChn,numClasses,movieIdx);
        [performanceSVMHbR(i2,i3),performanceSVMHbRFold(i2,i3,:)] = trainClassifierSVMHRF(cumsumHbRM,mlList,numChn,numClasses,movieIdx);
        [performanceSVMHbT(i2,i3),performanceSVMHbTFold(i2,i3,:)] = trainClassifierSVMHRF(cumsumHbTM,mlList,numChn,numClasses,movieIdx);

        [performanceBaggingHbO(i2,i3),performanceBaggingHbOFold(i2,i3,:)] = trainClassifierBaggingTreesSingleHRF(cumsumHbOM,mlList,numChn,numClasses,movieIdx);
        [performanceBaggingHbR(i2,i3),performanceBaggingHbRFold(i2,i3,:)] = trainClassifierBaggingTreesSingleHRF(cumsumHbRM,mlList,numChn,numClasses,movieIdx);
        [performanceBaggingHbT(i2,i3),performanceBaggingHbTFold(i2,i3,:)] = trainClassifierBaggingTreesSingleHRF(cumsumHbTM,mlList,numChn,numClasses,movieIdx);

        [performanceBoostedHbO(i2,i3),performanceBoostedHbOFold(i2,i3,:)] = trainClassifierBoostedTreesSingleHRF(cumsumHbOM,mlList,numChn,numClasses,movieIdx);
        [performanceBoostedHbR(i2,i3),performanceBoostedHbRFold(i2,i3,:)] = trainClassifierBoostedTreesSingleHRF(cumsumHbRM,mlList,numChn,numClasses,movieIdx);
        [performanceBoostedHbT(i2,i3),performanceBoostedHbTFold(i2,i3,:)] = trainClassifierBoostedTreesSingleHRF(cumsumHbTM,mlList,numChn,numClasses,movieIdx);
    end
end

if numClasses == 2
    %savePerfFN = 'performanceLinearDiscriminantUpdated_LR.mat';
    savePerfFN = 'performanceLinearDiscriminantUpdated_SSBeta_LR.mat';
else
    %savePerfFN = 'performanceLinearDiscriminantUpdated.mat';
    savePerfFN = 'performanceLinearDiscriminantUpdated_SSBeta.mat';
end
save([processedDataDir filesep savePerfFN]);

conf = 0.95;

for i = 1:length(timePt)
    
    temp = performanceLDALedoitHbOFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDALedoitHbOMean(i) = mean(temp(:));
    performanceLDALedoitHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceLDALedoitHbRFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDALedoitHbRMean(i) = mean(temp(:));
    performanceLDALedoitHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceLDALedoitHbTFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDALedoitHbTMean(i) = mean(temp(:));
    performanceLDALedoitHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % CERNN
    temp = performanceLDACERNNHbOFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDACERNNHbOMean(i) = mean(temp(:));
    performanceLDACERNNHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceLDACERNNHbRFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDACERNNHbRMean(i) = mean(temp(:));
    performanceLDACERNNHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceLDACERNNHbTFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDACERNNHbTMean(i) = mean(temp(:));
    performanceLDACERNNHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % Log
    temp = performanceLogRegHbOFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLogRegHbOMean(i) = mean(temp(:));
    performanceLogRegHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceLogRegHbRFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLogRegHbRMean(i) = mean(temp(:));
    performanceLogRegHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceLogRegHbTFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLogRegHbTMean(i) = mean(temp(:));
    performanceLogRegHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % SVM
    temp = performanceSVMHbOFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceSVMHbOMean(i) = mean(temp(:));
    performanceSVMHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceSVMHbRFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceSVMHbRMean(i) = mean(temp(:));
    performanceSVMHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceSVMHbTFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceSVMHbTMean(i) = mean(temp(:));
    performanceSVMHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % Bagging
    temp = performanceBaggingHbOFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceBaggingHbOMean(i) = mean(temp(:));
    performanceBaggingHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceBaggingHbRFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceBaggingHbRMean(i) = mean(temp(:));
    performanceBaggingHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceBaggingHbTFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceBaggingHbTMean(i) = mean(temp(:));
    performanceBaggingHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % Boosting
    temp = performanceBoostedHbOFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceBoostedHbOMean(i) = mean(temp(:));
    performanceBoostedHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceBoostedHbRFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceBoostedHbRMean(i) = mean(temp(:));
    performanceBoostedHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performanceBoostedHbTFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceBoostedHbTMean(i) = mean(temp(:));
    performanceBoostedHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
end

numClassifiers = 6;
cmap = jet(numClassifiers);

subplot(1,3,1);hold on;
% for j = 1:size(performanceArrMultiHbO,1)
%     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
%     plot(timePt./fs-2,performanceArrMultiHbO(j,:),'Color',cmap(j,:));
% end
% plot(timePt./fs-2,performanceArrMultiHbO_95,'Color',cmap(1,:));
% plot(timePt./fs-2,performanceArrMultiHbO_85,'Color',cmap(2,:));
% plot(timePt./fs-2,performanceArrMultiHbO_75,'Color',cmap(3,:));
% plot(timePt./fs-2,performanceArrMultiHbO_65,'Color',cmap(4,:));
% plot(timePt./fs-2,performanceArrMultiHbO_55,'Color',cmap(5,:));
% plot(timePt./fs-2,performanceLDALedoitHbO,'Color',cmap(6,:));
% plot(timePt./fs-2,performanceLDACERNNHbO,'Color',cmap(7,:));
% % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
% % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
% plot(timePt./fs-2,performanceLogRegHbO,'Color',cmap(8,:));
% plot(timePt./fs-2,performanceSVMHbO,'Color',cmap(9,:));

errorbar(timePt./fs-2,performanceLDALedoitHbOMean,performanceLDALedoitHbOCI(:,2),'Color',cmap(1,:));
errorbar(timePt./fs-2,performanceLDACERNNHbOMean,performanceLDACERNNHbOCI(:,2),'Color',cmap(2,:));
% plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
% plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
errorbar(timePt./fs-2,performanceLogRegHbOMean,performanceLogRegHbOCI(:,2),'Color',cmap(3,:));
errorbar(timePt./fs-2,performanceSVMHbOMean,performanceSVMHbOCI(:,2),'Color',cmap(4,:));
errorbar(timePt./fs-2,performanceBaggingHbOMean,performanceBaggingHbOCI(:,2),'Color',cmap(5,:));
errorbar(timePt./fs-2,performanceBoostedHbOMean,performanceBoostedHbOCI(:,2),'Color',cmap(6,:));

ylim(yLimAxis);
xlim([0 7]);
title(sprintf('Sbj %s: Δ[HbO] Multi',num2str(sbjNum)));
% legend({'PCA 95%','PCA 85%','PCA 75%','PCA 65%','PCA 55%','LDAShrink',...
%     'QDAShrink','LDAGD','LogReg','SVM'});
% legend({'PCA 95%','PCA 85%','PCA 75%','PCA 65%','PCA 55%','LDAShrink',...
%     'LDA CERNN','LogReg','SVM'});
legend({'LDAShrink','LDA CERNN','LogReg','SVM','Bagging','Boosted'});
xlabel('Time [s]');ylabel('Accuracy');
hold off;

subplot(1,3,2);hold on;
% for j = 1:size(performanceArrMultiHbR,1)
%     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
%     plot(timePt./fs-2,performanceArrMultiHbR(j,:),'Color',cmap(j,:));
% end
% plot(timePt./fs-2,performanceArrMultiHbR_95,'Color',cmap(1,:));
% plot(timePt./fs-2,performanceArrMultiHbR_85,'Color',cmap(2,:));
% plot(timePt./fs-2,performanceArrMultiHbR_75,'Color',cmap(3,:));
% plot(timePt./fs-2,performanceArrMultiHbR_65,'Color',cmap(4,:));
% plot(timePt./fs-2,performanceArrMultiHbR_55,'Color',cmap(5,:));
% plot(timePt./fs-2,performanceLDALedoitHbR,'Color',cmap(6,:));
% plot(timePt./fs-2,performanceLDACERNNHbR,'Color',cmap(7,:));
% % plot(timePt./fs-2,performanceQDALedoitHbR,'Color',cmap(7,:));
% % plot(timePt./fs-2,performanceLDAGDHbR,'Color',cmap(8,:));
% plot(timePt./fs-2,performanceLogRegHbR,'Color',cmap(8,:));
% plot(timePt./fs-2,performanceSVMHbR,'Color',cmap(9,:));

errorbar(timePt./fs-2,performanceLDALedoitHbRMean,performanceLDALedoitHbRCI(:,2),'Color',cmap(1,:));
errorbar(timePt./fs-2,performanceLDACERNNHbRMean,performanceLDACERNNHbRCI(:,2),'Color',cmap(2,:));
% plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
% plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
errorbar(timePt./fs-2,performanceLogRegHbRMean,performanceLogRegHbRCI(:,2),'Color',cmap(3,:));
errorbar(timePt./fs-2,performanceSVMHbRMean,performanceSVMHbRCI(:,2),'Color',cmap(4,:));
errorbar(timePt./fs-2,performanceBaggingHbRMean,performanceBaggingHbRCI(:,2),'Color',cmap(5,:));
errorbar(timePt./fs-2,performanceBoostedHbRMean,performanceBoostedHbRCI(:,2),'Color',cmap(6,:));

ylim(yLimAxis);
xlim([0 7]);
title(sprintf('Sbj %s: Δ[HbR] Multi',num2str(sbjNum)));
%legend('S1D1','S1D2','S1D3','S1D4');
xlabel('Time [s]');ylabel('Accuracy');
hold off;

subplot(1,3,3);hold on;
% for j = 1:size(performanceArrMultiHbT,1)
%     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
%     plot(timePt./fs-2,performanceArrMultiHbT(j,:),'Color',cmap(j,:));
% end
% plot(timePt./fs-2,performanceArrMultiHbT_95,'Color',cmap(1,:));
% plot(timePt./fs-2,performanceArrMultiHbT_85,'Color',cmap(2,:));
% plot(timePt./fs-2,performanceArrMultiHbT_75,'Color',cmap(3,:));
% plot(timePt./fs-2,performanceArrMultiHbT_65,'Color',cmap(4,:));
% plot(timePt./fs-2,performanceArrMultiHbT_55,'Color',cmap(5,:));
% plot(timePt./fs-2,performanceLDALedoitHbT,'Color',cmap(6,:));
% plot(timePt./fs-2,performanceLDACERNNHbT,'Color',cmap(7,:));
% % plot(timePt./fs-2,performanceQDALedoitHbT,'Color',cmap(7,:));
% % plot(timePt./fs-2,performanceLDAGDHbT,'Color',cmap(8,:));
% plot(timePt./fs-2,performanceLogRegHbT,'Color',cmap(8,:));
% plot(timePt./fs-2,performanceSVMHbT,'Color',cmap(9,:));

errorbar(timePt./fs-2,performanceLDALedoitHbTMean,performanceLDALedoitHbTCI(:,2),'Color',cmap(1,:));
errorbar(timePt./fs-2,performanceLDACERNNHbTMean,performanceLDACERNNHbTCI(:,2),'Color',cmap(2,:));
% plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
% plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
errorbar(timePt./fs-2,performanceLogRegHbTMean,performanceLogRegHbTCI(:,2),'Color',cmap(3,:));
errorbar(timePt./fs-2,performanceSVMHbTMean,performanceSVMHbTCI(:,2),'Color',cmap(4,:));
errorbar(timePt./fs-2,performanceBaggingHbTMean,performanceBaggingHbTCI(:,2),'Color',cmap(5,:));
errorbar(timePt./fs-2,performanceBoostedHbTMean,performanceBoostedHbTCI(:,2),'Color',cmap(6,:));

ylim(yLimAxis);
xlim([0 7]);
title(sprintf('Sbj %s: Δ[HbT] Multi',num2str(sbjNum)));
%legend('S1D1','S1D2','S1D3','S1D4');
xlabel('Time [s]');ylabel('Accuracy');
hold off;

% subplot(1,4,4);hold on;
% for j = 1:size(performanceArrMultiHbT,1)
%     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
%     plot(timePt./fs-2,zeros(1,length(timePt./fs-2)),'Color',cmap(j,:));
% end
% %legend(chnName(1:30));
% annotation('textbox', [0.8, 0.05, 0.1, 0.1], 'String', ['Score: ' num2str(behScore) '%'],'FontSize',8);
% annotation('textbox', [0.71, 0.2, 0.1, 0.1],'String','calcCumSum\_MultiOnly\_DiffStartT\_AllChn\_AllBasis.m','FontSize',8);
% title(sprintf('All Channels. Diff Start T: %s-Class',num2str(numClasses)));
% hold off;

if numClasses == 2
    %fn = sprintf('PerformanceCumsumVsTime_DiffStartT_LR_AllChns_DiffClassifiers');
    fn = sprintf('PerformanceCumsumVsTime_DiffStartT_SSBeta_LR_AllChns_DiffClassifiers');
else
    %fn = sprintf('PerformanceCumsumVsTime_DiffStartT_AllChns_DiffClassifiers');
    fn = sprintf('PerformanceCumsumVsTime_DiffStartT_SSBeta_AllChns_DiffClassifiers');
end

if saveOp == 1
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
end

end