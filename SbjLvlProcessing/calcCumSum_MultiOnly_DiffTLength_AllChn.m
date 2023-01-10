% For sbj 08 and after.
%
% After running all sbjs, run grpPlotPerfVsDiffTLen.m
% 
% Add repts and CI. Used for publication
% Also test different "islands", subset of probe for portability
%
% Update to use trials from GLM CV SSBeta pipeline

function calcCumSum_MultiOnly_DiffTLength_AllChn(sbjNum,numClasses,saveOp,rejTrOp)
processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
rawDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
figSaveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' ...
    num2str(sbjNum) '\Figures\Classification'];

fs = 50;
timeLen = [0.1 0.2 0.5 1 1.5 2 3 4 5].*fs;
startT = 4*fs;
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

idxLFEF = 1:10;
idxRFEF = 11:20;
idxLSTG = 21;
idxRSTG = 22;
idxLIPS = [23 24 27 28];
idxRIPS = [25 26 29 30];

chnName = getChnName(numChn);

%figure(1);hold on;
figure('units','normalized','outerposition',[0 0 1 1]);hold on;

% Too lazy to modularize this
% Basis 1
load([processedDataDir filesep 'intermediateOutputsCombined_Basis1.mat']);
% Each array is channels x time x trial
if numClasses == 2
    if rejTrOp
        load([processedDataDir filesep 'singleTrialsUpdated_LR_Basis1.mat'],'singleTrialHRFHbOM',...
            'singleTrialHRFHbRM',...
            'singleTrialHRFHbTM','indexMoviesTest');
    else
        load([processedDataDir filesep 'singleTrialsUpdated_LR_Basis1.mat'],'singleTrialHRFHbOM',...
            'singleTrialHRFHbRM',...
            'singleTrialHRFHbTM','indexMoviesTest');
    end
    
    yLimAxis = [0 1];
else
    load([processedDataDir filesep 'singleTrialsUpdated_Basis1.mat'],'singleTrialHRFHbOM',...
        'singleTrialHRFHbRM',...
        'singleTrialHRFHbTM','indexMoviesTest');
    yLimAxis = [0 1];
end

if strcmp(sbjNum,'15')
    indexMoviesTest = indexMoviesTest(2:end,:);
end

% Offset % Chn x Time x Trials
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

%isSDNew = 1;

% after selectRS, cut down from 36 to 30 chns. Remove SS. Both steps can be
% implemented simultaneously using bit logic.
[singleTrialHRFHbOM,mlList] = selectRS(singleTrialHRFHbOM,1,mlList);
[singleTrialHRFHbRM,~] = selectRS(singleTrialHRFHbRM,1);
[singleTrialHRFHbTM,~] = selectRS(singleTrialHRFHbTM,1);

behFN = [rawDataDir filesep 'responses_' sbjNum];

if isSDNew
    [singleTrialHRFHbOM,movieIdx,behScore] = keepCorrectTrials(sbjNum,singleTrialHRFHbOM,behFN,numClasses,indexMoviesTest);
    [singleTrialHRFHbRM] = keepCorrectTrials(sbjNum,singleTrialHRFHbRM,behFN,numClasses,indexMoviesTest);
    [singleTrialHRFHbTM] = keepCorrectTrials(sbjNum,singleTrialHRFHbTM,behFN,numClasses,indexMoviesTest);
else
    [singleTrialHRFHbOM,movieIdx,behScore] = keepCorrectTrials_Special(sbjNum,singleTrialHRFHbOM,behFN,numClasses,indexMoviesTest);
    [singleTrialHRFHbRM] = keepCorrectTrials_Special(sbjNum,singleTrialHRFHbRM,behFN,numClasses,indexMoviesTest);
    [singleTrialHRFHbTM] = keepCorrectTrials_Special(sbjNum,singleTrialHRFHbTM,behFN,numClasses,indexMoviesTest);
end

performanceLDALedoitHbO = zeros(length(timeLen),nrept);
performanceLDALedoitHbR = zeros(length(timeLen),nrept);
performanceLDALedoitHbT = zeros(length(timeLen),nrept);

performance_LFEF_HbO = zeros(length(timeLen),nrept);
performance_LFEF_HbR = zeros(length(timeLen),nrept);
performance_LFEF_HbT = zeros(length(timeLen),nrept);

performance_RFEF_HbO = zeros(length(timeLen),nrept);
performance_RFEF_HbR = zeros(length(timeLen),nrept);
performance_RFEF_HbT = zeros(length(timeLen),nrept);

performance_LSTG_HbO = zeros(length(timeLen),nrept);
performance_LSTG_HbR = zeros(length(timeLen),nrept);
performance_LSTG_HbT = zeros(length(timeLen),nrept);

performance_RSTG_HbO = zeros(length(timeLen),nrept);
performance_RSTG_HbR = zeros(length(timeLen),nrept);
performance_RSTG_HbT = zeros(length(timeLen),nrept);

performance_LIPS_HbO = zeros(length(timeLen),nrept);
performance_LIPS_HbR = zeros(length(timeLen),nrept);
performance_LIPS_HbT = zeros(length(timeLen),nrept);

performance_RIPS_HbO = zeros(length(timeLen),nrept);
performance_RIPS_HbR = zeros(length(timeLen),nrept);
performance_RIPS_HbT = zeros(length(timeLen),nrept);

% performanceLDACERNNHbO = zeros(length(timeLen),nrept);
% performanceLDACERNNHbR = zeros(length(timeLen),nrept);
% performanceLDACERNNHbT = zeros(length(timeLen),nrept);
% 
% performanceLogRegHbO = zeros(length(timeLen),nrept);
% performanceLogRegHbR = zeros(length(timeLen),nrept);
% performanceLogRegHbT = zeros(length(timeLen),nrept);
% 
% performanceSVMHbO = zeros(length(timeLen),nrept);
% performanceSVMHbR = zeros(length(timeLen),nrept);
% performanceSVMHbT = zeros(length(timeLen),nrept);
% 
% performanceBaggingHbO = zeros(length(timeLen),nrept);
% performanceBaggingHbR = zeros(length(timeLen),nrept);
% performanceBaggingHbT = zeros(length(timeLen),nrept);
% 
% performanceBoostedHbO = zeros(length(timeLen),nrept);
% performanceBoostedHbR = zeros(length(timeLen),nrept);
% performanceBoostedHbT = zeros(length(timeLen),nrept);

performanceLDALedoitHbOFold = zeros(length(timeLen),nrept,kFold);
performanceLDALedoitHbRFold = zeros(length(timeLen),nrept,kFold);
performanceLDALedoitHbTFold = zeros(length(timeLen),nrept,kFold);

performance_LFEF_HbOFold = zeros(length(timeLen),nrept,kFold);
performance_LFEF_HbRFold = zeros(length(timeLen),nrept,kFold);
performance_LFEF_HbTFold = zeros(length(timeLen),nrept,kFold);

performance_RFEF_HbOFold = zeros(length(timeLen),nrept,kFold);
performance_RFEF_HbRFold = zeros(length(timeLen),nrept,kFold);
performance_RFEF_HbTFold = zeros(length(timeLen),nrept,kFold);

performance_LSTG_HbOFold = zeros(length(timeLen),nrept,kFold);
performance_LSTG_HbRFold = zeros(length(timeLen),nrept,kFold);
performance_LSTG_HbTFold = zeros(length(timeLen),nrept,kFold);

performance_RSTG_HbOFold = zeros(length(timeLen),nrept,kFold);
performance_RSTG_HbRFold = zeros(length(timeLen),nrept,kFold);
performance_RSTG_HbTFold = zeros(length(timeLen),nrept,kFold);

performance_LIPS_HbOFold = zeros(length(timeLen),nrept,kFold);
performance_LIPS_HbRFold = zeros(length(timeLen),nrept,kFold);
performance_LIPS_HbTFold = zeros(length(timeLen),nrept,kFold);

performance_RIPS_HbOFold = zeros(length(timeLen),nrept,kFold);
performance_RIPS_HbRFold = zeros(length(timeLen),nrept,kFold);
performance_RIPS_HbTFold = zeros(length(timeLen),nrept,kFold);

% performanceLDACERNNHbOFold = zeros(length(timeLen),nrept,kFold);
% performanceLDACERNNHbRFold = zeros(length(timeLen),nrept,kFold);
% performanceLDACERNNHbTFold = zeros(length(timeLen),nrept,kFold);
% 
% performanceLogRegHbOFold = zeros(length(timeLen),nrept,kFold);
% performanceLogRegHbRFold = zeros(length(timeLen),nrept,kFold);
% performanceLogRegHbTFold = zeros(length(timeLen),nrept,kFold);
% 
% performanceSVMHbOFold = zeros(length(timeLen),nrept,kFold);
% performanceSVMHbRFold = zeros(length(timeLen),nrept,kFold);
% performanceSVMHbTFold = zeros(length(timeLen),nrept,kFold);
% 
% performanceBaggingHbOFold = zeros(length(timeLen),nrept,kFold);
% performanceBaggingHbRFold = zeros(length(timeLen),nrept,kFold);
% performanceBaggingHbTFold = zeros(length(timeLen),nrept,kFold);
% 
% performanceBoostedHbOFold = zeros(length(timeLen),nrept,kFold);
% performanceBoostedHbRFold = zeros(length(timeLen),nrept,kFold);
% performanceBoostedHbTFold = zeros(length(timeLen),nrept,kFold);

performanceLDALedoitHbOCI = zeros(length(timeLen),2);
performanceLDALedoitHbRCI = zeros(length(timeLen),2);
performanceLDALedoitHbTCI = zeros(length(timeLen),2);

performance_LFEF_HbOCI = zeros(length(timeLen),2);
performance_LFEF_HbRCI = zeros(length(timeLen),2);
performance_LFEF_HbTCI = zeros(length(timeLen),2);

performance_RFEF_HbOCI = zeros(length(timeLen),2);
performance_RFEF_HbRCI = zeros(length(timeLen),2);
performance_RFEF_HbTCI = zeros(length(timeLen),2);

performance_LSTG_HbOCI = zeros(length(timeLen),2);
performance_LSTG_HbRCI = zeros(length(timeLen),2);
performance_LSTG_HbTCI = zeros(length(timeLen),2);

performance_RSTG_HbOCI = zeros(length(timeLen),2);
performance_RSTG_HbRCI = zeros(length(timeLen),2);
performance_RSTG_HbTCI = zeros(length(timeLen),2);

performance_LIPS_HbOCI = zeros(length(timeLen),2);
performance_LIPS_HbRCI = zeros(length(timeLen),2);
performance_LIPS_HbTCI = zeros(length(timeLen),2);

performance_RIPS_HbOCI = zeros(length(timeLen),2);
performance_RIPS_HbRCI = zeros(length(timeLen),2);
performance_RIPS_HbTCI = zeros(length(timeLen),2);

% performanceLDACERNNHbOCI = zeros(length(timeLen),2);
% performanceLDACERNNHbRCI = zeros(length(timeLen),2);
% performanceLDACERNNHbTCI = zeros(length(timeLen),2);
% 
% performanceLogRegHbOCI = zeros(length(timeLen),2);
% performanceLogRegHbRCI = zeros(length(timeLen),2);
% performanceLogRegHbTCI = zeros(length(timeLen),2);
% 
% performanceSVMHbOCI = zeros(length(timeLen),2);
% performanceSVMHbRCI = zeros(length(timeLen),2);
% performanceSVMHbTCI = zeros(length(timeLen),2);
% 
% performanceBaggingHbOCI = zeros(length(timeLen),2);
% performanceBaggingHbRCI = zeros(length(timeLen),2);
% performanceBaggingHbTCI = zeros(length(timeLen),2);
% 
% performanceBoostedHbOCI = zeros(length(timeLen),2);
% performanceBoostedHbRCI = zeros(length(timeLen),2);
% performanceBoostedHbTCI = zeros(length(timeLen),2);

% mean
performanceLDALedoitHbOMean = zeros(length(timeLen),1);
performanceLDALedoitHbRMean = zeros(length(timeLen),1);
performanceLDALedoitHbTMean = zeros(length(timeLen),1);

performance_LFEF_HbOMean = zeros(length(timeLen),1);
performance_LFEF_HbRMean = zeros(length(timeLen),1);
performance_LFEF_HbTMean = zeros(length(timeLen),1);

performance_RFEF_HbOMean = zeros(length(timeLen),1);
performance_RFEF_HbRMean = zeros(length(timeLen),1);
performance_RFEF_HbTMean = zeros(length(timeLen),1);

performance_LSTG_HbOMean = zeros(length(timeLen),1);
performance_LSTG_HbRMean = zeros(length(timeLen),1);
performance_LSTG_HbTMean = zeros(length(timeLen),1);

performance_RSTG_HbOMean = zeros(length(timeLen),1);
performance_RSTG_HbRMean = zeros(length(timeLen),1);
performance_RSTG_HbTMean = zeros(length(timeLen),1);

performance_LIPS_HbOMean = zeros(length(timeLen),1);
performance_LIPS_HbRMean = zeros(length(timeLen),1);
performance_LIPS_HbTMean = zeros(length(timeLen),1);

performance_RIPS_HbOMean = zeros(length(timeLen),1);
performance_RIPS_HbRMean = zeros(length(timeLen),1);
performance_RIPS_HbTMean = zeros(length(timeLen),1);

% performanceLDACERNNHbOMean = zeros(length(timeLen),1);
% performanceLDACERNNHbRMean = zeros(length(timeLen),1);
% performanceLDACERNNHbTMean = zeros(length(timeLen),1);
% 
% performanceLogRegHbOMean = zeros(length(timeLen),1);
% performanceLogRegHbRMean = zeros(length(timeLen),1);
% performanceLogRegHbTMean = zeros(length(timeLen),1);
% 
% performanceSVMHbOMean = zeros(length(timeLen),1);
% performanceSVMHbRMean = zeros(length(timeLen),1);
% performanceSVMHbTMean = zeros(length(timeLen),1);
% 
% performanceBaggingHbOMean = zeros(length(timeLen),1);
% performanceBaggingHbRMean = zeros(length(timeLen),1);
% performanceBaggingHbTMean = zeros(length(timeLen),1);
% 
% performanceBoostedHbOMean = zeros(length(timeLen),1);
% performanceBoostedHbRMean = zeros(length(timeLen),1);
% performanceBoostedHbTMean = zeros(length(timeLen),1);

for i2 = 2:length(timeLen)
    
    temp = cumsum(singleTrialHRFHbOM(:,startT:startT+timeLen(i2),:),2);
    cumsumHbOM = squeeze(temp(:,end,:));
    cumsumHbOM_LFEF = cumsumHbOM(idxLFEF,:,:);
    cumsumHbOM_RFEF = cumsumHbOM(idxRFEF,:,:);
    cumsumHbOM_LSTG = cumsumHbOM(idxLSTG,:,:);
    cumsumHbOM_RSTG = cumsumHbOM(idxRSTG,:,:);
    cumsumHbOM_LIPS = cumsumHbOM(idxLIPS,:,:);
    cumsumHbOM_RIPS = cumsumHbOM(idxRIPS,:,:);
    
%     temp = cumsum(singleTrialHRFHbRS(:,timePt(i):timePt(i)+lenT,:),2);
%     cumsumHbRS = squeeze(temp(:,end,:));
    
    temp = cumsum(singleTrialHRFHbRM(:,startT:startT+timeLen(i2),:),2);
    cumsumHbRM = squeeze(temp(:,end,:));
    cumsumHbRM_LFEF = cumsumHbRM(idxLFEF,:,:);
    cumsumHbRM_RFEF = cumsumHbRM(idxRFEF,:,:);
    cumsumHbRM_LSTG = cumsumHbRM(idxLSTG,:,:);
    cumsumHbRM_RSTG = cumsumHbRM(idxRSTG,:,:);
    cumsumHbRM_LIPS = cumsumHbRM(idxLIPS,:,:);
    cumsumHbRM_RIPS = cumsumHbRM(idxRIPS,:,:);
    
%     temp = cumsum(singleTrialHRFHbTS(:,startT:timePt(i),:),2);
%     cumsumHbTS = squeeze(temp(:,end,:));
    
    temp = cumsum(singleTrialHRFHbTM(:,startT:startT+timeLen(i2),:),2);
    cumsumHbTM = squeeze(temp(:,end,:));
    cumsumHbTM_LFEF = cumsumHbTM(idxLFEF,:,:);
    cumsumHbTM_RFEF = cumsumHbTM(idxRFEF,:,:);
    cumsumHbTM_LSTG = cumsumHbTM(idxLSTG,:,:);
    cumsumHbTM_RSTG = cumsumHbTM(idxRSTG,:,:);
    cumsumHbTM_LIPS = cumsumHbTM(idxLIPS,:,:);
    cumsumHbTM_RIPS = cumsumHbTM(idxRIPS,:,:);

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
        
        [performance_LFEF_HbO(i2,i3), performance_LFEF_HbOFold(i2,i3,:)] = train_RLDA_Ledoit(cumsumHbOM_LFEF,mlList,size(cumsumHbOM_LFEF,1),numClasses,movieIdx);
        [performance_LFEF_HbR(i2,i3), performance_LFEF_HbRFold(i2,i3,:)] = train_RLDA_Ledoit(cumsumHbRM_LFEF,mlList,size(cumsumHbRM_LFEF,1),numClasses,movieIdx);
        [performance_LFEF_HbT(i2,i3), performance_LFEF_HbTFold(i2,i3,:)] = train_RLDA_Ledoit(cumsumHbTM_LFEF,mlList,size(cumsumHbTM_LFEF,1),numClasses,movieIdx);
        
        [performance_RFEF_HbO(i2,i3), performance_RFEF_HbOFold(i2,i3,:)] = train_RLDA_Ledoit(cumsumHbOM_RFEF,mlList,size(cumsumHbOM_RFEF,1),numClasses,movieIdx);
        [performance_RFEF_HbR(i2,i3), performance_RFEF_HbRFold(i2,i3,:)] = train_RLDA_Ledoit(cumsumHbRM_RFEF,mlList,size(cumsumHbRM_RFEF,1),numClasses,movieIdx);
        [performance_RFEF_HbT(i2,i3), performance_RFEF_HbTFold(i2,i3,:)] = train_RLDA_Ledoit(cumsumHbTM_RFEF,mlList,size(cumsumHbTM_RFEF,1),numClasses,movieIdx);
        
        [performance_LSTG_HbO(i2,i3), performance_LSTG_HbOFold(i2,i3,:)] = train_RLDA_Ledoit(cumsumHbOM_LSTG,mlList,size(cumsumHbOM_LSTG,1),numClasses,movieIdx);
        [performance_LSTG_HbR(i2,i3), performance_LSTG_HbRFold(i2,i3,:)] = train_RLDA_Ledoit(cumsumHbRM_LSTG,mlList,size(cumsumHbRM_LSTG,1),numClasses,movieIdx);
        [performance_LSTG_HbT(i2,i3), performance_LSTG_HbTFold(i2,i3,:)] = train_RLDA_Ledoit(cumsumHbTM_LSTG,mlList,size(cumsumHbTM_LSTG,1),numClasses,movieIdx);
        
        [performance_RSTG_HbO(i2,i3), performance_RSTG_HbOFold(i2,i3,:)] = train_RLDA_Ledoit(cumsumHbOM_RSTG,mlList,size(cumsumHbOM_RSTG,1),numClasses,movieIdx);
        [performance_RSTG_HbR(i2,i3), performance_RSTG_HbRFold(i2,i3,:)] = train_RLDA_Ledoit(cumsumHbRM_RSTG,mlList,size(cumsumHbRM_RSTG,1),numClasses,movieIdx);
        [performance_RSTG_HbT(i2,i3), performance_RSTG_HbTFold(i2,i3,:)] = train_RLDA_Ledoit(cumsumHbTM_RSTG,mlList,size(cumsumHbTM_RSTG,1),numClasses,movieIdx);
        
        [performance_LIPS_HbO(i2,i3), performance_LIPS_HbOFold(i2,i3,:)] = train_RLDA_Ledoit(cumsumHbOM_LIPS,mlList,size(cumsumHbOM_LIPS,1),numClasses,movieIdx);
        [performance_LIPS_HbR(i2,i3), performance_LIPS_HbRFold(i2,i3,:)] = train_RLDA_Ledoit(cumsumHbRM_LIPS,mlList,size(cumsumHbRM_LIPS,1),numClasses,movieIdx);
        [performance_LIPS_HbT(i2,i3), performance_LIPS_HbTFold(i2,i3,:)] = train_RLDA_Ledoit(cumsumHbTM_LIPS,mlList,size(cumsumHbTM_LIPS,1),numClasses,movieIdx);
        
        [performance_RIPS_HbO(i2,i3), performance_RIPS_HbOFold(i2,i3,:)] = train_RLDA_Ledoit(cumsumHbOM_RIPS,mlList,size(cumsumHbOM_RIPS,1),numClasses,movieIdx);
        [performance_RIPS_HbR(i2,i3), performance_RIPS_HbRFold(i2,i3,:)] = train_RLDA_Ledoit(cumsumHbRM_RIPS,mlList,size(cumsumHbRM_RIPS,1),numClasses,movieIdx);
        [performance_RIPS_HbT(i2,i3), performance_RIPS_HbTFold(i2,i3,:)] = train_RLDA_Ledoit(cumsumHbTM_RIPS,mlList,size(cumsumHbTM_RIPS,1),numClasses,movieIdx);

%         [performanceLDACERNNHbO(i2,i3),performanceLDACERNNHbOFold(i2,i3,:)] = train_RLDA_CERNN(cumsumHbOM,mlList,numChn,numClasses,movieIdx);
%         [performanceLDACERNNHbR(i2,i3),performanceLDACERNNHbRFold(i2,i3,:)] = train_RLDA_CERNN(cumsumHbRM,mlList,numChn,numClasses,movieIdx);
%         [performanceLDACERNNHbT(i2,i3),performanceLDACERNNHbTFold(i2,i3,:)] = train_RLDA_CERNN(cumsumHbTM,mlList,numChn,numClasses,movieIdx);
% 
%     %     performanceQDALedoitHbO(i) = train_RQDA_Ledoit(cumsumHbOM,mlList,numChn,numClasses,movieIdx);
%     %     performanceQDALedoitHbR(i) = train_RQDA_Ledoit(cumsumHbRM,mlList,numChn,numClasses,movieIdx);
%     %     performanceQDALedoitHbT(i) = train_RQDA_Ledoit(cumsumHbTM,mlList,numChn,numClasses,movieIdx);
% 
%     %     performanceLDAGDHbO(i) = train_RLDA_Matlab(cumsumHbOM,mlList,numChn,numClasses,movieIdx);
%     %     performanceLDAGDHbR(i) = train_RLDA_Matlab(cumsumHbRM,mlList,numChn,numClasses,movieIdx);
%     %     performanceLDAGDHbT(i) = train_RLDA_Matlab(cumsumHbTM,mlList,numChn,numClasses,movieIdx);
% 
%         [performanceLogRegHbO(i2,i3),performanceLogRegHbOFold(i2,i3,:)] = trainClassifierLogisticRegressionHRF(cumsumHbOM,mlList,numChn,numClasses,movieIdx);
%         [performanceLogRegHbR(i2,i3),performanceLogRegHbRFold(i2,i3,:)] = trainClassifierLogisticRegressionHRF(cumsumHbRM,mlList,numChn,numClasses,movieIdx);
%         [performanceLogRegHbT(i2,i3),performanceLogRegHbTFold(i2,i3,:)] = trainClassifierLogisticRegressionHRF(cumsumHbTM,mlList,numChn,numClasses,movieIdx);
% 
% 
%     %     [performanceArrMultiHbO(i),predRespMHbO] = trainClassifierLinearDiscriminantSingleHRFAllChns(cumsumHbOM,mlActAuto{1},numClasses,movieIdx);
%     %     [performanceArrMultiHbR(i),predRespMHbR] = trainClassifierLinearDiscriminantSingleHRFAllChns(cumsumHbRM,mlActAuto{1},numClasses,movieIdx);
%     %     [performanceArrMultiHbT(i),predRespMHbT] = trainClassifierLinearDiscriminantSingleHRFAllChns(cumsumHbTM,mlActAuto{1},numClasses,movieIdx);
% 
%         [performanceSVMHbO(i2,i3),performanceSVMHbOFold(i2,i3,:)] = trainClassifierSVMHRF(cumsumHbOM,mlList,numChn,numClasses,movieIdx);
%         [performanceSVMHbR(i2,i3),performanceSVMHbRFold(i2,i3,:)] = trainClassifierSVMHRF(cumsumHbRM,mlList,numChn,numClasses,movieIdx);
%         [performanceSVMHbT(i2,i3),performanceSVMHbTFold(i2,i3,:)] = trainClassifierSVMHRF(cumsumHbTM,mlList,numChn,numClasses,movieIdx);
% 
%         [performanceBaggingHbO(i2,i3),performanceBaggingHbOFold(i2,i3,:)] = trainClassifierBaggingTreesSingleHRF(cumsumHbOM,mlList,numChn,numClasses,movieIdx);
%         [performanceBaggingHbR(i2,i3),performanceBaggingHbRFold(i2,i3,:)] = trainClassifierBaggingTreesSingleHRF(cumsumHbRM,mlList,numChn,numClasses,movieIdx);
%         [performanceBaggingHbT(i2,i3),performanceBaggingHbTFold(i2,i3,:)] = trainClassifierBaggingTreesSingleHRF(cumsumHbTM,mlList,numChn,numClasses,movieIdx);
% 
%         [performanceBoostedHbO(i2,i3),performanceBoostedHbOFold(i2,i3,:)] = trainClassifierBoostedTreesSingleHRF(cumsumHbOM,mlList,numChn,numClasses,movieIdx);
%         [performanceBoostedHbR(i2,i3),performanceBoostedHbRFold(i2,i3,:)] = trainClassifierBoostedTreesSingleHRF(cumsumHbRM,mlList,numChn,numClasses,movieIdx);
%         [performanceBoostedHbT(i2,i3),performanceBoostedHbTFold(i2,i3,:)] = trainClassifierBoostedTreesSingleHRF(cumsumHbTM,mlList,numChn,numClasses,movieIdx);
    end
end

if numClasses == 2
    savePerfFN = 'PerformanceCumsumVsTime_DiffTLen_StartMovie_LR_AllChns.mat';
else
    savePerfFN = 'PerformanceCumsumVsTime_DiffTLen_StartMovie_AllChns.mat';
end
save([processedDataDir filesep savePerfFN]);

conf = 0.95;

for i = 1:length(timeLen)
    
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
    
    % LFEF
    temp = performance_LFEF_HbOFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performance_LFEF_HbOMean(i) = mean(temp(:));
    performance_LFEF_HbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performance_LFEF_HbRFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performance_LFEF_HbRMean(i) = mean(temp(:));
    performance_LFEF_HbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performance_LFEF_HbTFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performance_LFEF_HbTMean(i) = mean(temp(:));
    performance_LFEF_HbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % RFEF
    temp = performance_RFEF_HbOFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performance_RFEF_HbOMean(i) = mean(temp(:));
    performance_RFEF_HbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performance_RFEF_HbRFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performance_RFEF_HbRMean(i) = mean(temp(:));
    performance_RFEF_HbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performance_RFEF_HbTFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performance_RFEF_HbTMean(i) = mean(temp(:));
    performance_RFEF_HbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    %LSTG
    temp = performance_LSTG_HbOFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performance_LSTG_HbOMean(i) = mean(temp(:));
    performance_LSTG_HbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performance_LSTG_HbRFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performance_LSTG_HbRMean(i) = mean(temp(:));
    performance_LSTG_HbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performance_LSTG_HbTFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performance_LSTG_HbTMean(i) = mean(temp(:));
    performance_LSTG_HbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % RSTG
    temp = performance_RSTG_HbOFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performance_RSTG_HbOMean(i) = mean(temp(:));
    performance_RSTG_HbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performance_RSTG_HbRFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performance_RSTG_HbRMean(i) = mean(temp(:));
    performance_RSTG_HbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performance_RSTG_HbTFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performance_RSTG_HbTMean(i) = mean(temp(:));
    performance_RSTG_HbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    % LIPS
    temp = performance_LIPS_HbOFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performance_LIPS_HbOMean(i) = mean(temp(:));
    performance_LIPS_HbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performance_LIPS_HbRFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performance_LIPS_HbRMean(i) = mean(temp(:));
    performance_LIPS_HbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performance_LIPS_HbTFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performance_LIPS_HbTMean(i) = mean(temp(:));
    performance_LIPS_HbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    %RIPS
    temp = performance_RIPS_HbOFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performance_RIPS_HbOMean(i) = mean(temp(:));
    performance_RIPS_HbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performance_RIPS_HbRFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performance_RIPS_HbRMean(i) = mean(temp(:));
    performance_RIPS_HbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
    temp = performance_RIPS_HbTFold(i,:,:);
    se = std(temp(:))/sqrt(nrept*kFold);
    performance_RIPS_HbTMean(i) = mean(temp(:));
    performance_RIPS_HbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
%     % CERNN
%     temp = performanceLDACERNNHbOFold(i,:,:);
%     se = std(temp(:))/sqrt(nrept*kFold);
%     performanceLDACERNNHbOMean(i) = mean(temp(:));
%     performanceLDACERNNHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
%     
%     temp = performanceLDACERNNHbRFold(i,:,:);
%     se = std(temp(:))/sqrt(nrept*kFold);
%     performanceLDACERNNHbRMean(i) = mean(temp(:));
%     performanceLDACERNNHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
%     
%     temp = performanceLDACERNNHbTFold(i,:,:);
%     se = std(temp(:))/sqrt(nrept*kFold);
%     performanceLDACERNNHbTMean(i) = mean(temp(:));
%     performanceLDACERNNHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
%     
%     % Log
%     temp = performanceLogRegHbOFold(i,:,:);
%     se = std(temp(:))/sqrt(nrept*kFold);
%     performanceLogRegHbOMean(i) = mean(temp(:));
%     performanceLogRegHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
%     
%     temp = performanceLogRegHbRFold(i,:,:);
%     se = std(temp(:))/sqrt(nrept*kFold);
%     performanceLogRegHbRMean(i) = mean(temp(:));
%     performanceLogRegHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
%     
%     temp = performanceLogRegHbTFold(i,:,:);
%     se = std(temp(:))/sqrt(nrept*kFold);
%     performanceLogRegHbTMean(i) = mean(temp(:));
%     performanceLogRegHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
%     
%     % SVM
%     temp = performanceSVMHbOFold(i,:,:);
%     se = std(temp(:))/sqrt(nrept*kFold);
%     performanceSVMHbOMean(i) = mean(temp(:));
%     performanceSVMHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
%     
%     temp = performanceSVMHbRFold(i,:,:);
%     se = std(temp(:))/sqrt(nrept*kFold);
%     performanceSVMHbRMean(i) = mean(temp(:));
%     performanceSVMHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
%     
%     temp = performanceSVMHbTFold(i,:,:);
%     se = std(temp(:))/sqrt(nrept*kFold);
%     performanceSVMHbTMean(i) = mean(temp(:));
%     performanceSVMHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
%     
%     % Bagging
%     temp = performanceBaggingHbOFold(i,:,:);
%     se = std(temp(:))/sqrt(nrept*kFold);
%     performanceBaggingHbOMean(i) = mean(temp(:));
%     performanceBaggingHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
%     
%     temp = performanceBaggingHbRFold(i,:,:);
%     se = std(temp(:))/sqrt(nrept*kFold);
%     performanceBaggingHbRMean(i) = mean(temp(:));
%     performanceBaggingHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
%     
%     temp = performanceBaggingHbTFold(i,:,:);
%     se = std(temp(:))/sqrt(nrept*kFold);
%     performanceBaggingHbTMean(i) = mean(temp(:));
%     performanceBaggingHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
%     
%     % Boosting
%     temp = performanceBoostedHbOFold(i,:,:);
%     se = std(temp(:))/sqrt(nrept*kFold);
%     performanceBoostedHbOMean(i) = mean(temp(:));
%     performanceBoostedHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
%     
%     temp = performanceBoostedHbRFold(i,:,:);
%     se = std(temp(:))/sqrt(nrept*kFold);
%     performanceBoostedHbRMean(i) = mean(temp(:));
%     performanceBoostedHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
%     
%     temp = performanceBoostedHbTFold(i,:,:);
%     se = std(temp(:))/sqrt(nrept*kFold);
%     performanceBoostedHbTMean(i) = mean(temp(:));
%     performanceBoostedHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    
end

numSubSets = 7;
cmap = jet(numSubSets);

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

errorbar(timeLen./fs,performanceLDALedoitHbOMean,performanceLDALedoitHbOCI(:,2),'Color',cmap(1,:));
errorbar(timeLen./fs,performance_LFEF_HbOMean,performance_LFEF_HbOCI(:,2),'Color',cmap(2,:));
errorbar(timeLen./fs,performance_RFEF_HbOMean,performance_RFEF_HbOCI(:,2),'Color',cmap(3,:));
errorbar(timeLen./fs,performance_LSTG_HbOMean,performance_LSTG_HbOCI(:,2),'Color',cmap(4,:));
errorbar(timeLen./fs,performance_RSTG_HbOMean,performance_RSTG_HbOCI(:,2),'Color',cmap(5,:));
errorbar(timeLen./fs,performance_LIPS_HbOMean,performance_LIPS_HbOCI(:,2),'Color',cmap(6,:));
errorbar(timeLen./fs,performance_RIPS_HbOMean,performance_RIPS_HbOCI(:,2),'Color',cmap(7,:));
% errorbar(timeLen./fs-2,performanceLDACERNNHbOMean,performanceLDACERNNHbOCI(:,2),'Color',cmap(2,:));
% % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
% % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
% errorbar(timeLen./fs-2,performanceLogRegHbOMean,performanceLogRegHbOCI(:,2),'Color',cmap(3,:));
% errorbar(timeLen./fs-2,performanceSVMHbOMean,performanceSVMHbOCI(:,2),'Color',cmap(4,:));
% errorbar(timeLen./fs-2,performanceBaggingHbOMean,performanceBaggingHbOCI(:,2),'Color',cmap(5,:));
% errorbar(timeLen./fs-2,performanceBoostedHbOMean,performanceBoostedHbOCI(:,2),'Color',cmap(6,:));

ylim(yLimAxis);
%xlim([0 7]);
title(sprintf('Sbj %s: Δ[HbO] Multi',num2str(sbjNum)));
% legend({'PCA 95%','PCA 85%','PCA 75%','PCA 65%','PCA 55%','LDAShrink',...
%     'QDAShrink','LDAGD','LogReg','SVM'});
% legend({'PCA 95%','PCA 85%','PCA 75%','PCA 65%','PCA 55%','LDAShrink',...
%     'LDA CERNN','LogReg','SVM'});
%legend({'LDAShrink','LDA CERNN','LogReg','SVM','Bagging','Boosted'});
%legend('LDAShrink');
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

errorbar(timeLen./fs,performanceLDALedoitHbRMean,performanceLDALedoitHbRCI(:,2),'Color',cmap(1,:));
errorbar(timeLen./fs,performance_LFEF_HbRMean,performance_LFEF_HbRCI(:,2),'Color',cmap(2,:));
errorbar(timeLen./fs,performance_RFEF_HbRMean,performance_RFEF_HbRCI(:,2),'Color',cmap(3,:));
errorbar(timeLen./fs,performance_LSTG_HbRMean,performance_LSTG_HbRCI(:,2),'Color',cmap(4,:));
errorbar(timeLen./fs,performance_RSTG_HbRMean,performance_RSTG_HbRCI(:,2),'Color',cmap(5,:));
errorbar(timeLen./fs,performance_LIPS_HbRMean,performance_LIPS_HbRCI(:,2),'Color',cmap(6,:));
errorbar(timeLen./fs,performance_RIPS_HbRMean,performance_RIPS_HbRCI(:,2),'Color',cmap(7,:));
% errorbar(timeLen./fs-2,performanceLDACERNNHbRMean,performanceLDACERNNHbRCI(:,2),'Color',cmap(2,:));
% % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
% % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
% errorbar(timeLen./fs-2,performanceLogRegHbRMean,performanceLogRegHbRCI(:,2),'Color',cmap(3,:));
% errorbar(timeLen./fs-2,performanceSVMHbRMean,performanceSVMHbRCI(:,2),'Color',cmap(4,:));
% errorbar(timeLen./fs-2,performanceBaggingHbRMean,performanceBaggingHbRCI(:,2),'Color',cmap(5,:));
% errorbar(timeLen./fs-2,performanceBoostedHbRMean,performanceBoostedHbRCI(:,2),'Color',cmap(6,:));

ylim(yLimAxis);
%xlim([0 7]);
title(sprintf('Sbj %s: Δ[HbR] Multi',num2str(sbjNum)));
%legend('S1D1','S1D2','S1D3','S1D4');
legend({'All','LFEF','RFEF','LSTG','RSTG','LIPS','RIPS'});
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

errorbar(timeLen./fs,performanceLDALedoitHbTMean,performanceLDALedoitHbTCI(:,2),'Color',cmap(1,:));
errorbar(timeLen./fs,performance_LFEF_HbTMean,performance_LFEF_HbTCI(:,2),'Color',cmap(2,:));
errorbar(timeLen./fs,performance_RFEF_HbTMean,performance_RFEF_HbTCI(:,2),'Color',cmap(3,:));
errorbar(timeLen./fs,performance_LSTG_HbTMean,performance_LSTG_HbTCI(:,2),'Color',cmap(4,:));
errorbar(timeLen./fs,performance_RSTG_HbTMean,performance_RSTG_HbTCI(:,2),'Color',cmap(5,:));
errorbar(timeLen./fs,performance_LIPS_HbTMean,performance_LIPS_HbTCI(:,2),'Color',cmap(6,:));
errorbar(timeLen./fs,performance_RIPS_HbTMean,performance_RIPS_HbTCI(:,2),'Color',cmap(7,:));
% errorbar(timeLen./fs-2,performanceLDACERNNHbTMean,performanceLDACERNNHbTCI(:,2),'Color',cmap(2,:));
% % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
% % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
% errorbar(timeLen./fs-2,performanceLogRegHbTMean,performanceLogRegHbTCI(:,2),'Color',cmap(3,:));
% errorbar(timeLen./fs-2,performanceSVMHbTMean,performanceSVMHbTCI(:,2),'Color',cmap(4,:));
% errorbar(timeLen./fs-2,performanceBaggingHbTMean,performanceBaggingHbTCI(:,2),'Color',cmap(5,:));
% errorbar(timeLen./fs-2,performanceBoostedHbTMean,performanceBoostedHbTCI(:,2),'Color',cmap(6,:));

ylim(yLimAxis);
%xlim([0 7]);
title(sprintf('Sbj %s: Δ[HbT] Multi',num2str(sbjNum)));
%legend('S1D1','S1D2','S1D3','S1D4');
%legend('LDAShrink');
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
    fn = sprintf('PerformanceCumsumVsTime_DiffTLen_StartMovie_LR_AllChns');
else
    fn = sprintf('PerformanceCumsumVsTime_DiffTLen_StartMovie_AllChns');
end

if saveOp == 1
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
end

end