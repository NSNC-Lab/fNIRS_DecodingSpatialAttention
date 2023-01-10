% CV where SS beta coefficients from training fold is passed to test fold
%
% Only for sbj 08 and 10.
% Here, focus on multi-movies condition only. Faster

function preprocessFNIRS06_CV_GLM_ssBeta(s,numClasses,rejTrOp,rejChnOp)

sbjNum = s.name;
rawDataFN = s.fName;
movieList = s.movieList;
behData = s.resp;
startT = s.startT;
endT = s.endT;

saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
% Convert nirs to snirf file format
% snirf1 is entire data
snirf1 = SnirfClass(load([rawDataFN{1} '.nirs'],'-mat'));
snirf3 = SnirfClass(load([rawDataFN{2} '.nirs'],'-mat'));
snirf4 = SnirfClass(load([rawDataFN{3} '.nirs'],'-mat'));
snirf1.Info()

% Extract aux and convert to stimclass
allS2 = find(snirf1.aux(1,1).dataTimeSeries>1);
aInd2 = find(allS2(2:end)-allS2(1:end-1)==1);
allS3 = find(snirf3.aux(1,1).dataTimeSeries>1);
aInd3 = find(allS3(2:end)-allS3(1:end-1)==1);
allS4 = find(snirf4.aux(1,1).dataTimeSeries>1);
aInd4 = find(allS4(2:end)-allS4(1:end-1)==1);
allS2(aInd2) = [];
allS2 = allS2./50;
allS3(aInd3) = [];
allS3 = allS3./50;
allS4(aInd4) = [];
allS4 = allS4./50;

% oldStimClass1Length2 = size(allS2,1);
% oldStimClass1Length3 = oldStimClass1Length2 + size(allS3,1);
% oldStimClass1Length4 = oldStimClass1Length3 + size(allS4,1);

if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
    cueOnsetIndex = 1:4:720;
else
    cueOnsetIndex = 1:4:1080;
end
% cueOnsetIndex2 = 1:oldStimClass1Length2;
% cueOnsetIndex3 = oldStimClass1Length2+1:oldStimClass1Length3;
% cueOnsetIndex4 = oldStimClass1Length3+1:oldStimClass1Length4;

% trialNum2 = sum(ismember(cueOnsetIndex2,cueOnsetIndex));
% trialNum3 = sum(ismember(cueOnsetIndex3,cueOnsetIndex));
% trialNum4 = sum(ismember(cueOnsetIndex4,cueOnsetIndex));

% lastTrial2 = trialNum2;
% lastTrial3 = lastTrial2 + trialNum3;
% lastTrial4 = lastTrial3 + trialNum4;

snirf2TLen = size(snirf1.data.time,1)/(50);
snirf3TLen = size(snirf3.data.time,1)/(50);
allS3 = allS3+snirf2TLen;
allS4 = allS4+snirf2TLen+snirf3TLen;

allS = [allS2; allS3; allS4];

allS = allS(cueOnsetIndex);

if startT ~= 1
    % allS is in sec
    allS = allS-startT;
end

fs = 50;
%timePt = (0:0.25*fs:5*fs)+2*fs;
timePt = (0:0.25*fs:7*fs);

% split triggers into 6 categories
load([saveDir filesep movieList '.mat'],'indexMoviesTest');
load([saveDir filesep behData '.mat'],'responsesA','responsesV','correctRespA','correctRespV');

indexMoviesTest = updateMovieList(allS,indexMoviesTest);

dAll = DataClass([snirf1.data.dataTimeSeries; snirf3.data.dataTimeSeries; snirf4.data.dataTimeSeries],...
    0:1/50:(size(snirf1.data.dataTimeSeries,1)+size(snirf3.data.dataTimeSeries,1)+size(snirf4.data.dataTimeSeries,1)-1)/50,snirf1.data.measurementList);

if endT ~= -1
    dAll.dataTimeSeries = dAll.dataTimeSeries(startT*fs:endT*fs,:);
    dAll.time = dAll.time(startT*fs:endT*fs);
else
    dAll.dataTimeSeries = dAll.dataTimeSeries(startT*fs:end,:);
    dAll.time = dAll.time(startT*fs:end);
end

data = dAll;
probe = snirf1.probe;
mlActMan = {};
tIncMan = {};
Aaux = [];
rcMap = [];

if rejChnOp
    mlActAuto = hmrR_PruneChannels(data,probe,mlActMan,tIncMan,[7000  10000000],1.5,[0  45]);
else
    mlActAuto = {ones(size(snirf1.data.measurementList,2),1)};
end

dod = hmrR_Intensity2OD(data);

% For sbj 12, params in 1st line for motion correct/artifact perform
% significantly better than params in 2nd line.
[dod] = hmrR_MotionCorrectPCArecurse(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,15,4,0.97,5,1);
%[dod,tInc,svs,nSV,tInc0] = hmrR_MotionCorrectPCArecurse(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,20,5,0.97,5,1);

tIncAuto = hmrR_MotionArtifact(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,15,5);
%tIncAuto = hmrR_MotionArtifact(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,20,4);

if rejTrOp == 1
    
    leftOnsetMAll = StimClass('leftMulti');
    rightOnsetMAll = StimClass('rightMulti');
    centerOnsetMAll = StimClass('centerMulti');
    
    leftMultiIndexAll = indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1;
    rightMultiIndexAll = indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1;
    centerMultiIndexAll = indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==1;
    
    AddStims(leftOnsetMAll, allS(leftMultiIndexAll));
    AddStims(rightOnsetMAll, allS(rightMultiIndexAll));
    AddStims(centerOnsetMAll, allS(centerMultiIndexAll));

    updateStates(leftOnsetMAll);
    updateStates(rightOnsetMAll);
    updateStates(centerOnsetMAll);
   
    stimAll(1,1) = leftOnsetMAll;
    stimAll(1,2) = rightOnsetMAll;
    stimAll(1,3) = centerOnsetMAll;
    
    [stimAll,~] = hmrR_StimRejection(dod,stimAll,tIncAuto,tIncMan,[-2  15]);
   
    [indexMoviesTest] = updateMovieListStates(indexMoviesTest,rejTrOp,stimAll);
    
else
    
    [indexMoviesTest] = updateMovieListStates(indexMoviesTest,rejTrOp);
    
end

% 90 trials. Multi-movies condition only.
trialIdx = (responsesA==correctRespA)&(responsesV==correctRespV)&(indexMoviesTest(:,5)==1)'&(indexMoviesTest(:,7)>0)';
if numClasses == 2
    trialIdx = trialIdx & (indexMoviesTest(:,2)==1|indexMoviesTest(:,2)==2)';
end
numCorrectTrials = sum(trialIdx);
indexMoviesMultiple = indexMoviesTest(trialIdx,:);
allS = allS(trialIdx);
%cueOnsetIndexMultiple = cueOnsetIndex(trialIdx);

% Initialize variables
kFolds = 5;
nReps = 10;

performanceLDALedoitHbO = zeros(nReps*kFolds,length(timePt));
performanceLDACERNNHbO = zeros(nReps*kFolds,length(timePt));
performanceLogRegHbO = zeros(nReps*kFolds,length(timePt));
performanceSVMHbO = zeros(nReps*kFolds,length(timePt));
performanceBaggingHbO = zeros(nReps*kFolds,length(timePt));
performanceBoostedHbO = zeros(nReps*kFolds,length(timePt));

performanceLDALedoitHbR = zeros(nReps*kFolds,length(timePt));
performanceLDACERNNHbR = zeros(nReps*kFolds,length(timePt));
performanceLogRegHbR = zeros(nReps*kFolds,length(timePt));
performanceSVMHbR = zeros(nReps*kFolds,length(timePt));
performanceBaggingHbR = zeros(nReps*kFolds,length(timePt));
performanceBoostedHbR = zeros(nReps*kFolds,length(timePt));

performanceLDALedoitHbT = zeros(nReps*kFolds,length(timePt));
performanceLDACERNNHbT = zeros(nReps*kFolds,length(timePt));
performanceLogRegHbT = zeros(nReps*kFolds,length(timePt));
performanceSVMHbT = zeros(nReps*kFolds,length(timePt));
performanceBaggingHbT = zeros(nReps*kFolds,length(timePt));
performanceBoostedHbT = zeros(nReps*kFolds,length(timePt));

for iRep = 1:nReps
    
    cvp = cvpartition(numCorrectTrials, 'KFold',kFolds);
    
    for iFold = 1:kFolds

        foldIdx = (iRep-1)*kFolds+iFold;
        
        % snirf 2. Cell of StimClass array.
%         leftOnsetSTr = StimClass('leftSingle');
%         rightOnsetSTr = StimClass('rightSingle');
%         centerOnsetSTr = StimClass('centerSingle');
        leftOnsetMTr = StimClass('leftMulti');
        rightOnsetMTr = StimClass('rightMulti');
        centerOnsetMTr = StimClass('centerMulti');

%         leftOnsetSTst = StimClass('leftSingle');
%         rightOnsetSTst = StimClass('rightSingle');
%         centerOnsetSTst = StimClass('centerSingle');
        leftOnsetMTst = StimClass('leftMulti');
        rightOnsetMTst = StimClass('rightMulti');
        centerOnsetMTst = StimClass('centerMulti');

%         indexMoviesTrainSingle = indexMoviesSingle(cvp.training(iFold),:);
        indexMoviesTrainMultiple = indexMoviesMultiple(cvp.training(iFold),:);

%         cueOnsetIndexSingleTr = cueOnsetIndexSingle(cvp.training(iFold));
%         cueOnsetIndexMultipleTr = cueOnsetIndexMultiple(cvp.training(iFold));

%         indexMoviesTestSingle = indexMoviesSingle(cvp.test(iFold),:);
        indexMoviesTestMultiple = indexMoviesMultiple(cvp.test(iFold),:);
%         cueOnsetIndexSingleTst = cueOnsetIndexSingle(cvp.test(iFold));
%         cueOnsetIndexMultipleTst = cueOnsetIndexMultiple(cvp.test(iFold));

        allSMultipleTr = allS(cvp.training(iFold));
        allSMultipleTst = allS(cvp.test(iFold));

        % index of training size (72 trials), boolean for left condition
%         leftSingleIndexTr = indexMoviesTrainSingle(:,2)==2 & indexMoviesTrainSingle(:,5)==0;
%         rightSingleIndexTr = indexMoviesTrainSingle(:,2)==1 & indexMoviesTrainSingle(:,5)==0;
%         centerSingleIndexTr = indexMoviesTrainSingle(:,2)==3 & indexMoviesTrainSingle(:,5)==0;
        leftMultiIndexTr = indexMoviesTrainMultiple(:,2)==2 & indexMoviesTrainMultiple(:,5)==1;
        rightMultiIndexTr = indexMoviesTrainMultiple(:,2)==1 & indexMoviesTrainMultiple(:,5)==1;
        centerMultiIndexTr = indexMoviesTrainMultiple(:,2)==3 & indexMoviesTrainMultiple(:,5)==1;

%         leftSingleIndexTst = indexMoviesTestSingle(:,2)==2 & indexMoviesTestSingle(:,5)==0;
%         rightSingleIndexTst = indexMoviesTestSingle(:,2)==1 & indexMoviesTestSingle(:,5)==0;
%         centerSingleIndexTst = indexMoviesTestSingle(:,2)==3 & indexMoviesTestSingle(:,5)==0;
        leftMultiIndexTst = indexMoviesTestMultiple(:,2)==2 & indexMoviesTestMultiple(:,5)==1;
        rightMultiIndexTst = indexMoviesTestMultiple(:,2)==1 & indexMoviesTestMultiple(:,5)==1;
        centerMultiIndexTst = indexMoviesTestMultiple(:,2)==3 & indexMoviesTestMultiple(:,5)==1;

%         AddStims(leftOnsetSTr, allS(cueOnsetIndexSingleTr(leftSingleIndexTr)));
%         AddStims(rightOnsetSTr, allS(cueOnsetIndexSingleTr(rightSingleIndexTr)));
%         AddStims(centerOnsetSTr, allS(cueOnsetIndexSingleTr(centerSingleIndexTr)));
        AddStims(leftOnsetMTr, allSMultipleTr(leftMultiIndexTr));
        AddStims(rightOnsetMTr, allSMultipleTr(rightMultiIndexTr));
        AddStims(centerOnsetMTr, allSMultipleTr(centerMultiIndexTr));

%         AddStims(leftOnsetSTst, allS(cueOnsetIndexSingleTst(leftSingleIndexTst)));
%         AddStims(rightOnsetSTst, allS(cueOnsetIndexSingleTst(rightSingleIndexTst)));
%         AddStims(centerOnsetSTst, allS(cueOnsetIndexSingleTst(centerSingleIndexTst)));
        AddStims(leftOnsetMTst, allSMultipleTst(leftMultiIndexTst));
        AddStims(rightOnsetMTst, allSMultipleTst(rightMultiIndexTst));
        AddStims(centerOnsetMTst, allSMultipleTst(centerMultiIndexTst));

%         updateStates(leftOnsetSTr);
%         updateStates(rightOnsetSTr);
%         updateStates(centerOnsetSTr);
        updateStates(leftOnsetMTr);
        updateStates(rightOnsetMTr);
        updateStates(centerOnsetMTr);

%         updateStates(leftOnsetSTst);
%         updateStates(rightOnsetSTst);
%         updateStates(centerOnsetSTst);
        updateStates(leftOnsetMTst);
        updateStates(rightOnsetMTst);
        updateStates(centerOnsetMTst);

        % I prefer this over SetStim so I can control index
%         stimTr(1,1) = leftOnsetSTr;
%         stimTr(1,2) = rightOnsetSTr;
%         stimTr(1,3) = centerOnsetSTr;
        stimTr(1,1) = leftOnsetMTr;
        stimTr(1,2) = rightOnsetMTr;
        if numClasses == 3
            stimTr(1,3) = centerOnsetMTr;
        end

%         stimTst(1,1) = leftOnsetSTst;
%         stimTst(1,2) = rightOnsetSTst;
%         stimTst(1,3) = centerOnsetSTst;
        stimTst(1,1) = leftOnsetMTst;
        stimTst(1,2) = rightOnsetMTst;
        if numClasses == 3
            stimTst(1,3) = centerOnsetMTst;
        end

        %obj = SnirfClass(data, stim, probe, aux);
        %snirfTr = SnirfClass(dAll,stimTr,snirf1.probe,snirf1.aux);

        %stimTr = snirf1.stim;

        [stimTr,~] = hmrR_StimRejection(dod,stimTr,tIncAuto,tIncMan,[-2  15]);

        %dodBPFilt = hmrR_BandpassFilt(dod,0.01,0.5);
        dodBPFilt = hmrR_BandpassFilt(dod,0.01,10);
        %dod = hmrR_BandpassFilt(dod,0,0.5);

        dc = hmrR_OD2Conc(dodBPFilt,probe,[1  1  1]);

        % GLM here
%         [~,~,~,dcNewTr,~,~,~,~,~,~,beta] = hmrR_GLM_ssBeta(dc,stimTr,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
%             [-2  15],1,1,[1.0 1.0 0.0 0.0 0.0 0.0],15,1,0,0);
        
%         [~,~,~,dcNewTr,~,~,~,~,~,~,betaSS] = hmrR_GLM_HbT(dc,stimTr,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
%             [-2  15],1,1,[1.0 1.0 0.0 0.0 0.0 0.0],15,1,0,0);
        [~,~,~,dcNewTr,~,~,~,~,~,~,betaSS] = hmrR_GLM_MN(dc,stimTr,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
            [-2  15],1,1,[1.0 1.0 0.0 0.0 0.0 0.0],15,1,0,0);

        % format beta var here
        ssBeta = betaSS{1}(end,:,:);

        [stimTst,~] = hmrR_StimRejection(dod,stimTst,tIncAuto,tIncMan,[-2  15]);

        % technically, don't need this. Can just use dcNewTr
        dcNewTst = hmrR_ssBeta_CV(data,dc, probe, mlActAuto, tIncAuto, squeeze(ssBeta));

        % [dcAvg,dcAvgStd,nTrials,dcNew,dcResid,dcSum2,beta,R,hmrstats,bvar,ssBeta] = hmrR_ssBeta_CV(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
        %     [-2  15],1,1,[1.0 1.0 0.0 0.0 0.0 0.0],15,1,0,0);

        % extract and format training and test trials
%         if numClasses == 2
%             idxS = indexMoviesTrainSingle(:,2)==2 | indexMoviesTrainSingle(:,2)==1;
%             idxM = indexMoviesTrainMultiple(:,2)==2 | indexMoviesTrainMultiple(:,2)==1;
%             allSSingleTr = allS(cueOnsetIndexSingleTr(idxS));
%             allSMultipleTr = allS(cueOnsetIndexMultipleTr(idxM));
%             
%             idxS = indexMoviesTestSingle(:,2)==2 | indexMoviesTestSingle(:,2)==1;
%             idxM = indexMoviesTestMultiple(:,2)==2 | indexMoviesTestMultiple(:,2)==1;
%             allSSingleTst = allS(cueOnsetIndexSingleTst(idxS));
%             allSMultipleTst = allS(cueOnsetIndexMultipleTst(idxM));
%         else
%             allSSingleTr = allS(cueOnsetIndexSingleTr);
%             allSMultipleTr = allS(cueOnsetIndexMultipleTr);
%             allSSingleTst = allS(cueOnsetIndexSingleTst);
%             allSMultipleTst = allS(cueOnsetIndexMultipleTst);
%         end

%         if rejTrOp == 1
%         
%             [allSMultipleTr,indexMoviesTrainMultiple] = ...
%                 updateStimAndMovieList(stimTr,allSMultipleTr,indexMoviesTrainMultiple);
% 
%             [allSMultipleTst,indexMoviesTestMultiple] = ...
%                 updateStimAndMovieList(stimTst,allSMultipleTst,indexMoviesTestMultiple);
%         
%         end

        [trials_HbOM_Tr, trials_HbRM_Tr, trials_HbTM_Tr] ...
            = createSingleTrialHRF_MultiOnly_NoSave(dcNewTr, allSMultipleTr);
        
        [trials_HbOM_Tst, trials_HbRM_Tst, trials_HbTM_Tst] ...
            = createSingleTrialHRF_MultiOnly_NoSave(dcNewTst, allSMultipleTst);

        % train classifier
        [performanceLDALedoitHbO(foldIdx,:), performanceLDACERNNHbO(foldIdx,:),...
            performanceLogRegHbO(foldIdx,:), performanceSVMHbO(foldIdx,:),...
            performanceBaggingHbO(foldIdx,:), performanceBoostedHbO(foldIdx,:)] = ...
            calcCumSum_DiffStartT_AllChn_CompClassifiers_Final_NoSave(sbjNum,trials_HbOM_Tr,trials_HbOM_Tst,...
            indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlActAuto);
        
        [performanceLDALedoitHbR(foldIdx,:), performanceLDACERNNHbR(foldIdx,:),...
            performanceLogRegHbR(foldIdx,:), performanceSVMHbR(foldIdx,:),...
            performanceBaggingHbR(foldIdx,:), performanceBoostedHbR(foldIdx,:)] = ...
            calcCumSum_DiffStartT_AllChn_CompClassifiers_Final_NoSave(sbjNum,trials_HbRM_Tr,trials_HbRM_Tst,...
            indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlActAuto);
        
        [performanceLDALedoitHbT(foldIdx,:), performanceLDACERNNHbT(foldIdx,:),...
            performanceLogRegHbT(foldIdx,:), performanceSVMHbT(foldIdx,:),...
            performanceBaggingHbT(foldIdx,:), performanceBoostedHbT(foldIdx,:)] = ...
            calcCumSum_DiffStartT_AllChn_CompClassifiers_Final_NoSave(sbjNum,trials_HbTM_Tr,trials_HbTM_Tst,...
            indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlActAuto);

    end
     
end

% save folds to file

if ~exist(processedDataDir,'dir')
    mkdir(processedDataDir);
end

if numClasses == 2
    if rejTrOp
        %fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_LR_RejTr_SNR_1.5_OrigFunc.mat'];
        fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_LR_RejTr_SNR_1.5_OrigFunc_10Hz.mat'];
    else
        fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_LR_SNR_1.5.mat'];
    end
else
    if rejTrOp
        fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_RejTr_SNR_1.5_10Hz.mat'];
    else
        fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_SNR_1.5.mat'];
    end
end

save(fileName,'performanceLDALedoitHbO','performanceLDACERNNHbO',...
    'performanceLogRegHbO','performanceSVMHbO','performanceBaggingHbO',...
    'performanceBoostedHbO','performanceLDALedoitHbR','performanceLDACERNNHbR',...
    'performanceLogRegHbR','performanceSVMHbR','performanceBaggingHbR',...
    'performanceBoostedHbR','performanceLDALedoitHbT','performanceLDACERNNHbT',...
    'performanceLogRegHbT','performanceSVMHbT','performanceBaggingHbT',...
    'performanceBoostedHbT');
    
end