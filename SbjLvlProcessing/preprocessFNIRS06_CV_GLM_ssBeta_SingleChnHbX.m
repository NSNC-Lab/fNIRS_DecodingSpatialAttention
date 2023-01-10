% CV where SS beta coefficients from training fold is passed to test fold
% 03/24/2022: Remove rejected trials from classification
%   Remove cueOnsetIndex, ineffective method. Just directly filter allS
% Only for sbj 12 and after

function preprocessFNIRS06_CV_GLM_ssBeta_SingleChnHbX(s,numClasses,rejTrOp,rejChnOp)

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

if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
    cueOnsetIndex = 1:4:720;
else
    cueOnsetIndex = 1:4:360;
end

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

load([saveDir filesep movieList '.mat'],'indexMoviesTest');
load([saveDir filesep behData '.mat'],'responsesA','responsesV','correctRespA','correctRespV');

indexMoviesTest = updateMovieList(allS,indexMoviesTest);

% dAll = DataClass([snirf1.data.dataTimeSeries; snirf3.data.dataTimeSeries; snirf4.data.dataTimeSeries],...
%     0:1/50:(size(snirf1.data.dataTimeSeries,1)+size(snirf3.data.dataTimeSeries,1)+size(snirf4.data.dataTimeSeries,1)-1)/50,snirf1.data.measurementList);

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

% Here define stim class solely for trials rejection.
%[stimTr,~] = hmrR_StimRejection(dod,stim,tIncAuto,tIncMan,[-2  15]);
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

% 90 trials
trialIdx = (responsesA==correctRespA)&(responsesV==correctRespV)&(indexMoviesTest(:,5)==1)'&(indexMoviesTest(:,7)>0)';
if numClasses == 2
    trialIdx = trialIdx & (indexMoviesTest(:,2)==1|indexMoviesTest(:,2)==2)';
end
if strcmp(sbjNum,'15')
    trialIdx(1) = 0;
end
numCorrectTrials = sum(trialIdx);
indexMoviesMultiple = indexMoviesTest(trialIdx,:);
allS = allS(trialIdx);
%cueOnsetIndexMultiple = cueOnsetIndex(trialIdx);

% Initialize variables
kFolds = 5;
nReps = 10;

% performanceLDALedoitHbO = zeros(nReps*kFolds,length(timePt));
% performanceLDACERNNHbO = zeros(nReps*kFolds,length(timePt));
% performanceLogRegHbO = zeros(nReps*kFolds,length(timePt));
% performanceSVMHbO = zeros(nReps*kFolds,length(timePt));
% performanceBaggingHbO = zeros(nReps*kFolds,length(timePt));
% performanceBoostedHbO = zeros(nReps*kFolds,length(timePt));
% 
% performanceLDALedoitHbR = zeros(nReps*kFolds,length(timePt));
% performanceLDACERNNHbR = zeros(nReps*kFolds,length(timePt));
% performanceLogRegHbR = zeros(nReps*kFolds,length(timePt));
% performanceSVMHbR = zeros(nReps*kFolds,length(timePt));
% performanceBaggingHbR = zeros(nReps*kFolds,length(timePt));
% performanceBoostedHbR = zeros(nReps*kFolds,length(timePt));
% 
% performanceLDALedoitHbT = zeros(nReps*kFolds,length(timePt));
% performanceLDACERNNHbT = zeros(nReps*kFolds,length(timePt));
% performanceLogRegHbT = zeros(nReps*kFolds,length(timePt));
% performanceSVMHbT = zeros(nReps*kFolds,length(timePt));
% performanceBaggingHbT = zeros(nReps*kFolds,length(timePt));
% performanceBoostedHbT = zeros(nReps*kFolds,length(timePt));

performanceLDAHbO = zeros(nReps*kFolds,30,length(timePt));
performanceLDAHbR = zeros(nReps*kFolds,30,length(timePt));
performanceLDAHbT = zeros(nReps*kFolds,30,length(timePt));

for iRep = 1:nReps
    
    cvp = cvpartition(numCorrectTrials, 'KFold',kFolds);
    
    for iFold = 1:kFolds

        foldIdx = (iRep-1)*kFolds+iFold;
        
        % snirf 2. Cell of StimClass array.
        leftOnsetMTr = StimClass('leftMulti');
        rightOnsetMTr = StimClass('rightMulti');
        centerOnsetMTr = StimClass('centerMulti');

        leftOnsetMTst = StimClass('leftMulti');
        rightOnsetMTst = StimClass('rightMulti');
        centerOnsetMTst = StimClass('centerMulti');
        
        indexMoviesTrainMultiple = indexMoviesMultiple(cvp.training(iFold),:);

        indexMoviesTestMultiple = indexMoviesMultiple(cvp.test(iFold),:);
        
        allSMultipleTr = allS(cvp.training(iFold));
        allSMultipleTst = allS(cvp.test(iFold));

        leftMultiIndexTr = indexMoviesTrainMultiple(:,2)==2 & indexMoviesTrainMultiple(:,5)==1;
        rightMultiIndexTr = indexMoviesTrainMultiple(:,2)==1 & indexMoviesTrainMultiple(:,5)==1;
        centerMultiIndexTr = indexMoviesTrainMultiple(:,2)==3 & indexMoviesTrainMultiple(:,5)==1;

        leftMultiIndexTst = indexMoviesTestMultiple(:,2)==2 & indexMoviesTestMultiple(:,5)==1;
        rightMultiIndexTst = indexMoviesTestMultiple(:,2)==1 & indexMoviesTestMultiple(:,5)==1;
        centerMultiIndexTst = indexMoviesTestMultiple(:,2)==3 & indexMoviesTestMultiple(:,5)==1;

        % Only correct trials
        AddStims(leftOnsetMTr, allSMultipleTr(leftMultiIndexTr));
        AddStims(rightOnsetMTr, allSMultipleTr(rightMultiIndexTr));
        AddStims(centerOnsetMTr, allSMultipleTr(centerMultiIndexTr));

        AddStims(leftOnsetMTst, allSMultipleTst(leftMultiIndexTst));
        AddStims(rightOnsetMTst, allSMultipleTst(rightMultiIndexTst));
        AddStims(centerOnsetMTst, allSMultipleTst(centerMultiIndexTst));

        updateStates(leftOnsetMTr);
        updateStates(rightOnsetMTr);
        updateStates(centerOnsetMTr);

        updateStates(leftOnsetMTst);
        updateStates(rightOnsetMTst);
        updateStates(centerOnsetMTst);

        % I prefer this over SetStim so I can control index
        stimTr(1,1) = leftOnsetMTr;
        stimTr(1,2) = rightOnsetMTr;
        if numClasses == 3
            stimTr(1,3) = centerOnsetMTr;
        end

        stimTst(1,1) = leftOnsetMTst;
        stimTst(1,2) = rightOnsetMTst;
        if numClasses == 3
            stimTst(1,3) = centerOnsetMTst;
        end

        [stimTr,~] = hmrR_StimRejection(dod,stimTr,tIncAuto,tIncMan,[-2  15]);

        dodBPFilt = hmrR_BandpassFilt(dod,0.01,0.5);

        dc = hmrR_OD2Conc(dodBPFilt,probe,[1  1  1]);

        % GLM here
        % hmrR_GLM_ssBeta
%         [~,~,~,dcNewTr,~,~,~,~,~,~,betaSS] = hmrR_GLM_HbT(dc,stimTr,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
%             [-2  15],1,1,[1.0 1.0 0.0 0.0 0.0 0.0],15,1,0,0);

        [~,~,~,dcNewTr,~,~,~,~,~,~,betaSS] = hmrR_GLM_MN(dc,stimTr,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
            [-2  15],1,1,[1.0 1.0 0.0 0.0 0.0 0.0],15,1,0,0);

        % format beta var here
        ssBeta = betaSS{1}(end,:,:);

        [stimTst,~] = hmrR_StimRejection(dod,stimTst,tIncAuto,tIncMan,[-2  15]);

        % actually this is not needed, can use stimTst on dcNewTr.
        dcNewTst = hmrR_ssBeta_CV(data,dc, probe, mlActAuto, tIncAuto, squeeze(ssBeta));

%         allSMultipleTr = indexMoviesTrainMultiple(:,6);
%         allSMultipleTst = indexMoviesTestMultiple(:,6);
        
        
        % reject trial, stimTr into allSMultipleTr and
        
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
%         [performanceLDAHbO(foldIdx,:,:)] = ...
%             calcCumSum_DiffStartT_SingleChn_NoSave(sbjNum,trials_HbOM_Tr,trials_HbOM_Tst,...
%             indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlActAuto);
%         
%         [performanceLDAHbR(foldIdx,:,:)] = ...
%             calcCumSum_DiffStartT_SingleChn_NoSave(sbjNum,trials_HbRM_Tr,trials_HbRM_Tst,...
%             indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlActAuto);
%         
%         [performanceLDAHbT(foldIdx,:,:)] = ...
%             calcCumSum_DiffStartT_SingleChn_NoSave(sbjNum,trials_HbTM_Tr,trials_HbTM_Tst,...
%             indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlActAuto);
        
        performanceLDAHbT(foldIdx,:,:) = calcCumSum_DiffStartT_SingleChnHbX_NoSave(sbjNum,...
            trials_HbOM_Tr,trials_HbOM_Tst,trials_HbRM_Tr,trials_HbRM_Tst,trials_HbTM_Tr,trials_HbTM_Tst,...
            indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlActAuto);

    end
    
end

% save folds to file

if ~exist(processedDataDir,'dir')
    mkdir(processedDataDir);
end

if numClasses == 2
    if rejTrOp
        fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_SingleChnHbX_LR_RejTr_SNR_1.5_OrigHomer3.mat'];
        %fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_LR_RejTr_SNR_1.5_NewHomer3.mat'];
    else
        fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_SingleChnHbX_LR_SNR_1.5.mat'];
    end
else
    if rejTrOp
        fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_SingleChnHbX_RejTr_SNR_1.5.mat'];
    else
        fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_SingleChnHbX_SNR_1.5.mat'];
    end
end

save(fileName,'performanceLDAHbO','performanceLDAHbR','performanceLDAHbT');
    
end