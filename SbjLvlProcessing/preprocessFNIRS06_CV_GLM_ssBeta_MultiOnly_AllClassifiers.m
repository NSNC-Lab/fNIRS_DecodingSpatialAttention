% CV where SS beta coefficients from training fold is passed to test fold

% Only for sbj 12 and after

% Add more classifiers to validate logistic regression and svm unusual
% results. Only test one time point to expedite process. Not for
% publication

% Trim start and end T. Apply to dataTimeSeries and dataTime.

function preprocessFNIRS06_CV_GLM_ssBeta_MultiOnly_AllClassifiers(s,numClasses)

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
snirf1.Info()

% Extract aux and convert to stimclass
allS2 = find(snirf1.aux(1,1).dataTimeSeries>1);
aInd2 = find(allS2(2:end)-allS2(1:end-1)==1);

allS2(aInd2) = [];
allS2 = allS2./50;

if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
    cueOnsetIndex = 1:4:720;
else
    cueOnsetIndex = 1:4:360;
end

allS = allS2;

if startT ~= 1
    % allS is in sec
    allS = allS-startT;
end

fs = 50;
timePt = (0:0.25*fs:5*fs)+2*fs;

load([saveDir filesep movieList '.mat'],'indexMoviesTest');
load([saveDir filesep behData '.mat'],'responsesA','responsesV','correctRespA','correctRespV');

% dAll = DataClass([snirf1.data.dataTimeSeries; snirf3.data.dataTimeSeries; snirf4.data.dataTimeSeries],...
%     0:1/50:(size(snirf1.data.dataTimeSeries,1)+size(snirf3.data.dataTimeSeries,1)+size(snirf4.data.dataTimeSeries,1)-1)/50,snirf1.data.measurementList);

dAll = DataClass(snirf1.data.dataTimeSeries,...
    0:1/50:(size(snirf1.data.dataTimeSeries,1)-1)/50,snirf1.data.measurementList);

%obj = SnirfClass(data, stim, probe, aux);

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
%dod = snirf1.dod;
%mlActAuto = snirf1.mlActAuto;

%tIncAuto = snirf1.tIncAuto;
%dc = snirf1.dc;
Aaux = [];
rcMap = [];

%mlActAuto = hmrR_PruneChannels(data,probe,mlActMan,tIncMan,[10000  10000000],2,[0  45]);
mlActAuto = {ones(size(snirf1.data.measurementList,2),1)};
dod = hmrR_Intensity2OD(data);

% For sbj 12, params in 1st line for motion correct/artifact perform
% significantly better than params in 2nd line.
[dod] = hmrR_MotionCorrectPCArecurse(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,15,4,0.97,5,1);
%[dod,tInc,svs,nSV,tInc0] = hmrR_MotionCorrectPCArecurse(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,20,5,0.97,5,1);

tIncAuto = hmrR_MotionArtifact(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,15,5);
%tIncAuto = hmrR_MotionArtifact(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,20,4);

% 90 trials
trialIdx = (responsesA==correctRespA)&(responsesV==correctRespV)&(indexMoviesTest(:,5)==1)';
if numClasses == 2
    trialIdx = trialIdx & (indexMoviesTest(:,2)==1|indexMoviesTest(:,2)==2)';
end
if strcmp(sbjNum,'15')
    trialIdx(1) = 0;
end
numCorrectTrials = sum(trialIdx);
indexMoviesMultiple = indexMoviesTest(trialIdx,:);
cueOnsetIndexMultiple = cueOnsetIndex(trialIdx);

% Initialize variables
KFolds = 5;
nReps = 10;

% performanceLDALedoitHbO = zeros(nReps*KFolds,length(timePt));
% performanceLDACERNNHbO = zeros(nReps*KFolds,length(timePt));
% performanceLogRegHbO = zeros(nReps*KFolds,length(timePt));
% performanceSVMHbO = zeros(nReps*KFolds,length(timePt));
% performanceBaggingHbO = zeros(nReps*KFolds,length(timePt));
% performanceBoostedHbO = zeros(nReps*KFolds,length(timePt));
% 
% performanceLDALedoitHbR = zeros(nReps*KFolds,length(timePt));
% performanceLDACERNNHbR = zeros(nReps*KFolds,length(timePt));
% performanceLogRegHbR = zeros(nReps*KFolds,length(timePt));
% performanceSVMHbR = zeros(nReps*KFolds,length(timePt));
% performanceBaggingHbR = zeros(nReps*KFolds,length(timePt));
% performanceBoostedHbR = zeros(nReps*KFolds,length(timePt));
% 
% performanceLDALedoitHbT = zeros(nReps*KFolds,length(timePt));
% performanceLDACERNNHbT = zeros(nReps*KFolds,length(timePt));
% performanceLogRegHbT = zeros(nReps*KFolds,length(timePt));
% performanceSVMHbT = zeros(nReps*KFolds,length(timePt));
% performanceBaggingHbT = zeros(nReps*KFolds,length(timePt));
% performanceBoostedHbT = zeros(nReps*KFolds,length(timePt));

% get number of correct trials and their index


for iRep = 1:nReps
    
    % rejected trials? Not 90 or 60.
    cvp = cvpartition(numCorrectTrials, 'KFold',KFolds);
    
    for iFold = 1:KFolds

        foldIdx = (iRep-1)*5+iFold;
        
        % snirf 2. Cell of StimClass array.
        leftOnsetMTr = StimClass('leftMulti');
        rightOnsetMTr = StimClass('rightMulti');
        centerOnsetMTr = StimClass('centerMulti');

        leftOnsetMTst = StimClass('leftMulti');
        rightOnsetMTst = StimClass('rightMulti');
        centerOnsetMTst = StimClass('centerMulti');

        % training size (72 trials)
        indexMoviesTrainMultiple = indexMoviesMultiple(cvp.training(iFold),:);

        % list, cueOnsetIndexSingleTr is training size (72 trials).
        % cueOnsetIndexSingle is single condition size (90 trials).
        cueOnsetIndexMultipleTr = cueOnsetIndexMultiple(cvp.training(iFold));

        indexMoviesTestMultiple = indexMoviesMultiple(cvp.test(iFold),:);
        cueOnsetIndexMultipleTst = cueOnsetIndexMultiple(cvp.test(iFold));
        %allSTst = allS(cvp.test(iFold));

        % index of training size (72 trials), boolean for left condition
        leftMultiIndexTr = indexMoviesTrainMultiple(:,2)==2 & indexMoviesTrainMultiple(:,5)==1;
        rightMultiIndexTr = indexMoviesTrainMultiple(:,2)==1 & indexMoviesTrainMultiple(:,5)==1;
        centerMultiIndexTr = indexMoviesTrainMultiple(:,2)==3 & indexMoviesTrainMultiple(:,5)==1;

        leftMultiIndexTst = indexMoviesTestMultiple(:,2)==2 & indexMoviesTestMultiple(:,5)==1;
        rightMultiIndexTst = indexMoviesTestMultiple(:,2)==1 & indexMoviesTestMultiple(:,5)==1;
        centerMultiIndexTst = indexMoviesTestMultiple(:,2)==3 & indexMoviesTestMultiple(:,5)==1;

        % Only correct trials
        AddStims(leftOnsetMTr, allS(cueOnsetIndexMultipleTr(leftMultiIndexTr)));
        AddStims(rightOnsetMTr, allS(cueOnsetIndexMultipleTr(rightMultiIndexTr)));
        AddStims(centerOnsetMTr, allS(cueOnsetIndexMultipleTr(centerMultiIndexTr)));

        AddStims(leftOnsetMTst, allS(cueOnsetIndexMultipleTst(leftMultiIndexTst)));
        AddStims(rightOnsetMTst, allS(cueOnsetIndexMultipleTst(rightMultiIndexTst)));
        AddStims(centerOnsetMTst, allS(cueOnsetIndexMultipleTst(centerMultiIndexTst)));

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
        [~,~,~,dcNewTr,~,~,~,~,~,~,beta] = hmrR_GLM_ssBeta(dc,stimTr,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
            [-2  15],1,1,[1.0 1.0 0.0 0.0 0.0 0.0],15,1,0,0);

        % format beta var here
        ssBeta = beta{1}(end,:,:);

        [stimTst,~] = hmrR_StimRejection(dod,stimTst,tIncAuto,tIncMan,[-2  15]);

        dcNewTst = hmrR_ssBeta_CV(data,dc, probe, mlActAuto, tIncAuto, squeeze(ssBeta));

        % extract and format training and test trials
        allSMultipleTr = allS(cueOnsetIndexMultipleTr);
        allSMultipleTst = allS(cueOnsetIndexMultipleTst);
        
%         if numClasses == 2
%             idxM = indexMoviesTrainMultiple(:,2)==2 | indexMoviesTrainMultiple(:,2)==1;
%             allSMultipleTr = allS(cueOnsetIndexMultipleTr(idxM));
%             
%             idxM = indexMoviesTestMultiple(:,2)==2 | indexMoviesTestMultiple(:,2)==1;
%             allSMultipleTst = allS(cueOnsetIndexMultipleTst(idxM));
%         else
%             allSMultipleTr = allS(cueOnsetIndexMultipleTr);
%             allSMultipleTst = allS(cueOnsetIndexMultipleTst);
%         end
        
        [trials_HbOM_Tr, trials_HbRM_Tr, trials_HbTM_Tr] ...
            = createSingleTrialHRF_MultiOnly_NoSave(dcNewTr, allSMultipleTr);
        
        [trials_HbOM_Tst, trials_HbRM_Tst, trials_HbTM_Tst] ...
            = createSingleTrialHRF_MultiOnly_NoSave(dcNewTst, allSMultipleTst);

        % train classifier
        performanceStrHbOTemp = calcCumSum_DiffStartT_AllClassifiers(...
            sbjNum,trials_HbOM_Tr,trials_HbOM_Tst,indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlActAuto);
        
        performanceStrHbRTemp = calcCumSum_DiffStartT_AllClassifiers(...
            sbjNum,trials_HbRM_Tr,trials_HbRM_Tst,indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlActAuto);
        
        performanceStrHbTTemp = calcCumSum_DiffStartT_AllClassifiers(...
            sbjNum,trials_HbTM_Tr,trials_HbTM_Tst,indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlActAuto);
        
        performanceStrHbO.performanceLDALedoit(foldIdx,:) = performanceStrHbOTemp.performanceLDALedoit;
        performanceStrHbO.performanceLDACERNN(foldIdx,:) = performanceStrHbOTemp.performanceLDACERNN;
        performanceStrHbO.performanceLogReg(foldIdx,:) = performanceStrHbOTemp.performanceLogReg;
        performanceStrHbO.performanceSVM(foldIdx,:) = performanceStrHbOTemp.performanceSVM;
        performanceStrHbO.performanceBagging(foldIdx,:) = performanceStrHbOTemp.performanceBagging;
        performanceStrHbO.performanceBoosted(foldIdx,:) = performanceStrHbOTemp.performanceBoosted;

        performanceStrHbO.performanceCoarseKNN(foldIdx,:) = performanceStrHbOTemp.performanceCoarseKNN;
        performanceStrHbO.performanceCoarseTree(foldIdx,:) = performanceStrHbOTemp.performanceCoarseTree;
        performanceStrHbO.performanceCosineKNN(foldIdx,:) = performanceStrHbOTemp.performanceCosineKNN;
        performanceStrHbO.performanceCubicKNN(foldIdx,:) = performanceStrHbOTemp.performanceCubicKNN;
        performanceStrHbO.performanceCubicSVM(foldIdx,:) = performanceStrHbOTemp.performanceCubicSVM;
        performanceStrHbO.performanceFineGaussianSVM(foldIdx,:) = performanceStrHbOTemp.performanceFineGaussianSVM;
        performanceStrHbO.performanceFineKNN(foldIdx,:) = performanceStrHbOTemp.performanceFineKNN;
        performanceStrHbO.performanceFineTree(foldIdx,:) = performanceStrHbOTemp.performanceFineTree;
        performanceStrHbO.performanceGaussianNaiveBayes(foldIdx,:) = performanceStrHbOTemp.performanceGaussianNaiveBayes;
        performanceStrHbO.performanceKernelNaiveBayes(foldIdx,:) = performanceStrHbOTemp.performanceKernelNaiveBayes;
        performanceStrHbO.performanceWeightedKNN(foldIdx,:) = performanceStrHbOTemp.performanceWeightedKNN;
        performanceStrHbO.performanceRandom(foldIdx,:) = performanceStrHbOTemp.performanceRandom;
        
        performanceStrHbR.performanceLDALedoit(foldIdx,:) = performanceStrHbRTemp.performanceLDALedoit;
        performanceStrHbR.performanceLDACERNN(foldIdx,:) = performanceStrHbRTemp.performanceLDACERNN;
        performanceStrHbR.performanceLogReg(foldIdx,:) = performanceStrHbRTemp.performanceLogReg;
        performanceStrHbR.performanceSVM(foldIdx,:) = performanceStrHbRTemp.performanceSVM;
        performanceStrHbR.performanceBagging(foldIdx,:) = performanceStrHbRTemp.performanceBagging;
        performanceStrHbR.performanceBoosted(foldIdx,:) = performanceStrHbRTemp.performanceBoosted;

        performanceStrHbR.performanceCoarseKNN(foldIdx,:) = performanceStrHbRTemp.performanceCoarseKNN;
        performanceStrHbR.performanceCoarseTree(foldIdx,:) = performanceStrHbRTemp.performanceCoarseTree;
        performanceStrHbR.performanceCosineKNN(foldIdx,:) = performanceStrHbRTemp.performanceCosineKNN;
        performanceStrHbR.performanceCubicKNN(foldIdx,:) = performanceStrHbRTemp.performanceCubicKNN;
        performanceStrHbR.performanceCubicSVM(foldIdx,:) = performanceStrHbRTemp.performanceCubicSVM;
        performanceStrHbR.performanceFineGaussianSVM(foldIdx,:) = performanceStrHbRTemp.performanceFineGaussianSVM;
        performanceStrHbR.performanceFineKNN(foldIdx,:) = performanceStrHbRTemp.performanceFineKNN;
        performanceStrHbR.performanceFineTree(foldIdx,:) = performanceStrHbRTemp.performanceFineTree;
        performanceStrHbR.performanceGaussianNaiveBayes(foldIdx,:) = performanceStrHbRTemp.performanceGaussianNaiveBayes;
        performanceStrHbR.performanceKernelNaiveBayes(foldIdx,:) = performanceStrHbRTemp.performanceKernelNaiveBayes;
        performanceStrHbR.performanceWeightedKNN(foldIdx,:) = performanceStrHbRTemp.performanceWeightedKNN;
        performanceStrHbR.performanceRandom(foldIdx,:) = performanceStrHbRTemp.performanceRandom;
        
        performanceStrHbT.performanceLDALedoit(foldIdx,:) = performanceStrHbTTemp.performanceLDALedoit;
        performanceStrHbT.performanceLDACERNN(foldIdx,:) = performanceStrHbTTemp.performanceLDACERNN;
        performanceStrHbT.performanceLogReg(foldIdx,:) = performanceStrHbTTemp.performanceLogReg;
        performanceStrHbT.performanceSVM(foldIdx,:) = performanceStrHbTTemp.performanceSVM;
        performanceStrHbT.performanceBagging(foldIdx,:) = performanceStrHbTTemp.performanceBagging;
        performanceStrHbT.performanceBoosted(foldIdx,:) = performanceStrHbTTemp.performanceBoosted;

        performanceStrHbT.performanceCoarseKNN(foldIdx,:) = performanceStrHbTTemp.performanceCoarseKNN;
        performanceStrHbT.performanceCoarseTree(foldIdx,:) = performanceStrHbTTemp.performanceCoarseTree;
        performanceStrHbT.performanceCosineKNN(foldIdx,:) = performanceStrHbTTemp.performanceCosineKNN;
        performanceStrHbT.performanceCubicKNN(foldIdx,:) = performanceStrHbTTemp.performanceCubicKNN;
        performanceStrHbT.performanceCubicSVM(foldIdx,:) = performanceStrHbTTemp.performanceCubicSVM;
        performanceStrHbT.performanceFineGaussianSVM(foldIdx,:) = performanceStrHbTTemp.performanceFineGaussianSVM;
        performanceStrHbT.performanceFineKNN(foldIdx,:) = performanceStrHbTTemp.performanceFineKNN;
        performanceStrHbT.performanceFineTree(foldIdx,:) = performanceStrHbTTemp.performanceFineTree;
        performanceStrHbT.performanceGaussianNaiveBayes(foldIdx,:) = performanceStrHbTTemp.performanceGaussianNaiveBayes;
        performanceStrHbT.performanceKernelNaiveBayes(foldIdx,:) = performanceStrHbTTemp.performanceKernelNaiveBayes;
        performanceStrHbT.performanceWeightedKNN(foldIdx,:) = performanceStrHbTTemp.performanceWeightedKNN;
        performanceStrHbT.performanceRandom(foldIdx,:) = performanceStrHbTTemp.performanceRandom;
        
        
        
        
        % train classifier
%         [performanceLDALedoitHbO(foldIdx,:), performanceLDACERNNHbO(foldIdx,:),...
%             performanceLogRegHbO(foldIdx,:), performanceSVMHbO(foldIdx,:),...
%             performanceBaggingHbO(foldIdx,:), performanceBoostedHbO(foldIdx,:)] = ...
%             calcCumSum_DiffStartT_AllChn_CompClassifiers_Final_NoSave(sbjNum,trials_HbOM_Tr,trials_HbOM_Tst,...
%             indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlActAuto);
%         
%         [performanceLDALedoitHbR(foldIdx,:), performanceLDACERNNHbR(foldIdx,:),...
%             performanceLogRegHbR(foldIdx,:), performanceSVMHbR(foldIdx,:),...
%             performanceBaggingHbR(foldIdx,:), performanceBoostedHbR(foldIdx,:)] = ...
%             calcCumSum_DiffStartT_AllChn_CompClassifiers_Final_NoSave(sbjNum,trials_HbRM_Tr,trials_HbRM_Tst,...
%             indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlActAuto);
%         
%         [performanceLDALedoitHbT(foldIdx,:), performanceLDACERNNHbT(foldIdx,:),...
%             performanceLogRegHbT(foldIdx,:), performanceSVMHbT(foldIdx,:),...
%             performanceBaggingHbT(foldIdx,:), performanceBoostedHbT(foldIdx,:)] = ...
%             calcCumSum_DiffStartT_AllChn_CompClassifiers_Final_NoSave(sbjNum,trials_HbTM_Tr,trials_HbTM_Tst,...
%             indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlActAuto);

    end
    
    
end

% save folds to file

if ~exist(processedDataDir,'dir')
    mkdir(processedDataDir);
end

if numClasses == 2
    fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_LR_AllClassifiers.mat'];
else
    fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_AllClassifiers.mat'];
end
save(fileName,'performanceStrHbO','performanceStrHbR','performanceStrHbT');
    
end