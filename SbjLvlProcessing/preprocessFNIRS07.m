% GLM Cross validation. 2 Rounds of GLM Fitting.
% For sbj 08 and 10.

function preprocessFNIRS07(sbjNum,rawDataFN,respData,numClasses)
saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
rawDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
figSaveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' ...
    num2str(sbjNum) '\Figures\Classification'];
behFN = [rawDataDir filesep 'responses_' sbjNum];

%numClasses = 3;
fs = 50;
timePt = (0:0.5*fs:6*fs)+2*fs;
trialTime = [-2 15];

if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
    isSDNew = 0;
else
    isSDNew = 1;
end

% Convert nirs to snirf file format
snirf2 = SnirfClass(load([rawDataFN{1} '.nirs'],'-mat'));
snirf3 = SnirfClass(load([rawDataFN{2} '.nirs'],'-mat'));
snirf4 = SnirfClass(load([rawDataFN{3} '.nirs'],'-mat'));
snirf2.Info()
snirfTemp = SnirfClass();

% Extract aux and convert to stimclass
allS2 = find(snirf2.aux(1,1).dataTimeSeries>1);
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

oldStimClass1Length2 = size(allS2,1);
oldStimClass1Length3 = oldStimClass1Length2 + size(allS3,1);
oldStimClass1Length4 = oldStimClass1Length3 + size(allS4,1);

cueOnsetIndex = 1:4:720;
cueOnsetIndex2 = 1:oldStimClass1Length2;
cueOnsetIndex3 = oldStimClass1Length2+1:oldStimClass1Length3;
cueOnsetIndex4 = oldStimClass1Length3+1:oldStimClass1Length4;

trialNum2 = sum(ismember(cueOnsetIndex2,cueOnsetIndex));
trialNum3 = sum(ismember(cueOnsetIndex3,cueOnsetIndex));
trialNum4 = sum(ismember(cueOnsetIndex4,cueOnsetIndex));

lastTrial2 = trialNum2;
lastTrial3 = lastTrial2 + trialNum3;
lastTrial4 = lastTrial3 + trialNum4;

snirf2TLen = size(snirf2.data.time,1)/(50);
snirf3TLen = size(snirf3.data.time,1)/(50);
allS3 = allS3+snirf2TLen;
allS4 = allS4+snirf2TLen+snirf3TLen;

allS = [allS2; allS3; allS4];

% split triggers into 6 categories
load([saveDir filesep respData '.mat'],'fixedMaskerList','indexMoviesTest','maskerMovies','numTrials','uniqueMovies');

%cv = cvpartition(size(allS,1),'KFold',KFolds);

indexMoviesTest = [indexMoviesTest (1:size(indexMoviesTest,1))'];

% snirf 2
leftOnsetS = StimClass('leftSingle');
rightOnsetS = StimClass('rightSingle');
centerOnsetS = StimClass('centerSingle');
leftOnsetM = StimClass('leftMulti');
rightOnsetM = StimClass('rightMulti');
centerOnsetM = StimClass('centerMulti');

leftSingleIndex = indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==0;
rightSingleIndex = indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==0;
centerSingleIndex = indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==0;
leftMultiIndex = indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1;
rightMultiIndex = indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1;
centerMultiIndex = indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==1;

% leftSingleIndex = leftSingleIndex(leftSingleIndex<= lastTrial2);
% rightSingleIndex = rightSingleIndex(rightSingleIndex<= lastTrial2);
% centerSingleIndex = centerSingleIndex(centerSingleIndex<= lastTrial2);
% leftMultiIndex = leftMultiIndex(leftMultiIndex<= lastTrial2);
% rightMultiIndex = rightMultiIndex(rightMultiIndex<= lastTrial2);
% centerMultiIndex = centerMultiIndex(centerMultiIndex<= lastTrial2);

AddStims(leftOnsetS, allS(cueOnsetIndex(leftSingleIndex)));
AddStims(rightOnsetS, allS(cueOnsetIndex(rightSingleIndex)));
AddStims(centerOnsetS, allS(cueOnsetIndex(centerSingleIndex)));
AddStims(leftOnsetM, allS(cueOnsetIndex(leftMultiIndex)));
AddStims(rightOnsetM, allS(cueOnsetIndex(rightMultiIndex)));
AddStims(centerOnsetM, allS(cueOnsetIndex(centerMultiIndex)));

updateStates(leftOnsetS);
updateStates(rightOnsetS);
updateStates(centerOnsetS);
updateStates(leftOnsetM);
updateStates(rightOnsetM);
updateStates(centerOnsetM);

snirf2.stim(1,1) = leftOnsetS;
snirf2.stim(1,2) = rightOnsetS;
snirf2.stim(1,3) = centerOnsetS;
snirf2.stim(1,4) = leftOnsetM;
snirf2.stim(1,5) = rightOnsetM;
snirf2.stim(1,6) = centerOnsetM;

% combine snirf2.data.dataTimeSeries and time
snirf2.data.dataTimeSeries = [snirf2.data.dataTimeSeries; snirf3.data.dataTimeSeries; snirf4.data.dataTimeSeries];
%snirf2.data.time = [snirf2.data.time;snirf3.data.time;snirf4.data.time];
snirf2.data.time = 0:1/50:(size(snirf2.data.dataTimeSeries,1)-1)/50;

data = snirf2.data;
probe = snirf2.probe;
mlActMan = {};
tIncMan = {};
%dod = snirf1.dod;
%mlActAuto = snirf1.mlActAuto;
stim = snirf2.stim;
%tIncAuto = snirf1.tIncAuto;
%dc = snirf1.dc;
Aaux = [];
rcMap = [];

%mlActAuto = hmrR_PruneChannels(data,probe,mlActMan,tIncMan,[10000  10000000],2,[0  45]);
mlActAuto = {ones(size(snirf2.data.measurementList,2),1)};

dod = hmrR_Intensity2OD(data);

[dod,tInc,svs,nSV,tInc0] = hmrR_MotionCorrectPCArecurse(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,15,4,0.97,5,1);

tIncAuto = hmrR_MotionArtifact(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,15,5);
 
[stim,tRange] = hmrR_StimRejection(dod,stim,tIncAuto,tIncMan,[-2  15]);

dod = hmrR_BandpassFilt(dod,0.01,0.5);
dc = hmrR_OD2Conc(dod,probe,[1  1  1]); 

chnNum = size(dod.dataTimeSeries,2)/2;
tLen = size(dod.dataTimeSeries,1);

zeroT = 2*fs;


% channel X HRF time X trials x trials
% dcTestAvg = zeros(chnNum,length(-2*fs:15*fs-1),size(indexMoviesTest,1),size(indexMoviesTest,1));
% dcTestAvgStd = zeros(chnNum,length(-2*fs:15*fs-1),size(indexMoviesTest,1),size(indexMoviesTest,1));
% % trials X trials X chnnels X whole data time
% dcTestNew  = zeros(size(indexMoviesTest,1),size(indexMoviesTest,1),chnNum,tLen);
% dcTestResid  = zeros(size(indexMoviesTest,1),size(indexMoviesTest,1),chnNum,tLen);
if strcmp(sbjNum,'15')
    indexMoviesTest = indexMoviesTest(2:end,:);
end

% Filter to only correct trials
load(behFN,'responsesV','responsesA','correctRespV','correctRespA');
idx = (responsesA==correctRespA).*(responsesV==correctRespV);

rejTr = find(~idx);

dod.dataTimeSeries((allS(cueOnsetIndex(rejTr),1)+trialTime(1))*fs:(allS(cueOnsetIndex(rejTr),1)+trialTime(2))*fs,:) = 0;

origIdxMoviesTest = indexMoviesTest;
indexMoviesTest = indexMoviesTest(logical(idx),:);

if numClasses==2
    idx = (indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1) | (indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1);
    indexMoviesTestCond = indexMoviesTest(idx,:);
    numTr = sum(logical(idx));
    %indexMoviesTestCond(:,6)=1:size(indexMoviesTestCond,1);
else
    indexMoviesTestCond = indexMoviesTest(indexMoviesTest(:,5)==indexMoviesTest(i1,5),:);
    numTr = sum(logical(idx));
    %indexMoviesTestCond(:,6)=1:size(indexMoviesTestCond,1);
end

% if numClasses==2
%     numTrOrig = size(indexMoviesTestCond,1)/3;
% else
%     numTrOrig = size(indexMoviesTestCond,1)/2;
% end

if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
    origNumChn = 42;
else
    origNumChn = 36;
end

% Array too big, out of memory issue
% Don't save, just input to classification directly.
dcTestAvg = zeros(origNumChn*3,length(trialTime(1)*fs:trialTime(2)*fs),numTr);
dcTestAvgHbOOrig = zeros(origNumChn,length(trialTime(1)*fs:trialTime(2)*fs),numTr);
dcTestAvgHbROrig = zeros(origNumChn,length(trialTime(1)*fs:trialTime(2)*fs),numTr);
dcTestAvgHbTOrig = zeros(origNumChn,length(trialTime(1)*fs:trialTime(2)*fs),numTr);

performanceArr_LDA_LW_HbO = zeros(numTr,length(timePt));
performanceArr_LDA_CERNN_HbO = zeros(numTr,length(timePt));
performanceArr_Logistics_HbO = zeros(numTr,length(timePt));
performanceArr_SVM_Lin_HbO = zeros(numTr,length(timePt));
performanceArr_Bagging_HbO = zeros(numTr,length(timePt));
performanceArr_Boosting_HbO = zeros(numTr,length(timePt));

performanceArr_LDA_LW_HbR = zeros(numTr,length(timePt));
performanceArr_LDA_CERNN_HbR = zeros(numTr,length(timePt));
performanceArr_Logistics_HbR = zeros(numTr,length(timePt));
performanceArr_SVM_Lin_HbR = zeros(numTr,length(timePt));
performanceArr_Bagging_HbR = zeros(numTr,length(timePt));
performanceArr_Boosting_HbR = zeros(numTr,length(timePt));

performanceArr_LDA_LW_HbT = zeros(numTr,length(timePt));
performanceArr_LDA_CERNN_HbT = zeros(numTr,length(timePt));
performanceArr_Logistics_HbT = zeros(numTr,length(timePt));
performanceArr_SVM_Lin_HbT = zeros(numTr,length(timePt));
performanceArr_Bagging_HbT = zeros(numTr,length(timePt));
performanceArr_Boosting_HbT = zeros(numTr,length(timePt));

param.tstIdx = zeros(1,numTr);
param.trnIdx = zeros(numTr,numTr-1);

% After removing incorrect trials, only multiple movies condition
indexMoviesTest = indexMoviesTest(idx,:);
cueOnsetIndexOrig = cueOnsetIndex;
cueOnsetIndex = cueOnsetIndex(idx);
%allS = allS(cueOnsetIndex);

% Start looping folds.
for i1 = 1:size(indexMoviesTest,1)
%for i1 = 1:1
    % These are index of indexMoviesTestCond. Sixth col of
    % indexMoviesTestCond is original trial #.
    param.tstIdx(i1)=find(indexMoviesTestCond(:,6) == indexMoviesTest(i1,6));
    param.trnIdx(i1,:) = setdiff(1:size(indexMoviesTestCond,1), find(indexMoviesTestCond(:,6) == indexMoviesTest(i1,6)));
    %param.trnIdx(i1,:) = setdiff(indexMoviesTestCond(:,6), find(indexMoviesTestCond(:,6) == indexMoviesTest(i1,6)));
    %param.trnIdx(i1,:) = setdiff(1:size(indexMoviesTestCond,1), find(indexMoviesTestCond(:,6) == indexMoviesTest(i1,6)));
    %param.trnIdx(i1,:) = setdiff(indexMoviesTestCond(:,6),i1);
    trNum = i1;
    
    snirfTest = SnirfClass();
    dodTr = dod;
    leftOnsetS = StimClass('leftSingle');
    rightOnsetS = StimClass('rightSingle');
    centerOnsetS = StimClass('centerSingle');
    leftOnsetM = StimClass('leftMulti');
    rightOnsetM = StimClass('rightMulti');
    centerOnsetM = StimClass('centerMulti');

    indexMoviesTestTemp = indexMoviesTest;
    allSTemp = allS;
    cueOnsetIndexTemp = cueOnsetIndex;
%     indexMoviesTestTemp(i1,:) = [];
%     allSTemp(i1,:) = [];
%     cueOnsetIndexTemp(i1) = [];
    
    leftSingleIndex = indexMoviesTestTemp(:,2)==2 & indexMoviesTestTemp(:,5)==0;
    rightSingleIndex = indexMoviesTestTemp(:,2)==1 & indexMoviesTestTemp(:,5)==0;
    centerSingleIndex = indexMoviesTestTemp(:,2)==3 & indexMoviesTestTemp(:,5)==0;
    leftMultiIndex = indexMoviesTestTemp(:,2)==2 & indexMoviesTestTemp(:,5)==1;
    rightMultiIndex = indexMoviesTestTemp(:,2)==1 & indexMoviesTestTemp(:,5)==1;
    centerMultiIndex = indexMoviesTestTemp(:,2)==3 & indexMoviesTestTemp(:,5)==1;
    
%     leftSingleIndex = indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==0 & cv.training(i)==1;
%     rightSingleIndex = indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==0 & cv.training(i)==1;
%     centerSingleIndex = indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==0 & cv.training(i)==1;
%     leftMultiIndex = indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1 & cv.training(i)==1;
%     rightMultiIndex = indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1 & cv.training(i)==1;
%     centerMultiIndex = indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==1 & cv.training(i)==1;
    
%     leftSingleIndexTest = indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==0 & cv.test(i)==1;
%     rightSingleIndexTest = indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==0 & cv.test(i)==1;
%     centerSingleIndexTest = indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==0 & cv.test(i)==1;
%     leftMultiIndexTest = indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1 & cv.test(i)==1;
%     rightMultiIndexTest = indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1 & cv.test(i)==1;
%     centerMultiIndexTest = indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==1 & cv.test(i)==1;
%     % Not needed, do single trial (leave one out)
%     trainingStim = allS(cv.training(i));
%     testStim = allS(cv.test(i));

    % leftSingleIndex = leftSingleIndex(leftSingleIndex<= lastTrial2);
    % rightSingleIndex = rightSingleIndex(rightSingleIndex<= lastTrial2);
    % centerSingleIndex = centerSingleIndex(centerSingleIndex<= lastTrial2);
    % leftMultiIndex = leftMultiIndex(leftMultiIndex<= lastTrial2);
    % rightMultiIndex = rightMultiIndex(rightMultiIndex<= lastTrial2);
    % centerMultiIndex = centerMultiIndex(centerMultiIndex<= lastTrial2);

    AddStims(leftOnsetS, allSTemp(cueOnsetIndexTemp(leftSingleIndex)));
    AddStims(rightOnsetS, allSTemp(cueOnsetIndexTemp(rightSingleIndex)));
    AddStims(centerOnsetS, allSTemp(cueOnsetIndexTemp(centerSingleIndex)));
    AddStims(leftOnsetM, allSTemp(cueOnsetIndexTemp(leftMultiIndex)));
    AddStims(rightOnsetM, allSTemp(cueOnsetIndexTemp(rightMultiIndex)));
    AddStims(centerOnsetM, allSTemp(cueOnsetIndexTemp(centerMultiIndex)));

    updateStates(leftOnsetS);
    updateStates(rightOnsetS);
    updateStates(centerOnsetS);
    updateStates(leftOnsetM);
    updateStates(rightOnsetM);
    updateStates(centerOnsetM);

    snirfTemp.stim(1,1) = leftOnsetS;
    snirfTemp.stim(1,2) = rightOnsetS;
    snirfTemp.stim(1,3) = centerOnsetS;
    snirfTemp.stim(1,4) = leftOnsetM;
    snirfTemp.stim(1,5) = rightOnsetM;
    snirfTemp.stim(1,6) = centerOnsetM;
    
    stimTr = snirfTemp.stim;
    
%     for i = 1:size(testStim,1)
%         dodTr.dataTimeSeries((testStim(i)-2)*fs:(testStim(i)+15)*fs,:) = 0;
%     end

%     for i = 1:size(rejTr,1)
%         %dodTr.dataTimeSeries((testStim(i)-2)*fs:(testStim(i)+15)*fs,:) = 0;
%         dodTr.dataTimeSeries((allSTemp(cueOnsetIndexTemp(i1),1)+trialTime(1))*fs:(allSTemp(cueOnsetIndexTemp(i1),1)+trialTime(2))*fs,:) = 0;
%     end
    %dodTr.dataTimeSeries((allSTemp(cueOnsetIndexOrig(rejTr),1)+trialTime(1))*fs:(allSTemp(cueOnsetIndexOrig(rejTr),1)+trialTime(2))*fs,:) = 0;
    %dodTr.dataTimeSeries((allSTemp(cueOnsetIndexTemp(rejTr),1)+trialTime(1))*fs:(allSTemp(cueOnsetIndexTemp(rejTr),1)+trialTime(2))*fs,:) = 0;

    dodTr.dataTimeSeries((allSTemp(cueOnsetIndexTemp(i1),1)+trialTime(1))*fs:(allSTemp(cueOnsetIndexTemp(i1),1)+trialTime(2))*fs,:) = 0;
    
    dcTr = hmrR_OD2Conc(dodTr,probe,[1  1  1]);

    % set segments of dc = 0 for testing block.

    [dcTrAvg,dcAvgStd,nTrials,dcNew,dcResid,dcSum2,beta,R,hmrstats] = ...
        hmrR_GLM(dcTr,stimTr,probe,mlActAuto,Aaux,tIncAuto,rcMap,trialTime,1,1,[1  1  0  0  0  0],15,1,3,0);

    iCounter = 0;
    
    for i2 = 1:size(indexMoviesTest,1)
    %for i2 = 1:1

        if indexMoviesTest(i2,5) == indexMoviesTest(i1,5)
            iCounter = iCounter + 1;
            indexMoviesTestTest = indexMoviesTest(i2,:);
            allSTempTest = allS;
            cueOnsetIndexTest = cueOnsetIndex;
            %indexMoviesTestTest = indexMoviesTestTest(i2,:);
            allSTempTest = allSTempTest(cueOnsetIndexTest(i2),:);

            if indexMoviesTestTest(:,2)==2 && indexMoviesTestTest(:,5)==0
                leftOnsetSTest = StimClass('leftSingle');
                AddStims(leftOnsetSTest,allSTempTest);
                updateStates(leftOnsetSTest);
                snirfTest.stim(1,1) = leftOnsetSTest;
                thisCond = 1;
            elseif indexMoviesTestTest(:,2)==1 && indexMoviesTestTest(:,5)==0
                rightOnsetSTest = StimClass('rightSingle');
                AddStims(rightOnsetSTest,allSTempTest);
                updateStates(rightOnsetSTest);
                snirfTest.stim(1,1) = rightOnsetSTest;
                thisCond = 2;
            elseif indexMoviesTestTest(:,2)==3 && indexMoviesTestTest(:,5)==0
                centerOnsetSTest = StimClass('centerSingle');
                AddStims(centerOnsetSTest,allSTempTest);
                updateStates(centerOnsetSTest);
                snirfTest.stim(1,1) = centerOnsetSTest;
                thisCond = 3;
            elseif indexMoviesTestTest(:,2)==2 && indexMoviesTestTest(:,5)==1
                leftOnsetMTest = StimClass('leftMulti');
                AddStims(leftOnsetMTest,allSTempTest);
                updateStates(leftOnsetMTest);
                snirfTest.stim(1,1) = leftOnsetMTest;
                thisCond = 4;
            elseif indexMoviesTestTest(:,2)==1 && indexMoviesTestTest(:,5)==1
                rightOnsetMTest = StimClass('rightMulti');
                AddStims(rightOnsetMTest,allSTempTest);
                updateStates(rightOnsetMTest);
                snirfTest.stim(1,1) = rightOnsetMTest;
                thisCond = 5;
            elseif indexMoviesTestTest(:,2)==3 && indexMoviesTestTest(:,5)==1
                centerOnsetMTest = StimClass('centerMulti');
                AddStims(centerOnsetMTest,allSTempTest);
                updateStates(centerOnsetMTest);
                snirfTest.stim(1,1) = centerOnsetMTest;
                thisCond = 6;
            end

            stimTest = snirfTest.stim;

            % time X chromophores X chns X cond
            % Example, when thisCond = 6
            % test = dcTrAvgReshape(:,2,2);
            % test2 = dcTrAvg.dataTimeSeries(:,5+126*5);
            dcTrAvgReshape = reshape(dcTrAvg.dataTimeSeries,size(dcTrAvg.dataTimeSeries,1),3,42,6);
            dcTrAvgReshape = dcTrAvgReshape(:,:,:,thisCond);

            dcTrial = DataClass(dc.dataTimeSeries((allSTemp(cueOnsetIndexTemp(i2),1)+trialTime(1))*fs:(allSTemp(cueOnsetIndexTemp(i2),1)+trialTime(2)+1)*fs,:),...
                dc.time(1,(allSTemp(cueOnsetIndexTemp(i2),1)+trialTime(1))*fs:(allSTemp(cueOnsetIndexTemp(i2),1)+trialTime(2)+1)*fs),...
                dc.measurementList);
            
            %dc.dataTimeSeries = dc.dataTimeSeries((allSTemp(cueOnsetIndexTemp(i2),1)-2)*fs:(allSTemp(cueOnsetIndexTemp(i2),1)+15)*fs,:);
            %dc.time = dc.time(1,(allSTemp(cueOnsetIndexTemp(i2),1)-2)*fs:(allSTemp(cueOnsetIndexTemp(i2),1)+15)*fs);
            tIncAutoTrial{1} = ones(size(dcTrial.time,2),1);
            %dc.dataTimeSeries((allSTemp(cueOnsetIndexTemp(i2),1)-2)*fs:(allSTemp(cueOnsetIndexTemp(i2),1)+15)*fs,:)= 0;

            %[dcAvgTemp,dcStdTemp,nTrials,dcNewTemp,dcResidTemp,dcSum2,beta,R,hmrstats] = ...
            dcAvgTemp = hmrR_GLM_CVTest(dcTrial,...
                stimTest,probe,mlActAuto,Aaux,tIncAutoTrial,rcMap,trialTime,1,5,dcTrAvgReshape,15,1,3,0);

            dcTestAvg(:,:,iCounter) = dcAvgTemp.dataTimeSeries';
            
            % Test this!!!
            hbOIdx = 1:3:size(dcTestAvg,1);
            % channels x time x trials
            dcTestAvgHbOOrig(:,:,iCounter) = dcTestAvg(hbOIdx,:,iCounter);
            %dcTestAvgHbO(2:end,:,iCounter) = dcTestAvg(hbOIdx,:,iCounter);
            %dcTestAvgHbO(1,:,iCounter) = indexMoviesTest(i2,2);
            
            hbRIdx = 2:3:size(dcTestAvg,1);
            dcTestAvgHbROrig(:,:,iCounter) = dcTestAvg(hbRIdx,:,iCounter);
            %dcTestAvgHbR(2:end,:,iCounter) = dcTestAvg(hbRIdx,:,iCounter);
            %dcTestAvgHbR(1,:,iCounter) = indexMoviesTest(i2,2);
            
            hbTIdx = 3:3:size(dcTestAvg,1);
            dcTestAvgHbTOrig(:,:,iCounter) = dcTestAvg(hbTIdx,:,iCounter);
            %dcTestAvgHbT(2:end,:,iCounter) = dcTestAvg(hbTIdx,:,iCounter);
            %dcTestAvgHbT(1,:,iCounter) = indexMoviesTest(i2,2);
        end
        
%         dcTestAvg(:,:,i1,i2) = dcAvgTemp;
%         dcTestAvgStd(:,:,i1,i2) = dcStdTemp;
%         dcTestNew(i1,i2,:) = dcNewTemp;
%         dcTestResid(i1,i2,:) = dcResidTemp;
        
    end
    % performanceArr = zeros(size(trials,1)-1,length(timePt));
    %[performanceArr(i1,:,:),~] = trainClassifierLinearDiscriminantfNIRSGLMCV(dcTestAvgHbO,timePt,mlActAuto{1},numClasses,param,trNum);
    
    % Offset
    singleTrialHRFHbOM = offsetTrials(dcTestAvgHbOOrig,zeroT);
    singleTrialHRFHbRM = offsetTrials(dcTestAvgHbROrig,zeroT);
    singleTrialHRFHbTM = offsetTrials(dcTestAvgHbTOrig,zeroT);

    % after convert2SD2, cut down from 42 chns to 36 chns. Remove old chns
    if ~isSDNew
        [singleTrialHRFHbOM,mlList] = convert2SD2(singleTrialHRFHbOM,mlActAuto{1});
        [singleTrialHRFHbRM,~] = convert2SD2(singleTrialHRFHbRM);
        [singleTrialHRFHbTM,~] = convert2SD2(singleTrialHRFHbTM);
    else
        mlList = mlActAuto{1};
    end

%     isSDNew = 1;

    % after selectRS, cut down from 36 to 30 chns. Remove SS
    [dcTestAvgHbO,mlList] = selectRS(singleTrialHRFHbOM,1,mlList);
    [dcTestAvgHbR,~] = selectRS(singleTrialHRFHbRM,1);
    [dcTestAvgHbT,~] = selectRS(singleTrialHRFHbTM,1);

%     [dcTestAvgHbO,movieIdx,behScore] = keepCorrectTrials_Param_Special(sbjNum,singleTrialHRFHbOM,behFN,numClasses,origIdxMoviesTest);
%     [dcTestAvgHbR] = keepCorrectTrials_Param_Special(sbjNum,singleTrialHRFHbRM,behFN,numClasses,origIdxMoviesTest);
%     [dcTestAvgHbT] = keepCorrectTrials_Param_Special(sbjNum,singleTrialHRFHbTM,behFN,numClasses,origIdxMoviesTest);
            
    numChn = 30;
            
    performanceArr_LDA_LW_HbO(i1,:) = train_RLDA_Ledoit_LOOCV(dcTestAvgHbO,mlActAuto{1},numChn,numClasses,timePt,indexMoviesTest,param,trNum);
    performanceArr_LDA_CERNN_HbO(i1,:) = train_RLDA_CERNN_LOOCV(dcTestAvgHbO,mlActAuto{1},numChn,numClasses,timePt,indexMoviesTest,param,trNum);
    performanceArr_Logistics_HbO(i1,:) = train_LogisticRegression_LOOCV(dcTestAvgHbO,mlActAuto{1},numChn,numClasses,timePt,indexMoviesTest,param,trNum);
    performanceArr_SVM_Lin_HbO(i1,:) = trainClassifier_LinearSVM_LOOCV(dcTestAvgHbO,mlActAuto{1},numChn,numClasses,timePt,indexMoviesTest,param,trNum);
    performanceArr_Bagging_HbO(i1,:) = trainClassifier_BaggingTrees_LOOCV(dcTestAvgHbO,mlActAuto{1},numChn,numClasses,timePt,indexMoviesTest,param,trNum);
    performanceArr_Boosting_HbO(i1,:) = trainClassifier_BoostingTrees_LOOCV(dcTestAvgHbO,mlActAuto{1},numChn,numClasses,timePt,indexMoviesTest,param,trNum);
    
    performanceArr_LDA_LW_HbR(i1,:) = train_RLDA_Ledoit_LOOCV(dcTestAvgHbR,mlActAuto{1},numChn,numClasses,timePt,indexMoviesTest,param,trNum);
    performanceArr_LDA_CERNN_HbR(i1,:) = train_RLDA_CERNN_LOOCV(dcTestAvgHbR,mlActAuto{1},numChn,numClasses,timePt,indexMoviesTest,param,trNum);
    performanceArr_Logistics_HbR(i1,:) = train_LogisticRegression_LOOCV(dcTestAvgHbR,mlActAuto{1},numChn,numClasses,timePt,indexMoviesTest,param,trNum);
    performanceArr_SVM_Lin_HbR(i1,:) = trainClassifier_LinearSVM_LOOCV(dcTestAvgHbR,mlActAuto{1},numChn,numClasses,timePt,indexMoviesTest,param,trNum);
    performanceArr_Bagging_HbR(i1,:) = trainClassifier_BaggingTrees_LOOCV(dcTestAvgHbR,mlActAuto{1},numChn,numClasses,timePt,indexMoviesTest,param,trNum);
    performanceArr_Boosting_HbR(i1,:) = trainClassifier_BoostingTrees_LOOCV(dcTestAvgHbR,mlActAuto{1},numChn,numClasses,timePt,indexMoviesTest,param,trNum);
    
    performanceArr_LDA_LW_HbT(i1,:) = train_RLDA_Ledoit_LOOCV(dcTestAvgHbT,mlActAuto{1},numChn,numClasses,timePt,indexMoviesTest,param,trNum);
    performanceArr_LDA_CERNN_HbT(i1,:) = train_RLDA_CERNN_LOOCV(dcTestAvgHbT,mlActAuto{1},numChn,numClasses,timePt,indexMoviesTest,param,trNum);
    performanceArr_Logistics_HbT(i1,:) = train_LogisticRegression_LOOCV(dcTestAvgHbT,mlActAuto{1},numChn,numClasses,timePt,indexMoviesTest,param,trNum);
    performanceArr_SVM_Lin_HbT(i1,:) = trainClassifier_LinearSVM_LOOCV(dcTestAvgHbT,mlActAuto{1},numChn,numClasses,timePt,indexMoviesTest,param,trNum);
    performanceArr_Bagging_HbT(i1,:) = trainClassifier_BaggingTrees_LOOCV(dcTestAvgHbT,mlActAuto{1},numChn,numClasses,timePt,indexMoviesTest,param,trNum);
    performanceArr_Boosting_HbT(i1,:) = trainClassifier_BoostingTrees_LOOCV(dcTestAvgHbT,mlActAuto{1},numChn,numClasses,timePt,indexMoviesTest,param,trNum);
    
%     snirf2.Save([saveDir filesep rawDataFN{1} 'Combined.snirf']);
%     fileName = [processedDataDir filesep 'intermediateOutputsCombined.mat'];
%     save(fileName,'dcTr','dcAvgTr','dcAvgStd','nTrials','dcNew',...
%         'dcResid','dcSum2','beta','R','hmrstats','dod','stim','tRange',...
%         'tIncAuto','mlActAuto','tIncAuto','snirf2');


end

fileName = [processedDataDir filesep 'GLM_CV_Performances.mat'];
save(fileName,'performanceArr_LDA_LW_HbO','performanceArr_LDA_CERNN_HbO',...
    'performanceArr_Logistics_HbO','performanceArr_SVM_Lin_HbO',...
    'performanceArr_Bagging_HbO','performanceArr_Boosting_HbO',...
    'performanceArr_LDA_LW_HbR','performanceArr_LDA_CERNN_HbR',...
    'performanceArr_Logistics_HbR','performanceArr_SVM_Lin_HbR',...
    'performanceArr_Bagging_HbR','performanceArr_Boosting_HbR',...
    'performanceArr_LDA_LW_HbT','performanceArr_LDA_CERNN_HbT',...
    'performanceArr_Logistics_HbT','performanceArr_SVM_Lin_HbT',...
    'performanceArr_Bagging_HbT','performanceArr_Boosting_HbT');

avgPerformance_LDA_LW_HbO = squeeze(mean(performanceArr_LDA_LW_HbO,1));
avgPerformance_LDA_CERNN_HbO = squeeze(mean(performanceArr_LDA_CERNN_HbO,1));
avgPerformance_Logistics_HbO = squeeze(mean(performanceArr_Logistics_HbO,1));
avgPerformance_SVM_Lin_HbO = squeeze(mean(performanceArr_SVM_Lin_HbO,1));
avgPerformance_Bagging_HbO = squeeze(mean(performanceArr_Bagging_HbO,1));
avgPerformance_Boosting_HbO = squeeze(mean(performanceArr_Boosting_HbO,1));

avgPerformance_LDA_LW_HbR = squeeze(mean(performanceArr_LDA_LW_HbR,1));
avgPerformance_LDA_CERNN_HbR = squeeze(mean(performanceArr_LDA_CERNN_HbR,1));
avgPerformance_Logistics_HbR = squeeze(mean(performanceArr_Logistics_HbR,1));
avgPerformance_SVM_Lin_HbR = squeeze(mean(performanceArr_SVM_Lin_HbR,1));
avgPerformance_Bagging_HbR = squeeze(mean(performanceArr_Bagging_HbR,1));
avgPerformance_Boosting_HbR = squeeze(mean(performanceArr_Boosting_HbR,1));

avgPerformance_LDA_LW_HbT = squeeze(mean(performanceArr_LDA_LW_HbT,1));
avgPerformance_LDA_CERNN_HbT = squeeze(mean(performanceArr_LDA_CERNN_HbT,1));
avgPerformance_Logistics_HbT = squeeze(mean(performanceArr_Logistics_HbT,1));
avgPerformance_SVM_Lin_HbT = squeeze(mean(performanceArr_SVM_Lin_HbT,1));
avgPerformance_Bagging_HbT = squeeze(mean(performanceArr_Bagging_HbT,1));
avgPerformance_Boosting_HbT = squeeze(mean(performanceArr_Boosting_HbT,1));

numClassifiers = 6;
cmap = jet(numClassifiers);

figure('units','normalized','outerposition',[0 0 1 1]);hold on;
subplot(1,3,1);
plot(timePt./fs-2,avgPerformance_LDA_LW_HbO,'Color',cmap(1,:));
plot(timePt./fs-2,avgPerformance_LDA_CERNN_HbO,'Color',cmap(2,:));
plot(timePt./fs-2,avgPerformance_Logistics_HbO,'Color',cmap(3,:));
plot(timePt./fs-2,avgPerformance_SVM_Lin_HbO,'Color',cmap(4,:));
plot(timePt./fs-2,avgPerformance_Bagging_HbO,'Color',cmap(5,:));
plot(timePt./fs-2,avgPerformance_Boosting_HbO,'Color',cmap(6,:));

ylim([0 1]);
title(sprintf('Sbj %s: Δ[HbO] Multi',num2str(sbjNum)));
legend({'LDAShrink','LDA CERNN','LogReg','SVM','Bagging','Boosted'});
xlabel('Time [s]');ylabel('Accuracy');
hold off;

subplot(1,3,2);
plot(timePt./fs-2,avgPerformance_LDA_LW_HbR,'Color',cmap(1,:));
plot(timePt./fs-2,avgPerformance_LDA_CERNN_HbR,'Color',cmap(2,:));
plot(timePt./fs-2,avgPerformance_Logistics_HbR,'Color',cmap(3,:));
plot(timePt./fs-2,avgPerformance_SVM_Lin_HbR,'Color',cmap(4,:));
plot(timePt./fs-2,avgPerformance_Bagging_HbR,'Color',cmap(5,:));
plot(timePt./fs-2,avgPerformance_Boosting_HbR,'Color',cmap(6,:));

ylim([0 1]);
title(sprintf('Sbj %s: Δ[HbR] Multi',num2str(sbjNum)));
legend({'LDAShrink','LDA CERNN','LogReg','SVM','Bagging','Boosted'});
xlabel('Time [s]');ylabel('Accuracy');
hold off;

subplot(1,3,3);
plot(timePt./fs-2,avgPerformance_LDA_LW_HbT,'Color',cmap(1,:));
plot(timePt./fs-2,avgPerformance_LDA_CERNN_HbT,'Color',cmap(2,:));
plot(timePt./fs-2,avgPerformance_Logistics_HbT,'Color',cmap(3,:));
plot(timePt./fs-2,avgPerformance_SVM_Lin_HbT,'Color',cmap(4,:));
plot(timePt./fs-2,avgPerformance_Bagging_HbT,'Color',cmap(5,:));
plot(timePt./fs-2,avgPerformance_Boosting_HbT,'Color',cmap(6,:));

ylim([0 1]);
title(sprintf('Sbj %s: Δ[HbT] Multi',num2str(sbjNum)));
legend({'LDAShrink','LDA CERNN','LogReg','SVM','Bagging','Boosted'});
xlabel('Time [s]');ylabel('Accuracy');
hold off;

if numClasses == 2
    fn = sprintf('PerformanceCumsumVsTime_GLM_CV_LR_AllChns_DiffClassifiers');
else
    fn = sprintf('PerformanceCumsumVsTime_GLM_CV_AllChns_DiffClassifiers');
end

%if saveOp == 1
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
%end

end
