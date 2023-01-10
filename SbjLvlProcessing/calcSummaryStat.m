% calculate R correlation coefficient with GLM fit HRF, and avg R, number
% of channels rejected, behavioral score, p and t-value score for all pairs
% For fNIRS
function calcSummaryStat(sbjNum,snrThresh)
rawDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
procDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
fn = ['summaryStat_' num2str(snrThresh) '.mat'];

load([rawDir filesep sbjNum '.mat'],'s');

% R is from intermediateOutputsCombined_Basis1.mat file
%load([procDir filesep 'intermediateOutputsCombined_Basis1.mat'],'R','dc','stim');

% Rerun preprocess step.
%sbjNum = s.name;
rawDataFN = s.fName;
movieList = s.movieList;
behData = s.resp;
startT = s.startT;
endT = s.endT;

fs = 50;

dataTimeSeries = [];
timeL = 0;
%allS = zeros(720,1);
allS = zeros(1080,1);
startI = 1;

for i = 1:length(rawDataFN)
    
    snirfTemp = SnirfClass(load( [rawDataFN{i} '.nirs'],'-mat'));
    allSTemp = find(snirfTemp.aux(1,1).dataTimeSeries>1);
    aIndTemp = find(allSTemp(2:end)-allSTemp(1:end-1)==1);
    allSTemp(aIndTemp) = [];
    stimLen = size(allSTemp,1);
    endI = startI + stimLen-1;
    allS(startI:endI) = allSTemp;
    startI = endI + 1;
    
    dataTimeSeries = [dataTimeSeries; snirfTemp.data.dataTimeSeries];
    timeL = timeL + size(snirfTemp.data.time,1);
    
end

ml = snirfTemp.data.measurementList;

if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
    cueOnsetIndex = 1:4:720;
else
    cueOnsetIndex = 1:4:360;
end

allS = allS(cueOnsetIndex);

if startT ~= 1
    % allS is in sec
    allS = allS-startT;
end

load([rawDir filesep movieList '.mat'],'indexMoviesTest');
load([rawDir filesep behData '.mat'],'responsesA','responsesV','correctRespA','correctRespV');

%indexMoviesTest = updateMovieList(allS,indexMoviesTest);

% Define stimclass here
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

% dAll = DataClass([snirf1.data.dataTimeSeries; snirf3.data.dataTimeSeries; snirf4.data.dataTimeSeries],...
%     0:1/50:(size(snirf1.data.dataTimeSeries,1)+size(snirf3.data.dataTimeSeries,1)+size(snirf4.data.dataTimeSeries,1)-1)/50,snirf1.data.measurementList);

dAll = DataClass(dataTimeSeries,0:1/50:((timeL-1)/50),ml);

%     dAll = DataClass(snirf1.data.dataTimeSeries,...
%         0:1/50:(size(snirf1.data.dataTimeSeries,1)-1)/50,snirf1.data.measurementList);

if endT ~= -1
    dAll.dataTimeSeries = dAll.dataTimeSeries(startT*fs:endT*fs,:);
    dAll.time = dAll.time(startT*fs:endT*fs);
else
    dAll.dataTimeSeries = dAll.dataTimeSeries(startT*fs:end,:);
    dAll.time = dAll.time(startT*fs:end);
end

data = dAll;
probe = snirfTemp.probe;
mlActMan = {};
tIncMan = {};
Aaux = [];
rcMap = [];

%mlActAuto = hmrR_PruneChannels(data,probe,mlActMan,tIncMan,[10000  10000000],2,[0  45]);
    %mlActAuto = hmrR_PruneChannels(data,probe,mlActMan,tIncMan,[10000  10000000],2,[0  45]);
mlActAuto = hmrR_PruneChannels(data,probe,mlActMan,tIncMan,[7000  10000000],snrThresh,[0  45]);

dod = hmrR_Intensity2OD(data);

% For sbj 12, params in 1st line for motion correct/artifact perform
% significantly better than params in 2nd line.
[dod] = hmrR_MotionCorrectPCArecurse(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,15,4,0.97,5,1);
%[dod,tInc,svs,nSV,tInc0] = hmrR_MotionCorrectPCArecurse(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,20,5,0.97,5,1);

tIncAuto = hmrR_MotionArtifact(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,15,5);

[stimAll,~] = hmrR_StimRejection(dod,stimAll,tIncAuto,tIncMan,[-2  15]);

ml690 = mlActAuto{1}(1:length(mlActAuto{1})/2);
ml870 = mlActAuto{1}(length(mlActAuto{1})/2+1:length(mlActAuto{1}));

numUnusedChns = length(mlActAuto{1})/2-sum(ml690&ml870);

if strcmp(sbjNum,'04')||strcmp(sbjNum,'07') || strcmp(sbjNum,'08') || strcmp(sbjNum,'10')
    isNew = 0;
else
    isNew = 1;
end

chnList = getChnList(isNew);

rejChns = chnList(~logical(ml690&ml870));

% Avg time of trial
timeTr = diff(allS);
avgTimeTr = mean(timeTr)/fs;

% Number of trials rejected
numRejTr = 0;
for i = 1:size(stimAll,2)
    numRejTr = numRejTr + sum(find(stimAll(i).states(:,2)<0));
end

%numRejTr = sum(find(stimAll.states(:,2)<0));

behScore = sum((responsesA==correctRespA)&(responsesV==correctRespV))/length(correctRespV);

save([procDir filesep fn],'mlActAuto','behScore',...
    'numUnusedChns','rejChns','avgTimeTr','numRejTr');

end