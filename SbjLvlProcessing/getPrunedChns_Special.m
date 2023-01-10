% for sbj 08 and 10
function [mlActAuto_Unused, allCues] = getPrunedChns_Special(s,snrThresh)

sbjNum = s.name;
rawDataFN = s.fName;
movieList = s.movieList;
%behData = s.resp;
startT = s.startT;
endT = s.endT;

fs = 50;

saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
% Convert nirs to snirf file format
snirf1 = SnirfClass(load([rawDataFN{1} '.nirs'],'-mat'));
snirf3 = SnirfClass(load([rawDataFN{2} '.nirs'],'-mat'));
snirf4 = SnirfClass(load([rawDataFN{3} '.nirs'],'-mat'));
snirf1.Info()

% Extract aux and convert to stimclass
allS1 = find(snirf1.aux(1,1).dataTimeSeries>1);
aInd1 = find(allS1(2:end)-allS1(1:end-1)==1);
allS3 = find(snirf3.aux(1,1).dataTimeSeries>1);
aInd3 = find(allS3(2:end)-allS3(1:end-1)==1);
allS4 = find(snirf4.aux(1,1).dataTimeSeries>1);
aInd4 = find(allS4(2:end)-allS4(1:end-1)==1);
allS1(aInd1) = [];
allS1 = allS1./50;
allS3(aInd3) = [];
allS3 = allS3./50;
allS4(aInd4) = [];
allS4 = allS4./50;

oldStimClass1Length1 = size(allS1,1);
oldStimClass1Length3 = oldStimClass1Length1 + size(allS3,1);
oldStimClass1Length4 = oldStimClass1Length3 + size(allS4,1);

if strcmp(sbjNum,'04')||strcmp(sbjNum,'07')
    cueOnsetIndex = 1:4:1080;
else
    cueOnsetIndex = 1:4:720;
end
%cueOnsetIndex = 1:4:360;
cueOnsetIndex1 = 1:oldStimClass1Length1;
cueOnsetIndex3 = oldStimClass1Length1+1:oldStimClass1Length3;
cueOnsetIndex4 = oldStimClass1Length3+1:oldStimClass1Length4;

trialNum1 = sum(ismember(cueOnsetIndex1,cueOnsetIndex));
trialNum3 = sum(ismember(cueOnsetIndex3,cueOnsetIndex));
trialNum4 = sum(ismember(cueOnsetIndex4,cueOnsetIndex));

lastTrial1 = trialNum1;
lastTrial3 = lastTrial1 + trialNum3;
lastTrial4 = lastTrial3 + trialNum4;

snirf1TLen = size(snirf1.data.time,1)/(50);
snirf3TLen = size(snirf3.data.time,1)/(50);
allS3 = allS3+snirf1TLen;
allS4 = allS4+snirf1TLen+snirf3TLen;

allS = [allS1; allS3; allS4];

if startT ~= 1
    % allS is in sec
    allS = allS-startT;
end

allCues = allS(cueOnsetIndex);
%allS = allS1;

% split triggers into 6 categories
load([saveDir filesep movieList '.mat'],'fixedMaskerList','indexMoviesTest','maskerMovies','numTrials','uniqueMovies');

% snirf 2
% leftOnsetS = StimClass('leftSingle');
% rightOnsetS = StimClass('rightSingle');
% centerOnsetS = StimClass('centerSingle');
leftOnsetM = StimClass('leftMulti');
rightOnsetM = StimClass('rightMulti');
centerOnsetM = StimClass('centerMulti');

% leftSingleIndex = indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==0;
% rightSingleIndex = indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==0;
% centerSingleIndex = indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==0;
leftMultiIndex = indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1;
rightMultiIndex = indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1;
centerMultiIndex = indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==1;

% leftSingleIndex = leftSingleIndex(leftSingleIndex<= lastTrial2);
% rightSingleIndex = rightSingleIndex(rightSingleIndex<= lastTrial2);
% centerSingleIndex = centerSingleIndex(centerSingleIndex<= lastTrial2);
% leftMultiIndex = leftMultiIndex(leftMultiIndex<= lastTrial2);
% rightMultiIndex = rightMultiIndex(rightMultiIndex<= lastTrial2);
% centerMultiIndex = centerMultiIndex(centerMultiIndex<= lastTrial2);

% AddStims(leftOnsetS, allS(cueOnsetIndex(leftSingleIndex)));
% AddStims(rightOnsetS, allS(cueOnsetIndex(rightSingleIndex)));
% AddStims(centerOnsetS, allS(cueOnsetIndex(centerSingleIndex)));
AddStims(leftOnsetM, allS(cueOnsetIndex(leftMultiIndex)));
AddStims(rightOnsetM, allS(cueOnsetIndex(rightMultiIndex)));
AddStims(centerOnsetM, allS(cueOnsetIndex(centerMultiIndex)));

% updateStates(leftOnsetS);
% updateStates(rightOnsetS);
% updateStates(centerOnsetS);
updateStates(leftOnsetM);
updateStates(rightOnsetM);
updateStates(centerOnsetM);

% snirf2.stim(1,1) = leftOnsetS;
% snirf2.stim(1,2) = rightOnsetS;
% snirf2.stim(1,3) = centerOnsetS;
snirf1.stim(1,1) = leftOnsetM;
snirf1.stim(1,2) = rightOnsetM;
snirf1.stim(1,3) = centerOnsetM;

% combine snirf2.data.dataTimeSeries and time
snirf1.data.dataTimeSeries = [snirf1.data.dataTimeSeries; snirf3.data.dataTimeSeries; snirf4.data.dataTimeSeries];
snirf1.data.time = 0:1/50:(size(snirf1.data.dataTimeSeries,1)-1)/50;

dAll = DataClass(snirf1.data.dataTimeSeries,...
    0:1/50:(size(snirf1.data.dataTimeSeries,1)-1)/50,snirf1.data.measurementList);

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
stim = snirf1.stim;
%tIncAuto = snirf1.tIncAuto;
%dc = snirf1.dc;
Aaux = [];
rcMap = [];

mlActAuto_Unused = hmrR_PruneChannels(data,probe,mlActMan,tIncMan,[7000  10000000],snrThresh,[0  45]);

end