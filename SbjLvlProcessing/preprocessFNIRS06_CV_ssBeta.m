% STATUS: experimental.
% Not for publication.
% 
% SYNTAX:
% preprocessFNIRS06_CV_ssBeta(sbjNum,rawDataFN,respData)
% 
% DESCRIPTION:
% Set SS Beta a priori to 1. No GLM fitting here.
% 
% RESTRICTION:
% Only for sbj 08 and 10.
% 
% INPUTS:
% sbjNum - subject ID. 2-digits string. Ex: '08'.
% rawDataFN - 1xN cell array of nirs data recordings. Use struct s.fName as
%   input.
% respData - file name of mat file containing different variables related
% to experiment:
%   fixedMaskerList - cell array: file names of fixed masker movies
%   indexMoviesTest - numTrials x 5 double array:
%       col 1: index of target movies in uniqueMovies
%       col 2: index of spatial location
%       col 3: boolean: masker is fixed or random
%       col 4: index of masker movies in fixedMaskerList(?)
%       col 5: boolean: condition is target-alone or target+maskers
%   maskerMovies - cell array: file names of random masker movies
%   numTrials - int: number of trials the subject took
%   uniqueMovies - cell array: file names of target movies
%
% RETURNED VARIABLES:
% None.
% 
% FILES SAVED:
% 1) save snirf as 'CombinedBasis1_SSBetaPrior.snirf'
% 2) save selected variables as
%   'intermediateOutputsCombined_Basis1_SSBetaPrior.mat'
% 
% PLOTTING:
% None.

function preprocessFNIRS06_CV_ssBeta(sbjNum,rawDataFN,respData)
saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
% Convert nirs to snirf file format
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

oldStimClass1Length2 = size(allS2,1);
oldStimClass1Length3 = oldStimClass1Length2 + size(allS3,1);
oldStimClass1Length4 = oldStimClass1Length3 + size(allS4,1);

if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
    cueOnsetIndex = 1:4:720;
else
    % For sbj 07 and before.
    cueOnsetIndex = 1:4:1080;
end
cueOnsetIndex2 = 1:oldStimClass1Length2;
cueOnsetIndex3 = oldStimClass1Length2+1:oldStimClass1Length3;
cueOnsetIndex4 = oldStimClass1Length3+1:oldStimClass1Length4;

trialNum2 = sum(ismember(cueOnsetIndex2,cueOnsetIndex));
trialNum3 = sum(ismember(cueOnsetIndex3,cueOnsetIndex));
trialNum4 = sum(ismember(cueOnsetIndex4,cueOnsetIndex));

lastTrial2 = trialNum2;
lastTrial3 = lastTrial2 + trialNum3;
lastTrial4 = lastTrial3 + trialNum4;

snirf2TLen = size(snirf1.data.time,1)/(50);
snirf3TLen = size(snirf3.data.time,1)/(50);
allS3 = allS3+snirf2TLen;
allS4 = allS4+snirf2TLen+snirf3TLen;

allS = [allS2; allS3; allS4];

% split triggers into 6 categories
load([saveDir filesep respData '.mat'],'indexMoviesTest');

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

snirf1.stim(1,1) = leftOnsetS;
snirf1.stim(1,2) = rightOnsetS;
snirf1.stim(1,3) = centerOnsetS;
snirf1.stim(1,4) = leftOnsetM;
snirf1.stim(1,5) = rightOnsetM;
snirf1.stim(1,6) = centerOnsetM;

% combine snirf2.data.dataTimeSeries and time
snirf1.data.dataTimeSeries = [snirf1.data.dataTimeSeries; snirf3.data.dataTimeSeries; snirf4.data.dataTimeSeries];
%snirf2.data.time = [snirf2.data.time;snirf3.data.time;snirf4.data.time];
snirf1.data.time = 0:1/50:(size(snirf1.data.dataTimeSeries,1)-1)/50;

data = snirf1.data;
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

%mlActAuto = hmrR_PruneChannels(data,probe,mlActMan,tIncMan,[10000  10000000],2,[0  45]);
mlActAuto = {ones(size(snirf1.data.measurementList,2),1)};
dod = hmrR_Intensity2OD(data);

% For sbj 12, params in 1st line for motion correct/artifact perform
% significantly better than params in 2nd line.
[dod,tInc,svs,nSV,tInc0] = hmrR_MotionCorrectPCArecurse(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,15,4,0.97,5,1);
%[dod,tInc,svs,nSV,tInc0] = hmrR_MotionCorrectPCArecurse(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,20,5,0.97,5,1);

if strcmp(sbjNum,'24')||strcmp(sbjNum,'25')
    tIncAuto = hmrR_MotionArtifact(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,20,5);
else
    tIncAuto = hmrR_MotionArtifact(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,15,5);
end
%tIncAuto = hmrR_MotionArtifact(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,20,4);

[stim,tRange] = hmrR_StimRejection(dod,stim,tIncAuto,tIncMan,[-2  15]);

dod = hmrR_BandpassFilt(dod,0.01,0.5);
%dod = hmrR_BandpassFilt(dod,0,0.5);

dc = hmrR_OD2Conc(dod,probe,[1  1  1]);

dcNew = hmrR_ssBeta_CV(data,dc, probe, mlActAuto, tIncAuto);

% [dcAvg,dcAvgStd,nTrials,dcNew,dcResid,dcSum2,beta,R,hmrstats,bvar,ssBeta] = hmrR_ssBeta_CV(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
%     [-2  15],1,1,[1.0 1.0 0.0 0.0 0.0 0.0],15,1,0,0);

if ~exist(processedDataDir,'dir')
    mkdir(processedDataDir);
end

snirf1.Save([saveDir filesep rawDataFN{1} 'CombinedBasis1_SSBetaPrior.snirf']);
fileName = [processedDataDir filesep 'intermediateOutputsCombined_Basis1_SSBetaPrior.mat'];
save(fileName,'dc','dcNew','dod','stim','tRange',...
    'tIncAuto','mlActAuto','tIncAuto','snirf1','allS');

end