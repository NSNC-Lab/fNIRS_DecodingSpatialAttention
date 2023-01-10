% this starts with nirs files so can run again to recreate snirf files.
% this combine nirs files so we can compute HRF from whole time series.
% Skip pruning step
% Different parameters

function preprocessFNIRS06_BVar_MultiOnly(sbjNum,rawDataFN,movieList)
saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
% Convert nirs to snirf file format
snirf1 = SnirfClass(load([rawDataFN{1} '.nirs'],'-mat'));
%snirf3 = SnirfClass(load([rawDataFN{2} '.nirs'],'-mat'));
%snirf4 = SnirfClass(load([rawDataFN{3} '.nirs'],'-mat'));
snirf1.Info()

% Extract aux and convert to stimclass
allS1 = find(snirf1.aux(1,1).dataTimeSeries>1);
aInd1 = find(allS1(2:end)-allS1(1:end-1)==1);
% allS3 = find(snirf3.aux(1,1).dataTimeSeries>1);
% aInd3 = find(allS3(2:end)-allS3(1:end-1)==1);
% allS4 = find(snirf4.aux(1,1).dataTimeSeries>1);
% aInd4 = find(allS4(2:end)-allS4(1:end-1)==1);
allS1(aInd1) = [];
allS1 = allS1./50;
% allS3(aInd3) = [];
% allS3 = allS3./50;
% allS4(aInd4) = [];
% allS4 = allS4./50;

oldStimClass1Length1 = size(allS1,1);
% oldStimClass1Length3 = oldStimClass1Length1 + size(allS3,1);
% oldStimClass1Length4 = oldStimClass1Length3 + size(allS4,1);

cueOnsetIndex = 1:4:720;
%cueOnsetIndex = 1:4:360;
cueOnsetIndex1 = 1:oldStimClass1Length1;
% cueOnsetIndex3 = oldStimClass1Length1+1:oldStimClass1Length3;
% cueOnsetIndex4 = oldStimClass1Length3+1:oldStimClass1Length4;

trialNum1 = sum(ismember(cueOnsetIndex1,cueOnsetIndex));
% trialNum3 = sum(ismember(cueOnsetIndex3,cueOnsetIndex));
% trialNum4 = sum(ismember(cueOnsetIndex4,cueOnsetIndex));

lastTrial1 = trialNum1;
% lastTrial3 = lastTrial1 + trialNum3;
% lastTrial4 = lastTrial3 + trialNum4;

snirf1TLen = size(snirf1.data.time,1)/(50);
% snirf3TLen = size(snirf3.data.time,1)/(50);
% allS3 = allS3+snirf1TLen;
% allS4 = allS4+snirf1TLen+snirf3TLen;

%allS = [allS1; allS3; allS4];
allS = allS1;

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
%snirf2.data.dataTimeSeries = [snirf2.data.dataTimeSeries; snirf3.data.dataTimeSeries; snirf4.data.dataTimeSeries];
%snirf2.data.time = [snirf2.data.time;snirf3.data.time;snirf4.data.time];
%snirf2.data.time = 0:1/50:(size(snirf2.data.dataTimeSeries,1)-1)/50;

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

%mlActAuto_Unused = hmrR_PruneChannels(data,probe,mlActMan,tIncMan,[10000  10000000],2,[0  45]);

mlActAuto = {ones(size(snirf1.data.measurementList,2),1)};
dod = hmrR_Intensity2OD(data);

% For sbj 12, params in 1st line for motion correct/artifact perform
% significantly better than params in 2nd line.
[dod,tInc,svs,nSV,tInc0] = hmrR_MotionCorrectPCArecurse(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,15,4,0.97,5,1);
%[dod,tInc,svs,nSV,tInc0] = hmrR_MotionCorrectPCArecurse(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,20,5,0.97,5,1);

tIncAuto = hmrR_MotionArtifact(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,15,5);
%tIncAuto = hmrR_MotionArtifact(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,20,4);

[stim,tRange] = hmrR_StimRejection(dod,stim,tIncAuto,tIncMan,[-2  15]);

dod = hmrR_BandpassFilt(dod,0.01,0.5);
%dod = hmrR_BandpassFilt(dod,0,0.5);

dc = hmrR_OD2Conc(dod,probe,[1  1  1]);

[dcAvg,dcAvgStd,nTrials,dcNew,dcResid,dcSum2,beta,R,hmrstats,bvar] = hmrR_GLM_bvar(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
    [-2  15],1,1,[1.0 1.0 0.0 0.0 0.0 0.0],15,1,0,0);

if ~exist(processedDataDir,'dir')
    mkdir(processedDataDir);
end

snirf1.Save([saveDir filesep rawDataFN{1} 'CombinedBasis1.snirf']);
fileName = [processedDataDir filesep 'intermediateOutputsCombined_Basis1.mat'];
save(fileName,'dc','dcAvg','dcAvgStd','nTrials','dcNew',...
    'dcResid','dcSum2','beta','R','bvar','hmrstats','dod','stim','tRange',...
    'tIncAuto','mlActAuto','tIncAuto','snirf1','allS');

% [dcAvg,dcAvgStd,nTrials,dcNew,dcResid,dcSum2,beta,R,hmrstats] = hmrR_GLM(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
%     [-2  15],1,2,[0.1 3.0 10.0 1.8 3.0 10.0],15,1,3,0);
% 
% snirf1.Save([saveDir filesep rawDataFN{1} 'CombinedBasis2.snirf']);
% fileName = [processedDataDir filesep 'intermediateOutputsCombined_Basis2.mat'];
% save(fileName,'dc','dcAvg','dcAvgStd','nTrials','dcNew',...
%     'dcResid','dcSum2','beta','R','hmrstats','dod','stim','tRange',...
%     'tIncAuto','mlActAuto','tIncAuto','snirf1','allS');
% 
% [dcAvg,dcAvgStd,nTrials,dcNew,dcResid,dcSum2,beta,R,hmrstats] = hmrR_GLM(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
%     [-2  15],1,3,[0.1 3.0 10.0 1.8 3.0 10.0],15,1,3,0);
% 
% snirf1.Save([saveDir filesep rawDataFN{1} 'CombinedBasis3.snirf']);
% fileName = [processedDataDir filesep 'intermediateOutputsCombined_Basis3.mat'];
% save(fileName,'dc','dcAvg','dcAvgStd','nTrials','dcNew',...
%     'dcResid','dcSum2','beta','R','hmrstats','dod','stim','tRange',...
%     'tIncAuto','mlActAuto','tIncAuto','snirf1','allS');
% 
% [dcAvg,dcAvgStd,nTrials,dcNew,dcResid,dcSum2,beta,R,hmrstats] = hmrR_GLM(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
%     [-2  15],1,4,[8.6 0.547 10.0 8.6 0.547 10.0],15,1,3,0);
% 
% snirf1.Save([saveDir filesep rawDataFN{1} 'CombinedBasis4.snirf']);
% fileName = [processedDataDir filesep 'intermediateOutputsCombined_Basis4.mat'];
% save(fileName,'dc','dcAvg','dcAvgStd','nTrials','dcNew',...
%     'dcResid','dcSum2','beta','R','hmrstats','dod','stim','tRange',...
%     'tIncAuto','mlActAuto','tIncAuto','snirf1','allS');

end