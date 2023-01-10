% this starts with nirs files so can run again to recreate snirf files.
% this combine nirs files so we can compute HRF from whole time series.
% All 4 basis. Skip pruning step.
% For sbj 03. For sbj 08 and 10, use preprocessFNIRS06.m. For 12 and after,
% use preprocessFNIRS05_MultiOnly

function preprocessFNIRS06_Sbj03(sbjNum,rawDataFN,respData)
saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
% Convert nirs to snirf file format
snirf1 = SnirfClass(load([rawDataFN{1} '.nirs'],'-mat'));
snirf1.Info()

% Extract aux and convert to stimclass
allS2 = find(snirf1.aux(1,1).dataTimeSeries>1);
aInd2 = find(allS2(2:end)-allS2(1:end-1)==1);
allS2(aInd2) = [];
allS2 = allS2./50;

oldStimClass1Length2 = size(allS2,1);

if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
    cueOnsetIndex = 1:4:720;
else
    cueOnsetIndex = 1:4:1080;
end
cueOnsetIndex2 = 1:oldStimClass1Length2;

trialNum2 = sum(ismember(cueOnsetIndex2,cueOnsetIndex));

lastTrial2 = trialNum2;

snirf2TLen = size(snirf1.data.time,1)/(50);

allS = [allS2];

idxCutoff = find(cueOnsetIndex>length(allS));
cueOnsetIndex(idxCutoff) = [];

% split triggers into 6 categories
load([saveDir filesep respData '.mat'],'fixedMaskerList','indexMoviesTest','maskerMovies','numTrials','uniqueMovies');

% snirf 2
leftOnsetS = StimClass('leftSingle');
rightOnsetS = StimClass('rightSingle');
centerOnsetS = StimClass('centerSingle');
leftOnsetM = StimClass('leftMulti');
rightOnsetM = StimClass('rightMulti');
centerOnsetM = StimClass('centerMulti');

indexMoviesTest = indexMoviesTest(1:trialNum2,:);

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
%snirf1.data.dataTimeSeries = [snirf1.data.dataTimeSeries; snirf3.data.dataTimeSeries; snirf4.data.dataTimeSeries];
%snirf2.data.time = [snirf2.data.time;snirf3.data.time;snirf4.data.time];
%snirf1.data.time = 0:1/50:(size(snirf1.data.dataTimeSeries,1)-1)/50;

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

tIncAuto = hmrR_MotionArtifact(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,15,5);
%tIncAuto = hmrR_MotionArtifact(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,20,4);

[stim,tRange] = hmrR_StimRejection(dod,stim,tIncAuto,tIncMan,[-2  15]);

dod = hmrR_BandpassFilt(dod,0.01,0.5);
%dod = hmrR_BandpassFilt(dod,0,0.5);

dc = hmrR_OD2Conc(dod,probe,[1  1  1]);

[dcAvg,dcAvgStd,nTrials,dcNew,dcResid,dcSum2,beta,R,hmrstats] = hmrR_GLM(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
    [-2  15],1,1,[1.0 1.0 0.0 0.0 0.0 0.0],15,1,0,0);

snirf1.Save([saveDir filesep rawDataFN{1} 'Combined_Basis1.snirf']);
fileName = [processedDataDir filesep 'intermediateOutputsCombined_Basis1.mat'];
save(fileName,'dc','dcAvg','dcAvgStd','nTrials','dcNew',...
    'dcResid','dcSum2','beta','R','hmrstats','dod','stim','tRange',...
    'tIncAuto','mlActAuto','tIncAuto','snirf1','allS');

% #2
[dcAvg,dcAvgStd,nTrials,dcNew,dcResid,dcSum2,beta,R,hmrstats] = hmrR_GLM(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
    [-2  15],1,2,[0.1 3.0 10.0 1.8 3.0 10.0],15,1,3,0);

snirf1.Save([saveDir filesep rawDataFN{1} 'Combined_Basis2.snirf']);
fileName = [processedDataDir filesep 'intermediateOutputsCombined_Basis2.mat'];
save(fileName,'dc','dcAvg','dcAvgStd','nTrials','dcNew',...
    'dcResid','dcSum2','beta','R','hmrstats','dod','stim','tRange',...
    'tIncAuto','mlActAuto','tIncAuto','snirf1','allS');

% #3
[dcAvg,dcAvgStd,nTrials,dcNew,dcResid,dcSum2,beta,R,hmrstats] = hmrR_GLM(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
    [-2  15],1,3,[0.1 3.0 10.0 1.8 3.0 10.0],15,1,3,0);

snirf1.Save([saveDir filesep rawDataFN{1} 'Combined_Basis3.snirf']);
fileName = [processedDataDir filesep 'intermediateOutputsCombined_Basis3.mat'];
save(fileName,'dc','dcAvg','dcAvgStd','nTrials','dcNew',...
    'dcResid','dcSum2','beta','R','hmrstats','dod','stim','tRange',...
    'tIncAuto','mlActAuto','tIncAuto','snirf1','allS');

% #4. For experiment 08 #4, when saving, file may be corrupt. error closing
% file.
[dcAvg,dcAvgStd,nTrials,dcNew,dcResid,dcSum2,beta,R,hmrstats] = hmrR_GLM(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
    [-2  15],1,4,[8.6 0.547 10.0 8.6 0.547 10.0],15,1,3,0);

% Full save
% snirf1.Save([saveDir filesep rawDataFN{1} 'Combined_Basis4.snirf']);
% fileName = [processedDataDir filesep 'intermediateOutputsCombined_Basis4.mat'];
% save(fileName,'dc','dcAvg','dcAvgStd','nTrials','dcNew',...
%     'dcResid','dcSum2','beta','R','hmrstats','dod','stim','tRange',...
%     'tIncAuto','mlActAuto','tIncAuto','snirf1','allS');

% Mini save, only for Experiment 08
snirf1.Save([saveDir filesep rawDataFN{1} 'Combined_Basis4.snirf']);
fileName = [processedDataDir filesep 'intermediateOutputsCombined_Basis4.mat'];
save(fileName,'dc','dcAvg','dcNew',...
    'beta','R','hmrstats','dod','stim','tRange',...
    'mlActAuto','snirf1','allS');


% 
% fileName = [processedDataDir filesep 'intermediateOutputsCombined_Basis4_Part2.mat'];
% save(fileName,'stim','tRange',...
%     'tIncAuto','mlActAuto','tIncAuto','snirf1','allS');

end