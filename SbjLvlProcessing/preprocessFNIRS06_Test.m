% 4/11/2022: For testing GLM. Multi only.
% this starts with nirs files so can run again to recreate snirf files.
% this combine nirs files so we can compute HRF from whole time series.
% Skip pruning step
% Different parameters

function preprocessFNIRS06_Test(sbjNum,rawDataFN,movieList)
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
% photon flux: # of photons/sec/area
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

%[data_yavg, data_yavgstd, nTrials, data_ynew, data_yresid, data_ysum2, ...
%    beta_blks, yR_blks, bvar, hmrstats, betaSS_blks, dof_blks, sse_blks] = ...
%    hmrR_GLM_HbT(data_y, stim, probe, mlActAuto, Aaux, tIncAuto, rcMap, trange, glmSolveMethod, idxBasis, paramsBasis, rhoSD_ssThresh, flagNuisanceRMethod, driftOrder, c_vector, betaList)

% betaTest = 3:8;
% betaList = 1:16;
% betaList(betaTest) = [];
betaList = 3:8;
numCond = size(stim,2);
condList = [1 1 1];
betaT.betaList = betaList;
betaT.condList = [1 0 0;0 1 0;0 0 1;1 -1 0;1 0 -1;0 1 -1];

% full model
% check beta_label var!
% complex number when idxBasis == 4. Add abs
% preallocate arrays for faster processing
% add R^2 and partial R^2 and partial F-statistics
% p val for t-statistics now computed with correct dof
[dcAvg,dcAvgStd,~,dcNew,~,~,beta,~,hmrstats,bvar,betaSS] = hmrR_GLM_MN(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
    [-2  15],1,1,[1.0 1.0 0.0 0.0 0.0 0.0],15,1,0,betaT.condList,0);

% [dcAvg,~,~,~,~,~,~,~,~] = hmrR_GLM(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
%         [-2  15],1,4,[8.6 0.547 10.0 8.6 0.547 10.0],15,1,3,0);

% reduced model
% betaTest = 3:8;
% betaList = 1:16;
% betaList(betaTest) = [];
%betaF.betaList = betaList;
% dof_reduced = zeros(1,numCond);
% sse_reduced = zeros(numCond,size(hmrstats.sse_full{1},2),size(hmrstats.sse_full{1},3));
for i = 1:numCond
    tempCondList = condList;
    tempCondList(i) = 0;
%     betaF.condList = tempCondList;
%     [~,~,~,~,~,~,~,~,~,~,~,dof_temp,sse_temp] = hmrR_GLM_HbT_V02(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
%         [-2  15],1,1,[1.0 1.0 0.0 0.0 0.0 0.0],15,1,0,0,tempCondList);

    % Basis #2
    [~,~,~,~,~,~,~,~,hmrstats2,bvar2, betaSS_blks2] = hmrR_GLM_MN(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
        [-2  15],1,2,[0.1 3.0 1.8 3.0 0.1 3.0],15,1,0,0,0,tempCondList);
    
    % Basis #3
    [~,~,~,~,~,~,~,~,hmrstats3,bvar3, betaSS_blks3] = hmrR_GLM_MN(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
        [-2  15],1,3,[0.1 3.0 1.8 3.0 0.1 3.0],15,1,3,0,0,tempCondList);
    
    % Basis #4
    [~,~,~,~,~,~,~,~,hmrstats4,bvar4, betaSS_blks4] = hmrR_GLM_MN(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
        [-2  15],1,4,[8.6 0.547 8.6 0.547 8.6 0.547],15,1,3,0,0,tempCondList);

%     dof_reduced(i) = hmrstats.dof_blks{1};
%     sse_reduced(i,:,:) = hmrstats.sse_blks{1};
end

% if ~exist(processedDataDir,'dir')
%     mkdir(processedDataDir);
% end
% 
% snirf1.Save([saveDir filesep rawDataFN{1} 'CombinedBasis1.snirf']);
% fileName = [processedDataDir filesep 'intermediateOutputsCombined_Basis1.mat'];
% save(fileName,'dc','dcAvg','dcAvgStd','dcNew','beta',...
%     'betaSS','hmrstats','dod','stim','tRange','dof_full','sse_full',...
%     'tIncAuto','mlActAuto','tIncAuto','snirf1','allS','bvar','dof_reduced',...
%     'sse_reduced');

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

%dcAvgReshape = group.subjs(1,1).runs(1).procStream.output.dcAvg.GetDataTimeSeries('reshape');
%tHRF = group.subjs(1,1).runs.procStream.output.dcAvg.GetTime();

% numConds = 6;
% numConc = 3;
% numChns = size(dcAvg.measurementList,2)/(numConds*numConc);
% sizeCond = size(dcAvg.measurementList,2)/numConds;
% sizeChn = sizeCond/numChns;
% time x (cond x chn x conc)
% Look at S1D1, S1D2, S1D3, & S1D4
% For validation first. Plot 4 channels for one condition in one figure
% for i = 1:numConds
% for i = 1:1
%     figure(); hold on;
%     % for j = 1:numChns
%     for j = 1:1
%         plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i-1) + sizeChn*(j-1) + 1));
%         plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i-1) + sizeChn*(j-1) + 2));
%         plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i-1) + sizeChn*(j-1) + 3)); hold off;
%     end
%     legend('HbO','HbR','HbT');
%     srcStr = num2str(dcAvg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).sourceIndex);
%     detectorStr = num2str(dcAvg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).detectorIndex);
%     %label = num2str(dcavg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).dataTypeLabel);
%     condStr = snirf1.stim(1,i).name;
%     %titleStr = [condStr 'S' srcStr 'D' detectorStr];
%     title(sprintf('%s S%sD%s',condStr,srcStr,detectorStr));
% end

% % For actual analysis, plot 6 conditions from one channel in one figure
% colorList = {[0, 0.447, 0.741],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560] 	};
% for i = 1:6
%     figure(); hold on;
%     % for j = 1:numChns
%     for j = 1:4
%         plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i-1) + sizeChn*(j-1) + 1),'-','Color',colorList{j});
%         plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i-1) + sizeChn*(j-1) + 2),':','Color',colorList{j});
%         plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i-1) + sizeChn*(j-1) + 3),'-.','Color',colorList{j});
%     end
%     legend('S1D1 HbO','S1D1 HbR','S1D1 HbT','S1D2 HbO','S1D2 HbR','S1D2 HbT','S1D3 HbO','S1D3 HbR','S1D3 HbT','S1D4 HbO','S1D4 HbR','S1D4 HbT');
%     srcStr = num2str(dcAvg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).sourceIndex);
%     detectorStr = num2str(dcAvg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).detectorIndex);
%     %label = num2str(dcavg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).dataTypeLabel);
%     condStr = snirf2.stim(1,i).name;
%     %title(sprintf('%s S%sD%s',condStr,srcStr,detectorStr));
%     title(sprintf('%s',condStr));
%     hold off;
% end


%save('Experiment03.mat','dod','dc','dodAvg','dcAvg','dodAvgStd','dcAvgStd','dodSum2','dcSum2','tHRF','nTrials','ch','grpAvgPass','misc');
% output is a ProcResultClass
% save('Experiment03.mat','output');
end