% split into single and multi groups.
% new. Use trials from dcNew var
% Use different method of concatenating to validate the code in
% calcCumSumUpdateDCNew.m
% work for all sbj 08-16

function calcSlopeValue(sbjNum,numClasses,saveOp)
% channels pruned stored in mlActAuto var
processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
rawDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
figSaveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' ...
    num2str(sbjNum) '\Figures\Classification'];
load([processedDataDir filesep 'intermediateOutputsCombined_Basis1.mat']);

fs = 50;
%timeLngth = 1.5*fs;
%timePt = (0:0.5*fs:12*fs)+2*fs;
startT = 2*fs;
timePt = (0:0.5*fs:7*fs)+2*fs;
zeroT = 2*fs;
numChn = 30;
chnName = getChnName(numChn);
cmap = jet(numChn);
if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
    isSDNew = 0;
else
    isSDNew = 1;
end

% Pick 4 channels
numConds = 6;
numConc = 3;
numChns = size(dcAvg.measurementList,2)/(numConds*numConc);
sizeCond = size(dcAvg.measurementList,2)/numConds;
sizeChn = sizeCond/numChns;

% Cell is array of 6 different conditions
% Each array is channels x time x trial
if numClasses == 2
    load([processedDataDir filesep 'singleTrialsUpdated_LR_Basis1.mat'],'singleTrialHRFHbOM',...
        'singleTrialHRFHbRM',...
        'singleTrialHRFHbTM','indexMoviesTest');
else
    load([processedDataDir filesep 'singleTrialsUpdated_Basis1.mat'],'singleTrialHRFHbOM',...
        'singleTrialHRFHbRM',...
        'singleTrialHRFHbTM','indexMoviesTest');
end

timePt = (0:0.5*fs:12*fs)+2*fs;

zeroT = 2*fs;

if strcmp(sbjNum,'15')
    indexMoviesTest = indexMoviesTest(2:end,:);
end

% Offset
for i1 = 1:size(singleTrialHRFHbOM,1)
    for i2 = 1:size(singleTrialHRFHbOM,3)
        offset = singleTrialHRFHbOM(i1,zeroT,i2);
        %offset = mean(singleTrials(i1,1:zeroT,i2));
        singleTrialHRFHbOM(i1,:,i2) = singleTrialHRFHbOM(i1,:,i2) - offset;
        
        offset = singleTrialHRFHbRM(i1,zeroT,i2);
        singleTrialHRFHbRM(i1,:,i2) = singleTrialHRFHbRM(i1,:,i2) - offset;
        
        offset = singleTrialHRFHbTM(i1,zeroT,i2);
        singleTrialHRFHbTM(i1,:,i2) = singleTrialHRFHbTM(i1,:,i2) - offset;
    end
end

% after convert2SD2, cut down from 42 chns to 36 chns. Remove old chns
if ~isSDNew
    [singleTrialHRFHbOM,mlList] = convert2SD2(singleTrialHRFHbOM,mlActAuto{1});
    [singleTrialHRFHbRM,~] = convert2SD2(singleTrialHRFHbRM);
    [singleTrialHRFHbTM,~] = convert2SD2(singleTrialHRFHbTM);
else
    mlList = mlActAuto{1};
end

%isSDNew = 1;

% after selectRS, cut down from 36 to 30 chns. Remove SS
[singleTrialHRFHbOM,mlList] = selectRS(singleTrialHRFHbOM,1,mlList);
[singleTrialHRFHbRM,~] = selectRS(singleTrialHRFHbRM,1);
[singleTrialHRFHbTM,~] = selectRS(singleTrialHRFHbTM,1);

behFN = [rawDataDir filesep 'responses_' sbjNum];

if ~isSDNew
    [singleTrialHRFHbOM,movieIdx,behScore] = keepCorrectTrials_Special(sbjNum,singleTrialHRFHbOM,behFN,numClasses,indexMoviesTest);
    [singleTrialHRFHbRM] = keepCorrectTrials_Special(sbjNum,singleTrialHRFHbRM,behFN,numClasses,indexMoviesTest);
    [singleTrialHRFHbTM] = keepCorrectTrials_Special(sbjNum,singleTrialHRFHbTM,behFN,numClasses,indexMoviesTest);
else
    [singleTrialHRFHbOM,movieIdx,behScore] = keepCorrectTrials(sbjNum,singleTrialHRFHbOM,behFN,numClasses,indexMoviesTest);
    [singleTrialHRFHbRM] = keepCorrectTrials(sbjNum,singleTrialHRFHbRM,behFN,numClasses,indexMoviesTest);
    [singleTrialHRFHbTM] = keepCorrectTrials(sbjNum,singleTrialHRFHbTM,behFN,numClasses,indexMoviesTest);
end

slopeValueHRFHbOM = zeros(size(singleTrialHRFHbOM,1),1,size(singleTrialHRFHbOM,3));
slopeValueHRFHbRM = zeros(size(singleTrialHRFHbRM,1),1,size(singleTrialHRFHbRM,3));
slopeValueHRFHbTM = zeros(size(singleTrialHRFHbTM,1),1,size(singleTrialHRFHbTM,3));

startT = 2*fs;

performanceArrMultiHbO = zeros(size(singleTrialHRFHbOM,1),length(timePt));
performanceArrMultiHbR = zeros(size(singleTrialHRFHbRM,1),length(timePt));
performanceArrMultiHbT = zeros(size(singleTrialHRFHbTM,1),length(timePt));

predRespMHbO = zeros(size(singleTrialHRFHbOM,1),size(singleTrialHRFHbOM,3),length(timePt));
predRespMHbR = zeros(size(singleTrialHRFHbRM,1),size(singleTrialHRFHbRM,3),length(timePt));
predRespMHbT = zeros(size(singleTrialHRFHbTM,1),size(singleTrialHRFHbTM,3),length(timePt));

for i = 1:length(timePt)
    for j = 1:size(singleTrialHRFHbTM,1)
        for k = 1:size(singleTrialHRFHbTM,3)
            X = [ones(length(startT:timePt(i)),1) (startT:timePt(i))'];
            
            Beta = ((X'*X)\X')*squeeze(singleTrialHRFHbOM(j,startT:timePt(i),k))';
            slopeValueHRFHbOM(j,1,k) = Beta(2);
            
            Beta = ((X'*X)\X')*squeeze(singleTrialHRFHbRM(j,startT:timePt(i),k))';
            slopeValueHRFHbRM(j,1,k) = Beta(2);
            
            Beta = ((X'*X)\X')*squeeze(singleTrialHRFHbTM(j,startT:timePt(i),k))';
            slopeValueHRFHbTM(j,1,k) = Beta(2);
        end
    end
    
    %slopeValueHRFHbORS = cat(2,slopeValueHRFHbOS,slopeValueHRFHbRS);
    %slopeValueHRFHbORS = permute(slopeValueHRFHbORS,[1 3 2]);
    %slopeValueHRFHbORM = cat(2,slopeValueHRFHbOM,slopeValueHRFHbRM);
    %slopeValueHRFHbORM = permute(slopeValueHRFHbORM,[1 3 2]);

    %[performanceArrSingleHbO(:,i),predRespSHbO(:,:,i)] = trainClassifierLinearDiscriminantfNIRS(slopeValueHRFHbOS,mlActAuto{1},numClasses,indexMoviesTestSingle);
    [performanceArrMultiHbO(:,i),predRespMHbO(:,:,i)] = trainClassifierLinearDiscriminantfNIRS(slopeValueHRFHbOM,mlActAuto{1},numClasses,movieIdx);

    %[performanceArrSingleHbR(:,i),predRespSHbR(:,:,i)] = trainClassifierLinearDiscriminantfNIRS(slopeValueHRFHbRS,mlActAuto{1},numClasses,indexMoviesTestSingle);
    [performanceArrMultiHbR(:,i),predRespMHbR(:,:,i)] = trainClassifierLinearDiscriminantfNIRS(slopeValueHRFHbRM,mlActAuto{1},numClasses,movieIdx);

    %[performanceArrSingleHbT(:,i),predRespSHbT(:,:,i)] = trainClassifierLinearDiscriminantfNIRS(slopeValueHRFHbTS,mlActAuto{1},numClasses,indexMoviesTestSingle);
    [performanceArrMultiHbT(:,i),predRespMHbT(:,:,i)] = trainClassifierLinearDiscriminantfNIRS(slopeValueHRFHbTM,mlActAuto{1},numClasses,movieIdx);

    %[performanceArrSingleHbOR(:,i),predRespSHbOR(:,:,i)] = trainClassifierLinearDiscriminantfNIRS(slopeValueHRFHbORS,mlActAuto{1},numClasses,indexMoviesTestSingle);
    %[performanceArrMultiHbOR(:,i),predRespMHbOR(:,:,i)] = trainClassifierLinearDiscriminantfNIRS(slopeValueHRFHbORM,mlActAuto{1},numClasses,indexMoviesTestMultiple);

end

% if numClasses == 2
%     savePerfFN = 'performanceSlopeLinearDiscriminantUpdated_LR.mat';
% else
%     savePerfFN = 'performanceSlopeLinearDiscriminantUpdated.mat';
% end
% save([processedDataDir filesep savePerfFN],'performanceArrMultiHbO','performanceArrSingleHbO',...
%     'performanceArrMultiHbR','performanceArrSingleHbR',...
%     'performanceArrMultiHbT','performanceArrSingleHbT',...
%     'performanceArrSingleHbOR','performanceArrMultiHbOR',...
%     'predRespSHbO','predRespMHbO','singleTrialHRFHbOS','singleTrialHRFHbOM',...
%     'singleTrialHRFHbRS','singleTrialHRFHbRM',...
%     'singleTrialHRFHbTS','singleTrialHRFHbTM',...
%     'predRespSHbR','predRespMHbR',...
%     'predRespSHbT','predRespMHbT',...
%     'predRespSHbOR','predRespMHbOR');

%plotPerformanceVsTime(sbjNum,savePerfFN,processedDataDir,'Slope');
%plotPerformanceVsTimeVarName(sbjNum,performanceArrSingleHbOR,savePerfFN,processedDataDir,'Single Slope',timePt);
%plotPerformanceVsTimeVarName(sbjNum,performanceArrMultiHbOR,savePerfFN,processedDataDir,'Multi Slope',timePt);
figure('units','normalized','outerposition',[0 0 1 1]);hold on;
subplot(1,4,1);hold on;
for j = 1:size(performanceArrMultiHbO,1)
    plot(timePt./fs-2,performanceArrMultiHbO(j,:),'Color',cmap(j,:));
end
ylim([0.3 1]);
xlim([0 7]);
title(sprintf('Sbj %s: Δ[HbO] Multi Basis 1',num2str(sbjNum)));
xlabel('Time [s]');ylabel('Accuracy');
hold off;

subplot(1,4,2);hold on;
for j = 1:size(performanceArrMultiHbR,1)
    plot(timePt./fs-2,performanceArrMultiHbR(j,:),'Color',cmap(j,:));
end
ylim([0.3 1]);
xlim([0 7]);
title(sprintf('Sbj %s: Δ[HbR] Multi Basis 1',num2str(sbjNum)));
xlabel('Time [s]');ylabel('Accuracy');
hold off;

subplot(1,4,3);hold on;
for j = 1:size(performanceArrMultiHbT,1)
    plot(timePt./fs-2,performanceArrMultiHbT(j,:),'Color',cmap(j,:));
end
ylim([0.3 1]);
xlim([0 7]);
title(sprintf('Sbj %s: Δ[HbT] Multi Basis 1',num2str(sbjNum)));
xlabel('Time [s]');ylabel('Accuracy');
hold off;

subplot(1,4,[4]);hold on;
for j = 1:size(performanceArrMultiHbT,1)
    %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
    plot(timePt./fs-2,zeros(1,length(timePt./fs-2)),'Color',cmap(j,:));
end
legend(chnName(1:30));
annotation('textbox', [0.8, 0.05, 0.1, 0.1], 'String', ['Score: ' num2str(behScore) '%'],'FontSize',8);
annotation('textbox', [0.71, 0.2, 0.1, 0.1],'String','calcCumSum\_MultiOnly\_AllBasis.m','FontSize',8);
hold off;

if saveOp
    fn = sprintf('PerformanceSlopeVsTime_LR');
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
end

% if numClasses == 2
%     savePerfFN = 'performanceLinearDiscriminantUpdated_LR.mat';
% else
%     savePerfFN = 'performanceLinearDiscriminantUpdated.mat';
% end
% save([processedDataDir filesep savePerfFN],'performanceArrMultiHbO',...
%     'performanceArrMultiHbR',...
%     'performanceArrMultiHbT',...
%     'predRespMHbO','singleTrialHRFHbOM',...
%     'singleTrialHRFHbRM',...
%     'singleTrialHRFHbTM',...
%     'predRespMHbR',...
%     'predRespMHbT');

end