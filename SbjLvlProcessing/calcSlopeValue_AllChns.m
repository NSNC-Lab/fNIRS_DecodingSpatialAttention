% split into single and multi groups.
% new. Use trials from dcNew var
% Use different method of concatenating to validate the code in
% calcCumSumUpdateDCNew.m
% work for all sbj 08-16

function calcSlopeValue_AllChns(sbjNum,numClasses,saveOp)
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

slopeValueHRFHbOM = zeros(size(singleTrialHRFHbOM,1),size(singleTrialHRFHbOM,3));
slopeValueHRFHbRM = zeros(size(singleTrialHRFHbRM,1),size(singleTrialHRFHbRM,3));
slopeValueHRFHbTM = zeros(size(singleTrialHRFHbTM,1),size(singleTrialHRFHbTM,3));

startT = 2*fs;

performanceLDALedoitHbO = zeros(length(timePt),1);
performanceLDALedoitHbR = zeros(length(timePt),1);
performanceLDALedoitHbT = zeros(length(timePt),1);

performanceLDACERNNHbO = zeros(length(timePt),1);
performanceLDACERNNHbR = zeros(length(timePt),1);
performanceLDACERNNHbT = zeros(length(timePt),1);

performanceLogRegHbO = zeros(length(timePt),1);
performanceLogRegHbR = zeros(length(timePt),1);
performanceLogRegHbT = zeros(length(timePt),1);

performanceSVMHbO = zeros(length(timePt),1);
performanceSVMHbR = zeros(length(timePt),1);
performanceSVMHbT = zeros(length(timePt),1);

performanceBaggingHbO = zeros(length(timePt),1);
performanceBaggingHbR = zeros(length(timePt),1);
performanceBaggingHbT = zeros(length(timePt),1);

performanceBoostedHbO = zeros(length(timePt),1);
performanceBoostedHbR = zeros(length(timePt),1);
performanceBoostedHbT = zeros(length(timePt),1);

predRespMHbO = zeros(size(singleTrialHRFHbOM,1),size(singleTrialHRFHbOM,3),length(timePt));
predRespMHbR = zeros(size(singleTrialHRFHbRM,1),size(singleTrialHRFHbRM,3),length(timePt));
predRespMHbT = zeros(size(singleTrialHRFHbTM,1),size(singleTrialHRFHbTM,3),length(timePt));

for i = 2:length(timePt)
    for j = 1:size(singleTrialHRFHbTM,1)
        for k = 1:size(singleTrialHRFHbTM,3)
            X = [ones(length(startT:timePt(i)),1) (startT:timePt(i))'];
            
            Beta = ((X'*X)\X')*squeeze(singleTrialHRFHbOM(j,startT:timePt(i),k))';
            slopeValueHRFHbOM(j,k) = Beta(2);
            
            Beta = ((X'*X)\X')*squeeze(singleTrialHRFHbRM(j,startT:timePt(i),k))';
            slopeValueHRFHbRM(j,k) = Beta(2);
            
            Beta = ((X'*X)\X')*squeeze(singleTrialHRFHbTM(j,startT:timePt(i),k))';
            slopeValueHRFHbTM(j,k) = Beta(2);
        end
    end

    performanceLDALedoitHbO(i) = train_RLDA_Ledoit(slopeValueHRFHbOM,mlList,numChn,numClasses,movieIdx);
    performanceLDALedoitHbR(i) = train_RLDA_Ledoit(slopeValueHRFHbRM,mlList,numChn,numClasses,movieIdx);
    performanceLDALedoitHbT(i) = train_RLDA_Ledoit(slopeValueHRFHbTM,mlList,numChn,numClasses,movieIdx);
    
%     performanceLDACERNNHbO(i) = train_RLDA_CERNN(slopeValueHRFHbOM,mlList,numChn,numClasses,movieIdx);
%     performanceLDACERNNHbR(i) = train_RLDA_CERNN(slopeValueHRFHbRM,mlList,numChn,numClasses,movieIdx);
%     performanceLDACERNNHbT(i) = train_RLDA_CERNN(slopeValueHRFHbTM,mlList,numChn,numClasses,movieIdx);
    
    performanceLogRegHbO(i) = trainClassifierLogisticRegressionHRF(slopeValueHRFHbOM,mlList,numChn,numClasses,movieIdx);
    performanceLogRegHbR(i) = trainClassifierLogisticRegressionHRF(slopeValueHRFHbRM,mlList,numChn,numClasses,movieIdx);
    performanceLogRegHbT(i) = trainClassifierLogisticRegressionHRF(slopeValueHRFHbTM,mlList,numChn,numClasses,movieIdx);
    
    performanceSVMHbO(i) = trainClassifierSVMHRF(slopeValueHRFHbOM,mlList,numChn,numClasses,movieIdx);
    performanceSVMHbR(i) = trainClassifierSVMHRF(slopeValueHRFHbRM,mlList,numChn,numClasses,movieIdx);
    performanceSVMHbT(i) = trainClassifierSVMHRF(slopeValueHRFHbTM,mlList,numChn,numClasses,movieIdx);
    
    performanceBaggingHbO(i) = trainClassifierBaggingTreesSingleHRF(slopeValueHRFHbOM,mlList,numChn,numClasses,movieIdx);
    performanceBaggingHbR(i) = trainClassifierBaggingTreesSingleHRF(slopeValueHRFHbRM,mlList,numChn,numClasses,movieIdx);
    performanceBaggingHbT(i) = trainClassifierBaggingTreesSingleHRF(slopeValueHRFHbTM,mlList,numChn,numClasses,movieIdx);
    
    performanceBoostedHbO(i) = trainClassifierBoostedTreesSingleHRF(slopeValueHRFHbOM,mlList,numChn,numClasses,movieIdx);
    performanceBoostedHbR(i) = trainClassifierBoostedTreesSingleHRF(slopeValueHRFHbRM,mlList,numChn,numClasses,movieIdx);
    performanceBoostedHbT(i) = trainClassifierBoostedTreesSingleHRF(slopeValueHRFHbTM,mlList,numChn,numClasses,movieIdx);
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

numClassifiers = 6;
cmap = jet(numClassifiers);

figure('units','normalized','outerposition',[0 0 1 1]);hold on;
subplot(1,3,1);hold on;
plot(timePt./fs-2,performanceLDALedoitHbO,'Color',cmap(1,:));
plot(timePt./fs-2,performanceLDACERNNHbO,'Color',cmap(2,:));
plot(timePt./fs-2,performanceLogRegHbO,'Color',cmap(3,:));
plot(timePt./fs-2,performanceSVMHbO,'Color',cmap(4,:));
plot(timePt./fs-2,performanceBaggingHbO,'Color',cmap(5,:));
plot(timePt./fs-2,performanceBoostedHbO,'Color',cmap(6,:));
ylim([0.3 1]);
xlim([0 7]);
title(sprintf('Sbj %s: Δ[HbO] Multi Basis 1',num2str(sbjNum)));
legend({'LDAShrink','LDA CERNN','LogReg','SVM','Bagging','Boosted'});
xlabel('Time [s]');ylabel('Accuracy');
hold off;

subplot(1,3,2);hold on;
plot(timePt./fs-2,performanceLDALedoitHbR,'Color',cmap(1,:));
plot(timePt./fs-2,performanceLDACERNNHbR,'Color',cmap(2,:));
plot(timePt./fs-2,performanceLogRegHbR,'Color',cmap(3,:));
plot(timePt./fs-2,performanceSVMHbR,'Color',cmap(4,:));
plot(timePt./fs-2,performanceBaggingHbR,'Color',cmap(5,:));
plot(timePt./fs-2,performanceBoostedHbR,'Color',cmap(6,:));
ylim([0.3 1]);
xlim([0 7]);
title(sprintf('Sbj %s: Δ[HbR] Multi Basis 1',num2str(sbjNum)));
xlabel('Time [s]');ylabel('Accuracy');
hold off;

subplot(1,3,3);hold on;
plot(timePt./fs-2,performanceLDALedoitHbT,'Color',cmap(1,:));
plot(timePt./fs-2,performanceLDACERNNHbT,'Color',cmap(2,:));
plot(timePt./fs-2,performanceLogRegHbT,'Color',cmap(3,:));
plot(timePt./fs-2,performanceSVMHbT,'Color',cmap(4,:));
plot(timePt./fs-2,performanceBaggingHbT,'Color',cmap(5,:));
plot(timePt./fs-2,performanceBoostedHbT,'Color',cmap(6,:));
ylim([0.3 1]);
xlim([0 7]);
title(sprintf('Sbj %s: Δ[HbT] Multi Basis 1',num2str(sbjNum)));
xlabel('Time [s]');ylabel('Accuracy');
hold off;

% subplot(1,4,[4]);hold on;
% for j = 1:size(performanceArrMultiHbT,1)
%     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
%     plot(timePt./fs-2,zeros(1,length(timePt./fs-2)),'Color',cmap(j,:));
% end
% legend(chnName(1:30));
% annotation('textbox', [0.8, 0.05, 0.1, 0.1], 'String', ['Score: ' num2str(behScore) '%'],'FontSize',8);
% annotation('textbox', [0.71, 0.2, 0.1, 0.1],'String','calcCumSum\_MultiOnly\_AllBasis.m','FontSize',8);
% hold off;

if saveOp
    fn = sprintf('PerformanceSlopeVsTime_AllChns_LR');
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