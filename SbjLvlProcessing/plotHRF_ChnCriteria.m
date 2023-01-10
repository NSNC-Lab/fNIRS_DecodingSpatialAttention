% % Work with sbj 08 and after
%
% Use preprocessFNIRS06_MultiOnly.m results, which use GLM regression on entire
% dataset (no cross-validation, no channel or trials rejection). Update preprocessFNIRS06.m and 
% preprocessFNIRS06_MultiOnly.m to reject channels and trials.
%
% Pick channel with highest CA/lowest p-value and plot that.
%
% for criteriaOp == 1: highest Ca
% for criteriaOp == 2: lowest p-val
% for criteriaOp == 3: FEF location based on Petit and Haxby 1999 fMRI
% study associated with saccade movement
%   Use S2D1 and S5D9
%   Or S1D1 and S5D9

function plotHRF_ChnCriteria(sbjNum,criteriaOp)

figSaveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum) '\Figures\HRFs'];
processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
rawDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
load([processedDataDir filesep 'intermediateOutputsCombined_Basis1.mat'],'dcAvg','dcAvgStd','dcNew','beta','bvar');
%load([rawDir filesep sbjNum '.mat'],'s');

tSize = size(dcNew.dataTimeSeries,1);

% sigAct is HbX x #Channels x #conditions
% sigDiff is HbX x #Channels x #pairs
% beta is also a function. be careful
[sigAct,sigDiff,~,pvalDiff] = calcZTest(beta{1},bvar,tSize,1);
% ssIdx = [7,22,24,26,29,32];
% sigAct(:,ssIdx,:) = [];
% sigDiff(:,ssIdx,:) = [];

sigFN = [processedDataDir filesep 'sig.mat'];
save(sigFN,'sigAct','sigDiff');

numConds = 3;
numConc = 3;
numChns = size(dcAvg.measurementList,2)/(numConds*numConc);
sizeCond = size(dcAvg.measurementList,2)/numConds;
sizeChn = sizeCond/numChns;
%condSubPlotIdx = [1 3 5 2 4 6];
condSubPlotIdx = [1 2 3];

colorIdx = [1 3 2 2 1 3];

lineStyle = {'-' '-' '-' '-' '--' '--'};

% srcIdxGrp = {[1 2 3 4] [5 6] [8 9 10 11] [12 13 14 15] [16 17 18 19] [20 21] [23] [25] ...
%     [27 28 29] [31 32 33] [35 36] [37 38] [40] [42]};

srcIdxGrp = {[1 2 3 4] [5 6] [8 9 10 11] [12 13 14 15] [16 17 18 19] [20 21] [23] [25] ...
    [27 28] [30 31] [33 34] [35 36]};

remapChnNum = [1 2 3 4 5 6 8 9 10 11 12 13 14 15 16 17 18 19 20 21 23 25 ...
    27 28 30 31 33 34 35 36];

fs = 50;
timePt = (0:0.5*fs:12*fs)+2*fs;

%timeIdx = find(timePt==(5*fs+2*fs));

% subplot is 8x6.
% chnName = {'S1D1 Middle Posterior FEF Left',...
%     'S1D2 Posterior to sPCS/tgPSC Left',...
%     'S1D3 Middle FEF Left',...
%     'S1D4 sPCS/tgPCS Left',...
%     'S2D1 Posterior FEF Left',...
%     'S2D2 Inferior to FEF Left',...
%     'S3D3 Anterior FEF Left',...
%     'S3D4 iPCS/tgPCS Left',...
%     'S3D5 Anterior to FEF Left',...
%     'S3D6 iPCS Left',...
%     'S4D9 Anterior FEF Right',...
%     'S4D10 iPCS/tgPCS Right',...
%     'S4D11 Anterior to FEF Right',...
%     'S4D12 iPCS Right',...
%     'S5D7 Middle Posterior FEF Right',...
%     'S5D8 Posterior to sPCS/tgPSC Right',...
%     'S5D9 Middle FEF Right',...
%     'S5D10 sPCS/tgPCS Right',...
%     'S6D7 Posterior to FEF Right',...
%     'S6D8 Inferor to FEF Right',...
%     'S7D16 Posterior STG/PT Left',...
%     'S8D17 Posterior STG/PT Right',...
%     'S9D13 IPS3/IPS2/SPL1 Left',...
%     'S9D14 IPS3/antIPS/IPS4 Left',...
%     'S10D13 IPS3/IPS2/SPL1 Right',...
%     'S10D15	IPS3/antIPS/IPS4 Right',...
%     'S11D13	Superior to IPS3/IPS2/SPL1 Left',...
%     'S11D14	IPS4 Left',...
%     'S12D13	Superior to IPS3/IPS2/SPL1 Right',...
%     'S12D15	IPS4 Right'};

chnName = {'FEF',...
    'sPCS/tgPSC',...
    'FEF/tgPSC/sPCS',...
    'sPCS/tgPCS',...
    'FEF',...
    'tgPSC',...
    'SS1',...
    'tgPCS/sPCS/FEF',...
    'iPCS/tgPCS',...
    'FEF',...
    'iPCS/FEF/cIFS/dIPFC',...11
    'tgPCS/sPCS/FEF',...
    'iPCS/tgPCS',...
    'FEF',...
    'iPCS/FEF/cIFS/dIPFC',...
    'FEF',...
    'sPCS/tgPSC',...
    'FEF/tgPCS/sPCS',...
    'sPCS/tgPCS',...
    'FEF',...
    'tgPCS',...
    'SS2',...
    'SMG/STG/PT',...
    'SS3',...
    'SMG/STG/PT',...
    'SS4',...
    'IPS3/IPS2/SPL1',...
    'IPS3/antIPS/IPS4',...
    'SS5',...
    'IPS3/IPS2/SPL1',...
    'IPS3/antIPS/IPS4',...
    'SS6',...
    'IPS3/IPS2/mSPL',...
    'IPS4/mSPL',...
    'IPS3/IPS2/mSPL',...
    'IPS4/mSPL'};

% {'S1D1 Middle Posterior FEF Left',... 21
%     'S1D2 Posterior to sPCS/tgPSC Left',... 20
%     'S1D3 Middle FEF Left',... 15
%     'S1D4 sPCS/tgPCS Left',... 14
%     'S2D1 Posterior FEF Left',... 27
%     'S2D2 Inferior to FEF Left',... 26
%     'S3D3 Anterior FEF Left',... 9
%     'S3D4 iPCS/tgPCS Left',... 8
%     'S3D5 Anterior to FEF Left',... 3
%     'S3D6 iPCS Left',... 2
%     'S4D9 Anterior FEF Right',... 10
%     'S4D10 iPCS/tgPCS Right',... 11
%     'S4D11 Anterior to FEF Right',... 4
%     'S4D12 iPCS Right',... 5
%     'S5D7 Middle Posterior FEF Right',... 22
%     'S5D8 Posterior to sPCS/tgPSC Right',... 23
%     'S5D9 Middle FEF Right',... 16
%     'S5D10 sPCS/tgPCS Right',... 17
%     'S6D7 Posterior to FEF Right',... 28
%     'S6D8 Inferor to FEF Right',... 29
%     'S7D16 Posterior STG/PT Left',... 31
%     'S8D17 Posterior STG/PT Right',... 36
%     'S9D13 IPS3/IPS2/SPL1 Left',... 45
%     'S9D14 IPS3/antIPS/IPS4 Left',... 44
%     'S10D13 IPS3/IPS2/SPL1 Right',... 46
%     'S10D15	IPS3/antIPS/IPS4 Right',... 47
%     'S11D13	Superior to IPS3/IPS2/SPL1 Left',... 39
%     'S11D14	IPS4 Left',... 38
%     'S12D13	Superior to IPS3/IPS2/SPL1 Right',... 40
%     'S12D15	IPS4 Right'}; 41

% chnSubPlotInd = [21 20 15 14 27 26 9 8 3 2 10 11 4 5 22 23 16 17 28 29 31 ...
%     36 45 44 46 47 39 38 40 41];

colorList = {[0, 0.447, 0.741],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],...
    [0.4940, 0.1840, 0.5560],	[0.4660, 0.6740, 0.1880] 	,[0.3010, 0.7450, 0.9330]};

timeSecPostSh = [dcAvg.time; flipud(dcAvg.time)];

% pick highest CA. Here rerun since it's short
if criteriaOp == 1
    
    numRept = 10;
    
    numClasses = 2;
    
    % from createSingleTrialHRF_Updated_DiffBasis.m or
    % createSingleTrialHRF_Updated_MultiOnly_DiffBasis.m

    % 90 trials
    load([processedDataDir filesep 'singleTrialsUpdated_Basis1.mat'],'singleTrialHRFHbOM',...
        'indexMoviesTest');
    
    load([processedDataDir filesep 'intermediateOutputsCombined_Basis1.mat'],'mlActAuto');

    if strcmp(sbjNum,'08') || strcmp(sbjNum,'10')
        isSDNew = 0;
    else
        isSDNew = 1;
    end

    zeroT = 2*fs;
    
    singleTrialHRFHbOM = offsetTrials(singleTrialHRFHbOM,zeroT);

    % after convert2SD2, cut down from 42 chns to 36 chns. Remove old
    % chns
    if ~isSDNew
        [singleTrialHRFHbOM,mlActAutoNew] = convert2SD2(singleTrialHRFHbOM,mlActAuto{1});
%         [singleTrialHRFHbRM,~] = convert2SD2(singleTrialHRFHbRM);
%         [singleTrialHRFHbTM,~] = convert2SD2(singleTrialHRFHbTM);
    else
        mlActAutoNew = mlActAuto{1};
    end

    % terrible coding practice but it works
    %isSDNew = 1;

    % after selectRS, cut down from 36 to 30 chns. Remove SS
    [singleTrialHRFHbOM,mlList] = selectRS(singleTrialHRFHbOM,1,mlActAutoNew);
%     [singleTrialHRFHbRM,~] = selectRS(singleTrialHRFHbRM,1);
%     [singleTrialHRFHbTM,~] = selectRS(singleTrialHRFHbTM,1);

    behFN = [rawDir filesep 'responses_' sbjNum];

    if ~isSDNew
%         multipleIndex = indexMoviesTest(:,5)==1;
%         indexMoviesTest = indexMoviesTest(multipleIndex,:);
        % here, trials are filtered (90 trials) but indexMoviesTest and
        % behFN aren't (180 trials)
        [singleTrialHRFHbOM,movieIdx] = keepCorrectTrials_Special(sbjNum,singleTrialHRFHbOM,behFN,numClasses,indexMoviesTest);
%         [singleTrialHRFHbRM] = keepCorrectTrials_Special(sbjNum,singleTrialHRFHbRM,behFN,numClasses,indexMoviesTest);
%         [singleTrialHRFHbTM] = keepCorrectTrials_Special(sbjNum,singleTrialHRFHbTM,behFN,numClasses,indexMoviesTest);
    else
        [singleTrialHRFHbOM,movieIdx] = keepCorrectTrials(sbjNum,singleTrialHRFHbOM,behFN,numClasses,indexMoviesTest);
%         [singleTrialHRFHbRM] = keepCorrectTrials(sbjNum,singleTrialHRFHbRM,behFN,numClasses,indexMoviesTest);
%         [singleTrialHRFHbTM] = keepCorrectTrials(sbjNum,singleTrialHRFHbTM,behFN,numClasses,indexMoviesTest);
    end

    tempVarPerfHbO = zeros(size(singleTrialHRFHbOM,1),numRept);
%     tempVarPerfHbR = zeros(length(chnName),numRept);
%     tempVarPerfHbT = zeros(length(chnName),numRept);

    timeStrtInt = 2;
    timeLngth = 1*fs;
    idxT = timePt == (timeStrtInt+2)*fs;
    %for i3=idxT:idx
%    for i3=1:length(timePt)

    thisTimePt = timePt(idxT);
    for i4=1:numRept

        tempVarPerfHbO(:,i4) = trainClassifierLinearDiscriminantfNIRS_Cumsum_DiffStart(singleTrialHRFHbOM,thisTimePt,timeLngth,mlList,numClasses,movieIdx);
%         tempVarPerfHbR(:,i4) = trainClassifierLinearDiscriminantfNIRS_Cumsum_DiffStart(singleTrialHRFHbRM,thisTimePt,timeLngth,mlList,numClasses,movieIdx);
%         tempVarPerfHbT(:,i4) = trainClassifierLinearDiscriminantfNIRS_Cumsum_DiffStart(singleTrialHRFHbTM,thisTimePt,timeLngth,mlList,numClasses,movieIdx);

    end
    
    varPerfHbO = squeeze(mean(tempVarPerfHbO,2));
        
    [~,chnNum] = max(varPerfHbO);
    
elseif criteriaOp == 2
    % pick lowest p-value
    [~,chnNum] = min(squeeze(pvalDiff(1,:,1)));
elseif criteriaOp == 3
    chnNum = 5;
elseif criteriaOp == 4
    chnNum = 18;
end

origChnNum = chnNum;
chnNum = remapChnNum(chnNum);

figure('units','normalized','outerposition',[0 0 1 1]); hold on;

subplot(1,2,1);

hold on;
xline(0);hold on;
xline(5);
yline(0);

for i2 = 1:3
    % ystd = y std err. y CI = (t val) * (y std err)
    inBet = [squeeze(dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(chnNum-1) + 1)...
        +dcAvgStd.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(chnNum-1) + 1)); ...
        flipud(squeeze(dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(chnNum-1) + 1)...
        -dcAvgStd.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(chnNum-1) + 1)))];
    h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');

    set(h(i2),'facealpha',0.1);
end

h2 = zeros(1,3);

% i2 is condition number
for i2 = 1:3
    %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),'-','Color',colorList{j-0});hold on;
    if sigAct(1,chnNum,i2)
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    h2(i2) = plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(chnNum-1) + 1),lineStyle{i2},...
        'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;
    
end

if sigDiff(1,chnNum,1)
    LRMark = '*';
else
    LRMark = ' ';
end

if sigDiff(1,chnNum,2)
    LCMark = '+';
else
    LCMark = ' ';
end

if sigDiff(1,chnNum,3)
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('Sbj %s. Chn %s HbO',sbjNum,chnName{origChnNum});
subtitleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
title(titleStr);
subtitle(subtitleStr);
xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])
xlabel('Time [s]');
ylabel('Δ Concentration x Area [M*mm]');
legend(h2,'Left Multiple','Center Multiple','Right Multiple');

% HbR
subplot(1,2,2);

hold on;
xline(0);hold on;
xline(5);
yline(0);

for i2 = 1:3
    inBet = [squeeze(dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(chnNum-1) + 2)...
        +dcAvgStd.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(chnNum-1) + 2)); ...
        flipud(squeeze(dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(chnNum-1) + 2)...
        -dcAvgStd.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(chnNum-1) + 2)))];
    h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');

    set(h(i2),'facealpha',0.1);
end

h2 = zeros(1,3);

% i2 is condition number
for i2 = 1:3
    %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),'-','Color',colorList{j-0});hold on;
    if sigAct(1,chnNum,i2)
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    h2(i2) = plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(chnNum-1) + 2),lineStyle{i2},...
        'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;
    
end

if sigDiff(2,chnNum,1)
    LRMark = '*';
else
    LRMark = ' ';
end

if sigDiff(2,chnNum,2)
    LCMark = '+';
else
    LCMark = ' ';
end

if sigDiff(2,chnNum,3)
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('Sbj %s. Chn %s HbR',sbjNum,chnName{origChnNum});
subtitleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
title(titleStr);
subtitle(subtitleStr);
xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])
xlabel('Time [s]');
ylabel('Δ Concentration x Area [M*mm]');
legend(h2,'Left Multiple','Center Multiple','Right Multiple');

hold off;
if criteriaOp == 1
    fn = sprintf('HRF_MaxCA');
elseif criteriaOp == 2
    fn = sprintf('HRF_MinPVal');
elseif criteriaOp == 3
    fn = sprintf('HRF_FEF_Left');
elseif criteriaOp == 4
    fn = sprintf('HRF_FEF_Right');
end

if ~exist(figSaveDir,'dir')
    mkdir(figSaveDir);
end

print(gcf,[figSaveDir filesep fn],'-dpng','-r250');

end