% STATUS: active
% 
% SYNTAX:
% plotAllHRF_Grp_ChnCriteria_Fig4_TTest(saveOp,opHBCorrected)
% 
% DESCRIPTION:
% Plot HRFs for all 3 chromophores and spatial locations, grouped by task
%   performance. plotAllHRF_Grp_ChnCriteria_Fig4 used Welch t-test. Here
%   use standard t-test.
% 
% RESTRICTION:
% None.
% 
% INPUTS:
% saveOp - int: option to save figure.
%       0 - don't save
%       1 - save
% opHBCorrected - int: option to use corrected p-value
%
% RETURNED VARIABLES:
% None.
% 
% FILES SAVED:
% save figure of HRFs.
% 
% PLOTTING:
% Figure of HRFs

function plotAllHRF_Grp_ChnCriteria_Fig4_TTest(saveOp,opHBCorrected)

goodSbjList = {'08','12','13','16','19','21','24'};
badSbjList = {'14','15','22','23','25'};
fn = sprintf('HbOR_HRF_AllSbjs');
strGrp = 'All';
figSaveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll\Figures\Classification';

leftChnNum = 1;
leftStr = 'Left';
rightChnNum = 18; % 18 is for S5D9, 17 is for S5D8
rightStr = 'Right';

fs = 50;
opPlotInd = 1;
t = -2:1/fs:15;

timeSecPostSh = [t fliplr(t)];

colorIdx = [1 3 2];

colorList = {[0, 0.447, 0.741],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],...
    [0.4940, 0.1840, 0.5560],	[0.4660, 0.6740, 0.1880] 	,[0.3010, 0.7450, 0.9330]};

colorListLight = {[0, 0.447, 0.741, 0.5],[0.8500, 0.3250, 0.0980, 0.5],[0.9290, 0.6940, 0.1250, 0.5],...
    [0.4940, 0.1840, 0.5560],	[0.4660, 0.6740, 0.1880] 	,[0.3010, 0.7450, 0.9330]};

lineStyle = {'-' '-' '-'};

chnName = {'FEF',...
    'sPCS/tgPSC',...
    'FEF/tgPSC/sPCS',...
    'sPCS/tgPCS',...
    'FEF',...
    'tgPSC',...
    'tgPCS/sPCS/FEF',...
    'iPCS/tgPCS',...
    'FEF',...
    'iPCS/FEF/cIFS/dIPFC',...10
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
    'SMG/STG/PT',...
    'SMG/STG/PT',...
    'IPS3/IPS2/SPL1',...
    'IPS3/antIPS/IPS4',...
    'IPS3/IPS2/SPL1',...
    'IPS3/antIPS/IPS4',...
    'IPS3/IPS2/mSPL',...
    'IPS4/mSPL',...
    'IPS3/IPS2/mSPL',...
    'IPS4/mSPL'};

srcIdxGrp = {[1 2 3 4] [5 6] [8 9 10 11] [12 13 14 15] [16 17 18 19] [20 21] [23] [25] ...
    [27 28] [30 31] [33 34] [35 36]};

totNumChns_Homer3 = 36*3*3;

indHbGood = zeros(size(t,2),totNumChns_Homer3,length(goodSbjList));
indHbBad = zeros(size(t,2),totNumChns_Homer3,length(badSbjList));

chnNum = 0;

% HbX x #Channels x #conditions)
numChrom = 3;
numChns = 36;
numConds = 3;
numCoef = 16;
yAvgGrpGood = zeros(length(goodSbjList),numChrom,numChns,numConds);
yAvgGrpBad = zeros(length(goodSbjList),numChrom,numChns,numConds);
betaGoodGrp = zeros(length(goodSbjList),numChrom,numChns,numConds);
betaVarGoodGrp = zeros(length(goodSbjList),numChrom,numChns,numConds);
betaBadGrp = zeros(length(badSbjList),numChrom,numChns,numConds);
betaVarBadGrp = zeros(length(badSbjList),numChrom,numChns,numConds);
sigActGood = zeros(3,numChns,numConds);
sigDiffGood = zeros(3,numChns,numConds);
sigActBad = zeros(3,numChns,numConds);
sigDiffBad = zeros(3,numChns,numConds);

betaIdxRng = 3:8;
numBeta = length(betaIdxRng);

[~, tRangeStart] = min(abs(t-1));
[~, tRangeEnd] = min(abs(t-6));
tRangeIdx = [tRangeStart:1:tRangeEnd];

for i = 1:length(goodSbjList)
    sbjNum = goodSbjList{i};
    saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
    load([saveDir filesep 'intermediateOutputsCombined_Basis1.mat'],'dcAvg','mlActAuto','beta','bvar','hmrstats');
    %load([saveDir filesep 'intermediateOutputsCombined_Basis1_10Hz.mat'],'dcAvg','mlActAuto','beta','bvar','hmrstats');
    
    betaVar = beta{1};
    
    if strcmp(sbjNum,'08')
        temp = hmrstats.pval_contrast(1:3,:,:);
        temp(:,[29,33,39,40,41,42],:) = [];
        sigActGood = sigActGood + (temp<0.05);
        temp = hmrstats.pval_contrast(4:6,:,:);
        temp(:,[29,33,39,40,41,42],:) = [];
        sigDiffGood = sigDiffGood + (temp<0.05);
    else
        sigActGood = sigActGood + (hmrstats.pval_contrast(1:3,:,:)<0.05);
        sigDiffGood = sigDiffGood + (hmrstats.pval_contrast(4:6,:,:)<0.05);
    end
    
    hrf_GLM = dcAvg.dataTimeSeries;
    mlList = dcAvg.measurementList;
    
    if strcmp(sbjNum,'08') || strcmp(sbjNum,'10')
        isSDNew = 0;
        condSubPlotIdx = [1 2 3];
    else
        isSDNew = 1;
        condSubPlotIdx = [1 2 3];
    end
    
    if ~isSDNew
        % after rmvSingleCond, cut down from 756 chns to 378 chns
        %[hrf_GLM,mlList] = rmvSingleCond(hrf_GLM,origNumConds,mlList);
        % after convert2SD2, cut down from 378 chns to 324 chns
        [hrf_GLM,mlList,mlActAuto] = convert2SD2_cnt(hrf_GLM,mlList,mlActAuto{1});
        %betaVar(:,:,:,[1 2 3]) = [];
        betaVar(:,:,[29,33,39,40,41,42],:) = [];
    else
        mlActAuto = mlActAuto{1};
    end
    
    sizeCond = size(mlList,2)/numConds;
    numChns = size(mlList,2)/(numConds*numChrom);
    sizeChn = sizeCond/numChns;
    
    % dcAvg.dataTimeSeries is time x channel
    for i1 = 1:length(srcIdxGrp)
        for j = 1:length(srcIdxGrp{i1})
            for i2 = 1:3
                indHbGood(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1,i) = ...
                    hrf_GLM(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1);
                indHbGood(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2,i) = ...
                    hrf_GLM(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2);
                indHbGood(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3,i) = ...
                    hrf_GLM(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3);
                yAvgGrpGood(i,1,srcIdxGrp{i1}(j),condSubPlotIdx(i2)) = mean(hrf_GLM(tRangeIdx,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),1);
                yAvgGrpGood(i,2,srcIdxGrp{i1}(j),condSubPlotIdx(i2)) = mean(hrf_GLM(tRangeIdx,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2),1);
                yAvgGrpGood(i,3,srcIdxGrp{i1}(j),condSubPlotIdx(i2)) = mean(hrf_GLM(tRangeIdx,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3),1);
            end
        end
    end
                
    % HbX x #Channels x #conditions)
    %betaAvg = squeeze(mean(beta(betaIdxRng,:,:,:),1));
    betaGoodGrp(i,:,:,:) = squeeze(mean(betaVar(betaIdxRng,:,:,:),1));
    
    if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
        bvarIdx = [1 2 3];
    else
        bvarIdx = [1 2 3];
    end
    
    for i2 = 1:numChns
        for i3 = 1:numConds
            for i4 = 1:numChrom
                betaVarGoodGrp(i,i4,i2,i3) = (1/(numBeta^2))*(sum(bvar((bvarIdx(i3)-1)*numCoef+betaIdxRng,i2,i4),1));
            end
        end
    end
    
end

sigActGood = sigActGood./length(goodSbjList);
sigDiffGood = sigDiffGood./length(goodSbjList);

% bad sbj list
for i = 1:length(badSbjList)
    sbjNum = badSbjList{i};
    saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
    load([saveDir filesep 'intermediateOutputsCombined_Basis1.mat'],'dcAvg','mlActAuto','beta','bvar','hmrstats');
    
    betaVar = beta{1};
    
    sigActBad = sigActBad + (hmrstats.pval_contrast(1:3,:,:)<0.05);
    sigDiffBad = sigDiffBad + (hmrstats.pval_contrast(4:6,:,:)<0.05);
    
    hrf_GLM = dcAvg.dataTimeSeries;
    mlList = dcAvg.measurementList;
    
    if strcmp(sbjNum,'08') || strcmp(sbjNum,'10')
        isSDNew = 0;
        condSubPlotIdx = [1 2 3];
    else
        isSDNew = 1;
        condSubPlotIdx = [1 2 3];
    end
    
    if ~isSDNew
        % after rmvSingleCond, cut down from 756 chns to 378 chns
        %[hrf_GLM,mlList] = rmvSingleCond(hrf_GLM,origNumConds,mlList);
        % after convert2SD2, cut down from 378 chns to 324 chns
        [hrf_GLM,mlList,mlActAuto] = convert2SD2_cnt(hrf_GLM,mlList,mlActAuto{1});
        %betaVar(:,:,:,[1 2 3]) = [];
        betaVar(:,:,[29,33,39,40,41,42],:) = [];
    else
        mlActAuto = mlActAuto{1};
    end
    
    sizeCond = size(mlList,2)/numConds;
    numChns = size(mlList,2)/(numConds*numChrom);
    sizeChn = sizeCond/numChns;
    
    % dcAvg.dataTimeSeries is time x channel
    for i1 = 1:length(srcIdxGrp)
        for j = 1:length(srcIdxGrp{i1})
            for i2 = 1:3
                indHbBad(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1,i) = ...
                    hrf_GLM(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1);
                indHbBad(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2,i) = ...
                    hrf_GLM(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2);
                indHbBad(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3,i) = ...
                    hrf_GLM(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3);
                yAvgGrpBad(i,1,srcIdxGrp{i1}(j),condSubPlotIdx(i2)) = mean(hrf_GLM(tRangeIdx,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),1);
                yAvgGrpBad(i,2,srcIdxGrp{i1}(j),condSubPlotIdx(i2)) = mean(hrf_GLM(tRangeIdx,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2),1);
                yAvgGrpBad(i,3,srcIdxGrp{i1}(j),condSubPlotIdx(i2)) = mean(hrf_GLM(tRangeIdx,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3),1);
            end
        end
    end
                
    % HbX x #Channels x #conditions)
    %betaAvg = squeeze(mean(beta(betaIdxRng,:,:,:),1));
    betaBadGrp(i,:,:,:) = squeeze(mean(betaVar(betaIdxRng,:,:,:),1));
    
    if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
        bvarIdx = [1 2 3];
    else
        bvarIdx = [1 2 3];
    end
    
    for i2 = 1:numChns
        for i3 = 1:numConds
            for i4 = 1:numChrom
                betaVarBadGrp(i,i4,i2,i3) = (1/(numBeta^2))*(sum(bvar((bvarIdx(i3)-1)*numCoef+betaIdxRng,i2,i4),1));
            end
        end
    end
    
end

sigActBad = sigActBad./length(badSbjList);
sigDiffBad = sigDiffBad./length(badSbjList);

% left, center, right cond
% HbO, HbR, HbT
sigActGoodLftTbl = squeeze(sigActGood(:,leftChnNum,:));
sigDiffGoodLftTbl = squeeze(sigDiffGood(:,leftChnNum,:));
sigActGoodRightTbl = squeeze(sigActGood(:,rightChnNum,:));
sigDiffGoodRightTbl = squeeze(sigDiffGood(:,rightChnNum,:));

sigActBadLftTbl = squeeze(sigActBad(:,leftChnNum,:));
sigDiffBadLftTbl = squeeze(sigDiffBad(:,leftChnNum,:));
sigActBadRightTbl = squeeze(sigActBad(:,rightChnNum,:));
sigDiffBadRightTbl = squeeze(sigDiffBad(:,rightChnNum,:));

numPairs = numConds*(numConds-1)/2;
pValActGood = zeros(numChrom,numChns,numConds);
pValActBad = zeros(numChrom,numChns,numConds);
pValDiffGood = zeros(numChrom,numChns,numConds);
pValDiffBad = zeros(numChrom,numChns,numConds);
pValActGoodHBC = zeros(numChrom,numChns,numConds);
pValActBadHBC = zeros(numChrom,numChns,numConds);
pValDiffGoodHBC = zeros(numChrom,numChns,numConds);
pValDiffBadHBC = zeros(numChrom,numChns,numConds);
sigActGood = zeros(numChrom,numChns,numConds);
sigDiffGood = zeros(numChrom,numChns,numPairs);
sigActBad = zeros(numChrom,numChns,numConds);
sigDiffBad = zeros(numChrom,numChns,numPairs);
m = 3;
listPairs = nchoosek([1;2;3],2);

% beta is #sbj x HbX x #Channels x # conditions
for i1 = 1:numChrom
    for i2 = 1:numChns
        for i3 = 1:numConds
            %[~,sigActGood(i1,i2,i3)] = ttest(betaGoodGrp(:,i1,i2,i3));
            %[~,sigActBad(i1,i2,i3)] = ttest(betaBadGrp(:,i1,i2,i3));
            [~,pValActGood(i1,i2,i3),~,stats{i1,i2,i3}] = ttest(yAvgGrpGood(:,i1,i2,i3));
            [~,pValActBad(i1,i2,i3)] = ttest(yAvgGrpBad(:,i1,i2,i3));
        end
    end
end

if opHBCorrected
for i1 = 1:numChrom
    for i2 = 1:numChns
        tempPVal = zeros(1,3);
        [~,idTemp] = sort(pValActGood(i1,i2,:));
        for iK = 1:m
            tempPVal(iK) = min(pValActGood(i1,i2,idTemp(1,1,iK))*(m-iK+1),1);
            pValActGoodHBC(i1,i2,idTemp(1,1,iK)) = max(tempPVal);
        end
        tempPVal = zeros(1,3);
        [~,idTemp] = sort(pValActBad(i1,i2,:));
        for iK = 1:m
            tempPVal(iK) = min(pValActBad(i1,i2,idTemp(1,1,iK))*(m-iK+1),1);
            pValActBadHBC(i1,i2,idTemp(1,1,iK)) = max(tempPVal);
        end
%         for i3 = 1:numConds
%             sigAct(i1,i2,i3) = pValActHBC(i1,i2,i3)<pCV;
%         end
    end
end
else
    pValActGoodHBC = pValActGood;
    pValActBadHBC = pValActBad;
end

for i1 = 1:numChrom
    for i2 = 1:numChns
        for i3 = 1:numPairs
            %[~,sigDiffGood(i1,i2,i3)] = ttest(betaGoodGrp(:,i1,i2,listPairs(i3,1)),betaGoodGrp(:,i1,i2,listPairs(i3,2)));
            %[~,sigDiffBad(i1,i2,i3)] = ttest(betaBadGrp(:,i1,i2,listPairs(i3,1)),betaBadGrp(:,i1,i2,listPairs(i3,2)));
            [~,pValDiffGood(i1,i2,i3)] = ttest2(yAvgGrpGood(:,i1,i2,listPairs(i3,1)),yAvgGrpGood(:,i1,i2,listPairs(i3,2)));
            [~,pValDiffBad(i1,i2,i3)] = ttest2(yAvgGrpBad(:,i1,i2,listPairs(i3,1)),yAvgGrpBad(:,i1,i2,listPairs(i3,2)));
        end
    end
end

if opHBCorrected
for i1 = 1:numChrom
    for i2 = 1:numChns
        tempPVal = zeros(1,3);
        [~,idTemp] = sort(pValDiffGood(i1,i2,:));
        for iK = 1:m
            tempPVal(iK) = min(pValDiffGood(i1,i2,idTemp(1,1,iK))*(m-iK+1),1);
            pValDiffGoodHBC(i1,i2,idTemp(1,1,iK)) = max(tempPVal);
        end
        tempPVal = zeros(1,3);
        [~,idTemp] = sort(pValDiffBad(i1,i2,:));
        for iK = 1:m
            tempPVal(iK) = min(pValDiffBad(i1,i2,idTemp(1,1,iK))*(m-iK+1),1);
            pValDiffBadHBC(i1,i2,idTemp(1,1,iK)) = max(tempPVal);
        end
%         for i3 = 1:numConds
%             sigAct(i1,i2,i3) = pValActHBC(i1,i2,i3)<pCV;
%         end
    end
end
else
    pValDiffGoodHBC = pValDiffGood;
    pValDiffBadHBC = pValDiffBad;
end

% [h,p,ci,stats] = ttest(x);
% [h,p,ci,stats] = ttest(x);
% [sigActGood,sigDiffGood] = calcWelchTTest(betaGoodGrp,betaVarGoodGrp);
% [sigActBad,sigDiffBad] = calcWelchTTest(betaBadGrp,betaVarBadGrp);

critPVal = 0.05;
sigFN = [saveDir filesep 'sig_' strGrp '.mat'];
save(sigFN,'sigActGood','sigDiffGood','sigActBad','sigDiffBad');

grpHbGood = squeeze(mean(indHbGood,3));
%grpHbStd = squeeze(std(indHb,0,3));
grpHbGoodStdErr = std(indHbGood,0,3)/(size(indHbGood,3)-1);
% symmetric in this case (not always symmetric)
grpHbGoodCIU = squeeze(grpHbGoodStdErr.*tinv(0.975,size(indHbGood,3)-1));
grpHbGoodCIL = squeeze(grpHbGoodStdErr.*tinv(0.025,size(indHbGood,3)-1));

grpHbBad = squeeze(mean(indHbBad,3));
grpHbBadStdErr = std(indHbBad,0,3)/(size(indHbBad,3)-1);
grpHbBadCIU = squeeze(grpHbBadStdErr.*tinv(0.975,size(indHbBad,3)-1));
grpHbBadCIL = squeeze(grpHbBadStdErr.*tinv(0.025,size(indHbBad,3)-1));

figure();hold on;

condSubPlotIdx = [1 2 3];

% Left HbO Good
subplot(3,4,1);hold on;
        
xline(0);hold on;
xline(5);
yline(0);
h = zeros(1,3);

for i2 = 1:3

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(leftChnNum-1) + 1;
    inBet = [squeeze(grpHbGood(:,chnIdx) + grpHbGoodCIU(:,chnIdx)); ...
        flipud(squeeze(grpHbGood(:,chnIdx) + grpHbGoodCIL(:,chnIdx)))];
    h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');
    
    set(h(i2),'facealpha',0.1);

    if pValActGoodHBC(1,leftChnNum,i2) < critPVal
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    plot(t,grpHbGood(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if pValDiffGoodHBC(1,leftChnNum,1)  < critPVal
    LRMark = '*';
else
    LRMark = ' ';
end

if pValDiffGoodHBC(1,leftChnNum,2)  < critPVal
    LCMark = '+';
else
    LCMark = ' ';
end

if pValDiffGoodHBC(1,leftChnNum,3)  < critPVal
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
%title(['HbO ' chnName(leftChnNum)]);
title(['HbO Left High']);
subtitle(titleStr);

% Left HbO Bad
subplot(3,4,3);hold on;
        
xline(0);hold on;
xline(5);
yline(0);
h = zeros(1,3);

for i2 = 1:3

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(leftChnNum-1) + 1;
    inBet = [squeeze(grpHbBad(:,chnIdx) + grpHbBadCIU(:,chnIdx)); ...
        flipud(squeeze(grpHbBad(:,chnIdx) + grpHbBadCIL(:,chnIdx)))];
    h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');
    
    set(h(i2),'facealpha',0.1);

    if pValActBadHBC(1,leftChnNum,i2)  < critPVal
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    plot(t,grpHbBad(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if pValDiffBadHBC(1,leftChnNum,1)  < critPVal
    LRMark = '*';
else
    LRMark = ' ';
end

if pValDiffBadHBC(1,leftChnNum,2)  < critPVal
    LCMark = '+';
else
    LCMark = ' ';
end

if pValDiffBadHBC(1,leftChnNum,3)  < critPVal
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
%title(['HbO ' chnName(leftChnNum)]);
title(['HbO Left Low']);
subtitle(titleStr);

% Left HbR Good
subplot(3,4,5);hold on;

xline(0);
xline(5);
yline(0);

for i2 = 1:3
    
    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(leftChnNum-1) + 2;
    
    inBet = [squeeze(grpHbGood(:,chnIdx)+grpHbGoodCIU(:,chnIdx)); ...
        flipud(squeeze(grpHbGood(:,chnIdx)+grpHbGoodCIL(:,chnIdx)))];
    h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');hold on;

    set(h(i2),'facealpha',0.1); hold on;
end

h2 = zeros(1,3);

for i2 = 1:3

    if pValActGoodHBC(2,leftChnNum,i2) < critPVal
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(leftChnNum-1) + 2;
    h2(i2) = plot(t,grpHbGood(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if pValDiffGoodHBC(2,leftChnNum,1) < critPVal
    LRMark = '*';
else
    LRMark = ' ';
end

if pValDiffGoodHBC(2,leftChnNum,2) < critPVal
    LCMark = '+';
else
    LCMark = ' ';
end

if pValDiffGoodHBC(2,leftChnNum,3) < critPVal
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
%title(['HbR ' chnName(leftChnNum)]);
title(['HbR Left High']);
subtitle(titleStr);

% Left HbR Bad
subplot(3,4,7);hold on;

xline(0);
xline(5);
yline(0);

for i2 = 1:3
    
    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(leftChnNum-1) + 2;
    
    inBet = [squeeze(grpHbBad(:,chnIdx)+grpHbBadCIU(:,chnIdx)); ...
        flipud(squeeze(grpHbBad(:,chnIdx)+grpHbBadCIL(:,chnIdx)))];
    h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');hold on;

    set(h(i2),'facealpha',0.1); hold on;
end

h2 = zeros(1,3);

for i2 = 1:3

    if pValActBadHBC(2,leftChnNum,i2) < critPVal
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(leftChnNum-1) + 2;
    h2(i2) = plot(t,grpHbBad(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if pValDiffBadHBC(2,leftChnNum,1) < critPVal
    LRMark = '*';
else
    LRMark = ' ';
end

if pValDiffBadHBC(2,leftChnNum,2) < critPVal
    LCMark = '+';
else
    LCMark = ' ';
end

if pValDiffBadHBC(2,leftChnNum,3) < critPVal
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
%title(['HbR ' chnName(leftChnNum)]);
title(['HbO Left Low']);
subtitle(titleStr);

% Left HbT Good
subplot(3,4,9);hold on;

xline(0);
xline(5);
yline(0);

for i2 = 1:3
    
    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(leftChnNum-1) + 3;
    
    inBet = [squeeze(grpHbGood(:,chnIdx)+grpHbGoodCIU(:,chnIdx)); ...
        flipud(squeeze(grpHbGood(:,chnIdx)+grpHbGoodCIL(:,chnIdx)))];
    h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');hold on;

    set(h(i2),'facealpha',0.1); hold on;
end

h2 = zeros(1,3);

for i2 = 1:3

    if pValActGoodHBC(3,leftChnNum,i2) < critPVal
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(leftChnNum-1) + 3;
    h2(i2) = plot(t,grpHbGood(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if pValDiffGoodHBC(3,leftChnNum,1) < critPVal
    LRMark = '*';
else
    LRMark = ' ';
end

if pValDiffGoodHBC(3,leftChnNum,2) < critPVal
    LCMark = '+';
else
    LCMark = ' ';
end

if pValDiffGoodHBC(3,leftChnNum,3) < critPVal
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
%title(['HbR ' chnName(leftChnNum)]);
title(['HbT Left High']);
subtitle(titleStr);

% Left HbT Bad
subplot(3,4,11);hold on;

xline(0);
xline(5);
yline(0);

for i2 = 1:3
    
    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(leftChnNum-1) + 3;
    
    inBet = [squeeze(grpHbBad(:,chnIdx)+grpHbBadCIU(:,chnIdx)); ...
        flipud(squeeze(grpHbBad(:,chnIdx)+grpHbBadCIL(:,chnIdx)))];
    h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');hold on;

    set(h(i2),'facealpha',0.1); hold on;
end

h2 = zeros(1,3);

for i2 = 1:3

    if pValActBadHBC(3,leftChnNum,i2) < critPVal
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(leftChnNum-1) + 3;
    h2(i2) = plot(t,grpHbBad(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if pValDiffBadHBC(3,leftChnNum,1) < critPVal
    LRMark = '*';
else
    LRMark = ' ';
end

if pValDiffBadHBC(3,leftChnNum,2) < critPVal
    LCMark = '+';
else
    LCMark = ' ';
end

if pValDiffBadHBC(3,leftChnNum,3) < critPVal
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
%title(['HbR ' chnName(leftChnNum)]);
title(['HbT Left Low']);
subtitle(titleStr);

% Right HbO Good
subplot(3,4,2);hold on;
        
xline(0);hold on;
xline(5);
yline(0);
h = zeros(1,3);

for i2 = 1:3

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(rightChnNum-1) + 1;
    inBet = [squeeze(grpHbGood(:,chnIdx) + grpHbGoodCIU(:,chnIdx)); ...
        flipud(squeeze(grpHbGood(:,chnIdx) + grpHbGoodCIL(:,chnIdx)))];
    h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');
    
    set(h(i2),'facealpha',0.1);

    if pValActGoodHBC(1,rightChnNum,i2) < critPVal
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    plot(t,grpHbGood(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if pValDiffGoodHBC(1,rightChnNum,1) < critPVal
    LRMark = '*';
else
    LRMark = ' ';
end

if pValDiffGoodHBC(1,rightChnNum,2) < critPVal
    LCMark = '+';
else
    LCMark = ' ';
end

if pValDiffGoodHBC(1,rightChnNum,3) < critPVal
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
%title(['HbO ' chnName(rightChnNum)]);
title(['HbO Right High']);
subtitle(titleStr);

% Right HbO Bad
subplot(3,4,4);hold on;
        
xline(0);hold on;
xline(5);
yline(0);
h = zeros(1,3);

for i2 = 1:3

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(rightChnNum-1) + 1;
    inBet = [squeeze(grpHbBad(:,chnIdx) + grpHbBadCIU(:,chnIdx)); ...
        flipud(squeeze(grpHbBad(:,chnIdx) + grpHbBadCIL(:,chnIdx)))];
    h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');
    
    set(h(i2),'facealpha',0.1);

    if pValActBadHBC(1,rightChnNum,i2) < critPVal
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    plot(t,grpHbBad(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if pValDiffBadHBC(1,rightChnNum,1) < critPVal
    LRMark = '*';
else
    LRMark = ' ';
end

if pValDiffBadHBC(1,rightChnNum,2) < critPVal
    LCMark = '+';
else
    LCMark = ' ';
end

if pValDiffBadHBC(1,rightChnNum,3) < critPVal
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
%title(['HbO ' chnName(rightChnNum)]);
title(['HbO Right Low']);
subtitle(titleStr);

% Right HbR Good
subplot(3,4,6);hold on;

xline(0);
xline(5);
yline(0);

for i2 = 1:3
    
    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(rightChnNum-1) + 2;
    
    inBet = [squeeze(grpHbGood(:,chnIdx)+grpHbGoodCIU(:,chnIdx)); ...
        flipud(squeeze(grpHbGood(:,chnIdx)+grpHbGoodCIL(:,chnIdx)))];
    h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');hold on;

    set(h(i2),'facealpha',0.1); hold on;
end

h2 = zeros(1,3);

for i2 = 1:3

    if pValActGoodHBC(2,rightChnNum,i2) < critPVal
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(rightChnNum-1) + 2;
    h2(i2) = plot(t,grpHbGood(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if pValDiffGoodHBC(2,rightChnNum,1) < critPVal
    LRMark = '*';
else
    LRMark = ' ';
end

if pValDiffGoodHBC(2,rightChnNum,2) < critPVal
    LCMark = '+';
else
    LCMark = ' ';
end

if pValDiffGoodHBC(2,rightChnNum,3) < critPVal
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
%title(['HbR ' chnName(rightChnNum)]);
title(['HbR Right High']);
subtitle(titleStr);

% Right HbR Bad
subplot(3,4,8);hold on;

xline(0);
xline(5);
yline(0);

for i2 = 1:3
    
    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(rightChnNum-1) + 2;
    
    inBet = [squeeze(grpHbBad(:,chnIdx)+grpHbBadCIU(:,chnIdx)); ...
        flipud(squeeze(grpHbBad(:,chnIdx)+grpHbBadCIL(:,chnIdx)))];
    h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');hold on;

    set(h(i2),'facealpha',0.1); hold on;
end

h2 = zeros(1,3);

for i2 = 1:3

    if pValActBadHBC(2,rightChnNum,i2) < critPVal
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(rightChnNum-1) + 2;
    h2(i2) = plot(t,grpHbBad(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if pValDiffBadHBC(2,rightChnNum,1) < critPVal
    LRMark = '*';
else
    LRMark = ' ';
end

if pValDiffBadHBC(2,rightChnNum,2) < critPVal
    LCMark = '+';
else
    LCMark = ' ';
end

if pValDiffBadHBC(2,rightChnNum,3) < critPVal
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
%title(['HbR ' chnName(rightChnNum)]);
title(['HbR Right Low']);
subtitle(titleStr);

% Right HbT Good
subplot(3,4,10);hold on;

xline(0);
xline(5);
yline(0);

for i2 = 1:3
    
    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(rightChnNum-1) + 3;
    
    inBet = [squeeze(grpHbGood(:,chnIdx)+grpHbGoodCIU(:,chnIdx)); ...
        flipud(squeeze(grpHbGood(:,chnIdx)+grpHbGoodCIL(:,chnIdx)))];
    h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');hold on;

    set(h(i2),'facealpha',0.1); hold on;
end

h2 = zeros(1,3);

for i2 = 1:3

    if pValActGoodHBC(3,rightChnNum,i2) < critPVal
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(rightChnNum-1) + 3;
    h2(i2) = plot(t,grpHbGood(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if pValDiffGoodHBC(3,rightChnNum,1) < critPVal
    LRMark = '*';
else
    LRMark = ' ';
end

if pValDiffGoodHBC(3,rightChnNum,2) < critPVal
    LCMark = '+';
else
    LCMark = ' ';
end

if pValDiffGoodHBC(3,rightChnNum,3) < critPVal
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
%title(['HbR ' chnName(rightChnNum)]);
title(['HbT Right High']);
subtitle(titleStr);

% Right HbT Bad
subplot(3,4,12);hold on;

xline(0);
xline(5);
yline(0);

for i2 = 1:3
    
    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(rightChnNum-1) + 3;
    
    inBet = [squeeze(grpHbBad(:,chnIdx)+grpHbBadCIU(:,chnIdx)); ...
        flipud(squeeze(grpHbBad(:,chnIdx)+grpHbBadCIL(:,chnIdx)))];
    h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');hold on;

    set(h(i2),'facealpha',0.1); hold on;
end

h2 = zeros(1,3);

for i2 = 1:3

    if pValActBadHBC(3,rightChnNum,i2) < critPVal
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(rightChnNum-1) + 3;
    h2(i2) = plot(t,grpHbBad(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if pValDiffBadHBC(3,rightChnNum,1) < critPVal
    LRMark = '*';
else
    LRMark = ' ';
end

if pValDiffBadHBC(3,rightChnNum,2) < critPVal
    LCMark = '+';
else
    LCMark = ' ';
end

if pValDiffBadHBC(3,rightChnNum,3) < critPVal
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
%title(['HbR ' chnName(rightChnNum)]);
title(['HbT Right Low']);
subtitle(titleStr);


%legend({'Left Single','Center Single','Right Multiple','Right Single','Left Multiple','Center Multiple'},'Location','northeast');
legend(h2,'Left Multiple','Center Multiple','Right Multiple');
%legend('Left Multiple','Center Multiple','Right Multiple');

sgtitle(sprintf('FEF By Behavioral Performance'));
hold off;

if saveOp
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% OLD CODE %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%