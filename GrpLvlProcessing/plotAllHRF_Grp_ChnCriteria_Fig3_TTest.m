% STATUS: active
% 
% SYNTAX:
% plotAllHRF_Grp_ChnCriteria_Fig3_TTest(saveOp,opHBCorrected)
% 
% DESCRIPTION:
% Plot HRFs for all 3 chromophores and spatial locations.
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

function plotAllHRF_Grp_ChnCriteria_Fig3_TTest(saveOp,opHBCorrected)

sbjList = {'08','12','13','14','15','16','19','21','22','23','24','25'};
fn = sprintf('HbOR_HRF_AllSbjs');
strGrp = 'All';
figSaveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll\Figures\Classification';

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

indHb = zeros(size(t,2),totNumChns_Homer3,length(sbjList));

chnNum = 0;

% HbX x #Channels x #conditions)
numChrom = 3;
numChns = 36;
numConds = 3;
numCoef = 16;
betaGrp = zeros(length(sbjList),numChrom,numChns,numConds);
yAvgGrp = zeros(length(sbjList),numChrom,numChns,numConds);
betaVarGrp = zeros(length(sbjList),numChrom,numChns,numConds);
critPVal = 0.05;

betaIdxRng = 3:10;
numBeta = length(betaIdxRng);

for i = 1:length(sbjList)
    sbjNum = sbjList{i};
    saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
    load([saveDir filesep 'intermediateOutputsCombined_Basis1.mat'],'dcAvg','mlActAuto','beta','bvar');

    %load([saveDir filesep 'intermediateOutputsCombined_Basis1_10Hz.mat'],'dcAvg','mlActAuto','beta','bvar');
    betaVar = beta{1};
    
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
    
    [~, tRangeStart] = min(abs(t-1));
    [~, tRangeEnd] = min(abs(t-6));
    tRangeIdx = [tRangeStart:1:tRangeEnd];
    
    % dcAvg.dataTimeSeries is time x channel
    for i1 = 1:length(srcIdxGrp)
        for j = 1:length(srcIdxGrp{i1})
            for i2 = 1:3
                indHb(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1,i) = ...
                    hrf_GLM(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1);
                indHb(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2,i) = ...
                    hrf_GLM(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2);
                indHb(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3,i) = ...
                    hrf_GLM(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3);
                yAvgGrp(i,1,srcIdxGrp{i1}(j),condSubPlotIdx(i2)) = mean(hrf_GLM(tRangeIdx,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),1);
                yAvgGrp(i,2,srcIdxGrp{i1}(j),condSubPlotIdx(i2)) = mean(hrf_GLM(tRangeIdx,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2),1);
                yAvgGrp(i,3,srcIdxGrp{i1}(j),condSubPlotIdx(i2)) = mean(hrf_GLM(tRangeIdx,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3),1);
            end
        end
    end
                
    % HbX x #Channels x #conditions)
    %betaAvg = squeeze(mean(beta(betaIdxRng,:,:,:),1));
    betaGrp(i,:,:,:) = squeeze(mean(betaVar(betaIdxRng,:,:,:),1));
    
    if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
        bvarIdx = [1 2 3];
    else
        bvarIdx = [1 2 3];
    end
    
    for i2 = 1:numChns
        for i3 = 1:numConds
            for i4 = 1:numChrom
                betaVarGrp(i,i4,i2,i3) = (1/(numBeta^2))*(sum(bvar((bvarIdx(i3)-1)*numCoef+betaIdxRng,i2,i4),1));
            end
        end
    end
    
end

numPairs = numConds*(numConds-1)/2;
pValAct = zeros(numChrom,numChns,numConds);
pValDiff = zeros(numChrom,numChns,numConds);
pValActHBC = zeros(numChrom,numChns,numConds);
pValDiffHBC = zeros(numChrom,numChns,numConds);
sigAct = zeros(numChrom,numChns,numConds);
sigDiff = zeros(numChrom,numChns,numPairs);
sigActBad = zeros(numChrom,numChns,numConds);
sigDiffBad = zeros(numChrom,numChns,numPairs);
m = 3;

listPairs = nchoosek([1;2;3],2);

% beta is #sbj x HbX x #Channels x # conditions
for i1 = 1:numChrom
    for i2 = 1:numChns
        for i3 = 1:numConds
            %[~,sigAct(i1,i2,i3)] = ttest(betaGrp(:,i1,i2,i3));
            [~,pValAct(i1,i2,i3),~,stats{i1,i2,i3}] = ttest(yAvgGrp(:,i1,i2,i3));
            %[~,sigActBad(i1,i2,i3)] = ttest(betaBadGrp(:,i1,i2,i3));
        end
    end
end

if opHBCorrected
for i1 = 1:numChrom
    for i2 = 1:numChns
        tempPVal = zeros(1,3);
        [~,idTemp] = sort(pValAct(i1,i2,:));
        for iK = 1:m
            tempPVal(iK) = min(pValAct(i1,i2,idTemp(1,1,iK))*(m-iK+1),1);
            pValActHBC(i1,i2,idTemp(1,1,iK)) = max(tempPVal);
        end
%         for i3 = 1:numConds
%             sigAct(i1,i2,i3) = pValActHBC(i1,i2,i3)<pCV;
%         end
    end
end
else
    pValActHBC = pValAct;
end

for i1 = 1:numChrom
    for i2 = 1:numChns
        for i3 = 1:numPairs
            %[~,sigDiff(i1,i2,i3)] = ttest(betaGrp(:,i1,i2,listPairs(i3,1)),betaGrp(:,i1,i2,listPairs(i3,2)));
            [~,pValDiff(i1,i2,i3)] = ttest2(yAvgGrp(:,i1,i2,listPairs(i3,1)),yAvgGrp(:,i1,i2,listPairs(i3,2)));
            %[~,sigDiffBad(i1,i2,i3)] = ttest(betaBadGrp(:,i1,i2,listPairs(i3,1)),betaBadGrp(:,i1,i2,listPairs(i3,2)));
        end
    end
end

if opHBCorrected
for i1 = 1:numChrom
    for i2 = 1:numChns
        tempPVal = zeros(1,3);
        [~,idTemp] = sort(pValDiff(i1,i2,:));
        for iK = 1:m
            tempPVal(iK) = min(pValDiff(i1,i2,idTemp(1,1,iK))*(m-iK+1),1);
            pValDiffHBC(i1,i2,idTemp(1,1,iK)) = max(tempPVal);
        end
%         for i3 = 1:numConds
%             sigDiff(i1,i2,i3) = pValDiffHBC(i1,i2,i3)<pCV;
%         end
    end
end
else
    pValDiffHBC = pValDiff;
end

% beta is #sbj x HbX x #Channels x # conditions
% [sigAct,sigDiff] = calcWelchTTest(betaGrp,betaVarGrp);

sigFN = [saveDir filesep 'sig_' strGrp '.mat'];
save(sigFN,'sigAct','sigDiff');

grpHb = squeeze(mean(indHb,3));
%grpHbStd = squeeze(std(indHb,0,3));
grpHbStdErr = std(indHb,0,3)/(size(indHb,3)-1);
% symmetric in this case (not always symmetric)
grpHbCIU = squeeze(grpHbStdErr.*tinv(0.975,size(indHb,3)-1));
grpHbCIL = squeeze(grpHbStdErr.*tinv(0.025,size(indHb,3)-1));

figure();hold on;

condSubPlotIdx = [1 2 3];

leftChnNum = 1;
leftStr = 'Left';
rightChnNum = 18; % S5D9
%rightChnNum = 17; % S5D8
rightStr = 'Right';

% Left HbO
subplot(3,2,1);hold on;
        
xline(0);hold on;
xline(5);
yline(0);
h = zeros(1,3);

for i2 = 1:3

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(leftChnNum-1) + 1;
    inBet = [squeeze(grpHb(:,chnIdx) + grpHbCIU(:,chnIdx)); ...
        flipud(squeeze(grpHb(:,chnIdx) + grpHbCIL(:,chnIdx)))];
    h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');
    
    set(h(i2),'facealpha',0.1);

    if pValActHBC(1,leftChnNum,i2) < critPVal
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    plot(t,grpHb(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
ylim([-8e-5 8e-5]);
%set(gca,'XTick',[], 'YTick', [])

if pValDiffHBC(1,leftChnNum,1) < critPVal
    LRMark = '*';
else
    LRMark = ' ';
end

if pValDiffHBC(1,leftChnNum,2) < critPVal
    LCMark = '+';
else
    LCMark = ' ';
end

if pValDiffHBC(1,leftChnNum,3) < critPVal
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
title(['HbO ' chnName(leftChnNum)]);
subtitle(titleStr);

% Left HbR
subplot(3,2,3);hold on;

xline(0);
xline(5);
yline(0);

for i2 = 1:3
    
    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(leftChnNum-1) + 2;
    
    inBet = [squeeze(grpHb(:,chnIdx)+grpHbCIU(:,chnIdx)); ...
        flipud(squeeze(grpHb(:,chnIdx)+grpHbCIL(:,chnIdx)))];
    h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');hold on;

    set(h(i2),'facealpha',0.1); hold on;
end

h2 = zeros(1,3);

for i2 = 1:3

    if pValActHBC(2,leftChnNum,i2) < critPVal
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(leftChnNum-1) + 2;
    h2(i2) = plot(t,grpHb(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
ylim([-8e-5 8e-5]);
%set(gca,'XTick',[], 'YTick', [])

if pValDiffHBC(2,leftChnNum,1) < critPVal
    LRMark = '*';
else
    LRMark = ' ';
end

if pValDiffHBC(2,leftChnNum,2) < critPVal
    LCMark = '+';
else
    LCMark = ' ';
end

if pValDiffHBC(2,leftChnNum,3) < critPVal
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
title(['HbR ' chnName(leftChnNum)]);
subtitle(titleStr);

% Left HbT
subplot(3,2,5);hold on;

xline(0);
xline(5);
yline(0);

for i2 = 1:3
    
    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(leftChnNum-1) + 3;
    
    inBet = [squeeze(grpHb(:,chnIdx)+grpHbCIU(:,chnIdx)); ...
        flipud(squeeze(grpHb(:,chnIdx)+grpHbCIL(:,chnIdx)))];
    h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');hold on;

    set(h(i2),'facealpha',0.1); hold on;
end

h2 = zeros(1,3);

for i2 = 1:3

    if pValActHBC(3,leftChnNum,i2) < critPVal
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(leftChnNum-1) + 3;
    h2(i2) = plot(t,grpHb(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
ylim([-8e-5 8e-5]);
%set(gca,'XTick',[], 'YTick', [])

if pValDiffHBC(3,leftChnNum,1) < critPVal
    LRMark = '*';
else
    LRMark = ' ';
end

if pValDiffHBC(3,leftChnNum,2) < critPVal
    LCMark = '+';
else
    LCMark = ' ';
end

if pValDiffHBC(3,leftChnNum,3) < critPVal
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
title(['HbT ' chnName(leftChnNum)]);
subtitle(titleStr);

% Right HbO
subplot(3,2,2);hold on;
        
xline(0);hold on;
xline(5);
yline(0);
h = zeros(1,3);

for i2 = 1:3

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(rightChnNum-1) + 1;
    inBet = [squeeze(grpHb(:,chnIdx) + grpHbCIU(:,chnIdx)); ...
        flipud(squeeze(grpHb(:,chnIdx) + grpHbCIL(:,chnIdx)))];
    h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');
    
    set(h(i2),'facealpha',0.1);

    if pValActHBC(1,rightChnNum,i2) < critPVal
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    plot(t,grpHb(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
ylim([-8e-5 8e-5]);
%set(gca,'XTick',[], 'YTick', [])

if pValDiffHBC(1,rightChnNum,1) < critPVal
    LRMark = '*';
else
    LRMark = ' ';
end

if pValDiffHBC(1,rightChnNum,2) < critPVal
    LCMark = '+';
else
    LCMark = ' ';
end

if pValDiffHBC(1,rightChnNum,3) < critPVal
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
title(['HbO ' chnName(rightChnNum)]);
subtitle(titleStr);

% Right HbR
subplot(3,2,4);hold on;

xline(0);
xline(5);
yline(0);

for i2 = 1:3
    
    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(rightChnNum-1) + 2;
    
    inBet = [squeeze(grpHb(:,chnIdx)+grpHbCIU(:,chnIdx)); ...
        flipud(squeeze(grpHb(:,chnIdx)+grpHbCIL(:,chnIdx)))];
    h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');hold on;

    set(h(i2),'facealpha',0.1); hold on;
end

h2 = zeros(1,3);

for i2 = 1:3

    if pValActHBC(2,rightChnNum,i2) < critPVal
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(rightChnNum-1) + 2;
    h2(i2) = plot(t,grpHb(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
ylim([-8e-5 8e-5]);
%set(gca,'XTick',[], 'YTick', [])

if pValDiffHBC(2,rightChnNum,1) < critPVal
    LRMark = '*';
else
    LRMark = ' ';
end

if pValDiffHBC(2,rightChnNum,2) < critPVal
    LCMark = '+';
else
    LCMark = ' ';
end

if pValDiffHBC(2,rightChnNum,3) < critPVal
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
title(['HbR ' chnName(rightChnNum)]);
subtitle(titleStr);

% Right HbT
subplot(3,2,6);hold on;

xline(0);
xline(5);
yline(0);

for i2 = 1:3
    
    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(rightChnNum-1) + 3;
    
    inBet = [squeeze(grpHb(:,chnIdx)+grpHbCIU(:,chnIdx)); ...
        flipud(squeeze(grpHb(:,chnIdx)+grpHbCIL(:,chnIdx)))];
    h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');hold on;

    set(h(i2),'facealpha',0.1); hold on;
end

h2 = zeros(1,3);

for i2 = 1:3

    if pValActHBC(3,rightChnNum,i2) < critPVal
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(rightChnNum-1) + 3;
    h2(i2) = plot(t,grpHb(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
ylim([-8e-5 8e-5]);
%set(gca,'XTick',[], 'YTick', [])

if pValDiffHBC(3,rightChnNum,1) < critPVal
    LRMark = '*';
else
    LRMark = ' ';
end

if pValDiffHBC(3,rightChnNum,2) < critPVal
    LCMark = '+';
else
    LCMark = ' ';
end

if pValDiffHBC(3,rightChnNum,3) < critPVal
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
title(['HbT ' chnName(rightChnNum)]);
subtitle(titleStr);

%legend({'Left Single','Center Single','Right Multiple','Right Single','Left Multiple','Center Multiple'},'Location','northeast');
legend(h2,'Left','Center','Right');
%legend('Left Multiple','Center Multiple','Right Multiple');

sgtitle(sprintf('FEF Grand Average',strGrp));
hold off;

if saveOp
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% OLD CODE %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if sbjGrp==1
%     % all
%     sbjList = {'08','12','13','14','15','16','19','21','22','23','24','25'};
%     %sbjList = {'12','13','14','15','16','19','21','22','23','24'};
%     fn = sprintf('HbOR_HRF_AllSbjs');
%     strGrp = 'All';
% elseif sbjGrp == 2
%     % high-performing
%     sbjList = {'08','12','13','16','19','21','24'};
%     fn = sprintf('HbOR_HRF_HighPerfSbjs');
%     strGrp = 'High-Performing';
% elseif sbjGrp == 0
%     % low-performing
%     sbjList = {'14','15','22','23','25'};
%     fn = sprintf('HbOR_HRF_LowPerfSbjs');
%     strGrp = 'Low-Performing';
% end
% 
% figSaveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll\Figures\Classification';
% 
% fs = 50;
% opPlotInd = 1;
% t = -2:1/fs:15;
% 
% timeSecPostSh = [t fliplr(t)];
% 
% colorIdx = [1 3 2];
% 
% colorList = {[0, 0.447, 0.741],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],...
%     [0.4940, 0.1840, 0.5560],	[0.4660, 0.6740, 0.1880] 	,[0.3010, 0.7450, 0.9330]};
% 
% colorListLight = {[0, 0.447, 0.741, 0.5],[0.8500, 0.3250, 0.0980, 0.5],[0.9290, 0.6940, 0.1250, 0.5],...
%     [0.4940, 0.1840, 0.5560],	[0.4660, 0.6740, 0.1880] 	,[0.3010, 0.7450, 0.9330]};
% 
% 
% 
% lineStyle = {'-' '-' '-'};
% 
% % chnName = {'S1D1 Middle Posterior FEF Left',...
% %     'S1D2 Posterior to sPCS/tgPSC Left',...
% %     'S1D3 Middle FEF Left',...
% %     'S1D4 sPCS/tgPCS Left',...
% %     'S2D1 Posterior FEF Left',...
% %     'S2D2 Inferior to FEF Left',...
% %     'S3D3 Anterior FEF Left',...
% %     'S3D4 iPCS/tgPCS Left',...
% %     'S3D5 Anterior to FEF Left',...
% %     'S3D6 iPCS Left',...
% %     'S4D9 Anterior FEF Right',...
% %     'S4D10 iPCS/tgPCS Right',...
% %     'S4D11 Anterior to FEF Right',...
% %     'S4D12 iPCS Right',...
% %     'S5D7 Middle Posterior FEF Right',...
% %     'S5D8 Posterior to sPCS/tgPSC Right',...
% %     'S5D9 Middle FEF Right',...
% %     'S5D10 sPCS/tgPCS Right',...
% %     'S6D7 Posterior to FEF Right',...
% %     'S6D8 Inferor to FEF Right',...
% %     'S7D16 Posterior STG/PT Left',...
% %     'S8D17 Posterior STG/PT Right',...
% %     'S9D13 IPS3/IPS2/SPL1 Left',...
% %     'S9D14 IPS3/antIPS/IPS4 Left',...
% %     'S10D13 IPS3/IPS2/SPL1 Right',...
% %     'S10D15 IPS3/antIPS/IPS4 Right',...
% %     'S11D13 Superior to IPS3/IPS2/SPL1 Left',...
% %     'S11D14 IPS4 Left',...
% %     'S12D13 Superior to IPS3/IPS2/SPL1 Right',...
% %     'S12D15 IPS4 Right'};
% 
% % chnName = {'FEF',...
% %     'sPCS/tgPSC',...
% %     'FEF',...
% %     'sPCS/tgPCS',...
% %     'FEF',...
% %     'FEF',...
% %     'FEF',...
% %     'iPCS/tgPCS',...
% %     'FEF',...
% %     'iPCS',...
% %     'FEF',...
% %     'iPCS/tgPCS',...
% %     'FEF',...
% %     'iPCS',...
% %     'FEF',...
% %     'sPCS/tgPSC',...
% %     'FEF',...
% %     'sPCS/tgPCS',...
% %     'FEF',...
% %     'FEF',...
% %     'STG/PT',...
% %     'STG/PT',...
% %     'IPS3/IPS2/SPL1',...
% %     'IPS3/antIPS/IPS4',...
% %     'IPS3/IPS2/SPL1',...
% %     'IPS3/antIPS/IPS4',...
% %     'IPS3/IPS2/SPL1',...
% %     'IPS4',...
% %     'IPS3/IPS2/SPL1',...
% %     'IPS4'};
% 
% chnName = {'FEF',...
%     'sPCS/tgPSC',...
%     'FEF/tgPSC/sPCS',...
%     'sPCS/tgPCS',...
%     'FEF',...
%     'tgPSC',...
%     'tgPCS/sPCS/FEF',...
%     'iPCS/tgPCS',...
%     'FEF',...
%     'iPCS/FEF/cIFS/dIPFC',...10
%     'tgPCS/sPCS/FEF',...
%     'iPCS/tgPCS',...
%     'FEF',...
%     'iPCS/FEF/cIFS/dIPFC',...
%     'FEF',...
%     'sPCS/tgPSC',...
%     'FEF/tgPCS/sPCS',...
%     'sPCS/tgPCS',...
%     'FEF',...
%     'tgPCS',...
%     'SMG/STG/PT',...
%     'SMG/STG/PT',...
%     'IPS3/IPS2/SPL1',...
%     'IPS3/antIPS/IPS4',...
%     'IPS3/IPS2/SPL1',...
%     'IPS3/antIPS/IPS4',...
%     'IPS3/IPS2/mSPL',...
%     'IPS4/mSPL',...
%     'IPS3/IPS2/mSPL',...
%     'IPS4/mSPL'};
% 
% srcIdxGrp = {[1 2 3 4] [5 6] [8 9 10 11] [12 13 14 15] [16 17 18 19] [20 21] [23] [25] ...
%     [27 28] [30 31] [33 34] [35 36]};
% 
% chnSubPlotInd = [21 20 15 14 27 26 9 8 3 2 10 11 4 5 22 23 16 17 28 29 31 ...
%     36 45 44 46 47 39 38 40 41];
% 
% numChn = length(chnName);
% origNumConds = 6;
% numConds = 3;
% numConc = 3;
% 
% totNumChns_Homer3 = 36*3*3;
% 
% indHb = zeros(size(t,2),totNumChns_Homer3,length(sbjList));
% %indHbR = zeros(size(t,2),numChn,length(sbjList));
% %indHbT = zeros(size(t,2),length(sbjList));
% 
% chnNum = 0;
% 
% % HbX x #Channels x #conditions)
% numChrom = 3;
% numChns = 36;
% numConds = 3;
% numCoef = 16;
% betaGrp = zeros(length(sbjList),numChrom,numChns,numConds);
% betaVarGrp = zeros(length(sbjList),numChrom,numChns,numConds);
% 
% betaIdxRng = 3:8;
% numBeta = length(betaIdxRng);
% 
% for i = 1:length(sbjList)
%     sbjNum = sbjList{i};
%     saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
%     load([saveDir filesep 'intermediateOutputsCombined_Basis1.mat'],'dcAvg','mlActAuto','beta','bvar');
%     
%     betaVar = beta{1};
%     
%     hrf_GLM = dcAvg.dataTimeSeries;
%     mlList = dcAvg.measurementList;
%     
%     if strcmp(sbjNum,'08') || strcmp(sbjNum,'10')
%         isSDNew = 0;
%         condSubPlotIdx = [2 3 1];
%     else
%         isSDNew = 1;
%         condSubPlotIdx = [1 2 3];
%     end
%     
%     if ~isSDNew
%         % after rmvSingleCond, cut down from 756 chns to 378 chns
%         %[hrf_GLM,mlList] = rmvSingleCond(hrf_GLM,origNumConds,mlList);
%         % after convert2SD2, cut down from 378 chns to 324 chns
%         [hrf_GLM,mlList,mlActAuto] = convert2SD2_cnt(hrf_GLM,mlList,mlActAuto{1});
%         %betaVar(:,:,:,[1 2 3]) = [];
%         betaVar(:,:,[29,33,39,40,41,42],:) = [];
%     else
%         mlActAuto = mlActAuto{1};
%     end
%     
%     sizeCond = size(mlList,2)/numConds;
%     numChns = size(mlList,2)/(numConds*numConc);
%     sizeChn = sizeCond/numChns;
%     
%     % dcAvg.dataTimeSeries is time x channel
%     for i1 = 1:length(srcIdxGrp)
%         for j = 1:length(srcIdxGrp{i1})
%             for i2 = 1:3
%                 indHb(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1,i) = ...
%                     hrf_GLM(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1);
%                 indHb(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2,i) = ...
%                     hrf_GLM(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2);
%                 indHb(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3,i) = ...
%                     hrf_GLM(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3);
%             end
%         end
%     end
%                 
%     % HbX x #Channels x #conditions)
%     %betaAvg = squeeze(mean(beta(betaIdxRng,:,:,:),1));
%     betaGrp(i,:,:,:) = squeeze(mean(betaVar(betaIdxRng,:,:,:),1));
%     
%     if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
%         bvarIdx = [1 2 3];
%     else
%         bvarIdx = [1 2 3];
%     end
%     
%     for i2 = 1:numChns
%         for i3 = 1:numConds
%             for i4 = 1:numChrom
%                 betaVarGrp(i,i4,i2,i3) = (1/(numBeta^2))*(sum(bvar((bvarIdx(i3)-1)*numCoef+betaIdxRng,i2,i4),1));
%             end
%         end
%     end
%     
% end
% 
% % beta is #sbj x HbX x #Channels x # conditions
% [sigAct,sigDiff] = calcWelchTTest(betaGrp,betaVarGrp);
% 
% sigFN = [saveDir filesep 'sig_Bad.mat'];
% save(sigFN,'sigAct','sigDiff');
% 
% grpHb = squeeze(mean(indHb,3));
% %grpHbStd = squeeze(std(indHb,0,3));
% grpHbStdErr = std(indHb,0,3)/(size(indHb,3)-1);
% % symmetric in this case (not always symmetric)
% grpHbCIU = squeeze(grpHbStdErr.*tinv(0.975,size(indHb,3)-1));
% grpHbCIL = squeeze(grpHbStdErr.*tinv(0.025,size(indHb,3)-1));
% 
% % grpHbWeighted = sum(w.*indHb,3)./sum(w.3);
% % grpHbStdWeighted = squeeze(std(indHb,0,3));
% 
% figure('units','normalized','outerposition',[0 0 1 1]); hold on;
% 
% condSubPlotIdx = [1 2 3];
% 
% % now plot
% if criteriaOp == 3
%     chnNum = 5;
%     fn = [fn 'LFEF'];
%     strCrit = 'Left';
% elseif criteriaOp == 4
%     chnNum = 18;
%     fn = [fn 'RFEF'];
%     strCrit = 'Right';
% end
% 
% subplot(2,1,1);hold on;
%         
% xline(0);hold on;
% xline(5);
% yline(0);
% h = zeros(1,3);
% 
% for i2 = 1:3
% 
%     chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(chnNum-1) + 1;
% 
% %             inBet = [squeeze(grpHb(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1)...
% %                 +grpHbStd(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1)); ...
% %                 fliplr(squeeze(grpHb(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1)...
% %                 -grpHbStd(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1)))];
%     inBet = [squeeze(grpHb(:,chnIdx) + grpHbCIU(:,chnIdx)); ...
%         flipud(squeeze(grpHb(:,chnIdx) + grpHbCIL(:,chnIdx)))];
%     h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');
%     
%     set(h(i2),'facealpha',0.1);
% 
%     if sigAct(1,chnNum,i2)
%         lineWidth = 2;
%     else
%         lineWidth = 0.5;
%     end
% 
%     plot(t,grpHb(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;
% 
% end
% 
% xlim([-2 15]);
% %set(gca,'XTick',[], 'YTick', [])
% 
% if sigDiff(1,chnNum,1)
%     LRMark = '*';
% else
%     LRMark = ' ';
% end
% 
% if sigDiff(1,chnNum,2)
%     LCMark = '+';
% else
%     LCMark = ' ';
% end
% 
% if sigDiff(1,chnNum,3)
%     RCMark = 'x';
% else
%     RCMark = ' ';
% end
% 
% titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
% title(['HbO ' chnName(chnNum)]);
% subtitle(titleStr);
% 
% %legend('Left Multiple','Center Multiple','Right Multiple','Location','Southeast');
%         
% subplot(2,1,2);hold on;
% 
% xline(0);
% xline(5);
% yline(0);
% 
% for i2 = 1:3
%     
%     chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(chnNum-1) + 2;
%     
%     inBet = [squeeze(grpHb(:,chnIdx)+grpHbCIU(:,chnIdx)); ...
%         flipud(squeeze(grpHb(:,chnIdx)+grpHbCIL(:,chnIdx)))];
%     h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');hold on;
% 
%     set(h(i2),'facealpha',0.1); hold on;
% end
% 
% h2 = zeros(1,3);
% 
% for i2 = 1:3
% 
%     if sigAct(2,chnNum,i2)
%         lineWidth = 2;
%     else
%         lineWidth = 0.5;
%     end
% 
%     chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(chnNum-1) + 2;
%     h2(i2) = plot(t,grpHb(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;
% 
% end
% 
% xlim([-2 15]);
% %set(gca,'XTick',[], 'YTick', [])
% 
% if sigDiff(2,chnNum,1)
%     LRMark = '*';
% else
%     LRMark = ' ';
% end
% 
% if sigDiff(2,chnNum,2)
%     LCMark = '+';
% else
%     LCMark = ' ';
% end
% 
% if sigDiff(2,chnNum,3)
%     RCMark = 'x';
% else
%     RCMark = ' ';
% end
% 
% titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
% title(['HbR ' chnName(chnNum)]);
% subtitle(titleStr);
% 
% 
% %legend({'Left Single','Center Single','Right Multiple','Right Single','Left Multiple','Center Multiple'},'Location','northeast');
% legend(h2,'Left Multiple','Center Multiple','Right Multiple');
% %legend('Left Multiple','Center Multiple','Right Multiple');

