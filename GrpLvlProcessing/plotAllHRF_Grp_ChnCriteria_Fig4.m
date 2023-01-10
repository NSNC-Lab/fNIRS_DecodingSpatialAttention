% Arithmetic average of HRFs
% Light line color, plot ind HRFs
% Different processing for sbj 08 and 10 because they use old SD design and
% has both single and multi conditions
% new codes, ignore grpPlotHRF

function plotAllHRF_Grp_ChnCriteria_Fig4(saveOp)

goodSbjList = {'08','12','13','16','19','21','24'};
badSbjList = {'14','15','22','23','25'};
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

indHbGood = zeros(size(t,2),totNumChns_Homer3,length(goodSbjList));
indHbBad = zeros(size(t,2),totNumChns_Homer3,length(badSbjList));

chnNum = 0;

% HbX x #Channels x #conditions)
numChrom = 3;
numChns = 36;
numConds = 3;
numCoef = 16;
betaGoodGrp = zeros(length(goodSbjList),numChrom,numChns,numConds);
betaVarGoodGrp = zeros(length(goodSbjList),numChrom,numChns,numConds);
betaBadGrp = zeros(length(badSbjList),numChrom,numChns,numConds);
betaVarBadGrp = zeros(length(badSbjList),numChrom,numChns,numConds);

betaIdxRng = 3:10;
numBeta = length(betaIdxRng);

for i = 1:length(goodSbjList)
    sbjNum = goodSbjList{i};
    saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
    load([saveDir filesep 'intermediateOutputsCombined_Basis1.mat'],'dcAvg','mlActAuto','beta','bvar');
    
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

% bad sbj list
for i = 1:length(badSbjList)
    sbjNum = badSbjList{i};
    saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
    load([saveDir filesep 'intermediateOutputsCombined_Basis1.mat'],'dcAvg','mlActAuto','beta','bvar');
    
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

% beta is #sbj x HbX x #Channels x # conditions
% [h,p,ci,stats] = ttest(x);
% [h,p,ci,stats] = ttest(x);
[sigActGood,sigDiffGood] = calcWelchTTest(betaGoodGrp,betaVarGoodGrp);
[sigActBad,sigDiffBad] = calcWelchTTest(betaBadGrp,betaVarBadGrp);

% sigFN = [saveDir filesep 'sig_' strGrp '.mat'];
% save(sigFN,'sigActGood','sigDiffGood','sigActBad','sigDiffBad');

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

leftChnNum = 1;
leftStr = 'Left';
rightChnNum = 18;
rightStr = 'Right';

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

    if sigActGood(1,leftChnNum,i2)
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    plot(t,grpHbGood(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if sigDiffGood(1,leftChnNum,1)
    LRMark = '*';
else
    LRMark = ' ';
end

if sigDiffGood(1,leftChnNum,2)
    LCMark = '+';
else
    LCMark = ' ';
end

if sigDiffGood(1,leftChnNum,3)
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

    if sigActBad(1,leftChnNum,i2)
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    plot(t,grpHbBad(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if sigDiffBad(1,leftChnNum,1)
    LRMark = '*';
else
    LRMark = ' ';
end

if sigDiffBad(1,leftChnNum,2)
    LCMark = '+';
else
    LCMark = ' ';
end

if sigDiffBad(1,leftChnNum,3)
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

    if sigActGood(2,leftChnNum,i2)
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(leftChnNum-1) + 2;
    h2(i2) = plot(t,grpHbGood(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if sigDiffGood(2,leftChnNum,1)
    LRMark = '*';
else
    LRMark = ' ';
end

if sigDiffGood(2,leftChnNum,2)
    LCMark = '+';
else
    LCMark = ' ';
end

if sigDiffGood(2,leftChnNum,3)
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

    if sigActBad(2,leftChnNum,i2)
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(leftChnNum-1) + 2;
    h2(i2) = plot(t,grpHbBad(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if sigDiffBad(2,leftChnNum,1)
    LRMark = '*';
else
    LRMark = ' ';
end

if sigDiffBad(2,leftChnNum,2)
    LCMark = '+';
else
    LCMark = ' ';
end

if sigDiffBad(2,leftChnNum,3)
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

    if sigActGood(3,leftChnNum,i2)
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(leftChnNum-1) + 3;
    h2(i2) = plot(t,grpHbGood(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if sigDiffGood(3,leftChnNum,1)
    LRMark = '*';
else
    LRMark = ' ';
end

if sigDiffGood(3,leftChnNum,2)
    LCMark = '+';
else
    LCMark = ' ';
end

if sigDiffGood(3,leftChnNum,3)
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
%title(['HbR ' chnName(leftChnNum)]);
title(['HbR Left High']);
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

    if sigActBad(3,leftChnNum,i2)
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(leftChnNum-1) + 3;
    h2(i2) = plot(t,grpHbBad(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if sigDiffBad(3,leftChnNum,1)
    LRMark = '*';
else
    LRMark = ' ';
end

if sigDiffBad(3,leftChnNum,2)
    LCMark = '+';
else
    LCMark = ' ';
end

if sigDiffBad(3,leftChnNum,3)
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
%title(['HbR ' chnName(leftChnNum)]);
title(['HbO Left Low']);
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

    if sigActGood(1,rightChnNum,i2)
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    plot(t,grpHbGood(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if sigDiffGood(1,rightChnNum,1)
    LRMark = '*';
else
    LRMark = ' ';
end

if sigDiffGood(1,rightChnNum,2)
    LCMark = '+';
else
    LCMark = ' ';
end

if sigDiffGood(1,rightChnNum,3)
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

    if sigActGood(1,rightChnNum,i2)
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    plot(t,grpHbBad(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if sigDiffBad(1,rightChnNum,1)
    LRMark = '*';
else
    LRMark = ' ';
end

if sigDiffBad(1,rightChnNum,2)
    LCMark = '+';
else
    LCMark = ' ';
end

if sigDiffBad(1,rightChnNum,3)
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

    if sigActGood(2,rightChnNum,i2)
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(rightChnNum-1) + 2;
    h2(i2) = plot(t,grpHbGood(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if sigDiffGood(2,rightChnNum,1)
    LRMark = '*';
else
    LRMark = ' ';
end

if sigDiffGood(2,rightChnNum,2)
    LCMark = '+';
else
    LCMark = ' ';
end

if sigDiffGood(2,rightChnNum,3)
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

    if sigActBad(2,rightChnNum,i2)
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(rightChnNum-1) + 2;
    h2(i2) = plot(t,grpHbBad(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if sigDiffBad(2,rightChnNum,1)
    LRMark = '*';
else
    LRMark = ' ';
end

if sigDiffBad(2,rightChnNum,2)
    LCMark = '+';
else
    LCMark = ' ';
end

if sigDiffBad(2,rightChnNum,3)
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

    if sigActGood(3,rightChnNum,i2)
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(rightChnNum-1) + 3;
    h2(i2) = plot(t,grpHbGood(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if sigDiffGood(3,rightChnNum,1)
    LRMark = '*';
else
    LRMark = ' ';
end

if sigDiffGood(3,rightChnNum,2)
    LCMark = '+';
else
    LCMark = ' ';
end

if sigDiffGood(3,rightChnNum,3)
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
%title(['HbR ' chnName(rightChnNum)]);
title(['HbR Right High']);
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

    if sigActBad(3,rightChnNum,i2)
        lineWidth = 2;
    else
        lineWidth = 0.5;
    end

    chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(rightChnNum-1) + 3;
    h2(i2) = plot(t,grpHbBad(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;

end

xlim([-2 15]);
%set(gca,'XTick',[], 'YTick', [])

if sigDiffBad(3,rightChnNum,1)
    LRMark = '*';
else
    LRMark = ' ';
end

if sigDiffBad(3,rightChnNum,2)
    LCMark = '+';
else
    LCMark = ' ';
end

if sigDiffBad(3,rightChnNum,3)
    RCMark = 'x';
else
    RCMark = ' ';
end

titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
%title(['HbR ' chnName(rightChnNum)]);
title(['HbR Right Low']);
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