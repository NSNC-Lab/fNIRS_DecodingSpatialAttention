% fNIRS version
function plotBands_fNIRS(sbjNum,saveOp)

numClasses = 3;
exportDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum) '\Figures\HRFs'];
saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
dir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];

load([saveDir filesep 'intermediateOutputsCombined_Basis1.mat']);
load([dir filesep sbjNum '.mat'],'s');

cLim = 0.5*10^-4;
%cLim = 2*10^-4;

numConds = 3;
numConc = 3;
numChns = size(dcAvg.measurementList,2)/(numConds*numConc);
sizeCond = size(dcAvg.measurementList,2)/numConds;
sizeChn = sizeCond/numChns;
%condSubPlotIdx = [1 3 5 2 4 6];
condSubPlotIdx = [1 2 3];

lineStyle = {'-' '-' '-'};

if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
    srcIdxGrp = {[1 2 3 4] [5 6] [8 9 10 11] [12 13 14 15] [16 17 18 19] [20 21] [23] [25] ...
    [27 28 29] [31 32 33] [35 36] [37 38] [40] [42]};
else
    srcIdxGrp = {[1 2 3 4] [5 6] [8 9 10 11] [12 13 14 15] [16 17 18 19] [20 21] [23] [25] ...
        [27 28] [30 31] [33 34] [35 36]};
end

% number of channels, can pull
if strcmp(sbjNum,'08') ||strcmp(sbjNum,'10')
    mlActAuto_Unused = getPrunedChns_Special(sbjNum,s.fName,s.movieList,2);
    
    chnName = categorical({'S1D1 Middle Posterior FEF Left',...
    'S1D2 Posterior to sPCS/tgPSC Left',...
    'S1D3 Middle FEF Left',...
    'S1D4 sPCS/tgPCS Left',...
    'S2D1 Posterior FEF Left',...
    'S2D2 Inferior to FEF Left',...
    'S3D3 Anterior FEF Left',...
    'S3D4 iPCS/tgPCS Left',...
    'S3D5 Anterior to FEF Left',...
    'S3D6 iPCS Left',...
    'S4D9 Anterior FEF Right',...
    'S4D10 iPCS/tgPCS Right',...
    'S4D11 Anterior to FEF Right',...
    'S4D12 iPCS Right',...
    'S5D7 Middle Posterior FEF Right',...
    'S5D8 Posterior to sPCS/tgPSC Right',...
    'S5D9 Middle FEF Right',...
    'S5D10 sPCS/tgPCS Right',...
    'S6D7 Posterior to FEF Right',...
    'S6D8 Inferor to FEF Right',...
    'S7D18 Posterior STG/PT Left',...
    'S8D19 Posterior STG/PT Right',...
    'S9D13 IPS3/IPS2/SPL1 Left',...
    'S9D14 IPS3/antIPS/IPS4 Left',...
    'S9D16 IPS3/latIPS/antIPS Left',...
    'S10D13 IPS3/IPS2/SPL1 Right',...
    'S10D15	IPS3/antIPS/IPS4 Right',...
    'S10D17 IPS3/latIPS/antIPS Right',...
    'S11D13	Superior to IPS3/IPS2/SPL1 Left',...
    'S11D14	IPS4 Left',...
    'S12D13	Superior to IPS3/IPS2/SPL1 Right',...
    'S12D15	IPS4 Right',...
    'S13D28 Ant Posterior STG/PT Left',...
    'S14D29 Ant Posterior STG/PT Right'});

    reorderIdx = [33 21 10 8 4 2 6 9 7 3 1 5 30 24 25 29 23 31 26 13 11 17 ...
        15 19 32 27 28 14 12 18 16 20 22 34];
else
    mlActAuto_Unused = getPrunedChns(sbjNum,s.fName,s.movieList,2);
    
    chnName = categorical({'S1D1 Middle Posterior FEF Left',...
    'S1D2 Posterior to sPCS/tgPSC Left',...
    'S1D3 Middle FEF Left',...
    'S1D4 sPCS/tgPCS Left',...
    'S2D1 Posterior FEF Left',...
    'S2D2 Inferior to FEF Left',...
    'S3D3 Anterior FEF Left',...
    'S3D4 iPCS/tgPCS Left',...
    'S3D5 Anterior to FEF Left',...
    'S3D6 iPCS Left',...
    'S4D9 Anterior FEF Right',...
    'S4D10 iPCS/tgPCS Right',...
    'S4D11 Anterior to FEF Right',...
    'S4D12 iPCS Right',...
    'S5D7 Middle Posterior FEF Right',...
    'S5D8 Posterior to sPCS/tgPSC Right',...
    'S5D9 Middle FEF Right',...
    'S5D10 sPCS/tgPCS Right',...
    'S6D7 Posterior to FEF Right',...
    'S6D8 Inferor to FEF Right',...
    'S7D16 Posterior STG/PT Left',...
    'S8D17 Posterior STG/PT Right',...
    'S9D13 IPS3/IPS2/SPL1 Left',...
    'S9D14 IPS3/antIPS/IPS4 Left',...
    'S10D13 IPS3/IPS2/SPL1 Right',...
    'S10D15	IPS3/antIPS/IPS4 Right',...
    'S11D13	Superior to IPS3/IPS2/SPL1 Left',...
    'S11D14	IPS4 Left',...
    'S12D13	Superior to IPS3/IPS2/SPL1 Right',...
    'S12D15	IPS4 Right'});

    reorderIdx = [21 10 8 4 2 6 9 7 3 1 5 28 24 27 23 29 25 13 11 17 15 19 30 ...
        26 14 12 18 16 20 22];
end

ml690 = mlActAuto_Unused{1}(1:length(mlActAuto_Unused{1})/2);
ml870 = mlActAuto_Unused{1}(length(mlActAuto_Unused{1})/2+1:length(mlActAuto_Unused{1}));
acceptedChns = (ml690==ml870);

fs = 50;
t = -2:1/fs:12;

maxPlot = zeros(1,length(srcIdxGrp));
minPlot = zeros(1,length(srcIdxGrp));

channelIdxLeft = sizeCond*(condSubPlotIdx(1)-1)+sizeChn*(horzcat(srcIdxGrp{:})-1)+1;
channelIdxRight = sizeCond*(condSubPlotIdx(2)-1)+sizeChn*(horzcat(srcIdxGrp{:})-1)+1;
channelIdxCenter = sizeCond*(condSubPlotIdx(3)-1)+sizeChn*(horzcat(srcIdxGrp{:})-1)+1;

leftTrials = dcAvg.dataTimeSeries(:,channelIdxLeft);
rightTrials = dcAvg.dataTimeSeries(:,channelIdxRight);
centerTrials = dcAvg.dataTimeSeries(:,channelIdxCenter);

for i1 = 1:length(srcIdxGrp)
    tempMax = zeros(1,3);
    tempMin = zeros(1,3);
    for i2 = 1:3
        tempMax(i2) = max(max(dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(:)-1) + 1)));
        tempMin(i2) = min(min(dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(:)-1) + 1)));
    end
    if max(tempMax) ~= 0
        maxPlot(i1) = max(tempMax);
    else
        maxPlot(i1) = 1e-4;
    end
    if min(tempMin) ~= 0
        minPlot(i1) = min(tempMin);
    else
        minPlot(i1) = -1e-4;
    end
end

R = R{1};

chnNum = 0;

chnName = chnName(reorderIdx);
leftTrials = leftTrials(:,reorderIdx);
rightTrials = rightTrials(:,reorderIdx);
centerTrials = centerTrials(:,reorderIdx);

% # trials x chn x time

figure('units','normalized','outerposition',[0 0 1 1]);subplot(1,3,1);
y = 1:length(chnName);
%imAlpha=ones(size(leftSingleAvgPreSorted));
%imAlpha(isnan(leftSingleAvgPreSorted))=0;
%imagesc(t,y,leftSingleAvgPreSorted,'AlphaData',imAlpha);
imagesc(t,y,leftTrials');
hold on;
%line([minDiffPrePart1*fsRatio, minDiffPrePart1*fsRatio],[1,64]);
set(gca,'YTick',y,'YTickLabel',chnName);
title('Left [Part 1]');
xlabel('time [ms]');
ylabel('\muV');
caxis manual
%caxis([lowerLim upperLim]);
caxis([-cLim cLim]);
colorbar;
hold off;

% Right Part 1
%figure(3);
subplot(1,3,3);
%figure('units','normalized','outerposition',[0 0 1 1]);
%t = 1:fsRatio:size(avgEEGRightPart1Sorted,2)*fsRatio;
%imAlpha=ones(size(rightSingleAvgPreSorted));
%imAlpha(isnan(rightSingleAvgPreSorted))=0;
%imagesc(t,y,rightSingleAvgPreSorted,'AlphaData',imAlpha);
imagesc(t,y,rightTrials');
hold on;
%line([minDiffPrePart1*fsRatio, minDiffPrePart1*fsRatio],[1,64]);
%set(gca,'YTick',y,'YTickLabel',chnName);
title('Right [Part 1]');
xlabel('time [ms]');
ylabel('\muV');
caxis manual
%caxis([lowerLim upperLim]);
caxis([-cLim cLim]);
colorbar;
hold off;

% Center Part 1
%figure(5);
subplot(1,3,2);
%figure('units','normalized','outerposition',[0 0 1 1]);
%t = 1:fsRatio:size(avgEEGCenterPart1Sorted,2)*fsRatio;
% imAlpha=ones(size(centerSingleAvgPreSorted));
% imAlpha(isnan(centerSingleAvgPreSorted))=0;
% imagesc(t,y,centerSingleAvgPreSorted,'AlphaData',imAlpha);
imagesc(t,y,centerTrials');
hold on;
%line([minDiffPrePart1*fsRatio, minDiffPrePart1*fsRatio],[1,64]);
%set(gca,'YTick',y,'YTickLabel',chnName);
title('Center [Part 1]');
xlabel('time [ms]');
ylabel('\muV');
caxis manual
%caxis([lowerLim upperLim]);
caxis([-cLim cLim]);
colorbar;
hold off;


sgtitle(sprintf('%s: Multiple Movie',s.name));

if saveOp == 1
    filename = [exportDir filesep 'Band_MultipleMovie.png'];
    print(gcf,filename,'-dpng','-r350');
end
