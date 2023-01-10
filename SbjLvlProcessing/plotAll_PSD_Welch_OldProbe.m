% For sbj 08 and 10
function plotAll_PSD_Welch_OldProbe(sbjNum)

figSaveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum) '\Figures\HRFs'];
saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
rawDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
load([saveDir filesep 'intermediateOutputsCombined_Basis1.mat']);
load([saveDir filesep 'performanceLinearDiscriminantUpdated_LR.mat']);
load([rawDir filesep sbjNum '.mat'],'s');

numConds = 3;
numConc = 3;
numChns = size(dcAvg.measurementList,2)/(numConds*numConc);
sizeCond = size(dcAvg.measurementList,2)/numConds;
sizeChn = sizeCond/numChns;
%condSubPlotIdx = [1 3 5 2 4 6];
condSubPlotIdx = [1 2 3];

colorIdx = [1 3 2 2 1 3];

lineStyle = {'-' '-' '-' '-' '--' '--'};

if strcmp(sbjNum,'07')||strcmp(sbjNum,'08') || strcmp(sbjNum,'10')
    srcIdxGrp = {[1 2 3 4] [5 6] [8 9 10 11] [12 13 14 15] [16 17 18 19] [20 21] [23] [25] ...
        [27 28 29] [31 32 33] [35 36] [37 38] [40] [42]};
else
    srcIdxGrp = {[1 2 3 4] [5 6] [8 9 10 11] [12 13 14 15] [16 17 18 19] [20 21] [23] [25] ...
        [27 28] [30 31] [33 34] [35 36]};
end

fs = 50;
timePt = (0:0.5*fs:12*fs)+2*fs;

timeIdx = find(timePt==(5*fs+2*fs));

% number of channels, can pull
if strcmp(sbjNum,'07')|| strcmp(sbjNum,'08') ||strcmp(sbjNum,'10')
    mlActAuto_Unused = getPrunedChns_Special(sbjNum,s.fName,s.movieList,2);
else
    mlActAuto_Unused = getPrunedChns(sbjNum,s.fName,s.movieList,2);
end

ml690 = mlActAuto_Unused{1}(1:length(mlActAuto_Unused{1})/2);
ml870 = mlActAuto_Unused{1}(length(mlActAuto_Unused{1})/2+1:length(mlActAuto_Unused{1}));
acceptedChns = (ml690==ml870);
%numUnusedChns = length(mlActAuto_Unused{1})/2-sum(ml690==ml870);

% subplot is 8x6.
chnName = {'S1D1 Middle Posterior FEF Left',...
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
    'S12D15	IPS4 Right'};

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

if strcmp(sbjNum,'07')||strcmp(sbjNum,'08') || strcmp(sbjNum,'10')
    chnSubPlotInd = [21 20 15 14 27 26 9 8 3 2 10 11 4 5 22 23 16 17 28 29 31 ...
        36 45 44 50 46 47 53 39 38 40 41 7 12];
else
    chnSubPlotInd = [21 20 15 14 27 26 9 8 3 2 10 11 4 5 22 23 16 17 28 29 31 ...
        36 45 44 46 47 39 38 40 41];
end

colorList = {[0, 0.447, 0.741],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],...
    [0.4940, 0.1840, 0.5560],	[0.4660, 0.6740, 0.1880] 	,[0.3010, 0.7450, 0.9330]};

maxPlot = zeros(1,length(srcIdxGrp));
minPlot = zeros(1,length(srcIdxGrp));

% pxx is freq x chns. w is freq x 1.
%[pxx,w] = calcWelch(s.name,s.fName,s.movieList);
[pxx,w] = calcWelch_OldProbe(s.name,s.fName,s.movieList);

for i1 = 1:length(srcIdxGrp)
    tempMax = zeros(1,3);
    tempMin = zeros(1,3);
    for i2 = 1:3
        tempMax(i2) = max(max(pxx(:,(srcIdxGrp{i1}(:)-1) + 1)));
        tempMin(i2) = min(min(pxx(:,(srcIdxGrp{i1}(:)-1) + 1)));
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

figure('units','normalized','outerposition',[0 0 1 1]); hold on;
for i1 = 1:length(srcIdxGrp)
    for j = 1:length(srcIdxGrp{i1})
        chnNum = chnNum + 1;
        
        subplot(9,6,chnSubPlotInd(chnNum));hold on;

            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),'-','Color',colorList{j-0});hold on;
        plot(w/pi, 10*log10(pxx(:,(srcIdxGrp{i1}(j)))),lineStyle{1},'Color',colorList{colorIdx(1)});hold on;
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2),':','Color',colorList{j-0});hold on;
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3),'-.','Color',colorList{j-0});hold on;

        %legend(chnNameGrp{i1},'Location','northeastoutside');
        %legend({'Left Single','Center Single','Right Multiple','Right Single','Left Multiple','Center Multiple'},'Location','northeast');
%         srcStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).sourceIndex);
%         detectorStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).detectorIndex);
%         ylabel('[M*mm]')
%         xlabel('Time [s]');
        if ~acceptedChns(chnNum)
            ax = gca;
            ax.XColor = 'red';
            ax.YColor = 'red';
        end
        %label = num2str(dcavg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).dataTypeLabel);
        %condStr = snirf1.stim(1,i2).name;
        %title(sprintf('%s',num2str(performanceArrSingleHbO(srcIdxGrp{i1}(j),timeIdx))));
%         title(sprintf('S%sD%s R: %s CAS: %s CAM: %s',...
%             srcStr,detectorStr,num2str(R(srcIdxGrp{i1}(j),1)),...
%             num2str(performanceArrSingleHbO(srcIdxGrp{i1}(j),timeIdx)),...
%             num2str(performanceArrMultiHbO(srcIdxGrp{i1}(j),timeIdx))));
        %subtitle(sprintf('HbO %s', chnName{chnNum}));
        
    end
end

subplot(9,6,54);hold on;
plot(0,0);
plot(0,0);
plot(0,0);
legend('Left Multiple','Center Multiple','Right Multiple');

sgtitle(sprintf('Sbj %s 690 nm',sbjNum));
hold off;

if ~exist(figSaveDir,'dir')
    mkdir(figSaveDir);
end

fn = sprintf('PSD_690nm');
print(gcf,[figSaveDir filesep fn],'-dpng','-r250');

chnNum = 0;

%%%%%%%%%%%%
figure('units','normalized','outerposition',[0 0 1 1]); hold on;
for i1 = 1:length(srcIdxGrp)
    for j = 1:length(srcIdxGrp{i1})
        chnNum = chnNum + 1;

        subplot(9,6,chnSubPlotInd(chnNum));hold on;

            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),'-','Color',colorList{j-0});hold on;
        plot(w/pi, 10*log10(pxx(:,(srcIdxGrp{i1}(j)+42))),lineStyle{1},'Color',colorList{colorIdx(1)});hold on;
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2),':','Color',colorList{j-0});hold on;
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3),'-.','Color',colorList{j-0});hold on;

        %legend(chnNameGrp{i1},'Location','northeastoutside');
%         legend({'Left Single','Center Single','Right Multiple','Right Single','Left Multiple','Center Multiple'},'Location','northeast');
%         srcStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).sourceIndex);
%         detectorStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).detectorIndex);
%         ylabel('[M*mm]')
%         xlabel('Time [s]');

        %label = num2str(dcavg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).dataTypeLabel);
        %condStr = snirf1.stim(1,i2).name;
%         title(sprintf('S%sD%s R: %s CAS: %s CAM: %s',...
%             srcStr,detectorStr,num2str(R(srcIdxGrp{i1}(j),2)),...
%             num2str(performanceArrSingleHbR(srcIdxGrp{i1}(j),timeIdx)),...
%             num2str(performanceArrMultiHbR(srcIdxGrp{i1}(j),timeIdx))));
%         subtitle(sprintf('HbR %s',chnName{chnNum}));
        
        %title(sprintf('%s',condStr));
        if ~acceptedChns(chnNum)
            ax = gca;
            ax.XColor = 'red';
            ax.YColor = 'red';
        end
        hold off;
%         fn = sprintf('HRF_S%sD%s',srcStr,detectorStr);
%         print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
%         close;
    end
end
sgtitle(sprintf('Sbj %s 830 nm',sbjNum));
subplot(9,6,54);hold on;
plot(0,0);
plot(0,0);
plot(0,0);
legend('Left Multiple','Center Multiple','Right Multiple');

%sgtitle(sprintf('Sbj %s HbO',sbjNum));
hold off;
fn = sprintf('PSD_830nm');
print(gcf,[figSaveDir filesep fn],'-dpng','-r250');

end