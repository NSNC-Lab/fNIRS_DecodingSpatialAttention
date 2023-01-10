% Plot HRF using outputs from scripting
% plotHRF.m plot HRF using outputs from group.mat from Homer3 GUI
% Completed
% This plot all channels from one source, divide into 6 subplots, one for
% each condition
function plotHRF3(sbjNum)

figSaveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum) '\Figures\HRFs'];
saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
load([saveDir filesep 'intermediateOutputsCombined.mat']);

numConds = 6;
numConc = 3;
numChns = size(dcAvg.measurementList,2)/(numConds*numConc);
sizeCond = size(dcAvg.measurementList,2)/numConds;
sizeChn = sizeCond/numChns;

condSubPlotIdx = [1 3 5 2 4 6];

colorIdx = [1 3 2 2 1 3];

lineStyle = {'-' '-' ':' '-' ':' ':'};

srcIdxGrp = {[1 2 3 4] [5 6] [8 9 10 11] [12 13 14 15] [16 17 18 19] [20 21] [23] [25] ...
    [27 28 29] [31 32 33] [35 36] [37 38] [40] [42]};


% chnNameGrp = {{'S1D1 HbO','S1D1 HbR','S1D1 HbT','S1D2 HbO','S1D2 HbR','S1D2 HbT','S1D3 HbO','S1D3 HbR','S1D3 HbT','S1D4 HbO','S1D4 HbR','S1D4 HbT'},...
%     {'S2D1 HbO','S2D1 HbR','S2D1 HbT','S2D2 HbO','S2D2 HbR','S2D2 HbT'},...
%     {'S3D3 HbO','S3D3 HbR','S3D3 HbT','S3D4 HbO','S3D4 HbR','S3D4 HbT','S3D5 HbO','S3D5 HbR','S3D5 HbT','S3D6 HbO','S3D6 HbR','S3D6 HbT'},...
%     {'S4D9 HbO','S4D9 HbR','S4D9 HbT','S4D10 HbO','S4D10 HbR','S4D10 HbT','S4D11 HbO','S4D11 HbR','S4D11 HbT','S4D12 HbO','S4D12 HbR','S4D12 HbT'},...
%     {'S5D7 HbO','S5D7 HbR','S5D7 HbT','S5D8 HbO','S5D8 HbR','S5D8 HbT','S5D9 HbO','S5D9 HbR','S5D9 HbT','S5D10 HbO','S5D10 HbR','S5D10 HbT'},...
%     {'S6D7 HbO','S6D7 HbR','S6D7 HbT','S6D8 HbO','S6D8 HbR','S6D8 HbT'},...
%     {'S7D18 HbO','S7D18 HbR','S7D18 HbT'},...
%     {'S8D19 HbO','S8D19 HbR','S8D19 HbT'},...
%     {'S9D13 HbO','S9D13 HbR','S9D13 HbT','S9D14 HbO','S9D14 HbR','S9D14 HbT','S9D16 HbO','S9D16 HbR','S9D16 HbT'},...
%     {'S10D13 HbO','S10D13 HbR','S10D13 HbT','S10D15 HbO','S10D15 HbR','S10D15 HbT','S10D17 HbO','S10D17 HbR','S10D17 HbT'},...
%     {'S11D13 HbO','S11D13 HbR','S11D13 HbT','S11D14 HbO','S11D14 HbR','S11D14 HbT'},...
%     {'S12D13 HbO','S12D13 HbR','S12D13 HbT','S12D15 HbO','S12D15 HbR','S12D15 HbT'},...
%     {'S13D28 HbO','S13D28 HbR','S13D28 HbT'},...
%     {'S14D29 HbO','S14D29 HbR','S14D29 HbT'}};

chnNameGrp = {{'S1D1 HbO','S1D1 HbR','S1D1 HbT','S1D2 HbO','S1D2 HbR','S1D2 HbT','S1D3 HbO','S1D3 HbR','S1D3 HbT','S1D4 HbO','S1D4 HbR','S1D4 HbT'},...
    {'S2D1 HbO','S2D1 HbR','S2D1 HbT','S2D2 HbO','S2D2 HbR','S2D2 HbT'},...
    {'S3D3 HbO','S3D3 HbR','S3D3 HbT','S3D4 HbO','S3D4 HbR','S3D4 HbT','S3D5 HbO','S3D5 HbR','S3D5 HbT','S3D6 HbO','S3D6 HbR','S3D6 HbT'},...
    {'S4D9 HbO','S4D9 HbR','S4D9 HbT','S4D10 HbO','S4D10 HbR','S4D10 HbT','S4D11 HbO','S4D11 HbR','S4D11 HbT','S4D12 HbO','S4D12 HbR','S4D12 HbT'},...
    {'S5D7 HbO','S5D7 HbR','S5D7 HbT','S5D8 HbO','S5D8 HbR','S5D8 HbT','S5D9 HbO','S5D9 HbR','S5D9 HbT','S5D10 HbO','S5D10 HbR','S5D10 HbT'},...
    {'S6D7 HbO','S6D7 HbR','S6D7 HbT','S6D8 HbO','S6D8 HbR','S6D8 HbT'},...
    {'S7D18 HbO','S7D18 HbR','S7D18 HbT'},...
    {'S8D19 HbO','S8D19 HbR','S8D19 HbT'},...
    {'S9D13 HbO','S9D13 HbR','S9D13 HbT','S9D14 HbO','S9D14 HbR','S9D14 HbT','S9D16 HbO','S9D16 HbR','S9D16 HbT'},...
    {'S10D13 HbO','S10D13 HbR','S10D13 HbT','S10D15 HbO','S10D15 HbR','S10D15 HbT','S10D17 HbO','S10D17 HbR','S10D17 HbT'},...
    {'S11D13 HbO','S11D13 HbR','S11D13 HbT','S11D14 HbO','S11D14 HbR','S11D14 HbT'},...
    {'S12D13 HbO','S12D13 HbR','S12D13 HbT','S12D15 HbO','S12D15 HbR','S12D15 HbT'},...
    {'S13D28 HbO','S13D28 HbR','S13D28 HbT'},...
    {'S14D29 HbO','S14D29 HbR','S14D29 HbT'}};

% #46-57

colorList = {[0, 0.447, 0.741],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],...
    [0.4940, 0.1840, 0.5560],	[0.4660, 0.6740, 0.1880] 	,[0.3010, 0.7450, 0.9330]};

maxPlot = zeros(1,length(srcIdxGrp));
minPlot = zeros(1,length(srcIdxGrp));

% find max/min of 6 groups for plotting
for i1 = 1:length(srcIdxGrp)
    tempMax = zeros(1,6);
    tempMin = zeros(1,6);
    for i2 = 1:6
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

for i1 = 1:length(srcIdxGrp)
    figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    for i2 = 1:6
        subplot(3,2,condSubPlotIdx(i2)); hold on;
        for j = 1:length(srcIdxGrp{i1}) % source 1
            plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),'-','Color',colorList{j-0});
            plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2),':','Color',colorList{j-0});
            plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3),'-.','Color',colorList{j-0});
        end
        ylim([minPlot(i1) maxPlot(i1)]);
        legend(chnNameGrp{i1},'Location','northeastoutside');
        %srcStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).sourceIndex);
        %detectorStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).detectorIndex);
        ylabel('[M*mm]')
        xlabel('Time [s]');
        xlim([-2 15]);
        condStr = snirf2.stim(1,condSubPlotIdx(i2)).name;
        title(sprintf('%s',condStr));
    end
    hold off;
    fn = sprintf('HRF_S%s',num2str(i1));
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
    close;
end
% 
% end

% for i1 = 1:length(srcIdxGrp)
%     for j = 1:length(srcIdxGrp{i1})
%         figure('units','normalized','outerposition',[0 0 1 1]); hold on;
%         for i2 = 1:6
%             %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),'-','Color',colorList{j-0});hold on;
%             plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),lineStyle{i2},'Color',colorList{colorIdx(i2)});hold on;
%             %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2),':','Color',colorList{j-0});hold on;
%             %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3),'-.','Color',colorList{j-0});hold on;
%         end
%         %legend(chnNameGrp{i1},'Location','northeastoutside');
%         legend({'Left Single','Center Single','Right Multiple','Right Single','Left Multiple','Center Multiple'},'Location','northeastoutside');
%         srcStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).sourceIndex);
%         detectorStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).detectorIndex);
%         ylabel('[M*mm]')
%         xlabel('Time [s]');
%         xlim([-2 15]);
%         %label = num2str(dcavg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).dataTypeLabel);
%         %condStr = snirf1.stim(1,i2).name;
%         title(sprintf('S%sD%s',srcStr,detectorStr));
%         %title(sprintf('%s',condStr));
%         hold off;
%         fn = sprintf('Cond_HRF_S%sD%s',srcStr,detectorStr);
%         print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
%         close;
%     end
% end

end
