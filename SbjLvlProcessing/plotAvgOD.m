% Plot delta OD using outputs from scripting (dod is OD)
% need to average first
% plotHRF.m plot HRF using outputs from group.mat from Homer3 GUI
figSaveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment03\Figures';
saveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment03';
load([saveDir filesep 'intermediateOutputs.mat']);

numConds = 6;
numConc = 3;
numChns = size(dcAvg.measurementList,2)/(numConds*numConc);
sizeCond = size(dcAvg.measurementList,2)/numConds;
sizeChn = sizeCond/numChns;
fs = 50;

condSubPlotIdx = [1 3 5 2 4 6];

srcIdxGrp = {[1 2 3 4] [5 6] [8 9 10 11] [12 13 14 15] [16 17 18 19] [20 21] [23] [25] ...
    [27 28 29] [31 32 33] [35 36] [37 38] [40] [42]};

% chnNameGrp = {{'S1D1','S1D2','S1D3','S1D4'},{'S2D1','S2D2'},{'S3D3','S3D4','S3D5','S3D6'},...
%     {'S4D9','S4D10','S4D11','S4D12'},{'S5D7','S5D8','S5D9','S5D10'},{'S6D7','S6D8'},...
%     {'S7D18'},{'S8D19'},{'S9D13','S9D14','S9D16'},{'S10D13','S10D15','S10D17'},...
%     {'S11D13','S11D14'},{'S12D13','S12D15'},{'S13D28'},{'S14D29'}};

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

% chnNameGrp = {{'','','S1D1 830nm','S1D1 670nm','','','S1D2 830nm','S1D2 670nm','','','S1D3 830nm','S1D3 670nm','','','S1D4 830nm','S1D4 670nm'},...
%     {'','','S2D1 830nm','S2D1 670 nm','','','S2D2 830nm','S2D2 670nm'},...
%     {'','','S3D3 830nm','S3D3 670nm','','','S3D4 830nm','S3D4 670nm','','','S3D5 830nm','S3D5 670nm','','','S3D6 830nm','S3D6 670nm'},...
%     {'','','S4D9 830nm','S4D9 670nm','','','S4D10 830nm','S4D10 670nm','','','S4D11 830nm','S4D11 670nm','','','S4D12 830nm','S4D12 670nm'},...
%     {'','','S5D7 830nm','S5D7 670nm','','','S5D8 830nm','S5D8 670nm','','','S5D9 830nm','S5D9 670nm','','','S5D10 830nm','S5D10 670nm'},...
%     {'','','S6D7 830nm','S6D7 670nm','','','S6D8 830nm','S6D8 670nm'},...
%     {'','','S7D18 830nm','S7D18 670nm'},...
%     {'','','S8D19 830nm','S8D19 670nm'},...
%     {'','','S9D13 830nm','S9D13 670nm','','','S9D14 830nm','S9D14 670nm','','','S9D16 830nm','S9D16 670nm'},...
%     {'','','S10D13 830nm','S10D13 670nm','','','S10D15 830nm','S10D15 670nm','','','S10D17 830nm','S10D17 670nm'},...
%     {'','','S11D13 830nm','S11D13 670nm','','','S11D14 830nm','S11D14 670nm'},...
%     {'','','S12D13 830nm','S12D13 670nm','','','S12D15 830nm','S12D15 670nm'},...
%     {'','','S13D28 830nm','S13D28 670nm'},...
%     {'','','S14D29 830nm','S14D29 670nm'}};

% chnNameGrp = {{'S1D1 830nm','S1D1 670nm','S1D2 830nm','S1D2 670nm','S1D3 830nm','S1D3 670nm','S1D4 830nm','S1D4 670nm'},...
%     {'S2D1 830nm','S2D1 670 nm','S2D2 830nm','S2D2 670nm'},...
%     {'S3D3 830nm','S3D3 670nm','S3D4 830nm','S3D4 670nm','S3D5 830nm','S3D5 670nm','S3D6 830nm','S3D6 670nm'},...
%     {'S4D9 830nm','S4D9 670nm','S4D10 830nm','S4D10 670nm','S4D11 830nm','S4D11 670nm','S4D12 830nm','S4D12 670nm'},...
%     {'S5D7 830nm','S5D7 670nm','S5D8 830nm','S5D8 670nm','S5D9 830nm','S5D9 670nm','S5D10 830nm','S5D10 670nm'},...
%     {'S6D7 830nm','S6D7 670nm','S6D8 830nm','S6D8 670nm'},...
%     {'S7D18 830nm','S7D18 670nm'},...
%     {'S8D19 830nm','S8D19 670nm'},...
%     {'S9D13 830nm','S9D13 670nm','S9D14 830nm','S9D14 670nm','S9D16 830nm','S9D16 670nm'},...
%     {'S10D13 830nm','S10D13 670nm','S10D15 830nm','S10D15 670nm','S10D17 830nm','S10D17 670nm'},...
%     {'S11D13 830nm','S11D13 670nm','S11D14 830nm','S11D14 670nm'},...
%     {'S12D13 830nm','S12D13 670nm','S12D15 830nm','S12D15 670nm'},...
%     {'S13D28 830nm','S13D28 670nm'},...
%     {'S14D29 830nm','S14D29 670nm'}};

chnNameGrp = {'Left Single','Right Single','Center Single','Left Multi','Right Multi','Center Multi'};

% #46-57

colorList1 = {[0, 0.447, 0.741],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560] 	};
colorList2 = {[0, 0.48, 0.741],[0.8500, 0.36, 0.0980],[0.9290, 0.64, 0.1250],[0.4940, 0.22, 0.5560] 	};

avgRaw830 = zeros(numChns,length(dcAvg.time));
avgRaw670 = zeros(numChns,length(dcAvg.time));
sdRaw830 = zeros(numChns,length(dcAvg.time));
sdRaw670 = zeros(numChns,length(dcAvg.time));
avgTrialHRFCell830 = cell(6,1);
avgTrialHRFCell670 = cell(6,1);
sdTrialHRFCell830 = cell(6,1);
sdTrialHRFCell670 = cell(6,1);

for i1 = 1:numConds
    % chns x trials x time
    chnXtrials830 = zeros(numChns,size(stim(1,i1).data,1),length(dcAvg.time));
    chnXtrials670 = zeros(numChns,size(stim(1,i1).data,1),length(dcAvg.time));
    for i2 = 1:size(stim(1,i1).data,1)
        startPt = stim(1,i1).data(i2,1)*fs-2*fs;
        endPt = stim(1,i1).data(i2,1)*fs+15*fs;
        for i3 = 1:numChns
            if endPt < size(dod.dataTimeSeries,1)
                chnXtrials830(i3,i2,:) = squeeze(dod.dataTimeSeries(startPt:endPt,i3));
                chnXtrials670(i3,i2,:) = squeeze(dod.dataTimeSeries(startPt:endPt,i3+42));
            end
        end
    end
    avgRaw830 = squeeze(mean(chnXtrials830,2));
    avgRaw670 = squeeze(mean(chnXtrials670,2));
    sdRaw830 = squeeze(std(chnXtrials830,0,2));
    sdRaw670 = squeeze(std(chnXtrials670,0,2));
    avgTrialHRFCell830{i1} = avgRaw830;
    avgTrialHRFCell670{i1} = avgRaw670;
    sdTrialHRFCell830{i1} = sdRaw830;
    sdTrialHRFCell670{i1} = sdRaw670;
end

% each figure has 6 different conditions for one channel.
for i1 = 1:length(srcIdxGrp)
    for j = 1:length(srcIdxGrp{i1}) % source 1
        figure('units','normalized','outerposition',[0 0 1 1]); hold on;
        for i2 = 1:numConds
        
            avgRaw830 = avgTrialHRFCell830{i2};
            avgRaw670 = avgTrialHRFCell670{i2};
            sdRaw830 = sdTrialHRFCell830{i2};
            sdRaw670 = sdTrialHRFCell670{i2};
        
            %subplot(3,2,condSubPlotIdx(i2)); hold on;
        
        %for j = 16:19 % 46-57, source 5
        %for j = 20:21 % source 6
        %for j = 5:6
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),'-','Color',colorList{j-0});
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2),':','Color',colorList{j-0});
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3),'-.','Color',colorList{j-0});
            
            curve1 = avgRaw830(srcIdxGrp{i1}(j),:) + sdRaw830(srcIdxGrp{i1}(j),:);
            curve2 = avgRaw830(srcIdxGrp{i1}(j),:) - sdRaw830(srcIdxGrp{i1}(j),:);
            x2 = [dcAvg.time', fliplr(dcAvg.time')];
            inBetween = [curve1, fliplr(curve2)];
            %h = fill(x2, inBetween, colorList1{j-0});
            %set(h,'facealpha',.5);
            %set(h,'edgealpha',.5);
            hold on;
            
            curve1 = avgRaw670(srcIdxGrp{i1}(j),:) + sdRaw670(srcIdxGrp{i1}(j),:);
            curve2 = avgRaw670(srcIdxGrp{i1}(j),:) - sdRaw670(srcIdxGrp{i1}(j),:);
            %x2 = [x, fliplr(x)];
            inBetween = [curve1, fliplr(curve2)];
            %h = fill(x2, inBetween, colorList2{j-0});
            %set(h,'facealpha',.5)
            %set(h,'edgealpha',.5);
            hold on;
            %plot(dcAvg.time, avgRaw830(srcIdxGrp{i1}(j),:),'-','Color',colorList1{j-0});hold on;
            plot(dcAvg.time, avgRaw830(srcIdxGrp{i1}(j),:));hold on;
            %plot(dcAvg.time, avgRaw670(srcIdxGrp{i1}(j),:),':','Color',colorList2{j-0});hold on;
        end
        
        %legend(chnNameGrp{i1},'Location','northeastoutside');
        legend(chnNameGrp,'Location','northeastoutside');
        %legend('S1D1 HbO','S1D1 HbR','S1D1 HbT','S1D2 HbO','S1D2 HbR','S1D2 HbT','S1D3 HbO','S1D3 HbR','S1D3 HbT','S1D4 HbO','S1D4 HbR','S1D4 HbT');
        %legend('S5D7 HbO','S5D7 HbR','S5D7 HbT','S5D8 HbO','S5D8 HbR','S5D8 HbT','S5D9 HbO','S5D9 HbR','S5D9 HbT','S5D10 HbO','S5D10 HbR','S5D10 HbT');
        %legend('S6D7 HbO','S6D7 HbR','S6D7 HbT','S6D8 HbO','S6D8 HbR','S6D8 HbT');
        %legend('S2D1 HbO','S2D1 HbR','S2D1 HbT','S2D2 HbO','S2D2 HbR','S2D2 HbT');
        srcStr = num2str(dcAvg.measurementList(1,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).sourceIndex);
        detectorStr = num2str(dcAvg.measurementList(1,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).detectorIndex);
        if (mlActAuto{1}(srcIdxGrp{i1}(j)) == 1 && mlActAuto{1}(srcIdxGrp{i1}(j)+42)== 1)
            pruned = '';
        else
            pruned = 'Pruned';
        end
        ylabel('AU');
        xlabel('Time [s]');
        xlim([-2 15]);
        %label = num2str(dcavg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).dataTypeLabel);
        condStr = snirf1.stim(1,i2).name;
        %title(sprintf('%s S%sD%s',condStr,srcStr,detectorStr));
        title(sprintf('OD S%sD%s %s',srcStr, detectorStr, pruned));
        fn = sprintf('OD_S%sD%s',num2str(i1),detectorStr);
        %print(gcf,[figSaveDir filesep fn],'-dpng','-r250');

    end
    hold off;
end


% figure has 6 subplots, each different conditions. Each subplot has all
% channels associated with one source
% for i1 = 1:length(srcIdxGrp)
%     figure('units','normalized','outerposition',[0 0 1 1]); hold on;
%     for i2 = 1:numConds
%         
%         avgRaw830 = avgTrialHRFCell830{i2};
%         avgRaw670 = avgTrialHRFCell670{i2};
%         sdRaw830 = sdTrialHRFCell830{i2};
%         sdRaw670 = sdTrialHRFCell670{i2};
%         
%         subplot(3,2,condSubPlotIdx(i2)); hold on;
%         for j = 1:length(srcIdxGrp{i1}) % source 1
%         %for j = 16:19 % 46-57, source 5
%         %for j = 20:21 % source 6
%         %for j = 5:6
%             %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),'-','Color',colorList{j-0});
%             %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2),':','Color',colorList{j-0});
%             %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3),'-.','Color',colorList{j-0});
%             
%             curve1 = avgRaw830(srcIdxGrp{i1}(j),:) + sdRaw830(srcIdxGrp{i1}(j),:);
%             curve2 = avgRaw830(srcIdxGrp{i1}(j),:) - sdRaw830(srcIdxGrp{i1}(j),:);
%             x2 = [dcAvg.time', fliplr(dcAvg.time')];
%             inBetween = [curve1, fliplr(curve2)];
%             %h = fill(x2, inBetween, colorList1{j-0});
%             %set(h,'facealpha',.5);
%             %set(h,'edgealpha',.5);
%             hold on;
%             
%             curve1 = avgRaw670(srcIdxGrp{i1}(j),:) + sdRaw670(srcIdxGrp{i1}(j),:);
%             curve2 = avgRaw670(srcIdxGrp{i1}(j),:) - sdRaw670(srcIdxGrp{i1}(j),:);
%             %x2 = [x, fliplr(x)];
%             inBetween = [curve1, fliplr(curve2)];
%             %h = fill(x2, inBetween, colorList2{j-0});
%             %set(h,'facealpha',.5)
%             %set(h,'edgealpha',.5);
%             hold on;
%             plot(dcAvg.time, avgRaw830(srcIdxGrp{i1}(j),:),'-','Color',colorList1{j-0});hold on;
%             plot(dcAvg.time, avgRaw670(srcIdxGrp{i1}(j),:),':','Color',colorList2{j-0});
%         end
%         
%         legend(chnNameGrp{i1},'Location','northeastoutside');
%         %legend('S1D1 HbO','S1D1 HbR','S1D1 HbT','S1D2 HbO','S1D2 HbR','S1D2 HbT','S1D3 HbO','S1D3 HbR','S1D3 HbT','S1D4 HbO','S1D4 HbR','S1D4 HbT');
%         %legend('S5D7 HbO','S5D7 HbR','S5D7 HbT','S5D8 HbO','S5D8 HbR','S5D8 HbT','S5D9 HbO','S5D9 HbR','S5D9 HbT','S5D10 HbO','S5D10 HbR','S5D10 HbT');
%         %legend('S6D7 HbO','S6D7 HbR','S6D7 HbT','S6D8 HbO','S6D8 HbR','S6D8 HbT');
%         %legend('S2D1 HbO','S2D1 HbR','S2D1 HbT','S2D2 HbO','S2D2 HbR','S2D2 HbT');
%         srcStr = num2str(dcAvg.measurementList(1,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).sourceIndex);
%         detectorStr = num2str(dcAvg.measurementList(1,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).detectorIndex);
%         ylabel('AU')
%         xlabel('Time [s]');
%         xlim([-2 15]);
%         %label = num2str(dcavg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).dataTypeLabel);
%         condStr = snirf1.stim(1,i2).name;
%         %title(sprintf('%s S%sD%s',condStr,srcStr,detectorStr));
%         title(sprintf('Raw %s',condStr));
%         fn = sprintf('Raw_S%s',num2str(i1));
%         %print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
% 
%     end
%     hold off;
% end