% For sbj 08 and 10
function plotAllHRF_Special(sbjNum,opCI,opTTest)

figSaveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum) '\Figures\HRFs'];
saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
rawDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
load([saveDir filesep 'intermediateOutputsCombined_Basis1.mat'],'dcAvg','dcAvgStd','dcNew','beta','bvar','dof_full',...
    'sse_full','sse_reduced','dof_full','dof_reduced','hmrstats');
load([rawDir filesep sbjNum '.mat'],'s');

tSize = size(dcNew.dataTimeSeries,1);
% sigAct is HbX x #Channels x #conditions
% sigDiff is HbX x #Channels x #pairs
% beta is also a function. be careful
%[sigAct,sigDiff] = calcZTest(beta{1},bvar,tSize,0);
if opTTest
    sigAct = hmrstats.pval_contrast(1:3,:,:)<0.05;
    sigDiff = hmrstats.pval_contrast(4:6,:,:)<0.05;
else
    [sigAct,sigDiff] = calcFTest(sse_full,sse_reduced,dof_full,dof_reduced);
end

ssIdx = [7,22,24,26,30,34,39,41];
sigAct(:,ssIdx,:) = [];
sigDiff(:,ssIdx,:) = [];

numConds = 3;
numConc = 3;
numChns = size(dcAvg.measurementList,2)/(numConds*numConc);
sizeCond = size(dcAvg.measurementList,2)/numConds;
sizeChn = sizeCond/numChns;
% Original ordeer is: Center Multiple, Center Single, Left M, Left S, RM,
% RS
% L R C
%condSubPlotIdx = [4 6 2];
condSubPlotIdx = [1 2 3];
sigSubPlotIdx = [1 2 3];
%sigSubPlotIdx = [4 6 2];

sigFN = [saveDir filesep 'sig.mat'];
save(sigFN,'sigAct','sigDiff');

colorIdx = [1 2 3 1 2 3];

lineStyle = {'-' '-' '-' '-' '-' '-'};

srcIdxGrp = {[1 2 3 4] [5 6] [8 9 10 11] [12 13 14 15] [16 17 18 19] [20 21] [23] [25] ...
    [27 28 29] [31 32 33] [35 36] [37 38] [40] [42]};

% srcIdxGrp = {[1 2 3 4] [5 6] [8 9 10 11] [12 13 14 15] [16 17 18 19] [20 21] [23] [25] ...
%     [27 28] [30 31] [33 34] [35 36]};

fs = 50;
timePt = (0:0.5*fs:12*fs)+2*fs;

timeIdx = find(timePt==(5*fs+2*fs));

% number of channels, can pull
if strcmp(sbjNum,'04') ||strcmp(sbjNum,'07') ||strcmp(sbjNum,'08') ||strcmp(sbjNum,'10')
    mlActAuto_Unused = getPrunedChns_Special(s,1.5);
else
    mlActAuto_Unused = getPrunedChns(s,1.5);
end

ml690 = mlActAuto_Unused{1}(1:length(mlActAuto_Unused{1})/2);
ml870 = mlActAuto_Unused{1}(length(mlActAuto_Unused{1})/2+1:length(mlActAuto_Unused{1}));
acceptedChns = (ml690==ml870);
%numUnusedChns = length(mlActAuto_Unused{1})/2-sum(ml690==ml870);

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
%     'S9D16 IPS3/latIPS/antIPS Left',...
%     'S10D13 IPS3/IPS2/SPL1 Right',...
%     'S10D15	IPS3/antIPS/IPS4 Right',...
%     'S10D17 IPS3/latIPS/antIPS Right',...
%     'S11D13	Superior to IPS3/IPS2/SPL1 Left',...
%     'S11D14	IPS4 Left',...
%     'S12D13	Superior to IPS3/IPS2/SPL1 Right',...
%     'S12D15	IPS4 Right',...
%     'S13D28 Ant Posterior STG/PT Left',...
%     'S14D29 Ant Posterior STG/PT Right'};

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
    'tgPCS',...20
    'SMG/STG/PT',...
    'SMG/STG/PT',...
    'IPS3/IPS2/SPL1',...
    'IPS3/antIPS/IPS4',...
    'IPS3/latIPS',...
    'IPS3/IPS2/SPL1',...
    'IPS3/antIPS/IPS4',...
    'IPS3/latIPS',...
    'IPS3/IPS2/mSPL',...
    'IPS4/mSPL',...
    'IPS3/IPS2/mSPL',...
    'IPS4/mSPL',...
    'STG/PT',...
    'STG/PT'};

% chnSubPlotInd = [21 20 15 14 27 26 9 8 3 2 10 11 4 5 22 23 16 17 28 29 31 ...
%     36 45 44 46 47 39 38 40 41];

chnSubPlotInd = [21 20 15 14 27 26 9 8 3 2 10 11 4 5 22 23 16 17 28 29 31 ...
    36 45 44 50 46 47 53 39 38 40 41 7 12];

colorList = {[0, 0.447, 0.741],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],...
    [0.4940, 0.1840, 0.5560],	[0.4660, 0.6740, 0.1880] 	,[0.3010, 0.7450, 0.9330]};

maxPlot = zeros(1,length(srcIdxGrp));
minPlot = zeros(1,length(srcIdxGrp));

timeSecPostSh = [dcAvg.time; flipud(dcAvg.time)];

%dof = dof_full{1};
dof = dof_full;

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

%R = R{1};

chnNum = 0;

figure('units','normalized','outerposition',[0 0 1 1]); hold on;
for i1 = 1:length(srcIdxGrp)
    for j = 1:length(srcIdxGrp{i1})
        chnNum = chnNum + 1;
        
        subplot(9,6,chnSubPlotInd(chnNum));hold on;
        
        xline(0);
        xline(5);
        yline(0);
        
        if opCI
            %nTrials = 30;
            for i2 = 1:3
                
                y_se = dcAvgStd.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1);
                % still biased estimator. should be nTrials-nB-1;
                % dof > 10000 so already close to normal dist
                y_ci = tinv(1-(1-0.95)/2,dof)*y_se;
                inBet = [squeeze(dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1)...
                    + y_ci); ...
                    flipud(squeeze(dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1)...
                    - y_ci))];
                h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');

                set(h(i2),'facealpha',0.1);
            end
        else
            for i2 = 1:3
                inBet = [squeeze(dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1)...
                    +dcAvgStd.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1)); ...
                    flipud(squeeze(dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1)...
                    -dcAvgStd.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1)))];
                h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');

                set(h(i2),'facealpha',0.1);
            end
        end
        
        
        for i2 = 1:3
            
            if sigAct(i2,chnNum,1)
                lineWidth = 2;
            else
                lineWidth = 0.5;
            end
            
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),'-','Color',colorList{j-0});hold on;
            plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),lineStyle{i2},...
                'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2),':','Color',colorList{j-0});hold on;
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3),'-.','Color',colorList{j-0});hold on;
        end
        %legend(chnNameGrp{i1},'Location','northeastoutside');
        %legend({'Left Single','Center Single','Right Multiple','Right Single','Left Multiple','Center Multiple'},'Location','northeast');
%         srcStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).sourceIndex);
%         detectorStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).detectorIndex);
%         ylabel('[M*mm]')
%         xlabel('Time [s]');

        %if sigDiff(1,chnNum,condSubPlotIdx(1))
        if sigDiff(1,chnNum,1)
            LRMark = '*';
        else
            LRMark = ' ';
        end
        
        %if sigDiff(1,chnNum,condSubPlotIdx(2))
        if sigDiff(2,chnNum,1)
            LCMark = '+';
        else
            LCMark = ' ';
        end
        
        %if sigDiff(1,chnNum,condSubPlotIdx(3))
        if sigDiff(3,chnNum,1)
            RCMark = 'x';
        else
            RCMark = ' ';
        end
                
        titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
        title(chnName(chnNum));
        subtitle(titleStr);
        xlim([-2 15]);
        set(gca,'XTick',[], 'YTick', [])

        xlim([-2 15]);
%         if ~acceptedChns(chnNum)
%             ax = gca;
%             ax.XColor = 'red';
%             ax.YColor = 'red';
%         end
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
%legend({'Left Single','Center Single','Right Multiple','Right Single','Left Multiple','Center Multiple'},'Location','northeast');
legend('Left Multiple','Center Multiple','Right Multiple');

sgtitle(sprintf('Sbj %s HbO',sbjNum));
hold off;
fn = sprintf('HbO_AllHRF');
print(gcf,[figSaveDir filesep fn],'-dpng','-r250');

chnNum = 0;

%%%%%%%%%%%%
figure('units','normalized','outerposition',[0 0 1 1]); hold on;
for i1 = 1:length(srcIdxGrp)
    for j = 1:length(srcIdxGrp{i1})
        chnNum = chnNum + 1;

        subplot(9,6,chnSubPlotInd(chnNum));hold on;
                
        xline(0);
        xline(5);
        yline(0);
        
        if opCI
            for i2 = 1:3
                y_se = dcAvgStd.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2);
                % still biased estimator. should be nTrials-nB-1;
                y_ci = tinv(1-(1-0.95)/2,dof)*y_se;
                inBet = [squeeze(dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2)...
                    + y_ci); ...
                    flipud(squeeze(dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2)...
                    - y_ci))];
                h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');

                set(h(i2),'facealpha',0.1);
            end
        else
            for i2 = 1:3
                inBet = [squeeze(dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2)...
                    +dcAvgStd.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2)); ...
                    flipud(squeeze(dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2)...
                    -dcAvgStd.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2)))];
                h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');

                set(h(i2),'facealpha',0.1);
            end
        end
        
        for i2 = 1:3
            
            if sigAct(i2,chnNum,2)
                lineWidth = 2;
            else
                lineWidth = 0.5;
            end
            
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),'-','Color',colorList{j-0});hold on;
            plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2),lineStyle{i2},...
                'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;
        
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2),':','Color',colorList{j-0});hold on;
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3),'-.','Color',colorList{j-0});hold on;
        end
        %legend(chnNameGrp{i1},'Location','northeastoutside');
%         legend({'Left Single','Center Single','Right Multiple','Right Single','Left Multiple','Center Multiple'},'Location','northeast');
%         srcStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).sourceIndex);
%         detectorStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).detectorIndex);
%         ylabel('[M*mm]')
%         xlabel('Time [s]');
        xlim([-2 15]);
        %label = num2str(dcavg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).dataTypeLabel);
        %condStr = snirf1.stim(1,i2).name;
%         title(sprintf('S%sD%s R: %s CAS: %s CAM: %s',...
%             srcStr,detectorStr,num2str(R(srcIdxGrp{i1}(j),2)),...
%             num2str(performanceArrSingleHbR(srcIdxGrp{i1}(j),timeIdx)),...
%             num2str(performanceArrMultiHbR(srcIdxGrp{i1}(j),timeIdx))));
%         subtitle(sprintf('HbR %s',chnName{chnNum}));
        
        %title(sprintf('%s',condStr));
        
        %if sigDiff(2,chnNum,condSubPlotIdx(1))
        if sigDiff(1,chnNum,2)
            LRMark = '*';
        else
            LRMark = ' ';
        end
        
        %if sigDiff(2,chnNum,condSubPlotIdx(2))
        if sigDiff(2,chnNum,2)
            LCMark = '+';
        else
            LCMark = ' ';
        end
        
        %if sigDiff(2,chnNum,condSubPlotIdx(3))
        if sigDiff(3,chnNum,2)
            RCMark = 'x';
        else
            RCMark = ' ';
        end
                
        titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
        title(chnName(chnNum));
        subtitle(titleStr);
        xlim([-2 15]);
        set(gca,'XTick',[], 'YTick', [])
        
%         if ~acceptedChns(chnNum)
%             ax = gca;
%             ax.XColor = 'red';
%             ax.YColor = 'red';
%         end
        hold off;
%         fn = sprintf('HRF_S%sD%s',srcStr,detectorStr);
%         print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
%         close;
    end
end
sgtitle(sprintf('Sbj %s HbR',sbjNum));
subplot(9,6,54);hold on;
plot(0,0);
plot(0,0);
plot(0,0);
%legend({'Left Single','Center Single','Right Multiple','Right Single','Left Multiple','Center Multiple'},'Location','northeast');
legend('Left Multiple','Center Multiple','Right Multiple');

%sgtitle(sprintf('Sbj %s HbO',sbjNum));
hold off;
fn = sprintf('HbR_AllHRF');
print(gcf,[figSaveDir filesep fn],'-dpng','-r250');

% HbT
chnNum = 0;
figure('units','normalized','outerposition',[0 0 1 1]); hold on;
%tiledlayout(8,6,'TileSpacing','tight','Padding','tight');hold on;
for i1 = 1:length(srcIdxGrp)
    for j = 1:length(srcIdxGrp{i1})
        chnNum = chnNum + 1;

        %subplot(8,6,chnSubPlotInd(chnNum));hold on;
        subplot(9,6,chnSubPlotInd(chnNum));hold on;
        %nexttile(chnSubPlotInd(chnNum));hold on;
        
        xline(0);
        xline(5);
        yline(0);
        if opCI
            %dof = 30;
            for i2 = 1:3
                y_se = dcAvgStd.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3);
                % still biased estimator. should be nTrials-nB-1;
                y_ci = tinv(1-(1-0.95)/2,dof)*y_se;
                inBet = [squeeze(dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3)...
                    + y_ci); ...
                    flipud(squeeze(dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3)...
                    - y_ci))];
                h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');

                set(h(i2),'facealpha',0.1);
            end
        else
            for i2 = 1:3
                inBet = [squeeze(dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3)...
                    +dcAvgStd.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3)); ...
                    flipud(squeeze(dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3)...
                    -dcAvgStd.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3)))];
                h(i2) = fill(timeSecPostSh, inBet, colorList{colorIdx(i2)},'LineStyle','none');

                set(h(i2),'facealpha',0.1);
            end
        end
        
        for i2 = 1:3
            if sigAct(3,chnNum,i2)
                lineWidth = 2;
            else
                lineWidth = 0.5;
            end
            
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),'-','Color',colorList{j-0});hold on;
            plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3),lineStyle{i2},...
                'Color',colorList{colorIdx(i2)},'LineWidth',lineWidth);hold on;
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2),':','Color',colorList{j-0});hold on;
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3),'-.','Color',colorList{j-0});hold on;
        end
        %legend(chnNameGrp{i1},'Location','northeastoutside');
%         legend({'Left Single','Center Single','Right Multiple','Right Single','Left Multiple','Center Multiple'},'Location','northeast');
%         srcStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).sourceIndex);
%         detectorStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).detectorIndex);
%         ylabel('[M*mm]')
%         xlabel('Time [s]');
        xlim([-2 15]);
        set(gca,'XTick',[], 'YTick', [])
        %label = num2str(dcavg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).dataTypeLabel);
        %condStr = snirf1.stim(1,i2).name;
%         title(sprintf('S%sD%s R: %s CAS: %s CAM: %s',...
%             srcStr,detectorStr,num2str(R(srcIdxGrp{i1}(j),2)),...
%             num2str(performanceArrSingleHbR(srcIdxGrp{i1}(j),timeIdx)),...
%             num2str(performanceArrMultiHbR(srcIdxGrp{i1}(j),timeIdx))));
%         subtitle(sprintf('HbR %s',chnName{chnNum}));
        
        %title(sprintf('%s',condStr));
        
        if sigDiff(3,chnNum,1)
            LRMark = '*';
        else
            LRMark = ' ';
        end
        
        if sigDiff(3,chnNum,2)
            LCMark = '+';
        else
            LCMark = ' ';
        end
        
        if sigDiff(3,chnNum,3)
            RCMark = 'x';
        else
            RCMark = ' ';
        end
                
        titleStr = sprintf('%s        %s        %s',LRMark,LCMark,RCMark);
        title(chnName(chnNum),'FontSize',8);
        subtitle(titleStr);
        
%         if ~acceptedChns(chnNum)
%             ax = gca;
%             ax.XColor = 'red';
%             ax.YColor = 'red';
%         end
        hold off;
%         fn = sprintf('HRF_S%sD%s',srcStr,detectorStr);
%         print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
%         close;
    end
end
%sgtitle(sprintf('Sbj %s HbT',sbjNum));
subplot(9,6,54);hold on;
%nexttile(48);hold on;
plot(0,0);
plot(0,0);
plot(0,0);
legend('Left','Center','Right');
%legend('Right','Left','Center');
title(sprintf('Sbj %s HbT',sbjNum));
%legend('Left Multiple','Right Multiple','Center Multiple');

%sgtitle(sprintf('Sbj %s HbO',sbjNum));
hold off;
fn = sprintf('HbT_AllHRF');
print(gcf,[figSaveDir filesep fn],'-dpng','-r250');


end