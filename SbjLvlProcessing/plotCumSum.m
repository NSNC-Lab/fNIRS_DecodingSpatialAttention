saveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment03';
figSaveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment03\Figures\Classification';
load([saveDir filesep 'intermediateOutputs.mat']);
load([saveDir filesep 'singleTrials.mat'],'singleTrialHRFCell');
load([saveDir filesep 'iniList_alh_202154_.mat']);
load([saveDir filesep 'performance.mat']);
%savePerfFN = 'performance.mat';

srcIdxGrp = {[1 2 3 4] [5 6] [8 9 10 11] [12 13 14 15] [16 17 18 19] [20 21] [23] [25] ...
    [27 28 29] [31 32 33] [35 36] [37 38] [40] [42]};

srcPairIdxGrp = {[1 5] [2 6] [3 4] [7 8] [9 10] [11 12] [13 14]};
    
chnName = {'S1D1','S1D2','S1D3','S1D4','S2D1','S2D2','S3D3','S3D4','S3D5','S3D6',...
    'S4D9','S4D10','S4D11','S4D12','S5D7','S5D8','S5D9','S5D10','S6D7','S6D8',...
    'S7D18','S8D19','S9D13','S9D14','S9D16','S10D13','S10D15','S10D17',...
    'S11D13','S11D14','S12D13','S12D15','S13D28','S14D29'};

chnNameGrp = {{'S1D1','S1D2','S1D3','S1D4'},{'S2D1','S2D2'},{'S3D3','S3D4','S3D5','S3D6'},...
    {'S4D9','S4D10','S4D11','S4D12'},{'S5D7','S5D8','S5D9','S5D10'},{'S6D7','S6D8'},...
    {'S7D18'},{'S8D19'},{'S9D13','S9D14','S9D16'},{'S10D13','S10D15','S10D17'},...
    {'S11D13','S11D14'},{'S12D13','S12D15'},{'S13D28'},{'S14D29'}};

chnInd = ([1 2 3 4 5 6 8 9 10 11 12 13 14 15 16 17 18 19 20 21 23 25 27 28 29 ...
    31 32 33 35 36 37 38 40 42]- 1) *3 + 1;

fs = 50;
timePt = (0:0.5*fs:12*fs)+2*fs;

numConds = 6;
numConc = 3;
numChns = size(dcAvg.measurementList,2)/(numConds*numConc);
sizeCond = size(dcAvg.measurementList,2)/numConds;
sizeChn = sizeCond/numChns;

colorList = {[0, 0.447, 0.741],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560] 	};

for i = 1:length(srcPairIdxGrp)
    figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    
    subplot(4,2,1);hold on;
    for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(1)})
        plot(timePt./fs-2,performanceArrSingle((srcIdxGrp{i}(j)-1)*3 + 1,:));
    end
    titleStr = sprintf('S%s Δ[HbO] Single Movie',num2str(srcPairIdxGrp{i}(1)));
    title(titleStr);
    legend(chnNameGrp{srcPairIdxGrp{i}(1)});
    xlabel('Time [s]');ylabel('Accuracy');
    ylim([0 1]);
    
    subplot(4,2,2);hold on;
    for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(2)})
        plot(timePt./fs-2,performanceArrSingle((srcIdxGrp{srcPairIdxGrp{i}(2)}(j)-1)*3 + 1,:));
    end
    titleStr = sprintf('S%s Δ[HbO] Single Movie',num2str(srcPairIdxGrp{i}(2)));
    title(titleStr);
    legend(chnNameGrp{srcPairIdxGrp{i}(2)});
    xlabel('Time [s]');ylabel('Accuracy');
    ylim([0 1]);
    
    subplot(4,2,3);hold on;
    for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(1)})
        plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(0)+(srcIdxGrp{srcPairIdxGrp{i}(1)}(j)-1)*3+1),'-','Color',colorList{j-0});
    end
    legend(chnNameGrp{srcPairIdxGrp{i}(1)});
    ylabel('[M*mm]')
    xlabel('Time [s]');
    xlim([-2 15]);
    title('Left Condition. Left Hemi.');
    ylim([-2.5e-4 2.5e-4]);
    
    subplot(4,2,4);hold on;
    for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(2)})
        plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(0)+(srcIdxGrp{srcPairIdxGrp{i}(2)}(j)-1)*3+1),'-','Color',colorList{j-0});
    end
    legend(chnNameGrp{srcPairIdxGrp{i}(2)});
    ylabel('[M*mm]')
    xlabel('Time [s]');
    xlim([-2 15]);
    title('Left Condition. Right Hemi.');
    ylim([-2.5e-4 2.5e-4]);
    
    subplot(4,2,5);hold on;
    for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(1)})
        plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(1)+(srcIdxGrp{srcPairIdxGrp{i}(1)}(j)-1)*3+1),'-','Color',colorList{j-0});
    end
    legend(chnNameGrp{srcPairIdxGrp{i}(1)});
    ylabel('[M*mm]')
    xlabel('Time [s]');
    xlim([-2 15]);
    title('Right Condition. Left Hemi.');
    ylim([-2.5e-4 2.5e-4]);
    
    subplot(4,2,6);hold on;
    for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(2)})
        plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(1)+(srcIdxGrp{srcPairIdxGrp{i}(2)}(j)-1)*3+1),'-','Color',colorList{j-0});
    end
    legend(chnNameGrp{srcPairIdxGrp{i}(2)});
    ylabel('[M*mm]')
    xlabel('Time [s]');
    xlim([-2 15]);
    title('Right Condition. Right Hemi.');
    ylim([-2.5e-4 2.5e-4]);
    
    subplot(4,2,7);hold on;
    for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(1)})
        plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(2)+(srcIdxGrp{srcPairIdxGrp{i}(1)}(j)-1)*3+1),'-','Color',colorList{j-0});
    end
    legend(chnNameGrp{srcPairIdxGrp{i}(1)});
    ylabel('[M*mm]')
    xlabel('Time [s]');
    xlim([-2 15]);
    title('Center Condition. Left Hemi.');
    ylim([-2.5e-4 2.5e-4]);
    
    subplot(4,2,8);hold on;
    for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(2)})
        plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(2)+(srcIdxGrp{srcPairIdxGrp{i}(2)}(j)-1)*3+1),'-','Color',colorList{j-0});
    end
    legend(chnNameGrp{srcPairIdxGrp{i}(2)});
    ylabel('[M*mm]')
    xlabel('Time [s]');
    xlim([-2 15]);
    title('Center Condition. Right Hemi.');
    ylim([-2.5e-4 2.5e-4]);

    fnStr = sprintf('Single_S%sHbO',num2str(i));
    print(gcf,[figSaveDir filesep fnStr],'-dpng','-r250');

end







for i = 1:length(srcPairIdxGrp)
    figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    
    subplot(4,2,1);hold on;
    for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(1)})
        plot(timePt./fs-2,performanceArrMulti((srcIdxGrp{i}(j)-1)*3 + 1,:));
    end
    titleStr = sprintf('S%s Δ[HbO] Multi Movie',num2str(srcPairIdxGrp{i}(1)));
    fnStr = sprintf('Single_S%sHbO',num2str(srcPairIdxGrp{1}(1)));
    title(titleStr);
    legend(chnNameGrp{srcPairIdxGrp{i}(1)});
    xlabel('Time [s]');ylabel('Accuracy');
    ylim([0 1]);
    
    subplot(4,2,2);hold on;
    for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(2)})
        plot(timePt./fs-2,performanceArrMulti((srcIdxGrp{srcPairIdxGrp{i}(2)}(j)-1)*3 + 1,:));
    end
    titleStr = sprintf('S%s Δ[HbO] Multi Movie',num2str(srcPairIdxGrp{i}(2)));
    title(titleStr);
    legend(chnNameGrp{srcPairIdxGrp{i}(2)});
    xlabel('Time [s]');ylabel('Accuracy');
    ylim([0 1]);
    
    subplot(4,2,3);hold on;
    for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(1)})
        plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(3)+(srcIdxGrp{srcPairIdxGrp{i}(1)}(j)-1)*3+1),'-','Color',colorList{j-0});
    end
    legend(chnNameGrp{srcPairIdxGrp{i}(1)});
    ylabel('[M*mm]')
    xlabel('Time [s]');
    xlim([-2 15]);
    title('Left Condition. Left Hemi.');
    ylim([-2.5e-4 2.5e-4]);
    
    subplot(4,2,4);hold on;
    for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(2)})
        plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(3)+(srcIdxGrp{srcPairIdxGrp{i}(2)}(j)-1)*3+1),'-','Color',colorList{j-0});
    end
    legend(chnNameGrp{srcPairIdxGrp{i}(2)});
    ylabel('[M*mm]')
    xlabel('Time [s]');
    xlim([-2 15]);
    title('Left Condition. Right Hemi.');
    ylim([-2.5e-4 2.5e-4]);
    
    subplot(4,2,5);hold on;
    for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(1)})
        plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(4)+(srcIdxGrp{srcPairIdxGrp{i}(1)}(j)-1)*3+1),'-','Color',colorList{j-0});
    end
    legend(chnNameGrp{srcPairIdxGrp{i}(1)});
    ylabel('[M*mm]')
    xlabel('Time [s]');
    xlim([-2 15]);
    title('Right Condition. Left Hemi.');
    ylim([-2.5e-4 2.5e-4]);
    
    subplot(4,2,6);hold on;
    for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(2)})
        plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(4)+(srcIdxGrp{srcPairIdxGrp{i}(2)}(j)-1)*3+1),'-','Color',colorList{j-0});
    end
    legend(chnNameGrp{srcPairIdxGrp{i}(2)});
    ylabel('[M*mm]')
    xlabel('Time [s]');
    xlim([-2 15]);
    title('Right Condition. Right Hemi.');
    ylim([-2.5e-4 2.5e-4]);
    
    subplot(4,2,7);hold on;
    for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(1)})
        plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(5)+(srcIdxGrp{srcPairIdxGrp{i}(1)}(j)-1)*3+1),'-','Color',colorList{j-0});
    end
    legend(chnNameGrp{srcPairIdxGrp{i}(1)});
    ylabel('[M*mm]')
    xlabel('Time [s]');
    xlim([-2 15]);
    title('Center Condition. Left Hemi.');
    ylim([-2.5e-4 2.5e-4]);
    
    subplot(4,2,8);hold on;
    for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(2)})
        plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(5)+(srcIdxGrp{srcPairIdxGrp{i}(2)}(j)-1)*3+1),'-','Color',colorList{j-0});
    end
    legend(chnNameGrp{srcPairIdxGrp{i}(2)});
    ylabel('[M*mm]')
    xlabel('Time [s]');
    xlim([-2 15]);
    title('Center Condition. Right Hemi.');
    ylim([-2.5e-4 2.5e-4]);

    fnStr = sprintf('Multi_S%sHbO',num2str(i));
    print(gcf,[figSaveDir filesep fnStr],'-dpng','-r250');

end