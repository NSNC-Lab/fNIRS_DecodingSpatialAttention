% split into single and multi groups.
% new. Use trials from dcNew var

function plotPerformanceVsTimeAllChns(sbjNum,fileName,fileDir,featureName)
% channels pruned stored in mlActAuto var
%saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
%saveDir2 = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
% this file created by 
%load([saveDir filesep 'intermediateOutputsCombined.mat']);
load([fileDir filesep fileName],'performanceArrMultiHbO','performanceArrSingleHbO',...
    'performanceArrMultiHbR','performanceArrSingleHbR',...
    'performanceArrMultiHbT','performanceArrSingleHbT',...
    'predRespSHbO','predRespMHbO','singleTrialHRFHbOS','singleTrialHRFHbOM',...
    'singleTrialHRFHbRS','singleTrialHRFHbRM',...
    'singleTrialHRFHbTS','singleTrialHRFHbTM',...
    'predRespSHbR','predRespMHbR',...
    'predRespSHbT','predRespMHbT');

fs = 50;

% % Pick 4 channels
% numConds = 6;
% numConc = 3;
% numChns = size(dcAvg.measurementList,2)/(numConds*numConc);
% sizeCond = size(dcAvg.measurementList,2)/numConds;
% sizeChn = sizeCond/numChns;

% Cell is array of 6 different conditions
% Each array is channels x time x trial
% created by plotSingleTrialHRF_Updated.m
%load([saveDir2 filesep 'iniList_edq_2021614_.mat']);

%numClasses = 3;

timePt = (0:0.5*fs:12*fs)+2*fs;
%timePt = 6*fs+2*fs;

zeroT = 2*fs;

% % plot performance of channels using one source for each fig
% for i = 1:length(srcIdxGrp)
%     figure('units','normalized','outerposition',[0 0 1 1]); hold on;
%     for j = 1:length(srcIdxGrp{i})
%         plot(timePt./fs-2,performanceArrSingle((srcIdxGrp{i}(j)-1)*3 + 1,:));
%     end
%     titleStr = sprintf('S%s Δ[HbO] Single Movie',num2str(i));
%     fnStr = sprintf('Single_S%sHbO',num2str(i));
%     title(titleStr);
%     legend(chnNameGrp{i});
%     xlabel('Time [s]');ylabel('Accuracy');
%     hold off;
%     print(gcf,fnStr,'-dpng','-r250');
% end
% 
% for i = 1:length(srcIdxGrp)
%     figure('units','normalized','outerposition',[0 0 1 1]); hold on;
%     for j = 1:length(srcIdxGrp{i})
%         plot(timePt./fs-2,performanceArrMulti((srcIdxGrp{i}(j)-1)*3 + 1,:));
%     end
%     titleStr = sprintf('S%s Δ[HbO] Multi Movie',num2str(i));
%     fnStr = sprintf('Multi_S%sHbO',num2str(i));
%     title(titleStr);
%     legend(chnNameGrp{i});
%     xlabel('Time [s]');ylabel('Accuracy');
%     hold off;
%     print(gcf,fnStr,'-dpng','-r250');
% end

% % each figure divided into 8 subplots. channels corresponding to left and
% right hemisphere. top row is classification performance. bottom row is
% HRF of all trials for each condition

% colorList = {[0, 0.447, 0.741],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560] 	};
% 
% for i = 1:length(srcPairIdxGrp)
%     figure('units','normalized','outerposition',[0 0 1 1]); hold on;
%     
%     subplot(4,2,1);hold on;
%     for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(1)})
%         plot(timePt./fs-2,performanceArrSingle((srcIdxGrp{i}(j)-1)*3 + 1,:));
%     end
%     titleStr = sprintf('S%s Δ[HbO] Single Movie',num2str(srcPairIdxGrp{1}(1)));
%     fnStr = sprintf('Single_S%sHbO',num2str(srcPairIdxGrp{1}(1)));
%     title(titleStr);
%     legend(chnNameGrp{i});
%     xlabel('Time [s]');ylabel('Accuracy');
%     ylim([0 1]);
%     
%     subplot(4,2,2);hold on;
%     for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(2)})
%         plot(timePt./fs-2,performanceArrSingle((srcIdxGrp{srcPairIdxGrp{i}(2)}(j)-1)*3 + 1,:));
%     end
%     titleStr = sprintf('S%s Δ[HbO] Single Movie',num2str(srcPairIdxGrp{1}(2)));
%     title(titleStr);
%     legend(chnNameGrp{i});
%     xlabel('Time [s]');ylabel('Accuracy');
%     ylim([0 1]);
%     
%     subplot(4,2,3);hold on;
%     for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(1)})
%         plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(0)+chnInd(srcIdxGrp{srcPairIdxGrp{i}(1)}(j))),'-','Color',colorList{j-0});
%     end
%     legend(chnNameGrp{srcPairIdxGrp{i}(1)});
%     ylabel('[M*mm]')
%     xlabel('Time [s]');
%     xlim([-2 15]);
%     title('Right Condition. Left Hemi.');
%     
%     subplot(4,2,4);hold on;
%     for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(2)})
%         plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(0)+chnInd(srcIdxGrp{srcPairIdxGrp{i}(2)}(j))),'-','Color',colorList{j-0});
%     end
%     legend(chnNameGrp{srcPairIdxGrp{i}(2)});
%     ylabel('[M*mm]')
%     xlabel('Time [s]');
%     xlim([-2 15]);
%     title('Right Condition. Right Hemi.');
%     
%     subplot(4,2,5);hold on;
%     for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(1)})
%         plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(1)+chnInd(srcIdxGrp{srcPairIdxGrp{i}(1)}(j))),'-','Color',colorList{j-0});
%     end
%     legend(chnNameGrp{srcPairIdxGrp{i}(1)});
%     ylabel('[M*mm]')
%     xlabel('Time [s]');
%     xlim([-2 15]);
%     title('Left Condition. Left Hemi.');
%     
%     subplot(4,2,6);hold on;
%     for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(2)})
%         plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(1)+chnInd(srcIdxGrp{srcPairIdxGrp{i}(2)}(j))),'-','Color',colorList{j-0});
%     end
%     legend(chnNameGrp{srcPairIdxGrp{i}(2)});
%     ylabel('[M*mm]')
%     xlabel('Time [s]');
%     xlim([-2 15]);
%     title('Left Condition. Right Hemi.');
%     
%     subplot(4,2,7);hold on;
%     for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(1)})
%         plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(2)+chnInd(srcIdxGrp{srcPairIdxGrp{i}(1)}(j))),'-','Color',colorList{j-0});
%     end
%     legend(chnNameGrp{srcPairIdxGrp{i}(1)});
%     ylabel('[M*mm]')
%     xlabel('Time [s]');
%     xlim([-2 15]);
%     title('Center Condition. Left Hemi.');
%     
%     subplot(4,2,8);hold on;
%     for j = 1:length(srcIdxGrp{srcPairIdxGrp{i}(2)})
%         plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(2)+chnInd(srcIdxGrp{srcPairIdxGrp{i}(2)}(j))),'-','Color',colorList{j-0});
%     end
%     legend(chnNameGrp{srcPairIdxGrp{i}(2)});
%     ylabel('[M*mm]')
%     xlabel('Time [s]');
%     xlim([-2 15]);
%     title('Center Condition. Right Hemi.');
% 
%     fnStr = sprintf('Single_S%sHbO',num2str(i));
%     %print(gcf,fnStr,'-dpng','-r250');
% 
% end
    
    
% plot performance as a function of ending cumulative sum time point    
figure();hold on;
plot(timePt./fs-2,performanceArrSingleHbO);
plot(timePt./fs-2,performanceArrMultiHbO);
plot(timePt./fs-2,performanceArrSingleHbR);
plot(timePt./fs-2,performanceArrMultiHbR);
plot(timePt./fs-2,performanceArrSingleHbT);
plot(timePt./fs-2,performanceArrMultiHbT);
ylim([0 1]);
title(sprintf('Sbj %s: Δ[HbO] %s Single Movie',num2str(sbjNum),featureName));
legend('HbO Single','HbO Multi','HbR Single','HbR Multi','HbT Single','HbT Multi');
xlabel('Time [s]');ylabel('Accuracy');
hold off;



% % plot performance as bar chart for each channel
% figure();
% X = categorical(chnName);
% X = reordercats(X,chnName);
% 
% bar(X,performanceArrSingle(chnInd,tInd));
% title('Performance for Single Condition at t=6.5s');
% 
% figure();
% X = categorical(chnName);
% X = reordercats(X,chnName);
% 
% bar(X,performanceArrMulti(chnInd,tInd));
% title('Performance for Multi Condition at t=6.5s');
end
