% Status: Inactive
% 
% split into single and multi groups.
% new. Use trials from dcNew var
% Plot bar chart by channel using ROI names.
% Use this one.

function calcCumSumUpdateDCNew_BarChart_ChnName(sbjNum,numClasses)
% channels pruned stored in mlActAuto var
processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
%saveDir2 = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment03';
load([processedDataDir filesep 'intermediateOutputsCombined.mat']);

%rejectedChn = find(mlActAuto{1}==0);
% need to save old stim and compare with new stim var
%rejectedTrials = find(

% Channels index grouped by src # Used for plotting
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
probe = snirf2.probe;
Aaux = [];
rcMap = [];

% Pick 4 channels
numConds = 6;
numConc = 3;
numChns = size(dcAvg.measurementList,2)/(numConds*numConc);
sizeCond = size(dcAvg.measurementList,2)/numConds;
sizeChn = sizeCond/numChns;

% Cell is array of 6 different conditions
% Each array is channels x time x trial
if numClasses == 2
    load([processedDataDir filesep 'singleTrialsUpdated_LR.mat'],'singleTrialHRFHbOS','singleTrialHRFHbOM',...
        'singleTrialHRFHbRS','singleTrialHRFHbRM',...
        'singleTrialHRFHbTS','singleTrialHRFHbTM','indexMoviesTest');
else
    load([processedDataDir filesep 'singleTrialsUpdated.mat'],'singleTrialHRFHbOS','singleTrialHRFHbOM',...
        'singleTrialHRFHbRS','singleTrialHRFHbRM',...
        'singleTrialHRFHbTS','singleTrialHRFHbTM','indexMoviesTest');
end
%load([saveDir2 filesep 'iniList_alh_202154_.mat']);

% HbO
%singleTrialsHbO = cat(1,ones(1,size(singleTrialHRFHbOS,2),size(singleTrialHRFHbOS,3)),singleTrialHRFHbOS);
%leftTrialsSHbO = cat(1,2*ones(1,size(singleTrialHRFHbOS,2),size(singleTrialHRFHbOS,3)),singleTrialHRFHbOS);
%centerTrialsSHbO = cat(1,3*ones(1,size(singleTrialHRFCellHbO{3},2),size(singleTrialHRFCellHbO{3},3)),singleTrialHRFCellHbO{3});
%singleTrialsHbO = cat(3,rightTrialsSHbO,leftTrialsSHbO);
%singleTrialsHbO = cat(3,singleTrialsHbO,centerTrialsSHbO);

%multiTrialsHbO = cat(1,ones(1,size(singleTrialHRFHbOS{5},2),size(singleTrialHRFHbOS{5},3)),singleTrialHRFHbOS{5});
%leftTrialsMHbO = cat(1,2*ones(1,size(singleTrialHRFHbOS{4},2),size(singleTrialHRFHbOS{4},3)),singleTrialHRFHbOS{4});
%centerTrialsMHbO = cat(1,3*ones(1,size(singleTrialHRFCellHbO{6},2),size(singleTrialHRFCellHbO{6},3)),singleTrialHRFCellHbO{6});
%multiTrialsHbO = cat(3,rightTrialsMHbO,leftTrialsMHbO);
%multiTrialsHbO = cat(3,multiTrialsHbO,centerTrialsMHbO);

%HbR
% rightTrialsSHbR = cat(1,ones(1,size(singleTrialHRFCellHbR{2},2),size(singleTrialHRFCellHbR{2},3)),singleTrialHRFCellHbR{2});
% leftTrialsSHbR = cat(1,2*ones(1,size(singleTrialHRFCellHbR{1},2),size(singleTrialHRFCellHbR{1},3)),singleTrialHRFCellHbR{1});
% %centerTrialsSHbR = cat(1,3*ones(1,size(singleTrialHRFCellHbR{3},2),size(singleTrialHRFCellHbR{3},3)),singleTrialHRFCellHbR{3});
% singleTrialsHbR = cat(3,rightTrialsSHbR,leftTrialsSHbR);
% %singleTrialsHbR = cat(3,singleTrialsHbR,centerTrialsSHbR);
% 
% rightTrialsMHbR = cat(1,ones(1,size(singleTrialHRFCellHbR{5},2),size(singleTrialHRFCellHbR{5},3)),singleTrialHRFCellHbR{5});
% leftTrialsMHbR = cat(1,2*ones(1,size(singleTrialHRFCellHbR{4},2),size(singleTrialHRFCellHbR{4},3)),singleTrialHRFCellHbR{4});
% %centerTrialsMHbR = cat(1,3*ones(1,size(singleTrialHRFCellHbR{6},2),size(singleTrialHRFCellHbR{6},3)),singleTrialHRFCellHbR{6});
% multiTrialsHbR = cat(3,rightTrialsMHbR,leftTrialsMHbR);
% %multiTrialsHbR = cat(3,multiTrialsHbR,centerTrialsMHbR);
% 
% %HbT
% rightTrialsSHbT = cat(1,ones(1,size(singleTrialHRFCellHbT{2},2),size(singleTrialHRFCellHbT{2},3)),singleTrialHRFCellHbT{2});
% leftTrialsSHbT = cat(1,2*ones(1,size(singleTrialHRFCellHbT{1},2),size(singleTrialHRFCellHbT{1},3)),singleTrialHRFCellHbT{1});
% %centerTrialsSHbT = cat(1,3*ones(1,size(singleTrialHRFCellHbT{3},2),size(singleTrialHRFCellHbT{3},3)),singleTrialHRFCellHbT{3});
% singleTrialsHbT = cat(3,rightTrialsSHbT,leftTrialsSHbT);
% %singleTrialsHbT = cat(3,singleTrialsHbT,centerTrialsSHbT);
% 
% rightTrialsMHbT = cat(1,ones(1,size(singleTrialHRFCellHbT{5},2),size(singleTrialHRFCellHbT{5},3)),singleTrialHRFCellHbT{5});
% leftTrialsMHbT = cat(1,2*ones(1,size(singleTrialHRFCellHbT{4},2),size(singleTrialHRFCellHbT{4},3)),singleTrialHRFCellHbT{4});
% %centerTrialsMHbT = cat(1,3*ones(1,size(singleTrialHRFCellHbT{6},2),size(singleTrialHRFCellHbT{6},3)),singleTrialHRFCellHbT{6});
% multiTrialsHbT = cat(3,rightTrialsMHbT,leftTrialsMHbT);
% %multiTrialsHbT = cat(3,multiTrialsHbT,centerTrialsMHbT);

timePt = (0:0.5*fs:12*fs)+2*fs;
%timePt = 4*fs+2*fs;

tInd = find(timePt==(6*fs+2*fs));

zeroT = 2*fs;

% offset
for i1 = 2:size(singleTrialHRFHbOS,1)
    for i2 = 1:size(singleTrialHRFHbOS,3)
        offset = singleTrialHRFHbOS(i1,zeroT,i2);
        singleTrialHRFHbOS(i1,:,i2) = singleTrialHRFHbOS(i1,:,i2) - offset;
        
        offset = singleTrialHRFHbRS(i1,zeroT,i2);
        singleTrialHRFHbRS(i1,:,i2) = singleTrialHRFHbRS(i1,:,i2) - offset;
        
        offset = singleTrialHRFHbTS(i1,zeroT,i2);
        singleTrialHRFHbTS(i1,:,i2) = singleTrialHRFHbTS(i1,:,i2) - offset;
    end
end

% % baseline correction
% for i1 = 2:size(singleTrials,1)
%     for i2 = 1:size(singleTrials,3)
%         offset = mean(singleTrials(i1,1:zeroT,i2));
%         singleTrials(i1,:,i2) = singleTrials(i1,:,i2) - offset;
%     end
% end

% %validate offset
% figure();
% plot(dcAvg.time,squeeze(singleTrials(3,:,:)));
% title('Single-Movie Condition. All Trials.');
% xlabel('Time [s]');
% ylabel('M*mm');

% Offset
for i1 = 2:size(singleTrialHRFHbOM,1)
    for i2 = 1:size(singleTrialHRFHbOM,3)
        offset = singleTrialHRFHbOM(i1,zeroT,i2);
        %offset = mean(singleTrials(i1,1:zeroT,i2));
        singleTrialHRFHbOM(i1,:,i2) = singleTrialHRFHbOM(i1,:,i2) - offset;
        
        offset = singleTrialHRFHbRM(i1,zeroT,i2);
        singleTrialHRFHbRM(i1,:,i2) = singleTrialHRFHbRM(i1,:,i2) - offset;
        
        offset = singleTrialHRFHbTM(i1,zeroT,i2);
        singleTrialHRFHbTM(i1,:,i2) = singleTrialHRFHbTM(i1,:,i2) - offset;
    end
end

if numClasses == 2
    singleIndex = (indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==0) | (indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==0);
    multipleIndex = (indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1) | (indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1);
else
    singleIndex = (indexMoviesTest(:,5)==0);
    multipleIndex = (indexMoviesTest(:,5)==1);
end
indexMoviesTestSingle = indexMoviesTest(singleIndex,:);
indexMoviesTestMultiple = indexMoviesTest(multipleIndex,:);

startT = 2*fs;

[performanceArrSingleHbO,predRespSHbO] = trainClassifierLinearDiscriminantfNIRS_Cumsum(singleTrialHRFHbOS,startT,timePt,mlActAuto{1},numClasses,indexMoviesTestSingle);
[performanceArrMultiHbO,predRespMHbO] = trainClassifierLinearDiscriminantfNIRS_Cumsum(singleTrialHRFHbOM,startT,timePt,mlActAuto{1},numClasses,indexMoviesTestMultiple);

[performanceArrSingleHbR,predRespSHbR] = trainClassifierLinearDiscriminantfNIRS_Cumsum(singleTrialHRFHbRS,startT,timePt,mlActAuto{1},numClasses,indexMoviesTestSingle);
[performanceArrMultiHbR,predRespMHbR] = trainClassifierLinearDiscriminantfNIRS_Cumsum(singleTrialHRFHbRM,startT,timePt,mlActAuto{1},numClasses,indexMoviesTestMultiple);

[performanceArrSingleHbT,predRespSHbT] = trainClassifierLinearDiscriminantfNIRS_Cumsum(singleTrialHRFHbTS,startT,timePt,mlActAuto{1},numClasses,indexMoviesTestSingle);
[performanceArrMultiHbT,predRespMHbT] = trainClassifierLinearDiscriminantfNIRS_Cumsum(singleTrialHRFHbTM,startT,timePt,mlActAuto{1},numClasses,indexMoviesTestMultiple);

if numClasses == 2
    savePerfFN = 'performanceLinearDiscriminantUpdated_LR.mat';
else
    savePerfFN = 'performanceLinearDiscriminantUpdated.mat';
end
save([processedDataDir filesep savePerfFN],'performanceArrMultiHbO','performanceArrSingleHbO',...
    'performanceArrMultiHbR','performanceArrSingleHbR',...
    'performanceArrMultiHbT','performanceArrSingleHbT',...
    'predRespSHbO','predRespMHbO','singleTrialHRFHbOS','singleTrialHRFHbOM',...
    'singleTrialHRFHbRS','singleTrialHRFHbRM',...
    'singleTrialHRFHbTS','singleTrialHRFHbTM',...
    'predRespSHbR','predRespMHbR',...
    'predRespSHbT','predRespMHbT');

%plotPerformanceCumsum(sbjNum,numClasses)
plotPerformanceBarChartROINames(sbjNum,performanceArrSingleHbO,tInd,'Cumulative Sum ??[HbO] Single',1)
plotPerformanceBarChartROINames(sbjNum,performanceArrSingleHbR,tInd,'Cumulative Sum ??[HbR] Single',1)
plotPerformanceBarChartROINames(sbjNum,performanceArrSingleHbT,tInd,'Cumulative Sum ??[HbT] Single',1)
plotPerformanceBarChartROINames(sbjNum,performanceArrMultiHbO,tInd,'Cumulative Sum ??[HbO] Multiple',1)
plotPerformanceBarChartROINames(sbjNum,performanceArrMultiHbR,tInd,'Cumulative Sum ??[HbR] Multiple',1)
plotPerformanceBarChartROINames(sbjNum,performanceArrMultiHbT,tInd,'Cumulative Sum ??[HbT] Multiple',1)

% % plot performance of channels using one source for each fig
% for i = 1:length(srcIdxGrp)
%     figure('units','normalized','outerposition',[0 0 1 1]); hold on;
%     for j = 1:length(srcIdxGrp{i})
%         plot(timePt./fs-2,performanceArrSingle((srcIdxGrp{i}(j)-1)*3 + 1,:));
%     end
%     titleStr = sprintf('S%s ??[HbO] Single Movie',num2str(i));
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
%     titleStr = sprintf('S%s ??[HbO] Multi Movie',num2str(i));
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
%     titleStr = sprintf('S%s ??[HbO] Single Movie',num2str(srcPairIdxGrp{1}(1)));
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
%     titleStr = sprintf('S%s ??[HbO] Single Movie',num2str(srcPairIdxGrp{1}(2)));
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
    
    
% % plot performance as a function of ending cumulative sum time point    
% figure();hold on;
% for j = 1:size(performanceArrSingleHbO,1)
%     %plot(timePt./fs-2,performanceArrSingle(sizeChn*(j-1) + 1,:));
%     plot(timePt./fs-2,performanceArrSingleHbO(j,:));
% end
% title(sprintf('Sbj %s: ??[HbO] Single Movie',num2str(sbjNum)));
% %legend('S1D1','S1D2','S1D3','S1D4');
% xlabel('Time [s]');ylabel('Accuracy');
% hold off;
% 
% figure();hold on;
% for j = 1:size(performanceArrMultiHbO,1)
%     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
%     plot(timePt./fs-2,performanceArrMultiHbO(j,:));
% end
% title(sprintf('Sbj %s: ??[HbO] Multi Movie',num2str(sbjNum)));
% %legend('S1D1','S1D2','S1D3','S1D4');
% xlabel('Time [s]');ylabel('Accuracy');
% hold off;

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