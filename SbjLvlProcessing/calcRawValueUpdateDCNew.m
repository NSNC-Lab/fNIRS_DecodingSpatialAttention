% split into single and multi groups.
% new. Use trials from dcNew var
% Use different method of concatenating to validate the code in
% calcCumSumUpdateDCNew.m

function calcRawValueUpdateDCNew(sbjNum)
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
load([processedDataDir filesep 'singleTrialsUpdated_LR.mat'],'singleTrialHRFHbOS','singleTrialHRFHbOM',...
    'singleTrialHRFHbRS','singleTrialHRFHbRM',...
    'singleTrialHRFHbTS','singleTrialHRFHbTM','indexMoviesTest');

numClasses = 2;

timePt = (0:0.5*fs:12*fs)+2*fs;
%timePt = 6*fs+2*fs;

%tInd = find(timePt==(6.5*fs+2*fs));

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

performanceArrSingleHbO = zeros(size(singleTrialHRFHbOS,1),length(timePt));
performanceArrMultiHbO = zeros(size(singleTrialHRFHbOM,1),length(timePt));
performanceArrSingleHbR = zeros(size(singleTrialHRFHbRS,1),length(timePt));
performanceArrMultiHbR = zeros(size(singleTrialHRFHbRM,1),length(timePt));
performanceArrSingleHbT = zeros(size(singleTrialHRFHbTS,1),length(timePt));
performanceArrMultiHbT = zeros(size(singleTrialHRFHbTM,1),length(timePt));

predRespSHbO = zeros(size(singleTrialHRFHbOS,1),size(singleTrialHRFHbOS,3),length(timePt));
predRespMHbO = zeros(size(singleTrialHRFHbOM,1),size(singleTrialHRFHbOM,3),length(timePt));
predRespSHbR = zeros(size(singleTrialHRFHbRS,1),size(singleTrialHRFHbRS,3),length(timePt));
predRespMHbR = zeros(size(singleTrialHRFHbRM,1),size(singleTrialHRFHbRM,3),length(timePt));
predRespSHbT = zeros(size(singleTrialHRFHbTS,1),size(singleTrialHRFHbTS,3),length(timePt));
predRespMHbT = zeros(size(singleTrialHRFHbTM,1),size(singleTrialHRFHbTM,3),length(timePt));

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

singleIndex = (indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==0) | (indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==0);
multipleIndex = (indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1) | (indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1);
indexMoviesTestSingle = indexMoviesTest(singleIndex,:);
indexMoviesTestMultiple = indexMoviesTest(multipleIndex,:);

for i = 1:length(timePt)
    rawValueHRFHbOS(:,1,:) = singleTrialHRFHbOS(:,timePt(i),:);
    rawValueHRFHbOM(:,1,:) = singleTrialHRFHbOM(:,timePt(i),:);
    rawValueHRFHbRS(:,1,:) = singleTrialHRFHbRS(:,timePt(i),:);
    rawValueHRFHbRM(:,1,:) = singleTrialHRFHbRM(:,timePt(i),:);
    rawValueHRFHbTS(:,1,:) = singleTrialHRFHbTS(:,timePt(i),:);
    rawValueHRFHbTM(:,1,:) = singleTrialHRFHbTM(:,timePt(i),:);

    [performanceArrSingleHbO(:,i),predRespSHbO(:,:,i)] = trainClassifierLinearDiscriminantfNIRS(rawValueHRFHbOS,mlActAuto{1},numClasses,indexMoviesTestSingle);
    [performanceArrMultiHbO(:,i),predRespMHbO(:,:,i)] = trainClassifierLinearDiscriminantfNIRS(rawValueHRFHbOM,mlActAuto{1},numClasses,indexMoviesTestMultiple);

    [performanceArrSingleHbR(:,i),predRespSHbR(:,:,i)] = trainClassifierLinearDiscriminantfNIRS(rawValueHRFHbRS,mlActAuto{1},numClasses,indexMoviesTestSingle);
    [performanceArrMultiHbR(:,i),predRespMHbR(:,:,i)] = trainClassifierLinearDiscriminantfNIRS(rawValueHRFHbRM,mlActAuto{1},numClasses,indexMoviesTestMultiple);

    [performanceArrSingleHbT(:,i),predRespSHbT(:,:,i)] = trainClassifierLinearDiscriminantfNIRS(rawValueHRFHbTS,mlActAuto{1},numClasses,indexMoviesTestSingle);
    [performanceArrMultiHbT(:,i),predRespMHbT(:,:,i)] = trainClassifierLinearDiscriminantfNIRS(rawValueHRFHbTM,mlActAuto{1},numClasses,indexMoviesTestMultiple);

end

savePerfFN = 'performanceRawLinearDiscriminantUpdated_LR.mat';
save([processedDataDir filesep savePerfFN],'performanceArrMultiHbO','performanceArrSingleHbO',...
    'performanceArrMultiHbR','performanceArrSingleHbR',...
    'performanceArrMultiHbT','performanceArrSingleHbT',...
    'predRespSHbO','predRespMHbO','singleTrialHRFHbOS','singleTrialHRFHbOM',...
    'singleTrialHRFHbRS','singleTrialHRFHbRM',...
    'singleTrialHRFHbTS','singleTrialHRFHbTM',...
    'predRespSHbR','predRespMHbR',...
    'predRespSHbT','predRespMHbT');

plotPerformanceVsTime(sbjNum,savePerfFN,processedDataDir,'Raw')

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
    
    
% % plot performance as a function of ending cumulative sum time point    
% figure();hold on;
% for j = 1:size(performanceArrSingleHbO,1)
%     %plot(timePt./fs-2,performanceArrSingle(sizeChn*(j-1) + 1,:));
%     plot(timePt./fs-2,performanceArrSingleHbO(j,:));
% end
% title(sprintf('Sbj %s: Δ[HbO] Single Movie',num2str(sbjNum)));
% %legend('S1D1','S1D2','S1D3','S1D4');
% xlabel('Time [s]');ylabel('Accuracy');
% hold off;
% 
% figure();hold on;
% for j = 1:size(performanceArrMultiHbO,1)
%     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
%     plot(timePt./fs-2,performanceArrMultiHbO(j,:));
% end
% title(sprintf('Sbj %s: Δ[HbO] Multi Movie',num2str(sbjNum)));
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