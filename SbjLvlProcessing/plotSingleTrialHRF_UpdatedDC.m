% different from redefineStimulusClassesUpdated02.
% removed prune channels step to observe auditory cortical regions.
% Updated to work with updated singleTrialHRFCellHbO
% plotSingleTrialHRF_UpdatedDC_New is no longer in use.
function plotSingleTrialHRF_UpdatedDC(sbjNum)
saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
saveDir2 = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
load([saveDir filesep 'intermediateOutputsCombined.mat']);
load([saveDir filesep 'singleTrialsUpdated.mat'],'singleTrialHRFHbOS','singleTrialHRFHbRS','singleTrialHRFHbTS',...
    'singleTrialHRFHbOM','singleTrialHRFHbRM','singleTrialHRFHbTM','indexMoviesTest');
% Convert nirs to snirf file format
%snirf1 = SnirfClass(load('20210504_1348_02.nirs','-mat'));
%snirf1.Save('20210504_1348_02.snirf');
%snirf1.Info()

%dcAvgReshape = group.subjs(1,1).runs(1).procStream.output.dcAvg.GetDataTimeSeries('reshape');
%tHRF = group.subjs(1,1).runs.procStream.output.dcAvg.GetTime();
fs = 50;

numConds = 6;
numConc = 3;
numChns = size(dcAvg.measurementList,2)/(numConds*numConc);
sizeCond = size(dcAvg.measurementList,2)/numConds;
sizeChn = sizeCond/numChns;
% time x (cond x chn x conc)
% Look at S1D1, S1D2, S1D3, & S1D4
% For validation first. Plot 4 channels for one condition in one figure
% for i = 1:numConds
% for i = 1:1
%     figure(); hold on;
%     % for j = 1:numChns
%     for j = 1:1
%         plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i-1) + sizeChn*(j-1) + 1));
%         plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i-1) + sizeChn*(j-1) + 2));
%         plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i-1) + sizeChn*(j-1) + 3)); hold off;
%     end
%     legend('HbO','HbR','HbT');
%     srcStr = num2str(dcAvg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).sourceIndex);
%     detectorStr = num2str(dcAvg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).detectorIndex);
%     %label = num2str(dcavg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).dataTypeLabel);
%     condStr = snirf1.stim(1,i).name;
%     %titleStr = [condStr 'S' srcStr 'D' detectorStr];
%     title(sprintf('%s S%sD%s',condStr,srcStr,detectorStr));
% end

% For actual analysis, plot 6 conditions from one channel in one figure
colorList = {[0, 0.447, 0.741],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560] 	};
% for i = 1:6
%     figure(); hold on;
%     % for j = 1:numChns
%     for j = 1:4
%         plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i-1) + sizeChn*(j-1) + 1),'-','Color',colorList{j});
%         plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i-1) + sizeChn*(j-1) + 2),':','Color',colorList{j});
%         plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i-1) + sizeChn*(j-1) + 3),'-.','Color',colorList{j});
%     end
%     legend('S1D1 HbO','S1D1 HbR','S1D1 HbT','S1D2 HbO','S1D2 HbR','S1D2 HbT','S1D3 HbO','S1D3 HbR','S1D3 HbT','S1D4 HbO','S1D4 HbR','S1D4 HbT');
%     srcStr = num2str(dcAvg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).sourceIndex);
%     detectorStr = num2str(dcAvg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).detectorIndex);
%     %label = num2str(dcavg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).dataTypeLabel);
%     condStr = snirf1.stim(1,i).name;
%     %title(sprintf('%s S%sD%s',condStr,srcStr,detectorStr));
%     title(sprintf('%s',condStr));
%     hold off;
% end


%save('Experiment03.mat','dod','dc','dodAvg','dcAvg','dodAvgStd','dcAvgStd','dodSum2','dcSum2','tHRF','nTrials','ch','grpAvgPass','misc');
% output is a ProcResultClass
% save('Experiment03.mat','output');% First, store list of channels pruned
% Second, store list of trials rejected. Skipped for now bc none rejected
% calculate HRF from single trial
% Use HRF from dcNew variable

% channels pruned stored in mlActAuto var

rejectedChn = find(mlActAuto{1}==0);
% need to save old stim and compare with new stim var
%rejectedTrials = find(

% Pick 4 channels
numChns = size(singleTrialHRFHbOS,1);
%indexMoviesTest = (indexMoviesTest(:,2)==2) | (indexMoviesTest(:,2)==1);
indexMoviesTestS = indexMoviesTest((indexMoviesTest(:,5)==0),:);
indexMoviesTestM = indexMoviesTest((indexMoviesTest(:,5)==1),:);
singleTrialIndex = zeros(size(indexMoviesTest,1)/6,3);
multipleTrialIndex = zeros(size(indexMoviesTest,1)/6,3);
% left
singleTrialIndex(:,1) = find((indexMoviesTestS(:,2)==2));
% right
singleTrialIndex(:,2) = find((indexMoviesTestS(:,2)==1));
% center
singleTrialIndex(:,3) = find((indexMoviesTestS(:,2)==3));
multipleTrialIndex(:,1) = find((indexMoviesTestM(:,2)==2));
multipleTrialIndex(:,2) = find((indexMoviesTestM(:,2)==1));
multipleTrialIndex(:,3) = find((indexMoviesTestM(:,2)==3));

zeroT = 2*fs;

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

for i1 = 2:size(singleTrialHRFHbOM,1)
    for i2 = 1:size(singleTrialHRFHbOM,3)
        offset = singleTrialHRFHbOM(i1,zeroT,i2);
        singleTrialHRFHbOM(i1,:,i2) = singleTrialHRFHbOM(i1,:,i2) - offset;
        
        offset = singleTrialHRFHbRM(i1,zeroT,i2);
        singleTrialHRFHbRM(i1,:,i2) = singleTrialHRFHbRM(i1,:,i2) - offset;
        
        offset = singleTrialHRFHbTM(i1,zeroT,i2);
        singleTrialHRFHbTM(i1,:,i2) = singleTrialHRFHbTM(i1,:,i2) - offset;
    end
end


%condSubPlotIdx = [1 3 5 2 4 6];

% plot all trials
for j = 1:numChns
%for j = 1:3
    figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    for i = 1:3
        subplot(2,3,i); hold on;
        avgHRF = squeeze(mean(singleTrialHRFHbOS(j,:,singleTrialIndex(:,i)),3));
        %for k = 1:size(singleTrialIndex,1)
        for k = fix(size(singleTrialIndex,1)/2):size(singleTrialIndex,1)
            % channels x time x trial
            plot(dcAvg.time,squeeze(singleTrialHRFHbOS(j,:,singleTrialIndex(k,i))),'--','Color',[0.8500 0.3250 0.0980 0.5]);hold on;
        end
        %plot(dcAvg.time,dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i)-1) + sizeChn*(j-1) + 1),'-','Color',[0 0.4470 0.7410 1]);hold on;
        plot(dcAvg.time,dcAvg.dataTimeSeries(:,sizeCond*(i-1) + sizeChn*(j-1) + 1),'-','Color',[0 0.4470 0.7410 1]);hold on;
        plot(dcAvg.time,avgHRF,'-','Color',[0.9290, 0.6940, 0.1250]);hold on;
        title(sprintf('%s S%sD%s',stim(1,i).name,...
            num2str(dcAvg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).sourceIndex),...
            num2str(dcAvg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).detectorIndex)));
        ylabel('[M * mm]');
        xlabel('Time [s]');
        xlim([-2 15]);
        ylim([-3e-4 3e-4]);
    end
    
    for i = 4:6
        subplot(2,3,i); hold on;
        avgHRF = squeeze(mean(singleTrialHRFHbOM(j,:,multipleTrialIndex(:,i-3)),3));
        %for k = 1:size(multipleTrialIndex,1)
        for k = fix(size(multipleTrialIndex,1)/2):size(multipleTrialIndex,1)
            % channels x time x trial
            plot(dcAvg.time,squeeze(singleTrialHRFHbOM(j,:,multipleTrialIndex(k,i-3))),'--','Color',[0.8500 0.3250 0.0980 0.5]);hold on;
        end
        %plot(dcAvg.time,dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i)-1) + sizeChn*(j-1) + 1),'-','Color',[0 0.4470 0.7410 1]);hold on;
        plot(dcAvg.time,dcAvg.dataTimeSeries(:,sizeCond*(i-1) + sizeChn*(j-1) + 1),'-','Color',[0 0.4470 0.7410 1]);hold on;
        plot(dcAvg.time,avgHRF,'-','Color',[0.9290, 0.6940, 0.1250]);hold on;
        title(sprintf('%s S%sD%s',stim(1,i).name,...
            num2str(dcAvg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).sourceIndex),...
            num2str(dcAvg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).detectorIndex)));
        ylabel('[M * mm]');
        xlabel('Time [s]');
        xlim([-2 15]);
        ylim([-3e-4 3e-4]);
    end
    hold off;
end

end
