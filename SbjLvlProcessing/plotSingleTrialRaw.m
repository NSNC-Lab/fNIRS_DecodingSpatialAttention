% different from redefineStimulusClassesUpdated02.
% removed prune channels step to observe auditory cortical regions.
% For sbj 12 and afterward
function plotSingleTrialRaw(sbjNum)

saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
rawDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
figSaveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum) '\Figures\HRFs'];
load([rawDir filesep sbjNum '.mat'],'s');

numConds = 3;
numConc = 3;

% Convert nirs to snirf file format
snirf1 = SnirfClass(load([s.fName{1} '.nirs'],'-mat'));
%snirf3 = SnirfClass(load([rawDataFN{2} '.nirs'],'-mat'));
%snirf4 = SnirfClass(load([rawDataFN{3} '.nirs'],'-mat'));
snirf1.Info()

% Extract aux and convert to stimclass
allS1 = find(snirf1.aux(1,1).dataTimeSeries>1);
aInd1 = find(allS1(2:end)-allS1(1:end-1)==1);
% allS3 = find(snirf3.aux(1,1).dataTimeSeries>1);
% aInd3 = find(allS3(2:end)-allS3(1:end-1)==1);
% allS4 = find(snirf4.aux(1,1).dataTimeSeries>1);
% aInd4 = find(allS4(2:end)-allS4(1:end-1)==1);
allS1(aInd1) = [];
allS1 = allS1./50;
% allS3(aInd3) = [];
% allS3 = allS3./50;
% allS4(aInd4) = [];
% allS4 = allS4./50;

oldStimClass1Length1 = size(allS1,1);
% oldStimClass1Length3 = oldStimClass1Length1 + size(allS3,1);
% oldStimClass1Length4 = oldStimClass1Length3 + size(allS4,1);

cueOnsetIndex = 1:4:720;
%cueOnsetIndex = 1:4:360;
cueOnsetIndex1 = 1:oldStimClass1Length1;
% cueOnsetIndex3 = oldStimClass1Length1+1:oldStimClass1Length3;
% cueOnsetIndex4 = oldStimClass1Length3+1:oldStimClass1Length4;

trialNum1 = sum(ismember(cueOnsetIndex1,cueOnsetIndex));
% trialNum3 = sum(ismember(cueOnsetIndex3,cueOnsetIndex));
% trialNum4 = sum(ismember(cueOnsetIndex4,cueOnsetIndex));

lastTrial1 = trialNum1;
% lastTrial3 = lastTrial1 + trialNum3;
% lastTrial4 = lastTrial3 + trialNum4;

snirf1TLen = size(snirf1.data.time,1)/(50);
% snirf3TLen = size(snirf3.data.time,1)/(50);
% allS3 = allS3+snirf1TLen;
% allS4 = allS4+snirf1TLen+snirf3TLen;

%allS = [allS1; allS3; allS4];
allS = allS1;

% split triggers into 6 categories
load([saveDir filesep s.movieList '.mat'],'fixedMaskerList','indexMoviesTest','maskerMovies','numTrials','uniqueMovies');

% snirf 2
% leftOnsetS = StimClass('leftSingle');
% rightOnsetS = StimClass('rightSingle');
% centerOnsetS = StimClass('centerSingle');
leftOnsetM = StimClass('leftMulti');
rightOnsetM = StimClass('rightMulti');
centerOnsetM = StimClass('centerMulti');

% leftSingleIndex = indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==0;
% rightSingleIndex = indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==0;
% centerSingleIndex = indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==0;
leftMultiIndex = indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1;
rightMultiIndex = indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1;
centerMultiIndex = indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==1;

% leftSingleIndex = leftSingleIndex(leftSingleIndex<= lastTrial2);
% rightSingleIndex = rightSingleIndex(rightSingleIndex<= lastTrial2);
% centerSingleIndex = centerSingleIndex(centerSingleIndex<= lastTrial2);
% leftMultiIndex = leftMultiIndex(leftMultiIndex<= lastTrial2);
% rightMultiIndex = rightMultiIndex(rightMultiIndex<= lastTrial2);
% centerMultiIndex = centerMultiIndex(centerMultiIndex<= lastTrial2);

% AddStims(leftOnsetS, allS(cueOnsetIndex(leftSingleIndex)));
% AddStims(rightOnsetS, allS(cueOnsetIndex(rightSingleIndex)));
% AddStims(centerOnsetS, allS(cueOnsetIndex(centerSingleIndex)));
AddStims(leftOnsetM, allS(cueOnsetIndex(leftMultiIndex)));
AddStims(rightOnsetM, allS(cueOnsetIndex(rightMultiIndex)));
AddStims(centerOnsetM, allS(cueOnsetIndex(centerMultiIndex)));

% updateStates(leftOnsetS);
% updateStates(rightOnsetS);
% updateStates(centerOnsetS);
updateStates(leftOnsetM);
updateStates(rightOnsetM);
updateStates(centerOnsetM);

% snirf2.stim(1,1) = leftOnsetS;
% snirf2.stim(1,2) = rightOnsetS;
% snirf2.stim(1,3) = centerOnsetS;
snirf1.stim(1,1) = leftOnsetM;
snirf1.stim(1,2) = rightOnsetM;
snirf1.stim(1,3) = centerOnsetM;

%numChns = size(dcAvg.measurementList,2)/(numConds*numConc);
numChns = size(snirf1.data.dataTimeSeries,2)/2;
%sizeCond = size(dcAvg.measurementList,2)/numConds;
%sizeChn = sizeCond/numChns;

snirf1.data.dataTimeSeries;

condSubPlotIdx = [2 1 3];
fs = 50;

timeAxis = -1:1/fs:7;

zeroT = 2*fs;

% plot all trials
for j = 1:numChns
    figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    for i = 1:length(snirf1.stim)
        
        %thisId = indexMoviesTest(:,2)==condSubPlotIdx(i);
        %thisTrials = squeeze(snirf1.data.dataTimeSeires(j,:,thisId));
        %for k = 1:size(thisTrials,2)
        %for k = 1:size(snirf1.stim(i).data,1)
        for k = 1:8
            subplot(8,3,(k-1)*3+i); hold on;
            % singleTrialHRFHbOM = zeros(size(dcAvg.dataTimeSeries,2)/chromXCond,size(dcAvg.time,1),size(allSMultiple,1)-subDim);
            %plot(dcAvg.time,squeeze(thisTrials(:,k)),'--','Color',[0.8500 0.3250 0.0980 0.5]);hold on;
            temp = snirf1.data.dataTimeSeries((snirf1.stim(i).data(k,1)-1)*fs:(snirf1.stim(i).data(k,1)+7)*fs,j);
            %temp = lowpass(temp - temp(zeroT),10,fs);
            temp = (temp - temp(zeroT));
            plot(timeAxis,temp,'Color',[0.8500 0.3250 0.0980 0.5]);
            if k==1
                title(sprintf('%s S%sD%s',snirf1.stim(1,i).name,...
                    num2str(snirf1.data.measurementList(1,j).sourceIndex),...
                    num2str(snirf1.data.measurementList(1,j).detectorIndex)));
            end
        end
        %plot(dcAvg.time,dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i)-1) + sizeChn*(j-1) + 1),'-','Color',[0 0.4470 0.7410 1]);hold on;
        %plot(dcAvg.time,dcAvg.dataTimeSeries(:,sizeCond*(i-1) + sizeChn*(j-1) + 1),'-','Color',[0 0.4470 0.7410 1]);hold on;
%         title(sprintf('%s S%sD%s',stim(1,i).name,...
%             num2str(dcAvg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).sourceIndex),...
%             num2str(dcAvg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).detectorIndex)));

        ylabel('[M * mm]');
        xlabel('Time [s]');
        %xlim([-2 7]);
        %ylim([-1e-3 1e-3]);
    end
    hold off;
    fn = sprintf('HRF_IndTrials_Raw_S%sD%s',...
        num2str(snirf1.data.measurementList(1,j).sourceIndex),...
        num2str(snirf1.data.measurementList(1,j).detectorIndex));
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
    close();
end
end
