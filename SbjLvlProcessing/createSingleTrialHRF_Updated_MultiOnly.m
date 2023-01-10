% Extract single HRF from dcNew variable, the output of GLM fitting function
% save to mat file.
function createSingleTrialHRF_Updated_MultiOnly(sbjNum,movieList,numClasses)
%saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
load([saveDir filesep movieList '.mat'],'indexMoviesTest');
%load([processedDataDir filesep 'intermediateOutputsCombined.mat']);
load([processedDataDir filesep 'intermediateOutputsCombined_Basis1.mat']);

fs  = 50;
% dcNew.dataTimeSeries: time x concentration x channels
dcReshape = dcNew.GetDataTimeSeries('reshape');

% get dc var from hmrR_OD2Conc func and split into ind trials
cueOnsetIndex = 1:4:360;
if numClasses == 2
    %singleIndex = (indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==0) | (indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==0);
    multipleIndex = (indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1) | (indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1);
else
    %singleIndex = (indexMoviesTest(:,5)==0);
    multipleIndex = (indexMoviesTest(:,5)==1);
end
%allSSingle = allS(cueOnsetIndex(singleIndex));
allSMultiple = allS(cueOnsetIndex(multipleIndex));

chromXCond = 3*3;

singleTrialHRFHbOM = zeros(size(dcAvg.dataTimeSeries,2)/chromXCond,size(dcAvg.time,1),size(allSMultiple,1));
singleTrialHRFHbRM = zeros(size(dcAvg.dataTimeSeries,2)/chromXCond,size(dcAvg.time,1),size(allSMultiple,1));
singleTrialHRFHbTM = zeros(size(dcAvg.dataTimeSeries,2)/chromXCond,size(dcAvg.time,1),size(allSMultiple,1));

tLen = size(dcAvg.dataTimeSeries,1) - 1;

for i = 1:size(allSMultiple,1)
    if (allSMultiple(i)+15)*fs < size(dcReshape,1)
%         singleTrialHRFHbOM(:,:,i) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,1,:))';
%         singleTrialHRFHbRM(:,:,i) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,2,:))';
%         singleTrialHRFHbTM(:,:,i) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,3,:))';
        endT = fix((allSMultiple(i)-2)*fs)+tLen;
        singleTrialHRFHbOM(:,:,i) = squeeze(dcReshape((fix((allSMultiple(i)-2)*fs)):endT,1,:))';
        singleTrialHRFHbRM(:,:,i) = squeeze(dcReshape((fix((allSMultiple(i)-2)*fs)):endT,2,:))';
        singleTrialHRFHbTM(:,:,i) = squeeze(dcReshape((fix((allSMultiple(i)-2)*fs)):endT,3,:))';
    end
end


% for i = 1:length(stim)
%     % channels x time x trial
%     thisSingleTrialHRFHbO = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),length(squeeze(stim(1,i).data(:,1))));
%     thisSingleTrialHRFHbR = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),length(squeeze(stim(1,i).data(:,1))));
%     thisSingleTrialHRFHbT = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),length(squeeze(stim(1,i).data(:,1))));
%     disp(['i: ', num2str(i)]);
%     for j = 1:length(squeeze(stim(1,i).data(:,1)))
%         disp(['j: ',num2str(j)]);
%         
%         if (stim(1,i).data(j,1)+15)*fs < size(dcReshape,1)
%             thisSingleTrialHRFHbO(:,:,j) = squeeze(dcReshape((stim(1,i).data(j,1)-2)*fs:(stim(1,i).data(j,1)+15)*fs,1,:))';
%             thisSingleTrialHRFHbR(:,:,j) = squeeze(dcReshape((stim(1,i).data(j,1)-2)*fs:(stim(1,i).data(j,1)+15)*fs,2,:))';
%             thisSingleTrialHRFHbT(:,:,j) = squeeze(dcReshape((stim(1,i).data(j,1)-2)*fs:(stim(1,i).data(j,1)+15)*fs,3,:))';
%         end
%     end
%     singleTrialHRFCellHbO{i} = thisSingleTrialHRFHbO;
%     singleTrialHRFCellHbR{i} = thisSingleTrialHRFHbR;
%     singleTrialHRFCellHbT{i} = thisSingleTrialHRFHbT;
% end

if numClasses == 2
    fn = [processedDataDir filesep 'singleTrialsUpdated_LR_Basis1.mat'];
else
    fn = [processedDataDir filesep 'singleTrialsUpdated_Basis1.mat'];
end
save(fn,'singleTrialHRFHbOM','singleTrialHRFHbRM','singleTrialHRFHbTM','indexMoviesTest');

% % For actual analysis, plot 6 conditions from one channel in one figure
% colorList = {[0, 0.447, 0.741],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560] 	};
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

end
