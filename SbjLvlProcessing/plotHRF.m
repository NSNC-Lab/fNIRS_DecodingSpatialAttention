% plot HRF for each condition using output from Homer3 GUI

snirf1 = SnirfLoad('20210504_1348_02.snirf');
saveFigDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment03\Figures';
load('groupResults.mat');
%group.subjs(1).Load();
% Dimension = time x hbo/hbr/hbt x chn x stimclass
dcAvg = group.subjs(1,1).runs(1).procStream.output.dcAvg.GetDataTimeSeries('reshape');
tHRF = group.subjs(1,1).runs.procStream.output.dcAvg.GetTime();

numChns = size(dcAvg,3);
numConds = size(dcAvg,4);
numData = size(dcAvg,2);

% for i = 1:numChns
%     for j = 1:numConds
%         for k = 1:numData
for i = 1:1
    for j = 1:1
        for k = 1:1
            figure();
            plot(tHRF,squeeze(dcAvg(:,k,i,j)));
            srcStr = num2str(group.subjs(1,1).runs.procStream.output.dcAvg.measurementList(numConds*(j-1) + numChns*i + numData*k).sourceIndex);
            detectorStr = num2str(group.subjs(1,1).runs.procStream.output.dcAvg.measurementList(numConds*j + numChns*i + numData*k).detectorIndex);
            label = num2str(group.subjs(1,1).runs.procStream.output.dcAvg.measurementList(numConds*j + numChns*i + numData*k).dataTypeLabel);
            condStr = group.subjs(1,1).runs.procStream.output.misc.stim(1,group.subjs(1,1).runs.procStream.output.dcAvg.measurementList(numConds*j + numChns*i + numData*k).dataTypeIndex).name;
            titleStr = [condStr 'S' srcStr 'D' detectorStr ' ' label ];
            %title(titleStr);
            %title(sprintf('PSD Subject %s %s',num2str(sbjNum),chnsListPart2{i+1}));
            title(sprintf('%s S%sD%s %s',condStr,srcStr,detectorStr,label));
            
            filename = [saveFigDir filesep titleStr '.png'];
            %print(gcf,filename,'-dpng','-r250');
        end
    end
end

