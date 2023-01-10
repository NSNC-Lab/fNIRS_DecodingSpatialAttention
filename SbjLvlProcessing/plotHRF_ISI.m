% find good subject and channel
% ISI bins: 27s, 27-30s, 30-33s, 33-36s, 36-39s, >39s
% For sbj 8 and 10

function plotHRF_ISI(sbjNum,rawDataFN,respData)

srcIdxGrp = {[1 2 3 4] [5 6] [8 9 10 11] [12 13 14 15] [16 17 18 19] [20 21] [23] [25] ...
    [27 28 29] [31 32 33] [35 36] [37 38] [40] [42]};

figSaveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum) '\Figures\HRFs'];

chnName = {'Middle Posterior FEF Left',...
    'Posterior to sPCS/tgPSC Left',...
    'Middle FEF Left',...
    'sPCS/tgPCS Left',...
    'Posterior FEF Left',...
    'Inferior to FEF Left',...
    'Anterior FEF Left',...
    'iPCS/tgPCS Left',...
    'Anterior to FEF Left',...
    'iPCS Left',...
    'Anterior FEF Right',...
    'iPCS/tgPCS Right',...
    'Anterior to FEF Right',...
    'iPCS Right',...
    'Middle Posterior FEF Right',...
    'Posterior to sPCS/tgPSC Right',...
    'Middle FEF Right',...
    'sPCS/tgPCS Right',...
    'Posterior to FEF Right',...
    'Inferor to FEF Right',...
    'Post Posterior STG/PT Left',...
    'Posterior STG/PT Right',...
    'IPS3/IPS2/SPL1 Left',...
    'IPS3/antIPS/IPS4 Left',...
    'IPS3/latIPS/antIPS Left',...
    'IPS3/IPS2/SPL1 Right',...
    'IPS3/antIPS/IPS4 Right',...
    'IPS3/latIPS/antIPS Right',...
    'Superior to IPS3/IPS2/SPL1 Left',...
    'IPS4 Left',...
    'Superior to IPS3/IPS2/SPL1 Right',...
    'IPS4 Right',...
    'Ant Posterior STG/PT Left',...
    'Ant Posterior STG/PT Right'};

saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
saveDir2 = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum) 'GUI'];
processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
% Convert nirs to snirf file format
snirf2 = SnirfClass(load([rawDataFN{1} '.nirs'],'-mat'));
snirf3 = SnirfClass(load([rawDataFN{2} '.nirs'],'-mat'));
snirf4 = SnirfClass(load([rawDataFN{3} '.nirs'],'-mat'));
snirf2.Info()

% Extract aux and convert to stimclass
allS2 = find(snirf2.aux(1,1).dataTimeSeries>1);
aInd2 = find(allS2(2:end)-allS2(1:end-1)==1);
allS3 = find(snirf3.aux(1,1).dataTimeSeries>1);
aInd3 = find(allS3(2:end)-allS3(1:end-1)==1);
allS4 = find(snirf4.aux(1,1).dataTimeSeries>1);
aInd4 = find(allS4(2:end)-allS4(1:end-1)==1);
allS2(aInd2) = [];
allS2 = allS2./50;
allS3(aInd3) = [];
allS3 = allS3./50;
allS4(aInd4) = [];
allS4 = allS4./50;

oldStimClass1Length2 = size(allS2,1);
oldStimClass1Length3 = oldStimClass1Length2 + size(allS3,1);
oldStimClass1Length4 = oldStimClass1Length3 + size(allS4,1);

%cueOnsetIndex = 1:4:1080;
cueOnsetIndex = 1:4:720;

cueStartIndex = 1:4:713;
cueEndIndex = 5:4:717;

cueOnsetIndex2 = 1:oldStimClass1Length2;
cueOnsetIndex3 = oldStimClass1Length2+1:oldStimClass1Length3;
cueOnsetIndex4 = oldStimClass1Length3+1:oldStimClass1Length4;

trialNum2 = sum(ismember(cueOnsetIndex2,cueOnsetIndex));
trialNum3 = sum(ismember(cueOnsetIndex3,cueOnsetIndex));
trialNum4 = sum(ismember(cueOnsetIndex4,cueOnsetIndex));

lastTrial2 = trialNum2;
lastTrial3 = lastTrial2 + trialNum3;
lastTrial4 = lastTrial3 + trialNum4;

snirf2TLen = size(snirf2.data.time,1)/(50);
snirf3TLen = size(snirf3.data.time,1)/(50);
allS3 = allS3+snirf2TLen;
allS4 = allS4+snirf2TLen+snirf3TLen;

allS = [allS2; allS3; allS4];
currT = allS(cueStartIndex);
nextT = allS(cueEndIndex);
ISI = (nextT-currT);
bin1 = ISI<28;
bin1 = [bin1; false];
%shortHRF = cueStartIndex(bin1);
bin2 = ISI>=28;
bin2 = [bin2; true];
%longHRF = cueStartIndex(bin2);

% split triggers into 6 categories
load([saveDir filesep respData '.mat'],'fixedMaskerList','indexMoviesTest','maskerMovies','numTrials','uniqueMovies');

% [8 4 12 6 2 10; 7 3 11 5 1 9]
condSubPlotIdx = [8 4 12 6 2 10; 7 3 11 5 1 9];
% orig idx: 
%condSubPlotIdx = [12 6 9 3 10 4 7 1 
condSubPlotIdx = [1 2 3 4 5 6; 7 8 9 10 11 12];
% legendTextShort = {'Left Single Short','Center Single Short','Right Single Short',...
%     'Left Multiple Short','Center Multiple Short','Right Multiple Short'};
% legendTextLong = {'Left Single Long','Center Single Long','Right Single Long',...
%     'Left Multiple Long','Center Multiple Long','Right Multiple Long'};
% snirf 2
% 1 left single short
% 2 right single short
% 3 center single short
% 4 left multiple short
% 5 right multiple short
% 6 center multiple short
% 7 left single long
% 8 right single long
% 9 center single long
% 10 left multiple long
% 11 right multiple long
% 12 center multiple long
leftOnsetSS = StimClass('leftSingleShort'); %8
rightOnsetSS = StimClass('rightSingleShort'); %12
centerOnsetSS = StimClass('centerSingleShort'); %4
leftOnsetMS = StimClass('leftMultiShort'); %6
rightOnsetMS = StimClass('rightMultiShort'); %10
centerOnsetMS = StimClass('centerMultiShort'); %2

leftOnsetSL = StimClass('leftSingleLong'); %7
rightOnsetSL = StimClass('rightSingleLong');% 11
centerOnsetSL = StimClass('centerSingleLong'); %3
leftOnsetML = StimClass('leftMultiLong'); %5
rightOnsetML = StimClass('rightMultiLong'); %9
centerOnsetML = StimClass('centerMultiLong'); %1

leftSingleIndex = indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==0;
rightSingleIndex = indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==0;
centerSingleIndex = indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==0;
leftMultiIndex = indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1;
rightMultiIndex = indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1;
centerMultiIndex = indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==1;

% leftSingleIndex = leftSingleIndex(leftSingleIndex<= lastTrial2);
% rightSingleIndex = rightSingleIndex(rightSingleIndex<= lastTrial2);
% centerSingleIndex = centerSingleIndex(centerSingleIndex<= lastTrial2);
% leftMultiIndex = leftMultiIndex(leftMultiIndex<= lastTrial2);
% rightMultiIndex = rightMultiIndex(rightMultiIndex<= lastTrial2);
% centerMultiIndex = centerMultiIndex(centerMultiIndex<= lastTrial2);

numTr = 0;

tempIdx = bin1.*leftSingleIndex;
numTr = numTr + length(find(tempIdx));
allSOnset = allS(cueOnsetIndex);
AddStims(leftOnsetSS, allSOnset(logical(tempIdx)));
tempIdx = bin1.*rightSingleIndex;
numTr = numTr + length(find(tempIdx));
AddStims(rightOnsetSS, allSOnset(logical(tempIdx)));
tempIdx = bin1.*centerSingleIndex;
numTr = numTr + length(find(tempIdx));
AddStims(centerOnsetSS, allSOnset(logical(tempIdx)));
tempIdx = bin1.*leftMultiIndex;
numTr = numTr + length(find(tempIdx));
AddStims(leftOnsetMS, allSOnset(logical(tempIdx)));
tempIdx = bin1.*rightMultiIndex;
numTr = numTr + length(find(tempIdx));
AddStims(rightOnsetMS, allSOnset(logical(tempIdx)));
tempIdx = bin1.*centerMultiIndex;
numTr = numTr + length(find(tempIdx));
AddStims(centerOnsetMS, allSOnset(logical(tempIdx)));

tempIdx = bin2.*leftSingleIndex;
numTr = numTr + length(find(tempIdx));
AddStims(leftOnsetSL, allSOnset(logical(tempIdx)));
tempIdx = bin2.*rightSingleIndex;
numTr = numTr + length(find(tempIdx));
AddStims(rightOnsetSL, allSOnset(logical(tempIdx)));
tempIdx = bin2.*centerSingleIndex;
numTr = numTr + length(find(tempIdx));
AddStims(centerOnsetSL, allSOnset(logical(tempIdx)));
tempIdx = bin2.*leftMultiIndex;
numTr = numTr + length(find(tempIdx));
AddStims(leftOnsetML, allSOnset(logical(tempIdx)));
tempIdx = bin2.*rightMultiIndex;
numTr = numTr + length(find(tempIdx));
AddStims(rightOnsetML, allSOnset(logical(tempIdx)));
tempIdx = bin2.*centerMultiIndex;
numTr = numTr + length(find(tempIdx));
AddStims(centerOnsetML, allSOnset(logical(tempIdx)));

updateStates(leftOnsetSS);
updateStates(rightOnsetSS);
updateStates(centerOnsetSS);
updateStates(leftOnsetMS);
updateStates(rightOnsetMS);
updateStates(centerOnsetMS);

updateStates(leftOnsetSL);
updateStates(rightOnsetSL);
updateStates(centerOnsetSL);
updateStates(leftOnsetML);
updateStates(rightOnsetML);
updateStates(centerOnsetML);

snirf2.stim(1,1) = leftOnsetSS;
snirf2.stim(1,2) = rightOnsetSS;
snirf2.stim(1,3) = centerOnsetSS;
snirf2.stim(1,4) = leftOnsetMS;
snirf2.stim(1,5) = rightOnsetMS;
snirf2.stim(1,6) = centerOnsetMS;

snirf2.stim(1,7) = leftOnsetSL;
snirf2.stim(1,8) = rightOnsetSL;
snirf2.stim(1,9) = centerOnsetSL;
snirf2.stim(1,10) = leftOnsetML;
snirf2.stim(1,11) = rightOnsetML;
snirf2.stim(1,12) = centerOnsetML;

% combine snirf2.data.dataTimeSeries and time
snirf2.data.dataTimeSeries = [snirf2.data.dataTimeSeries; snirf3.data.dataTimeSeries; snirf4.data.dataTimeSeries];
%snirf2.data.time = [snirf2.data.time;snirf3.data.time;snirf4.data.time];
snirf2.data.time = 0:1/50:(size(snirf2.data.dataTimeSeries,1)-1)/50;

data = snirf2.data;
probe = snirf2.probe;
mlActMan = {};
tIncMan = {};
%dod = snirf1.dod;
%mlActAuto = snirf1.mlActAuto;
stim = snirf2.stim;
%tIncAuto = snirf1.tIncAuto;
%dc = snirf1.dc;
Aaux = [];
rcMap = [];

mlActAuto = hmrR_PruneChannels(data,probe,mlActMan,tIncMan,[10000  10000000],2,[0  45]);

dod = hmrR_Intensity2OD(data);

%[dod,tInc,svs,nSV,tInc0] = hmrR_MotionCorrectPCArecurse(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,15,4,0.97,5,1);
[dod,tInc,svs,nSV,tInc0] = hmrR_MotionCorrectPCArecurse(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,20,5,0.97,5,1);

%tIncAuto = hmrR_MotionArtifact(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,20,4);
tIncAuto = hmrR_MotionArtifact(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,50,5);

[stim,tRange] = hmrR_StimRejection(dod,stim,tIncAuto,tIncMan,[-2  15]);

dod = hmrR_BandpassFilt(dod,0,0.5);
%dod = hmrR_BandpassFilt(dod,0.01,0.5);

dc = hmrR_OD2Conc(dod,probe,[1  1  1]);

[dcAvg,dcAvgStd,nTrials,dcNew,dcResid,dcSum2,beta,R,hmrstats] = hmrR_GLM(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
    [-2  15],1,1,[1.0 1.0 0.0 0.0 0.0 0.0],15,1,3,0);

snirf2.Save([saveDir2 filesep 'ISI.snirf']);
% fileName = [processedDataDir filesep 'intermediateOutputsISI.mat'];
% save(fileName,'dc','dcAvg','dcAvgStd','nTrials','dcNew',...
%     'dcResid','dcSum2','beta','R','hmrstats','dod','stim','tRange',...
%     'tIncAuto','mlActAuto','tIncAuto','snirf2','allS');

numConds = 12;
numConc = 3;
numChns = size(dcAvg.measurementList,2)/(numConds*numConc);
sizeCond = size(dcAvg.measurementList,2)/numConds;
sizeChn = sizeCond/numChns;
R = R{1};
chnNum = 0;

colorList = {[0, 0.447, 0.741],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],...
    [0.4940, 0.1840, 0.5560],	[0.4660, 0.6740, 0.1880] 	,[0.3010, 0.7450, 0.9330]};
lineStyle = {'-' '-' '-' '--' '--' '--'};
colorIdx = [1 3 2 1 3 2];

for i1 = 1:length(srcIdxGrp)
    for j = 1:length(srcIdxGrp{i1})
        chnNum = chnNum + 1;
        figure('units','normalized','outerposition',[0 0 1 1]); hold on;
        subplot(2,2,1);hold on;
        for i2 = 1:6
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),'-','Color',colorList{j-0});hold on;
            plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(1,i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),lineStyle{i2},'Color',colorList{colorIdx(i2)});hold on;
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2),':','Color',colorList{j-0});hold on;
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3),'-.','Color',colorList{j-0});hold on;
        end
        %legend(chnNameGrp{i1},'Location','northeastoutside');
        %legend({'Left Single','Center Single','Right Multiple','Right Single','Left Multiple','Center Multiple'},'Location','northeast');
        srcStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(1,i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).sourceIndex);
        detectorStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).detectorIndex);
        ylabel('[M*mm]')
        xlabel('Time [s]');
        xlim([-2 15]);
        %label = num2str(dcavg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).dataTypeLabel);
        %condStr = snirf1.stim(1,i2).name;
        title(sprintf('S%sD%s Short R: %s CAS: %s CAM: %s',...
            srcStr,detectorStr,num2str(R(srcIdxGrp{i1}(j),1))));
        subtitle(sprintf('HbO %s', chnName{chnNum}));
        
        subplot(2,2,2);hold on;
        for i2 = 1:6
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),'-','Color',colorList{j-0});hold on;
            plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(2,i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),lineStyle{i2},'Color',colorList{colorIdx(i2)});hold on;
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2),':','Color',colorList{j-0});hold on;
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3),'-.','Color',colorList{j-0});hold on;
        end
        %legend(chnNameGrp{i1},'Location','northeastoutside');
        %legend({'Left Single','Center Single','Right Multiple','Right Single','Left Multiple','Center Multiple'},'Location','northeast');
        srcStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(2,i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).sourceIndex);
        detectorStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).detectorIndex);
        ylabel('[M*mm]')
        xlabel('Time [s]');
        xlim([-2 15]);
        %label = num2str(dcavg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).dataTypeLabel);
        %condStr = snirf1.stim(1,i2).name;
        title(sprintf('S%sD%s Long R: %s CAS: %s CAM: %s',...
            srcStr,detectorStr,num2str(R(srcIdxGrp{i1}(j),2))));
        subtitle(sprintf('HbO %s',chnName{chnNum}));
        
        %title(sprintf('%s',condStr));
        %hold off;
%        fn = sprintf('HRF_S%sD%s',srcStr,detectorStr);
%        print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
%        close;
        
        subplot(2,2,3);hold on;
        for i2 = 1:6
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),'-','Color',colorList{j-0});hold on;
            plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(1,i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2),lineStyle{i2},'Color',colorList{colorIdx(i2)});hold on;
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2),':','Color',colorList{j-0});hold on;
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3),'-.','Color',colorList{j-0});hold on;
        end
        %legend(chnNameGrp{i1},'Location','northeastoutside');
        %legend({'Left Single','Center Single','Right Multiple','Right Single','Left Multiple','Center Multiple'},'Location','northeast');
        srcStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(1,i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).sourceIndex);
        detectorStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).detectorIndex);
        ylabel('[M*mm]')
        xlabel('Time [s]');
        xlim([-2 15]);
        %label = num2str(dcavg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).dataTypeLabel);
        %condStr = snirf1.stim(1,i2).name;
        title(sprintf('S%sD%s Short R: %s CAS: %s CAM: %s',...
            srcStr,detectorStr,num2str(R(srcIdxGrp{i1}(j),1))));
        subtitle(sprintf('HbR %s', chnName{chnNum}));
        
        subplot(2,2,4);hold on;
        for i2 = 1:6
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(i2-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),'-','Color',colorList{j-0});hold on;
            plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2),lineStyle{i2},'Color',colorList{colorIdx(i2)});hold on;
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2),':','Color',colorList{j-0});hold on;
            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 3),'-.','Color',colorList{j-0});hold on;
        end
        %legend(chnNameGrp{i1},'Location','northeastoutside');
        legend({'Left Single','Center Single','Right Single','Left Multiple','Center Multiple','Right Multiple'},'Location','northeast');
        srcStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).sourceIndex);
        detectorStr = num2str(dcAvg.measurementList(1,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1).detectorIndex);
        ylabel('[M*mm]')
        xlabel('Time [s]');
        xlim([-2 15]);
        %label = num2str(dcavg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).dataTypeLabel);
        %condStr = snirf1.stim(1,i2).name;
        title(sprintf('S%sD%s Long R: %s CAS: %s CAM: %s',...
            srcStr,detectorStr,num2str(R(srcIdxGrp{i1}(j),2))));
        subtitle(sprintf('HbR %s',chnName{chnNum}));
        
        %title(sprintf('%s',condStr));
        hold off;
        fn = sprintf('HRF_S%sD%s_ISI',srcStr,detectorStr);
        print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
        close;
    end
end











%dcReshape = dcNew.GetDataTimeSeries('reshape');




% tempIdx = bin1.*leftSingleIndex;
% for i = 1:size(allSOnset(logical(tempIdx)),1)
%     if (allSSingle(i)+15)*fs < size(dcReshape,1)
%         tempAr = allSOnset(logical(tempIdx));
%         trial_HbO_SS_Left(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,1,:))';
%         trial_HbR_SS_Left(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,2,:))';
%         trial_HbT_SS_Left(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,3,:))';
%     end
% end
% 
% tempIdx = bin1.*rightSingleIndex;
% for i = 1:size(allSOnset(logical(tempIdx)),1)
%     if (allSSingle(i)+15)*fs < size(dcReshape,1)
%         tempAr = allSOnset(logical(tempIdx));
%         trial_HbO_SS_Right(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,1,:))';
%         trial_HbR_SS_Right(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,2,:))';
%         trial_HbT_SS_Right(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,3,:))';
%     end
% end
% 
% tempIdx = bin1.*centerSingleIndex;
% for i = 1:size(allSOnset(logical(tempIdx)),1)
%     if (allSSingle(i)+15)*fs < size(dcReshape,1)
%         tempAr = allSOnset(logical(tempIdx));
%         trial_HbO_SS_Center(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,1,:))';
%         trial_HbR_SS_Center(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,2,:))';
%         trial_HbT_SS_Center(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,3,:))';
%     end
% end
% 
% tempIdx = bin1.*leftMultipleIndex;
% for i = 1:size(allSOnset(logical(tempIdx)),1)
%     if (allSSingle(i)+15)*fs < size(dcReshape,1)
%         tempAr = allSOnset(logical(tempIdx));
%         trial_HbO_MS_Left(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,1,:))';
%         trial_HbR_MS_Left(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,2,:))';
%         trial_HbT_MS_Left(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,3,:))';
%         
%     end
% end
% tempIdx = bin1.*rightMultipleIndex;
% for i = 1:size(allSOnset(logical(tempIdx)),1)
%     if (allSSingle(i)+15)*fs < size(dcReshape,1)
%         tempAr = allSOnset(logical(tempIdx));
%         trial_HbO_MS_Right(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,1,:))';
%         trial_HbR_MS_Right(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,2,:))';
%         trial_HbT_MS_Right(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,3,:))';
%     end
% end
% tempIdx = bin1.*centerMultipleIndex;
% for i = 1:size(allSOnset(logical(tempIdx)),1)
%     if (allSSingle(i)+15)*fs < size(dcReshape,1)
%         tempAr = allSOnset(logical(tempIdx));
%         trial_HbO_MS_Center(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,1,:))';
%         trial_HbR_MS_Center(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,2,:))';
%         trial_HbT_MS_Center(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,3,:))';
%     end
% end
%         
% 
% tempIdx = bin2.*leftSingleIndex;
% for i = 1:size(allSOnset(logical(tempIdx)),1)
%     if (allSSingle(i)+15)*fs < size(dcReshape,1)
%         tempAr = allSOnset(logical(tempIdx));
%         trial_HbO_SL_Left(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,1,:))';
%         trial_HbR_SL_Left(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,2,:))';
%         trial_HbT_SL_Left(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,3,:))';
%     end
% end
% 
% tempIdx = bin2.*rightSingleIndex;
% for i = 1:size(allSOnset(logical(tempIdx)),1)
%     if (allSSingle(i)+15)*fs < size(dcReshape,1)
%         tempAr = allSOnset(logical(tempIdx));
%         trial_HbO_SL_Right(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,1,:))';
%         trial_HbR_SL_Right(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,2,:))';
%         trial_HbT_SL_Right(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,3,:))';
%     end
% end
% 
% tempIdx = bin2.*centerSingleIndex;
% for i = 1:size(allSOnset(logical(tempIdx)),1)
%     if (allSSingle(i)+15)*fs < size(dcReshape,1)
%         tempAr = allSOnset(logical(tempIdx));
%         trial_HbO_SL_Center(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,1,:))';
%         trial_HbR_SL_Center(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,2,:))';
%         trial_HbT_SL_Center(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,3,:))';
%     end
% end
% 
% tempIdx = bin2.*leftMultipleIndex;
% for i = 1:size(allSOnset(logical(tempIdx)),1)
%     if (allSSingle(i)+15)*fs < size(dcReshape,1)
%         tempAr = allSOnset(logical(tempIdx));
%         trial_HbO_ML_Left(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,1,:))';
%         trial_HbR_ML_Left(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,2,:))';
%         trial_HbT_ML_Left(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,3,:))';
%         
%     end
% end
% tempIdx = bin2.*rightMultipleIndex;
% for i = 1:size(allSOnset(logical(tempIdx)),1)
%     if (allSSingle(i)+15)*fs < size(dcReshape,1)
%         tempAr = allSOnset(logical(tempIdx));
%         trial_HbO_ML_Right(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,1,:))';
%         trial_HbR_ML_Right(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,2,:))';
%         trial_HbT_ML_Right(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,3,:))';
%     end
% end
% tempIdx = bin2.*centerMultipleIndex;
% for i = 1:size(allSOnset(logical(tempIdx)),1)
%     if (allSSingle(i)+15)*fs < size(dcReshape,1)
%         tempAr = allSOnset(logical(tempIdx));
%         trial_HbO_ML_Center(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,1,:))';
%         trial_HbR_ML_Center(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,2,:))';
%         trial_HbT_ML_Center(:,:,i) = squeeze(dcReshape((fix(tempAr(i))-2)*fs:(fix(tempAr(i))+15)*fs,3,:))';
%     end
% end
        



