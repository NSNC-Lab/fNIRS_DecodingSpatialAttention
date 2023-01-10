% first stim mark is start of cue
% second stim mark is start of movie
% third stim mark is start of question
% fourth stim mark is end of question

function plotHistISI(sbjNum)
    
rawDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
load([rawDataDir filesep sbjNum '.mat'],'s');
fName = s.fName;
snirf1 = SnirfClass(load([rawDataDir filesep fName{1} '.nirs'],'-mat'));

% for sbj 10 and earlier
if strcmp(sbjNum,'08') || strcmp(sbjNum,'10')
    
    %cueOnsetIndex = 1:4:1080;
    cueStartIndex = 1:4:713;
    cueEndIndex = 5:4:717;

    snirf2 = SnirfClass(load([rawDataDir filesep fName{2} '.nirs'],'-mat'));
    snirf3 = SnirfClass(load([rawDataDir filesep fName{3} '.nirs'],'-mat'));
    
    allS = find(snirf1.aux(1,1).dataTimeSeries>1);
    aInd1 = find(allS(2:end)-allS(1:end-1)==1);
    allS2 = find(snirf2.aux(1,1).dataTimeSeries>1);
    aInd2 = find(allS2(2:end)-allS2(1:end-1)==1);
    allS3 = find(snirf3.aux(1,1).dataTimeSeries>1);
    aInd3 = find(allS3(2:end)-allS3(1:end-1)==1);
    
    allS(aInd1) = [];
    allS = allS./50;
    allS2(aInd2) = [];
    allS2 = allS2./50;
    allS3(aInd3) = [];
    allS3 = allS3./50;
    
    snirf1TLen = size(snirf1.data.time,1)/(50);
    snirf2TLen = size(snirf2.data.time,1)/(50);
    allS2 = allS2+snirf1TLen;
    allS3 = allS3+snirf1TLen+snirf2TLen;

    allS = [allS; allS2; allS3];
    
    currT = allS(cueStartIndex);
    nextT = allS(cueEndIndex);
    
%     % plot histogram
%     histogram(endT-startT);
%     title(['Blank Screen Interval for Sbj: ' sbjNum]);
%     xlabel('s');
    
    % plot histogram
    histogram(nextT-currT);
    title(['ISI for Sbj: ' sbjNum]);
    xlabel('s');
    subtitle(sprintf('Mean %s Median %s SD %s',num2str(mean(nextT-currT)),...
        num2str(median(nextT-currT)),num2str(std(nextT-currT))));
% for sbj 12 and later
else
    %cueOnsetIndex = 1:4:1080;
    cueStartIndex = 1:4:353;
    cueEndIndex = 5:4:357;

    snirf1 = SnirfClass(load([fName{1} '.nirs'],'-mat'));

    allS = find(snirf1.aux(1,1).dataTimeSeries>1);
    aInd1 = find(allS(2:end)-allS(1:end-1)==1);
    % allS3 = find(snirf3.aux(1,1).dataTimeSeries>1);
    % aInd3 = find(allS3(2:end)-allS3(1:end-1)==1);
    % allS4 = find(snirf4.aux(1,1).dataTimeSeries>1);
    % aInd4 = find(allS4(2:end)-allS4(1:end-1)==1);
    allS(aInd1) = [];
    allS = allS./50;
    
    currT = allS(cueStartIndex);
    nextT = allS(cueEndIndex);
    
%     % plot histogram
%     histogram(endT-startT);
%     title(['Blank Screen for Sbj: ' sbjNum]);
%     xlabel('s');
    
    % plot histogram
    histogram(nextT-currT);
    title(['ISI for Sbj: ' sbjNum]);
    xlabel('s');
    subtitle(sprintf('Mean %s Median %s SD %s',num2str(mean(nextT-currT)),...
        num2str(median(nextT-currT)),num2str(std(nextT-currT))));
end