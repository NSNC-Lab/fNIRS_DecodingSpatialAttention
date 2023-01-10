% number of different conditions
numCond = 6;
% number of source-detector channels
numChn = 42;
% HbO, HbR, HbT
numConc = 3;
sizeOneCond = numChn*numConc;

procDataDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment08';
homer3GroupDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\GroupCAMapping\20210907';
load([procDataDir filesep 'performanceLinearDiscriminantUpdated_LR.mat']);
%load([homer3GroupDir filesep 'GroupCAMapping.mat']);
load([homer3GroupDir filesep '20210907.mat']);


fs = 50;
timePt = (0:0.5*fs:12*fs)+2*fs;
tInd = find(timePt==(5*fs+2*fs));

for i1 = 1:numChn
    for i2 = 1:numCond
        %trio = [(i1-1)*numConc+1 (i1-1)*numConc+2 (i1-1)*numConc+3] + (i2-1)*sizeOneCond;
        HbOIdx = (i1-1)*numConc+1 + (i2-1)*sizeOneCond;
        HbRIdx = (i1-1)*numConc+2 + (i2-1)*sizeOneCond;
        HbTIdx = (i1-1)*numConc+3 + (i2-1)*sizeOneCond;
        output.dcAvg.dataTimeSeries(:,HbOIdx) = performanceArrMultiHbO(i1,tInd);
        output.dcAvg.dataTimeSeries(:,HbRIdx) = performanceArrMultiHbR(i1,tInd);
        output.dcAvg.dataTimeSeries(:,HbTIdx) = performanceArrMultiHbT(i1,tInd);
    end
end

save([homer3GroupDir filesep '20210907.mat'],'output');
