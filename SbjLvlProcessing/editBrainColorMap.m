
sbjGrp = 1;

if sbjGrp == 1
    sbjGrpStr = 'HighPerf';
elseif sbjGrp == 2
    sbjGrpStr = 'All';
end

% For GrpAsSbj
% number of different conditions
numCond = 3;
% number of source-detector channels
numChn = 36;
numChnWOSS = 30;
% HbO, HbR, HbT
numConc = 3;
sizeOneCond = numChn*numConc;


procDataDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment08';
dataSaveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll\';
homer3GroupDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\GroupCA_20220114\GrpAsSbj';
% grpBarChartROINames_04.m
%load([dataSaveDir filesep 'grpIndChnCA_' sbjGrpStr],'varPerfHbO','varPerfHbR','varPerfHbT');

% grpBarChartROINames_RCoeff
%load([dataSaveDir filesep 'grpIndChnR_' sbjGrpStr],'rValHbT');
load([dataSaveDir filesep 'grpIndChnCA_' sbjGrpStr],'varPerfHbT');
%load([homer3GroupDir filesep 'GroupCAMapping.mat']);
load([homer3GroupDir filesep 'GrpAsSbj.mat']);
% this is run lvl: 20220427_1251_03.mat
% this is sbj lvl: Subj_20220427_1251_03.mat
% this is grp lvl: GroupCA_20220114.mat
% these will load ProcResultClass

 % ssIdx = [7,22,24,26,29,32];
chnNumConv = [1 2 3 4 5 6 8 9 10 11 12 13 14 15 16 17 18 19 20 21 23 25 27 ...
    28 30 31 33 34 35 36];
chnNumSSConv = [7 22 24 26 29 32];

fs = 50;
timePt = (0:0.5*fs:12*fs)+2*fs;
tInd = find(timePt==(5*fs+2*fs));

for i1 = 1:numChnWOSS
    chnNumWSS = chnNumConv(i1);
    for i2 = 1:numCond
        %trio = [(i1-1)*numConc+1 (i1-1)*numConc+2 (i1-1)*numConc+3] + (i2-1)*sizeOneCond;
        HbOIdx = (chnNumWSS-1)*numConc+1 + (i2-1)*sizeOneCond;
        HbRIdx = (chnNumWSS-1)*numConc+2 + (i2-1)*sizeOneCond;
        HbTIdx = (chnNumWSS-1)*numConc+3 + (i2-1)*sizeOneCond;
        % varPerfHbO only has 30 channels
        output.dcAvg.dataTimeSeries(:,HbOIdx) = max(0.43,varPerfHbT(i1));
        output.dcAvg.dataTimeSeries(:,HbRIdx) = max(0.43,varPerfHbT(i1));
        output.dcAvg.dataTimeSeries(:,HbTIdx) = max(0.43,varPerfHbT(i1));
    end
end

for i1 = 1:length(chnNumSSConv)
    chnNumSS = chnNumSSConv(i1);
    for i2 = 1:numCond
        HbOIdx = (chnNumSS-1)*numConc+1 + (i2-1)*sizeOneCond;
        HbRIdx = (chnNumSS-1)*numConc+2 + (i2-1)*sizeOneCond;
        HbTIdx = (chnNumSS-1)*numConc+3 + (i2-1)*sizeOneCond;
        output.dcAvg.dataTimeSeries(:,HbOIdx) = 0.43;
        output.dcAvg.dataTimeSeries(:,HbRIdx) = 0.43;
        output.dcAvg.dataTimeSeries(:,HbTIdx) = 0.43;
    end
end

save([homer3GroupDir filesep 'GrpAsSbj.mat'],'output');
