% Extract single HRF from dcNew variable, the output of GLM fitting function
% save to mat file.
% For experiment 03 and 06
function createSingleTrialHRF_Updated_DiffBasis_Sbj03(sbjNum,respData,numClasses)
%saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
load([saveDir filesep respData '.mat'],'fixedMaskerList','indexMoviesTest','maskerMovies','numTrials','uniqueMovies');

% Gaussian Basis: #1
load([processedDataDir filesep 'intermediateOutputsCombined_Basis1.mat']);

fs  = 50;
% dcNew.dataTimeSeries: time x concentration x channels
dcReshape = dcNew.GetDataTimeSeries('reshape');

% get dc var from hmrR_OD2Conc func and split into ind trials
cueOnsetIndex = 1:4:1080;
idxCutoff = find(cueOnsetIndex>length(allS));
cueOnsetIndex(idxCutoff) = [];

oldStimClass1Length2 = size(allS,1);
cueOnsetIndex2 = 1:oldStimClass1Length2;
trialNum2 = sum(ismember(cueOnsetIndex2,cueOnsetIndex));
indexMoviesTest = indexMoviesTest(1:trialNum2,:);

if numClasses == 2
    singleIndex = (indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==0) | (indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==0);
    multipleIndex = (indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1) | (indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1);
else
    singleIndex = (indexMoviesTest(:,5)==0);
    multipleIndex = (indexMoviesTest(:,5)==1);
end
allSSingle = allS(cueOnsetIndex(singleIndex));
allSMultiple = allS(cueOnsetIndex(multipleIndex));

% channels x time x trial
singleTrialHRFHbOS = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSSingle,1));
singleTrialHRFHbRS = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSSingle,1));
singleTrialHRFHbTS = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSSingle,1));

singleTrialHRFHbOM = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSMultiple,1));
singleTrialHRFHbRM = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSMultiple,1));
singleTrialHRFHbTM = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSMultiple,1));

for i = 1:size(allSSingle,1)
    if (allSSingle(i)+15)*fs < size(dcReshape,1)
        singleTrialHRFHbOS(:,:,i) = squeeze(dcReshape((fix(allSSingle(i))-2)*fs:(fix(allSSingle(i))+15)*fs,1,:))';
        singleTrialHRFHbRS(:,:,i) = squeeze(dcReshape((fix(allSSingle(i))-2)*fs:(fix(allSSingle(i))+15)*fs,2,:))';
        singleTrialHRFHbTS(:,:,i) = squeeze(dcReshape((fix(allSSingle(i))-2)*fs:(fix(allSSingle(i))+15)*fs,3,:))';
    end
end

for i = 1:size(allSMultiple,1)
    if (allSMultiple(i)+15)*fs < size(dcReshape,1)
        singleTrialHRFHbOM(:,:,i) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,1,:))';
        singleTrialHRFHbRM(:,:,i) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,2,:))';
        singleTrialHRFHbTM(:,:,i) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,3,:))';
    end
end

if numClasses == 2
    fn = [processedDataDir filesep 'singleTrialsUpdated_LR_Basis1.mat'];
else
    fn = [processedDataDir filesep 'singleTrialsUpdated_Basis1.mat'];
end
save(fn,'singleTrialHRFHbOS','singleTrialHRFHbRS','singleTrialHRFHbTS',...
    'singleTrialHRFHbOM','singleTrialHRFHbRM','singleTrialHRFHbTM','indexMoviesTest');

% Basis 2

% Gaussian Basis: #1
load([processedDataDir filesep 'intermediateOutputsCombined_Basis2.mat']);

fs  = 50;
% dcNew.dataTimeSeries: time x concentration x channels
dcReshape = dcNew.GetDataTimeSeries('reshape');

% get dc var from hmrR_OD2Conc func and split into ind trials
cueOnsetIndex = 1:4:1080;
idxCutoff = find(cueOnsetIndex>length(allS));
cueOnsetIndex(idxCutoff) = [];

oldStimClass1Length2 = size(allS,1);
cueOnsetIndex2 = 1:oldStimClass1Length2;
trialNum2 = sum(ismember(cueOnsetIndex2,cueOnsetIndex));
indexMoviesTest = indexMoviesTest(1:trialNum2,:);

if numClasses == 2
    singleIndex = (indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==0) | (indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==0);
    multipleIndex = (indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1) | (indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1);
else
    singleIndex = (indexMoviesTest(:,5)==0);
    multipleIndex = (indexMoviesTest(:,5)==1);
end
allSSingle = allS(cueOnsetIndex(singleIndex));
allSMultiple = allS(cueOnsetIndex(multipleIndex));

% channels x time x trial
singleTrialHRFHbOS = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSSingle,1));
singleTrialHRFHbRS = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSSingle,1));
singleTrialHRFHbTS = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSSingle,1));

singleTrialHRFHbOM = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSMultiple,1));
singleTrialHRFHbRM = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSMultiple,1));
singleTrialHRFHbTM = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSMultiple,1));

for i = 1:size(allSSingle,1)
    if (allSSingle(i)+15)*fs < size(dcReshape,1)
        singleTrialHRFHbOS(:,:,i) = squeeze(dcReshape((fix(allSSingle(i))-2)*fs:(fix(allSSingle(i))+15)*fs,1,:))';
        singleTrialHRFHbRS(:,:,i) = squeeze(dcReshape((fix(allSSingle(i))-2)*fs:(fix(allSSingle(i))+15)*fs,2,:))';
        singleTrialHRFHbTS(:,:,i) = squeeze(dcReshape((fix(allSSingle(i))-2)*fs:(fix(allSSingle(i))+15)*fs,3,:))';
    end
end

for i = 1:size(allSMultiple,1)
    if (allSMultiple(i)+15)*fs < size(dcReshape,1)
        singleTrialHRFHbOM(:,:,i) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,1,:))';
        singleTrialHRFHbRM(:,:,i) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,2,:))';
        singleTrialHRFHbTM(:,:,i) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,3,:))';
    end
end

if numClasses == 2
    fn = [processedDataDir filesep 'singleTrialsUpdated_LR_Basis2.mat'];
else
    fn = [processedDataDir filesep 'singleTrialsUpdated_Basis2.mat'];
end
save(fn,'singleTrialHRFHbOS','singleTrialHRFHbRS','singleTrialHRFHbTS',...
    'singleTrialHRFHbOM','singleTrialHRFHbRM','singleTrialHRFHbTM','indexMoviesTest');


% Basis 3
% Gaussian Basis: #1
load([processedDataDir filesep 'intermediateOutputsCombined_Basis3.mat']);

fs  = 50;
% dcNew.dataTimeSeries: time x concentration x channels
dcReshape = dcNew.GetDataTimeSeries('reshape');

% get dc var from hmrR_OD2Conc func and split into ind trials
cueOnsetIndex = 1:4:1080;
idxCutoff = find(cueOnsetIndex>length(allS));
cueOnsetIndex(idxCutoff) = [];

oldStimClass1Length2 = size(allS,1);
cueOnsetIndex2 = 1:oldStimClass1Length2;
trialNum2 = sum(ismember(cueOnsetIndex2,cueOnsetIndex));
indexMoviesTest = indexMoviesTest(1:trialNum2,:);

if numClasses == 2
    singleIndex = (indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==0) | (indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==0);
    multipleIndex = (indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1) | (indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1);
else
    singleIndex = (indexMoviesTest(:,5)==0);
    multipleIndex = (indexMoviesTest(:,5)==1);
end
allSSingle = allS(cueOnsetIndex(singleIndex));
allSMultiple = allS(cueOnsetIndex(multipleIndex));

% channels x time x trial
singleTrialHRFHbOS = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSSingle,1));
singleTrialHRFHbRS = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSSingle,1));
singleTrialHRFHbTS = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSSingle,1));

singleTrialHRFHbOM = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSMultiple,1));
singleTrialHRFHbRM = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSMultiple,1));
singleTrialHRFHbTM = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSMultiple,1));

for i = 1:size(allSSingle,1)
    if (allSSingle(i)+15)*fs < size(dcReshape,1)
        singleTrialHRFHbOS(:,:,i) = squeeze(dcReshape((fix(allSSingle(i))-2)*fs:(fix(allSSingle(i))+15)*fs,1,:))';
        singleTrialHRFHbRS(:,:,i) = squeeze(dcReshape((fix(allSSingle(i))-2)*fs:(fix(allSSingle(i))+15)*fs,2,:))';
        singleTrialHRFHbTS(:,:,i) = squeeze(dcReshape((fix(allSSingle(i))-2)*fs:(fix(allSSingle(i))+15)*fs,3,:))';
    end
end

for i = 1:size(allSMultiple,1)
    if (allSMultiple(i)+15)*fs < size(dcReshape,1)
        singleTrialHRFHbOM(:,:,i) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,1,:))';
        singleTrialHRFHbRM(:,:,i) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,2,:))';
        singleTrialHRFHbTM(:,:,i) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,3,:))';
    end
end

if numClasses == 2
    fn = [processedDataDir filesep 'singleTrialsUpdated_LR_Basis3.mat'];
else
    fn = [processedDataDir filesep 'singleTrialsUpdated_Basis3.mat'];
end
save(fn,'singleTrialHRFHbOS','singleTrialHRFHbRS','singleTrialHRFHbTS',...
    'singleTrialHRFHbOM','singleTrialHRFHbRM','singleTrialHRFHbTM','indexMoviesTest');


% Basis 4

% Gaussian Basis: #1
load([processedDataDir filesep 'intermediateOutputsCombined_Basis4.mat']);

fs  = 50;
% dcNew.dataTimeSeries: time x concentration x channels
dcReshape = dcNew.GetDataTimeSeries('reshape');

% get dc var from hmrR_OD2Conc func and split into ind trials
cueOnsetIndex = 1:4:1080;
idxCutoff = find(cueOnsetIndex>length(allS));
cueOnsetIndex(idxCutoff) = [];

oldStimClass1Length2 = size(allS,1);
cueOnsetIndex2 = 1:oldStimClass1Length2;
trialNum2 = sum(ismember(cueOnsetIndex2,cueOnsetIndex));
indexMoviesTest = indexMoviesTest(1:trialNum2,:);

if numClasses == 2
    singleIndex = (indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==0) | (indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==0);
    multipleIndex = (indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1) | (indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1);
else
    singleIndex = (indexMoviesTest(:,5)==0);
    multipleIndex = (indexMoviesTest(:,5)==1);
end
allSSingle = allS(cueOnsetIndex(singleIndex));
allSMultiple = allS(cueOnsetIndex(multipleIndex));

% channels x time x trial
singleTrialHRFHbOS = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSSingle,1));
singleTrialHRFHbRS = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSSingle,1));
singleTrialHRFHbTS = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSSingle,1));

singleTrialHRFHbOM = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSMultiple,1));
singleTrialHRFHbRM = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSMultiple,1));
singleTrialHRFHbTM = zeros(size(dcAvg.dataTimeSeries,2)/18,size(dcAvg.time,1),size(allSMultiple,1));

for i = 1:size(allSSingle,1)
    if (allSSingle(i)+15)*fs < size(dcReshape,1)
        singleTrialHRFHbOS(:,:,i) = squeeze(dcReshape((fix(allSSingle(i))-2)*fs:(fix(allSSingle(i))+15)*fs,1,:))';
        singleTrialHRFHbRS(:,:,i) = squeeze(dcReshape((fix(allSSingle(i))-2)*fs:(fix(allSSingle(i))+15)*fs,2,:))';
        singleTrialHRFHbTS(:,:,i) = squeeze(dcReshape((fix(allSSingle(i))-2)*fs:(fix(allSSingle(i))+15)*fs,3,:))';
    end
end

for i = 1:size(allSMultiple,1)
    if (allSMultiple(i)+15)*fs < size(dcReshape,1)
        singleTrialHRFHbOM(:,:,i) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,1,:))';
        singleTrialHRFHbRM(:,:,i) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,2,:))';
        singleTrialHRFHbTM(:,:,i) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,3,:))';
    end
end

if numClasses == 2
    fn = [processedDataDir filesep 'singleTrialsUpdated_LR_Basis4.mat'];
else
    fn = [processedDataDir filesep 'singleTrialsUpdated_Basis4.mat'];
end
save(fn,'singleTrialHRFHbOS','singleTrialHRFHbRS','singleTrialHRFHbTS',...
    'singleTrialHRFHbOM','singleTrialHRFHbRM','singleTrialHRFHbTM','indexMoviesTest');


end
