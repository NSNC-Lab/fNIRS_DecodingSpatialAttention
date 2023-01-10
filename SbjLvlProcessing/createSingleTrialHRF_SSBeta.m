% Extract single HRF from dcNew variable, the output of GLM fitting function
% save to mat file.
% For experiment 08 and 10

function createSingleTrialHRF_SSBeta(sbjNum,respData,numClasses)
%saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
load([saveDir filesep respData '.mat'],'indexMoviesTest');
load([processedDataDir filesep 'intermediateOutputsCombined_Basis1_SSBetaPrior.mat'],...
    'dcNew','allS');

fs  = 50;
% dcNew.dataTimeSeries: time x concentration x channels
dcReshape = dcNew.GetDataTimeSeries('reshape');
t = -2:1/fs:15;
tLen = length(t);

% get dc var from hmrR_OD2Conc func and split into ind trials
% cueOnsetIndex = 1:4:1080;
cueOnsetIndex = 1:4:720;
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
singleTrialHRFHbOS = zeros(size(dcReshape,3),size(t,2),size(allSSingle,1));
singleTrialHRFHbRS = zeros(size(dcReshape,3),size(t,2),size(allSSingle,1));
singleTrialHRFHbTS = zeros(size(dcReshape,3),size(t,2),size(allSSingle,1));

singleTrialHRFHbOM = zeros(size(dcReshape,3),size(t,2),size(allSMultiple,1));
singleTrialHRFHbRM = zeros(size(dcReshape,3),size(t,2),size(allSMultiple,1));
singleTrialHRFHbTM = zeros(size(dcReshape,3),size(t,2),size(allSMultiple,1));

for i = 1:size(allSSingle,1)
    if (allSSingle(i)+15)*fs < size(dcReshape,1)
%         singleTrialHRFHbOS(:,:,i) = squeeze(dcReshape((fix(allSSingle(i))-2)*fs:(fix(allSSingle(i))+15)*fs,1,:))';
%         singleTrialHRFHbRS(:,:,i) = squeeze(dcReshape((fix(allSSingle(i))-2)*fs:(fix(allSSingle(i))+15)*fs,2,:))';
%         singleTrialHRFHbTS(:,:,i) = squeeze(dcReshape((fix(allSSingle(i))-2)*fs:(fix(allSSingle(i))+15)*fs,3,:))';
        endT = fix((allSSingle(i)-2)*fs)+tLen;
        singleTrialHRFHbOS(:,:,i) = squeeze(dcReshape((fix((allSSingle(i)-2)*fs)):endT,1,:))';
        singleTrialHRFHbRS(:,:,i) = squeeze(dcReshape((fix((allSSingle(i)-2)*fs)):endT,2,:))';
        singleTrialHRFHbTS(:,:,i) = squeeze(dcReshape((fix((allSSingle(i)-2)*fs)):endT,3,:))';
    end
end

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

if numClasses == 2
    fn = [processedDataDir filesep 'singleTrialsUpdated_SSBeta_LR.mat'];
else
    fn = [processedDataDir filesep 'singleTrialsUpdated_SSBeta.mat'];
end
save(fn,'singleTrialHRFHbOM','singleTrialHRFHbRM','singleTrialHRFHbTM',...
    'singleTrialHRFHbOS','singleTrialHRFHbRS','singleTrialHRFHbTS','indexMoviesTest');

end
