% Extract single HRF from dcNew variable, the output of GLM fitting function
% save to mat file.
% All 4 basis, for sbj 12 and after. For 08 and 10, use createSingleTrialHRF_Updated_DiffBasis(sbjNum,respData,numClasses)
% For subject 15, exclude first trial because first trial starts at 0.08s,
% cannot get -2s pre-onset
function createSingleTrialHRF_Updated_MultiOnly_DiffBasis(sbjNum,movieList,opAllBasis)
%saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
load([saveDir filesep movieList '.mat'],'indexMoviesTest');

if strcmp(sbjNum,'15')
    startTrial = 2;
    subDim = 1;
else
    startTrial = 1;
    subDim = 0;
end

% Gaussian Basis: #1
load([processedDataDir filesep 'intermediateOutputsCombined_Basis1.mat'],'dcNew','dcAvg','allS');

fs  = 50;
% dcNew.dataTimeSeries: time x concentration x channels
dcReshape = dcNew.GetDataTimeSeries('reshape');

% get dc var from hmrR_OD2Conc func and split into ind trials
cueOnsetIndex = 1:4:360;
allSMultiple = allS(cueOnsetIndex);

% channels x time x trial
% if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
%     chromXCond = 3*6;
% else
%     chromXCond = 3*3;
% end
chromXCond = 3*3;

singleTrialHRFHbOM = zeros(size(dcAvg.dataTimeSeries,2)/chromXCond,size(dcAvg.time,1),size(allSMultiple,1)-subDim);
singleTrialHRFHbRM = zeros(size(dcAvg.dataTimeSeries,2)/chromXCond,size(dcAvg.time,1),size(allSMultiple,1)-subDim);
singleTrialHRFHbTM = zeros(size(dcAvg.dataTimeSeries,2)/chromXCond,size(dcAvg.time,1),size(allSMultiple,1)-subDim);

for i = startTrial:size(allSMultiple,1)
    if (allSMultiple(i)+15)*fs < size(dcReshape,1)
        singleTrialHRFHbOM(:,:,i-subDim) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,1,:))';
        singleTrialHRFHbRM(:,:,i-subDim) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,2,:))';
        singleTrialHRFHbTM(:,:,i-subDim) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,3,:))';
    end
end

fn = [processedDataDir filesep 'singleTrialsUpdated_Basis1.mat'];

save(fn,'singleTrialHRFHbOM','singleTrialHRFHbRM','singleTrialHRFHbTM','indexMoviesTest');

if opAllBasis == 1

    % modified gamma function convolved with a square-wave of duration T. #2
    load([processedDataDir filesep 'intermediateOutputsCombined_Basis2.mat']);

    fs  = 50;
    % dcNew.dataTimeSeries: time x concentration x channels
    dcReshape = dcNew.GetDataTimeSeries('reshape');

    % get dc var from hmrR_OD2Conc func and split into ind trials
    cueOnsetIndex = 1:4:360;
    allSMultiple = allS(cueOnsetIndex);

    % % channels x time x trial
    % if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
    %     chromXCond = 3*6;
    % else
    %     chromXCond = 3*3;
    % end
    chromXCond = 3*3;

    singleTrialHRFHbOM = zeros(size(dcAvg.dataTimeSeries,2)/chromXCond,size(dcAvg.time,1),size(allSMultiple,1)-subDim);
    singleTrialHRFHbRM = zeros(size(dcAvg.dataTimeSeries,2)/chromXCond,size(dcAvg.time,1),size(allSMultiple,1)-subDim);
    singleTrialHRFHbTM = zeros(size(dcAvg.dataTimeSeries,2)/chromXCond,size(dcAvg.time,1),size(allSMultiple,1)-subDim);

    for i = startTrial:size(allSMultiple,1)
        if (allSMultiple(i)+15)*fs < size(dcReshape,1)
            singleTrialHRFHbOM(:,:,i-subDim) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,1,:))';
            singleTrialHRFHbRM(:,:,i-subDim) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,2,:))';
            singleTrialHRFHbTM(:,:,i-subDim) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,3,:))';
        end
    end

    fn = [processedDataDir filesep 'singleTrialsUpdated_Basis2.mat'];

    save(fn,'singleTrialHRFHbOM','singleTrialHRFHbRM','singleTrialHRFHbTM','indexMoviesTest');

    % a modified gamma function and its derivative convolved with a square-wave of duration T.  #3
    load([processedDataDir filesep 'intermediateOutputsCombined_Basis3.mat']);

    fs  = 50;
    % dcNew.dataTimeSeries: time x concentration x channels
    dcReshape = dcNew.GetDataTimeSeries('reshape');

    % get dc var from hmrR_OD2Conc func and split into ind trials
    cueOnsetIndex = 1:4:360;
    allSMultiple = allS(cueOnsetIndex);

    % % channels x time x trial
    % if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
    %     chromXCond = 3*6;
    % else
    %     chromXCond = 3*3;
    % end
    chromXCond = 3*3;

    singleTrialHRFHbOM = zeros(size(dcAvg.dataTimeSeries,2)/chromXCond,size(dcAvg.time,1),size(allSMultiple,1)-subDim);
    singleTrialHRFHbRM = zeros(size(dcAvg.dataTimeSeries,2)/chromXCond,size(dcAvg.time,1),size(allSMultiple,1)-subDim);
    singleTrialHRFHbTM = zeros(size(dcAvg.dataTimeSeries,2)/chromXCond,size(dcAvg.time,1),size(allSMultiple,1)-subDim);

    for i = startTrial:size(allSMultiple,1)
        if (allSMultiple(i)+15)*fs < size(dcReshape,1)
            singleTrialHRFHbOM(:,:,i-subDim) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,1,:))';
            singleTrialHRFHbRM(:,:,i-subDim) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,2,:))';
            singleTrialHRFHbTM(:,:,i-subDim) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,3,:))';
        end
    end

    fn = [processedDataDir filesep 'singleTrialsUpdated_Basis3.mat'];

    save(fn,'singleTrialHRFHbOM','singleTrialHRFHbRM','singleTrialHRFHbTM','indexMoviesTest');

    %GAM function from 3dDeconvolve AFNI convolved with a square-wave of
    %duration T. #4
    load([processedDataDir filesep 'intermediateOutputsCombined_Basis4.mat']);

    fs  = 50;
    % dcNew.dataTimeSeries: time x concentration x channels
    dcReshape = dcNew.GetDataTimeSeries('reshape');

    % get dc var from hmrR_OD2Conc func and split into ind trials
    cueOnsetIndex = 1:4:360;
    allSMultiple = allS(cueOnsetIndex);

    % % channels x time x trial
    % if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
    %     chromXCond = 3*6;
    % else
    %     chromXCond = 3*3;
    % end
    chromXCond = 3*3;

    singleTrialHRFHbOM = zeros(size(dcAvg.dataTimeSeries,2)/chromXCond,size(dcAvg.time,1),size(allSMultiple,1)-subDim);
    singleTrialHRFHbRM = zeros(size(dcAvg.dataTimeSeries,2)/chromXCond,size(dcAvg.time,1),size(allSMultiple,1)-subDim);
    singleTrialHRFHbTM = zeros(size(dcAvg.dataTimeSeries,2)/chromXCond,size(dcAvg.time,1),size(allSMultiple,1)-subDim);

    for i = startTrial:size(allSMultiple,1)
        if (allSMultiple(i)+15)*fs < size(dcReshape,1)
            singleTrialHRFHbOM(:,:,i-subDim) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,1,:))';
            singleTrialHRFHbRM(:,:,i-subDim) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,2,:))';
            singleTrialHRFHbTM(:,:,i-subDim) = squeeze(dcReshape((fix(allSMultiple(i))-2)*fs:(fix(allSMultiple(i))+15)*fs,3,:))';
        end
    end

    fn = [processedDataDir filesep 'singleTrialsUpdated_Basis4.mat'];

    save(fn,'singleTrialHRFHbOM','singleTrialHRFHbRM','singleTrialHRFHbTM','indexMoviesTest');

end

end
