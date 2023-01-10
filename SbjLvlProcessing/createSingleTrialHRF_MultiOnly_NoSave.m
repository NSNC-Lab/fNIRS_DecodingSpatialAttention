% Extract single HRF from dcNew variable, the output of GLM fitting function
% save to mat file.
% For sbj 12 and after

% For allSMultiple, numClasses is predetermined.
%
% Fixed rounding error!

function [singleTrialHRFHbOM, singleTrialHRFHbRM, singleTrialHRFHbTM] ...
    = createSingleTrialHRF_MultiOnly_NoSave(trialsCont, allSMultiple)

fs  = 50;
% dcNew.dataTimeSeries: time x concentration x channels
dcReshape = trialsCont.GetDataTimeSeries('reshape');
%t = dcNew.GetT;
t = -2:1/fs:15;
tLen = length(t)-1;

% if strcmp(sbjNum,'15')
%     startTrial = 2;
%     subDim = 1;
% else
%     startTrial = 1;
%     subDim = 0;
% end

% channels x time x trial
singleTrialHRFHbOM = zeros(size(dcReshape,3),size(t,2),size(allSMultiple,1));
singleTrialHRFHbRM = zeros(size(dcReshape,3),size(t,2),size(allSMultiple,1));
singleTrialHRFHbTM = zeros(size(dcReshape,3),size(t,2),size(allSMultiple,1));

for i = 1:size(allSMultiple,1)
    if (allSMultiple(i)+15)*fs < size(dcReshape,1)
        endT = fix((allSMultiple(i)-2)*fs)+tLen;
        singleTrialHRFHbOM(:,:,i) = squeeze(dcReshape((fix((allSMultiple(i)-2)*fs)):endT,1,:))';
        singleTrialHRFHbRM(:,:,i) = squeeze(dcReshape((fix((allSMultiple(i)-2)*fs)):endT,2,:))';
        singleTrialHRFHbTM(:,:,i) = squeeze(dcReshape((fix((allSMultiple(i)-2)*fs)):endT,3,:))';

    end
end

end
