% Extract single HRF from dcNew variable, the output of GLM fitting function
% save to mat file.
% For experiment 08 and 10

% For allSSingle and allSMultiple, numClasses is predetermined.
% stimTr
%     stimTr(1,1) = leftOnsetSTr;
%     stimTr(1,2) = rightOnsetSTr;
%     stimTr(1,3) = centerOnsetSTr;
%     stimTr(1,4) = leftOnsetMTr;
%     stimTr(1,5) = rightOnsetMTr;
%     stimTr(1,6) = centerOnsetMTr;
%
% Fixed rounding error!

function [singleTrialHRFHbOS, singleTrialHRFHbRS,singleTrialHRFHbTS,...
    singleTrialHRFHbOM, singleTrialHRFHbRM, singleTrialHRFHbTM] = createSingleTrialHRF_NoSave(trialsCont, allSSingle, allSMultiple)

fs  = 50;
% dcNew.dataTimeSeries: time x concentration x channels
dcReshape = trialsCont.GetDataTimeSeries('reshape');
%t = trialsCont.GetT();
t = -2:1/fs:15;
tLen = length(t)-1;

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
        endT = fix((allSMultiple(i)-2)*fs)+tLen;
        singleTrialHRFHbOM(:,:,i) = squeeze(dcReshape((fix((allSMultiple(i)-2)*fs)):endT,1,:))';
        singleTrialHRFHbRM(:,:,i) = squeeze(dcReshape((fix((allSMultiple(i)-2)*fs)):endT,2,:))';
        singleTrialHRFHbTM(:,:,i) = squeeze(dcReshape((fix((allSMultiple(i)-2)*fs)):endT,3,:))';
    end
end

end
