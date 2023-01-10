% trials is channels x time x trials
% for continuous recording where data is time x (chn x cond x chromophore)

% for new SD, should cut down from 324 to 270 channels

function [data,varargout] = selectRS_cnt(data,isSDNew,mlList,mlActAuto)
    
numConds = 3;

if isSDNew
    % yes you can automate this list but too lazy
    ssIdx = [7,22,24,26,29,32];
else
    ssIdx = [7,22,24,26,30,34,39,41];
end
    
reshapeSSIdx = reshapeIdx(ssIdx,numConds,mlList);

data(:,reshapeSSIdx) = [];

switch nargin
    case 4
        ssIdxML = [ssIdx ssIdx+36];
        mlList(reshapeSSIdx) = [];
        mlActAuto(ssIdxML,:) = [];
        varargout{1} = mlList;
        varargout{2} = mlActAuto;
    otherwise
        varargout{1}= {};
end    

end