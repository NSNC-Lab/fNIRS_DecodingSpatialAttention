% for continuous recording where data is time x (chn x cond x chromophore)
% this is only for sbj 08 and 10. To convert old SD to new SD
% should cut down from 378 chns to 324 chns

function [data,varargout] = convert2SD2_cnt(data,mlList,mlActAuto)

numConds = 3;

% list of both LS and SS chns to remove, 42 to 36 chns
origIdxToRemove = [29,33,39,40,41,42];

reshapeIdxToRmv = reshapeIdx(origIdxToRemove,numConds,mlList);

data(:,reshapeIdxToRmv) = [];

switch nargin
    case 3
        origIdxToRemove = [origIdxToRemove origIdxToRemove+42];
        mlList(reshapeIdxToRmv) = [];
        mlActAuto(origIdxToRemove,:) = [];
        varargout{1} = mlList;
        varargout{2} = mlActAuto;
    otherwise
    varargout{1} = [];
end


end