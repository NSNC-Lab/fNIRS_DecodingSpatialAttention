% Should cut down from 756 to 378
% mlList/mlActAuto isn't affected

function [hrf_GLM,mlList] = rmvSingleCond(hrf_GLM,numConds,mlList)

condList = [1 3 5];

[reshapeIdxToRmv] = reshapeIdxByCond(condList,numConds,mlList);

hrf_GLM(:,reshapeIdxToRmv) = [];
mlList(:,reshapeIdxToRmv) = [];

end