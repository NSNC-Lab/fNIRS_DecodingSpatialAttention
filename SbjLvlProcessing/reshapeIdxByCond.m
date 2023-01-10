% reshape index from (chn x cond x conc) to just chn
% mlList is list of all (chn x cond x cond) and can be found in dcAvg.measurementList
% give list of conds, find all chn indx and conc belonging to cond list

function [reshapeIdx] = reshapeIdxByCond(condList,numConds,mlList)


numConc = 3;
numChns = size(mlList,2)/(numConds*numConc);
sizeCond = size(mlList,2)/numConds;
sizeChn = sizeCond/numChns;

reshapeIdx = zeros(sizeCond*length(condList),1);

for i = 1:length(condList)
    startIdx = (i-1)*sizeCond+1;
    endIdx  = i*sizeCond;
    reshapeIdx(startIdx:endIdx) = (condList(i)-1)*sizeCond+1:condList(i)*sizeCond;
end

end