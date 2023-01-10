% reshape index from (chn x cond x conc) to just chn
% mlList is list of all (chn x cond x cond) and can be found in dcAvg.measurementList
% given list of chn index, find all conditions and conc belonging to chn
% list.
% mlList increase conc first, next increase chn, next increase cond

function [reshapeIdx] = reshapeIdx(origIdx,numConds,mlList)

numConc = 3;
numChns = size(mlList,2)/(numConds*numConc);
sizeCond = size(mlList,2)/numConds;
sizeChn = sizeCond/numChns;
chromophoreIdx = [1 2 3];
reshapeIdx = zeros(length(chromophoreIdx)*length(origIdx)*numConds,1);
startCounter = 1;
for i1= 1:numConds
    for i2 = 1:length(origIdx)
    %         startIdx = (i1-1)*length(origIdx)+1;
    %         endIdx = i1*length(origIdx);
        startIdxML = sizeCond*(i1-1) + sizeChn.*(origIdx(i2)-1) + 1;
        endIdxML = sizeCond*(i1-1) + sizeChn.*(origIdx(i2)-1) + 3;
%         startIdx = sizeCond*(i1-1) + sizeChn.*((i2)-1) + 1;
%         endIdx = sizeCond*(i1-1) + sizeChn.*((i2)-1) + 3;
        reshapeIdx(startCounter:startCounter+2) = startIdxML:endIdxML;
        startCounter = startCounter+3;
    end
end

end