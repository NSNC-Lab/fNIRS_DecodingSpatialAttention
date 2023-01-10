% trials is channels x time x trials
function [trials,varargout] = selectRS(trials,isSDNew,mlActAuto)
    

if isSDNew
    ssIdx = [7,22,24,26,29,32];
else
    ssIdx = [7,22,24,26,30,34,39,41];
end
    
trials(ssIdx,:,:) = [];

switch nargin
    case 3
        ssIdxML = [ssIdx ssIdx+36];
        mlActAuto(ssIdxML,:) = [];
        varargout{1} = mlActAuto;
    otherwise
        varargout{1}= {};
end    

end