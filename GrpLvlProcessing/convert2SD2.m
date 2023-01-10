function [trials,varargout] = convert2SD2(trials,mlActAuto)
    origIdxToRemove = [29,33,39,40,41,42];
    trials(origIdxToRemove,:,:) = [];

    switch nargin
        case 2
            origIdxToRemove = [origIdxToRemove origIdxToRemove+42];
            mlActAuto(origIdxToRemove,:) = [];
            varargout{1} = mlActAuto;
        otherwise
        varargout{1} = [];
    end

end