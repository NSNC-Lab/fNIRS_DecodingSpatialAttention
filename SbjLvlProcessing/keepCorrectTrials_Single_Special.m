% if trials var already filtered out trials by condition/locations, then
% just filter out movieIdx to same criteria.
% for fNIRS only
% For sbj 08 and 10

function [trials,varargout] = keepCorrectTrials_Single_Special(sbjNum,trials,behFN,numClasses,movieIdx)
    load(behFN,'responsesV','responsesA','correctRespV','correctRespA');
    %load(origIdxFN,'indexMoviesTest');
    
    % trials only have multiple condition and left & right loc
    if numClasses == 2
        singleIndex = (movieIdx(:,2)==2 & movieIdx(:,5)==0) | (movieIdx(:,2)==1 & movieIdx(:,5)==0);
    else
        singleIndex = movieIdx(:,5)==1;
    end
    
    if strcmp(sbjNum,'15')
        responsesA = responsesA(2:end);
        responsesV = responsesV(2:end);
        correctRespA = correctRespA(2:end);
        correctRespV = correctRespV(2:end);
    end
    
    responsesA = responsesA(singleIndex);
    responsesV = responsesV(singleIndex);
    correctRespA = correctRespA(singleIndex);
    correctRespV = correctRespV(singleIndex);
    movieIdxNew = movieIdx(singleIndex,:);
    %trials = trials(:,:,multipleIndex);
    
    idx = (responsesA==correctRespA).*(responsesV==correctRespV);
    
    trials = trials(:,:,logical(idx));
    
    behScore = sum((responsesA==correctRespA).*(responsesV==correctRespV))/length(responsesA);
    
    switch nargout
        case 2
            movieIdxNew = movieIdxNew(logical(idx),:);
            varargout{1} = movieIdxNew;
        case 3
            movieIdxNew = movieIdxNew(logical(idx),:);
            varargout{1} = movieIdxNew;
            varargout{2} = behScore;
        otherwise
            varargout{1} = [];
            varargout{2} = [];
    end

end