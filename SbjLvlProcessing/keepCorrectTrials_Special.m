% if trials var already filtered out trials by condition/locations, then
% just filter out movieIdx to same criteria.
% for fNIRS only
% For sbj 08 and 10
% trials is already 90 trials, movieIdx has 180 trials.

function [trials,varargout] = keepCorrectTrials_Special(sbjNum,trials,behFN,numClasses,movieIdx)
    load(behFN,'responsesV','responsesA','correctRespV','correctRespA');
    %load(origIdxFN,'indexMoviesTest');
    
    if strcmp(sbjNum,'15')
        responsesA = responsesA(2:end);
        responsesV = responsesV(2:end);
        correctRespA = correctRespA(2:end);
        correctRespV = correctRespV(2:end);
        movieIdx = movieIdx(2:end,:);
        trials = trials(:,:,2:end);
    end
    
    % trials only have multiple condition and left & right loc
    if numClasses == 2
        multipleIndex = (movieIdx(:,2)==2 & movieIdx(:,5)==1) | (movieIdx(:,2)==1 & movieIdx(:,5)==1);
    else
        multipleIndex = movieIdx(:,5)==1;
    end
    
    responsesA = responsesA(multipleIndex);
    responsesV = responsesV(multipleIndex);
    correctRespA = correctRespA(multipleIndex);
    correctRespV = correctRespV(multipleIndex);
    movieIdxNew = movieIdx(multipleIndex,:);
 
    multipleIndex = (responsesA==correctRespA)' & (responsesV==correctRespV)';
    
    trials = trials(:,:,multipleIndex);
    movieIdxNew = movieIdxNew(multipleIndex,:);
    
    behScore = sum((responsesA==correctRespA)&(responsesV==correctRespV))/length(responsesA);
    
    switch nargout
        case 2
            %movieIdxNew = movieIdxNew(logical(idx),:);
            varargout{1} = movieIdxNew;
        case 3
            %movieIdxNew = movieIdxNew(logical(idx),:);
            varargout{1} = movieIdxNew;
            varargout{2} = behScore;
        otherwise
            varargout{1} = [];
            varargout{2} = [];
    end

end