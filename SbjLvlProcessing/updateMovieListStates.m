% allS is trials x 1. doubles
% stim.state is 

function [movieList] = updateMovieListStates(movieList,rejTrOp,stim)
    
    movieList(:,7) = ones(1,size(movieList,1));
    
    if rejTrOp
        
        numStimClasses = length(stim);
    
        for stimClassI = 1:numStimClasses

            thisStimClass = stim(stimClassI);

            stimStates = thisStimClass.states;

            for stimI = 1:size(stimStates,1)
                if stimStates(stimI,2) < 0
                    idx = movieList(:,6)==stimStates(stimI,1);
                    movieList(idx,7) = stimStates(stimI,2);
                end
            end
        end
    
    end

end