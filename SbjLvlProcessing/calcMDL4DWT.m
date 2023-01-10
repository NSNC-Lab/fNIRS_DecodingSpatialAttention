% Compute the minimum description length
% In fNIRS case, we use symmetric extension.
% Signal content at border are unimportant, fortunately.

% X is design matrix.
% Y is fNIRS signal

function mdlScore = calcMDL4DWT(Y,X)
    % other choice is 'per'
    dwtmode('sym');
    

    

    mdlScore = 1;
    
    function codeLengthUniversalPrior = calcLU(n)
        c = 2.865064;

        codeLengthUniversalPrior = log2(c) + max(log2(n),0) + max(log2(log2(n)),0) + ...
            max(log2(log2(log2(n))),0) + ...
            max(log2(log2(log2(log2(n)))),0) + ...
            max(log2(log2(log2(log2(log2(n))))),0);
    end
    
    function modCodeLength = calcModL(j)
        
        mj = log2(j)-1;
        
        for n = 0:mj
            modCodeLength = 
        end
    
    end
end