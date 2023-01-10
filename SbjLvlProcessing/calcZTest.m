% beta is (#coefficients x HbX x #Channels x #conditions)
% bvar is #coefficients x # Channels x #HbX
% bvar is square of standard error, not SD
% sigAct is HbX x #Channels x #conditions
% sigDiff is HbX x #Channels x #pairs
% For sbj 12, # of coefficients is 16.
% use welch t-test
% Average beta 
% Use Wald's test for single parameter estimation.
% Modify GLM function to return bvar
% tSize = size(dcNew.dataTimeSeries,1);
% I'm unsure if the formulas are correct...

% Update to add statistical test for HbT using error propagation/Var(X+Y) =
% Var(X)+ Var(Y) when X and Y are independent, or more generally, Var(X+Y) =
% Var(X) + Var(Y) + 2*Cov(X,Y).

% Review stat here in Appendix B & C & D: NIRS-SPM: Statistical parametric mapping for near-infrared spectroscopy
% Jong Chul Ye, et al. 2009

function [sigAct,sigDiff,pvalAct,pvalDiff] = calcZTest(beta,bvar,tSize,isSDNew)

    numCoef = size(beta,1);
    numChrom = size(beta,2);
    numChns = size(beta,3);
    numConds = size(beta,4);
    numPairs = numConds*(numConds-1)/2;
    
    sigAct = zeros(numChrom,numChns,numConds);
    sigDiff = zeros(numChrom,numChns,numPairs);
    
    betaIdxRng = 3:8;
    % HbX x #Channels x #conditions)
    betaAvg = squeeze(mean(beta(betaIdxRng,:,:,:),1));
    betaAvgVar = zeros(size(betaAvg));
    
    if ~isSDNew
        %condList = [1;2;3;4;5;6];
        % only care about multi
        listPairs = [4 6;4 2;6 2];
        numPairs = 3;
        sigDiff = zeros(numChrom,numChns,numPairs);
    else
        condList = [1;2;3];
        listPairs = nchoosek(condList,2);
    end
    
    numBeta = length(betaIdxRng);
    
    pvalThreshold = 0.05;
    
    tval = zeros(size(beta));
    % numChrom x numChns x numConds
    pvalAct = zeros(numChrom,numChns,numConds);
    %pvalAct = zeros(size(beta));
    % numChrom x numChns x numPairs
    pvalDiff = zeros(numChrom,numChns,numPairs);

    % 2-tailed t-test. t-stat estimated using Wald test
    for i2 = 1:numChns
        for i3 = 1:numConds
            for i4 = 1:numChrom
                % uncertainty is measured in STD
                betaAvgVar(i4,i2,i3) = (1/(numBeta^2))*(sum(bvar((i3-1)*numCoef+betaIdxRng,i2,i4),1));
                tval(i4,i2,i3) =  betaAvg(i4,i2,i3)./sqrt(betaAvgVar(i4,i2,i3));
                % two-tailed t-test.
                %pval(i4,i2,i3) = (1-tcdf(abs(tval(i4,i2,i3)),tSize-1))*2;
                % given tSize is massive, pretty close to Gaussian
                pvalAct(i4,i2,i3) = tcdf(abs(tval(i4,i2,i3)),tSize-1,'upper')*2;
                % one-tailed t-test
                %pval(i1,i4,i2,i3) = tcdf(abs(tval(i1,i4,i2,i3)),tSize-1,'upper');
                sigAct(i4,i2,i3) = (pvalAct(i4,i2,i3) < pvalThreshold);
            end
        end
    end
    
    % z-test
    for i2 = 1:numChns
        for i3 = 1:size(listPairs,1)
            for i4 = 1:numChrom
                numerator = (betaAvg(i4,i2,listPairs(i3,1))-betaAvg(i4,i2,listPairs(i3,2)));
                denominator = sqrt(betaAvgVar(i4,i2,listPairs(i3,1))+betaAvgVar(i4,i2,listPairs(i3,2)));
                zScore = numerator/denominator;
                % two-tailed using standard normal distribution
                pvalDiff(i4,i2,i3) = 2*normcdf(-abs(zScore));
                sigDiff(i4,i2,i3) = pvalDiff(i4,i2,i3) < pvalThreshold;
            end
        end
    end
end