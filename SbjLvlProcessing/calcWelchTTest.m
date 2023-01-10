% beta is #sbj x HbX x #Channels x #conditions
% bvar is #sbj x #HxB x #Channels x #conditions
% sigAct is HbX x #Channels x #conditions
% sigDiff is HbX x #Channels x #pairs 

% For sbj 12, # of coefficients is 16.
% use welch t-test
% Average beta 
% Use Wald's test for single parameter estimation.
% Given each observation, each with uncertainty, you can do weighted t-test
% However, there's barely any Google-able resources that cover it and
% Wikipedia is nada on it.

% Helm-Bonferonni correction for 3 comparisons.

function [sigAct,sigDiff] = calcWelchTTest(betaVar,bvar)

    numSbj = size(betaVar,1);
    numChrom = size(betaVar,2);
    numChns = size(betaVar,3);
    numConds = size(betaVar,4);
    numPairs = numConds*(numConds-1)/2;
    m = 3;
    
    sigAct = zeros(numChrom,numChns,numConds);
    sigDiff = zeros(numChrom,numChns,numPairs);
    pValAct = zeros(numChrom,numChns,numConds);
    pValDiff = zeros(numChrom,numChns,numConds);
    pValActHBC = zeros(numChrom,numChns,numConds);
    pValDiffHBC = zeros(numChrom,numChns,numConds);
    
    listPairs = nchoosek([1;2;3],2);
    
    pCV = 0.05;
    s = zeros(numChns,numConds,numChrom);
    se = zeros(numChns,numConds,numChrom);
    
    for i2 = 1:numChns
        for i3 = 1:numConds
            for i4 = 1:numChrom
                % use default 0.05 alpha level, 2-tailed. Both are default
                
                snum = 0;
                sdem = 0;
                t1 = 0;
                for i = 1:numSbj
                    % betas are not independent at all, but still assume
                    % nonetheless.
                    % weighted average % beta is #sbj x HbX x #Channels x #conditions
                    snum = snum + betaVar(i,i4,i2,i3)/bvar(i,i4,i2,i3);
                    sdem = sdem + 1/bvar(i,i4,i2,i3);
                end
                
                s(i2,i3,i4) = snum/sdem;
                
                for i = 1:numSbj
                    % weighted SE
                    t1 = t1 + (betaVar(i,i4,i2,i3)-s(i2,i3,i4))^2/bvar(i,i4,i2,i3);
                end
                
                se(i2,i3,i4) = sqrt(t1/((numSbj-1)*sdem));
                
                tStat = s(i2,i3,i4)/se(i2,i3,i4);
                dof = numSbj-1;
                
                p = 2*tcdf(abs(tStat),dof,'upper');
                pValAct(i4,i2,i3) = p;
                %sigAct(i4,i2,i3) = p<pCV;
            end
        end
    end
    
    for i2 = 1:numChns
        for i4 = 1:numChrom
            tempPVal = zeros(1,3);
            [~,idTemp] = sort(pValAct(i4,:,i2));
            for iK = 1:m
                tempPVal(iK) = pValAct(i4,idTemp(iK),i2)*(m-iK+1);
                pValActHBC(i4,i2,idTemp(iK)) = max(tempPVal);
            end
            for i3 = 1:numConds
                sigAct(i4,i2,i3) = pValActHBC(i4,i2,i3)<pCV;
            end
        end
    end
    
    % Welch's test
    for i2 = 1:numChns
        for i3 = 1:numPairs
            for i4 = 1:numChrom
%                 % use default 0.05 alpha level, 2-tailed. Both are default
%                 sigDiff(i4,i2,i3) = ttest2(squeeze(betaVar(:,i4,i2,listPairs(i3,1))),...
%                     squeeze(betaVar(:,i4,i2,listPairs(i3,2))),'Vartype','unequal');
                num = s(i2,listPairs(i3,1),i4)-s(i2,listPairs(i3,2),i4);
                dem = sqrt(se(i2,listPairs(i3,1),i4)^2 + se(i2,listPairs(i3,2),i4)^2);
                tStat2 = num/dem;
                % Welch-Satterhwaite equation
                dof_WS = dem^4/(se(i2,listPairs(i3,1),i4)^4/(numSbj-1)+se(i2,listPairs(i3,2),i4)^4/(numSbj-1));
                
                p = 2*tcdf(abs(tStat2),dof_WS,'upper');
                pValDiff(i4,i2,i3) = p;
                sigDiff(i4,i2,i3) = p<pCV;
            end
        end
    end
    
    for i2 = 1:numChns
        for i4 = 1:numChrom
            tempPVal = zeros(1,3);
            [~,idTemp] = sort(pValDiff(i4,:,i2));
            for iK = 1:m
                tempPVal(iK) = pValDiff(i4,idTemp(iK),i2)*(m-iK+1);
                pValDiffHBC(i4,idTemp(iK),i2) = max(tempPVal);
            end
            for i3 = 1:numConds
                sigDiff(i4,i2,i3) = pValDiffHBC(i4,i2,i3)<pCV;
            end
        end
    end

end