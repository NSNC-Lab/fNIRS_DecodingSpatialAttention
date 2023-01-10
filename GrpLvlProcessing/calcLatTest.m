% lateralization test using welch's t-test for group.
% beta is (#coefficients x HbX x #Channels x #conditions)
% bvar is #coefficients x # Channels x #HbX
% bvar is square of standard error, not SD

function calcLatTest(sbjList)

    % Good news, same for both old and new probe
    leftChn = [1:6];
    rightChn = [15:20];
    
    numChn = size(leftChn,2);
    numCond = 3;
    
    betaIdxRng = 4:8;
    numCoeffUsed = size(betaIdxRng,2);
    
%     sigActAbsLeftHbO = zeros(length(sbjList),numChn,numCond);
%     sigActAbsLeftHbR = zeros(length(sbjList),numChn,numCond);
%     sigActAbsRightHbO = zeros(length(sbjList),numChn,numCond);
%     sigActAbsRightHbR = zeros(length(sbjList),numChn,numCond);
    
    betaTotLeftHbO = zeros(length(sbjList),numChn,numCond);
    betaTotLeftHbR = zeros(length(sbjList),numChn,numCond);
    betaTotRightHbO = zeros(length(sbjList),numChn,numCond);
    betaTotRightHbR = zeros(length(sbjList),numChn,numCond);
    
    bvarTotLeftHbO = zeros(length(sbjList),numChn,numCond);
    bvarTotLeftHbR = zeros(length(sbjList),numChn,numCond);
    bvarTotRightHbO = zeros(length(sbjList),numChn,numCond);
    bvarTotRightHbR = zeros(length(sbjList),numChn,numCond);
    %sigDiffAbs = zeros(2,numChn,numPairs);
    
%     sigActGrpAbs = zeros(2,numSubsets,numCond);
%     sigDiffGrpAbs = zeros(2,numSubsets,numCond);
%     
%     idxDel = [25 28 33 34];
    
    numBeta = numChn*length(sbjList);

    for i = 1:length(sbjList)
        
        sbjNum = sbjList{i};
        saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
        load([saveDir filesep 'intermediateOutputsCombined_Basis1.mat'],'dcAvg','dcAvgStd','dcNew','beta','R','bvar','hmrstats');
        
%         sigFN = [saveDir filesep 'sig.mat'];
        % sigAct is HbX x #Channels x #conditions
        % sigDiff is HbX x #Channels x #pairs
        
        beta = beta{1};
        
        condIdx = [1 2 3];
        
        numCoef = size(beta,1);
                
        if strcmp(sbjNum,'08')
        
            % left center right
            condIdx = [1 2 3];
            
        
        end
        
        %(#coefficients x HbX x #Channels x #conditions)
        betaAvg = squeeze(mean(beta(betaIdxRng,:,:,:),1));
        betaTotLeftHbO(i,:,:) = betaAvg(1,leftChn,condIdx);
        betaTotLeftHbR(i,:,:) = betaAvg(2,leftChn,condIdx);
        betaTotRightHbO(i,:,:) = betaAvg(1,rightChn,condIdx);
        betaTotRightHbR(i,:,:) = betaAvg(2,rightChn,condIdx);
        
        for i3 = 1:3
            %#coefficients x # Channels x #HbX
            bvarTotLeftHbO(i,:,i3) = (1/(numCoeffUsed^2))*sum(bvar((condIdx(i3)-1)*numCoef+betaIdxRng,leftChn,1),1);
            bvarTotLeftHbR(i,:,i3) = (1/(numCoeffUsed^2))*sum(bvar((condIdx(i3)-1)*numCoef+betaIdxRng,leftChn,2),1);
            bvarTotRightHbO(i,:,i3) = (1/(numCoeffUsed^2))*sum(bvar((condIdx(i3)-1)*numCoef+betaIdxRng,rightChn,1),1);
            bvarTotRightHbR(i,:,i3) = (1/(numCoeffUsed^2))*sum(bvar((condIdx(i3)-1)*numCoef+betaIdxRng,rightChn,2),1);
        end
        
    end
    
    betaMeanLeftHbO = squeeze(mean(betaTotLeftHbO,[1 2]));
    betaMeanLeftHbR = squeeze(mean(betaTotLeftHbR,[1 2]));
    betaMeanRightHbO = squeeze(mean(betaTotRightHbO,[1 2]));
    betaMeanRightHbR = squeeze(mean(betaTotRightHbR,[1 2]));
    
    betaVarLeftHbO = (1/(numBeta^2))*squeeze(sum(bvarTotLeftHbO,[1 2]));
    betaVarLeftHbR = (1/(numBeta^2))*squeeze(sum(bvarTotLeftHbR,[1 2]));
    betaVarRightHbO = (1/(numBeta^2))*squeeze(sum(bvarTotRightHbO,[1 2]));
    betaVarRightHbR = (1/(numBeta^2))*squeeze(sum(bvarTotRightHbR,[1 2]));
    
    tvalHbO = (betaMeanLeftHbO - betaMeanRightHbO)./sqrt(betaVarLeftHbO ...
        + betaVarRightHbO);
    
    tvalHbR = (betaMeanLeftHbR - betaMeanRightHbR)./sqrt(betaVarLeftHbR ...
        + betaVarRightHbR);
    
    dofHbO = (betaVarLeftHbO + betaVarRightHbO).^2./(betaVarLeftHbO.^2/numBeta...
        + betaVarRightHbO.^2/numBeta);
    
    dofHbR = (betaVarLeftHbR + betaVarRightHbR).^2./(betaVarLeftHbR.^2/numBeta...
        + betaVarRightHbR.^2/numBeta);
    
    pvalHbO = tcdf(abs(tvalHbO),dofHbO,'upper').*2;
    pvalHbR = tcdf(abs(tvalHbR),dofHbR,'upper').*2;

end
