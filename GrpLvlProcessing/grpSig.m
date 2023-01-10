% compute sig table: absolute & normalized

function grpSig(sbjList)
    
    saveDirGp = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll';

    numChn = 30;
    numCond = 3;
    numPairs = 3;
    numSubsets = 6;
        
    idxLFEF = 1:10;
    idxRFEF = 11:20;
    idxLSTG = 21;
    idxRSTG = 22;
    idxLIPS = [23 24 27 28];
    idxRIPS = [25 26 29 30];
    
    subIdx = {idxLFEF, idxRFEF, idxLSTG, idxRSTG, idxLIPS, idxRIPS};

    sigActAbs = zeros(2,numChn,numCond);
    sigDiffAbs = zeros(2,numChn,numPairs);
    
    sigActGrpAbs = zeros(2,numSubsets,numCond);
    sigDiffGrpAbs = zeros(2,numSubsets,numCond);
    
    idxDel = [25 28 33 34];

    for i = 1:length(sbjList)
        
        sbjNum = sbjList{i};
        saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
        
        sigFN = [saveDir filesep 'sig.mat'];
        load(sigFN,'sigAct','sigDiff');
        
        if strcmp(sbjNum,'08')
            
            sigAct(:,idxDel,:) = [];
            sigDiff(:,idxDel,:) = [];
            
            sigAct(:,:,[1 3 5]) = [];
            %sigDiff(:,:,[1 3 5]) = [];
            
            sigAct = sigAct(:,:,[2 1 3]);
            %sigDiff = sigDiff(:,:,[4 6 2]);
            
        end
        
        sigActAbs = sigActAbs + sigAct;
        sigDiffAbs = sigDiffAbs + sigDiff;
        
    end
    
    for i = 1:length(subIdx)
        
        sigActGrpAbs(:,i,:) = sum(sigActAbs(:,subIdx{i},:),2);
        sigDiffGrpAbs(:,i,:) = sum(sigDiffAbs(:,subIdx{i},:),2);
        
    end
    
    sigActGrpNorm = sigActGrpAbs./numChn;
    sigDiffGrpNorm = sigDiffGrpAbs./numChn;
    
    fn = [saveDirGp filesep 'sigTot.mat'];
    save(fn,'sigActGrpAbs','sigDiffGrpAbs');

end