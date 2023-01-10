% CV where SS beta coefficients from training fold is passed to test fold
% 03/24/2022: Remove rejected trials from classification
%   Remove cueOnsetIndex, ineffective method. Just directly filter allS
% Only for sbj 12 and after
%
% For figure 7 and 8. Different window length and island.
%
% LDA Ledoit only
%
% Features Interaction

function preprocessFNIRS06_FInteraction(s,numClasses,rejTrOp,rejChnOp,fefOnlyOp,rerunOp,plotOp)

sbjNum = s.name;
rawDataFN = s.fName;
movieList = s.movieList;
behData = s.resp;
startT = s.startT;
endT = s.endT;

saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
    processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];

fs = 50;
timeLen = [0.1 0.2 0.5 1 1.5 2 3 4 5].*fs;
%timeLen = 1*fs;
    
nrept = 10;
kFold = 5;
zeroT = 2*fs;
    
if rerunOp
    
    % Convert nirs to snirf file format
    % snirf1 is entire data
    snirf1 = SnirfClass(load([rawDataFN{1} '.nirs'],'-mat'));
    snirf3 = SnirfClass(load([rawDataFN{2} '.nirs'],'-mat'));
    snirf4 = SnirfClass(load([rawDataFN{3} '.nirs'],'-mat'));
    snirf1.Info()

    % Extract aux and convert to stimclass
    allS2 = find(snirf1.aux(1,1).dataTimeSeries>1);
    aInd2 = find(allS2(2:end)-allS2(1:end-1)==1);
    allS3 = find(snirf3.aux(1,1).dataTimeSeries>1);
    aInd3 = find(allS3(2:end)-allS3(1:end-1)==1);
    allS4 = find(snirf4.aux(1,1).dataTimeSeries>1);
    aInd4 = find(allS4(2:end)-allS4(1:end-1)==1);
    allS2(aInd2) = [];
    allS2 = allS2./50;
    allS3(aInd3) = [];
    allS3 = allS3./50;
    allS4(aInd4) = [];
    allS4 = allS4./50;
    
    snirf2TLen = size(snirf1.data.time,1)/(50);
    snirf3TLen = size(snirf3.data.time,1)/(50);
    allS3 = allS3+snirf2TLen;
    allS4 = allS4+snirf2TLen+snirf3TLen;

    allS = [allS2; allS3; allS4];

    if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
        cueOnsetIndex = 1:4:720;
    else
        cueOnsetIndex = 1:4:360;
    end

    allS = allS(cueOnsetIndex);

    if startT ~= 1
        % allS is in sec
        allS = allS-startT;
    end

    load([saveDir filesep movieList '.mat'],'indexMoviesTest');
    load([saveDir filesep behData '.mat'],'responsesA','responsesV','correctRespA','correctRespV');

    indexMoviesTest = updateMovieList(allS,indexMoviesTest);

    % dAll = DataClass([snirf1.data.dataTimeSeries; snirf3.data.dataTimeSeries; snirf4.data.dataTimeSeries],...
    %     0:1/50:(size(snirf1.data.dataTimeSeries,1)+size(snirf3.data.dataTimeSeries,1)+size(snirf4.data.dataTimeSeries,1)-1)/50,snirf1.data.measurementList);

    dAll = DataClass([snirf1.data.dataTimeSeries; snirf3.data.dataTimeSeries; snirf4.data.dataTimeSeries],...
        0:1/50:(size(snirf1.data.dataTimeSeries,1)+size(snirf3.data.dataTimeSeries,1)+size(snirf4.data.dataTimeSeries,1)-1)/50,snirf1.data.measurementList);


    if endT ~= -1
        dAll.dataTimeSeries = dAll.dataTimeSeries(startT*fs:endT*fs,:);
        dAll.time = dAll.time(startT*fs:endT*fs);
    else
        dAll.dataTimeSeries = dAll.dataTimeSeries(startT*fs:end,:);
        dAll.time = dAll.time(startT*fs:end);
    end

    data = dAll;
    probe = snirf1.probe;
    mlActMan = {};
    tIncMan = {};
    Aaux = [];
    rcMap = [];

    if rejChnOp
        % Consider replacing this one with Welsh periodogram method.
        %mlActAuto = hmrR_PruneChannels(data,probe,mlActMan,tIncMan,[10000  10000000],2,[0  45]);
        mlActAuto = hmrR_PruneChannels(data,probe,mlActMan,tIncMan,[7000  10000000],1.5,[0  45]);
    else
        mlActAuto = {ones(size(snirf1.data.measurementList,2),1)};
    end
    
    dod = hmrR_Intensity2OD(data);

    % For sbj 12, params in 1st line for motion correct/artifact perform
    % significantly better than params in 2nd line.
    [dod] = hmrR_MotionCorrectPCArecurse(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,15,4,0.97,5,1);
    %[dod,tInc,svs,nSV,tInc0] = hmrR_MotionCorrectPCArecurse(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,20,5,0.97,5,1);

    tIncAuto = hmrR_MotionArtifact(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,15,5);
    %tIncAuto = hmrR_MotionArtifact(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,20,4);

    % Here define stim class solely for trials rejection.
    %[stimTr,~] = hmrR_StimRejection(dod,stim,tIncAuto,tIncMan,[-2  15]);

    % Index of 36 channel
    if fefOnlyOp
        idxChns = [12:21];
        idxIPSChns = [27 28 30 31 33 34 35 36];
        idxIPSChns_3FLO = [27 28 30 33 34];
        % 12, 16, 18
        %featsLO = [1:10 13 14 15 17 19 20 21];
        %featsLO = [13 14 15 17 19 20 21];
        %featsLO = [12 13 16 18 19 20 21];
        featsLO = [13 15 16 17 19 20 21];
    else
        idxChns = [1:10 12:21 23 25 27 28 31 32 35 36 37 38];
        featsLO = [1:10 13 14 15 17 19 20 21 23 25 27 28 31 32 35 36 37 38];
    end
    
    startTFeat = 2*fs;

    if rejTrOp == 1

        leftOnsetMAll = StimClass('leftMulti');
        rightOnsetMAll = StimClass('rightMulti');
        centerOnsetMAll = StimClass('centerMulti');

        leftMultiIndexAll = indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1;
        rightMultiIndexAll = indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1;
        centerMultiIndexAll = indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==1;

        AddStims(leftOnsetMAll, allS(leftMultiIndexAll));
        AddStims(rightOnsetMAll, allS(rightMultiIndexAll));
        AddStims(centerOnsetMAll, allS(centerMultiIndexAll));

        updateStates(leftOnsetMAll);
        updateStates(rightOnsetMAll);
        updateStates(centerOnsetMAll);

        stimAll(1,1) = leftOnsetMAll;
        stimAll(1,2) = rightOnsetMAll;
        stimAll(1,3) = centerOnsetMAll;

        [stimAll,~] = hmrR_StimRejection(dod,stimAll,tIncAuto,tIncMan,[-2  15]);

        [indexMoviesTest] = updateMovieListStates(indexMoviesTest,rejTrOp,stimAll);

    else

        [indexMoviesTest] = updateMovieListStates(indexMoviesTest,rejTrOp);

    end

    % 90 trials
    trialIdx = (responsesA==correctRespA)&(responsesV==correctRespV)&(indexMoviesTest(:,5)==1)'&(indexMoviesTest(:,7)>0)';
    if numClasses == 2
        trialIdx = trialIdx & (indexMoviesTest(:,2)==1|indexMoviesTest(:,2)==2)';
    end
    if strcmp(sbjNum,'15')
        trialIdx(1) = 0;
    end
    numCorrectTrials = sum(trialIdx);
    indexMoviesMultiple = indexMoviesTest(trialIdx,:);
    allS = allS(trialIdx);
    %cueOnsetIndexMultiple = cueOnsetIndex(trialIdx);

    % all channel
    %performanceLDALedoitHbO_LOFO = zeros(nrept*kFold,size(idxFEF,2));
    %performanceLDALedoitHbR_LOFO = zeros(nrept*kFold,size(idxFEF,2));
    performanceLDALedoitHbT_All = zeros(nrept*kFold,length(timeLen));
    performanceLDALedoitHbT_3FLO = zeros(nrept*kFold,length(timeLen));
    performanceLDALedoitHbT_IPS = zeros(nrept*kFold,length(timeLen));
    performanceLDALedoitHbT_IPS_3FLO = zeros(nrept*kFold,length(timeLen));

    mlListAll = mlActAuto{1};
    numChnAll = length(mlListAll)/2;
    mlList = mlListAll([idxChns idxChns+numChnAll]);
    numChn = length(mlList)/2;
    % all chns
    mlList_3FLO = mlListAll([featsLO featsLO+numChnAll]);
    %mlList_3FLO([featsLO featsLO+numChnAll]) = [];
    if fefOnlyOp

        mlListIPS = mlListAll([idxIPSChns idxIPSChns+numChnAll]);
        mlListIPS_3FLO = mlListAll([idxIPSChns_3FLO idxIPSChns_3FLO+numChnAll]);
        numChnIPS = length(mlListIPS)/2;
        %numChnIPS_3FLO = length(mlListIPS_3FLO)/2;
        
    end


    for iRep = 1:nrept

        cvp = cvpartition(numCorrectTrials, 'KFold',kFold);

        for iFold = 1:kFold

            foldIdx = (iRep-1)*kFold+iFold;

            % snirf 2. Cell of StimClass array.
            leftOnsetMTr = StimClass('leftMulti');
            rightOnsetMTr = StimClass('rightMulti');
            centerOnsetMTr = StimClass('centerMulti');

            leftOnsetMTst = StimClass('leftMulti');
            rightOnsetMTst = StimClass('rightMulti');
            centerOnsetMTst = StimClass('centerMulti');

            indexMoviesTrainMultiple = indexMoviesMultiple(cvp.training(iFold),:);

            indexMoviesTestMultiple = indexMoviesMultiple(cvp.test(iFold),:);

            allSMultipleTr = allS(cvp.training(iFold));
            allSMultipleTst = allS(cvp.test(iFold));

            leftMultiIndexTr = indexMoviesTrainMultiple(:,2)==2 & indexMoviesTrainMultiple(:,5)==1;
            rightMultiIndexTr = indexMoviesTrainMultiple(:,2)==1 & indexMoviesTrainMultiple(:,5)==1;
            centerMultiIndexTr = indexMoviesTrainMultiple(:,2)==3 & indexMoviesTrainMultiple(:,5)==1;

            leftMultiIndexTst = indexMoviesTestMultiple(:,2)==2 & indexMoviesTestMultiple(:,5)==1;
            rightMultiIndexTst = indexMoviesTestMultiple(:,2)==1 & indexMoviesTestMultiple(:,5)==1;
            centerMultiIndexTst = indexMoviesTestMultiple(:,2)==3 & indexMoviesTestMultiple(:,5)==1;

            % Only correct trials
            AddStims(leftOnsetMTr, allSMultipleTr(leftMultiIndexTr));
            AddStims(rightOnsetMTr, allSMultipleTr(rightMultiIndexTr));
            AddStims(centerOnsetMTr, allSMultipleTr(centerMultiIndexTr));

            AddStims(leftOnsetMTst, allSMultipleTst(leftMultiIndexTst));
            AddStims(rightOnsetMTst, allSMultipleTst(rightMultiIndexTst));
            AddStims(centerOnsetMTst, allSMultipleTst(centerMultiIndexTst));

            updateStates(leftOnsetMTr);
            updateStates(rightOnsetMTr);
            updateStates(centerOnsetMTr);

            updateStates(leftOnsetMTst);
            updateStates(rightOnsetMTst);
            updateStates(centerOnsetMTst);

            % I prefer this over SetStim so I can control index
            stimTr(1,1) = leftOnsetMTr;
            stimTr(1,2) = rightOnsetMTr;
            if numClasses == 3
                stimTr(1,3) = centerOnsetMTr;
            end

            stimTst(1,1) = leftOnsetMTst;
            stimTst(1,2) = rightOnsetMTst;
            if numClasses == 3
                stimTst(1,3) = centerOnsetMTst;
            end

            [stimTr,~] = hmrR_StimRejection(dod,stimTr,tIncAuto,tIncMan,[-2  15]);

            %dodBPFilt = hmrR_BandpassFilt(dod,0.01,0.5);
            dodBPFilt = hmrR_BandpassFilt(dod,0.01,10);

            dc = hmrR_OD2Conc(dodBPFilt,probe,[1  1  1]);

            % GLM here
            % original beta var only returns coefficient for temporal basis
%             [~,~,~,dcNewTr,~,~,betaOrig,~,~,~,beta] = hmrR_GLM_ssBeta(dc,stimTr,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
%                 [-2  15],1,1,[1.0 1.0 0.0 0.0 0.0 0.0],15,1,0,0);

            [~,~,~,dcNewTr,~,~,~,~,~,~,betaSS] = hmrR_GLM_MN(dc,stimTr,probe,mlActAuto,Aaux,tIncAuto,rcMap,...
                [-2  15],1,1,[1.0 1.0 0.0 0.0 0.0 0.0],15,1,0,0);

            % format beta var here
            ssBeta = betaSS{1}(end,:,:);

            [stimTst,~] = hmrR_StimRejection(dod,stimTst,tIncAuto,tIncMan,[-2  15]);

            % actually this is not needed, can use stimTst on dcNewTr.
            dcNewTst = hmrR_ssBeta_CV(data,dc, probe, mlActAuto, tIncAuto, squeeze(ssBeta));

    %         allSMultipleTr = indexMoviesTrainMultiple(:,6);
    %         allSMultipleTst = indexMoviesTestMultiple(:,6);


            % reject trial, stimTr into allSMultipleTr and

    %         if rejTrOp == 1
    %         
    %             [allSMultipleTr,indexMoviesTrainMultiple] = ...
    %                 updateStimAndMovieList(stimTr,allSMultipleTr,indexMoviesTrainMultiple);
    % 
    %             [allSMultipleTst,indexMoviesTestMultiple] = ...
    %                 updateStimAndMovieList(stimTst,allSMultipleTst,indexMoviesTestMultiple);
    %         
    %         end

            [~, ~, trials_HbTM_Tr] ...
                = createSingleTrialHRF_MultiOnly_NoSave(dcNewTr, allSMultipleTr);

            [~, ~, trials_HbTM_Tst] ...
                = createSingleTrialHRF_MultiOnly_NoSave(dcNewTst, allSMultipleTst);
            
            trials_HbTM_Tr = offsetTrials(trials_HbTM_Tr,zeroT);
            trials_HbTM_Tst = offsetTrials(trials_HbTM_Tst,zeroT);
            
            trials_HbOM_Tr_FEF = trials_HbTM_Tr(idxChns,:,:);
            trials_HbOM_Tst_FEF = trials_HbTM_Tst(idxChns,:,:);
            
            if fefOnlyOp
                
                trials_HbOM_Tr_IPS = trials_HbTM_Tr(idxIPSChns,:,:);
                trials_HbOM_Tst_IPS = trials_HbTM_Tst(idxIPSChns,:,:);
                
                trials_HbOM_Tr_IPS_3FLO = trials_HbTM_Tr(idxIPSChns_3FLO,:,:);
                trials_HbOM_Tst_IPS_3FLO = trials_HbTM_Tst(idxIPSChns_3FLO,:,:);
                
            end

            % 3-Features Left Out
            trials_HbOM_Tr_3FLO = trials_HbTM_Tr(featsLO,:,:);
            trials_HbOM_Tst_3FLO = trials_HbTM_Tst(featsLO,:,:);
            
            for iTimeIdx = 1:length(timeLen)
                
                thisTLen = timeLen(iTimeIdx);
                
                temp = cumsum(trials_HbOM_Tr_FEF(:,startTFeat:startTFeat+thisTLen,:),2);
                cumsumTr = squeeze(temp(:,end,:));

                temp = cumsum(trials_HbOM_Tst_FEF(:,startTFeat:startTFeat+thisTLen,:),2);
                cumsumTst = squeeze(temp(:,end,:));

                performanceLDALedoitHbT_All(foldIdx,iTimeIdx) = train_RLDA_Ledoit_TrTst(cumsumTr,cumsumTst,mlList,numChn,numClasses,indexMoviesTrainMultiple,indexMoviesTestMultiple);

                temp = cumsum(trials_HbOM_Tr_3FLO(:,startTFeat:startTFeat+thisTLen,:),2);
                cumsumTr_3FLO = squeeze(temp(:,end,:));

                temp = cumsum(trials_HbOM_Tst_3FLO(:,startTFeat:startTFeat+thisTLen,:),2);
                cumsumTst_3FLO = squeeze(temp(:,end,:));

                performanceLDALedoitHbT_3FLO(foldIdx,iTimeIdx) = train_RLDA_Ledoit_TrTst(cumsumTr_3FLO,cumsumTst_3FLO,mlList_3FLO,numChn-3,numClasses,indexMoviesTrainMultiple,indexMoviesTestMultiple);

                if fefOnlyOp
                    
                    temp = cumsum(trials_HbOM_Tr_IPS(:,startTFeat:startTFeat+thisTLen,:),2);
                    cumsumTrIPS = squeeze(temp(:,end,:));

                    temp = cumsum(trials_HbOM_Tst_IPS(:,startTFeat:startTFeat+thisTLen,:),2);
                    cumsumTstIPS = squeeze(temp(:,end,:));
                    
                    performanceLDALedoitHbT_IPS(foldIdx,iTimeIdx) = train_RLDA_Ledoit_TrTst(cumsumTrIPS,cumsumTstIPS,mlListIPS,numChnIPS,numClasses,indexMoviesTrainMultiple,indexMoviesTestMultiple);
                    
                    temp = cumsum(trials_HbOM_Tr_IPS_3FLO(:,startTFeat:startTFeat+thisTLen,:),2);
                    cumsumTrIPS_3FLO = squeeze(temp(:,end,:));
                    
                    temp = cumsum(trials_HbOM_Tst_IPS_3FLO(:,startTFeat:startTFeat+thisTLen,:),2);
                    cumsumTstIPS_3FLO = squeeze(temp(:,end,:));
                    
                    performanceLDALedoitHbT_IPS_3FLO(foldIdx,iTimeIdx) = train_RLDA_Ledoit_TrTst(cumsumTrIPS_3FLO,cumsumTstIPS_3FLO,mlListIPS_3FLO,numChnIPS-3,numClasses,indexMoviesTrainMultiple,indexMoviesTestMultiple);
                    
                end
                
            end
        end


    end

    % save folds to file

    if ~exist(processedDataDir,'dir')
        mkdir(processedDataDir);
    end

    if numClasses == 2
        if fefOnlyOp
            fileName = [processedDataDir filesep 'performance_3FLO_FEF_LR_10Hz.mat'];
        else
            fileName = [processedDataDir filesep 'performance_3FLO_LR_10Hz.mat'];
        end
    else
        if fefOnlyOp
            fileName = [processedDataDir filesep 'performance_3FLO_FEF_10Hz.mat'];
        else
            fileName = [processedDataDir filesep 'performance_3FLO_10Hz.mat'];
        end
    end

    save(fileName,'performanceLDALedoitHbT_All','performanceLDALedoitHbT_3FLO',...
        'performanceLDALedoitHbT_IPS','performanceLDALedoitHbT_IPS_3FLO');
    
else
    if numClasses == 2
        if fefOnlyOp
            fileName = [processedDataDir filesep 'performance_LOFO_LR.mat'];
        else
            fileName = [processedDataDir filesep 'performance_LOFO_LR_AllChns.mat'];
        end
    else
        if fefOnlyOp
            fileName = [processedDataDir filesep 'performance_LOFO.mat'];
        else
            fileName = [processedDataDir filesep 'performance_LOFO_AllChns.mat'];
        end
    end
    
    load(fileName,'performanceLDALedoitHbT_All','performanceLDALedoitHbT_LOFO');

end
   

% Not updated
if plotOp
    
    figSaveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' ...
        num2str(sbjNum) '\Figures\Classification'];

    yLimAxis = [0 1];
    
    performanceLDALedoitHbOCI = zeros(length(timeLen),2);
    performanceLDALedoitHbRCI = zeros(length(timeLen),2);
    performanceLDALedoitHbTCI = zeros(length(timeLen),2);

    performance_LFEF_HbOCI = zeros(length(timeLen),2);
    performance_LFEF_HbRCI = zeros(length(timeLen),2);
    performance_LFEF_HbTCI = zeros(length(timeLen),2);

    performance_RFEF_HbOCI = zeros(length(timeLen),2);
    performance_RFEF_HbRCI = zeros(length(timeLen),2);
    performance_RFEF_HbTCI = zeros(length(timeLen),2);

    performance_LSTG_HbOCI = zeros(length(timeLen),2);
    performance_LSTG_HbRCI = zeros(length(timeLen),2);
    performance_LSTG_HbTCI = zeros(length(timeLen),2);

    performance_RSTG_HbOCI = zeros(length(timeLen),2);
    performance_RSTG_HbRCI = zeros(length(timeLen),2);
    performance_RSTG_HbTCI = zeros(length(timeLen),2);

    performance_LIPS_HbOCI = zeros(length(timeLen),2);
    performance_LIPS_HbRCI = zeros(length(timeLen),2);
    performance_LIPS_HbTCI = zeros(length(timeLen),2);

    performance_RIPS_HbOCI = zeros(length(timeLen),2);
    performance_RIPS_HbRCI = zeros(length(timeLen),2);
    performance_RIPS_HbTCI = zeros(length(timeLen),2);

    % performanceLDACERNNHbOCI = zeros(length(timeLen),2);
    % performanceLDACERNNHbRCI = zeros(length(timeLen),2);
    % performanceLDACERNNHbTCI = zeros(length(timeLen),2);
    % 
    % performanceLogRegHbOCI = zeros(length(timeLen),2);
    % performanceLogRegHbRCI = zeros(length(timeLen),2);
    % performanceLogRegHbTCI = zeros(length(timeLen),2);
    % 
    % performanceSVMHbOCI = zeros(length(timeLen),2);
    % performanceSVMHbRCI = zeros(length(timeLen),2);
    % performanceSVMHbTCI = zeros(length(timeLen),2);
    % 
    % performanceBaggingHbOCI = zeros(length(timeLen),2);
    % performanceBaggingHbRCI = zeros(length(timeLen),2);
    % performanceBaggingHbTCI = zeros(length(timeLen),2);
    % 
    % performanceBoostedHbOCI = zeros(length(timeLen),2);
    % performanceBoostedHbRCI = zeros(length(timeLen),2);
    % performanceBoostedHbTCI = zeros(length(timeLen),2);

    % mean
    performanceLDALedoitHbOMean = zeros(length(timeLen),1);
    performanceLDALedoitHbRMean = zeros(length(timeLen),1);
    performanceLDALedoitHbTMean = zeros(length(timeLen),1);

    performance_LFEF_HbOMean = zeros(length(timeLen),1);
    performance_LFEF_HbRMean = zeros(length(timeLen),1);
    performance_LFEF_HbTMean = zeros(length(timeLen),1);

    performance_RFEF_HbOMean = zeros(length(timeLen),1);
    performance_RFEF_HbRMean = zeros(length(timeLen),1);
    performance_RFEF_HbTMean = zeros(length(timeLen),1);

    performance_LSTG_HbOMean = zeros(length(timeLen),1);
    performance_LSTG_HbRMean = zeros(length(timeLen),1);
    performance_LSTG_HbTMean = zeros(length(timeLen),1);

    performance_RSTG_HbOMean = zeros(length(timeLen),1);
    performance_RSTG_HbRMean = zeros(length(timeLen),1);
    performance_RSTG_HbTMean = zeros(length(timeLen),1);

    performance_LIPS_HbOMean = zeros(length(timeLen),1);
    performance_LIPS_HbRMean = zeros(length(timeLen),1);
    performance_LIPS_HbTMean = zeros(length(timeLen),1);

    performance_RIPS_HbOMean = zeros(length(timeLen),1);
    performance_RIPS_HbRMean = zeros(length(timeLen),1);
    performance_RIPS_HbTMean = zeros(length(timeLen),1);
    
    conf = 0.95;

    for i = 1:length(timeLen)

        temp = performanceLDALedoitHbO_LOFO(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performanceLDALedoitHbOMean(i) = mean(temp(:));
        performanceLDALedoitHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        temp = performanceLDALedoitHbR_LOFO(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performanceLDALedoitHbRMean(i) = mean(temp(:));
        performanceLDALedoitHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        temp = performanceLDALedoitHbT_LOFO(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performanceLDALedoitHbTMean(i) = mean(temp(:));
        performanceLDALedoitHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        % LFEF
        temp = performanceLDALedoitHbO_LFEF(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_LFEF_HbOMean(i) = mean(temp(:));
        performance_LFEF_HbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        temp = performanceLDALedoitHbR_LFEF(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_LFEF_HbRMean(i) = mean(temp(:));
        performance_LFEF_HbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        temp = performanceLDALedoitHbT_LFEF(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_LFEF_HbTMean(i) = mean(temp(:));
        performance_LFEF_HbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        % RFEF
        temp = performanceLDALedoitHbO_RFEF(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_RFEF_HbOMean(i) = mean(temp(:));
        performance_RFEF_HbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        temp = performanceLDALedoitHbR_RFEF(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_RFEF_HbRMean(i) = mean(temp(:));
        performance_RFEF_HbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        temp = performanceLDALedoitHbT_RFEF(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_RFEF_HbTMean(i) = mean(temp(:));
        performance_RFEF_HbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        %LSTG
        temp = performanceLDALedoitHbO_LSTG(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_LSTG_HbOMean(i) = mean(temp(:));
        performance_LSTG_HbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        temp = performanceLDALedoitHbR_LSTG(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_LSTG_HbRMean(i) = mean(temp(:));
        performance_LSTG_HbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        temp = performanceLDALedoitHbT_LSTG(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_LSTG_HbTMean(i) = mean(temp(:));
        performance_LSTG_HbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        % RSTG
        temp = performanceLDALedoitHbO_RSTG(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_RSTG_HbOMean(i) = mean(temp(:));
        performance_RSTG_HbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        temp = performanceLDALedoitHbR_RSTG(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_RSTG_HbRMean(i) = mean(temp(:));
        performance_RSTG_HbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        temp = performanceLDALedoitHbT_RSTG(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_RSTG_HbTMean(i) = mean(temp(:));
        performance_RSTG_HbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        % LIPS
        temp = performanceLDALedoitHbO_LIPS(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_LIPS_HbOMean(i) = mean(temp(:));
        performance_LIPS_HbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        temp = performanceLDALedoitHbR_LIPS(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_LIPS_HbRMean(i) = mean(temp(:));
        performance_LIPS_HbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        temp = performanceLDALedoitHbT_LIPS(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_LIPS_HbTMean(i) = mean(temp(:));
        performance_LIPS_HbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        %RIPS
        temp = performanceLDALedoitHbO_RIPS(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_RIPS_HbOMean(i) = mean(temp(:));
        performance_RIPS_HbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        temp = performanceLDALedoitHbR_RIPS(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_RIPS_HbRMean(i) = mean(temp(:));
        performance_RIPS_HbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        temp = performanceLDALedoitHbT_RIPS(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_RIPS_HbTMean(i) = mean(temp(:));
        performance_RIPS_HbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);

    %     % CERNN
    %     temp = performanceLDACERNNHbOFold(i,:,:);
    %     se = std(temp(:))/sqrt(nrept*kFold);
    %     performanceLDACERNNHbOMean(i) = mean(temp(:));
    %     performanceLDACERNNHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    %     
    %     temp = performanceLDACERNNHbRFold(i,:,:);
    %     se = std(temp(:))/sqrt(nrept*kFold);
    %     performanceLDACERNNHbRMean(i) = mean(temp(:));
    %     performanceLDACERNNHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    %     
    %     temp = performanceLDACERNNHbTFold(i,:,:);
    %     se = std(temp(:))/sqrt(nrept*kFold);
    %     performanceLDACERNNHbTMean(i) = mean(temp(:));
    %     performanceLDACERNNHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    %     
    %     % Log
    %     temp = performanceLogRegHbOFold(i,:,:);
    %     se = std(temp(:))/sqrt(nrept*kFold);
    %     performanceLogRegHbOMean(i) = mean(temp(:));
    %     performanceLogRegHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    %     
    %     temp = performanceLogRegHbRFold(i,:,:);
    %     se = std(temp(:))/sqrt(nrept*kFold);
    %     performanceLogRegHbRMean(i) = mean(temp(:));
    %     performanceLogRegHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    %     
    %     temp = performanceLogRegHbTFold(i,:,:);
    %     se = std(temp(:))/sqrt(nrept*kFold);
    %     performanceLogRegHbTMean(i) = mean(temp(:));
    %     performanceLogRegHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    %     
    %     % SVM
    %     temp = performanceSVMHbOFold(i,:,:);
    %     se = std(temp(:))/sqrt(nrept*kFold);
    %     performanceSVMHbOMean(i) = mean(temp(:));
    %     performanceSVMHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    %     
    %     temp = performanceSVMHbRFold(i,:,:);
    %     se = std(temp(:))/sqrt(nrept*kFold);
    %     performanceSVMHbRMean(i) = mean(temp(:));
    %     performanceSVMHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    %     
    %     temp = performanceSVMHbTFold(i,:,:);
    %     se = std(temp(:))/sqrt(nrept*kFold);
    %     performanceSVMHbTMean(i) = mean(temp(:));
    %     performanceSVMHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    %     
    %     % Bagging
    %     temp = performanceBaggingHbOFold(i,:,:);
    %     se = std(temp(:))/sqrt(nrept*kFold);
    %     performanceBaggingHbOMean(i) = mean(temp(:));
    %     performanceBaggingHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    %     
    %     temp = performanceBaggingHbRFold(i,:,:);
    %     se = std(temp(:))/sqrt(nrept*kFold);
    %     performanceBaggingHbRMean(i) = mean(temp(:));
    %     performanceBaggingHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    %     
    %     temp = performanceBaggingHbTFold(i,:,:);
    %     se = std(temp(:))/sqrt(nrept*kFold);
    %     performanceBaggingHbTMean(i) = mean(temp(:));
    %     performanceBaggingHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    %     
    %     % Boosting
    %     temp = performanceBoostedHbOFold(i,:,:);
    %     se = std(temp(:))/sqrt(nrept*kFold);
    %     performanceBoostedHbOMean(i) = mean(temp(:));
    %     performanceBoostedHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    %     
    %     temp = performanceBoostedHbRFold(i,:,:);
    %     se = std(temp(:))/sqrt(nrept*kFold);
    %     performanceBoostedHbRMean(i) = mean(temp(:));
    %     performanceBoostedHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);
    %     
    %     temp = performanceBoostedHbTFold(i,:,:);
    %     se = std(temp(:))/sqrt(nrept*kFold);
    %     performanceBoostedHbTMean(i) = mean(temp(:));
    %     performanceBoostedHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);

    end

    numSubSets = 7;
    cmap = jet(numSubSets);

    figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    subplot(1,3,1);hold on;
    % for j = 1:size(performanceArrMultiHbO,1)
    %     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
    %     plot(timePt./fs-2,performanceArrMultiHbO(j,:),'Color',cmap(j,:));
    % end
    % plot(timePt./fs-2,performanceArrMultiHbO_95,'Color',cmap(1,:));
    % plot(timePt./fs-2,performanceArrMultiHbO_85,'Color',cmap(2,:));
    % plot(timePt./fs-2,performanceArrMultiHbO_75,'Color',cmap(3,:));
    % plot(timePt./fs-2,performanceArrMultiHbO_65,'Color',cmap(4,:));
    % plot(timePt./fs-2,performanceArrMultiHbO_55,'Color',cmap(5,:));
    % plot(timePt./fs-2,performanceLDALedoitHbO,'Color',cmap(6,:));
    % plot(timePt./fs-2,performanceLDACERNNHbO,'Color',cmap(7,:));
    % % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
    % % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
    % plot(timePt./fs-2,performanceLogRegHbO,'Color',cmap(8,:));
    % plot(timePt./fs-2,performanceSVMHbO,'Color',cmap(9,:));

    errorbar(timeLen./fs,performanceLDALedoitHbOMean,performanceLDALedoitHbOCI(:,2),'Color',cmap(1,:));
    errorbar(timeLen./fs,performance_LFEF_HbOMean,performance_LFEF_HbOCI(:,2),'Color',cmap(2,:));
    errorbar(timeLen./fs,performance_RFEF_HbOMean,performance_RFEF_HbOCI(:,2),'Color',cmap(3,:));
    errorbar(timeLen./fs,performance_LSTG_HbOMean,performance_LSTG_HbOCI(:,2),'Color',cmap(4,:));
    errorbar(timeLen./fs,performance_RSTG_HbOMean,performance_RSTG_HbOCI(:,2),'Color',cmap(5,:));
    errorbar(timeLen./fs,performance_LIPS_HbOMean,performance_LIPS_HbOCI(:,2),'Color',cmap(6,:));
    errorbar(timeLen./fs,performance_RIPS_HbOMean,performance_RIPS_HbOCI(:,2),'Color',cmap(7,:));
    % errorbar(timeLen./fs-2,performanceLDACERNNHbOMean,performanceLDACERNNHbOCI(:,2),'Color',cmap(2,:));
    % % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
    % % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
    % errorbar(timeLen./fs-2,performanceLogRegHbOMean,performanceLogRegHbOCI(:,2),'Color',cmap(3,:));
    % errorbar(timeLen./fs-2,performanceSVMHbOMean,performanceSVMHbOCI(:,2),'Color',cmap(4,:));
    % errorbar(timeLen./fs-2,performanceBaggingHbOMean,performanceBaggingHbOCI(:,2),'Color',cmap(5,:));
    % errorbar(timeLen./fs-2,performanceBoostedHbOMean,performanceBoostedHbOCI(:,2),'Color',cmap(6,:));

    ylim(yLimAxis);
    %xlim([0 7]);
    title(sprintf('Sbj %s: Δ[HbO] Multi',num2str(sbjNum)));
    % legend({'PCA 95%','PCA 85%','PCA 75%','PCA 65%','PCA 55%','LDAShrink',...
    %     'QDAShrink','LDAGD','LogReg','SVM'});
    % legend({'PCA 95%','PCA 85%','PCA 75%','PCA 65%','PCA 55%','LDAShrink',...
    %     'LDA CERNN','LogReg','SVM'});
    %legend({'LDAShrink','LDA CERNN','LogReg','SVM','Bagging','Boosted'});
    %legend('LDAShrink');
    xlabel('Time [s]');ylabel('Accuracy');
    hold off;

    subplot(1,3,2);hold on;
    % for j = 1:size(performanceArrMultiHbR,1)
    %     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
    %     plot(timePt./fs-2,performanceArrMultiHbR(j,:),'Color',cmap(j,:));
    % end
    % plot(timePt./fs-2,performanceArrMultiHbR_95,'Color',cmap(1,:));
    % plot(timePt./fs-2,performanceArrMultiHbR_85,'Color',cmap(2,:));
    % plot(timePt./fs-2,performanceArrMultiHbR_75,'Color',cmap(3,:));
    % plot(timePt./fs-2,performanceArrMultiHbR_65,'Color',cmap(4,:));
    % plot(timePt./fs-2,performanceArrMultiHbR_55,'Color',cmap(5,:));
    % plot(timePt./fs-2,performanceLDALedoitHbR,'Color',cmap(6,:));
    % plot(timePt./fs-2,performanceLDACERNNHbR,'Color',cmap(7,:));
    % % plot(timePt./fs-2,performanceQDALedoitHbR,'Color',cmap(7,:));
    % % plot(timePt./fs-2,performanceLDAGDHbR,'Color',cmap(8,:));
    % plot(timePt./fs-2,performanceLogRegHbR,'Color',cmap(8,:));
    % plot(timePt./fs-2,performanceSVMHbR,'Color',cmap(9,:));

    errorbar(timeLen./fs,performanceLDALedoitHbRMean,performanceLDALedoitHbRCI(:,2),'Color',cmap(1,:));
    errorbar(timeLen./fs,performance_LFEF_HbRMean,performance_LFEF_HbRCI(:,2),'Color',cmap(2,:));
    errorbar(timeLen./fs,performance_RFEF_HbRMean,performance_RFEF_HbRCI(:,2),'Color',cmap(3,:));
    errorbar(timeLen./fs,performance_LSTG_HbRMean,performance_LSTG_HbRCI(:,2),'Color',cmap(4,:));
    errorbar(timeLen./fs,performance_RSTG_HbRMean,performance_RSTG_HbRCI(:,2),'Color',cmap(5,:));
    errorbar(timeLen./fs,performance_LIPS_HbRMean,performance_LIPS_HbRCI(:,2),'Color',cmap(6,:));
    errorbar(timeLen./fs,performance_RIPS_HbRMean,performance_RIPS_HbRCI(:,2),'Color',cmap(7,:));
    % errorbar(timeLen./fs-2,performanceLDACERNNHbRMean,performanceLDACERNNHbRCI(:,2),'Color',cmap(2,:));
    % % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
    % % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
    % errorbar(timeLen./fs-2,performanceLogRegHbRMean,performanceLogRegHbRCI(:,2),'Color',cmap(3,:));
    % errorbar(timeLen./fs-2,performanceSVMHbRMean,performanceSVMHbRCI(:,2),'Color',cmap(4,:));
    % errorbar(timeLen./fs-2,performanceBaggingHbRMean,performanceBaggingHbRCI(:,2),'Color',cmap(5,:));
    % errorbar(timeLen./fs-2,performanceBoostedHbRMean,performanceBoostedHbRCI(:,2),'Color',cmap(6,:));

    ylim(yLimAxis);
    %xlim([0 7]);
    title(sprintf('Sbj %s: Δ[HbR] Multi',num2str(sbjNum)));
    %legend('S1D1','S1D2','S1D3','S1D4');
    legend({'All','LFEF','RFEF','LSTG','RSTG','LIPS','RIPS'},'Location','southwest');
    xlabel('Time [s]');ylabel('Accuracy');
    hold off;

    subplot(1,3,3);hold on;
    % for j = 1:size(performanceArrMultiHbT,1)
    %     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
    %     plot(timePt./fs-2,performanceArrMultiHbT(j,:),'Color',cmap(j,:));
    % end
    % plot(timePt./fs-2,performanceArrMultiHbT_95,'Color',cmap(1,:));
    % plot(timePt./fs-2,performanceArrMultiHbT_85,'Color',cmap(2,:));
    % plot(timePt./fs-2,performanceArrMultiHbT_75,'Color',cmap(3,:));
    % plot(timePt./fs-2,performanceArrMultiHbT_65,'Color',cmap(4,:));
    % plot(timePt./fs-2,performanceArrMultiHbT_55,'Color',cmap(5,:));
    % plot(timePt./fs-2,performanceLDALedoitHbT,'Color',cmap(6,:));
    % plot(timePt./fs-2,performanceLDACERNNHbT,'Color',cmap(7,:));
    % % plot(timePt./fs-2,performanceQDALedoitHbT,'Color',cmap(7,:));
    % % plot(timePt./fs-2,performanceLDAGDHbT,'Color',cmap(8,:));
    % plot(timePt./fs-2,performanceLogRegHbT,'Color',cmap(8,:));
    % plot(timePt./fs-2,performanceSVMHbT,'Color',cmap(9,:));

    errorbar(timeLen./fs,performanceLDALedoitHbTMean,performanceLDALedoitHbTCI(:,2),'Color',cmap(1,:));
    errorbar(timeLen./fs,performance_LFEF_HbTMean,performance_LFEF_HbTCI(:,2),'Color',cmap(2,:));
    errorbar(timeLen./fs,performance_RFEF_HbTMean,performance_RFEF_HbTCI(:,2),'Color',cmap(3,:));
    errorbar(timeLen./fs,performance_LSTG_HbTMean,performance_LSTG_HbTCI(:,2),'Color',cmap(4,:));
    errorbar(timeLen./fs,performance_RSTG_HbTMean,performance_RSTG_HbTCI(:,2),'Color',cmap(5,:));
    errorbar(timeLen./fs,performance_LIPS_HbTMean,performance_LIPS_HbTCI(:,2),'Color',cmap(6,:));
    errorbar(timeLen./fs,performance_RIPS_HbTMean,performance_RIPS_HbTCI(:,2),'Color',cmap(7,:));
    % errorbar(timeLen./fs-2,performanceLDACERNNHbTMean,performanceLDACERNNHbTCI(:,2),'Color',cmap(2,:));
    % % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
    % % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
    % errorbar(timeLen./fs-2,performanceLogRegHbTMean,performanceLogRegHbTCI(:,2),'Color',cmap(3,:));
    % errorbar(timeLen./fs-2,performanceSVMHbTMean,performanceSVMHbTCI(:,2),'Color',cmap(4,:));
    % errorbar(timeLen./fs-2,performanceBaggingHbTMean,performanceBaggingHbTCI(:,2),'Color',cmap(5,:));
    % errorbar(timeLen./fs-2,performanceBoostedHbTMean,performanceBoostedHbTCI(:,2),'Color',cmap(6,:));

    ylim(yLimAxis);
    %xlim([0 7]);
    title(sprintf('Sbj %s: Δ[HbT] Multi',num2str(sbjNum)));
    %legend('S1D1','S1D2','S1D3','S1D4');
    %legend('LDAShrink');
    xlabel('Time [s]');ylabel('Accuracy');
    hold off;

    % subplot(1,4,4);hold on;
    % for j = 1:size(performanceArrMultiHbT,1)
    %     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
    %     plot(timePt./fs-2,zeros(1,length(timePt./fs-2)),'Color',cmap(j,:));
    % end
    % %legend(chnName(1:30));
    % annotation('textbox', [0.8, 0.05, 0.1, 0.1], 'String', ['Score: ' num2str(behScore) '%'],'FontSize',8);
    % annotation('textbox', [0.71, 0.2, 0.1, 0.1],'String','calcCumSum\_MultiOnly\_DiffStartT\_AllChn\_AllBasis.m','FontSize',8);
    % title(sprintf('All Channels. Diff Start T: %s-Class',num2str(numClasses)));
    % hold off;

    if numClasses == 2
        fn = sprintf('PerformanceCumsumVsTime_DiffTLen_StartMovie_LR_AllChns');
    else
        fn = sprintf('PerformanceCumsumVsTime_DiffTLen_StartMovie_AllChns');
    end

    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');

end

end