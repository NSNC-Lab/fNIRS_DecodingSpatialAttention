% CV where SS beta coefficients from training fold is passed to test fold
% 03/24/2022: Remove rejected trials from classification
%   Remove cueOnsetIndex, ineffective method. Just directly filter allS
% Only for sbj 12 and after
%
% For figure 7 and 8. Different window length and island.
%
% LDA Ledoit only

% Compare posterior vs anterior FEF
% Plot coding block not updated. Too lazy.

function preprocessFNIRS06_CV_GLM_ssBeta_PostVsAnt_MultiOnly(s,numClasses,rejTrOp,rejChnOp,rerunOp,plotOp)

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
    
nrept = 10;
kFold = 5;
    
if rerunOp
    
    % Convert nirs to snirf file format
    % snirf1 is entire data
    snirf1 = SnirfClass(load([rawDataFN{1} '.nirs'],'-mat'));
    snirf1.Info()

    % Extract aux and convert to stimclass
    allS2 = find(snirf1.aux(1,1).dataTimeSeries>1);
    aInd2 = find(allS2(2:end)-allS2(1:end-1)==1);

    allS2(aInd2) = [];
    allS2 = allS2./50;

    if strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
        cueOnsetIndex = 1:4:720;
    else
        cueOnsetIndex = 1:4:360;
    end

    allS = allS2(cueOnsetIndex);

    if startT ~= 1
        % allS is in sec
        allS = allS-startT;
    end

    load([saveDir filesep movieList '.mat'],'indexMoviesTest');
    load([saveDir filesep behData '.mat'],'responsesA','responsesV','correctRespA','correctRespV');

    indexMoviesTest = updateMovieList(allS,indexMoviesTest);

    % dAll = DataClass([snirf1.data.dataTimeSeries; snirf3.data.dataTimeSeries; snirf4.data.dataTimeSeries],...
    %     0:1/50:(size(snirf1.data.dataTimeSeries,1)+size(snirf3.data.dataTimeSeries,1)+size(snirf4.data.dataTimeSeries,1)-1)/50,snirf1.data.measurementList);

    dAll = DataClass(snirf1.data.dataTimeSeries,...
        0:1/50:(size(snirf1.data.dataTimeSeries,1)-1)/50,snirf1.data.measurementList);

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
    % try both splitting options
    % Option #1
    idxPFEF_1 = [1 2 4 5 6 16 17 19 20 21];
    idxAFEF_1 = [3 8 9 10 11 12 13 14 15 18];
    
    % Option #2
    idxPFEF_2 = [1 2 3 5 6 16 17 18 20 21];
    idxAFEF_2 = [4 8 9 10 11 12 13 14 15 19];
    
%     idxLFEF = 1:10;
%     idxRFEF = 12:21;
%     idxLSTG = 23;
%     idxRSTG = 25;
%     idxLIPS = [27 28 33 34];
%     idxRIPS = [30 31 35 36];

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
    % Just load existing data. No need to re-run. Too time-consuming.
%     performanceLDALedoitHbO_All = zeros(nrept*kFold,length(timeLen));
%     performanceLDALedoitHbR_All = zeros(nrept*kFold,length(timeLen));
%     performanceLDALedoitHbT_All = zeros(nrept*kFold,length(timeLen));

    % Left FEF
%     performanceLDALedoitHbO_LFEF = zeros(nrept*kFold,length(timeLen));
%     performanceLDALedoitHbR_LFEF = zeros(nrept*kFold,length(timeLen));
%     performanceLDALedoitHbT_LFEF = zeros(nrept*kFold,length(timeLen));

    % Right FEF
%     performanceLDALedoitHbO_PFEF_1 = zeros(nrept*kFold,length(timeLen));
%     performanceLDALedoitHbR_PFEF_1 = zeros(nrept*kFold,length(timeLen));
    performanceLDALedoitHbT_PFEF_1 = zeros(nrept*kFold,length(timeLen));

    % Left STG
%     performanceLDALedoitHbO_AFEF_1 = zeros(nrept*kFold,length(timeLen));
%     performanceLDALedoitHbR_AFEF_1 = zeros(nrept*kFold,length(timeLen));
    performanceLDALedoitHbT_AFEF_1 = zeros(nrept*kFold,length(timeLen));

    % Right STG
% %     performanceLDALedoitHbO_PFEF_2 = zeros(nrept*kFold,length(timeLen));
% %     performanceLDALedoitHbR_PFEF_2 = zeros(nrept*kFold,length(timeLen));
    performanceLDALedoitHbT_PFEF_2 = zeros(nrept*kFold,length(timeLen));

    % Left IPS
%     performanceLDALedoitHbO_AFEF_2 = zeros(nrept*kFold,length(timeLen));
%     performanceLDALedoitHbR_AFEF_2 = zeros(nrept*kFold,length(timeLen));
    performanceLDALedoitHbT_AFEF_2 = zeros(nrept*kFold,length(timeLen));

    mlList = mlActAuto{1};
    numChn = length(mlList)/2;
    mlList_PFEF_1 = mlList([idxPFEF_1 idxPFEF_1+numChn]);
    mlList_AFEF_1 = mlList([idxAFEF_1 idxAFEF_1+numChn]);
    mlList_PFEF_2 = mlList([idxPFEF_2 idxPFEF_2+numChn]);
    mlList_AFEF_2 = mlList([idxAFEF_2 idxAFEF_2+numChn]);
%     mlList_LIPS = mlList([idxLIPS idxLIPS+numChn]);
%     mlList_RIPS = mlList([idxRIPS idxRIPS+numChn]);

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

            dodBPFilt = hmrR_BandpassFilt(dod,0.01,0.5);

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

            [trials_HbOM_Tr, trials_HbRM_Tr, trials_HbTM_Tr] ...
                = createSingleTrialHRF_MultiOnly_NoSave(dcNewTr, allSMultipleTr);

            [trials_HbOM_Tst, trials_HbRM_Tst, trials_HbTM_Tst] ...
                = createSingleTrialHRF_MultiOnly_NoSave(dcNewTst, allSMultipleTst);

            % Here, split into subsets of channels.
%             trials_HbOM_Tr_LFEF = trials_HbOM_Tr(idxLFEF,:,:);
%             trials_HbOM_Tr_RFEF = trials_HbOM_Tr(idxRFEF,:,:);
%             trials_HbOM_Tr_LSTG = trials_HbOM_Tr(idxLSTG,:,:);
%             trials_HbOM_Tr_RSTG = trials_HbOM_Tr(idxRSTG,:,:);
%             trials_HbOM_Tr_LIPS = trials_HbOM_Tr(idxLIPS,:,:);
%             trials_HbOM_Tr_RIPS = trials_HbOM_Tr(idxRIPS,:,:);
% 
%             trials_HbRM_Tr_LFEF = trials_HbRM_Tr(idxLFEF,:,:);
%             trials_HbRM_Tr_RFEF = trials_HbRM_Tr(idxRFEF,:,:);
%             trials_HbRM_Tr_LSTG = trials_HbRM_Tr(idxLSTG,:,:);
%             trials_HbRM_Tr_RSTG = trials_HbRM_Tr(idxRSTG,:,:);
%             trials_HbRM_Tr_LIPS = trials_HbRM_Tr(idxLIPS,:,:);
%             trials_HbRM_Tr_RIPS = trials_HbRM_Tr(idxRIPS,:,:);

            trials_HbTM_Tr_PFEF_1 = trials_HbTM_Tr(idxPFEF_1,:,:);
            trials_HbTM_Tr_AFEF_1 = trials_HbTM_Tr(idxAFEF_1,:,:);
            trials_HbTM_Tr_PFEF_2 = trials_HbTM_Tr(idxPFEF_2,:,:);
            trials_HbTM_Tr_AFEF_2 = trials_HbTM_Tr(idxAFEF_2,:,:);

            % Test
%             trials_HbOM_Tst_LFEF = trials_HbOM_Tst(idxLFEF,:,:);
%             trials_HbOM_Tst_RFEF = trials_HbOM_Tst(idxRFEF,:,:);
%             trials_HbOM_Tst_LSTG = trials_HbOM_Tst(idxLSTG,:,:);
%             trials_HbOM_Tst_RSTG = trials_HbOM_Tst(idxRSTG,:,:);
%             trials_HbOM_Tst_LIPS = trials_HbOM_Tst(idxLIPS,:,:);
%             trials_HbOM_Tst_RIPS = trials_HbOM_Tst(idxRIPS,:,:);
% 
%             trials_HbRM_Tst_LFEF = trials_HbRM_Tst(idxLFEF,:,:);
%             trials_HbRM_Tst_RFEF = trials_HbRM_Tst(idxRFEF,:,:);
%             trials_HbRM_Tst_LSTG = trials_HbRM_Tst(idxLSTG,:,:);
%             trials_HbRM_Tst_RSTG = trials_HbRM_Tst(idxRSTG,:,:);
%             trials_HbRM_Tst_LIPS = trials_HbRM_Tst(idxLIPS,:,:);
%             trials_HbRM_Tst_RIPS = trials_HbRM_Tst(idxRIPS,:,:);

            trials_HbTM_Tst_PFEF_1 = trials_HbTM_Tst(idxPFEF_1,:,:);
            trials_HbTM_Tst_AFEF_1 = trials_HbTM_Tst(idxAFEF_1,:,:);
            trials_HbTM_Tst_PFEF_2 = trials_HbTM_Tst(idxPFEF_2,:,:);
            trials_HbTM_Tst_AFEF_2 = trials_HbTM_Tst(idxAFEF_2,:,:);

            % train classifier
%             performanceLDALedoitHbO_All(foldIdx,:) = ...
%                 calcCumSum_DiffTLen_AllChn_CompClassifiers_Final_NoSave(sbjNum,trials_HbOM_Tr,trials_HbOM_Tst,...
%                 indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlActAuto);
% 
%             performanceLDALedoitHbR_All(foldIdx,:) = ...
%                 calcCumSum_DiffTLen_AllChn_CompClassifiers_Final_NoSave(sbjNum,trials_HbRM_Tr,trials_HbRM_Tst,...
%                 indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlActAuto);
% 
%             performanceLDALedoitHbT_All(foldIdx,:) = ...
%                 calcCumSum_DiffTLen_AllChn_CompClassifiers_Final_NoSave(sbjNum,trials_HbTM_Tr,trials_HbTM_Tst,...
%                 indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlActAuto);
% 
%             % Left FEF
%             performanceLDALedoitHbO_LFEF(foldIdx,:) = ...
%                 calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbOM_Tr_LFEF,trials_HbOM_Tst_LFEF,...
%                 indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_LFEF);
% 
%             performanceLDALedoitHbR_LFEF(foldIdx,:) = ...
%                 calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbRM_Tr_LFEF,trials_HbRM_Tst_LFEF,...
%                 indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_LFEF);
% 
%             performanceLDALedoitHbT_LFEF(foldIdx,:) = ...
%                 calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbTM_Tr_LFEF,trials_HbTM_Tst_LFEF,...
%                 indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_LFEF);

            % Posterior FEF #1
%             performanceLDALedoitHbO_PFEF_1(foldIdx,:) = ...
%                 calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbOM_Tr_RFEF,trials_HbOM_Tst_RFEF,...
%                 indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_RFEF);
% 
%             performanceLDALedoitHbR_PFEF_1(foldIdx,:) = ...
%                 calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbRM_Tr_RFEF,trials_HbRM_Tst_RFEF,...
%                 indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_RFEF);

            performanceLDALedoitHbT_PFEF_1(foldIdx,:) = ...
                calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbTM_Tr_PFEF_1,trials_HbTM_Tst_PFEF_1,...
                indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_PFEF_1);

            % Anterior FEF #1
%             performanceLDALedoitHbO_AFEF_1(foldIdx,:) = ...
%                 calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbOM_Tr_LSTG,trials_HbOM_Tst_LSTG,...
%                 indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_LSTG);
% 
%             performanceLDALedoitHbR_AFEF_1(foldIdx,:) = ...
%                 calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbRM_Tr_LSTG,trials_HbRM_Tst_LSTG,...
%                 indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_LSTG);

            performanceLDALedoitHbT_AFEF_1(foldIdx,:) = ...
                calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbTM_Tr_AFEF_1,trials_HbTM_Tst_AFEF_1,...
                indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_AFEF_1);

            % Posterior FEF #2
%             performanceLDALedoitHbO_RSTG(foldIdx,:) = ...
%                 calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbOM_Tr_RSTG,trials_HbOM_Tst_RSTG,...
%                 indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_RSTG);
% 
%             performanceLDALedoitHbR_RSTG(foldIdx,:) = ...
%                 calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbRM_Tr_RSTG,trials_HbRM_Tst_RSTG,...
%                 indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_RSTG);

            performanceLDALedoitHbT_PFEF_2(foldIdx,:) = ...
                calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbTM_Tr_PFEF_2,trials_HbTM_Tst_PFEF_2,...
                indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_PFEF_2);

            % Anterior FEF #2
%             performanceLDALedoitHbO_LIPS(foldIdx,:) = ...
%                 calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbOM_Tr_LIPS,trials_HbOM_Tst_LIPS,...
%                 indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_LIPS);
% 
%             performanceLDALedoitHbR_LIPS(foldIdx,:) = ...
%                 calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbRM_Tr_LIPS,trials_HbRM_Tst_LIPS,...
%                 indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_LIPS);

            performanceLDALedoitHbT_AFEF_2(foldIdx,:) = ...
                calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbTM_Tr_AFEF_2,trials_HbTM_Tst_AFEF_2,...
                indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_AFEF_2);

        end


    end

    % save folds to file

    if ~exist(processedDataDir,'dir')
        mkdir(processedDataDir);
    end

    if numClasses == 2
        if rejTrOp
            fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_LR_PostVsAnt_RejTr.mat'];
        else
            fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_LR_PostVsAnt_.mat'];
        end
    else
        if rejTrOp
            fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_RejTr_PostVsAnt_.mat'];
        else
            fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_PostVsAnt_.mat'];
        end
    end

    save(fileName,'performanceLDALedoitHbT_PFEF_1','performanceLDALedoitHbT_AFEF_1',...
        'performanceLDALedoitHbT_PFEF_2','performanceLDALedoitHbT_AFEF_2');
    
else
    if numClasses == 2
        if rejTrOp
            fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_LR_PostVsAnt_RejTr.mat'];
        else
            fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_LR.mat'];
        end
    else
        if rejTrOp
            fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta_RejTr_PostVsAnt.mat'];
        else
            fileName = [processedDataDir filesep 'performance_GLM_CV_SSBeta.mat'];
        end
    end
    
    load(fileName,'performanceLDALedoitHbT_PFEF_1','performanceLDALedoitHbT_AFEF_1',...
        'performanceLDALedoitHbT_PFEF_2','performanceLDALedoitHbT_AFEF_2');

end
    
if plotOp
    
    figSaveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' ...
        num2str(sbjNum) '\Figures\Classification'];

    yLimAxis = [0 1];

%     performance_LFEF_HbOCI = zeros(length(timeLen),2);
%     performance_LFEF_HbRCI = zeros(length(timeLen),2);
    performance_PFEF_1_HbTCI = zeros(length(timeLen),2);
    performance_AFEF_1_HbTCI = zeros(length(timeLen),2);
    performance_PFEF_2_HbTCI = zeros(length(timeLen),2);
    performance_AFEF_2_HbTCI = zeros(length(timeLen),2);

    % mean
    performance_PFEF_1_HbTMean = zeros(length(timeLen),1);
    performance_AFEF_1_HbTMean = zeros(length(timeLen),1);
    performance_PFEF_2_HbTMean = zeros(length(timeLen),1);
    performance_AFEF_2_HbTMean = zeros(length(timeLen),1);

    conf = 0.95;

    for i = 1:length(timeLen)

        % LFEF
        temp = performanceLDALedoitHbT_PFEF_1(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_PFEF_1_HbTMean(i) = mean(temp(:));
        performance_PFEF_1_HbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
        
        temp = performanceLDALedoitHbT_AFEF_1(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_AFEF_1_HbTMean(i) = mean(temp(:));
        performance_AFEF_1_HbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
        
        temp = performanceLDALedoitHbT_PFEF_2(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_PFEF_2_HbTMean(i) = mean(temp(:));
        performance_PFEF_2_HbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);
        
        temp = performanceLDALedoitHbT_AFEF_2(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performance_AFEF_2_HbTMean(i) = mean(temp(:));
        performance_AFEF_2_HbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);

    end

    numSubSets = 4;
    cmap = jet(numSubSets);

    figure();hold on;
    errorbar(timeLen./fs,performance_PFEF_1_HbTMean,performance_PFEF_1_HbTCI(:,2),'Color',cmap(1,:));
    errorbar(timeLen./fs,performance_AFEF_1_HbTMean,performance_AFEF_1_HbTCI(:,2),'Color',cmap(2,:));
    errorbar(timeLen./fs,performance_PFEF_2_HbTMean,performance_PFEF_2_HbTCI(:,2),'Color',cmap(3,:));
    errorbar(timeLen./fs,performance_AFEF_2_HbTMean,performance_AFEF_2_HbTCI(:,2),'Color',cmap(4,:));
    
    ylim(yLimAxis);
    %xlim([0 7]);
    title(sprintf('Sbj %s: Δ[HbT] Multi',num2str(sbjNum)));
    % legend({'PCA 95%','PCA 85%','PCA 75%','PCA 65%','PCA 55%','LDAShrink',...
    %     'QDAShrink','LDAGD','LogReg','SVM'});
    % legend({'PCA 95%','PCA 85%','PCA 75%','PCA 65%','PCA 55%','LDAShrink',...
    %     'LDA CERNN','LogReg','SVM'});
    %legend({'LDAShrink','LDA CERNN','LogReg','SVM','Bagging','Boosted'});
    legend({'PFEF 1','AFEF 1','PFEF 2','AFEF 2'});
    xlabel('Time [s]');ylabel('Accuracy');
    hold off;

    if numClasses == 2
        fn = sprintf('PerformanceCumsumVsTime_PostVsAnt_StartMovie_LR_AllChns');
    else
        fn = sprintf('PerformanceCumsumVsTime_PostVsAnt_StartMovie_AllChns');
    end

    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');

end

end