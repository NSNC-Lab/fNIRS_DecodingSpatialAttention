% option to re-run training and classifying or pull from saved data
% sbjList is cell array of experiment "number". Ex: sbjList = {'01','02'};
% sbjList = {'08','10','12','13','14','15'};tInd = 8;featureName = 'Cumulative Sum';numRept = 5;numClasses = 2;saveOp = 0;rerunOp = 1;
% Group channels by "islands" of probe.

function grpPlotPerfVsTime_IndChn(goodSbj,numRept,numClasses,rerunOp,saveOp)

figSaveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll\Figures\Classification';
figDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll';

    if goodSbj == 1
        sbjList = {'08','12','13','16','19','21'};
        
        if numClasses == 2
            fn = sprintf('GrpPlot_CVPerfByGrpROINames_Good_2Class');
            dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
                filesep 'sbjListPerf_IndChn_Good_2Class.mat'];
        else
            fn = sprintf('GrpPlot_CVPerfByGrpROINames_Good_3Class');
            dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
                filesep 'sbjListPerf_IndChn_Good_3Class.mat'];
        end
    elseif goodSbj == 0
        sbjList = {'14','15'};

        if numClasses == 2
            dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
                filesep 'sbjListPerf_IndChn_Bad_2Class.mat'];
            fn = sprintf('GrpPlot_CVPerfByGrpROINames_Bad_2Class');
        else
            dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
                filesep 'sbjListPerf_IndChn_Bad_3Class.mat'];
            fn = sprintf('GrpPlot_CVPerfByGrpROINames_Bad_3Class');
        end
    else
        sbjList = {'08','12','14','15','16','17'};
        dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
            filesep 'sbjListPerf_IndChn_All.mat'];
    end

fs = 50;
timePt = (0:0.25*fs:12*fs)+2*fs;
zeroT = 2*fs;
timeLngth = 1*fs;
lenT = 1*fs;
numChn = 30;
timeStrtInt = 2;

P = 0:0.01:1;
if numClasses == 2
    N = 2;
else
    N = 3;
end

itf = log2(N) + P.*log2(P) + (1-P).*log2((1-P)/(N-1));

if numClasses == 3
    itf(1:33) = -itf(1:33);
end

% old list
% chnName = {'Middle Posterior FEF Left',...
%     'Posterior to sPCS/tgPSC Left',...
%     'Middle FEF Left',...
%     'sPCS/tgPCS Left',...
%     'Posterior FEF Left',...
%     'Inferior to FEF Left',...
%     'Anterior FEF Left',...
%     'iPCS/tgPCS Left',...
%     'Anterior to FEF Left',...
%     'iPCS Left',...
%     'Anterior FEF Right',...
%     'iPCS/tgPCS Right',...
%     'Anterior to FEF Right',...
%     'iPCS Right',...
%     'Middle Posterior FEF Right',...
%     'Posterior to sPCS/tgPSC Right',...
%     'Middle FEF Right',...
%     'sPCS/tgPCS Right',...
%     'Posterior to FEF Right',...
%     'Inferor to FEF Right',...
%     'Post Posterior STG/PT Left',...
%     'Posterior STG/PT Right',...
%     'IPS3/IPS2/SPL1 Left',...
%     'IPS3/antIPS/IPS4 Left',...
%     'IPS3/latIPS/antIPS Left',...
%     'IPS3/IPS2/SPL1 Right',...
%     'IPS3/antIPS/IPS4 Right',...
%     'IPS3/latIPS/antIPS Right',...
%     'Superior to IPS3/IPS2/SPL1 Left',...
%     'IPS4 Left',...
%     'Superior to IPS3/IPS2/SPL1 Right',...
%     'IPS4 Right',...
%     'Ant Posterior STG/PT Left',...
%     'Ant Posterior STG/PT Right'};

% new list
chnName = {'S1D1 Middle Posterior FEF Left',...
    'S1D2 Posterior to sPCS/tgPSC Left',...
    'S1D3 Middle FEF Left',...
    'S1D4 sPCS/tgPCS Left',...
    'S2D1 Posterior FEF Left',...
    'S2D2 Inferior to FEF Left',...
    'S3D3 Anterior FEF Left',...
    'S3D4 iPCS/tgPCS Left',...
    'S3D5 Anterior to FEF Left',...
    'S3D6 iPCS Left',...
    'S4D9 Anterior FEF Right',...
    'S4D10 iPCS/tgPCS Right',...
    'S4D11 Anterior to FEF Right',...
    'S4D12 iPCS Right',...
    'S5D7 Middle Posterior FEF Right',...
    'S5D8 Posterior to sPCS/tgPSC Right',...
    'S5D9 Middle FEF Right',...
    'S5D10 sPCS/tgPCS Right',...
    'S6D7 Posterior to FEF Right',...
    'S6D8 Inferor to FEF Right',...
    'S7D16 Posterior STG/PT Left',...
    'S8D17 Posterior STG/PT Right',...
    'S9D13 IPS3/IPS2/SPL1 Left',...
    'S9D14 IPS3/antIPS/IPS4 Left',...
    'S10D13 IPS3/IPS2/SPL1 Right',...
    'S10D15 IPS3/antIPS/IPS4 Right',...
    'S11D13 Superior to IPS3/IPS2/SPL1 Left',...
    'S11D14 IPS4 Left',...
    'S12D13 Superior to IPS3/IPS2/SPL1 Right',...
    'S12D15 IPS4 Right'};

% index for chnName
chnGrp = {[1 2 3 4 5 6 7 8 9 10] [11 12 13 14 15 16 17 18 19 20]...
    21 22 [23 24 27 28] [25 26 29 30]};

idxLFEF = 1:10;
idxRFEF = 11:20;
idxLSTG = 21;
idxRSTG = 22;
idxLIPS = [23 24 27 28];
idxRIPS = [25 26 29 30];

chnGrpName = {'L FEF/tgPSC','R FEF/tgPSC','L STG','R STG','L IPS','R IPS'};

varPerfHbO = zeros(length(sbjList),length(chnName),length(timePt));
varPerfHbR = zeros(length(sbjList),length(chnName),length(timePt));
varPerfHbT = zeros(length(sbjList),length(chnName),length(timePt));

varPerfBestHbO = zeros(length(sbjList),6,length(timePt));
varPerfBestHbR = zeros(length(sbjList),6,length(timePt));
varPerfBestHbT = zeros(length(sbjList),6,length(timePt));

varPerfSTDHbO = zeros(length(sbjList),length(chnName));
varPerfSTDHbR = zeros(length(sbjList),length(chnName));
varPerfSTDHbT = zeros(length(sbjList),length(chnName));

varPerfHbTGrp = zeros(length(sbjList),size(chnGrp,2));

% tempVarPerfHbO = zeros(length(sbjList),length(chnName),numRept);
% tempVarPerfHbR = zeros(length(sbjList),length(chnName),numRept);
% tempVarPerfHbT = zeros(length(sbjList),length(chnName),numRept);
        
if rerunOp == 0
    % This has never been tested!!! forgot about it!!!
%     for i = 1:length(sbjList)
        
%         sbjNum = sbjList{i};
%         
%         processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
%         
%         fn = [processedDataDir filesep 'perfLinearDiscriminant_LR_Basis1'];
% 
%         load(fn,'performanceArrMultiHbO','performanceArrMultiHbR',...
%             'performanceArrMultiHbT');
% %         % performanceArr = numChns x time
% %         varPerfHbO(i,:) = performanceArrMultiHbO(:,tInd)';
% %         varPerfHbR(i,:) = performanceArrMultiHbR(:,tInd)';
% %         varPerfHbT(i,:) = performanceArrMultiHbT(:,tInd)';
%         load('grpBarChartROINames.mat');
        

%     end

    %fn = [figDir filesep 'grp_Names.mat'];
    load(dataFN);
else
    for i = 1:length(sbjList)
        
        sbjNum = sbjList{i};
        
        disp(sbjNum);
        
        processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
        
        rawDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
        
        if numClasses == 2
            load([processedDataDir filesep 'singleTrialsUpdated_LR_Basis1.mat'],'singleTrialHRFHbOM',...
                'singleTrialHRFHbRM',...
                'singleTrialHRFHbTM','indexMoviesTest');
        else
            load([processedDataDir filesep 'singleTrialsUpdated_Basis1.mat'],'singleTrialHRFHbOM',...
                'singleTrialHRFHbRM',...
                'singleTrialHRFHbTM','indexMoviesTest');
        end
                
        processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
        load([processedDataDir filesep 'intermediateOutputsCombined_Basis1.mat'],'mlActAuto');
        
        if strcmp(sbjNum,'08')
            isSDNew = 0;
        elseif strcmp(sbjNum,'10')
            isSDNew = 0;
        elseif strcmp(sbjNum,'12')
            isSDNew = 1;
        else
            isSDNew = 1;
        end
        
        if strcmp(sbjNum,'15')
            indexMoviesTest = indexMoviesTest(2:end,:);
        end
        
        singleTrialHRFHbTM = offsetTrials(singleTrialHRFHbTM,zeroT);
        
        % after convert2SD2, cut down from 42 chns to 36 chns. Remove old
        % chns
        if ~isSDNew
            [singleTrialHRFHbOM,mlActAutoNew] = convert2SD2(singleTrialHRFHbOM,mlActAuto{1});
            [singleTrialHRFHbRM,~] = convert2SD2(singleTrialHRFHbRM);
            [singleTrialHRFHbTM,~] = convert2SD2(singleTrialHRFHbTM);
        else
            mlActAutoNew = mlActAuto{1};
        end
        
        % terrible coding practice but it works
        %isSDNew = 1;
        
        % after selectRS, cut down from 36 to 30 chns. Remove SS
        [singleTrialHRFHbOM,mlList] = selectRS(singleTrialHRFHbOM,1,mlActAutoNew);
        [singleTrialHRFHbRM,~] = selectRS(singleTrialHRFHbRM,1);
        [singleTrialHRFHbTM,~] = selectRS(singleTrialHRFHbTM,1);
        
        behFN = [rawDataDir filesep 'responses_' sbjNum];
        origIdxFN = [rawDataDir filesep 'movieList_' sbjNum];

        if ~isSDNew
            [singleTrialHRFHbOM,movieIdx] = keepCorrectTrials_Special(sbjNum,singleTrialHRFHbOM,behFN,numClasses,indexMoviesTest);
            [singleTrialHRFHbRM] = keepCorrectTrials_Special(sbjNum,singleTrialHRFHbRM,behFN,numClasses,indexMoviesTest);
            [singleTrialHRFHbTM] = keepCorrectTrials_Special(sbjNum,singleTrialHRFHbTM,behFN,numClasses,indexMoviesTest);
        else
            [singleTrialHRFHbOM,movieIdx] = keepCorrectTrials(sbjNum,singleTrialHRFHbOM,behFN,numClasses,indexMoviesTest);
            [singleTrialHRFHbRM] = keepCorrectTrials(sbjNum,singleTrialHRFHbRM,behFN,numClasses,indexMoviesTest);
            [singleTrialHRFHbTM] = keepCorrectTrials(sbjNum,singleTrialHRFHbTM,behFN,numClasses,indexMoviesTest);
        end

%         if numClasses == 2
%             multipleIndex = (indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1) | (indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1);
%         else
%             multipleIndex = (indexMoviesTest(:,5)==1);
%         end
%         indexMoviesTestMultiple = indexMoviesTest(multipleIndex,:);

%         % offset
%         for i1 = 2:size(singleTrialHRFHbOM,1)
%             for i2 = 1:size(singleTrialHRFHbOM,3)
%                 offset = singleTrialHRFHbOM(i1,zeroT,i2);
%                 %offset = mean(singleTrials(i1,1:zeroT,i2));
%                 singleTrialHRFHbOM(i1,:,i2) = singleTrialHRFHbOM(i1,:,i2) - offset;
% 
%                 offset = singleTrialHRFHbRM(i1,zeroT,i2);
%                 singleTrialHRFHbRM(i1,:,i2) = singleTrialHRFHbRM(i1,:,i2) - offset;
% 
%                 offset = singleTrialHRFHbTM(i1,zeroT,i2);
%                 singleTrialHRFHbTM(i1,:,i2) = singleTrialHRFHbTM(i1,:,i2) - offset;
%             end
%         end
        
        tempVarPerfHbO = zeros(length(chnName),length(timePt),numRept);
        tempVarPerfHbR = zeros(length(chnName),length(timePt),numRept);
        tempVarPerfHbT = zeros(length(chnName),length(timePt),numRept);
        
        idxT = find(timePt == (timeStrtInt+3)*fs);
%         for i3=idx:idx
        for i3=1:length(timePt)
%             temp = cumsum(singleTrialHRFHbTM(:,timePt(i3):timePt(i3)+lenT,:),2);
%             cumsumHbTM = squeeze(temp(:,end,:));
            thisTimePt = timePt(i3);
            for i4=1:numRept
            
                %performanceArr is channels x time
%                 %(trials,timePt,timeLngth,mlActAuto,numClasses,indexMoviesTest)
%                 [performanceArrMultiHbO,~] = trainClassifierLinearDiscriminantfNIRS_Cumsum_DiffStart(singleTrialHRFHbOM,timePt,timeLngth,mlActAutoNew,numClasses,movieIdx);
%                 [performanceArrMultiHbR,~] = trainClassifierLinearDiscriminantfNIRS_Cumsum_DiffStart(singleTrialHRFHbRM,timePt,timeLngth,mlActAutoNew,numClasses,movieIdx);
%                 [performanceArrMultiHbT,~] = trainClassifierLinearDiscriminantfNIRS_Cumsum_DiffStart(singleTrialHRFHbTM,timePt,timeLngth,mlActAutoNew,numClasses,movieIdx);
% 
%                 %we're not computing CI here.
%                 tempVarPerfHbO(:,i4) = performanceArrMultiHbO(:,tInd);
%                 tempVarPerfHbR(:,i4) = performanceArrMultiHbR(:,tInd);
%                 tempVarPerfHbT(:,i4) = performanceArrMultiHbT(:,tInd);
                
                tempVarPerfHbO(:,i3,i4) = trainClassifierLinearDiscriminantfNIRS_Cumsum_DiffStart(singleTrialHRFHbOM,thisTimePt,timeLngth,mlList,numClasses,movieIdx);
                tempVarPerfHbR(:,i3,i4) = trainClassifierLinearDiscriminantfNIRS_Cumsum_DiffStart(singleTrialHRFHbRM,thisTimePt,timeLngth,mlList,numClasses,movieIdx);
                tempVarPerfHbT(:,i3,i4) = trainClassifierLinearDiscriminantfNIRS_Cumsum_DiffStart(singleTrialHRFHbTM,thisTimePt,timeLngth,mlList,numClasses,movieIdx);
            
            end
        end
        
        varPerfHbO(i,:,:) = squeeze(mean(tempVarPerfHbO,3));
        varPerfHbR(i,:,:) = squeeze(mean(tempVarPerfHbR,3));
        varPerfHbT(i,:,:) = squeeze(mean(tempVarPerfHbT,3));
        
        [~,bestChn] = max(varPerfHbO(i,idxLFEF,idxT),[],2);
        varPerfBestHbO(i,1,:) = varPerfHbO(i,idxLFEF(bestChn),:);
        
        [~,bestChn] = max(varPerfHbR(i,idxLFEF,idxT),[],2);
        varPerfBestHbR(i,1,:) = varPerfHbR(i,idxLFEF(bestChn),:);
        
        [~,bestChn] = max(varPerfHbT(i,idxLFEF,idxT),[],2);
        varPerfBestHbT(i,1,:) = varPerfHbT(i,idxLFEF(bestChn),:);
        
        [~,bestChn] = max(varPerfHbO(i,idxRFEF,idxT),[],2);
        varPerfBestHbO(i,2,:) = varPerfHbO(i,idxRFEF(bestChn),:);
        
        [~,bestChn] = max(varPerfHbR(i,idxRFEF,idxT),[],2);
        varPerfBestHbR(i,2,:) = varPerfHbR(i,idxRFEF(bestChn),:);
        
        [~,bestChn] = max(varPerfHbT(i,idxRFEF,idxT),[],2);
        varPerfBestHbT(i,2,:) = varPerfHbT(i,idxRFEF(bestChn),:);
        
        [~,bestChn] = max(varPerfHbO(i,idxLSTG,idxT),[],2);
        varPerfBestHbO(i,3,:) = varPerfHbO(i,idxLSTG(bestChn),:);
        
        [~,bestChn] = max(varPerfHbR(i,idxLSTG,idxT),[],2);
        varPerfBestHbR(i,3,:) = varPerfHbR(i,idxLSTG(bestChn),:);
        
        [~,bestChn] = max(varPerfHbT(i,idxLSTG,idxT),[],2);
        varPerfBestHbT(i,3,:) = varPerfHbT(i,idxLSTG(bestChn),:);
        
        [~,bestChn] = max(varPerfHbO(i,idxRSTG,idxT),[],2);
        varPerfBestHbO(i,4,:) = varPerfHbO(i,bestChn,:);
        
        [~,bestChn] = max(varPerfHbR(i,idxRSTG,idxT),[],2);
        varPerfBestHbR(i,4,:) = varPerfHbR(i,idxRSTG(bestChn),:);
        
        [~,bestChn] = max(varPerfHbT(i,idxRSTG,idxT),[],2);
        varPerfBestHbT(i,4,:) = varPerfHbT(i,idxRSTG(bestChn),:);
        
        [~,bestChn] = max(varPerfHbO(i,idxLIPS,idxT),[],2);
        varPerfBestHbO(i,5,:) = varPerfHbO(i,idxLIPS(bestChn),:);
        
        [~,bestChn] = max(varPerfHbR(i,idxLIPS,idxT),[],2);
        varPerfBestHbR(i,5,:) = varPerfHbR(i,idxLIPS(bestChn),:);
        
        [~,bestChn] = max(varPerfHbT(i,idxLIPS,idxT),[],2);
        varPerfBestHbT(i,5,:) = varPerfHbT(i,idxLIPS(bestChn),:);
        
        [~,bestChn] = max(varPerfHbO(i,idxRIPS,idxT),[],2);
        varPerfBestHbO(i,6,:) = varPerfHbO(i,idxRIPS(bestChn),:);
        
        [~,bestChn] = max(varPerfHbR(i,idxRIPS,idxT),[],2);
        varPerfBestHbR(i,6,:) = varPerfHbR(i,idxRIPS(bestChn),:);
        
        [~,bestChn] = max(varPerfHbT(i,idxRIPS,idxT),[],2);
        varPerfBestHbT(i,6,:) = varPerfHbT(i,idxRIPS(bestChn),:);
        
%         varPerfSTDHbO(i,:,:) = squeeze(std(tempVarPerfHbO,0,3)');
%         varPerfSTDHbR(i,:,:) = squeeze(std(tempVarPerfHbR,0,3)');
%         varPerfSTDHbT(i,:,:) = squeeze(std(tempVarPerfHbT,0,3)');
    end
    
    %fn = [figDir filesep 'grp_Names'];
    save(dataFN);
    
end

% for i = 1:size(chnGrp,2)
%     varPerfHbTGrp(:,i) = max(varPerfHbT(:,chnGrp{i}),[],2);
% end

grpMeanPerfHbO = squeeze(mean(varPerfBestHbO,1));
grpMeanPerfHbR = squeeze(mean(varPerfBestHbR,1));
grpMeanPerfHbT = squeeze(mean(varPerfBestHbT,1));

grpSTDPerfHbO = squeeze(std(varPerfBestHbO,[],1));
grpSTDPerfHbR = squeeze(std(varPerfBestHbR,[],1));
grpSTDPerfHbT = squeeze(std(varPerfBestHbT,[],1));

cmap = jet(6);

figure('units','normalized','outerposition',[0 0 1 1]); hold on;

subplot(1,3,1);
for j = 1:6
    errorbar(timePt./fs-2,grpMeanPerfHbO(j,:),grpSTDPerfHbO(j,:),'Color',cmap(j,:));hold on;
end
yyaxis left;
ylim([0 1]);
ylabel('Accuracy');
yyaxis right;
if numClasses == 2
    ylim([-max(itf) max(itf)]);
    yticks([-max(itf) 0 max(itf)]);
    yticklabels({num2str(max(itf)) num2str(0) num2str(max(itf))});
else
    ylim([min(itf) max(itf)]);
    yticks([min(itf) 0 max(itf)]);
    yticks(linspace(min(itf), max(itf),6));
    tempIdx = linspace(1,100,6);
    yticklabels([-itf(2) -itf(round(tempIdx(2))) ...
        itf(round(tempIdx(3))) itf(round(tempIdx(4))) itf(round(tempIdx(5))) ...
        itf(round(tempIdx(6)))]);
    %yticklabels({num2str(max(itf)) num2str(0) num2str(max(itf))});  
end
ylabel('ITR bits/trial');
xlim([0 7]);
title(sprintf('Group Avg of Best Chn: Δ[HbO] Multi'));
legend(chnGrpName);
xlabel('Time [s]');
hold off;

subplot(1,3,2);
for j = 1:6
    errorbar(timePt./fs-2,grpMeanPerfHbR(j,:),grpSTDPerfHbR(j,:),'Color',cmap(j,:));hold on;
end
yyaxis left;
ylim([0 1]);
ylabel('Accuracy');
yyaxis right;
if numClasses == 2
    ylim([-max(itf) max(itf)]);
    yticks([-max(itf) 0 max(itf)]);
    yticklabels({num2str(max(itf)) num2str(0) num2str(max(itf))});
else
    ylim([min(itf) max(itf)]);
    yticks([min(itf) 0 max(itf)]);
    yticks(linspace(min(itf), max(itf),6));
    tempIdx = linspace(1,100,6);
    yticklabels([-itf(2) -itf(round(tempIdx(2))) ...
        itf(round(tempIdx(3))) itf(round(tempIdx(4))) itf(round(tempIdx(5))) ...
        itf(round(tempIdx(6)))]);
    %yticklabels({num2str(max(itf)) num2str(0) num2str(max(itf))});  
end
ylabel('ITR bits/trial');
xlim([0 7]);
title(sprintf('Group Avg of Best Chn: Δ[HbR] Multi'));
legend(chnGrpName);
xlabel('Time [s]');
hold off;

subplot(1,3,3);
for j = 1:6
    errorbar(timePt./fs-2,grpMeanPerfHbT(j,:),grpSTDPerfHbT(j,:),'Color',cmap(j,:));hold on;
end
yyaxis left;
ylim([0 1]);
ylabel('Accuracy');
yyaxis right;
if numClasses == 2
    ylim([-max(itf) max(itf)]);
    yticks([-max(itf) 0 max(itf)]);
    yticklabels({num2str(max(itf)) num2str(0) num2str(max(itf))});
else
    ylim([min(itf) max(itf)]);
    yticks([min(itf) 0 max(itf)]);
    yticks(linspace(min(itf), max(itf),6));
    tempIdx = linspace(1,100,6);
    yticklabels([-itf(2) -itf(round(tempIdx(2))) ...
        itf(round(tempIdx(3))) itf(round(tempIdx(4))) itf(round(tempIdx(5))) ...
        itf(round(tempIdx(6)))]);
    %yticklabels({num2str(max(itf)) num2str(0) num2str(max(itf))});  
end
ylabel('ITR bits/trial');
xlim([0 7]);
title(sprintf('Group Avg of Best Chn: Δ[HbT] Multi'));
legend(chnGrpName);
xlabel('Time [s]');
hold off;

if saveOp == 1
    %fn = sprintf('Sbj%s%ss%s',num2str(sbjNum),featureName,num2str(timePt(tInd)));
    
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
end

end
