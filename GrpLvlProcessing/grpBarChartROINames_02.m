% split into single and multi groups.
% new. Use trials from dcNew var
% option to re-run training and classifying or pull from saved data
% sbjList is cell array of experiment "number". Ex: sbjList = {'01','02'};
% sbjList = {'08','10','12','13','14','15'};tInd = 8;featureName = 'Cumulative Sum';numRept = 5;numClasses = 2;saveOp = 0;rerunOp = 1;
% 

function grpBarChartROINames_02(sbjGrp,numClasses,rejChnOp,rerunOp,saveOp)

numRept = 5;

if sbjGrp == 1
    sbjList = {'08','12','13','16','19','21','24'};
    sbjGrpStr = 'HighPerf';
elseif sbjGrp == 2
    %sbjList = {'08','10''12','13','14','15','16','17','18','19','20','21','22','23','24','25'};
    sbjList = {'08','12','13','14','15','16','19','21','22','23','24','25'};
    sbjGrpStr = 'AllSbjs';
elseif sbjGrp == 0
    sbjList = {'14','15','22','23','25'};
    sbjGrpStr = 'LowPerf';
end

figSaveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll\Figures\Classification';

fs = 50;
timePt = (-2*fs:0.25*fs:5*fs);
zeroT = 2*fs;
timeLngth = 1*fs;
lenT = 1*fs;
numChn = 30;
timeStrtInt = 3;

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

% 

varPerfHbO = zeros(length(sbjList),length(chnName));
varPerfHbR = zeros(length(sbjList),length(chnName));
varPerfHbT = zeros(length(sbjList),length(chnName));

varPerfSTDHbO = zeros(length(sbjList),length(chnName));
varPerfSTDHbR = zeros(length(sbjList),length(chnName));
varPerfSTDHbT = zeros(length(sbjList),length(chnName));

% tempVarPerfHbO = zeros(length(sbjList),length(chnName),numRept);
% tempVarPerfHbR = zeros(length(sbjList),length(chnName),numRept);
% tempVarPerfHbT = zeros(length(sbjList),length(chnName),numRept);
    
chnList = zeros(length(sbjList),30);
chnSbjCnt = zeros(1,30);

if rerunOp == 0
    for i = 1:length(sbjList)

        sbjNum = sbjList{i};

        processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];

        if numClasses==2
            % from preprocessFNIRS06_CV_GLM_ssBeta_SingleChn_MultiOnly.m
            fn = [processedDataDir filesep 'performance_GLM_CV_SSBeta_SingleChn_LR_RejTr_SNR_1.5_OrigHomer3.mat'];
        elseif numClasses == 3
            fn = [processedDataDir filesep 'performance_GLM_CV_SSBeta_SingleChn_RejTr_SNR_1.5_OrigHomer3.mat'];
        end

        fnStat = [processedDataDir filesep 'summaryStat_1.5.mat'];

    %         load(fn,'performanceArrMultiHbO','performanceArrMultiHbR',...
    %             'performanceArrMultiHbT');
    %         % performanceArr = numChns x time
    %         varPerfHbO(i,:) = performanceArrMultiHbO(:,tInd)';
    %         varPerfHbR(i,:) = performanceArrMultiHbR(:,tInd)';
    %         varPerfHbT(i,:) = performanceArrMultiHbT(:,tInd)';
        % only regular channels. Folds x channel x time
        load(fn,'performanceLDAHbO','performanceLDAHbR','performanceLDAHbT');
        load(fnStat,'mlActAuto');

        ml690 = mlActAuto{1}(1:length(mlActAuto{1})/2);
        ml870 = mlActAuto{1}(length(mlActAuto{1})/2+1:length(mlActAuto{1}));

        % convert old probe (42 chns) to new probe (36 chns)
        if strcmp(sbjNum,'08')
            origIdxToRemove = [29,33,39,40,41,42];
            ml690(origIdxToRemove) = [];
            ml870(origIdxToRemove) = [];
        end

        % filter out SS channels (from 36 total to 30 total)
        ssIdx = [7,22,24,26,29,32];
        ml690(ssIdx) = [];
        ml870(ssIdx) = [];

        chnList(i,:) = ml690&ml870;

        thisT = timeStrtInt*fs;
        tidx = abs(timePt-thisT)<10*eps("double");

        %tempHbO = squeeze(sum(performanceLDAHbO,1))/size(performanceLDAHbO,1);
        tempHbO = sum(squeeze(performanceLDAHbO(:,:,tidx)),1)/size(performanceLDAHbO,1);
        tempHbR = sum(squeeze(performanceLDAHbR(:,:,tidx)),1)/size(performanceLDAHbR,1);
        tempHbT = sum(squeeze(performanceLDAHbT(:,:,tidx)),1)/size(performanceLDAHbT,1);

        for i2 = 1:size(chnList,2)
            if rejChnOp && chnList(i,i2)

                chnSbjCnt(i2) = chnSbjCnt(i2) + 1;
                varPerfHbO(i,i2) = tempHbO(i2);
                varPerfHbR(i,i2) = tempHbR(i2);
                varPerfHbT(i,i2) = tempHbT(i2);

            elseif ~rejChnOp
                chnSbjCnt(i2) = chnSbjCnt(i2) + 1;
                varPerfHbO(i,i2) = tempHbO(i2);
                varPerfHbR(i,i2) = tempHbR(i2);
                varPerfHbT(i,i2) = tempHbT(i2);
            end
        end
    end
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
            load([processedDataDir filesep 'singleTrialsUpdated.mat'],'singleTrialHRFHbOM',...
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
        
        tempVarPerfHbO = zeros(length(chnName),numRept);
        tempVarPerfHbR = zeros(length(chnName),numRept);
        tempVarPerfHbT = zeros(length(chnName),numRept);
        
        idx = find(timePt == (timeStrtInt+2)*fs);
        for i3=idx:idx
        %for i3=1:length(timePt)
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
                
                %tempVarPerfHbT(i3,i4) = train_RLDA_Ledoit(cumsumHbTM,mlList,numChn,numClasses,movieIdx);
                tempVarPerfHbT(:,i4) = trainClassifierLinearDiscriminantfNIRS_Cumsum_DiffStart(singleTrialHRFHbTM,thisTimePt,timeLngth,mlList,numClasses,movieIdx);
            
            end
        end
        
%         varPerfHbO(i,:) = mean(tempVarPerfHbO,2)';
%         varPerfHbR(i,:) = mean(tempVarPerfHbR,2)';
        varPerfHbT(i,:) = mean(tempVarPerfHbT,2)';
        
%         varPerfSTDHbO(i,:) = std(tempVarPerfHbO,0,2)';
%         varPerfSTDHbR(i,:) = std(tempVarPerfHbR,0,2)';
        varPerfSTDHbT(i,:) = std(tempVarPerfHbT,0,2)';
    end
end

% grpMeanPerfHbO = mean(nonzeros(varPerfHbO));
% grpMeanPerfHbR = mean(nonzeros(varPerfHbR));
grpMeanPerfHbT = mean(nonzeros(varPerfHbT));

for i = 1:length(chnName)
%     grpMeanPerfHbO(i) = mean(nonzeros(varPerfHbO(:,i)));
%     grpMeanPerfHbR(i) = mean(nonzeros(varPerfHbR(:,i)));
    grpMeanPerfHbT(i) = mean(nonzeros(varPerfHbT(:,i)));
end

%figure();
X = categorical(chnName);
X = reordercats(X,chnName);

save('grpBarChartROINames.mat');

colorCodes = loadDefaultColors(1);

% figure(1);hold on;
% for i1 = 1:length(chnName)
%     tempY = varPerfHbO(:,i1);
%     tempX = i1*ones(1,length(sbjList));
%     for i2 = 1:length(sbjList)
%         plot(tempX(i2),tempY(i2),'o','MarkerFaceColor',colorCodes(i2,:),'MarkerEdgeColor',colorCodes(i2,:),'MarkerSize',5);
%     end
%     plot(i1,grpMeanPerfHbO(1,i1),'*','Color',colorCodes(length(sbjList)+1,:),'MarkerSize',10);
% end
% ax = gca;
% ax.XTick = 1:1:length(chnName);
% ax.XTickLabels = chnName;
% xtickangle(45);
% title('HbO: Sum. 1.5s interval. At 4s');
% hold off;
% 
% figure(2);hold on;
% for i1 = 1:length(chnName)
%     tempY = varPerfHbR(:,i1);
%     tempX = i1*ones(1,length(sbjList));
%     for i2 = 1:length(sbjList)
%         plot(tempX(i2),tempY(i2),'o','MarkerFaceColor',colorCodes(i2,:),'MarkerEdgeColor',colorCodes(i2,:),'MarkerSize',5);
%     end
%     plot(i1,grpMeanPerfHbR(1,i1),'*','Color',colorCodes(length(sbjList)+1,:),'MarkerSize',10);
% end
% ax = gca;
% ax.XTick = 1:1:length(chnName);
% ax.XTickLabels = chnName;
% xtickangle(45);
% title('HbR: Sum. 1.5s interval. At 4s');
% hold off;

figure(3);hold on;
set(gcf,'Position',get(0,'Screensize'));
for i1 = 1:length(chnName)
    tempY = varPerfHbT(:,i1);
    tempX = i1*ones(1,length(sbjList));
    for i2 = 1:length(sbjList)
        %plot(tempX(i2),tempY(i2),'o','MarkerFaceColor',colorCodes(i2,:),'MarkerEdgeColor',colorCodes(i2,:),'MarkerSize',5);
        plot(tempX(i2),tempY(i2),'o','MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410],'MarkerSize',5);
    end
    plot(i1,grpMeanPerfHbT(1,i1),'*','Color',[0.8500 0.3250 0.0980],'MarkerSize',10);
end
ylim([0.2 1]);
ax = gca;
ax.XTick = 1:1:length(chnName);
ax.XTickLabels = chnName;
xtickangle(45);
title(sprintf('HbT: Sum. %ss interval. At %ss',num2str(timeLngth/fs),num2str(timeStrtInt)));
hold off;

if saveOp == 1
    %fn = sprintf('Sbj%s%ss%s',num2str(sbjNum),featureName,num2str(timePt(tInd)));
    fn = sprintf('Grp_CVPerfByROINames_HbT');
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
end


end

