% split into single and multi groups.
% new. Use trials from dcNew var
% option to re-run training and classifying or pull from saved data
% sbjList is cell array of experiment "number". Ex: sbjList = {'01','02'};
% sbjList = {'08','10','12','13','14','15'};tInd = 8;featureName = 'Cumulative Sum';numRept = 5;numClasses = 2;saveOp = 0;rerunOp = 1;
% 

function grpBarChartROINames(sbjList,tInd,featureName,numRept,numClasses,rerunOp,saveOp)

figSaveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll\Figures\Classification';

fs = 50;
timePt = (0:0.5*fs:12*fs)+2*fs;
zeroT = 2*fs;
timeLngth = 1*fs;

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
        
if rerunOp == 0
    % This has never been tested!!! forgot about it!!!
    for i = 1:length(sbjList)
        
        sbjNum = sbjList{i};
        
        processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
        
        fn = [processedDataDir filesep 'perfLinearDiscriminant_LR_Basis4'];
%         load(fn,'performanceArrMultiHbO','performanceArrMultiHbR',...
%             'performanceArrMultiHbT');
%         % performanceArr = numChns x time
%         varPerfHbO(i,:) = performanceArrMultiHbO(:,tInd)';
%         varPerfHbR(i,:) = performanceArrMultiHbR(:,tInd)';
%         varPerfHbT(i,:) = performanceArrMultiHbT(:,tInd)';
        load('grpBarChartROINames.mat');
        

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
        [singleTrialHRFHbOM,mlActAutoNew] = selectRS(singleTrialHRFHbOM,1,mlActAutoNew);
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

        % offset
        for i1 = 2:size(singleTrialHRFHbOM,1)
            for i2 = 1:size(singleTrialHRFHbOM,3)
                offset = singleTrialHRFHbOM(i1,zeroT,i2);
                %offset = mean(singleTrials(i1,1:zeroT,i2));
                singleTrialHRFHbOM(i1,:,i2) = singleTrialHRFHbOM(i1,:,i2) - offset;

                offset = singleTrialHRFHbRM(i1,zeroT,i2);
                singleTrialHRFHbRM(i1,:,i2) = singleTrialHRFHbRM(i1,:,i2) - offset;

                offset = singleTrialHRFHbTM(i1,zeroT,i2);
                singleTrialHRFHbTM(i1,:,i2) = singleTrialHRFHbTM(i1,:,i2) - offset;
            end
        end
        
        tempVarPerfHbO = zeros(length(chnName),numRept);
        tempVarPerfHbR = zeros(length(chnName),numRept);
        tempVarPerfHbT = zeros(length(chnName),numRept);
        
        for i3=1:numRept
            
            %performanceArr is channels x time
            % (trials,timePt,timeLngth,mlActAuto,numClasses,indexMoviesTest)
            [performanceArrMultiHbO,~] = trainClassifierLinearDiscriminantfNIRS_Cumsum_DiffStart(singleTrialHRFHbOM,timePt,timeLngth,mlActAutoNew,numClasses,movieIdx);
            [performanceArrMultiHbR,~] = trainClassifierLinearDiscriminantfNIRS_Cumsum_DiffStart(singleTrialHRFHbRM,timePt,timeLngth,mlActAutoNew,numClasses,movieIdx);
            [performanceArrMultiHbT,~] = trainClassifierLinearDiscriminantfNIRS_Cumsum_DiffStart(singleTrialHRFHbTM,timePt,timeLngth,mlActAutoNew,numClasses,movieIdx);
                        
            % we're not computing CI here.
            tempVarPerfHbO(:,i3) = performanceArrMultiHbO(:,tInd);
            tempVarPerfHbR(:,i3) = performanceArrMultiHbR(:,tInd);
            tempVarPerfHbT(:,i3) = performanceArrMultiHbT(:,tInd);
            
        end
        
        varPerfHbO(i,:) = mean(tempVarPerfHbO,2)';
        varPerfHbR(i,:) = mean(tempVarPerfHbR,2)';
        varPerfHbT(i,:) = mean(tempVarPerfHbT,2)';
        
        varPerfSTDHbO(i,:) = std(tempVarPerfHbO,0,2)';
        varPerfSTDHbR(i,:) = std(tempVarPerfHbR,0,2)';
        varPerfSTDHbT(i,:) = std(tempVarPerfHbT,0,2)';
    end
end

grpMeanPerfHbO = mean(nonzeros(varPerfHbO));
grpMeanPerfHbR = mean(nonzeros(varPerfHbR));
grpMeanPerfHbT = mean(nonzeros(varPerfHbT));

for i = 1:length(chnName)
    grpMeanPerfHbO(i) = mean(nonzeros(varPerfHbO(:,i)));
    grpMeanPerfHbR(i) = mean(nonzeros(varPerfHbR(:,i)));
    grpMeanPerfHbT(i) = mean(nonzeros(varPerfHbT(:,i)));
end

figure();
X = categorical(chnName);
X = reordercats(X,chnName);

save('grpBarChartROINames.mat');

colorCodes = loadDefaultColors();

figure(1);hold on;
for i1 = 1:length(chnName)
    tempY = varPerfHbO(:,i1);
    tempX = i1*ones(1,length(sbjList));
    for i2 = 1:length(sbjList)
        plot(tempX(i2),tempY(i2),'o','MarkerFaceColor',colorCodes(i2,:),'MarkerEdgeColor',colorCodes(i2,:),'MarkerSize',5);
    end
    plot(i1,grpMeanPerfHbO(1,i1),'*','Color',colorCodes(length(sbjList)+1,:),'MarkerSize',10);
end
ax = gca;
ax.XTick = 1:1:length(chnName);
ax.XTickLabels = chnName;
xtickangle(45);
title('HbO: Sum. 1.5s interval. At 4s');
hold off;

figure(2);hold on;
for i1 = 1:length(chnName)
    tempY = varPerfHbR(:,i1);
    tempX = i1*ones(1,length(sbjList));
    for i2 = 1:length(sbjList)
        plot(tempX(i2),tempY(i2),'o','MarkerFaceColor',colorCodes(i2,:),'MarkerEdgeColor',colorCodes(i2,:),'MarkerSize',5);
    end
    plot(i1,grpMeanPerfHbR(1,i1),'*','Color',colorCodes(length(sbjList)+1,:),'MarkerSize',10);
end
ax = gca;
ax.XTick = 1:1:length(chnName);
ax.XTickLabels = chnName;
xtickangle(45);
title('HbR: Sum. 1.5s interval. At 4s');
hold off;

figure(3);hold on;
for i1 = 1:length(chnName)
    tempY = varPerfHbT(:,i1);
    tempX = i1*ones(1,length(sbjList));
    for i2 = 1:length(sbjList)
        plot(tempX(i2),tempY(i2),'o','MarkerFaceColor',colorCodes(i2,:),'MarkerEdgeColor',colorCodes(i2,:),'MarkerSize',5);
    end
    plot(i1,grpMeanPerfHbT(1,i1),'*','Color',colorCodes(length(sbjList)+1,:),'MarkerSize',10);
end
ax = gca;
ax.XTick = 1:1:length(chnName);
ax.XTickLabels = chnName;
xtickangle(45);
title('HbT: Sum. 1.5s interval. At 4s');
hold off;

if saveOp == 1
    fn = sprintf('Sbj%s%ss%s',num2str(sbjNum),featureName,num2str(timePt(tInd)));
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
end


end

