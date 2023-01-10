% split into single and multi groups.
% new. Use trials from dcNew var
% option to re-run training and classifying or pull from saved data
% sbjList is cell array of experiment "number". Ex: sbjList = {'01','02'};
% sbjList = {'08','10','12','13','14','15'};tInd = 8;featureName = 'Cumulative Sum';numRept = 5;numClasses = 2;saveOp = 0;rerunOp = 1;
% For Figure Draft 4. Now exclude rejected channel from calculating avg CA
% performance.
% Add sig test

function grp_Plot_Subset_LOFO_CV(sbjGrp,numClasses,fefOnlyOp,saveOp)

fs = 50;
%timeLen = [0.06 0.1 0.2 0.5 1 1.5 2 3 4 5].*fs;
timePt = (0:0.25*fs:7*fs);
deltaLen = [5 11 25 49];

if sbjGrp == 1
    sbjList = {'08','12','13','16','19','21','24'};
    sbjGrpStr = 'High';
elseif sbjGrp == 2
    sbjList = {'08','12','13','14','15','16','19','21','22','23','24','25'};
    sbjGrpStr = 'All';
elseif sbjGrp == 0
    sbjList = {'14','15','22','23','25'};
    sbjGrpStr = 'Low';
end

figSaveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll\Figures\Classification';
dataSaveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll';

fs = 50;
%timePt = (-2*fs:0.25*fs:5*fs);

% w/o SS, just use first 20 chns.
if fefOnlyOp
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
        'S6D8 Inferor to FEF Right'};
else
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
end

performance_All = zeros(10*length(sbjList),length(timePt));
performance_Subset  = zeros(10*length(sbjList),length(timePt));
performance_Delta = zeros(10*length(sbjList),length(deltaLen));
performance_Delta1  = zeros(10*length(sbjList),length(deltaLen));

chnList = zeros(length(sbjList),length(chnName));
pVal = zeros(1,length(chnName));
pvalStar = cell(1,length(chnName));

%thisT = 2*fs;
%tidx = abs(timePt-thisT)<10*eps("double");
        
for i = 1:length(sbjList)

    sbjNum = sbjList{i};

    processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];

    if numClasses == 2
        if fefOnlyOp
            fileName = [processedDataDir filesep 'performance_LOFO_Sbst_LR_NRept.mat'];
        else
            fileName = [processedDataDir filesep 'performance_LOFO_Sbst_LR_NRept.mat'];
        end
    end
    
    fnStat = [processedDataDir filesep 'summaryStat_1.5.mat'];

    % only regular channels. Folds x channel x time
    load(fileName,'performanceLDALedoitHbT_All_FS','performanceLDALedoitHbT_Subset');
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
    if fefOnlyOp
        ml690(21:end) = [];
        ml870(21:end) = [];
    end

    chnList(i,:) = ml690&ml870;
    
    startIdx = (i-1)*10+1;
    endIdx = i*10;
    
    performance_All(startIdx:endIdx,:) = performanceLDALedoitHbT_All_FS;
    performance_Subset(startIdx:endIdx,:) = performanceLDALedoitHbT_Subset;
    
end

% if fefOnlyOp
%     fn = 'grp_HighPerf_3FLO_FEF';
% else
%     fn = 'grp_HighPerf_3FLO';
% end

%save([dataSaveDir filesep fn],'allFEF_CA','allFEF_3FLO_CA','allIPS_CA','allIPS_3FLO_CA');

performance_All_Mean = mean(performance_All,1);
performance_Subset_Mean = mean(performance_Subset,1);
% performance_Delta_Mean = mean(performance_Delta,1);
% performance_Delta1_Mean = mean(performance_Delta1,1);

performanceLDALedoitHbT_All_CI = zeros(length(timePt),2);
performanceLDALedoitHbT_Subset_CI = zeros(length(timePt),2);
% performanceLDALedoitHbT_Delta_CI = zeros(length(deltaLen),2);
% performanceLDALedoitHbT_Delta1_CI = zeros(length(deltaLen),2);

conf = 0.95;nrept = 10;kFold = 5;

for i = 1:length(timePt)
    
    temp = performance_All(:,i);
    se = std(temp(:))/sqrt(nrept);
    performanceLDALedoitHbT_All_CI(i,:) = calcCITDist(nrept,se,conf);

    temp = performance_Subset(:,i);
    se = std(temp(:))/sqrt(nrept);
    performanceLDALedoitHbT_Subset_CI(i,:) = calcCITDist(nrept,se,conf);
    
end

% for i = 1:length(deltaLen)
% 
%     temp = performance_Delta(:,i);
%     se = std(temp(:))/sqrt(nrept*kFold);
%     performanceLDALedoitHbT_Delta_CI(i,:) = calcCITDist(nrept*kFold,se,conf);
% 
%     temp = performance_Delta1(:,i);
%     se = std(temp(:))/sqrt(nrept*kFold);
%     performanceLDALedoitHbT_Delta1_CI(i,:) = calcCITDist(nrept*kFold,se,conf);
%     
% end

figure();
% Absolute
%subplot(1,2,1);hold on;
errorbar(timePt./fs-2,performance_All_Mean,performanceLDALedoitHbT_All_CI(:,2));hold on;
errorbar(timePt./fs-2,performance_Subset_Mean,performanceLDALedoitHbT_Subset_CI(:,2));hold on;
% errorbar(deltaLen./fs,performance_Delta_Mean,performanceLDALedoitHbT_Delta_CI(:,2));hold on;
% errorbar(deltaLen./fs,performance_Delta1_Mean,performanceLDALedoitHbT_Delta1_CI(:,2));hold on;
legend('All','Subset');

if saveOp == 1
    
    if fefOnlyOp
        fn = sprintf('Grp_Subset_HbT');
    else
        fn = sprintf('Grp_Subset_HbT');
    end
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
    
end

end

