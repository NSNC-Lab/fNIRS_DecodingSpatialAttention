% split into single and multi groups.
% new. Use trials from dcNew var
% option to re-run training and classifying or pull from saved data
% sbjList is cell array of experiment "number". Ex: sbjList = {'01','02'};
% sbjList = {'08','10','12','13','14','15'};tInd = 8;featureName = 'Cumulative Sum';numRept = 5;numClasses = 2;saveOp = 0;rerunOp = 1;
% For Figure Draft 4. Now exclude rejected channel from calculating avg CA
% performance.
% Add sig test

function grp_PlotFInteraction_Ind(sbjGrp,numClasses,fefOnlyOp,saveOp)

fs = 50;
timeLen = [0.1 0.2 0.5 1 1.5 2 3 4 5].*fs;

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

diff_FEF_HbT = zeros(3,length(sbjList),length(timeLen));
diff_IPS_HbT = zeros(3,length(sbjList),length(timeLen));

allFEF_CA = zeros(length(sbjList),length(timeLen));
allFEF_Ind_CA  = zeros(3,length(sbjList),length(timeLen));
allIPS_CA = zeros(length(sbjList),length(timeLen));
allIPS_Ind_CA  = zeros(3,length(sbjList),length(timeLen));

chnList = zeros(length(sbjList),length(chnName));
pVal = zeros(1,length(chnName));
pvalStar = cell(1,length(chnName));

%thisT = 2*fs;
%tidx = abs(timePt-thisT)<10*eps("double");
        
for i = 1:length(sbjList)

    sbjNum = sbjList{i};

    processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];

    if numClasses == 2
        fn = [processedDataDir filesep 'performance_3FC_FEF_LR.mat'];
    elseif numClasses == 3
        fn = [processedDataDir filesep 'performance_3FC_FEF.mat'];
    end
    
    fnStat = [processedDataDir filesep 'summaryStat_1.5.mat'];

    % only regular channels. Folds x channel x time
    load(fn,'performanceLDALedoitHbT_FEF_Comb','performanceLDALedoitHbT_FEF_Ind',...
        'performanceLDALedoitHbT_IPS_Comb','performanceLDALedoitHbT_IPS_Ind');
    load(fnStat,'mlActAuto');
    
%     ml690 = mlActAuto{1}(1:length(mlActAuto{1})/2);
%     ml870 = mlActAuto{1}(length(mlActAuto{1})/2+1:length(mlActAuto{1}));
% 
%     % convert old probe (42 chns) to new probe (36 chns)
%     if strcmp(sbjNum,'08')
%         origIdxToRemove = [29,33,39,40,41,42];
%         ml690(origIdxToRemove) = [];
%         ml870(origIdxToRemove) = [];
%     end
% 
%     % filter out SS channels (from 36 total to 30 total)
%     ssIdx = [7,22,24,26,29,32];
%     ml690(ssIdx) = [];
%     ml870(ssIdx) = [];
%     if fefOnlyOp
%         ml690(21:end) = [];
%         ml870(21:end) = [];
%     end
% 
%     chnList(i,:) = ml690&ml870;
    
    allFEF_CA(i,:) = mean(performanceLDALedoitHbT_FEF_Comb,1);
    allFEF_Ind_CA(:,i,:) = squeeze(mean(performanceLDALedoitHbT_FEF_Ind,2));
    allIPS_CA(i,:) = mean(performanceLDALedoitHbT_IPS_Comb,1);
    allIPS_Ind_CA(:,i,:) = squeeze(mean(performanceLDALedoitHbT_IPS_Ind,2));
    
    % Show both absolute and diff (relative)
    for i2 = 1:3
        diff_FEF_HbT(i2,i,:) = squeeze(allFEF_CA(i,:)) - squeeze(allFEF_Ind_CA(i2,i,:))';
        diff_IPS_HbT(i2,i,:) = squeeze(allIPS_CA(i,:)) - squeeze(allIPS_Ind_CA(i2,i,:))';
    end
    
end

grp_Diff_FEF_Hbt_Avg = mean(diff_FEF_HbT,2)*100;
grp_Diff_IPS_Hbt_Avg = mean(diff_IPS_HbT,2)*100;

fn = 'grp_HighPerf_FInteraction_Ind';

save([dataSaveDir filesep fn],'allFEF_CA','allFEF_Ind_CA','allIPS_CA','allIPS_Ind_CA');

performanceLDALedoitHbTAll_CI = zeros(length(timeLen),2);
performanceLDALedoitHbTAll_Ind_CI = zeros(3,length(timeLen),2);
performanceLDALedoitHbTIPS_CI = zeros(length(timeLen),2);
performanceLDALedoitHbTIPS_Ind_CI = zeros(3,length(timeLen),2);

conf = 0.95;nrept = 10;kFold = 5;

for i = 1:length(timeLen)
    
    temp = performanceLDALedoitHbT_FEF_Comb(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDALedoitHbTAll_CI(i,:) = calcCITDist(nrept*kFold,se,conf);

    temp = performanceLDALedoitHbT_IPS_Comb(:,i);
    se = std(temp(:))/sqrt(nrept*kFold);
    performanceLDALedoitHbTIPS_CI(i,:) = calcCITDist(nrept*kFold,se,conf);

    for i2 = 1:3
        
        temp = performanceLDALedoitHbT_FEF_Ind(i2,:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performanceLDALedoitHbTAll_Ind_CI(i2,i,:) = calcCITDist(nrept*kFold,se,conf);

        temp = performanceLDALedoitHbT_IPS_Ind(i2,:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performanceLDALedoitHbTIPS_Ind_CI(i2,i,:) = calcCITDist(nrept*kFold,se,conf);
        
    end
    
end

% for i = 1:size(timeLen,2)
%     
%     x_temp = diff_LOFO_HbT(find(diff_LOFO_HbT(:,i)),i);
%     se_temp = std(diff_LOFO_HbT(find(diff_LOFO_HbT(:,i)),i))/sqrt(chnSbjCnt(i));
%     
%     t_temp = (mean(x_temp))/se_temp;
%     
%     pVal(i) = 1-tcdf(t_temp,chnSbjCnt(i)-1);
%     
%     if pVal(i)<0.35 && pVal(i) >= 0.25
%         pvalStar{i} = "*";
%     elseif pVal(i)<0.25 && pVal(i) >= 0.15
%         pvalStar{i} = "**";
%     elseif pVal(i)<0.15
%         pvalStar{i} = "***";
%     end
%     
% end

%dataSaveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll';

% save
%save([dataSaveDir filesep 'grpIndChnR_' sbjGrpStr],'grpVarDiff');

figure();
% Absolute
%subplot(1,2,1);hold on;

errorbar(timeLen./fs,mean(allFEF_CA,1),performanceLDALedoitHbTAll_CI(:,2));hold on;
errorbar(timeLen./fs,mean(squeeze(allFEF_Ind_CA(1,:,:)),1),performanceLDALedoitHbTAll_Ind_CI(1,:,2));hold on;
errorbar(timeLen./fs,mean(squeeze(allFEF_Ind_CA(2,:,:)),1),performanceLDALedoitHbTAll_Ind_CI(2,:,2));hold on;
errorbar(timeLen./fs,mean(squeeze(allFEF_Ind_CA(3,:,:)),1),performanceLDALedoitHbTAll_Ind_CI(3,:,2));hold on;

errorbar(timeLen./fs,mean(allIPS_CA,1),performanceLDALedoitHbTIPS_CI(:,2));hold on;
errorbar(timeLen./fs,mean(squeeze(allIPS_Ind_CA(1,:,:)),1),performanceLDALedoitHbTIPS_Ind_CI(1,:,2));hold on;
errorbar(timeLen./fs,mean(squeeze(allIPS_Ind_CA(2,:,:)),1),performanceLDALedoitHbTIPS_Ind_CI(2,:,2));hold on;
errorbar(timeLen./fs,mean(squeeze(allIPS_Ind_CA(3,:,:)),1),performanceLDALedoitHbTIPS_Ind_CI(3,:,2));hold on;

legend('FEF All','S4D9','S5D7','S5D9','IPS All','S10D15','S12D13','S12D15');

% Relative (Diff)
%subplot(1,2,2);hold on;
%bar([1:1:9],[grp_Diff_FEF_Hbt_Avg; grp_Diff_IPS_Hbt_Avg]);
%set(gca,'XTick',(1:1:2));
%set(gca, 'XTickLabel', {'FEF','IPS'});

if saveOp == 1
    
    if fefOnlyOp
        fn = sprintf('Grp_3FLO_FEF_ByROINames_HbT');
    else
        fn = sprintf('Grp_3FLO_AllChns_ByROINames_HbT');
    end
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
    
end

end

