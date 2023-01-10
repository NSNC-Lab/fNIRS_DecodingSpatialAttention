% split into single and multi groups.
% new. Use trials from dcNew var
% option to re-run training and classifying or pull from saved data
% sbjList is cell array of experiment "number". Ex: sbjList = {'01','02'};
% sbjList = {'08','10','12','13','14','15'};tInd = 8;featureName = 'Cumulative Sum';numRept = 5;numClasses = 2;saveOp = 0;rerunOp = 1;
% For Figure Draft 4. Now exclude rejected channel from calculating avg CA
% performance.
% Add sig test

function grpBarChartROINames_LOFO(sbjGrp,numClasses,fefOnlyOp,saveOp)

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

diff_LOFO_HbT = zeros(length(sbjList),length(chnName));
diff_LOFO_HbT_Tot = zeros(1,length(chnName));

chnList = zeros(length(sbjList),length(chnName));
chnSbjCnt = zeros(1,length(chnName));
pVal = zeros(1,length(chnName));
pvalStar = cell(1,length(chnName));

%thisT = 2*fs;
%tidx = abs(timePt-thisT)<10*eps("double");
        
for i = 1:length(sbjList)

    sbjNum = sbjList{i};

    processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];

    if numClasses == 2
        if fefOnlyOp
            fn = [processedDataDir filesep 'performance_LOFO_FEF_LR_10Hz.mat'];
        else
            fn = [processedDataDir filesep 'performance_LOFO_LR_10Hz.mat'];
        end
    elseif numClasses == 3
        fn = [processedDataDir filesep 'performance_LOFO.mat'];
    end
    
    fnStat = [processedDataDir filesep 'summaryStat_1.5.mat'];

    % only regular channels. Folds x channel x time
    load(fn,'performanceLDALedoitHbT_All','performanceLDALedoitHbT_LOFO');
    load(fnStat,'mlActAuto','behScore');
    
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
    
    allFEF_CA = mean(performanceLDALedoitHbT_All);
    
    diff_LOFO_HbT(i,:) = allFEF_CA - mean(performanceLDALedoitHbT_LOFO,1);
    
    for i2 = 1:size(chnList,2)
        if chnList(i,i2)
            temp = allFEF_CA - mean(performanceLDALedoitHbT_LOFO(:,i2),1);
            diff_LOFO_HbT_Tot(i2) = diff_LOFO_HbT_Tot(i2) + temp;
            chnSbjCnt(i2) = chnSbjCnt(i2) + 1;
        end
    end
    
end

grpVarDiff = mean(diff_LOFO_HbT,1);
    
diff_LOFO_HbT_Tot_Avg = diff_LOFO_HbT_Tot./chnSbjCnt*100;

if fefOnlyOp
    fn = 'grp_HighPerf_LOFO_FEF';
else
    fn = 'grp_HighPerf_LOFO';
end

save([dataSaveDir filesep fn],'diff_LOFO_HbT','diff_LOFO_HbT_Tot','chnSbjCnt');

for i = 1:size(chnList,2)
    
    x_temp = diff_LOFO_HbT(find(diff_LOFO_HbT(:,i)),i);
    se_temp = std(diff_LOFO_HbT(find(diff_LOFO_HbT(:,i)),i))/sqrt(chnSbjCnt(i));
    
    t_temp = (mean(x_temp))/se_temp;
    
    pVal(i) = 1-tcdf(t_temp,chnSbjCnt(i)-1);
    
    if pVal(i)<0.35 && pVal(i) >= 0.25
        pvalStar{i} = "*";
    elseif pVal(i)<0.25 && pVal(i) >= 0.15
        pvalStar{i} = "**";
    elseif pVal(i)<0.15
        pvalStar{i} = "***";
    end
    
end

%dataSaveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll';

% save
%save([dataSaveDir filesep 'grpIndChnR_' sbjGrpStr],'grpVarDiff');

[grpVarDiffSorted, origInd] = sort(diff_LOFO_HbT_Tot_Avg,'descend');

XSorted = chnName(origInd);

figure('units','normalized','outerposition',[0 0 1 1]); hold on;
hB = bar(grpVarDiffSorted);
hAx=gca; 
hT=[];
for i=1:length(hB.XData)  % iterate over number of bar objects
    if ~isempty(pvalStar{i}) && hB.YData(i)>0
        hT=[hT text(hB.XData(i) + hB.XOffset, hB.YData(i), pvalStar{i}, ...
            'VerticalAlignment','bottom','horizontalalign','center')];
    elseif ~isempty(pvalStar{i}) && hB.YData(i)<0
        hT=[hT text(hB.XData(i) + hB.XOffset, hB.YData(i)-0.1, pvalStar{i}, ...
            'VerticalAlignment','bottom','horizontalalign','center')];
    end
end
set(gca,'XTick',(1:1:length(chnName)));
set(gca, 'XTickLabel', XSorted);
xtickangle(45);
ylabel('Percentage Difference in Classification Accuracy [%]');
set(gca,'FontSize',14);

if saveOp == 1
    
    if fefOnlyOp
        fn = sprintf('Grp_LOFO_FEF_ByROINames_HbT');
    else
        fn = sprintf('Grp_LOFO_AllChns_ByROINames_HbT');
    end
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
    
end

end

