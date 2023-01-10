% split into single and multi groups.
% new. Use trials from dcNew var
% option to re-run training and classifying or pull from saved data
% sbjList is cell array of experiment "number". Ex: sbjList = {'01','02'};
% sbjList = {'08','10','12','13','14','15'};tInd = 8;featureName = 'Cumulative Sum';numRept = 5;numClasses = 2;saveOp = 0;rerunOp = 1;
% For Figure Draft 4. Now exclude rejected channel from calculating avg CA
% performance.

function grpBarChartROINames_04(sbjGrp,numClasses,rejChnOp,saveOp)

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
dataSaveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll\';

fs = 50;
timePt = (-2*fs:0.25*fs:5*fs);
zeroT = 2*fs;
timeLngth = 1*fs;
lenT = 1*fs;
numChn = 30;
timeStrtInt = 2;

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

varPerfHbO = zeros(1,length(chnName));
varPerfHbR = zeros(1,length(chnName));
varPerfHbT = zeros(1,length(chnName));

% varPerfSTDHbO = zeros(length(sbjList),length(chnName));
% varPerfSTDHbR = zeros(length(sbjList),length(chnName));
% varPerfSTDHbT = zeros(length(sbjList),length(chnName));

chnList = zeros(length(sbjList),30);
chnSbjCnt = zeros(1,30);
        
for i = 1:length(sbjList)

    sbjNum = sbjList{i};

    processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];

    if numClasses==2
        % from preprocessFNIRS06_CV_GLM_ssBeta_SingleChn_MultiOnly.m
        %fn = [processedDataDir filesep 'performance_GLM_CV_SSBeta_SingleChnHbX_LR_RejTr_SNR_1.5_OrigHomer3.mat'];
        fn = [processedDataDir filesep 'performance_GLM_CV_SSBeta_SingleChn_LR_RejTr_SNR_1.5_OrigHomer3_10Hz.mat'];
    elseif numClasses == 3
        fn = [processedDataDir filesep 'performance_GLM_CV_SSBeta_SingleChnHbX_RejTr_SNR_1.5_OrigHomer3.mat'];
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
    
    thisT = 3*fs;
    tidx = abs(timePt-thisT)<10*eps("double");

    %tempHbO = squeeze(sum(performanceLDAHbO,1))/size(performanceLDAHbO,1);
    tempHbO = sum(squeeze(performanceLDAHbO(:,:,tidx)),1)/size(performanceLDAHbO,1);
    tempHbR = sum(squeeze(performanceLDAHbR(:,:,tidx)),1)/size(performanceLDAHbR,1);
    tempHbT = sum(squeeze(performanceLDAHbT(:,:,tidx)),1)/size(performanceLDAHbT,1);
    
    for i2 = 1:size(chnList,2)
        if rejChnOp && chnList(i,i2)
            
            chnSbjCnt(i2) = chnSbjCnt(i2) + 1;
            varPerfHbO(i2) = varPerfHbO(i2)+tempHbO(i2);
            varPerfHbR(i2) = varPerfHbR(i2)+tempHbR(i2);
            varPerfHbT(i2) = varPerfHbT(i2)+tempHbT(i2);
           
        elseif ~rejChnOp
            chnSbjCnt(i2) = chnSbjCnt(i2) + 1;
            varPerfHbO(i2) = varPerfHbO(i2)+tempHbO(i2);
            varPerfHbR(i2) = varPerfHbR(i2)+tempHbR(i2);
            varPerfHbT(i2) = varPerfHbT(i2)+tempHbT(i2);
        end
    end
    
end

varPerfHbO = varPerfHbO./chnSbjCnt;
varPerfHbR = varPerfHbR./chnSbjCnt;
varPerfHbT = varPerfHbT./chnSbjCnt;

save([dataSaveDir filesep 'grpIndChnCA_' sbjGrpStr],'varPerfHbO','varPerfHbR','varPerfHbT');

% grpMeanPerfHbO = mean(nonzeros(varPerfHbO));
% grpMeanPerfHbR = mean(nonzeros(varPerfHbR));
% varPerfHbT = mean(nonzeros(varPerfHbT));

% for i = 1:length(chnName)
% %     grpMeanPerfHbO(i) = mean(nonzeros(varPerfHbO(:,i)));
% %     grpMeanPerfHbR(i) = mean(nonzeros(varPerfHbR(:,i)));
%     varPerfHbT(i) = mean(nonzeros(varPerfHbT(:,i)));
% end

diffHbT_HbO = varPerfHbT - varPerfHbO;
diffHbT_HbR = varPerfHbT - varPerfHbR;


%figure();
[varPerfHbOSorted, origIndHbO] = sort(varPerfHbO,'descend');
[varPerfHbRSorted, origIndHbR] = sort(varPerfHbR,'descend');
[varPerfHbTSorted, origIndHbT] = sort(varPerfHbT,'descend');

XSortedO = chnName(origIndHbO);
XSortedR = chnName(origIndHbR);
XSortedT = chnName(origIndHbT);

figure('units','normalized','outerposition',[0 0 1 1]); hold on;
bar(varPerfHbTSorted);
set(gca,'XTick',[1:1:30]);
set(gca, 'XTickLabel', XSortedT);
xtickangle(45);
ylabel('Accuracy [Decimal]');
set(gca,'FontSize',14);
title('HbT: 3-4s window');

if saveOp == 1
    %fn = sprintf('Sbj%s%ss%s',num2str(sbjNum),featureName,num2str(timePt(tInd)));
    fn = sprintf('Grp_CVPerfByROINames_HbT');
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
end

figure('units','normalized','outerposition',[0 0 1 1]); hold on;
bar(varPerfHbOSorted);
set(gca,'XTick',[1:1:30]);
set(gca, 'XTickLabel', XSortedO);
xtickangle(45);
ylabel('Accuracy [Decimal]');
set(gca,'FontSize',14);
title('HbO: 3-4s window');

if saveOp == 1
    %fn = sprintf('Sbj%s%ss%s',num2str(sbjNum),featureName,num2str(timePt(tInd)));
    fn = sprintf('Grp_CVPerfByROINames_HbO');
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
end

figure('units','normalized','outerposition',[0 0 1 1]); hold on;
bar(varPerfHbRSorted);
set(gca,'XTick',[1:1:30]);
set(gca, 'XTickLabel', XSortedR);
xtickangle(45);
ylabel('Accuracy [Decimal]');
set(gca,'FontSize',14);
title('HbR: 3-4s window');

% 
% 
% figure(3);hold on;
% set(gcf,'Position',get(0,'Screensize'));
% for i1 = 1:length(chnName)
%     tempY = varPerfHbT(:,i1);
%     tempX = i1*ones(1,length(sbjList));
%     for i2 = 1:length(sbjList)
%         plot(tempX(i2),tempY(i2),'o','MarkerFaceColor',colorCodes(i2,:),'MarkerEdgeColor',colorCodes(i2,:),'MarkerSize',5);
%     end
%     plot(i1,varPerfHbT(1,i1),'*','Color',colorCodes(length(sbjList)+1,:),'MarkerSize',10);
% end
% ax = gca;
% ax.XTick = 1:1:length(chnName);
% ax.XTickLabels = chnName;
% xtickangle(45);
% title(sprintf('HbT: Sum. %ss interval. At %ss',num2str(timeLngth/fs),num2str(timeStrtInt)));
% hold off;
% 
if saveOp == 1
    %fn = sprintf('Sbj%s%ss%s',num2str(sbjNum),featureName,num2str(timePt(tInd)));
    fn = sprintf('Grp_CVPerfByROINames_HbR');
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
end

figure('units','normalized','outerposition',[0 0 1 1]); hold on;
bar(diffHbT_HbO);
set(gca,'XTick',[1:1:30]);
set(gca, 'XTickLabel', chnName);
xtickangle(45);
ylabel('Accuracy [Decimal]');
set(gca,'FontSize',14);
title('HbT-HbO: 3-4s window');

figure('units','normalized','outerposition',[0 0 1 1]); hold on;
bar(diffHbT_HbR);
set(gca,'XTick',[1:1:30]);
set(gca, 'XTickLabel', chnName);
xtickangle(45);
ylabel('Accuracy [Decimal]');
set(gca,'FontSize',14);
title('HbT-HbR: 3-4s window');


end

