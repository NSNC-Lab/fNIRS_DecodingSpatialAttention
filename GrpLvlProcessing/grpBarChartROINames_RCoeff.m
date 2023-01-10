% split into single and multi groups.
% new. Use trials from dcNew var
% option to re-run training and classifying or pull from saved data
% sbjList is cell array of experiment "number". Ex: sbjList = {'01','02'};
% sbjList = {'08','10','12','13','14','15'};tInd = 8;featureName = 'Cumulative Sum';numRept = 5;numClasses = 2;saveOp = 0;rerunOp = 1;
% For Figure Draft 4. Now exclude rejected channel from calculating avg CA
% performance.

function grpBarChartROINames_RCoeff(sbjGrp,numClasses,rejChnOp,opPValFilt,saveOp)

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

varPerfHbO = zeros(length(sbjList),length(chnName));
varPerfHbR = zeros(length(sbjList),length(chnName));
varPerfHbT = zeros(length(sbjList),length(chnName));

behArr = zeros(1,length(sbjList));

% varPerfSTDHbO = zeros(length(sbjList),length(chnName));
% varPerfSTDHbR = zeros(length(sbjList),length(chnName));
% varPerfSTDHbT = zeros(length(sbjList),length(chnName));

chnList = zeros(length(sbjList),30);
chnSbjCnt = zeros(1,30);
rValHbT = zeros(1,30);
pValHbT = zeros(1,30);

thisT = 3*fs;
tidx = abs(timePt-thisT)<10*eps("double");
        
for i = 1:length(sbjList)

    sbjNum = sbjList{i};

    processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];

    if numClasses==2
        %fn = [processedDataDir filesep 'performance_GLM_CV_SSBeta_SingleChn_LR_RejTr_SNR_1.5_OrigHomer3.mat'];
        fn = [processedDataDir filesep 'performance_GLM_CV_SSBeta_SingleChn_LR_RejTr_SNR_1.5_OrigHomer3_10Hz.mat'];
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
    load(fnStat,'mlActAuto','behScore');
    
    behArr(i) = behScore;

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

%varPerfHbO = varPerfHbO./chnSbjCnt;
%varPerfHbR = varPerfHbR./chnSbjCnt;
%varPerfHbT = varPerfHbT./chnSbjCnt;

for i = 1:size(chnList,2)

    %disp(i);
    tempIdx = varPerfHbT(:,i)~=0;
    temp = varPerfHbT(tempIdx,i);
    tempBeh = behArr(tempIdx);
    
    % replace this with p value using student's t distribution
    [tempR,tempP] = corrcoef(temp,tempBeh');
    %tempP = 
    
    if size(tempR,1)==1 && size(tempR,2)==1
        rValHbT(i) = 2;
        pValHbT(i) = 2;
    else
        rValHbT(i) = tempR(1,2);
        pValHbT(i) = tempP(1,2);
    end

end

dataSaveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll';

% save
save([dataSaveDir filesep 'grpIndChnR_' sbjGrpStr],'rValHbT');

% grpMeanPerfHbO = mean(nonzeros(varPerfHbO));
% grpMeanPerfHbR = mean(nonzeros(varPerfHbR));
% varPerfHbT = mean(nonzeros(varPerfHbT));

% for i = 1:length(chnName)
% %     grpMeanPerfHbO(i) = mean(nonzeros(varPerfHbO(:,i)));
% %     grpMeanPerfHbR(i) = mean(nonzeros(varPerfHbR(:,i)));
%     varPerfHbT(i) = mean(nonzeros(varPerfHbT(:,i)));
% end

if opPValFilt
    pCritIdx = pValHbT<0.4;
%pValHbT = pValHbT(pCritIdx);
    rValHbT = rValHbT(pCritIdx);
    chnName = chnName(pCritIdx);
    chnSbjCnt = chnSbjCnt(pCritIdx);
end

[rValHbT, origInd] = sort(rValHbT,'descend');

XSorted = chnName(origInd);
chnSbjCnt = chnSbjCnt(origInd);

figure('units','normalized','outerposition',[0 0 1 1]); hold on;
hB = bar(rValHbT);
hAx=gca; 
hT=[];
for i=1:length(hB.XData)  % iterate over number of bar objects
    hT=[hT text(hB.XData(i) + hB.XOffset, hB.YData(i), num2str(chnSbjCnt(i)), ...
        'VerticalAlignment','bottom','horizontalalign','center')];
end
set(gca,'XTick',(1:1:30));
set(gca, 'XTickLabel', XSorted);
xtickangle(45);
ylabel('R Correlation Coefficient');
set(gca,'FontSize',14);

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
    
    if opPValFilt
        fn = sprintf('Grp_RValByROINames_PValFilt_HbT');
    else
        fn = sprintf('Grp_RValByROINames_HbT');
    end
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
    
end

% tempIdx = varPerfHbT(:,17)~=0;
% temp = varPerfHbT(tempIdx,17);
% tempBeh = behArr(tempIdx);
% 
% figure();
% mdl = fitlm(tempBeh',temp);
% x = anova(mdl);
% h = plot(mdl);


end

