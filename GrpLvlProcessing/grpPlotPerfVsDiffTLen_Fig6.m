% For figure 6 in Draft 3
% 2 subplots (1,2). HbT. Show only high performing group.
% Only show LDA classifier.
% 2-class classification at left panel.
% 3-class classification at right panel.
% look at 7-color scheme for readability.
% 08/01/2022, add CI options
%   1 for counting each fold as one sample point
%   2 for counting each sbj as one sample point

function grpPlotPerfVsDiffTLen_Fig6(goodSbj,numClasses,saveOp,itrOp,ciOp)

    if goodSbj == 1
        %sbjList = {'08','10','12','16','17','18','19'};
        %sbjList = {'08','12','13','16','17','19'};
        % eligible if less than 20 channels rejected
        sbjList = {'08','12','13','16','19','21','24'};
        %sbjList = {'08','12','13','16','19','21','24','18','20'};
        dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
            filesep 'sbjListPerf_DiffTLen_Good.mat'];
    elseif goodSbj == 0
        sbjList = {'14','15','22','23','25'};
        dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
            filesep 'sbjListPerf_DiffTLen_Bad.mat'];
    else
        %sbjList = {'08','10','12','13','14','15','16','17','18','19'};
        sbjList = {'08','12','13','14','15','16','19','21','22','23','24','25'};
        %sbjList = {'08','12','13','14','15','16','17','18','19','20','21','22','23','24','25'};
        dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
            filesep 'sbjListPerf_DiffTLen_All.mat'];
    end

    figSaveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll\Figures';

    fs = 50;
    %timePt = (0:0.25*fs:5*fs)+2*fs;
    timeLen = [0.1 0.2 0.5 1 1.5 2 3 4 5].*fs;
    zeroT = 2*fs;
    numSubsets = 7;
    
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

    sbjListPerfHbO = zeros(numSubsets,5*10*length(sbjList),length(timeLen));
    sbjListPerfHbR = zeros(numSubsets,5*10*length(sbjList),length(timeLen));
    sbjFoldListPerfHbT = zeros(numSubsets,5*10*length(sbjList),length(timeLen));
    
    sbjListPerfHbT = zeros(numSubsets,length(sbjList),length(timeLen));
    

%     sbjListPerfSTDHbO = zeros(numClassifiers,length(sbjList),length(timePt));
%     sbjListPerfSTDHbR = zeros(numClassifiers,length(sbjList),length(timePt));
%     sbjListPerfSTDHbT = zeros(numClassifiers,length(sbjList),length(timePt));

    for i = 1:length(sbjList)
        sbjNum = sbjList{i};

        processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];

        if numClasses == 2
            %savePerfFN = 'performanceLinearDiscriminantUpdated_DiffTLen_StartMovie_LR.mat';
            %savePerfFN = 'PerformanceCumsumVsTime_DiffTLen_StartMovie_LR_AllChns.mat';
            %savePerfFN = 'performance_GLM_CV_SSBeta_LR_RejTr_SNR_1.5.mat';
            %savePerfFN = 'performance_GLM_CV_SSBeta_LR_RejTr.mat';
            %savePerfFN = 'performance_GLM_CV_SSBeta_LR_DiffTLen_RejTr';
            savePerfFN = 'performance_GLM_CV_SSBeta_LR_DiffTLen_RejTr_10Hz';
            %savePerfFN = 'performance_GLM_CV_SSBeta_LR_RejTr_SNR_1.5_OrigHomer3.mat';
        else
            %savePerfFN = 'performanceLinearDiscriminantUpdated_DiffTLen_StartMovie.mat';
            %savePerfFN = 'PerformanceCumsumVsTime_DiffTLen_StartMovie_AllChns.mat';
            savePerfFN = 'performance_GLM_CV_SSBeta_RejTr_DiffTLen_.mat';
        end
%         load([processedDataDir filesep savePerfFN],'performanceLDALedoitHbOFold',...
%             'performanceLDALedoitHbRFold','performanceLDALedoitHbTFold',...
%             'performanceLDACERNNHbOFold','performanceLDACERNNHbRFold',...
%             'performanceLDACERNNHbTFold',...
%             'performanceLogRegHbOFold','performanceLogRegHbRFold',...
%             'performanceLogRegHbTFold',...
%             'performanceSVMHbOFold','performanceSVMHbRFold',...
%             'performanceSVMHbTFold',...
%             'performanceBaggingHbOFold','performanceBaggingHbRFold',...
%             'performanceBaggingHbTFold',...
%             'performanceBoostedHbOFold','performanceBoostedHbRFold',...
% %             'performanceBoostedHbTFold');
%         load([processedDataDir filesep savePerfFN],'performanceLDALedoitHbOFold',...
%             'performanceLDALedoitHbRFold','performanceLDALedoitHbTFold');
        load([processedDataDir filesep savePerfFN],...
            'performanceLDALedoitHbT_All',...
            'performanceLDALedoitHbT_LFEF',...
            'performanceLDALedoitHbT_RFEF',...
            'performanceLDALedoitHbT_LSTG',...
            'performanceLDALedoitHbT_RSTG',...
            'performanceLDALedoitHbT_LIPS',...
            'performanceLDALedoitHbT_RIPS');

    %     performanceLDALedoitHbOCI = zeros(length(timePt),2);
    %     performanceLDALedoitHbRCI = zeros(length(timePt),2);
    %     performanceLDALedoitHbTCI = zeros(length(timePt),2);
    % 
    %     performanceLDACERNNHbOCI = zeros(length(timePt),2);
    %     performanceLDACERNNHbRCI = zeros(length(timePt),2);
    %     performanceLDACERNNHbTCI = zeros(length(timePt),2);
    % 
    %     performanceLogRegHbOCI = zeros(length(timePt),2);
    %     performanceLogRegHbRCI = zeros(length(timePt),2);
    %     performanceLogRegHbTCI = zeros(length(timePt),2);
    % 
    %     performanceSVMHbOCI = zeros(length(timePt),2);
    %     performanceSVMHbRCI = zeros(length(timePt),2);
    %     performanceSVMHbTCI = zeros(length(timePt),2);
    % 
    %     performanceBaggingHbOCI = zeros(length(timePt),2);
    %     performanceBaggingHbRCI = zeros(length(timePt),2);
    %     performanceBaggingHbTCI = zeros(length(timePt),2);
    % 
    %     performanceBoostedHbOCI = zeros(length(timePt),2);
    %     performanceBoostedHbRCI = zeros(length(timePt),2);
    %     performanceBoostedHbTCI = zeros(length(timePt),2);


%         performanceLDALedoitHbOFoldRS = reshape(performanceLDALedoitHbOFold,...
%             [size(performanceLDALedoitHbOFold,1), size(performanceLDALedoitHbOFold,2)*size(performanceLDALedoitHbOFold,3)]);
%         performanceLDALedoitHbRFoldRS = reshape(performanceLDALedoitHbRFold,...
%             [size(performanceLDALedoitHbRFold,1), size(performanceLDALedoitHbRFold,2)*size(performanceLDALedoitHbRFold,3)]);
%         performanceLDALedoitHbTFoldRS = reshape(performanceLDALedoitHbTFold,...
%             [size(performanceLDALedoitHbTFold,1), size(performanceLDALedoitHbTFold,2)*size(performanceLDALedoitHbTFold,3)]);
% 
%         % LFEF
%         performance_LFEF_HbOFoldRS = reshape(performance_LFEF_HbOFold,...
%             [size(performance_LFEF_HbOFold,1), size(performance_LFEF_HbOFold,2)*size(performance_LFEF_HbOFold,3)]);
%         performance_LFEF_HbRFoldRS = reshape(performance_LFEF_HbRFold,...
%             [size(performance_LFEF_HbRFold,1), size(performance_LFEF_HbRFold,2)*size(performance_LFEF_HbRFold,3)]);
%         performance_LFEF_HbTFoldRS = reshape(performance_LFEF_HbTFold,...
%             [size(performance_LFEF_HbTFold,1), size(performance_LFEF_HbTFold,2)*size(performance_LFEF_HbTFold,3)]);
% 
%         % RFEF
%         performance_RFEF_HbOFoldRS = reshape(performance_RFEF_HbOFold,...
%             [size(performance_RFEF_HbOFold,1), size(performance_RFEF_HbOFold,2)*size(performance_RFEF_HbOFold,3)]);
%         performance_RFEF_HbRFoldRS = reshape(performance_RFEF_HbRFold,...
%             [size(performance_RFEF_HbRFold,1), size(performance_RFEF_HbRFold,2)*size(performance_RFEF_HbRFold,3)]);
%         performance_RFEF_HbTFoldRS = reshape(performance_RFEF_HbTFold,...
%             [size(performance_RFEF_HbTFold,1), size(performance_RFEF_HbTFold,2)*size(performance_RFEF_HbTFold,3)]);
%         
%         % LSTG
%         performance_LSTG_HbOFoldRS = reshape(performance_LSTG_HbOFold,...
%             [size(performance_LSTG_HbOFold,1), size(performance_LSTG_HbOFold,2)*size(performance_LSTG_HbOFold,3)]);
%         performance_LSTG_HbRFoldRS = reshape(performance_LSTG_HbRFold,...
%             [size(performance_LSTG_HbRFold,1), size(performance_LSTG_HbRFold,2)*size(performance_LSTG_HbRFold,3)]);
%         performance_LSTG_HbTFoldRS = reshape(performance_LSTG_HbTFold,...
%             [size(performance_LSTG_HbTFold,1), size(performance_LSTG_HbTFold,2)*size(performance_LSTG_HbTFold,3)]);
%         
%         % RSTG
%         performance_RSTG_HbOFoldRS = reshape(performance_RSTG_HbOFold,...
%             [size(performance_RSTG_HbOFold,1), size(performance_RSTG_HbOFold,2)*size(performance_RSTG_HbOFold,3)]);
%         performance_RSTG_HbRFoldRS = reshape(performance_RSTG_HbRFold,...
%             [size(performance_RSTG_HbRFold,1), size(performance_RSTG_HbRFold,2)*size(performance_RSTG_HbRFold,3)]);
%         performance_RSTG_HbTFoldRS = reshape(performance_RSTG_HbTFold,...
%             [size(performance_RSTG_HbTFold,1), size(performance_RSTG_HbTFold,2)*size(performance_RSTG_HbTFold,3)]);
%         
%         % LIPS
%         performance_LIPS_HbOFoldRS = reshape(performance_LIPS_HbOFold,...
%             [size(performance_LIPS_HbOFold,1), size(performance_LIPS_HbOFold,2)*size(performance_LIPS_HbOFold,3)]);
%         performance_LIPS_HbRFoldRS = reshape(performance_LIPS_HbRFold,...
%             [size(performance_LIPS_HbRFold,1), size(performance_LIPS_HbRFold,2)*size(performance_LIPS_HbRFold,3)]);
%         performance_LIPS_HbTFoldRS = reshape(performance_LIPS_HbTFold,...
%             [size(performance_LIPS_HbTFold,1), size(performance_LIPS_HbTFold,2)*size(performance_LIPS_HbTFold,3)]);
%         
%         % RIPS
%         performance_RIPS_HbOFoldRS = reshape(performance_RIPS_HbOFold,...
%             [size(performance_RIPS_HbOFold,1), size(performance_RIPS_HbOFold,2)*size(performance_RIPS_HbOFold,3)]);
%         performance_RIPS_HbRFoldRS = reshape(performance_RIPS_HbRFold,...
%             [size(performance_RIPS_HbRFold,1), size(performance_RIPS_HbRFold,2)*size(performance_RIPS_HbRFold,3)]);
%         performance_RIPS_HbTFoldRS = reshape(performance_RIPS_HbTFold,...
%             [size(performance_RIPS_HbTFold,1), size(performance_RIPS_HbTFold,2)*size(performance_RIPS_HbTFold,3)]);
        
        % sanity check. Passed!
    %     temp1 = performanceLDALedoitHbOFold(1,:,:);
    %     test1 = mean(temp1(:));
    %     temp2 = performanceLDALedoitHbOFold(2,:,:);
    %     test2 = mean(temp2(:));
    %     temp3 = performanceLDALedoitHbOFold(3,:,:);
    %     test3 = mean(temp3(:));
        % performanceLDALedoitHbO_All = zeros(nrept*kFold,length(timeLen));
        
        startIdx = (i-1)*50+1;
        endIdx = i*50;
        
        %sbjListPerfHbO(1,startIdx:endIdx,:) = mean(performanceLDALedoitHbO_All,1);
        %sbjListPerfHbR(1,startIdx:endIdx,:) = mean(performanceLDALedoitHbR_All,1);
        sbjFoldListPerfHbT(1,startIdx:endIdx,:) = performanceLDALedoitHbT_All;
        
        %sbjListPerfHbO(2,i,:) = mean(performanceLDALedoitHbO_LFEF,1);
        %sbjListPerfHbR(2,i,:) = mean(performanceLDALedoitHbR_LFEF,1);
        sbjFoldListPerfHbT(2,startIdx:endIdx,:) = performanceLDALedoitHbT_LFEF;
        
        %sbjListPerfHbO(3,i,:) = mean(performanceLDALedoitHbO_RFEF,1);
        %sbjListPerfHbR(3,i,:) = mean(performanceLDALedoitHbR_RFEF,1);
        sbjFoldListPerfHbT(3,startIdx:endIdx,:) = performanceLDALedoitHbT_RFEF;
        
        %sbjListPerfHbO(4,i,:) = mean(performanceLDALedoitHbO_LSTG,1);
        %sbjListPerfHbR(4,i,:) = mean(performanceLDALedoitHbR_LSTG,1);
        sbjFoldListPerfHbT(4,startIdx:endIdx,:) = performanceLDALedoitHbT_LSTG;
        
        %sbjListPerfHbO(5,i,:) = mean(performanceLDALedoitHbO_RSTG,1);
        %sbjListPerfHbR(5,i,:) = mean(performanceLDALedoitHbR_RSTG,1);
        sbjFoldListPerfHbT(5,startIdx:endIdx,:) = performanceLDALedoitHbT_RSTG;
        
        %sbjListPerfHbO(6,i,:) = mean(performanceLDALedoitHbO_LIPS,1);
        %sbjListPerfHbR(6,i,:) = mean(performanceLDALedoitHbR_LIPS,1);
        sbjFoldListPerfHbT(6,startIdx:endIdx,:) = performanceLDALedoitHbT_LIPS;
        
        %sbjListPerfHbO(7,i,:) = mean(performanceLDALedoitHbO_RIPS,1);
        %sbjListPerfHbR(7,i,:) = mean(performanceLDALedoitHbR_RIPS,1);
        sbjFoldListPerfHbT(7,startIdx:endIdx,:) = performanceLDALedoitHbT_RIPS;
        
        sbjListPerfHbT(1,i,:) = mean(performanceLDALedoitHbT_All,1);
        sbjListPerfHbT(2,i,:) = mean(performanceLDALedoitHbT_LFEF,1);
        sbjListPerfHbT(3,i,:) = mean(performanceLDALedoitHbT_RFEF,1);
        sbjListPerfHbT(4,i,:) = mean(performanceLDALedoitHbT_LSTG,1);
        sbjListPerfHbT(5,i,:) = mean(performanceLDALedoitHbT_RSTG,1);
        sbjListPerfHbT(6,i,:) = mean(performanceLDALedoitHbT_LIPS,1);
        sbjListPerfHbT(7,i,:) = mean(performanceLDALedoitHbT_RIPS,1);

    end
    
    save(dataFN,'sbjListPerfHbO','sbjListPerfHbR','sbjListPerfHbT');

    %grpPerfHbO = sum(sbjListPerfHbO,2) ./ sum(sbjListPerfHbO~=0,2);
    %grpPerfHbR = sum(sbjListPerfHbR,2) ./ sum(sbjListPerfHbR~=0,2);
    grpPerfHbT = sum(sbjFoldListPerfHbT,2) ./ sum(sbjFoldListPerfHbT~=0,2);
    
    if ciOp == 1
        grpPerfHbTSE = squeeze(std(sbjFoldListPerfHbT,[],2))./sqrt(sum(sbjFoldListPerfHbT~=0,2));
    elseif ciOp == 2
        grpPerfHbTSE = squeeze(std(sbjListPerfHbT,[],2))./squeeze(sqrt(sum(sbjListPerfHbT~=0,2)));
    end
    
    conf = 0.95;

    nSz = size(grpPerfHbTSE,1);
    grpPerfCIHbT = zeros(2,numSubsets,length(timeLen));
    
    for i = 1:numSubsets
        for idxTime = 1:length(timeLen)

            grpPerfCIHbT(:,i,idxTime) = calcCITDist(nSz,grpPerfHbTSE(i,idxTime),conf);

        end
    end

    cmap = loadDefaultColors(4);
    cmapROI = [0 0 0; cmap(1,:); cmap(1,:); cmap(3,:); cmap(3,:); cmap(2,:); cmap(2,:)];
    lineStyleROI = {'-','-','--','-','--','-','--'};
    yLimAxis = [0.3 1];
    
    [h,p,ci,stats] = ttest(sbjListPerfHbT(1,:,4),0.5);
    disp(p);
    disp(stats);
    
    figure();hold on;
    subsetIdx = [2 3 6 7];
    %for i = 1:numSubsets-2
    %for i = 1:numSubsets-4
    errorbar(squeeze(timeLen./fs),squeeze(grpPerfHbT(1,1,:)),abs(squeeze(grpPerfCIHbT(1,1,:))),abs(squeeze(grpPerfCIHbT(2,1,:))),...
        'Color',cmapROI(1,:),'LineStyle',lineStyleROI{1});
    for i = 1:length(subsetIdx)
        plot(timeLen./fs,squeeze(grpPerfHbT(subsetIdx(i),1,:)),...
            'Color',cmapROI(subsetIdx(i),:),'LineStyle',lineStyleROI{subsetIdx(i)});
    end

    if itrOp
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
    else
        ylim([0.3 1]);
        ylabel('Accuracy');
    end
    if numClasses == 2
        %yline(.5,'--');
    else
        %yline(0.333,'--');
    end
    grid;
    xticks([0.2 0.5 1 1.5 2 3 4 5]);
    %title(sprintf('CA vs Probe Subset: Δ[HbT]'));
    %legend({'LDAShrink','LDA CERNN','LogReg','SVM','Bagging','Boosted'},'Location','southeast');
    %legend({'All','LFEF','RFEF','LSTG','RSTG','LIPS','RIPS'},'Location','Southeast');
    %legend({'All','LFEF','RFEF','LSTG','RSTG'},'Location','Southeast');
    %legend({'All','LFEF','RFEF'},'Location','Southeast');
    legend({'All (30)','LFEF (10)','RFEF (10)','LIPS (4)','RIPS (4)'},'Location','Northwest');
    xlabel('Window Length [s]');
    hold off;
    
    %figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    figure();subplot(1,3,1);hold on;
    %subplot(1,2,1);hold on;
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

%     errorbar(timePt./fs-2,performanceLDALedoitHbOMean,performanceLDALedoitHbOCI(:,2),'Color',cmap(1,:));
%     errorbar(timePt./fs-2,performanceLDACERNNHbOMean,performanceLDACERNNHbOCI(:,2),'Color',cmap(2,:));
%     % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
%     % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
%     errorbar(timePt./fs-2,performanceLogRegHbOMean,performanceLogRegHbOCI(:,2),'Color',cmap(3,:));
%     errorbar(timePt./fs-2,performanceSVMHbOMean,performanceSVMHbOCI(:,2),'Color',cmap(4,:));
%     errorbar(timePt./fs-2,performanceBaggingHbOMean,performanceBaggingHbOCI(:,2),'Color',cmap(5,:));
%     errorbar(timePt./fs-2,performanceBoostedHbOMean,performanceBoostedHbOCI(:,2),'Color',cmap(6,:));

    subsetIdx = [1 2 3 6 7];
    %for i = 1:numSubsets-2
    %for i = 1:numSubsets-4
%     for i = 1:length(subsetIdx)
%         errorbar(timeLen./fs,grpPerfHbT(subsetIdx(i),:),abs(grpPerfCIHbT(1,:)),abs(grpPerfCIHbT(2,:)),...
%             'Color',cmapROI(subsetIdx(i),:),'LineStyle',lineStyleROI{subsetIdx(i)});
%     end
    for i = 1:length(subsetIdx)
        plot(timeLen./fs,squeeze(grpPerfHbT(subsetIdx(i),1,:)),...
            'Color',cmapROI(subsetIdx(i),:),'LineStyle',lineStyleROI{subsetIdx(i)});
    end

    if itrOp
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
    else
        ylim([0.3 1]);
        ylabel('Accuracy');
    end
    if numClasses == 2
        %yline(.5,'--');
    else
        %yline(0.333,'--');
    end
    grid;
    xticks([0.2 0.5 1 1.5 2 3 4 5]);
    %title(sprintf('CA vs Probe Subset: Δ[HbT]'));
    %legend({'LDAShrink','LDA CERNN','LogReg','SVM','Bagging','Boosted'},'Location','southeast');
    %legend({'All','LFEF','RFEF','LSTG','RSTG','LIPS','RIPS'},'Location','Southeast');
    %legend({'All','LFEF','RFEF','LSTG','RSTG'},'Location','Southeast');
    %legend({'All','LFEF','RFEF'},'Location','Southeast');
    legend({'All','LFEF','RFEF','LIPS','RIPS'},'Location','Southeast');
    xlabel('Window Length [s]');
    hold off;

    subplot(1,3,2);hold on;
    %subplot(1,2,1);hold on;
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

%     errorbar(timePt./fs-2,performanceLDALedoitHbOMean,performanceLDALedoitHbOCI(:,2),'Color',cmap(1,:));
%     errorbar(timePt./fs-2,performanceLDACERNNHbOMean,performanceLDACERNNHbOCI(:,2),'Color',cmap(2,:));
%     % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
%     % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
%     errorbar(timePt./fs-2,performanceLogRegHbOMean,performanceLogRegHbOCI(:,2),'Color',cmap(3,:));
%     errorbar(timePt./fs-2,performanceSVMHbOMean,performanceSVMHbOCI(:,2),'Color',cmap(4,:));
%     errorbar(timePt./fs-2,performanceBaggingHbOMean,performanceBaggingHbOCI(:,2),'Color',cmap(5,:));
%     errorbar(timePt./fs-2,performanceBoostedHbOMean,performanceBoostedHbOCI(:,2),'Color',cmap(6,:));

    %for i = 1:numSubsets-2
    %for i = 1:numSubsets-4
    subsetIdx = [1 2 6];
    for i = 1:length(subsetIdx)
        errorbar(timeLen./fs,squeeze(grpPerfHbT(subsetIdx(i),1,:)),abs(squeeze(grpPerfCIHbT(1,subsetIdx(i),:))),abs(squeeze(grpPerfCIHbT(2,subsetIdx(i),:))),...
            'Color',cmapROI(subsetIdx(i),:),'LineStyle',lineStyleROI{subsetIdx(i)});
    end
    %b = bar(timeLen./fs,grpPerfHbT(subsetIdx,:),4,'FaceColor','flat');
    %b = bar(1:9,grpPerfHbT(subsetIdx,:),'FaceColor','flat');
        
%     for i = 1:length(subsetIdx)
%         b(i).CData = cmapROI(subsetIdx(i),:);
%     end

    if itrOp
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
    else
        ylim([0.3 1]);
        ylabel('Accuracy');
    end
    if numClasses == 2
        %yline(.5,'--');
    else
        %yline(0.333,'--');
    end
    grid;
    xticks([0.2 0.5 1 1.5 2 3 4 5]);
    %title(sprintf('CA vs Probe Subset: Δ[HbT]'));
    %legend({'LDAShrink','LDA CERNN','LogReg','SVM','Bagging','Boosted'},'Location','southeast');
    %legend({'All','LFEF','RFEF','LSTG','RSTG','LIPS','RIPS'},'Location','Southeast');
    %legend({'All','LFEF','RFEF','LSTG','RSTG'},'Location','Southeast');
    %legend({'All','LFEF','RFEF'},'Location','Southeast');
    %legend({'All','LFEF','RFEF','LIPS','RIPS'},'Location','Southeast');
    xlabel('Window Length [s]');
    hold off;
    
    subplot(1,3,3);hold on;
    %subplot(1,2,1);hold on;
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

%     errorbar(timePt./fs-2,performanceLDALedoitHbOMean,performanceLDALedoitHbOCI(:,2),'Color',cmap(1,:));
%     errorbar(timePt./fs-2,performanceLDACERNNHbOMean,performanceLDACERNNHbOCI(:,2),'Color',cmap(2,:));
%     % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
%     % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
%     errorbar(timePt./fs-2,performanceLogRegHbOMean,performanceLogRegHbOCI(:,2),'Color',cmap(3,:));
%     errorbar(timePt./fs-2,performanceSVMHbOMean,performanceSVMHbOCI(:,2),'Color',cmap(4,:));
%     errorbar(timePt./fs-2,performanceBaggingHbOMean,performanceBaggingHbOCI(:,2),'Color',cmap(5,:));
%     errorbar(timePt./fs-2,performanceBoostedHbOMean,performanceBoostedHbOCI(:,2),'Color',cmap(6,:));

    %for i = 1:numSubsets-2
    %for i = 1:numSubsets-4
    subsetIdx = [1 3 7];
    for i = 1:length(subsetIdx)
        errorbar(timeLen./fs,squeeze(grpPerfHbT(subsetIdx(i),1,:)),abs(squeeze(grpPerfCIHbT(1,subsetIdx(i),:))),abs(squeeze(grpPerfCIHbT(2,subsetIdx(i),:))),...
            'Color',cmapROI(subsetIdx(i),:),'LineStyle',lineStyleROI{subsetIdx(i)});
    end
    %b = bar(timeLen./fs,grpPerfHbT(subsetIdx,:),4,'FaceColor','flat');
    %b = bar(1:9,grpPerfHbT(subsetIdx,:),'FaceColor','flat');
        
%     for i = 1:length(subsetIdx)
%         b(i).CData = cmapROI(subsetIdx(i),:);
%     end

    if itrOp
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
    else
        ylim([0.3 1]);
        ylabel('Accuracy');
    end
    if numClasses == 2
        %yline(.5,'--');
    else
        %yline(0.333,'--');
    end
    grid;
    xticks([0.2 0.5 1 1.5 2 3 4 5]);
    %title(sprintf('CA vs Probe Subset: Δ[HbT]'));
    %legend({'LDAShrink','LDA CERNN','LogReg','SVM','Bagging','Boosted'},'Location','southeast');
    %legend({'All','LFEF','RFEF','LSTG','RSTG','LIPS','RIPS'},'Location','Southeast');
    %legend({'All','LFEF','RFEF','LSTG','RSTG'},'Location','Southeast');
    %legend({'All','LFEF','RFEF'},'Location','Southeast');
    %legend({'All','LFEF','RFEF','LIPS','RIPS'},'Location','Southeast');
    xlabel('Window Length [s]');
    hold off;
    
    %subplot(1,3,2);hold on;
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

%     errorbar(timePt./fs-2,performanceLDALedoitHbRMean,performanceLDALedoitHbRCI(:,2),'Color',cmap(1,:));
%     errorbar(timePt./fs-2,performanceLDACERNNHbRMean,performanceLDACERNNHbRCI(:,2),'Color',cmap(2,:));
%     % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
%     % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
%     errorbar(timePt./fs-2,performanceLogRegHbRMean,performanceLogRegHbRCI(:,2),'Color',cmap(3,:));
%     errorbar(timePt./fs-2,performanceSVMHbRMean,performanceSVMHbRCI(:,2),'Color',cmap(4,:));
%     errorbar(timePt./fs-2,performanceBaggingHbRMean,performanceBaggingHbRCI(:,2),'Color',cmap(5,:));
%     errorbar(timePt./fs-2,performanceBoostedHbRMean,performanceBoostedHbRCI(:,2),'Color',cmap(6,:));

%     for i = 1:numSubsets
%         errorbar(timeLen./fs,grpPerfHbR(i,:),grpPerfSTDHbR(i,:),'Color',cmap(i,:));
%     end
% 
%     if itrOp
%         yyaxis left;
%         ylim([0 1]);
%         ylabel('Accuracy');
%         yyaxis right;
%         if numClasses == 2
%             ylim([-max(itf) max(itf)]);
%             yticks([-max(itf) 0 max(itf)]);
%             yticklabels({num2str(max(itf)) num2str(0) num2str(max(itf))});
%         else
%             ylim([min(itf) max(itf)]);
%             yticks([min(itf) 0 max(itf)]);
%             yticks(linspace(min(itf), max(itf),6));
%             tempIdx = linspace(1,100,6);
%             yticklabels([-itf(2) -itf(round(tempIdx(2))) ...
%                 itf(round(tempIdx(3))) itf(round(tempIdx(4))) itf(round(tempIdx(5))) ...
%                 itf(round(tempIdx(6)))]);
%             %yticklabels({num2str(max(itf)) num2str(0) num2str(max(itf))});  
%         end
%         ylabel('ITR bits/trial');
%     else
%         ylim([0 1]);
%         ylabel('Accuracy');
%     end
%     if numClasses == 2
%         yline(.5,'--');
%     else
%         yline(0.333,'--');
%     end
%     %xlim([0 7]);
%     title(sprintf('Group Avg: Δ[HbR] Multi'));
%     %legend('S1D1','S1D2','S1D3','S1D4');
%     xlabel('Window Length [s]');
%     hold off;
% 
%     subplot(1,3,3);hold on;
%     % for j = 1:size(performanceArrMultiHbT,1)
%     %     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
%     %     plot(timePt./fs-2,performanceArrMultiHbT(j,:),'Color',cmap(j,:));
%     % end
%     % plot(timePt./fs-2,performanceArrMultiHbT_95,'Color',cmap(1,:));
%     % plot(timePt./fs-2,performanceArrMultiHbT_85,'Color',cmap(2,:));
%     % plot(timePt./fs-2,performanceArrMultiHbT_75,'Color',cmap(3,:));
%     % plot(timePt./fs-2,performanceArrMultiHbT_65,'Color',cmap(4,:));
%     % plot(timePt./fs-2,performanceArrMultiHbT_55,'Color',cmap(5,:));
%     % plot(timePt./fs-2,performanceLDALedoitHbT,'Color',cmap(6,:));
%     % plot(timePt./fs-2,performanceLDACERNNHbT,'Color',cmap(7,:));
%     % % plot(timePt./fs-2,performanceQDALedoitHbT,'Color',cmap(7,:));
%     % % plot(timePt./fs-2,performanceLDAGDHbT,'Color',cmap(8,:));
%     % plot(timePt./fs-2,performanceLogRegHbT,'Color',cmap(8,:));
%     % plot(timePt./fs-2,performanceSVMHbT,'Color',cmap(9,:));
% 
% %     errorbar(timePt./fs-2,performanceLDALedoitHbTMean,performanceLDALedoitHbTCI(:,2),'Color',cmap(1,:));
% %     errorbar(timePt./fs-2,performanceLDACERNNHbTMean,performanceLDACERNNHbTCI(:,2),'Color',cmap(2,:));
% %     % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
% %     % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
% %     errorbar(timePt./fs-2,performanceLogRegHbTMean,performanceLogRegHbTCI(:,2),'Color',cmap(3,:));
% %     errorbar(timePt./fs-2,performanceSVMHbTMean,performanceSVMHbTCI(:,2),'Color',cmap(4,:));
% %     errorbar(timePt./fs-2,performanceBaggingHbTMean,performanceBaggingHbTCI(:,2),'Color',cmap(5,:));
% %     errorbar(timePt./fs-2,performanceBoostedHbTMean,performanceBoostedHbTCI(:,2),'Color',cmap(6,:));
% 
%     for i = 1:numSubsets
%         errorbar(timeLen./fs,grpPerfHbT(i,:),grpPerfSTDHbT(i,:),'Color',cmap(i,:));
%     end
% 
%     if itrOp
%         yyaxis left;
%         ylim([0 1]);
%         ylabel('Accuracy');
%         yyaxis right;
%         if numClasses == 2
%             ylim([-max(itf) max(itf)]);
%             yticks([-max(itf) 0 max(itf)]);
%             yticklabels({num2str(max(itf)) num2str(0) num2str(max(itf))});
%         else
%             ylim([min(itf) max(itf)]);
%             yticks([min(itf) 0 max(itf)]);
%             yticks(linspace(min(itf), max(itf),6));
%             tempIdx = linspace(1,100,6);
%             yticklabels([-itf(2) -itf(round(tempIdx(2))) ...
%                 itf(round(tempIdx(3))) itf(round(tempIdx(4))) itf(round(tempIdx(5))) ...
%                 itf(round(tempIdx(6)))]);
%             %yticklabels({num2str(max(itf)) num2str(0) num2str(max(itf))});  
%         end
%         ylabel('ITR bits/trial');
%     else
%         ylim([0 1]);
%         ylabel('Accuracy');
%     end
%     if numClasses == 2
%         yline(.5,'--');
%     else
%         yline(0.333,'--');
%     end
%     %xlim([0 7]);
%     title(sprintf('Group Avg: Δ[HbT] Multi'));
%     %legend('S1D1','S1D2','S1D3','S1D4');
%     xlabel('Window Length [s]');
%     hold off;

%     if goodSbj==1 && numClasses == 2
%         fn = sprintf('GrpPerformanceCumsumVsTime_DiffTLen_LR_AllChns_DiffClassifiers_GoodSbjs');
%     elseif goodSbj==1 && numClasses == 3
%         fn = sprintf('GrpPerformanceCumsumVsTime_DiffTLen_AllChns_DiffClassifiers_GoodSbjs');
%     elseif goodSbj==0 && numClasses == 2
%         fn = sprintf('GrpPerformanceCumsumVsTime_DiffTLen_LR_AllChns_DiffClassifiers_BadSbjs');
%     elseif goodSbj==0 && numClasses == 3
%         fn = sprintf('GrpPerformanceCumsumVsTime_DiffTLen_AllChns_DiffClassifiers_BadSbjs');
%     else
%         fn = sprintf('GrpPerformanceCumsumVsTime_DiffStartT_AllChns_DiffClassifiers_AllSbjs');
%     end
    
    fn = sprintf('Fig6_%sClass',num2str(numClasses));

    if saveOp == 1
        print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
    end

end