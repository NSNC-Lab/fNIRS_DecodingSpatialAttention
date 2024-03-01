% For figure 6 in Draft 3
% 2 subplots (1,2). HbT. Show only high performing group.
% Only show LDA classifier.
% 2-class classification at left panel.
% 3-class classification at right panel.
% look at 7-color scheme for readability.
% 08/01/2022, add CI options
%   1 for counting each fold as one sample point
%   2 for counting each sbj as one sample point

% STATUS: active
% 
% SYNTAX:
% grpPlotPerfVsDiffTLen_Fig6(goodSbj,numClasses,saveOp,itrOp,ciOp)
% 
% DESCRIPTION:
% Plot decoding performance (accuracies) as a function of decision window
%   length.
% 
% RESTRICTION:
% None.
% 
% INPUTS:
% goodSbj - int: group of subjects. Encoding as follow:
%       1 - high-performance group
%       2 - low-performance group
%       0 - all subjects
% numClasses - int: number of classes to classify
%       2 - classify between left & right
%       3 - classify between left, right & center
% saveOp - int: option to save
%       0 - don't save
%       1- save
% itrOp - int: option to display information transfer rate (itr)
%       0 - don't display itr.
%       1- display itr.
% ciOp - int: option to display error bar as confidence intervals:
%       0 - don't display error bars as confidence intervals
%       1 - display error bars as confidence intervals: counting each fold 
%           as one sample point
%       2 - display error bars as confidence intervals: counting each sbj
%           as one sample point
%
% RETURNED VARIABLES:
% None.
% 
% FILES SAVED:Matthew_nir_data_allroi
% 1) save figure of decoding performance as a function of time course of trial
% 2) save decoding performances of all subjects in one array variable.
% 
% PLOTTING:
% Line plots of decoding performance as a function of time course of trial
%grpPlotPerfVsDiffTLen_Fig6_allroi(1,2,0,0,2)

function grpPlotPerfVsDiffTLen_Fig6_allroi(goodSbj,numClasses,saveOp,itrOp,ciOp)

    if goodSbj == 1;
        % eligible if less than 20 channels rejected and audio and video behavioral perf
        % greater than 50% : subj 10 has low behav perf (less than 50%)
        sbjList = {'08','12','13','16','19','21','24'};
        dataFN = ['H:\My Drive\fNIRS\sudan_final_all_roi_data_paper\GroupResultAll' ...
            filesep 'sbjListPerf_DiffTLen_Good.mat'];
    elseif goodSbj == 0
        sbjList = {'14','15','22','23','25'};
        dataFN = ['H:\My Drive\fNIRS\sudan_final_all_roi_data_paper\GroupResultAll' ...
            filesep 'sbjListPerf_DiffTLen_Bad.mat'];
    else
        sbjList = {'08','12','13','14','15','16','19','21','22','23','24','25'};
        dataFN = ['H:\My Drive\fNIRS\sudan_final_all_roi_data_paper\GroupResultAll' ...
            filesep 'sbjListPerf_DiffTLen_All.mat'];
    end

    figSaveDir = 'H:\My Drive\fNIRS\sudan_final_all_roi_data_paper\GroupResultAll\Figures';

    fs = 50;
    timeLen = [0.1 0.2 0.5 1 1.5 2 3 4 5].*fs;
    zeroT = 2*fs;
    numSubsets = 2; % was 7 for 7 rois
  
    sbjListPerfHbO = zeros(numSubsets,5*10*length(sbjList),length(timeLen));
    sbjListPerfHbR = zeros(numSubsets,5*10*length(sbjList),length(timeLen));
    sbjFoldListPerfHbT = zeros(numSubsets,5*10*length(sbjList),length(timeLen));
    
    sbjListPerfHbT = zeros(numSubsets,length(sbjList),length(timeLen));
    
    for i = 1:length(sbjList)
        sbjNum = sbjList{i};

        processedDataDir = ['H:\My Drive\fNIRS\sudan_final_all_roi_data_paper\Experiment' num2str(sbjNum)];

        if numClasses == 2
            savePerfFN = 'performance_GLM_CV_SSBeta_LR_DiffTLen_RejTr_10Hz';
        else
            savePerfFN = 'performance_GLM_CV_SSBeta_RejTr_DiffTLen_.mat';
        end

        load([processedDataDir filesep savePerfFN],...
            'performanceLDALedoitHbT_All',...
            'performanceLDALedoitHbT_LRFEF');
        
        startIdx = (i-1)*50+1;
        endIdx = i*50;
        
        sbjFoldListPerfHbT(1,startIdx:endIdx,:) = performanceLDALedoitHbT_All;
        
        %sbjFoldListPerfHbT(2,startIdx:endIdx,:) = performanceLDALedoitHbT_LRFEF;
        
        sbjListPerfHbT(1,i,:) = mean(performanceLDALedoitHbT_All,1);

      %sbjListPerfHbT(2,i,:) = mean(performanceLDALedoitHbT_LRFEF,1);

    end
    
    save(dataFN,'sbjListPerfHbO','sbjListPerfHbR','sbjListPerfHbT');


    grpPerfHbT = sum(sbjFoldListPerfHbT,2) ./ sum(sbjFoldListPerfHbT~=0,2);
    
    if ciOp == 1
        grpPerfHbTSE = squeeze(std(sbjFoldListPerfHbT,[],2))./sqrt(sum(sbjFoldListPerfHbT~=0,2));
    elseif ciOp == 2
        grpPerfHbTSE = squeeze(std(sbjListPerfHbT,[],2))./squeeze(sqrt(sum(sbjListPerfHbT~=0,2)));
    end
    
    conf = 0.95;

    %nSz = size(grpPerfHbTSE,1);
    nSz = length(sbjList);
    grpPerfCIHbT = zeros(2,numSubsets,length(timeLen));
    
    for i = 1:numSubsets
        for idxTime = 1:length(timeLen)

            grpPerfCIHbT(:,i,idxTime) = calcCITDist(nSz,grpPerfHbTSE(i,idxTime),conf);

        end
    end
    
    [h,p,ci,stats] = ttest(sbjListPerfHbT(1,:,4),0.5);
    disp(p);
    disp(ci);
    disp(stats);
    xSig = tinv(0.975,stats.df)*stats.sd/sqrt(size(sbjListPerfHbT,2)) + 0.5;
    disp(xSig);
    
    %subsetIdx = [2];

    figure();hold on;
    
    % for ALL
    errorbar(squeeze(timeLen./fs),squeeze(grpPerfHbT(1,1,:)),abs(squeeze(grpPerfCIHbT(1,1,:))),abs(squeeze(grpPerfCIHbT(2,1,:))),...
        'Color',[0 0 0],'LineWidth',1.5)

    % for LRFEF
        %errorbar(squeeze(timeLen./fs),squeeze(grpPerfHbT(2,1,:)),abs(squeeze(grpPerfCIHbT(2,1,:))),abs(squeeze(grpPerfCIHbT(2,2,:))),...
         %'Color',[0 0 0],'LineWidth',1.5)


    %for i = 1:length(subsetIdx)
        %plot(timeLen./fs,squeeze(grpPerfHbT(subsetIdx(i),1,:)))
    %end

    yline(xSig,'Color',[243 21 185]./255,'LineStyle',':','LineWidth',1.8);

    x = 0:8; 

    %y_lower = repmat(ci(1), size(squeeze(timeLen./fs)));
    %y_upper = repmat(ci(2), size(squeeze(timeLen./fs)));

    %fill([x, fliplr(x)], [y_lower, fliplr(y_upper)], 'r', 'FaceAlpha', 0.1,'LineStyle','none');
    
    grid;
    xticks([0.2 0.5 1 1.5 2 3 4 5]);

    legend({'All ROIs', 'Chance level'},'Location','Northwest');
    xlabel('Window Length [s]');
    ylim([0.3 1])
    ylabel('Accuracy')
    hold off;
    set(gca, 'FontName', 'Arial', 'FontSize', 18);
    title("Average CV accuracy across the subjects")

    
end