% plot correlation between time of trial and CA performance at 2-3s
% window

% do the same thing for var of beta coeff.

function plot_Corr_CAvsAvgTime(numClasses)

    sbjList = {'08','12','13','14','15','16','19','21','22','23','24','25'};

    figSaveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll\Figures\Classification';
    figDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll';

    %fn = [figDir filesep 'grp_Names'];
%     fn = [figDir filesep 'sbjListPerf_DiffTLen_All.mat'];
%     load(fn,'sbjListPerfHbO','sbjListPerfHbR','sbjListPerfHbT');
    
    fs = 50;
    %timePt = (0:0.25*fs:5*fs)+2*fs;
    timePt = (-2*fs:0.25*fs:5*fs);
    
    avgTime = zeros(1,length(sbjList));
    
    for i = 1:length(sbjList)
        tempDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' ...
            sbjList{i} filesep 'summaryStat_1.5.mat'];
        
        load(tempDir,'avgTimeTr');
        avgTime(i) = avgTimeTr; 
    end
    
    % from grpPlotPerfVsTime.m
    if numClasses == 2
        dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
            filesep 'sbjListPerfAll_GLM_CV_SSBeta_LR_RejTr.mat'];
%         dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
%             filesep 'sbjListPerfAll_GLM_CV_SSBeta_LR_RejTr_10Hz.mat'];
    else
        dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
            filesep 'sbjListPerfAll_GLM_CV_SSBeta_RejTr.mat'];
%         dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
%             filesep 'sbjListPerfAll_GLM_CV_SSBeta_RejTr_10Hz.mat'];
    end
    load(dataFN,'sbjListPerfHbO','sbjListPerfHbR','sbjListPerfHbT');

    idxT = timePt == (2+2)*fs;
    
    %figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    figure();
    mdl = fitlm(avgTime,sbjListPerfHbT(:,idxT));
    tVal = mdl.Coefficients.tStat(2);
    x = anova(mdl);
    h = plot(mdl);
    h(1).Marker = 'o';
    h(1).MarkerSize = 8;
    h(1).MarkerFaceColor = 'black';
    h(1).MarkerEdgeColor = 'black';
    h(2).LineWidth = 1.5;
    h(3).LineWidth = 1.5;
    h(4).LineWidth = 1.5;
    xlabel('Average Time of Trial [s]');
    ylabel('Classification Accuracy');
    title('Average Time vs HbT CA at 4s');
%     annotation('textbox',[0.2 0.15 0.1 0.1], 'String', ['p Value: ' num2str(x.pValue(1)) ...
%         newline 'R Squared: ' num2str(mdl.Rsquared.Ordinary) newline ...
%         'R Squared Adjusted: ' num2str(mdl.Rsquared.Adjusted)]);
    annotation('textbox',[0.2 0.15 0.1 0.1], 'String', sprintf('p Value: %0.4f\n R Squared: %0.2f\n R Squared Adjusted: %0.2f', x.pValue(1) ...
        ,mdl.Rsquared.Ordinary, mdl.Rsquared.Adjusted));
    
    %scatter(behScore,sbjListPerfHbO(1,:,idxT));

end