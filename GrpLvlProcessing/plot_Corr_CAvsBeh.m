% plot correlation between behavioral score and CA performance at 2-3s
% window

% do the same thing for var of beta coeff.

function plot_Corr_CAvsBeh(numClasses)

    sbjList = {'08','12','13','14','15','16','19','21','22','23','24','25'};

    figSaveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll\Figures\Classification';
    figDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll';

    %fn = [figDir filesep 'grp_Names'];
    %load(fn,'sbjListPerfHbO','sbjListPerfHbR','sbjListPerfHbT');
    
    fs = 50;
    timePt = (-2*fs:0.25*fs:5*fs);
    
    behScore = zeros(1,length(sbjList));
    
    for i = 1:length(sbjList)
        tempDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' ...
            sbjList{i} filesep 'responses_' sbjList{i}];
        
        load(tempDir,'responsesA','responsesV','correctRespA','correctRespV');
        behScore(i) = sum((responsesA==correctRespA)&(responsesV==correctRespV))/length(responsesA); 
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
    load(dataFN,'sbjListPerfHbT');

    idxT = timePt == (3)*fs;
    
    %figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    figure();
    
    mdl = fitlm(behScore,sbjListPerfHbT(:,idxT));
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
    xlabel('Behavioral Score');
    ylabel('Classification Accuracy');
    title('Behavioral Score vs HbT CA at 3-4s');
%     annotation('textbox',[0.2 0.15 0.1 0.1], 'String', sprintf('p Value: %0.4f\n R Squared: %0.2f\n R Squared Adjusted: %0.2f', x.pValue(1) ...
%         ,mdl.Rsquared.Ordinary, mdl.Rsquared.Adjusted));
    annotation('textbox',[0.2 0.15 0.1 0.1], 'String', sprintf('p Value: %0.4f\nR Squared Adjusted: %0.2f', x.pValue(1) ...
        , mdl.Rsquared.Adjusted));
    
    %scatter(behScore,sbjListPerfHbO(1,:,idxT));

end