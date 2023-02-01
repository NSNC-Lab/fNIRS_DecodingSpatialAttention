% script to run all. Can take 1-3 days. Parallel computing version is
% available at ...
% Final figures require photo editing.
% Figure number refers to preprint version at 
% https://www.biorxiv.org/content/10.1101/2022.09.06.506821v2.full.pdf

prompt = 'Warning. It'' take 1-3 days to run the whole thing. Are you sure? (Y/N)';
x = input(prompt);

if strcmp(x,'Y')

    sbjList = {'12','13','14','15','16','19','21','22','23','24','25'};
    
    rawDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS';
    
    for i = 1:length(sbjList)
        load([rawDir filesep 'Experiment' sbjList{i} filesep sbjList{i} '.mat'],'s');
        preprocessFNIRS06_CV_GLM_ssBeta_MultiOnly_Version02(s,2,1,1,0.5);
        preprocessFNIRS06_CV_GLM_ssBeta_MultiOnly_Version02(s,3,1,1,0.5);
        preprocessFNIRS06_CV_GLM_ssBeta_DiffTLen_MultiOnly(s,2,1,1,1,0);
        preprocessFNIRS06_CV_GLM_ssBeta_SingleChn_MultiOnly(s,2,1,1);
        preprocessFNIRS06_MultiOnly(s.name,s.fName,s.movieList,0);
        preprocessFNIRS06_LOFO_MultiOnly(s,2,1,1,1,1,0);
    end
    
    load([rawDir filesep 'Experiment' sbjList{i} filesep '08.mat'],'s');
    preprocessFNIRS06_CV_GLM_ssBeta_Version02(s,2,1,1,0.5);
    preprocessFNIRS06_CV_GLM_ssBeta_Version02(s,3,1,1,0.5);
    preprocessFNIRS06_CV_GLM_ssBeta_DiffTLen(s,2,1,1,1,0);
    preprocessFNIRS06_CV_GLM_ssBeta_SingleChn(s,2,1,1);
    preprocessFNIRS06(s.name,s.fName,s.movieList,0);
    preprocessFNIRS06_LOFO(s,2,1,1,1,1,0);
    
    % Figure 3
    % Requires some photo editing
    plotAllHRF_Grp_ChnCriteria_Fig3_TTest(1,0)
    
    % Figure 4
    % Requires some photo editing
    plotAllHRF_Grp_ChnCriteria_Fig4_TTest(1,0)
    
    % Figure 5
    grpPlotPerfVsTime_Fig5_BarPlot(2,1,2);
    grpPlotPerfVsTime_Fig5_BarPlot(3,1,2);
    
    % For original version of Figure 5
    % grpPlotPerfVsTime_Fig5(2,1,0)
    % grpPlotPerfVsTime_Fig5(3,1,0)
    
    % Figure 6
    % This is needed for merging data from multiple subjects and saving it.
    grpPlotPerfVsTime(1,2,1,1,0,0)
    plot_Corr_CAvsBeh(2);
    plot_Corr_CAvsBeh(3);
    
    % Figure 7
    grpLinePlotROINames_04(1,2,1,1);
    
    % Figure 8
    grpPlotPerfVsDiffTLen_Fig6(1,2,1,0,1);
end