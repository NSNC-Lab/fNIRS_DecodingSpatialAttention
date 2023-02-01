% script to run all. Fast.
% Generate figures for the manuscript.
% Final figures require photo editing.
% Figure number refers to preprint version at 
% https://www.biorxiv.org/content/10.1101/2022.09.06.506821v2.full.pdf

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
