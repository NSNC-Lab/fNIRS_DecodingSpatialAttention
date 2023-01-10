function runPlotPerfBarChart(sbjNum,numClasses,feat)

processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];

% pick one, duh. set feat = 0 if don't want to re-run.
switch feat
    case 1
        % cumulative sum from 0, ind channel
        calcCumSum_MultiOnly_AllBasis(sbjNum,numClasses);
    case 2
        % cumulative sum over different interval, ind channel
        calcCumSum_DiffStartT_MultiOnly_AllBasis(sbjNum,numClasses);
    case 3
        % cumulative sum from 0, all channel
        calcCumSum_MultiOnly_AllChn_AllBasis(sbjNum,numClasses);
    case 4
        % cumulative sum over different interval, ind channel
        calcCumSum_MultiOnly_DiffStartT_AllChn_AllBasis(sbjNum,numClasses);
    case 5
        % cumulative sum starting from 1s with slope from -2 to -1s.
        calcCumSumWPreStimSlope_MultiOnly_AllBasis(sbjNum,numClasses)
end

if numClasses == 2
    savePerfFN = 'performanceLinearDiscriminantUpdated_LR.mat';
else
    savePerfFN = 'performanceLinearDiscriminantUpdated.mat';
end

load([processedDataDir filesep savePerfFN],'performanceArrMultiHbO','performanceArrMultiHbR');

saveOp = 1;
sortOp = 1;
maxTOp = 1;
tInd = 1;

plotPerformanceBarChartROINames(sbjNum,performanceArrMultiHbO,tInd,'HbO_CumsumSlope',saveOp,sortOp,maxTOp);
plotPerformanceBarChartROINames(sbjNum,performanceArrMultiHbR,tInd,'HbR_CumsumSlope',saveOp,sortOp,maxTOp);

end