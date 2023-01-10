% split into single and multi groups.
% new. Use trials from dcNew var
% For sbj 12 and after

function plotPerformanceBarChartROINames(sbjNum,perfVar,tInd,featureName,saveOp,sortOp,maxTOp)

figSaveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum) '\Figures\Classification'];
rawDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];

load([rawDir filesep sbjNum '.mat'],'s');

fs = 50;
timePt = (0:0.5*fs:12*fs)+2*fs;

if strcmp(sbjNum,'08') || strcmp(sbjNum,'10')
    isNew = 0;
    SSidx = [7,22,24,26,30,34,39,41];
else
    isNew = 1;
    SSidx = [7,22,24,26,29,32];
end

% number of channels, can pull
if strcmp(sbjNum,'08') ||strcmp(sbjNum,'10')
    mlActAuto_Unused = getPrunedChns_Special(sbjNum,s.fName,s.movieList,1.5);
else
    mlActAuto_Unused = getPrunedChns(sbjNum,s.fName,s.movieList,1.5);
end

ml690 = mlActAuto_Unused{1}(1:length(mlActAuto_Unused{1})/2);
ml870 = mlActAuto_Unused{1}(length(mlActAuto_Unused{1})/2+1:length(mlActAuto_Unused{1}));
acceptedChns = (ml690&ml870);

[chnName, LSChnList] = getChnList(isNew);

%perfVar(SSidx,:) = [];

% chnName = {'S1D1','S1D2','S1D3','S1D4','S2D1','S2D2','S3D3','S3D4','S3D5','S3D6',...
%     'S4D9','S4D10','S4D11','S4D12','S5D7','S5D8','S5D9','S5D10','S6D7','S6D8',...
%     'S7D18','S8D19','S9D13','S9D14','S9D16','S10D13','S10D15','S10D17',...
%     'S11D13','S11D14','S12D13','S12D15','S13D28','S14D29'};

% plot performance as bar chart for each channel
figure('units','normalized','outerposition',[0 0 1 1])
%X = categorical(LSChnList);
%X = reordercats(X,LSChnList);

if maxTOp
    perfVar = max(perfVar,[],2);
    tInd = 1;
    titleStr = sprintf('Sbj %s: %s Max Perf',num2str(sbjNum),featureName);
else
    perfVar = perfVar(:,tInd);
    titleStr = sprintf('Sbj %s: %s at %ss',num2str(sbjNum),featureName,num2str(timePt(tInd)));
end

if sortOp
    [perfVar, origInd] = sort(perfVar,'descend');

    XSorted = LSChnList(origInd);
end

XSorted = categorical(XSorted);
XSorted = reordercats(XSorted,string(XSorted));

%bar(XSorted,perfVar);
hold on;
b = bar(XSorted,perfVar,'FaceColor','flat');

for i = 1:length(LSChnList)
    % acceptedChns has 36 chns, fixed it down to 30 chns
    if ~acceptedChns(i)
        b.CData(origInd(i),:) = [0 0.8 0.8];
    end
end

%set(gca,'XTick',[1:1:30]);
set(gca, 'XTickLabel', XSorted);
xtickangle(45);
ylabel('Accuracy [Decimal]');
title(titleStr);
ylim([0 1]);

hold off;
if saveOp == 1
    %fn = sprintf('Sbj%s%ss%s',num2str(sbjNum),featureName,num2str(timePt(tInd)));
    fn = sprintf('BarChart_%s_Basis4_Max',featureName);
    print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
end

end
