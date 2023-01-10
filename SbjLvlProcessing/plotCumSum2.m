% plot performance as bar chart for each channel

saveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment03';
saveDir2 = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment03';
load([saveDir filesep 'intermediateOutputs.mat']);
load([saveDir filesep 'singleTrialsUpdated.mat'],'singleTrialHRFCell');
load([saveDir2 filesep 'iniList_alh_202154_.mat']);
load([saveDir filesep 'performanceLinearDiscriminantUpdated.mat']);
%load([saveDir filesep 'performanceLR.mat']);
%savePerfFN = 'performance.mat';

srcIdxGrp = {[1 2 3 4] [5 6] [8 9 10 11] [12 13 14 15] [16 17 18 19] [20 21] [23] [25] ...
    [27 28 29] [31 32 33] [35 36] [37 38] [40] [42]};

srcPairIdxGrp = {[1 5] [2 6] [3 4] [7 8] [9 10] [11 12] [13 14]};
    
chnName = {'S1D1','S1D2','S1D3','S1D4','S2D1','S2D2','S3D3','S3D4','S3D5','S3D6',...
    'S4D9','S4D10','S4D11','S4D12','S5D7','S5D8','S5D9','S5D10','S6D7','S6D8',...
    'S7D18','S8D19','S9D13','S9D14','S9D16','S10D13','S10D15','S10D17',...
    'S11D13','S11D14','S12D13','S12D15','S13D28','S14D29'};

chnNameGrp = {{'S1D1','S1D2','S1D3','S1D4'},{'S2D1','S2D2'},{'S3D3','S3D4','S3D5','S3D6'},...
    {'S4D9','S4D10','S4D11','S4D12'},{'S5D7','S5D8','S5D9','S5D10'},{'S6D7','S6D8'},...
    {'S7D18'},{'S8D19'},{'S9D13','S9D14','S9D16'},{'S10D13','S10D15','S10D17'},...
    {'S11D13','S11D14'},{'S12D13','S12D15'},{'S13D28'},{'S14D29'}};

% chnInd = ([1 2 3 4 5 6 8 9 10 11 12 13 14 15 16 17 18 19 20 21 23 25 27 28 29 ...
%     31 32 33 35 36 37 38 40 42]- 1) *3 + 3;

chnInd = [1 2 3 4 5 6 8 9 10 11 12 13 14 15 16 17 18 19 20 21 23 25 27 28 29 ...
    31 32 33 35 36 37 38 40 42];


fs = 50;
timePt = (0:0.5*fs:12*fs)+2*fs;

numConds = 6;
numConc = 3;
numChns = size(dcAvg.measurementList,2)/(numConds*numConc);
sizeCond = size(dcAvg.measurementList,2)/numConds;
sizeChn = sizeCond/numChns;

colorList = {[0, 0.447, 0.741],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560] 	};

tInd = find(timePt==(6.5*fs+2*fs));

figure();
X = categorical(chnName);
X = reordercats(X,chnName);

[resultSorted, origInd] = sort(performanceArrSingle(chnInd,tInd),'descend');

XSorted = X(origInd);

%bar(X,performanceArrSingle(chnInd,tInd));
bar(resultSorted);
set(gca, 'XTick',1:1:length(XSorted),'XTickLabel', XSorted);
xtickangle(45)
title('3-Classes Performance for Single Condition at t=6.5s HbO');
ylim([0 1]);
%yline(0.5);
yline(0.3333);

figure();
X = categorical(chnName);
X = reordercats(X,chnName);

[resultSorted, origInd] = sort(performanceArrMulti(chnInd,tInd),'descend');

XSorted = X(origInd);

bar(resultSorted);
set(gca, 'XTick',1:1:length(XSorted),'XTickLabel', XSorted);
xtickangle(45)
title('3-Classes Performance for Multi Condition at t=6.5s HbO');
ylim([0 1]);
%yline(0.5);
yline(0.3333);

% plot confusion matrices. One for best performing, one for chance level
% and one for non-0 worst performing.

% Best is chn #8, chance level is chn #29, worst #4
figure();
responseASingle = squeeze(singleTrials(1,1,:));
CBest = confusionmat(responseASingle,squeeze(predRespS(3,tInd,:)));
confusionchart(CBest);
title('Single. S1D3 HbO: Best-Performing');

figure();
CChance = confusionmat(responseASingle,squeeze(predRespS(29,tInd,:)));
confusionchart(CChance);
title('Single. S9D16 HbO: Chance-Level');

figure();
CWorst = confusionmat(responseASingle,squeeze(predRespS(14,tInd,:)));
confusionchart(CWorst);
title('Single. S4D12 HbO: Worst-Performing');

% Best is chn #8, chance level is chn #29, worst #4
figure();
responseAMulti = squeeze(multiTrials(1,1,:));
CBest = confusionmat(responseAMulti,squeeze(predRespM(3,tInd,:)));
confusionchart(CBest);
title('Multi. S1D3 HbO: Best-Performing');

% figure();
% CChance = confusionmat(responseAMulti,squeeze(predRespM(29,tInd,:)));
% confusionchart(CChance);
% title('Multi. S9D16 HbO: Chance-Level');

figure();
CWorst = confusionmat(responseAMulti,squeeze(predRespM(10,tInd,:)));
confusionchart(CWorst);
title('Multi. S3D5 HbO: Worst-Performing');
