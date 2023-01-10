function [snr,snr_Mean,snr_Std,cardiacSize,cardiac_Mean,cardiac_Std] = calc_AvgSNR_Welch(sbjNum)

figSaveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum) '\Figures\HRFs'];
saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
rawDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
load([saveDir filesep 'intermediateOutputsCombined_Basis1.mat']);
load([saveDir filesep 'performanceLinearDiscriminantUpdated_LR.mat']);
load([rawDir filesep sbjNum '.mat'],'s');

numConds = 3;
numConc = 3;
numChns = size(dcAvg.measurementList,2)/(numConds*numConc);
sizeCond = size(dcAvg.measurementList,2)/numConds;
sizeChn = sizeCond/numChns;
%condSubPlotIdx = [1 3 5 2 4 6];
condSubPlotIdx = [1 2 3];

colorIdx = [1 3 2 2 1 3];

lineStyle = {'-' '-' '-' '-' '--' '--'};

% srcIdxGrp = {[1 2 3 4] [5 6] [8 9 10 11] [12 13 14 15] [16 17 18 19] [20 21] [23] [25] ...
%     [27 28 29] [31 32 33] [35 36] [37 38] [40] [42]};

srcIdxGrp = {[1 2 3 4] [5 6] [8 9 10 11] [12 13 14 15] [16 17 18 19] [20 21] [23] [25] ...
    [27 28] [30 31] [33 34] [35 36]};

fs = 50;
timePt = (0:0.5*fs:12*fs)+2*fs;

timeIdx = find(timePt==(5*fs+2*fs));

% number of channels, can pull
if strcmp(sbjNum,'08') ||strcmp(sbjNum,'10')
    mlActAuto_Unused = getPrunedChns_Special(sbjNum,s.fName,s.movieList,2);
else
    mlActAuto_Unused = getPrunedChns(sbjNum,s.fName,s.movieList,2);
end

ml690 = mlActAuto_Unused{1}(1:length(mlActAuto_Unused{1})/2);
ml870 = mlActAuto_Unused{1}(length(mlActAuto_Unused{1})/2+1:length(mlActAuto_Unused{1}));
acceptedChns = (ml690==ml870);
%numUnusedChns = length(mlActAuto_Unused{1})/2-sum(ml690==ml870);

% subplot is 8x6.
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
    'S10D15	IPS3/antIPS/IPS4 Right',...
    'S11D13	Superior to IPS3/IPS2/SPL1 Left',...
    'S11D14	IPS4 Left',...
    'S12D13	Superior to IPS3/IPS2/SPL1 Right',...
    'S12D15	IPS4 Right'};

% {'S1D1 Middle Posterior FEF Left',... 21
%     'S1D2 Posterior to sPCS/tgPSC Left',... 20
%     'S1D3 Middle FEF Left',... 15
%     'S1D4 sPCS/tgPCS Left',... 14
%     'S2D1 Posterior FEF Left',... 27
%     'S2D2 Inferior to FEF Left',... 26
%     'S3D3 Anterior FEF Left',... 9
%     'S3D4 iPCS/tgPCS Left',... 8
%     'S3D5 Anterior to FEF Left',... 3
%     'S3D6 iPCS Left',... 2
%     'S4D9 Anterior FEF Right',... 10
%     'S4D10 iPCS/tgPCS Right',... 11
%     'S4D11 Anterior to FEF Right',... 4
%     'S4D12 iPCS Right',... 5
%     'S5D7 Middle Posterior FEF Right',... 22
%     'S5D8 Posterior to sPCS/tgPSC Right',... 23
%     'S5D9 Middle FEF Right',... 16
%     'S5D10 sPCS/tgPCS Right',... 17
%     'S6D7 Posterior to FEF Right',... 28
%     'S6D8 Inferor to FEF Right',... 29
%     'S7D16 Posterior STG/PT Left',... 31
%     'S8D17 Posterior STG/PT Right',... 36
%     'S9D13 IPS3/IPS2/SPL1 Left',... 45
%     'S9D14 IPS3/antIPS/IPS4 Left',... 44
%     'S10D13 IPS3/IPS2/SPL1 Right',... 46
%     'S10D15	IPS3/antIPS/IPS4 Right',... 47
%     'S11D13	Superior to IPS3/IPS2/SPL1 Left',... 39
%     'S11D14	IPS4 Left',... 38
%     'S12D13	Superior to IPS3/IPS2/SPL1 Right',... 40
%     'S12D15	IPS4 Right'}; 41

chnSubPlotInd = [21 20 15 14 27 26 9 8 3 2 10 11 4 5 22 23 16 17 28 29 31 ...
    36 45 44 46 47 39 38 40 41];

colorList = {[0, 0.447, 0.741],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],...
    [0.4940, 0.1840, 0.5560],	[0.4660, 0.6740, 0.1880] 	,[0.3010, 0.7450, 0.9330]};

maxPlot = zeros(1,length(srcIdxGrp));
minPlot = zeros(1,length(srcIdxGrp));

% pxx is freq x chns. w is freq x 1.
[pxx,w] = calcWelch(s.name,s.fName,s.movieList);

% estimated SNR
% convert rad/sample to Hz
% omega * fs * (1/(2*pi)) = Hz
% Hz * 2*pi/fs = omega
numBins = size(pxx,1);
upperHz = 20;
upperNormFreq = upperHz * 2/fs;
upperNormFreqIdx = round(upperNormFreq * numBins);
lowerHz = 0.5;
lowerNormFreq = lowerHz * 2/fs;
lowerNormFreqIdx = round(lowerNormFreq * numBins);

% heart rate of 45 beats/min
lowC = 0.75;
lowerCNormFreq = lowC * 2/fs;
lowerCNormFreqIdx = round(lowerCNormFreq * numBins);
% heart rate of 72 beats/min
highC = 1.2;
upperCNormFreq = lowC * 2/fs;
upperCNormFreqIdx = round(upperCNormFreq * numBins);



abs_pxx =abs(pxx);

snr = 10*log10(sum(abs_pxx(1:lowerNormFreqIdx,:),1)./sum(abs_pxx(upperNormFreqIdx:end,:),1));
snr_Mean = mean(snr);
snr_Std = std(snr);

plot(snr);

cardiacSize = 10*log10(sum(abs_pxx(lowerCNormFreqIdx:upperCNormFreqIdx,:),1));
cardiac_Mean = mean(cardiacSize);
cardiac_Std = std(cardiacSize);

end