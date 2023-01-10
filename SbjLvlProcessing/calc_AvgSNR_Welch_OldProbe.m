function [snr,snr_Mean,snr_Std,cardiacSize,cardiac_Mean,cardiac_Std] = calc_AvgSNR_Welch_OldProbe(sbjNum)

rawDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
load([rawDir filesep sbjNum '.mat'],'s');

% pxx is freq x chns. w is freq x 1.
[pxx,w] = calcWelch_OldProbe(s.name,s.fName,s.movieList);

fs = 50;

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

snr = 10*log10(sum(abs_pxx(1:lowerNormFreqIdx,:))./sum(abs_pxx(upperNormFreqIdx:end,:)));
snr_Mean = mean(snr);
snr_Std = std(snr);

plot(snr);

cardiacSize = 10*log10(sum(abs_pxx(lowerCNormFreqIdx:upperCNormFreqIdx,:),1));
cardiac_Mean = mean(cardiacSize);
cardiac_Std = std(cardiacSize);

end