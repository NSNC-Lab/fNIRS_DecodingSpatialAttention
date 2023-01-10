% plot raw, raw bandpassed under 0.5 Hz, raw bandpassed under 3 Hz, and raw
% median filtered. Explore other filtering options, linear and nonlinear filters.

% plot delta optical density unfiltered, bandpassed filtered under 0.5
% Hz, under 3 Hz, median filtered, and gaussian-filter.

% Butterworth, filter order 3.

% we will not compare motion correction algorithms here. Beyond the scope.

function plotSignal_MultiOnly(sbjNum,rawDataFN)

saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
% Convert nirs to snirf file format
snirf1 = SnirfClass(load([rawDataFN{1} '.nirs'],'-mat'));
%snirf3 = SnirfClass(load([rawDataFN{2} '.nirs'],'-mat'));
%snirf4 = SnirfClass(load([rawDataFN{3} '.nirs'],'-mat'));
snirf1.Info();

numChn = 30;
chnName = getChnName(numChn);

data = snirf1.data;
probe = snirf1.probe;
mlActMan = {};
tIncMan = {};
%dod = snirf1.dod;
%mlActAuto = snirf1.mlActAuto;
%stim = snirf1.stim;
%tIncAuto = snirf1.tIncAuto;
%dc = snirf1.dc;
%Aaux = [];
%rcMap = [];

%mlActAuto_Unused = hmrR_PruneChannels(data,probe,mlActMan,tIncMan,[10000  10000000],2,[0  45]);

mlActAuto = {ones(size(snirf1.data.measurementList,2),1)};
dod = hmrR_Intensity2OD(data);

% For sbj 12, params in 1st line for motion correct/artifact perform
% significantly better than params in 2nd line.
dod = hmrR_MotionCorrectPCArecurse(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,15,4,0.97,5,1);
%[dod,tInc,svs,nSV,tInc0] = hmrR_MotionCorrectPCArecurse(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,20,5,0.97,5,1);

dod_T = dod.GetTime();

dod_BP_p5_DC = hmrR_BandpassFilt(dod,0.01,0.5);
dod_BP_3_DC = hmrR_BandpassFilt(dod,0.01,3);

dod_raw = dod.GetDataTimeSeries();

dod_mf_data = medfilt1(dod_raw,11);
dod_mf_100_data = medfilt1(dod_raw,100);

gaussFilt = gausswin(100,2.5);
dod_gf_data = filter(gaussFilt,1,dod_raw);

dod_bp_p5 = dod_BP_p5_DC.GetDataTimeSeries();
dod_bp_3 = dod_BP_3_DC.GetDataTimeSeries();

dod_MedianFilt = DataClass(dod_mf_data,dod_T,snirf1.data.measurementList);

dod_MF_BP_p5_DC = hmrR_BandpassFilt(dod_MedianFilt,0.01,0.5);
dod_MF_BP_p5 = dod_MF_BP_p5_DC.GetDataTimeSeries();
%dod_GaussianFilt = DataClass(dod_gf_data,dod.data.dataTimeSeries,snirf1.data.measurementList);

%dod = hmrR_BandpassFilt(dod,0,0.5);
%time x concentration x channels
ssIdx = [7,22,24,26,29,32];

dod_raw = reshape(dod_raw,[size(dod_raw,1),size(dod_raw,2)/2,2]);
dod_bp_p5 = reshape(dod_bp_p5,[size(dod_bp_p5,1),size(dod_bp_p5,2)/2,2]);
dod_bp_3 = reshape(dod_bp_3,[size(dod_bp_3,1),size(dod_bp_3,2)/2,2]);
dod_mf_data = reshape(dod_mf_data,[size(dod_mf_data,1),size(dod_mf_data,2)/2,2]);
dod_mf_100_data = reshape(dod_mf_100_data,[size(dod_mf_100_data,1),size(dod_mf_100_data,2)/2,2]);
dod_gf_data = reshape(dod_gf_data,[size(dod_gf_data,1),size(dod_gf_data,2)/2,2]);
dod_MF_BP_p5 = reshape(dod_MF_BP_p5,[size(dod_MF_BP_p5,1),size(dod_MF_BP_p5,2)/2,2]);

fs = 50;
xTempS = [50 65]*fs;
xTempFS = xTempS(1):xTempS(2);
tTemp = dod_T(xTempS(1):xTempS(2));

dod_raw(:,ssIdx,:) = [];
dod_bp_p5(:,ssIdx,:) = [];
dod_bp_3(:,ssIdx,:) = [];
dod_mf_data(:,ssIdx,:) = [];
dod_mf_100_data(:,ssIdx,:) = [];
dod_gf_data(:,ssIdx,:) = [];
dod_MF_BP_p5(:,ssIdx,:) = [];

dod_raw = dod_raw(xTempFS,:,:);
dod_bp_p5 = dod_bp_p5(xTempFS,:,:);
dod_bp_3 = dod_bp_3(xTempFS,:,:);
dod_mf_data = dod_mf_data(xTempFS,:,:);
dod_mf_100_data = dod_mf_100_data(xTempFS,:,:);
dod_gf_data = dod_gf_data(xTempFS,:,:);
dod_MF_BP_p5 = dod_MF_BP_p5(xTempFS,:,:);

for iChn = 1:size(dod_raw,2)
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    yTemp = [dod_raw(:,iChn,2) dod_bp_p5(:,iChn,2) dod_bp_3(:,iChn,2)...
        dod_mf_data(:,iChn,2) dod_mf_100_data(:,iChn,2) dod_gf_data(:,iChn,2) dod_MF_BP_p5(:,iChn,2)];
    
    yLabels = {'Raw','BP 0.5 Hz','BP 3 Hz','Median 11th','Median 100th','Gaussian','Median & BP 0.5Hz'};
    
    stackedplot(tTemp,yTemp,"Title",chnName(iChn),"DisplayLabels",yLabels);
    
end