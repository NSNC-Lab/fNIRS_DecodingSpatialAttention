% Arithmetic average of HRFs
% Light line color, plot ind HRFs
% Different processing for sbj 08 and 10 because they use old SD design and
% has both single and multi conditions
% old codes

function grpPlotHRF(sbjList)

figSaveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll\Figures\Classification';

fs = 50;
opPlotInd = 1;
t = -2:1/fs:15;

colorIdx = [1 3 2];

colorList = {[0, 0.447, 0.741],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],...
    [0.4940, 0.1840, 0.5560],	[0.4660, 0.6740, 0.1880] 	,[0.3010, 0.7450, 0.9330]};

colorListLight = {[0, 0.447, 0.741, 0.5],[0.8500, 0.3250, 0.0980, 0.5],[0.9290, 0.6940, 0.1250, 0.5],...
    [0.4940, 0.1840, 0.5560],	[0.4660, 0.6740, 0.1880] 	,[0.3010, 0.7450, 0.9330]};

condSubPlotIdx = [1 2 3];

lineStyle = {'-' '-' '-'};

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
    'S10D15 IPS3/antIPS/IPS4 Right',...
    'S11D13 Superior to IPS3/IPS2/SPL1 Left',...
    'S11D14 IPS4 Left',...
    'S12D13 Superior to IPS3/IPS2/SPL1 Right',...
    'S12D15 IPS4 Right'};

srcIdxGrp = {[1 2 3 4] [5 6] [8 9 10 11] [12 13 14 15] [16 17 18 19] [20 21] [23] [25] ...
    [27 28] [30 31] [33 34] [35 36]};

numChn = length(chnName);
origNumConds = 6;
numConds = 3;
numConc = 3;

totNumChns_Homer3 = 36*3*3;

indHb = zeros(size(t,2),totNumChns_Homer3,length(sbjList));
%indHbR = zeros(size(t,2),numChn,length(sbjList));
%indHbT = zeros(size(t,2),length(sbjList));

chnNum = 0;

for i = 1:length(sbjList)
    sbjNum = sbjList{i};
    saveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
    load([saveDir filesep 'intermediateOutputsCombined_Basis1.mat'],'dcAvg','mlActAuto');
    
    hrf_GLM = dcAvg.dataTimeSeries;
    mlList = dcAvg.measurementList;
    
    if strcmp(sbjNum,'08') || strcmp(sbjNum,'10')
        isSDNew = 0;
    else
        isSDNew = 1;
    end
    
    % after rmvSingleCond, cut down from 756 chns to 378 chns
    
    
    % after convert2SD2, cut down from 378 chns to 324 chns
    if ~isSDNew
        [hrf_GLM,mlList] = rmvSingleCond(hrf_GLM,origNumConds,mlList);
        [hrf_GLM,mlList,mlActAuto] = convert2SD2_cnt(hrf_GLM,mlList,mlActAuto{1});
    else
        mlActAuto = mlActAuto{1};
    end
    
%     opRmvSS = 1;
% 
%     % after selectRS, cut down from 324 to 270 channels
%     [hrf_GLM,mlList,mlActAuto] = selectRS_cnt(hrf_GLM,opRmvSS,mlList,mlActAuto);
    
    sizeCond = size(mlList,2)/numConds;
    numChns = size(mlList,2)/(numConds*numConc);
    sizeChn = sizeCond/numChns;
    
    % dcAvg.dataTimeSeries is time x channel
    for i1 = 1:length(srcIdxGrp)
        for j = 1:length(srcIdxGrp{i1})
            for i2 = 1:3
                % dcAvg.dataTimeSeries is time x channels
                % plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),lineStyle{i2},'Color',colorList{colorIdx(i2)});hold on;
                indHb(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1,i) = ...
                    hrf_GLM(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1);
                indHb(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2,i) = ...
                    hrf_GLM(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2);
%                 indHbR(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2,i) = ...
%                     dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2);
            end
        end
    end
end

% grpHbO = indHbO./length(sbjList);
% grpHbR = indHbR./length(sbjList);
grpHb = squeeze(mean(indHb,3));
%grpHbR = squeeze(mean(indHbR,1));

% now plot
for i1 = 1:length(srcIdxGrp)
    for j = 1:length(srcIdxGrp{i1})
        figure('units','normalized','outerposition',[0 0 1 1]); hold on;
        chnNum = chnNum + 1;
        subplot(1,2,1);hold on;
        for i2 = 1:3
            % dcAvg.dataTimeSeries is time x channels
            if opPlotInd
                chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1;
                for k= 1:length(sbjList)
                    %sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1
                    
                    plot(t,indHb(:,chnIdx,k),'Color',colorListLight{colorIdx(i2)},'LineStyle','--');hold on;
                end
            end

            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 1),lineStyle{i2},'Color',colorList{colorIdx(i2)});hold on;
            plot(t,grpHb(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)});hold on;

            legend({'Left','Right','Center'});
            ylabel('[M*mm]')
            xlabel('Time [s]');
            xlim([-2 15]);
            %label = num2str(dcavg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).dataTypeLabel);
            %condStr = snirf1.stim(1,i2).name;
            title(sprintf('HbO %s',chnName{chnNum}));
        end

        % HbR
        subplot(1,2,2);hold on;
        for i2 = 1:3
            % dcAvg.dataTimeSeries is time x channels
            if opPlotInd
                chnIdx = sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2;
                for k = 1:length(sbjList)
                    plot(t,indHb(:,chnIdx,k),'Color',colorListLight{colorIdx(i2)},'LineStyle','--');hold on;
                end
            end

            %plot(dcAvg.time, dcAvg.dataTimeSeries(:,sizeCond*(condSubPlotIdx(i2)-1) + sizeChn*(srcIdxGrp{i1}(j)-1) + 2),lineStyle{i2},'Color',colorList{colorIdx(i2)});hold on;
            plot(t,grpHb(:,chnIdx),lineStyle{i2},'Color',colorList{colorIdx(i2)});hold on;

            legend({'Left','Right','Center'});
            ylabel('[M*mm]')
            xlabel('Time [s]');
            xlim([-2 15]);
            %label = num2str(dcavg.measurementList(1,sizeCond*(i-1) + sizeChn*(j-1) + 1).dataTypeLabel);
            %condStr = snirf1.stim(1,i2).name;
            title(sprintf('HbR %s',chnName{chnNum}));
        end
        srcStr = num2str(mlList(1,chnIdx).sourceIndex);
        detectorStr = num2str(mlList(1,chnIdx).detectorIndex);
        hold off;
        fn = sprintf('HRF_S%sD%s',srcStr,detectorStr);
        print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
        close;
    end
end



end