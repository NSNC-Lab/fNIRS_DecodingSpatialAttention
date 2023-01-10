% spaceOp
% 1: original feature space
% 2: PCA transformed space
% 3: log transformed space
function plotHistogramFeats(sbjNum,numClasses,spaceOp,saveOp)

processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];
rawDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
figSaveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' ...
    num2str(sbjNum) '\Figures\Classification'];

fs = 50;
timePt = (0:0.25*fs:5*fs)+2*fs;
zeroT = 2*fs;
numChn = 30;
cmap = jet(numChn);
if strcmp(sbjNum,'07')||strcmp(sbjNum,'08')||strcmp(sbjNum,'10')
    isSDNew = 0;
else
    isSDNew = 1;
end

chnName = getChnName(numChn);

% Too lazy to modularize this
% Basis 1
load([processedDataDir filesep 'intermediateOutputsCombined_Basis1.mat']);
% Each array is channels x time x trial
if numClasses == 2
    load([processedDataDir filesep 'singleTrialsUpdated_LR_Basis1.mat'],'singleTrialHRFHbOM',...
        'singleTrialHRFHbRM',...
        'singleTrialHRFHbTM','indexMoviesTest');
    
    yLimAxis = [0 1];
else
    load([processedDataDir filesep 'singleTrialsUpdated_Basis1.mat'],'singleTrialHRFHbOM',...
        'singleTrialHRFHbRM',...
        'singleTrialHRFHbTM','indexMoviesTest');
    yLimAxis = [0 1];
end

if strcmp(sbjNum,'15')
    indexMoviesTest = indexMoviesTest(2:end,:);
end

% Offset
singleTrialHRFHbOM = offsetTrials(singleTrialHRFHbOM,zeroT);
singleTrialHRFHbRM = offsetTrials(singleTrialHRFHbRM,zeroT);
singleTrialHRFHbTM = offsetTrials(singleTrialHRFHbTM,zeroT);

% after convert2SD2, cut down from 42 chns to 36 chns. Remove old chns
if ~isSDNew
    [singleTrialHRFHbOM,mlList] = convert2SD2(singleTrialHRFHbOM,mlActAuto{1});
    [singleTrialHRFHbRM,~] = convert2SD2(singleTrialHRFHbRM);
    [singleTrialHRFHbTM,~] = convert2SD2(singleTrialHRFHbTM);
else
    mlList = mlActAuto{1};
end

%isSDNew = 1;

% after selectRS, cut down from 36 to 30 chns. Remove SS
[singleTrialHRFHbOM,mlList] = selectRS(singleTrialHRFHbOM,1,mlList);
[singleTrialHRFHbRM,~] = selectRS(singleTrialHRFHbRM,1);
[singleTrialHRFHbTM,~] = selectRS(singleTrialHRFHbTM,1);

behFN = [rawDataDir filesep 'responses_' sbjNum];

if ~isSDNew
    [singleTrialHRFHbOM,movieIdx] = keepCorrectTrials_Special(sbjNum,singleTrialHRFHbOM,behFN,numClasses,indexMoviesTest);
    [singleTrialHRFHbRM] = keepCorrectTrials_Special(sbjNum,singleTrialHRFHbRM,behFN,numClasses,indexMoviesTest);
    [singleTrialHRFHbTM] = keepCorrectTrials_Special(sbjNum,singleTrialHRFHbTM,behFN,numClasses,indexMoviesTest);
else
    [singleTrialHRFHbOM,movieIdx,behScore] = keepCorrectTrials(sbjNum,singleTrialHRFHbOM,behFN,numClasses,indexMoviesTest);
    [singleTrialHRFHbRM] = keepCorrectTrials(sbjNum,singleTrialHRFHbRM,behFN,numClasses,indexMoviesTest);
    [singleTrialHRFHbTM] = keepCorrectTrials(sbjNum,singleTrialHRFHbTM,behFN,numClasses,indexMoviesTest);
end

lenT = 1.5*fs;

timeIdx = find(timePt==(1.5+2)*fs);

%for i = 2:length(timePt)
for i = timeIdx:timeIdx
    
%     temp = cumsum(singleTrialHRFHbOM(:,timePt(i):timePt(i)+lenT,:),2);
%     cumsumHbOM = squeeze(temp(:,end,:));
%     
% %     temp = cumsum(singleTrialHRFHbRS(:,timePt(i):timePt(i)+lenT,:),2);
% %     cumsumHbRS = squeeze(temp(:,end,:));
%     
%     temp = cumsum(singleTrialHRFHbRM(:,timePt(i):timePt(i)+lenT,:),2);
%     cumsumHbRM = squeeze(temp(:,end,:));
%     
% %     temp = cumsum(singleTrialHRFHbTS(:,startT:timePt(i),:),2);
% %     cumsumHbTS = squeeze(temp(:,end,:));
    
    temp = cumsum(singleTrialHRFHbTM(:,timePt(i):timePt(i)+lenT,:),2);
    cumsumHbTM = squeeze(temp(:,end,:));
    if spaceOp == 1
        strTitle = 'Original Feature Space';
    elseif spaceOp == 2
        [coeff,cumsumHbTM] = pca(cumsumHbTM');
        cumsumHbTM = cumsumHbTM';
        strTitle = 'PCA-Transformed Feature Space';
    elseif spaceOp == 3
%         cumsumHbTM = log(cumsumHbTM);
%         
%         
%         logRegModel = fitclinear(...
%             cumsumHbTM, ...
%             trainingResponse, ...
%             'Learner', 'logistic', ...
%             'Regularization', 'ridge', ...
%             'Solver', 'sgd', ...
%             'ClassNames', label);
%         
%         strTitle = 'Log-Transformed Feature Space';
    end
    leftIdx = find(movieIdx(:,2)==1);
    rightIdx = find(movieIdx(:,2)==2);
    % channels x trials
    cumsumHbTM_Left = cumsumHbTM(:,leftIdx);
    cumsumHbTM_Right = cumsumHbTM(:,rightIdx);

end

% Hard coded
startChn1 = 1;
endChn1 = 16;
startChn2 = 17;
endChn2 = 30;

% Part 1
figure('units','normalized','outerposition',[0 0 1 1]); hold on;
%for i = startChn1:endChn1
if spaceOp==1||spaceOp==2
    for i = 1:numChn
        subplot(5,6,i);
        binWidth = max(max(cumsumHbTM_Left(i,:)),max(cumsumHbTM_Right(i,:)))/20;
        hold on;
        histogram(cumsumHbTM_Left(i,:),'BinWidth',binWidth','FaceAlpha',0.5,'FaceColor',[0, 0.4470, 0.7410]);
        histogram(cumsumHbTM_Right(i,:),'BinWidth',binWidth','FaceAlpha',0.5,'FaceColor',[0.8500, 0.3250, 0.0980]);
        if spaceOp == 1
            title(chnName{i});
        elseif spaceOp == 2
            title(sprintf('Component %s', num2str(i)));
        elseif spaceOp == 3
            title(chnName{i});
        end
        if i == numChn
            legend('Left','Right');
        end
        %xlim([-0.03 0.03]);
        hold off;
    end
    sgtitle(sprintf('Sbj %s : %s',sbjNum,strTitle));
    if saveOp

        fn = sprintf('Sbj %s Histogram_%s',...
            sbjNum,strTitle);
        print(gcf,[figSaveDir filesep fn],'-dpng','-r250');

    end
else
    for i = 1:numChn
        minCumsum = min(cumsumHbTM(i,:));
        maxCumsum = max(cumsumHbTM(i,:));
    
    
    
        %pData = 1/(1+exp(logRegModel.Bias+logRegModel.Beta.*
    end 
end

% % Part 2
% figure('units','normalized','outerposition',[0 0 1 1]); hold on;
% for i = startChn2:endChn2
%     subplot(4,4,i-16);
%     binWidth = max(max(cumsumHbTM_Left(i,:)),max(cumsumHbTM_Right(i,:)))/20;
%     hold on;
%     histogram(cumsumHbTM_Left,'BinWidth',binWidth','FaceAlpha',0.5,'FaceColor',[0, 0.4470, 0.7410]);
%     histogram(cumsumHbTM_Right,'BinWidth',binWidth','FaceAlpha',0.5,'FaceColor',[0.8500, 0.3250, 0.0980]);
%     %histogram(meanAlphaLeft(:,jj),'BinWidth',binWidth,'FaceAlpha',0.5,'FaceColor',[0, 0.4470, 0.7410]);
%     %hold on;
%     %histogram(meanAlphaRight(:,jj),'BinWidth',binWidth,'FaceAlpha',0.5,'FaceColor',[0.8500, 0.3250, 0.0980]);
%     %histogram(meanAlphaCenter(:,jj),'BinWidth',binWidth,'FaceAlpha',0.5,'FaceColor',[0.9290, 0.6940, 0.1250]);
%     title(chnName{i});
%     legend('Left','Right');
%     sgtitle(strTitle);
%     hold off;
% end

end

