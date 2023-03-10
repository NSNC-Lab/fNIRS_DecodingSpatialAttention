% STATUS: active
% 
% SYNTAX:
% grpPlotPerfVsTime(goodSbj,numClasses,saveOp,rejTrOp,itrOp,errBarOp)
% 
% DESCRIPTION:
% Plot decoding performance (accuracies) as a function of time course of a
%   trial for all 3 chromophores.
% 
% RESTRICTION:
% None.
% 
% INPUTS:
% goodSbj - int: group of subjects. Encoding as follow:
%       1 - high-performance group
%       2 - low-performance group
%       0 - all subjects
% numClasses - int: number of classes to classify
%       2 - classify between left & right
%       3 - classify between left, right & center
% saveOp - int: option to save
%       0 - don't save
%       1- save
% rejTrOp - int: option to reject trials
%       0 - don't reject trials
%       1 - reject trials
% itrOp - int: option to display information transfer rate (itr)
%       0 - don't display itr.
%       1- display itr.
% errBarOp - int: option to display error bar:
%       0 - don't display error bars
%       1 - display error bars as confidence intervals
%
% RETURNED VARIABLES:
% None.
% 
% FILES SAVED:
% 1) save figure of decoding performance as a function of time course of trial
% 2) save decoding performances of all subjects in one array variable.
% 
% PLOTTING:
% Line plots of decoding performance as a function of time course of trial

function grpPlotPerfVsTime(goodSbj,numClasses,saveOp,rejTrOp,itrOp,errBarOp)

    if goodSbj == 1
        %sbjList = {'08','12','16','17','18','19'};
        % high-performing
        sbjList = {'08','12','13','16','19','21','24'};
        if rejTrOp
            if numClasses == 2
                dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
                    filesep 'sbjListPerfGood_GLM_CV_SSBeta_LR_RejTr.mat'];
            else
                dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
                    filesep 'sbjListPerfGood_GLM_CV_SSBeta_RejTr.mat'];
            end
        else
            if numClasses == 2
                dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
                    filesep 'sbjListPerfGood_GLM_CV_SSBeta_LR.mat'];
            else
                dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
                    filesep 'sbjListPerfGood_GLM_CV_SSBeta.mat'];
            end
        end
    elseif goodSbj == 0
        % low performing
        sbjList = {'14','15','22','23'};
        if rejTrOp
            if numClasses == 2
                dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
                    filesep 'sbjListPerfBad_GLM_CV_SSBeta_LR_RejTr.mat'];
            else
                dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
                    filesep 'sbjListPerfBad_GLM_CV_SSBeta_RejTr.mat'];
            end
        else
            if numClasses == 2
                dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
                    filesep 'sbjListPerfBad_GLM_CV_SSBeta_LR.mat'];
            else
                dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
                    filesep 'sbjListPerfBad_SSBeta.mat'];
            end
        end
    else
        if numClasses == 2
            % all
            sbjList = {'08','12','13','14','15','16','19','21','22','23','24','25'};
            if rejTrOp
                dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
                    filesep 'sbjListPerfAll_GLM_CV_SSBeta_LR_RejTr.mat'];
            else
                dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
                    filesep 'sbjListPerfAll_SSBeta.mat'];
            end
        else
            sbjList = {'08','12','13','14','15','16','19','21','22','23','24','25'};
            if rejTrOp
                dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
                    filesep 'sbjListPerfAll_GLM_CV_SSBeta_RejTr.mat'];
            else
                dataFN = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll' ...
                    filesep 'sbjListPerfAll_GLM_CV_SSBeta.mat'];
            end
        end
    end

    figSaveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll\Figures';

    fs = 50;
    %timePt = (0:0.25*fs:5*fs)+2*fs;
    %timePt = (0:0.25*fs:7*fs);
    timePt = -2*fs:0.25*fs:5*fs;
    %zeroT = 2*fs;
    if numClasses == 2
        numClassifiers = 6;
    else
        numClassifiers = 5;
    end
    
    P = 0:0.01:1;
    if numClasses == 2
        N = 2;
    else
        N = 3;
    end

    itf = log2(N) + P.*log2(P) + (1-P).*log2((1-P)/(N-1));
    
    if numClasses == 3
        itf(1:33) = -itf(1:33);
    end

    sbjListPerfHbO = zeros(numClassifiers,length(sbjList),length(timePt));
    sbjListPerfHbR = zeros(numClassifiers,length(sbjList),length(timePt));
    sbjListPerfHbT = zeros(numClassifiers,length(sbjList),length(timePt));
    
    grpPerfHbOAllFolds = zeros(numClassifiers,length(timePt),5*10*length(sbjList));
    grpPerfHbRAllFolds = zeros(numClassifiers,length(timePt),5*10*length(sbjList));
    grpPerfHbTAllFolds = zeros(numClassifiers,length(timePt),5*10*length(sbjList));

%     sbjListPerfSTDHbO = zeros(numClassifiers,length(sbjList),length(timePt));
%     sbjListPerfSTDHbR = zeros(numClassifiers,length(sbjList),length(timePt));
%     sbjListPerfSTDHbT = zeros(numClassifiers,length(sbjList),length(timePt));

    for i = 1:length(sbjList)
        sbjNum = sbjList{i};

        processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];

        if numClasses == 2
            if rejTrOp
                % preprocessFNIRS06_CV_GLM_ssBeta_MultiOnly
                %savePerfFN = 'performance_GLM_CV_SSBeta_LR_RejTr_SNR_1.5_OrigHomer3.mat';
                %savePerfFN = 'performance_GLM_CV_SSBeta_LR_RejTr_SNR_1.5_OrigHomer3.mat';
                %savePerfFN = 'performance_GLM_CV_SSBeta_LR_RejTr_SNR_1.5_OrigHomer3_10Hz.mat';
                savePerfFN = 'performance_GLM_CV_SSBeta_LR_RejTr_SNR_1.5_Version02.mat';
            else
                savePerfFN = 'performanceLinearDiscriminantUpdated_SSBeta_LR.mat';
            end
        else
            if rejTrOp
                %savePerfFN = 'performance_GLM_CV_SSBeta_RejTr_SNR_1.5_OrigFunc.mat';
                %savePerfFN = 'performance_GLM_CV_SSBeta_RejTr_SNR_1.5.mat';
                savePerfFN = 'performance_GLM_CV_SSBeta_RejTr_SNR_1.5_10Hz.mat';
            else
                savePerfFN = 'performanceLinearDiscriminantUpdated_SSBeta.mat';
            end
        end
        
%         load([processedDataDir filesep savePerfFN],'performanceLDALedoitHbOFold',...
%             'performanceLDALedoitHbRFold','performanceLDALedoitHbTFold',...
%             'performanceLDACERNNHbOFold','performanceLDACERNNHbRFold',...
%             'performanceLDACERNNHbTFold',...
%             'performanceLogRegHbOFold','performanceLogRegHbRFold',...
%             'performanceLogRegHbTFold',...
%             'performanceSVMHbOFold','performanceSVMHbRFold',...
%             'performanceSVMHbTFold',...
%             'performanceBaggingHbOFold','performanceBaggingHbRFold',...
%             'performanceBaggingHbTFold',...
%             'performanceBoostedHbOFold','performanceBoostedHbRFold',...
%             'performanceBoostedHbTFold');

        load([processedDataDir filesep savePerfFN],'performanceLDALedoitHbO',...
            'performanceLDALedoitHbR','performanceLDALedoitHbT',...
            'performanceLDACERNNHbO','performanceLDACERNNHbR',...
            'performanceLDACERNNHbT',...
            'performanceLogRegHbO','performanceLogRegHbR',...
            'performanceLogRegHbT',...
            'performanceSVMHbO','performanceSVMHbR',...
            'performanceSVMHbT',...
            'performanceBaggingHbO','performanceBaggingHbR',...
            'performanceBaggingHbT',...
            'performanceBoostedHbO','performanceBoostedHbR',...
            'performanceBoostedHbT');

    %     performanceLDALedoitHbOCI = zeros(length(timePt),2);
    %     performanceLDALedoitHbRCI = zeros(length(timePt),2);
    %     performanceLDALedoitHbTCI = zeros(length(timePt),2);
    % 
    %     performanceLDACERNNHbOCI = zeros(length(timePt),2);
    %     performanceLDACERNNHbRCI = zeros(length(timePt),2);
    %     performanceLDACERNNHbTCI = zeros(length(timePt),2);
    % 
    %     performanceLogRegHbOCI = zeros(length(timePt),2);
    %     performanceLogRegHbRCI = zeros(length(timePt),2);
    %     performanceLogRegHbTCI = zeros(length(timePt),2);
    % 
    %     performanceSVMHbOCI = zeros(length(timePt),2);
    %     performanceSVMHbRCI = zeros(length(timePt),2);
    %     performanceSVMHbTCI = zeros(length(timePt),2);
    % 
    %     performanceBaggingHbOCI = zeros(length(timePt),2);
    %     performanceBaggingHbRCI = zeros(length(timePt),2);
    %     performanceBaggingHbTCI = zeros(length(timePt),2);
    % 
    %     performanceBoostedHbOCI = zeros(length(timePt),2);
    %     performanceBoostedHbRCI = zeros(length(timePt),2);
    %     performanceBoostedHbTCI = zeros(length(timePt),2);

% I think old format is no longer needed but keep just in case
%         performanceLDALedoitHbOFoldRS = reshape(performanceLDALedoitHbOFold,...
%             [size(performanceLDALedoitHbOFold,1), size(performanceLDALedoitHbOFold,2)*size(performanceLDALedoitHbOFold,3)]);
%         performanceLDALedoitHbRFoldRS = reshape(performanceLDALedoitHbRFold,...
%             [size(performanceLDALedoitHbRFold,1), size(performanceLDALedoitHbRFold,2)*size(performanceLDALedoitHbRFold,3)]);
%         performanceLDALedoitHbTFoldRS = reshape(performanceLDALedoitHbTFold,...
%             [size(performanceLDALedoitHbTFold,1), size(performanceLDALedoitHbTFold,2)*size(performanceLDALedoitHbTFold,3)]);
% 
%         performanceLDACERNNHbOFoldRS = reshape(performanceLDACERNNHbOFold,...
%             [size(performanceLDACERNNHbOFold,1), size(performanceLDACERNNHbOFold,2)*size(performanceLDACERNNHbOFold,3)]);
%         performanceLDACERNNHbRFoldRS = reshape(performanceLDACERNNHbRFold,...
%             [size(performanceLDACERNNHbRFold,1), size(performanceLDACERNNHbRFold,2)*size(performanceLDACERNNHbRFold,3)]);
%         performanceLDACERNNHbTFoldRS = reshape(performanceLDACERNNHbTFold,...
%             [size(performanceLDACERNNHbTFold,1), size(performanceLDACERNNHbTFold,2)*size(performanceLDACERNNHbTFold,3)]);
% 
%         performanceLogRegHbOFoldRS = reshape(performanceLogRegHbOFold,...
%             [size(performanceLogRegHbOFold,1), size(performanceLogRegHbOFold,2)*size(performanceLogRegHbOFold,3)]);
%         performanceLogRegHbRFoldRS = reshape(performanceLogRegHbRFold,...
%             [size(performanceLogRegHbRFold,1), size(performanceLogRegHbRFold,2)*size(performanceLogRegHbRFold,3)]);
%         performanceLogRegHbTFoldRS = reshape(performanceLogRegHbTFold,...
%             [size(performanceLogRegHbTFold,1), size(performanceLogRegHbTFold,2)*size(performanceLogRegHbTFold,3)]);
% 
%         performanceSVMHbOFoldRS = reshape(performanceSVMHbOFold,...
%             [size(performanceSVMHbOFold,1), size(performanceSVMHbOFold,2)*size(performanceSVMHbOFold,3)]);
%         performanceSVMHbRFoldRS = reshape(performanceSVMHbRFold,...
%             [size(performanceSVMHbRFold,1), size(performanceSVMHbRFold,2)*size(performanceSVMHbRFold,3)]);
%         performanceSVMHbTFoldRS = reshape(performanceSVMHbTFold,...
%             [size(performanceSVMHbTFold,1), size(performanceSVMHbTFold,2)*size(performanceSVMHbTFold,3)]);
% 
%         performanceBaggingHbOFoldRS = reshape(performanceBaggingHbOFold,...
%             [size(performanceBaggingHbOFold,1), size(performanceBaggingHbOFold,2)*size(performanceBaggingHbOFold,3)]);
%         performanceBaggingHbRFoldRS = reshape(performanceBaggingHbRFold,...
%             [size(performanceBaggingHbRFold,1), size(performanceBaggingHbRFold,2)*size(performanceBaggingHbRFold,3)]);
%         performanceBaggingHbTFoldRS = reshape(performanceBaggingHbTFold,...
%             [size(performanceBaggingHbTFold,1), size(performanceBaggingHbTFold,2)*size(performanceBaggingHbTFold,3)]);
% 
%         performanceBoostedHbOFoldRS = reshape(performanceBoostedHbOFold,...
%             [size(performanceBoostedHbOFold,1), size(performanceBoostedHbOFold,2)*size(performanceBoostedHbOFold,3)]);
%         performanceBoostedHbRFoldRS = reshape(performanceBoostedHbRFold,...
%             [size(performanceBoostedHbRFold,1), size(performanceBoostedHbRFold,2)*size(performanceBoostedHbRFold,3)]);
%         performanceBoostedHbTFoldRS = reshape(performanceBoostedHbTFold,...
%             [size(performanceBoostedHbTFold,1), size(performanceBoostedHbTFold,2)*size(performanceBoostedHbTFold,3)]);

        % sanity check. Passed!
    %     temp1 = performanceLDALedoitHbOFold(1,:,:);
    %     test1 = mean(temp1(:));
    %     temp2 = performanceLDALedoitHbOFold(2,:,:);
    %     test2 = mean(temp2(:));
    %     temp3 = performanceLDALedoitHbOFold(3,:,:);
    %     test3 = mean(temp3(:));
%         sbjListPerfHbO(1,i,:) = mean(performanceLDALedoitHbOFoldRS,2);
%         sbjListPerfHbR(1,i,:) = mean(performanceLDALedoitHbRFoldRS,2);
%         sbjListPerfHbT(1,i,:) = mean(performanceLDALedoitHbTFoldRS,2);
% 
%         sbjListPerfHbO(2,i,:) = mean(performanceLDACERNNHbOFoldRS,2);
%         sbjListPerfHbR(2,i,:) = mean(performanceLDACERNNHbRFoldRS,2);
%         sbjListPerfHbT(2,i,:) = mean(performanceLDACERNNHbTFoldRS,2);
% 
%         sbjListPerfHbO(3,i,:) = mean(performanceLogRegHbOFoldRS,2);
%         sbjListPerfHbR(3,i,:) = mean(performanceLogRegHbRFoldRS,2);
%         sbjListPerfHbT(3,i,:) = mean(performanceLogRegHbTFoldRS,2);
% 
%         sbjListPerfHbO(4,i,:) = mean(performanceSVMHbOFoldRS,2);
%         sbjListPerfHbR(4,i,:) = mean(performanceSVMHbRFoldRS,2);
%         sbjListPerfHbT(4,i,:) = mean(performanceSVMHbTFoldRS,2);
% 
%         sbjListPerfHbO(5,i,:) = mean(performanceBaggingHbOFoldRS,2);
%         sbjListPerfHbR(5,i,:) = mean(performanceBaggingHbRFoldRS,2);
%         sbjListPerfHbT(5,i,:) = mean(performanceBaggingHbTFoldRS,2);
% 
%         sbjListPerfHbO(6,i,:) = mean(performanceBoostedHbOFoldRS,2);
%         sbjListPerfHbR(6,i,:) = mean(performanceBoostedHbRFoldRS,2);
%         sbjListPerfHbT(6,i,:) = mean(performanceBoostedHbTFoldRS,2);

        sbjListPerfHbO(1,i,:) = mean(performanceLDALedoitHbO,1);
        sbjListPerfHbR(1,i,:) = mean(performanceLDALedoitHbR,1);
        sbjListPerfHbT(1,i,:) = mean(performanceLDALedoitHbT,1);

        sbjListPerfHbO(2,i,:) = mean(performanceLDACERNNHbO,1);
        sbjListPerfHbR(2,i,:) = mean(performanceLDACERNNHbR,1);
        sbjListPerfHbT(2,i,:) = mean(performanceLDACERNNHbT,1);

        sbjListPerfHbO(3,i,:) = mean(performanceLogRegHbO,1);
        sbjListPerfHbR(3,i,:) = mean(performanceLogRegHbR,1);
        sbjListPerfHbT(3,i,:) = mean(performanceLogRegHbT,1);

        sbjListPerfHbO(4,i,:) = mean(performanceSVMHbO,1);
        sbjListPerfHbR(4,i,:) = mean(performanceSVMHbR,1);
        sbjListPerfHbT(4,i,:) = mean(performanceSVMHbT,1);

        sbjListPerfHbO(5,i,:) = mean(performanceBaggingHbO,1);
        sbjListPerfHbR(5,i,:) = mean(performanceBaggingHbR,1);
        sbjListPerfHbT(5,i,:) = mean(performanceBaggingHbT,1);

        sbjListPerfHbO(6,i,:) = mean(performanceBoostedHbO,1);
        sbjListPerfHbR(6,i,:) = mean(performanceBoostedHbR,1);
        sbjListPerfHbT(6,i,:) = mean(performanceBoostedHbT,1);
        
        startIdx = (i-1)*50+1;
        endIdx = i*50;
        
%         grpPerfHbOAllFolds(1,:,startIdx:endIdx) = performanceLDALedoitHbOFoldRS;
%         grpPerfHbRAllFolds(1,:,startIdx:endIdx) = performanceLDALedoitHbRFoldRS;
%         grpPerfHbTAllFolds(1,:,startIdx:endIdx) = performanceLDALedoitHbTFoldRS;
%         
%         grpPerfHbOAllFolds(2,:,startIdx:endIdx) = performanceLDACERNNHbOFoldRS;
%         grpPerfHbRAllFolds(2,:,startIdx:endIdx) = performanceLDACERNNHbRFoldRS;
%         grpPerfHbTAllFolds(2,:,startIdx:endIdx) = performanceLDACERNNHbTFoldRS;
%         
%         grpPerfHbOAllFolds(3,:,startIdx:endIdx) = performanceLogRegHbOFoldRS;
%         grpPerfHbRAllFolds(3,:,startIdx:endIdx) = performanceLogRegHbRFoldRS;
%         grpPerfHbTAllFolds(3,:,startIdx:endIdx) = performanceLogRegHbTFoldRS;
%         
%         grpPerfHbOAllFolds(4,:,startIdx:endIdx) = performanceSVMHbOFoldRS;
%         grpPerfHbRAllFolds(4,:,startIdx:endIdx) = performanceSVMHbRFoldRS;
%         grpPerfHbTAllFolds(4,:,startIdx:endIdx) = performanceSVMHbTFoldRS;
%         
%         grpPerfHbOAllFolds(5,:,startIdx:endIdx) = performanceBaggingHbOFoldRS;
%         grpPerfHbRAllFolds(5,:,startIdx:endIdx) = performanceBaggingHbRFoldRS;
%         grpPerfHbTAllFolds(5,:,startIdx:endIdx) = performanceBaggingHbTFoldRS;
%         
%         grpPerfHbOAllFolds(6,:,startIdx:endIdx) = performanceBoostedHbOFoldRS;
%         grpPerfHbRAllFolds(6,:,startIdx:endIdx) = performanceBoostedHbRFoldRS;
%         grpPerfHbTAllFolds(6,:,startIdx:endIdx) = performanceBoostedHbTFoldRS;

        grpPerfHbOAllFolds(1,:,startIdx:endIdx) = performanceLDALedoitHbO';
        grpPerfHbRAllFolds(1,:,startIdx:endIdx) = performanceLDALedoitHbR';
        grpPerfHbTAllFolds(1,:,startIdx:endIdx) = performanceLDALedoitHbT';
        
        grpPerfHbOAllFolds(2,:,startIdx:endIdx) = performanceLDACERNNHbO';
        grpPerfHbRAllFolds(2,:,startIdx:endIdx) = performanceLDACERNNHbR';
        grpPerfHbTAllFolds(2,:,startIdx:endIdx) = performanceLDACERNNHbT';
        
        grpPerfHbOAllFolds(3,:,startIdx:endIdx) = performanceLogRegHbO';
        grpPerfHbRAllFolds(3,:,startIdx:endIdx) = performanceLogRegHbR';
        grpPerfHbTAllFolds(3,:,startIdx:endIdx) = performanceLogRegHbT';
        
        grpPerfHbOAllFolds(4,:,startIdx:endIdx) = performanceSVMHbO';
        grpPerfHbRAllFolds(4,:,startIdx:endIdx) = performanceSVMHbR';
        grpPerfHbTAllFolds(4,:,startIdx:endIdx) = performanceSVMHbT';
        
        grpPerfHbOAllFolds(5,:,startIdx:endIdx) = performanceBaggingHbO';
        grpPerfHbRAllFolds(5,:,startIdx:endIdx) = performanceBaggingHbR';
        grpPerfHbTAllFolds(5,:,startIdx:endIdx) = performanceBaggingHbT';
        
        grpPerfHbOAllFolds(6,:,startIdx:endIdx) = performanceBoostedHbO';
        grpPerfHbRAllFolds(6,:,startIdx:endIdx) = performanceBoostedHbR';
        grpPerfHbTAllFolds(6,:,startIdx:endIdx) = performanceBoostedHbT';

    %     %STD
    %     sbjListPerfSTDHbO(1,i,:) = std(performanceLDALedoitHbOFoldRS,[],2);
    %     sbjListPerfSTDHbR(1,i,:) = std(performanceLDALedoitHbRFoldRS,[],2);
    %     sbjListPerfSTDHbT(1,i,:) = std(performanceLDALedoitHbTFoldRS,[],2);
    %     
    %     sbjListPerfSTDHbO(2,i,:) = std(performanceLDACERNNHbOFoldRS,[],2);
    %     sbjListPerfSTDHbR(2,i,:) = std(performanceLDACERNNHbRFoldRS,[],2);
    %     sbjListPerfSTDHbT(2,i,:) = std(performanceLDACERNNHbTFoldRS,[],2);
    %     
    %     sbjListPerfSTDHbO(3,i,:) = std(performanceLogRegHbOFoldRS,[],2);
    %     sbjListPerfSTDHbR(3,i,:) = std(performanceLogRegHbRFoldRS,[],2);
    %     sbjListPerfSTDHbT(3,i,:) = std(performanceLogRegHbTFoldRS,[],2);
    %     
    %     sbjListPerfSTDHbO(4,i,:) = std(performanceSVMHbOFoldRS,[],2);
    %     sbjListPerfSTDHbR(4,i,:) = std(performanceSVMHbRFoldRS,[],2);
    %     sbjListPerfSTDHbT(4,i,:) = std(performanceSVMHbTFoldRS,[],2);
    %     
    %     sbjListPerfSTDHbO(5,i,:) = std(performanceBaggingHbOFoldRS,[],2);
    %     sbjListPerfSTDHbR(5,i,:) = std(performanceBaggingHbRFoldRS,[],2);
    %     sbjListPerfSTDHbT(5,i,:) = std(performanceBaggingHbTFoldRS,[],2);
    %     
    %     sbjListPerfSTDHbO(6,i,:) = std(performanceBoostedHbOFoldRS,[],2);
    %     sbjListPerfSTDHbR(6,i,:) = std(performanceBoostedHbRFoldRS,[],2);
    %     sbjListPerfSTDHbT(6,i,:) = std(performanceBoostedHbTFoldRS,[],2);

    end
    
    save(dataFN,'sbjListPerfHbO','sbjListPerfHbR','sbjListPerfHbT');

    grpPerfHbO = squeeze(mean(sbjListPerfHbO,2));
    grpPerfHbR = squeeze(mean(sbjListPerfHbR,2));
    grpPerfHbT = squeeze(mean(sbjListPerfHbT,2));

%     grpPerfSTDHbO = squeeze(std(sbjListPerfHbO,[],2));
%     grpPerfSTDHbR = squeeze(std(sbjListPerfHbR,[],2));
%     grpPerfSTDHbT = squeeze(std(sbjListPerfHbT,[],2));
    
    % grpPerfHbOAllFolds is 6 x 21 x 200
    
    grpPerfHbOAllFoldsSE = squeeze(std(grpPerfHbOAllFolds,[],3))./sqrt(size(grpPerfHbOAllFolds,3));
    grpPerfHbRAllFoldsSE = squeeze(std(grpPerfHbRAllFolds,[],3))./sqrt(size(grpPerfHbOAllFolds,3));
    grpPerfHbTAllFoldsSE = squeeze(std(grpPerfHbTAllFolds,[],3))./sqrt(size(grpPerfHbOAllFolds,3));

    conf = 0.95;
    
    nSz = size(grpPerfHbOAllFolds,3)*ones(1,6);
    
    grpPerfCIHbO = zeros(6,2,length(timePt));
    grpPerfCIHbR = zeros(6,2,length(timePt));
    grpPerfCIHbT = zeros(6,2,length(timePt));
    
    for i2 = 1:length(timePt)
    
        grpPerfCIHbO(:,:,i2) = calcCITDist(nSz,grpPerfHbOAllFoldsSE(:,i2),conf);
        grpPerfCIHbR(:,:,i2) = calcCITDist(nSz,grpPerfHbRAllFoldsSE(:,i2),conf);
        grpPerfCIHbT(:,:,i2) = calcCITDist(nSz,grpPerfHbTAllFoldsSE(:,i2),conf);
    
    end
    
%     grpPerfCIHbONeg = grpPerfHbO-squeeze(grpPerfCIHbO(:,1,:));
%     grpPerfCIHbOPos = -grpPerfHbO+squeeze(grpPerfCIHbO(:,2,:));
%     
%     grpPerfCIHbRNeg = grpPerfHbR-squeeze(grpPerfCIHbR(:,1,:));
%     grpPerfCIHbRPos = -grpPerfHbR+squeeze(grpPerfCIHbR(:,2,:));
%     
%     grpPerfCIHbTNeg = grpPerfHbT-squeeze(grpPerfCIHbT(:,1,:));
%     grpPerfCIHbTPos = -grpPerfHbT+squeeze(grpPerfCIHbT(:,2,:));
    
    %cmap = jet(numClassifiers);
    cmap = loadDefaultColors(1);
    yLimAxis = [0.2 1];
    
    figure('units','normalized','outerposition',[0 0 0.5 0.6]); hold on;
    %figure();

    subplot(1,3,1);hold on;
    % for j = 1:size(performanceArrMultiHbO,1)
    %     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
    %     plot(timePt./fs-2,performanceArrMultiHbO(j,:),'Color',cmap(j,:));
    % end
    % plot(timePt./fs-2,performanceArrMultiHbO_95,'Color',cmap(1,:));
    % plot(timePt./fs-2,performanceArrMultiHbO_85,'Color',cmap(2,:));
    % plot(timePt./fs-2,performanceArrMultiHbO_75,'Color',cmap(3,:));
    % plot(timePt./fs-2,performanceArrMultiHbO_65,'Color',cmap(4,:));
    % plot(timePt./fs-2,performanceArrMultiHbO_55,'Color',cmap(5,:));
    % plot(timePt./fs-2,performanceLDALedoitHbO,'Color',cmap(6,:));
    % plot(timePt./fs-2,performanceLDACERNNHbO,'Color',cmap(7,:));
    % % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
    % % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
    % plot(timePt./fs-2,performanceLogRegHbO,'Color',cmap(8,:));
    % plot(timePt./fs-2,performanceSVMHbO,'Color',cmap(9,:));

%     errorbar(timePt./fs-2,performanceLDALedoitHbOMean,performanceLDALedoitHbOCI(:,2),'Color',cmap(1,:));
%     errorbar(timePt./fs-2,performanceLDACERNNHbOMean,performanceLDACERNNHbOCI(:,2),'Color',cmap(2,:));
%     % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
%     % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
%     errorbar(timePt./fs-2,performanceLogRegHbOMean,performanceLogRegHbOCI(:,2),'Color',cmap(3,:));
%     errorbar(timePt./fs-2,performanceSVMHbOMean,performanceSVMHbOCI(:,2),'Color',cmap(4,:));
%     errorbar(timePt./fs-2,performanceBaggingHbOMean,performanceBaggingHbOCI(:,2),'Color',cmap(5,:));
%     errorbar(timePt./fs-2,performanceBoostedHbOMean,performanceBoostedHbOCI(:,2),'Color',cmap(6,:));

    for i = 1:numClassifiers
        %errorbar(timePt./fs-2,grpPerfHbO(i,:),grpPerfSTDHbO(i,:),'Color',cmap(i,:));
        %errorbar(timePt./fs-2,grpPerfHbO(i,:),grpPerfCIHbO(i,:,:),'Color',cmap(i,:));
        %errorbar(timePt./fs-2,grpPerfHbO(i,:),abs(squeeze(grpPerfCIHbO(i,1,:))),abs(squeeze(grpPerfCIHbO(i,2,:))),'Color',cmap(i,:));
        if errBarOp
            errorbar(timePt./fs,grpPerfHbO(i,:),abs(squeeze(grpPerfCIHbO(i,1,:))),abs(squeeze(grpPerfCIHbO(i,2,:))),'Color',cmap(i,:));
        else
            plot(timePt./fs,grpPerfHbO(i,:),'Color',cmap(i,:));
        end
    end

    if itrOp
        yyaxis left;
        ylim(yLimAxis);
        ylabel('Accuracy');
        yyaxis right;
        if numClasses == 2
            ylim([-max(itf) max(itf)]);
            yticks([-max(itf) 0 max(itf)]);
            yticklabels({num2str(max(itf)) num2str(0) num2str(max(itf))});
        else
            ylim([min(itf) max(itf)]);
            yticks([min(itf) 0 max(itf)]);
            yticks(linspace(min(itf), max(itf),6));
            tempIdx = linspace(1,100,6);
            yticklabels([-itf(2) -itf(round(tempIdx(2))) ...
                itf(round(tempIdx(3))) itf(round(tempIdx(4))) itf(round(tempIdx(5))) ...
                itf(round(tempIdx(6)))]);
            %yticklabels({num2str(max(itf)) num2str(0) num2str(max(itf))});  
        end
        ylabel('ITR bits/trial');
    else
        ylabel('Accuracy');
        ylim(yLimAxis);
    end
    if numClasses == 2
        yline(.5,'--');
    else
        yline(0.333,'--');
    end
    xlim([-2 7]);
    title(sprintf('Group Avg: ??[HbO] Multi'));
    if numClasses == 2
        legend({'LDAShrink','LDA CERNN','Cosine KNN','Cubic-SVM','Random Forest','Boosted'},'Location','southeast');
    else
        legend({'LDAShrink','LDA CERNN','Cosine KNN','Cubic-SVM','Random Forest'},'Location','southeast');
    end
    xlabel('Time [s]');
    hold off;

    subplot(1,3,2);hold on;
    % for j = 1:size(performanceArrMultiHbR,1)
    %     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
    %     plot(timePt./fs-2,performanceArrMultiHbR(j,:),'Color',cmap(j,:));
    % end
    % plot(timePt./fs-2,performanceArrMultiHbR_95,'Color',cmap(1,:));
    % plot(timePt./fs-2,performanceArrMultiHbR_85,'Color',cmap(2,:));
    % plot(timePt./fs-2,performanceArrMultiHbR_75,'Color',cmap(3,:));
    % plot(timePt./fs-2,performanceArrMultiHbR_65,'Color',cmap(4,:));
    % plot(timePt./fs-2,performanceArrMultiHbR_55,'Color',cmap(5,:));
    % plot(timePt./fs-2,performanceLDALedoitHbR,'Color',cmap(6,:));
    % plot(timePt./fs-2,performanceLDACERNNHbR,'Color',cmap(7,:));
    % % plot(timePt./fs-2,performanceQDALedoitHbR,'Color',cmap(7,:));
    % % plot(timePt./fs-2,performanceLDAGDHbR,'Color',cmap(8,:));
    % plot(timePt./fs-2,performanceLogRegHbR,'Color',cmap(8,:));
    % plot(timePt./fs-2,performanceSVMHbR,'Color',cmap(9,:));

%     errorbar(timePt./fs-2,performanceLDALedoitHbRMean,performanceLDALedoitHbRCI(:,2),'Color',cmap(1,:));
%     errorbar(timePt./fs-2,performanceLDACERNNHbRMean,performanceLDACERNNHbRCI(:,2),'Color',cmap(2,:));
%     % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
%     % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
%     errorbar(timePt./fs-2,performanceLogRegHbRMean,performanceLogRegHbRCI(:,2),'Color',cmap(3,:));
%     errorbar(timePt./fs-2,performanceSVMHbRMean,performanceSVMHbRCI(:,2),'Color',cmap(4,:));
%     errorbar(timePt./fs-2,performanceBaggingHbRMean,performanceBaggingHbRCI(:,2),'Color',cmap(5,:));
%     errorbar(timePt./fs-2,performanceBoostedHbRMean,performanceBoostedHbRCI(:,2),'Color',cmap(6,:));

    for i = 1:numClassifiers
        %errorbar(timePt./fs-2,grpPerfHbR(i,:),grpPerfSTDHbR(i,:),'Color',cmap(i,:));
        %errorbar(timePt./fs-2,grpPerfHbR(i,:),abs(squeeze(grpPerfCIHbR(i,1,:))),abs(squeeze(grpPerfCIHbR(i,2,:))),'Color',cmap(i,:));
        if errBarOp
            errorbar(timePt./fs,grpPerfHbR(i,:),abs(squeeze(grpPerfCIHbO(i,1,:))),abs(squeeze(grpPerfCIHbO(i,2,:))),'Color',cmap(i,:));
        else
            plot(timePt./fs,grpPerfHbR(i,:),'Color',cmap(i,:));
        end
    end

    if itrOp
        yyaxis left;
        ylim(yLimAxis);
        ylabel('Accuracy');
        yyaxis right;
        if numClasses == 2
            ylim([-max(itf) max(itf)]);
            yticks([-max(itf) 0 max(itf)]);
            yticklabels({num2str(max(itf)) num2str(0) num2str(max(itf))});
        else
            ylim([min(itf) max(itf)]);
            yticks([min(itf) 0 max(itf)]);
            yticks(linspace(min(itf), max(itf),6));
            tempIdx = linspace(1,100,6);
            yticklabels([-itf(2) -itf(round(tempIdx(2))) ...
                itf(round(tempIdx(3))) itf(round(tempIdx(4))) itf(round(tempIdx(5))) ...
                itf(round(tempIdx(6)))]);
            %yticklabels({num2str(max(itf)) num2str(0) num2str(max(itf))});  
        end
        ylabel('ITR bits/trial');
    else
        ylabel('Accuracy');
        ylim(yLimAxis);
    end
    xlim([-2 7]);
    if numClasses == 2
        yline(.5,'--');
    else
        yline(0.333,'--');
    end
    title(sprintf('Group Avg: ??[HbR] Multi'));
    %legend('S1D1','S1D2','S1D3','S1D4');
    xlabel('Time [s]');
    hold off;

    subplot(1,3,3);hold on;
    % for j = 1:size(performanceArrMultiHbT,1)
    %     %plot(timePt./fs-2,performanceArrMulti(sizeChn*(j-1) + 1,:));
    %     plot(timePt./fs-2,performanceArrMultiHbT(j,:),'Color',cmap(j,:));
    % end
    % plot(timePt./fs-2,performanceArrMultiHbT_95,'Color',cmap(1,:));
    % plot(timePt./fs-2,performanceArrMultiHbT_85,'Color',cmap(2,:));
    % plot(timePt./fs-2,performanceArrMultiHbT_75,'Color',cmap(3,:));
    % plot(timePt./fs-2,performanceArrMultiHbT_65,'Color',cmap(4,:));
    % plot(timePt./fs-2,performanceArrMultiHbT_55,'Color',cmap(5,:));
    % plot(timePt./fs-2,performanceLDALedoitHbT,'Color',cmap(6,:));
    % plot(timePt./fs-2,performanceLDACERNNHbT,'Color',cmap(7,:));
    % % plot(timePt./fs-2,performanceQDALedoitHbT,'Color',cmap(7,:));
    % % plot(timePt./fs-2,performanceLDAGDHbT,'Color',cmap(8,:));
    % plot(timePt./fs-2,performanceLogRegHbT,'Color',cmap(8,:));
    % plot(timePt./fs-2,performanceSVMHbT,'Color',cmap(9,:));

%     errorbar(timePt./fs-2,performanceLDALedoitHbTMean,performanceLDALedoitHbTCI(:,2),'Color',cmap(1,:));
%     errorbar(timePt./fs-2,performanceLDACERNNHbTMean,performanceLDACERNNHbTCI(:,2),'Color',cmap(2,:));
%     % plot(timePt./fs-2,performanceQDALedoitHbO,'Color',cmap(7,:));
%     % plot(timePt./fs-2,performanceLDAGDHbO,'Color',cmap(8,:));
%     errorbar(timePt./fs-2,performanceLogRegHbTMean,performanceLogRegHbTCI(:,2),'Color',cmap(3,:));
%     errorbar(timePt./fs-2,performanceSVMHbTMean,performanceSVMHbTCI(:,2),'Color',cmap(4,:));
%     errorbar(timePt./fs-2,performanceBaggingHbTMean,performanceBaggingHbTCI(:,2),'Color',cmap(5,:));
%     errorbar(timePt./fs-2,performanceBoostedHbTMean,performanceBoostedHbTCI(:,2),'Color',cmap(6,:));

    for i = 1:numClassifiers
        %errorbar(timePt./fs-2,grpPerfHbT(i,:),grpPerfSTDHbT(i,:),'Color',cmap(i,:));
        %errorbar(timePt./fs-2,grpPerfHbT(i,:),abs(squeeze(grpPerfCIHbT(i,1,:))),abs(squeeze(grpPerfCIHbT(i,2,:))),'Color',cmap(i,:));
        if errBarOp
            errorbar(timePt./fs,grpPerfHbT(i,:),abs(squeeze(grpPerfCIHbO(i,1,:))),abs(squeeze(grpPerfCIHbO(i,2,:))),'Color',cmap(i,:));
        else
            plot(timePt./fs,grpPerfHbT(i,:),'Color',cmap(i,:));
        end
    end

    if itrOp
        yyaxis left;
        ylim(yLimAxis);
        ylabel('Accuracy');
        yyaxis right;
        if numClasses == 2
            ylim([-max(itf) max(itf)]);
            yticks([-max(itf) 0 max(itf)]);
            yticklabels({num2str(max(itf)) num2str(0) num2str(max(itf))});
        else
            ylim([min(itf) max(itf)]);
            yticks([min(itf) 0 max(itf)]);
            yticks(linspace(min(itf), max(itf),6));
            tempIdx = linspace(1,100,6);
            yticklabels([-itf(2) -itf(round(tempIdx(2))) ...
                itf(round(tempIdx(3))) itf(round(tempIdx(4))) itf(round(tempIdx(5))) ...
                itf(round(tempIdx(6)))]);
            %yticklabels({num2str(max(itf)) num2str(0) num2str(max(itf))});  
        end
        ylabel('ITR bits/trial');
    else
        ylabel('Accuracy');
        ylim(yLimAxis);
    end
    xlim([-2 7]);
    if numClasses == 2
        yline(.5,'--');
    else
        yline(0.333,'--');
    end
    title(sprintf('Group Avg: ??[HbT] Multi'));
    %legend('S1D1','S1D2','S1D3','S1D4');
    xlabel('Time [s]');
    hold off;


    if goodSbj==1 && numClasses == 2
        if rejTrOp
            fn = sprintf('GrpPerformanceCumsumVsTime_DiffStartT_LR_AllChns_DiffClassifiers_GoodSbjs_RejTr');
        else
            fn = sprintf('GrpPerformanceCumsumVsTime_DiffStartT_LR_AllChns_DiffClassifiers_GoodSbjs');
        end
        %fn = sprintf('GrpPerformanceCumsumVsTime_DiffStartT_LR_AllChns_DiffClassifiers_GoodSbjs_SSBeta_NoChncLvl');
    elseif goodSbj==1 && numClasses == 3
        if rejTrOp
            fn = sprintf('GrpPerformanceCumsumVsTime_DiffStartT_AllChns_DiffClassifiers_GoodSbjs_SSBeta_RejTr');
        %fn = sprintf('GrpPerformanceCumsumVsTime_DiffStartT_AllChns_DiffClassifiers_GoodSbjs');
        else
            fn = sprintf('GrpPerformanceCumsumVsTime_DiffStartT_AllChns_DiffClassifiers_GoodSbjs_SSBeta');
        end
    elseif goodSbj==0 && numClasses == 2
        if rejTrOp
        %fn = sprintf('GrpPerformanceCumsumVsTime_DiffStartT_LR_AllChns_DiffClassifiers_BadSbjs');
            fn = sprintf('GrpPerformanceCumsumVsTime_DiffStartT_LR_AllChns_DiffClassifiers_BadSbjs_SSBeta_RejTr');
        else
            fn = sprintf('GrpPerformanceCumsumVsTime_DiffStartT_LR_AllChns_DiffClassifiers_BadSbjs_SSBeta');
        end
    elseif goodSbj==0 && numClasses == 3
        if rejTrOp
        %fn = sprintf('GrpPerformanceCumsumVsTime_DiffStartT_AllChns_DiffClassifiers_BadSbjs');
            fn = sprintf('GrpPerformanceCumsumVsTime_DiffStartT_AllChns_DiffClassifiers_BadSbjs_SSBeta_RejTr');
        else
            fn = sprintf('GrpPerformanceCumsumVsTime_DiffStartT_AllChns_DiffClassifiers_BadSbjs_SSBeta');
        end
    else
        if rejTrOp
            fn = sprintf('GrpPerformanceCumsumVsTime_DiffStartT_AllChns_DiffClassifiers_AllSbjs_SSBeta_RejTr');
        else
        %fn = sprintf('GrpPerformanceCumsumVsTime_DiffStartT_AllChns_DiffClassifiers_AllSbjs');
            fn = sprintf('GrpPerformanceCumsumVsTime_DiffStartT_AllChns_DiffClassifiers_AllSbjs_SSBeta');
        end
    end

    if saveOp == 1
        print(gcf,[figSaveDir filesep fn],'-dpng','-r250');
    end

end