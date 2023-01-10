% This is single-channel classification using kfolds
% timePt needs to be in sample points
% folds are partitioned outside of this function.
% Classify all channels including pruned channels.
% for preprocessFNIRS06_CV_GLM_ssBeta_SingleChn_MultiOnly


function [validationScore] = trainLDA_SingleChannelHbX_TrTst(...
    cumsumTr,cumsumTst,numClasses,movieListTrain,movieListTest)
    % OLD trials is channels x time x trial
    % trials are channels x trial (p x n)

    
    trainingResponse = movieListTrain(:,2);
    testResponse = movieListTest(:,2);
    %predResp = zeros(size(cumsumTr,1),size(cumsumTr,2));
    %validationScore = zeros(size(cumsumTr,1),1);
    
    if numClasses == 2
        label = [1;2];
    else
        label = [1;2;3];
    end
    
    classificationDiscriminant = fitcdiscr(...
        cumsumTr', ...
        trainingResponse, ...
        'DiscrimType', 'diagLinear', ...
        'Gamma', 0, ...
        'FillCoeffs', 'off', ...
        'ClassNames', label);

    foldPredictions = predict(classificationDiscriminant,cumsumTst');
    correctPredictions = squeeze(foldPredictions == testResponse);
    isMissing = isnan(testResponse);
    correctPredictions = correctPredictions(~isMissing);
    validationScore = sum(correctPredictions)/length(correctPredictions);
        
end