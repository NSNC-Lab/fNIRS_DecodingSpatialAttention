% Architecture based on de Taillez, for stimulus reconstruction
% Since we're using EEG along with audio waveform, trim EEG segments to
% only active audio part.
% Extract envelope part of audio and downsample to same fs rate as EEG.


whichCond = 1;
sbjNum = 'pk39';
numClasses = 3;
nFolds = 5;
audioDir = 'C:\Users\mn0mn\Documents\Dissertation\MovieClips';
dataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedData'];
dir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedData'];
movieListFile = [dir filesep 'iniList_pk39_20191125_.mat'];
load(movieListFile,'indexMoviesTest');
[trialsPostICA, chnsList] = partitionERPEye(sbjNum,1,0,whichCond,dataDir);
numTrials = size(trialsPostICA.allTrials,1);

inputSize = size(trialsPostICA.leftTrials,2);
% all trials numTrials x numCh x time
% want numTrials x 1 cell array
% where each element is numCh x time
xCellArr = cell(numTrials,1);

% FIRST, PLOT IT!!!
for i = 1:numTrials
    [xCellArr{i,1}] = smoothdata(squeeze(trialsPostICA.allTrials(i,:,:))',"movmedian",10)';
end

yCellArr = categorical(squeeze(indexMoviesTest(:,2)));

c = cvpartition(length(yCellArr),"KFold",nFolds);

numHiddenUnits = 10;

acc = zeros(1,5);

% layers = [ ...
%     sequenceInputLayer(inputSize)
%     %batchNormalizationLayer()
%     bilstmLayer(numHiddenUnits,'OutputMode','last')
%     fullyConnectedLayer(numClasses)
%     softmaxLayer()
%     classificationLayer()];


filterSize1 = [61 61];
numFilters1 = 3;
filterSize2 = [61 2];
numFilters2 = 1;
% fcSize1 = [246 200];
% fcSize2 = [200 200];
% fcSize3 = [200 100];
% fcSize4 = [100 1];
fcSize1 = 200;
fcSize2 = 200;
fcSize3 = 100;
fcSize4 = 1;


% Same as Fig 4 of Ciccarelli 2019 papeer
% Follow Dry EEG parameters
layers = [ ...
    sequenceInputLayer(inputSize+1)
    batchNormalizationLayer()
    convolution2dLayer(filterSize1,numFilters1)
    % alpha parameter
    eluLayer(1)
    % can specify max pooling layer to be 1d in parameter
    maxPooling2dLayer()
    convolution2dLayer(filterSize2,numFilters2)
    batchNormalizationLayer()
    fullyConnectedLayer(fcSize1)
    eluLayer(1)
    % default is 0.5 probability
    dropoutLayer()
    batchNormalizationLayer()
    fullyConnectedLayer(fcSize2)
    eluLayer(1)
    % default is 0.5 probability
    dropoutLayer()
    batchNormalizationLayer()
    fullyConnectedLayer(fcSize3)
    eluLayer(1)
    % default is 0.5 probability
    dropoutLayer()
    fullyConnectedLayer(fcSize4)];
    
    

maxEpochs = 100;
miniBatchSize = 25;



for i = 1: 1
    
    
    idxTrain = training(c,i);
    idxTest = test(c,i);
    
    options = trainingOptions('adam', ...
        'ExecutionEnvironment','cpu', ...
        'GradientThreshold',1, ...
        'MaxEpochs',maxEpochs, ...
        'MiniBatchSize',miniBatchSize, ...
        'SequenceLength','longest', ...
        'Shuffle','never', ...
        'Verbose',0, ...
        'Plots','training-progress', ...
        'ValidationFrequency',10, ...
        'L2Regularization',1.0000e-01, ...
        'ValidationData',{xCellArr(idxTest),yCellArr(idxTest)});
    
    net = trainNetwork(xCellArr(idxTrain),yCellArr(idxTrain),layers,options);
    
    YPred = classify(net,xCellArr{idxTest}, ...
        'MiniBatchSize',miniBatchSize, ...
        'SequenceLength','longest');

    acc(i) = sum(YPred == yCellArr(idxTest))./numel(yCellArr(idxTest));
end

avgAcc = mean(acc);