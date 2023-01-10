% compute the maximum likelihood estimate of global trend
% For testing
% load('intermediateOutputsCombined.mat');
% calcDWTMatrix(dc.dataTimeSeries);
function theta = calcMLETrend(Y,X)
    theta = 1;

    W = calcDWTMatrix(X,J0);
    
    
    n0 = 2^(-J0+1)*N;
    eyeN0 = eyes(n0);
    bottomLeftA = zeros(N-n0,n0);
    
    leftA = [eyeN0;bottomLeftA];
    A = [leftA W];
    
    % 
    noiseCov = computeEstimatedNoiseCov(X);
    
    
end