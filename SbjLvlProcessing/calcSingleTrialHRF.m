% Calculate and save single-trial HRF bc it's taking forever
saveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment03';
load([saveDir filesep 'intermediateOutputs.mat']);

rejectedChn = find(mlActAuto{1}==0);
% need to save old stim and compare with new stim var
%rejectedTrials = find(

probe = snirf1.probe;
Aaux = [];
rcMap = [];

singleTrialHRFCell = cell(6,1);

% get dc var from hmrR_OD2Conc func and split into ind trials
for i = 1:length(stim)
    % channels x time x trial
    thisSingleTrialHRF = zeros(size(dcAvg.dataTimeSeries,2)/6,size(dcAvg.time,1),length(squeeze(stim(1,i).data(:,1))));
    disp(['i: ', num2str(i)]);
    for j = 1:length(squeeze(stim(1,i).data(:,1)))
        disp(['j: ',num2str(j)]);
        % Create dataclass object. Must have dataTimeSeries, time, and
        % measurement properties.
        %thisDC = DataClass(dc.dataTimeSeries(stim(1,i).data(j,1)-2:stim(1,i).data(j,1)+15,:),dc.time(stim(1,i).data(j,1)-2:stim(1,i).data(j,1)+15,1),dc.measurementList);
        
        % Create StimClass object
        thisStim = StimClass('thisStim');
        AddStims(thisStim, stim(1,i).data(j,1));
        updateStates(thisStim);
        
        %k = 
        %thisMLActAuto = {mlActAuto{1}(k)};
        %thisMLActAuto = {[1]};
        
        %[dcAvg,dcAvgStd,nTrials,dcNew,dcResid,dcSum2,beta,R,hmrstats] = hmrR_GLM(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,[-2  15],1,1,[1  1  0  0  0  0],15,1,3,0);
        [thisdcAvg,thisdcAvgStd,thisnTrials,thisdcNew,thisdcResid,thisdcSum2,thisbeta,thisR,thishmrstats] = hmrR_GLM(dc,thisStim,probe,mlActAuto,Aaux,tIncAuto,rcMap,[-2  15],1,1,[1  1  0  0  0  0],15,1,3,0);
        
        if (isempty(thisdcAvg.dataTimeSeries))
            thisSingleTrialHRF(:,:,j) = zeros(size(dcAvg.dataTimeSeries,2)/6,size(dcAvg.time,1));
        else
            thisSingleTrialHRF(:,:,j) = thisdcAvg.dataTimeSeries';
        end
        
    end
    singleTrialHRFCell{i} = thisSingleTrialHRF;
end

saveDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment03';
save([saveDir filesep 'singleTrials.mat'],'singleTrialHRFCell');