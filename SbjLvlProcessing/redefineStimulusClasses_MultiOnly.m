% redefine stimulus classes into 6 diff. classes
% to be runned after generating stimuli marks from aux data in Homer3 GUI
% Combined multiple snirf files into one.
% strictly for ExperimentXXGUI folder
% Only 1 snirf file

function redefineStimulusClasses_MultiOnly(sbjNum,rawDataFN,respData)

dir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum) 'GUI'];
load([dir filesep respData '.mat'],'fixedMaskerList','indexMoviesTest','maskerMovies','numTrials','uniqueMovies');

snirf2 = SnirfLoad([dir filesep rawDataFN{1} '.snirf']);

oldStimClass2 = snirf2.stim(1,2);

allS = [oldStimClass2.data];
% total number of stimuli
oldStimClass2Length = size(oldStimClass2.data,1);

cueOnsetIndex = 1:4:360;
cueOnsetIndex2 = 1:oldStimClass2Length;

% total number of trials
trialNum2 = sum(ismember(cueOnsetIndex2,cueOnsetIndex));

% snirf 2
leftOnsetM = StimClass('leftMultiple');
rightOnsetM = StimClass('rightMultiple');
centerOnsetM = StimClass('centerMultiple');

leftMultiIndex = indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1;
rightMultiIndex = indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1;
centerMultiIndex = indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==1;

% leftSingleIndex = leftSingleIndex(leftSingleIndex<= lastTrial2);
% rightSingleIndex = rightSingleIndex(rightSingleIndex<= lastTrial2);
% centerSingleIndex = centerSingleIndex(centerSingleIndex<= lastTrial2);
% leftMultiIndex = leftMultiIndex(leftMultiIndex<= lastTrial2);
% rightMultiIndex = rightMultiIndex(rightMultiIndex<= lastTrial2);
% centerMultiIndex = centerMultiIndex(centerMultiIndex<= lastTrial2);

AddStims(leftOnsetM, allS(cueOnsetIndex(leftMultiIndex)));
AddStims(rightOnsetM, allS(cueOnsetIndex(rightMultiIndex)));
AddStims(centerOnsetM, allS(cueOnsetIndex(centerMultiIndex)));

updateStates(leftOnsetM);
updateStates(rightOnsetM);
updateStates(centerOnsetM);

snirf2.stim(1,1) = leftOnsetM;
snirf2.stim(1,2) = rightOnsetM;
snirf2.stim(1,3) = centerOnsetM;

%snirf2.data.dataTimeSeries = [snirf2.data.dataTimeSeries;snirf3.data.dataTimeSeries;snirf4.data.dataTimeSeries];
%snirf2.data.time = [snirf2.data.time;snirf3.data.time;snirf4.data.time];
%snirf2.data.time = 0:1/50:(size(snirf2.data.dataTimeSeries,1)-1)/50;

snirf2.Save([rawDataFN{1} '.snirf']);

end