% redefine stimulus classes into 6 diff. classes
% to be runned after generating stimuli marks from aux data in Homer3 GUI
% Combined multiple snirf files into one.
% strictly for ExperimentXXGUI folder

function redefineStimulusClassesCombined(sbjNum,rawDataFN,respData)

dir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum) 'GUI'];
load([dir filesep respData '.mat'],'fixedMaskerList','indexMoviesTest','maskerMovies','numTrials','uniqueMovies');

snirf2 = SnirfLoad([dir filesep rawDataFN{1} '.snirf']);
snirf3 = SnirfLoad([dir filesep rawDataFN{2} '.snirf']);
snirf4 = SnirfLoad([dir filesep rawDataFN{3} '.snirf']);

oldStimClass2 = snirf2.stim(1,2);
oldStimClass3 = snirf3.stim(1,2);
oldStimClass4 = snirf4.stim(1,2);

snirf2TLen = size(snirf2.data.time,1)/(50);
snirf3TLen = size(snirf3.data.time,1)/(50);
%oldStimClass3.data(:,1) = oldStimClass3.data(:,1)-snirf2TLen;
oldStimClass4.data(:,1) = oldStimClass4.data(:,1)+snirf2TLen+snirf3TLen;

allS = [oldStimClass2.data; oldStimClass3.data; oldStimClass4.data];
% total number of stimuli
oldStimClass2Length = size(oldStimClass2.data,1);
oldStimClass3Length = oldStimClass2Length + size(oldStimClass3.data,1);
oldStimClass4Length = oldStimClass3Length + size(oldStimClass4.data,1);

cueOnsetIndex = 1:4:1080;
cueOnsetIndex2 = 1:oldStimClass2Length;
cueOnsetIndex3 = oldStimClass2Length+1:oldStimClass3Length;
cueOnsetIndex4 = oldStimClass3Length+1:oldStimClass4Length;

% total number of trials
trialNum2 = sum(ismember(cueOnsetIndex2,cueOnsetIndex));
trialNum3 = sum(ismember(cueOnsetIndex3,cueOnsetIndex));
trialNum4 = sum(ismember(cueOnsetIndex4,cueOnsetIndex));

lastTrial2 = trialNum2;
lastTrial3 = lastTrial2 + trialNum3;
lastTrial4 = lastTrial3 + trialNum4;

% snirf 2
leftOnsetS = StimClass('leftSingle');
rightOnsetS = StimClass('rightSingle');
centerOnsetS = StimClass('centerSingle');
leftOnsetM = StimClass('leftMulti');
rightOnsetM = StimClass('rightMulti');
centerOnsetM = StimClass('centerMulti');

leftSingleIndex = indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==0;
rightSingleIndex = indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==0;
centerSingleIndex = indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==0;
leftMultiIndex = indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1;
rightMultiIndex = indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1;
centerMultiIndex = indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==1;

% leftSingleIndex = leftSingleIndex(leftSingleIndex<= lastTrial2);
% rightSingleIndex = rightSingleIndex(rightSingleIndex<= lastTrial2);
% centerSingleIndex = centerSingleIndex(centerSingleIndex<= lastTrial2);
% leftMultiIndex = leftMultiIndex(leftMultiIndex<= lastTrial2);
% rightMultiIndex = rightMultiIndex(rightMultiIndex<= lastTrial2);
% centerMultiIndex = centerMultiIndex(centerMultiIndex<= lastTrial2);

AddStims(leftOnsetS, allS(cueOnsetIndex(leftSingleIndex)));
AddStims(rightOnsetS, allS(cueOnsetIndex(rightSingleIndex)));
AddStims(centerOnsetS, allS(cueOnsetIndex(centerSingleIndex)));
AddStims(leftOnsetM, allS(cueOnsetIndex(leftMultiIndex)));
AddStims(rightOnsetM, allS(cueOnsetIndex(rightMultiIndex)));
AddStims(centerOnsetM, allS(cueOnsetIndex(centerMultiIndex)));

updateStates(leftOnsetS);
updateStates(rightOnsetS);
updateStates(centerOnsetS);
updateStates(leftOnsetM);
updateStates(rightOnsetM);
updateStates(centerOnsetM);

snirf2.stim(1,1) = leftOnsetS;
snirf2.stim(1,2) = rightOnsetS;
snirf2.stim(1,3) = centerOnsetS;
snirf2.stim(1,4) = leftOnsetM;
snirf2.stim(1,5) = rightOnsetM;
snirf2.stim(1,6) = centerOnsetM;

snirf2.data.dataTimeSeries = [snirf2.data.dataTimeSeries;snirf3.data.dataTimeSeries;snirf4.data.dataTimeSeries];
%snirf2.data.time = [snirf2.data.time;snirf3.data.time;snirf4.data.time];
snirf2.data.time = 0:1/50:(size(snirf2.data.dataTimeSeries,1)-1)/50;

snirf2.Save([rawDataFN{1} 'Combined.snirf']);

end