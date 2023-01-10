% redefine stimulus classes into 6 diff. classes
% to be runned after generating stimuli marks from aux data in Homer3 GUI
% strictly for ExperimentXXGUI folder
% Old, no longer needed

function redefineStimulusClasses(sbjNum,rawDataFN,respData)

dir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum) 'GUI'];
load([dir filesep respData '.mat'],'fixedMaskerList','indexMoviesTest','maskerMovies','numTrials','uniqueMovies');

snirf2 = SnirfLoad([dir filesep rawDataFN{1} '.snirf']);
snirf3 = SnirfLoad([dir filesep rawDataFN{2} '.snirf']);
snirf4 = SnirfLoad([dir filesep rawDataFN{3} '.snirf']);

oldStimClass2 = snirf2.stim(1,2);
oldStimClass3 = snirf3.stim(1,2);
oldStimClass4 = snirf4.stim(1,2);

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

leftSingleIndex = find(indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==0);
rightSingleIndex = find(indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==0);
centerSingleIndex = find(indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==0);
leftMultiIndex = find(indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1);
rightMultiIndex = find(indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1);
centerMultiIndex = find(indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==1);

leftSingleIndex = leftSingleIndex(leftSingleIndex<= lastTrial2);
rightSingleIndex = rightSingleIndex(rightSingleIndex<= lastTrial2);
centerSingleIndex = centerSingleIndex(centerSingleIndex<= lastTrial2);
leftMultiIndex = leftMultiIndex(leftMultiIndex<= lastTrial2);
rightMultiIndex = rightMultiIndex(rightMultiIndex<= lastTrial2);
centerMultiIndex = centerMultiIndex(centerMultiIndex<= lastTrial2);

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

snirf2.Save([rawDataFN{1} '.snirf']);

% snirf 3
leftOnsetS = StimClass('leftSingle');
rightOnsetS = StimClass('rightSingle');
centerOnsetS = StimClass('centerSingle');
leftOnsetM = StimClass('leftMulti');
rightOnsetM = StimClass('rightMulti');
centerOnsetM = StimClass('centerMulti');

leftSingleIndex = find(indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==0);
rightSingleIndex = find(indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==0);
centerSingleIndex = find(indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==0);
leftMultiIndex = find(indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1);
rightMultiIndex = find(indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1);
centerMultiIndex = find(indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==1);

leftSingleIndex = leftSingleIndex((lastTrial2<leftSingleIndex) & (leftSingleIndex <= lastTrial3));
rightSingleIndex = rightSingleIndex((lastTrial2<rightSingleIndex) & (rightSingleIndex <= lastTrial3));
centerSingleIndex = centerSingleIndex((lastTrial2<centerSingleIndex) & (centerSingleIndex <= lastTrial3));
leftMultiIndex = leftMultiIndex((lastTrial2<leftMultiIndex) & (leftMultiIndex <= lastTrial3));
rightMultiIndex = rightMultiIndex((lastTrial2<rightMultiIndex) & (rightMultiIndex <= lastTrial3));
centerMultiIndex = centerMultiIndex((lastTrial2<centerMultiIndex) & (centerMultiIndex <= lastTrial3));

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

snirf3.stim(1,1) = leftOnsetS;
snirf3.stim(1,2) = rightOnsetS;
snirf3.stim(1,3) = centerOnsetS;
snirf3.stim(1,4) = leftOnsetM;
snirf3.stim(1,5) = rightOnsetM;
snirf3.stim(1,6) = centerOnsetM;

snirf3.Save([rawDataFN{2} '.snirf']);

% snirf 4
leftOnsetS = StimClass('leftSingle');
rightOnsetS = StimClass('rightSingle');
centerOnsetS = StimClass('centerSingle');
leftOnsetM = StimClass('leftMulti');
rightOnsetM = StimClass('rightMulti');
centerOnsetM = StimClass('centerMulti');

leftSingleIndex = find(indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==0);
rightSingleIndex = find(indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==0);
centerSingleIndex = find(indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==0);
leftMultiIndex = find(indexMoviesTest(:,2)==2 & indexMoviesTest(:,5)==1);
rightMultiIndex = find(indexMoviesTest(:,2)==1 & indexMoviesTest(:,5)==1);
centerMultiIndex = find(indexMoviesTest(:,2)==3 & indexMoviesTest(:,5)==1);

leftSingleIndex = leftSingleIndex((lastTrial3 < leftSingleIndex) & (leftSingleIndex <= lastTrial4));
rightSingleIndex = rightSingleIndex((lastTrial3 < rightSingleIndex) & (rightSingleIndex <= lastTrial4));
centerSingleIndex = centerSingleIndex((lastTrial3 < centerSingleIndex) & (centerSingleIndex <= lastTrial4));
leftMultiIndex = leftMultiIndex((lastTrial3 < leftMultiIndex) & (leftMultiIndex <= lastTrial4));
rightMultiIndex = rightMultiIndex((lastTrial3 < rightMultiIndex) & (rightMultiIndex <= lastTrial4));
centerMultiIndex = centerMultiIndex((lastTrial3 < centerMultiIndex) & (centerMultiIndex <= lastTrial4));

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

snirf4.stim(1,1) = leftOnsetS;
snirf4.stim(1,2) = rightOnsetS;
snirf4.stim(1,3) = centerOnsetS;
snirf4.stim(1,4) = leftOnsetM;
snirf4.stim(1,5) = rightOnsetM;
snirf4.stim(1,6) = centerOnsetM;

snirf4.Save([rawDataFN{3} '.snirf']);
