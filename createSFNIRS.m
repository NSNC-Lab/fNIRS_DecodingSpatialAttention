function createSFNIRS(sbjNum,fName,restFN,respFN,movieList,startT,endT)

rawDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];

s.name = sbjNum;
s.fName = fName; % rawDataFN
s.movieList = movieList;
s.restFN = restFN;
s.resp = respFN;
s.startT = startT;
% if endT == -1, end
s.endT = endT;

save([rawDataDir filesep sbjNum '.mat'],'s');

end

