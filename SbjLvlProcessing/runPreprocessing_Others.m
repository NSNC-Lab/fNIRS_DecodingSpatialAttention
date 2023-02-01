% script to run all. Can take up to 2-4 hours.
% Preprocess nirs data for Figure 3, 4, 7 and 8.

prompt = 'Warning. It'' take 1-3 days to run the whole thing. Are you sure? (Y/N)';
x = input(prompt);

if strcmp(x,'Y')

    sbjList = {'12','13','14','15','16','19','21','22','23','24','25'};
    
    rawDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS';
    
    for i = 1:length(sbjList)
        load([rawDir filesep 'Experiment' sbjList{i} filesep sbjList{i} '.mat'],'s');
        preprocessFNIRS06_CV_GLM_ssBeta_DiffTLen_MultiOnly(s,2,1,1,1,0);
        preprocessFNIRS06_CV_GLM_ssBeta_SingleChn_MultiOnly(s,2,1,1);
        preprocessFNIRS06_MultiOnly(s.name,s.fName,s.movieList,0);
        preprocessFNIRS06_LOFO_MultiOnly(s,2,1,1,1,1,0);
    end
    
    load([rawDir filesep 'Experiment' sbjList{i} filesep '08.mat'],'s');
    preprocessFNIRS06_CV_GLM_ssBeta_DiffTLen(s,2,1,1,1,0);
    preprocessFNIRS06_CV_GLM_ssBeta_SingleChn(s,2,1,1);
    preprocessFNIRS06(s.name,s.fName,s.movieList,0);
    preprocessFNIRS06_LOFO(s,2,1,1,1,1,0);
    
end