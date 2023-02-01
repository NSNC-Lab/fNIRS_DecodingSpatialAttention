% script to run all. Can take up to 1 full day. Parallel computing version is
% available at ...
% Preprocess nirs data for 2-class results for Figure 5 and 6.

prompt = 'Warning. It'' take up to 1 day to run the whole thing. Are you sure? (Y/N)';
x = input(prompt);

if strcmp(x,'Y')

    sbjList = {'12','13','14','15','16','19','21','22','23','24','25'};
    
    rawDir = 'C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS';
    
    for i = 1:length(sbjList)
        load([rawDir filesep 'Experiment' sbjList{i} filesep sbjList{i} '.mat'],'s');
        preprocessFNIRS06_CV_GLM_ssBeta_MultiOnly_Version02(s,2,1,1,0.5);
    end
    
    load([rawDir filesep 'Experiment' sbjList{i} filesep '08.mat'],'s');
    preprocessFNIRS06_CV_GLM_ssBeta_Version02(s,2,1,1,0.5);
    
end