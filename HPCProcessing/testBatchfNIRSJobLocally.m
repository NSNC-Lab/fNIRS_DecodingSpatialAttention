opIsLocal = 1;

sbjList = {'12','13','14','15','16','19','21','22','23','24'};
for i = 1:length(sbjList)
    sbjNum = sbjList{i};
    %rawDataDir = ['/projectnb2/binaural/mhn/RawDatafNIRS/Experiment' num2str(sbjNum)];
    rawDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\RawDatafNIRS\Experiment' num2str(sbjNum)];
    sLoad = load([rawDataDir filesep num2str(sbjNum)]);
    s = sLoad.s;
    preprocessFNIRS06_CV_GLM_ssBeta_MultiOnly_Parallel(s,2,1,1,0.5,opIsLocal);
    %preprocessFNIRS06_CV_GLM_ssBeta_MultiOnly_Parallel(s,3,1,1,0.5,opIsLocal);
    %preprocessFNIRS06_CV_GLM_ssBeta_DiffTLen_MultiOnly(s,2,1,1,1,0,opIsLocal);
    %preprocessFNIRS06_CV_GLM_ssBeta_DiffTLen_MultiOnly(s,3,1,1,1,0,opIsLocal);
    %createSingleTrialHRF_SSBeta(s.name,s.movieList,2);
    %calcCumSum_DiffStartT_AllChn_CompClassifiers_Final(s.name,2,1);
    close all;
end