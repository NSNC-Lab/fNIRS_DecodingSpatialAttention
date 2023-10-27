% sudan run 12-25 one shot
clear all
clc
sbjList = {'12','13','14','15','16','19','21','22','23','24','25'};
for i = 1:length(sbjList)
    sbj = sbjList{i};
    preprocessFNIRS06_CV_GLM_ssBeta_DiffTLen_MultiOnly_Sudan(sbj, 2, 1, 1, 1, 1, 1);
end

