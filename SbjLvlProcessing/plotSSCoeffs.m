% plot average of ss coefficients across all subjects for each channel


sbjList = {'12','13','14','15','16','17'};

figSaveDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\GroupResultAll\Figures'];

% hardcoded
ssBetaTotal = zeros(49,36,2);
ssBetaList = zeros(length(sbjList),36,2);


for i = 1:length(sbjList)
    sbjNum = sbjList{i};
    
    processedDataDir = ['C:\Users\mn0mn\Documents\ResearchProjects\spatailAttentionProject\ProcessedDatafNIRS\Experiment' num2str(sbjNum)];

    fileName = [processedDataDir filesep 'ssBeta.mat'];
    load(fileName,'ssBeta');
    
    ssBetaTotal = ssBetaTotal + ssBeta{1};
    ssBetaList(i,:,:) = ssBeta{1}(end,:,:);
    
end

ssBetaTotal = ssBetaTotal./length(sbjList);
ssBetaSTD = std(ssBetaList,[],1);

chnXAxis = 1:size(ssBeta{1},2);
subplot(1,2,1);hold on;
%scatter(chnXAxis,ssBeta{1}(end,:,1));
errorbar(chnXAxis,ssBetaTotal(end,:,1),squeeze(ssBetaSTD(1,:,1)),'o');
title('SS Coefficients for HbO');

subplot(1,2,2);hold on;
%scatter(chnXAxis,ssBeta{1}(end,:,2));
errorbar(chnXAxis,ssBetaTotal(end,:,2),squeeze(ssBetaSTD(1,:,2)),'o');
title('SS Coefficients for HbR');

hold off;
fn = sprintf('ss_Coefficients');
print(gcf,[figSaveDir filesep fn],'-dpng','-r250');