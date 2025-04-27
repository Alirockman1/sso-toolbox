%% Cleanup
fclose all;
close all;
clear all;
clc;
more off;
diary off;


%% debugging
rng(4);


%% Documentation / Archive
rngState = rng;
saveFolder = save_diary_files(mfilename);
goldenRatio = (1+sqrt(5))/2;
figureSize = [goldenRatio 1]*8.5;


%% plot results
nDivisions = 7;
file = sprintf('ContinuousTruss-%d-Division.mat',nDivisions);
load(file,'algoData');

cornerBoxRemovalIndex = 1:10;
planarTrimmingIndex = 11:20;

cornerBoxRemovalMeasure = zeros(1,length(cornerBoxRemovalIndex));
cornerBoxRemovalSamples = zeros(1,length(cornerBoxRemovalIndex));
planarTrimmingMeasure = zeros(1,length(planarTrimmingIndex));
planarTrimmingSamples = zeros(1,length(planarTrimmingIndex));

iCornerBoxRemoval = 1;
iPlanarTrimming = 1;
for i = 1:length(algoData)
    if(ismember(i,cornerBoxRemovalIndex))
        cornerBoxRemovalMeasure(iCornerBoxRemoval) = algoData(i).TotalMeasureAfterTrimNormalized(end);
        cornerBoxRemovalSamples(iCornerBoxRemoval) = algoData(i).NumberEvaluatedSamples(end);
        iCornerBoxRemoval = iCornerBoxRemoval + 1;
    else
        planarTrimmingMeasure(iPlanarTrimming) = algoData(i).TotalMeasureAfterTrimNormalized(end);
        planarTrimmingSamples(iPlanarTrimming) = algoData(i).NumberEvaluatedSamples(end);
        iPlanarTrimming = iPlanarTrimming + 1;
    end
end

figure;
plot(cornerBoxRemovalSamples,cornerBoxRemovalMeasure,'-o','color',color_palette_tol('cyan'),'LineWidth',2);
hold all;
plot(planarTrimmingSamples,planarTrimmingMeasure,'-o','color',color_palette_tol('purple'),'LineWidth',2);
grid minor;
xlabel('Number of Evaluated Samples Per Iteration','Interpreter','latex');
ylabel('Normalized Measure $$V/V_{ds}$$','Interpreter','latex');
legend({'Corner Box Removal','Planar Trimming'},'Interpreter','latex','location','east');
title(sprintf('Two Bar Continuous Truss - %d Divisions',nDivisions),'FontSize',14,'Interpreter','latex');
xlim([min(cornerBoxRemovalSamples),max(cornerBoxRemovalSamples)]);
save_print_figure(gcf,[saveFolder,sprintf('TwoBarContinuousResultsPlot-%d-Divisions',nDivisions)],'Size',figureSize,'PrintFormat',{'png','pdf'});


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;
