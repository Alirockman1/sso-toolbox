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
nBar = 6;
nDimension = 3;
file = sprintf('%02d-Bar-%dD-Batch.mat',nBar,nDimension);
load(file,'algoDataBox','algoDataComponent');

% load box shaped results  
boxShapedIndex = 1:10;
boxShapedMeasure = zeros(1,length(boxShapedIndex));
boxShapedSamples = zeros(1,length(boxShapedIndex));

for i = 1:length(algoDataBox)
    boxShapedMeasure(i) = algoDataBox(i).MeasureAfterTrimNormalized(end);
    boxShapedSamples(i) = algoDataBox(i).NumberEvaluatedSamples(end);
end

% load component shaped results
cornerBoxRemovalIndex = 1:10;
planarTrimmingIndex = 11:20;
cornerBoxRemovalMeasure = zeros(1,length(cornerBoxRemovalIndex));
cornerBoxRemovalSamples = zeros(1,length(cornerBoxRemovalIndex));
planarTrimmingMeasure = zeros(1,length(planarTrimmingIndex));
planarTrimmingSamples = zeros(1,length(planarTrimmingIndex));

iCornerBoxRemoval = 1;
iPlanarTrimming = 1;
for i = 1:length(algoDataComponent)
    if(ismember(i,cornerBoxRemovalIndex))
        cornerBoxRemovalMeasure(iCornerBoxRemoval) = algoDataComponent(i).TotalMeasureAfterTrimNormalized(end);
        cornerBoxRemovalSamples(iCornerBoxRemoval) = algoDataComponent(i).NumberEvaluatedSamples(end);
        iCornerBoxRemoval = iCornerBoxRemoval + 1;
    else
        planarTrimmingMeasure(iPlanarTrimming) = algoDataComponent(i).TotalMeasureAfterTrimNormalized(end);
        planarTrimmingSamples(iPlanarTrimming) = algoDataComponent(i).NumberEvaluatedSamples(end);
        iPlanarTrimming = iPlanarTrimming + 1;
    end
end

% plot comparison
figure;
plot(cornerBoxRemovalSamples,cornerBoxRemovalMeasure,'-o','color',color_palette_tol('cyan'),'LineWidth',2);
hold all;
plot(planarTrimmingSamples,planarTrimmingMeasure,'-o','color',color_palette_tol('purple'),'LineWidth',2);
plot(boxShapedSamples,boxShapedMeasure,'-o','color','k','LineWidth',2);
grid minor;
xlabel('Number of Evaluated Samples Per Iteration','Interpreter','latex');
ylabel('Normalized Measure $$V/V_{ds}$$','Interpreter','latex');
legend({'Corner Box Removal','Planar Trimming','Box-shaped'},'Interpreter','latex','location','east');
title(sprintf('%d Bar Truss - %d Design Variables',nBar,nDimension),'FontSize',14,'Interpreter','latex');
xlim([min(cornerBoxRemovalSamples),max(cornerBoxRemovalSamples)]);
set(gca, 'YScale', 'log');
save_print_figure(gcf,[saveFolder,sprintf('%02d-Bar-%dD-ResultsPlot',nBar,nDimension)],'Size',figureSize,'PrintFormat',{'png','pdf'});


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;
