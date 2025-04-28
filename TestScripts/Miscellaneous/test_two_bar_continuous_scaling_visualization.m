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
nDivisionsTotal = 1:7;

boxMeasureFinal = nan(1,length(nDivisionsTotal));
cornerBoxRemovalMeasureFinal = nan(1,length(nDivisionsTotal));
planarTrimmingMeasureFinal = nan(1,length(nDivisionsTotal));

for i = 1:length(nDivisionsTotal)
    nDivisions = nDivisionsTotal(i);

    % get component data
    file = sprintf('ContinuousTruss-%d-Division.mat',nDivisions);
    load(file,'algoData');

    cornerBoxRemovalIndex = 1:10;
    planarTrimmingIndex = 11:20;

    cornerBoxRemovalMeasure = nan(1,length(cornerBoxRemovalIndex));
    planarTrimmingMeasure = nan(1,length(planarTrimmingIndex));

    iCornerBoxRemoval = 1;
    iPlanarTrimming = 1;
    for j = 1:length(algoData)
        if(ismember(j,cornerBoxRemovalIndex))
            cornerBoxRemovalMeasure(iCornerBoxRemoval) = algoData(j).TotalMeasureAfterTrimNormalized(end);
            iCornerBoxRemoval = iCornerBoxRemoval + 1;
        else
            planarTrimmingMeasure(iPlanarTrimming) = algoData(j).TotalMeasureAfterTrimNormalized(end);
            iPlanarTrimming = iPlanarTrimming + 1;
        end
    end

    cornerBoxRemovalMeasureFinal(i) = max(cornerBoxRemovalMeasure);
    planarTrimmingMeasureFinal(i) = max(planarTrimmingMeasure);


    % get box data
    file = sprintf('ContinuousTruss-%d-Division-Box.mat',nDivisions);
    load(file,'algoDataBox');
    boxMeasureFinal(i) = algoDataBox.MeasureAfterTrimNormalized(end);
end

figure;
plot(nDivisionsTotal,cornerBoxRemovalMeasureFinal./boxMeasureFinal,'-o','color',color_palette_tol('cyan'),'LineWidth',2);
hold all;
plot(nDivisionsTotal,planarTrimmingMeasureFinal./boxMeasureFinal,'-o','color',color_palette_tol('purple'),'LineWidth',2);
grid minor;
xlabel('Number of Divisions','Interpreter','latex');
ylabel('Normalized Measure $$V/V_{box}$$','Interpreter','latex');
legend({'Corner Box Removal','Planar Trimming'},'Interpreter','latex','location','east');
title('Two Bar Continuous Truss - Tendency to Infinity','FontSize',14,'Interpreter','latex');
xlim([min(nDivisionsTotal),max(nDivisionsTotal)]);
set(gca, 'YScale', 'log');
save_print_figure(gcf,[saveFolder,'TwoBarContinuousResultsPlot-Tendency'],'Size',figureSize,'PrintFormat',{'png','pdf'});


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;
