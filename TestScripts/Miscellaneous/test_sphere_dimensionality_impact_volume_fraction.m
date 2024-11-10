%% Cleanup
close all;
fclose all;
clear all;
clc;
more off;
diary off;

%% debugging
rng(7);


%% Documentation / Archive
rngState = rng;
saveFolder = save_diary_files(mfilename);
goldenRatio = (1+sqrt(5))/2;
figureSize = [goldenRatio 1]*8.5;


%% 
R = 1;
nDimensionMaxPlot = 20;

nDimensionMaxTest = 2*nDimensionMaxPlot;
allDimension = 1:nDimensionMaxTest;

volumeCoefficientComplete = pi.^(allDimension'./2)./gamma(allDimension'./2+1);
volumeComplete = volumeCoefficientComplete*(R.^i);

radiusComponent = nan(nDimensionMaxTest,nDimensionMaxTest);
volumeFraction = nan(nDimensionMaxTest,nDimensionMaxTest);
for i=1:nDimensionMaxTest
    nComponentPossible = divisors(i)';
    radiusComponent(nComponentPossible,i) = R./sqrt(nComponentPossible);
    nDimensionComponent = i./nComponentPossible;

    volumeCoefficientComponent = volumeCoefficientComplete(nDimensionComponent);
    volumeComponent = volumeCoefficientComponent.*(radiusComponent(nComponentPossible,i).^nDimensionComponent);
    volumeFraction(nComponentPossible,i) = (volumeComponent.^nComponentPossible)./volumeComplete(i);
end

figure;
hold all;
for i=1:nDimensionMaxPlot
    validEntries = ~isnan(volumeFraction(i,:));

    currentFraction = volumeFraction(i,validEntries);
    currentValidDimension = allDimension(validEntries);

    if(i==1)
        legendDescription = 'Complete Solution Space';
    else
        legendDescription = sprintf('%d Components',i);
    end

    plot(currentValidDimension,currentFraction,'-*','DisplayName',legendDescription);
end
xlim([1 nDimensionMaxPlot]);
set(gca, 'YScale', 'log');
legend('location','southwest');
grid minor;
save_print_figure(gcf,[saveFolder,'VolumeLostToComplete'],'Size',figureSize);

volumeIncrease = volumeFraction./min(volumeFraction,[],1);
figure;
hold all;
for i=1:nDimensionMaxPlot
    validEntries = ~isnan(volumeIncrease(i,:));

    currentIncrease = volumeIncrease(i,validEntries);
    currentValidDimension = allDimension(validEntries);

    if(i==1)
        legendDescription = 'Complete Solution Space';
    else
        legendDescription = sprintf('%d Components',i);
    end

    plot(currentValidDimension,currentIncrease,'-*','DisplayName',legendDescription);
end
xlim([1 nDimensionMaxPlot]);
xlabel('Problem Dimension','Interpreter','latex');
ylabel('Volume Increase ($$V/V_{box}$$)','Interpreter','latex');
set(gca, 'YScale', 'log', 'FontSize', 14);
% legend('location','northwest','FontSize',8);
grid minor;
save_print_figure(gcf,[saveFolder,'VolumeGainedToBox'],'Size',figureSize,'PrintFormat',{'png','pdf'});


%% reference for two dimensions
% convergence2Component = 1:300;
% figure;
% plot(gamma(convergence2Component/2+1)./((gamma(convergence2Component/4+1)).^2.*2.^(convergence2Component./2)));
% grid minor;

