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
nDimensionMax = 50;
allDimension = 1:nDimensionMax;

volumeCoefficientComplete = pi.^(allDimension'./2)./gamma(allDimension'./2+1);
volumeComplete = volumeCoefficientComplete*(R.^i);

radiusComponent = nan(nDimensionMax,nDimensionMax);
volumeFraction = nan(nDimensionMax,nDimensionMax);
for i=1:nDimensionMax
    nComponentPossible = divisors(i)';
    radiusComponent(nComponentPossible,i) = R./sqrt(nComponentPossible);
    nDimensionComponent = i./nComponentPossible;

    volumeCoefficientComponent = volumeCoefficientComplete(nDimensionComponent);
    volumeComponent = volumeCoefficientComponent.*(radiusComponent(nComponentPossible,i).^nDimensionComponent);
    volumeFraction(nComponentPossible,i) = (volumeComponent.^nComponentPossible)./volumeComplete(i);
end

figure;
hold all;
for i=1:nDimensionMax
    validEntries = ~isnan(volumeFraction(i,:));

    currentFraction = volumeFraction(i,validEntries);
    currentValidDimension = allDimension(validEntries);

    plot(currentValidDimension,currentFraction,'-*','DisplayName',sprintf('%d Components',i));
end
set(gca, 'YScale', 'log');
legend;
grid minor;

volumeIncrease = volumeFraction./min(volumeFraction,[],1);
figure;
hold all;
for i=1:nDimensionMax
    validEntries = ~isnan(volumeIncrease(i,:));

    currentIncrease = volumeIncrease(i,validEntries);
    currentValidDimension = allDimension(validEntries);

    plot(currentValidDimension,currentIncrease,'-*','DisplayName',sprintf('%d Components',i));
end
set(gca, 'YScale', 'log');
legend;
grid minor;

%% reference for two dimensions
convergence2Component = 1:300;
figure;
plot(gamma(convergence2Component/2+1)./((gamma(convergence2Component/4+1)).^2.*2.^(convergence2Component./2)));
grid minor;
