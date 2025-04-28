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


%% Considered material properties
% Material properties arrays
youngsModulusMaterial = [200e9, 70e9, 116e9, 45e9]; % [Pa] Steel, Aluminum, Titanium, Magnesium


%% system definition
systemFunction = @truss_two_bar_continuous_area_only;


%% design space definition
systemParameter = [1,1,45e9,45e9,1000];
nDivisions = 7;

minimumArea = 0.001;
maximumArea = 1;

nIter = 150;
nSample = 200;
trimmingPasses = 'single';
trimmingOrder = 'lth';

performanceLowerLimitFactor = 0.9;
performanceUpperLimitFactor = 1.1;

evolutionIterationIndices = [];


%%
designSpaceLowerBound = [repmat(minimumArea,1,nDivisions),repmat(minimumArea,1,nDivisions)];
designSpaceUpperBound = [repmat(maximumArea,1,nDivisions),repmat(maximumArea,1,nDivisions)];
componentIndex = {[1:nDivisions],[nDivisions+1:2*nDivisions]};
initialDesign = mean([designSpaceLowerBound;designSpaceUpperBound],1);


%% compute initial design
bottomUpMapping = BottomUpMappingFunction(systemFunction,systemParameter);
performanceInitial = bottomUpMapping.response(initialDesign);

performanceLowerLimit = performanceLowerLimitFactor*performanceInitial;
performanceUpperLimit = performanceUpperLimitFactor*performanceInitial;


%% solution space computation
designEvaluator = DesignEvaluatorBottomUpMapping(...
    bottomUpMapping,...
    performanceLowerLimit,...
    performanceUpperLimit);

timeElapsedBox = tic;
optionsBox = sso_stochastic_options('box',...
    'NumberSamplesPerIterationExploration',nSample,...
    'NumberSamplesPerIterationConsolidation',nSample,...
    'FixIterNumberExploration',true,...
    'FixIterNumberConsolidation',true,...
    'MaxIterExploration',nIter,...
    'MaxIterConsolidation',nIter,...
    'GrowthRate',0.1^nDivisions,...
    'UseAdaptiveGrowthRate',true,...
    'ApplyLeanness','never',...
    'TrimmingOperationOptions',{'PassesCriterion',trimmingPasses},...
    'TrimmingOrderOptions',{'OrderPreference',trimmingOrder});
[boxSolutionSpace,optimizationDataBox] = sso_box_stochastic(designEvaluator,...
    initialDesign,designSpaceLowerBound,designSpaceUpperBound,optionsBox);
toc(timeElapsedBox)


%% plot results and metrics
% Individual results - Planar Trimming
resultsFolder = [saveFolder,'ResultsBox/'];
mkdir(resultsFolder);
algoDataBox = postprocess_sso_box_stochastic(optimizationDataBox);
plot_sso_box_stochastic_metrics(algoDataBox,...
    'SaveFolder',resultsFolder,...
    'CloseFigureAfterSaving',true,...
    'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

