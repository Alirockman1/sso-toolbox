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
nDivisions = 3;

minimumArea = 0.001;
maximumArea = 1;

trimmingPasses = 'single';
trimmingOrder = 'score';

performanceLowerLimitFactor = 0.9;
performanceUpperLimitFactor = 1.1;

evolutionIterationIndices = [];


%%
designSpaceLowerBound = [repmat(minimumArea,1,nDivisions),repmat(minimumArea,1,nDivisions)];
designSpaceUpperBound = [repmat(maximumArea,1,nDivisions),repmat(maximumArea,1,nDivisions)];
componentIndex = {[1:nDivisions]',[nDivisions+1:2*nDivisions]'};
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

fixedOptions = {...
    'MaxIter',300,...
    'FixIterNumberExploration',true,...
    'FixIterNumberConsolidation',true,...
    'CandidateSpaceOptions',{'NormalizeGrowthDirection',true,'CheckRedundantTrimmingGrowth',false,'CheckRedundantTrimmingUpdate',true,'CheckDuplicatePointsGrowth',false,'CheckDuplicatePointsUpdate',true,'MeasureEstimationFactor',2},...
    'TrimmingMethodOptions',{'NormalizeVariables',true},...
    'GrowthRate',0.1^nDivisions,...
    'UseAdaptiveGrowthRate',true,...
    'ApplyLeanness','never',...
    'MaximumNumberPaddingSamples',2000,...
    'MinimumNumberPaddingSamples',1000,...
    'UsePaddingSamplesInTrimming',true,...
    'UsePreviousEvaluatedSamplesConsolidation',false,...
    'ShapeSamplesUsefulExploration',false,...
    'ShapeSamplesUsefulConsolidation',true,...
    'TrimmingOperationOptions',{'PassesCriterion',trimmingPasses},...
    'TrimmingOrderOptions',{'OrderPreference',trimmingOrder},...
    'LoggingLevel','all'...
    };

batchOptions = batch_analysis_read_table('BatchTestTruss.xlsx');


%% run batch analysis
[solutionSpace,optimizationData,algoData,batchOptions] = batch_sso_stochastic_analysis(...
    batchOptions,...
    designEvaluator,...
    initialDesign,...
    designSpaceLowerBound,...
    designSpaceUpperBound,...
    componentIndex,...
    'FixedOptions',fixedOptions);


%% plot data
plot_batch_sso_stochastic_analysis_component_metrics(batchOptions,algoData,'SaveFolder',saveFolder);


%% Save and Stop Transcripting
save([saveFolder,'Data.mat'],'-v7.3');
diary off;

