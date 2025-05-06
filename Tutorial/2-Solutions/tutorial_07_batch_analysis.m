%TUTORIAL_07_BATCH_ANALYSIS Solving a problem with multiple sets of options
%   TUTORIAL_07_BATCH_ANALYSIS solves the same SSO problem using different
%   sets of options; this can help analyze the performance of the algorithm,
%   or just aid in finding the best options for a given problem.

%% cleanup
close all;
fclose all;
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
defaultFigureSize = [goldenRatio 1]*8.5;


%% base 
systemFunction = @tutorial_01_euclidean_distance_3d;
systemParameter = [0,0,0];
bottomUpMapping = BottomUpMappingFunction(systemFunction,systemParameter);

% performance requirements
performanceLowerLimit = -inf;
performanceUpperLimit = 5;
designEvaluator = DesignEvaluatorBottomUpMapping(bottomUpMapping,performanceLowerLimit,performanceUpperLimit);

% design space 
designSpaceLowerBound = [-6 -6 -6];
designSpaceUpperBound = [ 6  6  6];
initialPoint = [0,0,0];


%% run batch analysis
[solutionSpace,problemData,iterData,algoData,batchOptions] = batch_sso_stochastic_analysis(...
    'BatchTestHollowSphere.xlsx',...
    designEvaluator,...
    initialPoint,...
    designSpaceLowerBound,...
    designSpaceUpperBound,...
    'FixedOptions',{'UseAdaptiveGrowthRate',false});


%% plot data
plot_batch_sso_stochastic_analysis_box_metrics(batchOptions,algoData,'SaveFolder',saveFolder);


%% stop transcript 
save([saveFolder,'Data.mat']);
diary off;

