%TEST_BATCH_COMPONENT_HOLLOW_SPHERE Batch performance test for hollow sphere
%   TEST_BATCH_COMPONENT_HOLLOW_SPHERE performs component solution space 
%   optimizations in batch, saving the performance data. Factors can be
%	attributed according to the corresponding XLSX file.
%
%   Copyright 2025 Eduardo Rodrigues Della Noce
%   SPDX-License-Identifier: Apache-2.0

%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
%   
%       http://www.apache.org/licenses/LICENSE-2.0
%   
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

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
figureSize = [goldenRatio 1]*8.5;


%% setup problem parameters
nDivision = 15; % number of design variables (discritization)
componentIndex = {1:3,4:6,7:9,10:12,13:15};

initialDistanceWidth = 4; % distance in longitudinal (x) direction
initialDistanceHeight = 2; % distance in vertical (y) direction
nInterpolationPoint = 5000; % number of interpolation points


%% Problem Setup
systemFunction = @bead_slide_time;
systemParameter(1) = initialDistanceWidth; % distance in longitudinal (x) direction
systemParameter(2) = initialDistanceHeight; % distance in vertical (y) direction
systemParameter(3) = nInterpolationPoint; % number of interpolation points

designSpaceLowerBound = -initialDistanceHeight*ones(1,nDivision);
designSpaceUpperBound = initialDistanceHeight*ones(1,nDivision);


%% Run Initial Optimization (Find Optimal Time)
% initial guess - linear interpolation
initialDesign = linspace(initialDistanceHeight,0,nDivision+2); % in initial interpolation consider end-points... 
initialDesign = initialDesign(2:end-1); % ... which are then removed
pathWidthDivision = linspace(0,initialDistanceWidth,nDivision+2);

% run optimization
bottomUpMapping = BottomUpMappingFunction(systemFunction,systemParameter);
[designOptimal,timeOptimal,optimizationOutput,systemResponseData] = design_optimize_quantities_of_interest(...
    bottomUpMapping,...
    initialDesign,...
    designSpaceLowerBound,...
    designSpaceUpperBound,...
    @(performanceMeasure)performanceMeasure,...
    'OptimizationMethodOptions',{'Display','iter-detailed'});
referenceTime = bottomUpMapping.response(initialDesign);


%% Solve  Problem Analytically (Reference)
radiusAngleOptimal = fsolve(@(rt)brachistochrone_solve_analytical(rt,systemParameter),[1 1]);
anglePathAnalytical = linspace(0,radiusAngleOptimal(2),nInterpolationPoint);
pathWidthAnalytical = radiusAngleOptimal(1)*(anglePathAnalytical-sin(anglePathAnalytical));
pathHeightAnalytical = radiusAngleOptimal(1)*(-1+cos(anglePathAnalytical)) + initialDistanceHeight;


%% Perform Solution Space Computation
performanceLowerLimit = -inf;
performanceUpperLimit = referenceTime;
designEvaluator = DesignEvaluatorBottomUpMapping(bottomUpMapping,performanceLowerLimit,performanceUpperLimit);


%% run batch analysis
[solutionSpace,optimizationData,algoData,batchOptions] = batch_sso_stochastic_analysis(...
    'BatchTestBrachistochrone.xlsx',...
    designEvaluator,...
    designOptimal,...
    designSpaceLowerBound,...
    designSpaceUpperBound,...
    componentIndex,...
    'FixedOptions',{...
        'UseAdaptiveGrowthRate',true,...
        'FixIterNumberExploration',true,...
        'FixIterNumberConsolidation',true,...
        'TrimmingOperationOptions',{'PassesCriterion','reduced'},...
        'TrimmingOrderOptions',{'OrderPreference','score'},...
        'LoggingLevel','all',...
        'CandidateSpaceConstructor',@CandidateSpaceConvexHull,...
        'TrimmingMethodFunction',@component_trimming_method_planar_trimming});


%% plot data
plot_batch_sso_stochastic_analysis_component_metrics(batchOptions,algoData,'SaveFolder',saveFolder);


%% stop transcript 
save([saveFolder,'Data.mat']);
diary off;

