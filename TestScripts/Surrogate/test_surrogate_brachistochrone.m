%test_surrogate_brachistochrone Surrogate model of the brachistochrone problem
%   test_surrogate_brachistochrone creates a surrogate model for the 
%   brachistochrone problem. This is done using an active set strategy.
%
%   Copyright 2024 Eduardo Rodrigues Della Noce
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


%% setup problem parameters
systemFunction = @bead_slide_time;
nDivision = 15; % number of design variables (discritization)
initialDistanceWidth = 4; % distance in longitudinal (x) direction
initialDistanceHeight = 2; % distance in vertical (y) direction
nInterpolationPoint = 5000; % number of interpolation points
rangeAcceptable = 0.3; % accept responses that are (input*100)% worse than optimal solution


%% Problem Setup
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


%% Solve  Problem Analytically (Reference)
radiusAngleOptimal = fsolve(@(rt)brachistochrone_solve_analytical(rt,systemParameter),[1 1]);
anglePathAnalytical = linspace(0,radiusAngleOptimal(2),nInterpolationPoint);
pathWidthAnalytical = radiusAngleOptimal(1)*(anglePathAnalytical-sin(anglePathAnalytical));
pathHeightAnalytical = radiusAngleOptimal(1)*(-1+cos(anglePathAnalytical)) + initialDistanceHeight;


%% Perform Solution Space Computation
performanceLowerLimit = -inf;
performanceUpperLimit = timeOptimal*(1+rangeAcceptable);


%% Training
options = active_learning_model_training_options(...
    'FixIter',true,...
    'MaxIter',10,...
    'NumberCandidatesInitial',10000,...
    'ExplorationToBoundaryRatio',0.5,...
    'NumberSamplesEvaluation',50,...
    'LoggingLevel','all',...
    'NewModelInitialRegion',designOptimal);

[fastForwardModel,problemData,iterData] = active_learning_model_training(...
    bottomUpMapping,...
    designSpaceLowerBound,...
    designSpaceUpperBound,...
    performanceLowerLimit,...
    performanceUpperLimit,...
    options);

method_description = sprintf('Training: Active Learning %d Samples',...
    size(fastForwardModel.DesignSampleTrain,1));


%% Metrics
algoData = postprocess_active_learning_model_training(problemData,iterData);
plot_active_learning_model_training_metrics(algoData,'SaveFolder',saveFolder);


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

