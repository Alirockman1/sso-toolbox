%TEST_COMPONENT_BRACHISTOCHRONE Brachistochrone / Bead Descent problem 
%   TEST_COMPONENT_BRACHISTOCHRONE uses a discretized version of the bead 
%   descent problem to test the computation of component solution spaces.
%   First, the analytical curve of least time (brachistochrone) is computed, 
%   followed by the box-shaped solution space (with a requirement around the 
%   time it takes for the slide from A to B happen) and the component solution 
%   spaces. The results are compared after all computations are done. 
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

%% cleanup
close all;
fclose all;
clear all;
more off;
diary off;
clc;


%% debugging
rng(4);


%% Documentation / Archive
RNGstate = rng;
saveFolder = save_diary_files(mfilename);
goldenRatio = (1+sqrt(5))/2;
figureSize = [goldenRatio 1]*8.5;


%% setup problem parameters
nDivision = 15; % number of design variables (discritization)
initialDistanceWidth = 4; % distance in longitudinal (x) direction
initialDistanceHeight = 2; % distance in vertical (y) direction
nInterpolationPoint = 5000; % number of interpolation points
rangeAcceptable = 0.3; % accept responses that are (input*100)% worse than optimal solution


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


%% Solve  Problem Analytically (Reference)
radiusAngleOptimal = fsolve(@(rt)brachistochrone_solve_analytical(rt,systemParameter),[1 1]);
anglePathAnalytical = linspace(0,radiusAngleOptimal(2),nInterpolationPoint);
pathWidthAnalytical = radiusAngleOptimal(1)*(anglePathAnalytical-sin(anglePathAnalytical));
pathHeightAnalytical = radiusAngleOptimal(1)*(-1+cos(anglePathAnalytical)) + initialDistanceHeight;


%% Perform Solution Space Computation
performanceLowerLimit = -inf;
performanceUpperLimit = timeOptimal*(1+rangeAcceptable);

optionsBox = sso_stochastic_options('box',...
    'UseAdaptiveGrowthRate',true,...
    'NumberSamplesPerIteration',100,...
    'GrowthRate',0.2,...
    'FixIterNumberExploration',true,...
    'FixIterNumberConsolidation',true,...
    'MaxIterExploration',50,...
    'MaxIterConsolidation',50,...
    'TrimmingOperationOptions',{'PassesCriterion','reduced'},...
    'TrimmingOrderOptions',{'OrderPreference','score'},...
    'LoggingLevel','all');

designEvaluator = DesignEvaluatorBottomUpMapping(bottomUpMapping,performanceLowerLimit,performanceUpperLimit);
[designBoxOptimal,problemDataBox,iterDataBox] = sso_box_stochastic(designEvaluator,...
    designOptimal,designSpaceLowerBound,designSpaceUpperBound,optionsBox);


%% plot
figure; 
hold on;
for i=1:nDivision
    handleSolutionSpacePlot = plot([pathWidthDivision(i+1) pathWidthDivision(i+1)],designBoxOptimal(:,i),'g-','linewidth',12);
end
handleNumericalOptimal = plot(pathWidthDivision,[initialDistanceHeight,designOptimal,0]','bo-','linewidth',1);
plot([0 systemParameter(1)],[0 0],'k--');
handleBrachistochrone = plot(pathWidthAnalytical,pathHeightAnalytical,'r','linewidth',1);
grid minor;
xlim([0 systemParameter(1)])
title(sprintf('Approximate Optimal Time: %gs | Time Requirement: %gs (%g%% Leniency)',...
    timeOptimal,performanceUpperLimit,rangeAcceptable*100));
xlabel('Longitudinal Distance [m]');
ylabel('Vertical Distance [m]');
legend([handleSolutionSpacePlot handleNumericalOptimal handleBrachistochrone],{'Solution Spaces', ...
    'Numerical Optimal Solution','Analytical Optimal Solution'})


%% same problem, component solution spaces
options = sso_stochastic_options('component',...
    'UseAdaptiveGrowthRate',true,...
    'TargetAcceptedRatioExploration',0.9,...
    'GrowthRate',0.2,...
    'FixIterNumberExploration',true,...
    'FixIterNumberConsolidation',true,...
    'MaxIterExploration',50,...
    'MaxIterConsolidation',50,...
    'NumberSamplesPerIterationExploration',300,...
    'NumberSamplesPerIterationConsolidation',300,...
    'CandidateSpaceConstructorExploration',@CandidateSpaceConvexHull,...
    'CandidateSpaceConstructorConsolidation',@CandidateSpaceConvexHull,...
    'TrimmingMethodFunction',@component_trimming_method_planar_trimming,...
    'UsePaddingSamplesInTrimming',true,...
    'LoggingLevel','all',...
    'TrimmingOperationOptions',{'PassesCriterion','reduced'},...
    'TrimmingOrderOptions',{'OrderPreference','score'});

component = {1:3,4:6,7:9,10:12,13:15};
        
[componentSolutionSpace,problemDataComponent,iterDataComponent] = sso_component_stochastic(designEvaluator,...
    designBoxOptimal,designSpaceLowerBound,designSpaceUpperBound,component,options);


%% visualization
for i=1:size(component,2)
    currentBox = designBoxOptimal(:,component{i});
    currentComponentSpace = componentSolutionSpace(i);
    
    figure;
    hold on;
    plot_design_box_3d(gcf,currentBox,'FaceAlpha',0.5,'FaceColor',[0.8 0.8 0.8]);
    componentSolutionSpace(i).plot_candidate_space(gcf,...
        'FaceColor','green','EdgeColor','green','FaceAlpha',0.1);
    grid minor;
    
    %
    muBox(i) = prod(currentBox(2,:)-currentBox(1,:));
    
    %
    muComponent(i) = componentSolutionSpace(i).Measure;
    
    %
    title(sprintf('Box Volume = %gm^3 ; Component Volume = %gm^3 ;',...
            muBox(i),muComponent(i)));
end
fprintf('\nTotal Box Volume: %g\nTotal Component Volume: %g\nTotal Increase: %gx\n',...
    prod(muBox),prod(muComponent),prod(muComponent)/prod(muBox))


%% performance metrics
algoData = postprocess_sso_component_stochastic(problemDataComponent,iterDataComponent);
plot_sso_component_stochastic_metrics(algoData,...
    'SaveFolder',saveFolder,...
    'SaveFigureOptions',{'Size',figureSize});


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

