%TEST_COMPONENT_HOLLOW_SPHERE Component solution spaces for a sphere problem 
%   TEST_COMPONENT_HOLLOW_SPHERE computes a component solution space for a 
%   design problem where the good designs are enclosed in a sphere.
%   The analytical solution for this problem is known and is plotted together
%   with the computed component solution spaces for comparison.
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


%% function call
%
systemFunction = @distance_to_center;
systemParameter = [0,0,0];
%                        x1 x2 x3
designSpaceLowerBound = [-6 -6 -6];
designSpaceUpperBound = [ 6  6  6];
componentIndex = {[1,3],[2]};
%
performanceLowerLimit = -inf;
performanceUpperLimit = 5;
%
initialDesign = [0,0,0];


%% Component Opt - Function
timeElapsedAlgorithm = tic;
options = sso_stochastic_options('component',...
    'NumberSamplesPerIterationExploration',100,...
    'NumberSamplesPerIterationConsolidation',100,...
    'FixIterNumberExploration',true,...
    'FixIterNumberConsolidation',true,...
    'MaxIterExploration',20,...
    'MaxIterConsolidation',20,...
    'CandidateSpaceConstructor',@CandidateSpacePlanarTrimming,...
    'CandidateSpaceOptions',{'NormalizeGrowthDirection',false,'CheckRedundantTrimmingGrowth',true,'CheckRedundantTrimmingUpdate',true,'CheckDuplicatePointsGrowth',true,'CheckDuplicatePointsUpdate',true},...
    'TrimmingMethodFunction',@component_trimming_method_planar_trimming,...
    'TrimmingMethodOptions',{'NormalizeVariables',true,'TrimmingSlack',0.5},...
    'GrowthRate',0.2,...
    'UseAdaptiveGrowthRate',true,...
    'ApplyLeanness','never',...
    'UsePaddingSamplesInTrimming',true,...
    'UsePreviousEvaluatedSamplesConsolidation',false,...
    'ShapeSamplesUsefulExploration',false,...
    'ShapeSamplesUsefulConsolidation',true,...
    'TrimmingOperationOptions',{'PassesCriterion','full'},...
    'TrimmingOrderOptions',{'OrderPreference','score'});

bottomUpMapping = BottomUpMappingFunction(systemFunction,systemParameter);
designEvaluator = DesignEvaluatorBottomUpMapping(...
    bottomUpMapping,...
    performanceLowerLimit,...
    performanceUpperLimit);
[componentSolutionSpace,optimizationData] = sso_component_stochastic(designEvaluator,...
    initialDesign,designSpaceLowerBound,designSpaceUpperBound,componentIndex,options);
toc(timeElapsedAlgorithm)


%% Plot Comparison
% analytical solution
radiusOptimalAnalytical = performanceUpperLimit*sqrt(2/3);
heightOptimalAnalytical = performanceUpperLimit*sqrt(1/3);
circleAngle = linspace(0,2*pi,100);
analyticalSolutionCircleX = radiusOptimalAnalytical*cos(circleAngle);
analyticalSolutionCircleY = radiusOptimalAnalytical*sin(circleAngle);

% component 1
designSampleComponent = componentSolutionSpace(1).DesignSampleDefinition;
isInsideComponent = componentSolutionSpace(1).IsInsideDefinition;
figure;
plot(designSampleComponent(isInsideComponent,1),designSampleComponent(isInsideComponent,2),'g.');
hold on;
plot(designSampleComponent(~isInsideComponent,1),designSampleComponent(~isInsideComponent,2),'.','color',[227 114 34]/255);
plot(analyticalSolutionCircleX,analyticalSolutionCircleY,'m-','linewidth',2.0);
componentSolutionSpace(1).plot_candidate_space(gcf,'FaceColor','k','EdgeColor','k','FaceAlpha',0.1);
xlabel('x_1');
ylabel('x_3');
axis([-5 +5 -5 +5]);
grid minor;
legend({'Inside Component Space','Outside Component Space','Analytical Solution','Decision Boundary'});
save_print_figure(gcf,[saveFolder,'Component1TrimmingPlot']);

% component 2
designSampleComponent = componentSolutionSpace(2).DesignSampleDefinition;
isInsideComponent = componentSolutionSpace(2).IsInsideDefinition;
figure;
plot(designSampleComponent(isInsideComponent,1),zeros(size(designSampleComponent(isInsideComponent,1))),'g.');
hold on;
plot(designSampleComponent(~isInsideComponent,1),zeros(size(designSampleComponent(~isInsideComponent,1))),'.','color',[227 114 34]/255);
plot([-heightOptimalAnalytical heightOptimalAnalytical],[0 0],'mx','linewidth',2.0);
xlabel('x_2');
axis([-5 +5 -1 +1]);
grid minor;
legend({'Inside Component Space','Outside Component Space','Analytical Solution'});
save_print_figure(gcf,[saveFolder,'Component2TrimmingPlot']);


%%
initialPlotOptions = {'EdgeColor','none','FaceColor',color_palette_tol('blue'),'FaceAlpha',0.2};
grownPlotOptions = {'EdgeColor','none','FaceColor',color_palette_tol('green'),'FaceAlpha',0.2};
trimmedPlotOptions = {'EdgeColor','none','FaceColor',color_palette_tol('red'),'FaceAlpha',0.2};
plot_sso_component_evolution(optimizationData,...
    'InitialPlotOptions',initialPlotOptions,...
    'GrownPlotOptions',grownPlotOptions,...
    'TrimmedPlotOptions',trimmedPlotOptions,...
    'IterationIndices',[17 18 19 20 21 22],...
    'ComponentIndices',1);


%% 
algoData = postprocess_sso_component_stochastic(optimizationData);
plot_sso_component_stochastic_metrics(algoData,'SaveFolder',saveFolder);


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

