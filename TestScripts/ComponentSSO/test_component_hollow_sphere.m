%test_component_hollow_sphere Component solution spaces for a sphere problem 
%   test_component_hollow_sphere computes a component solution spaces with 
%   for a sphere problem.
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
RNGstate = rng;
savefolder = save_diary_files(mfilename);
goldenratio = (1+sqrt(5))/2;
figure_size = [goldenratio 1]*8.5;


%% function call
%
systemFunction = @distance_to_center;
systemParameter = [0,0,0];
%        x1  x2  x3
designSpaceLowerBound = [-6 -6 -6];
designSpaceUpperBound = [ 6  6  6];
Components = {[1,3],[2]};
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
    'CandidateSpaceConstructor',@CandidateSpaceCornerBoxRemoval,...
    'TrimmingMethodFunction',@component_trimming_method_corner_box_removal,...
    ... 'TrimmingMethodOptions',{'ReferenceDesigns','boundary-center'},...
    ... 'TrimmingMethodOptions',{'CornersToTest','away'},...
    ... 'TrimmingComponentChoiceOptions',{'WeightedCostType','ComponentDimension'},...
    'GrowthRate',0.2,...
    'UseAdaptiveGrowthRate',true,...
    'ApplyLeanness','never',...
    'UsePaddingSamplesInTrimming',true,...
    'UsePreviousEvaluatedSamplesConsolidation',true,...
    'ShapeSamplesUsefulExploration',false,...
    'ShapeSamplesUsefulConsolidation',false,...
    'TrimmingOperationOptions',{'PassesCriterion','full'},...
    'TrimmingOrderOptions',{'OrderPreference','score'},...
    'NumberPaddingSamples',1000);

bottomUpMapping = BottomUpMappingFunction(systemFunction,systemParameter);
designEvaluator = DesignEvaluatorBottomUpMapping(...
    bottomUpMapping,...
    performanceLowerLimit,...
    performanceUpperLimit);
[componentSolutionSpace,problemData,iterData] = sso_component_stochastic(designEvaluator,...
    initialDesign,designSpaceLowerBound,designSpaceUpperBound,Components,options);
toc(timeElapsedAlgorithm)


%% Plot Comparison
% analytical solution
r_maxtheory = performanceUpperLimit*sqrt(2/3);
h_maxtheory = performanceUpperLimit*sqrt(1/3);
circleAngle = linspace(0,2*pi,100);
anasol_1_x = r_maxtheory*cos(circleAngle);
anasol_1_y = r_maxtheory*sin(circleAngle);

% component 1
x_css = componentSolutionSpace(1).DesignSampleDefinition;
y_css = componentSolutionSpace(1).IsInsideDefinition;
figure;
plot(x_css(y_css,1),x_css(y_css,2),'g.');
hold on;
plot(x_css(~y_css,1),x_css(~y_css,2),'.','color',[227 114 34]/255);
plot(anasol_1_x,anasol_1_y,'m-','linewidth',2.0);
componentSolutionSpace(1).plot_candidate_space(gcf,'FaceColor','k','EdgeColor','k','FaceAlpha',0.1);
xlabel('x_1');
ylabel('x_3');
axis([-5 +5 -5 +5]);
grid minor;
legend({'Inside Component Space','Outside Component Space','Analytical Solution','Decision Boundary'});
save_print_figure(gcf,[savefolder,'Component1TrimmingPlot']);

% component 2
x_css = componentSolutionSpace(2).DesignSampleDefinition;
y_css = componentSolutionSpace(2).IsInsideDefinition;
figure;
plot(x_css(y_css,1),zeros(size(x_css(y_css,1))),'g.');
hold on;
plot(x_css(~y_css,1),zeros(size(x_css(~y_css,1))),'.','color',[227 114 34]/255);
plot([-h_maxtheory h_maxtheory],[0 0],'mx','linewidth',2.0);
xlabel('x_2');
axis([-5 +5 -1 +1]);
grid minor;
legend({'Inside Component Space','Outside Component Space','Analytical Solution'});
save_print_figure(gcf,[savefolder,'Component2TrimmingPlot']);


%% 
algoData = postprocess_sso_component_stochastic(problemData,iterData);
plot_sso_component_stochastic_metrics(algoData,'SaveFolder',savefolder);


%% Save and Stop Transcripting
save([savefolder,'Data.mat']);
diary off;

