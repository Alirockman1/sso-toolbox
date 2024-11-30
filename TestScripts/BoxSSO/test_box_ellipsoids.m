%test_box_ellipsoids Test requirement spaces computation with two ellipsoids
%   test_box_ellipsoids uses a bottom-up mapping with two defined ellipsoids, 
%   one which describes a region where performance is satisfied and another 
%   which describes where designs are physically feasible, to test the 
%   computation of requirement spaces with the stochastic algorithm.
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
close all;
fclose all;
clear all;
clc;
more off;
diary off;


%% debugging
rng(6);


%% Documentation / Archive
rngState = rng;
saveFolder = save_diary_files(mfilename);
goldenRatio = (1+sqrt(5))/2;
figureSize = [goldenRatio 1]*8.5;


%% function call
systemFunction = @two_ellipses_requirement_space;

designSpaceLowerBound = [-2 -2];
designSpaceUpperBound = [14 12];
initialPoint = [10 6];

performanceLowerLimit = -inf;
performanceUpperLimit = 1;
physicalFeasibilityLowerLimit = -inf;
physicalFeasibilityUpperLimit = 1;

requirementSpacesType = 'Omega2';
options = sso_stochastic_options('box',...
    'RequirementSpacesType',requirementSpacesType,...
    'ApplyLeanness','end-only',...
    'NumberSamplesPerIteration',200,...
    'FixIterNumberExploration',true,...
    'FixIterNumberConsolidation',true,...
    'UseAdaptiveGrowthRate',true);

bottomUpMapping = BottomUpMappingFunction(systemFunction);
designEvaluator = DesignEvaluatorBottomUpMapping(...
    bottomUpMapping,...
    performanceLowerLimit,...
    performanceUpperLimit,...
    'PhysicalFeasibilityLowerLimit',physicalFeasibilityLowerLimit,...
    'PhysicalFeasibilityUpperLimit',physicalFeasibilityUpperLimit);

[designBox,problemData,iterData] = sso_box_stochastic(designEvaluator,initialPoint,...
    designSpaceLowerBound,designSpaceUpperBound,options);


%% Plot Solution
plot_ellipse_2d(figure,[11,5],[7,2],1.8,'PatchOptions',{'FaceColor','none','EdgeColor','g'});
hold all;
plot_ellipse_2d(gcf,[5,6],[7,2],0.2,'PatchOptions',{'FaceColor','none','EdgeColor','b'});
grid minor;
hold on;
plot_design_box_2d(gcf,designBox);
plot(initialPoint(1),initialPoint(2),'rx');
xlabel('x_1');
ylabel('x_2');
save_print_figure(gcf,[saveFolder,'BasicSolutionBox'],'Size',figureSize);


%% Performance Metrics
algoData = postprocess_sso_box_stochastic(problemData,iterData);
plot_sso_box_stochastic_metrics(algoData,...
    'SaveFolder',saveFolder,...
    'SaveFigureOptions',{'Size',figureSize});


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

