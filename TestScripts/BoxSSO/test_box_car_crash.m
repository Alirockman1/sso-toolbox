%TEST_BOX_CAR_CRASH 2D Car Crash example for box-shaped solution spaces
%   TEST_BOX_CAR_CRASH allows to test the box SSO algorithm with a problem with
%   a known analytical solution, namely the simplifed 2D car crash problem.
%   Both solution and algorithm performance metrics are plotted at the end.
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


%% problem setup
bottomUpMapping = @car_crash_2d;
%               m d1c d2c   v0
parameter = [2000 0.3 0.3 15.6];
%                       E_rem    a_max order
performanceLowerLimit = [-inf        0  -inf];
performanceUpperLimit = [   0  32*9.81     0];

bottomUpMapping = BottomUpMappingFunction(bottomUpMapping,parameter);
designEvaluator = DesignEvaluatorBottomUpMapping(bottomUpMapping,performanceLowerLimit,performanceUpperLimit);


%% design space
%                            F1     F2   
designSpaceLowerBound = [   1e5    1e5];
designSpaceUpperBound = [   1e6    1e6];
startingPoint         = [4.05e5 4.05e5];


%% optimization call
options = sso_stochastic_options('box',...
    'NumberSamplesPerIteration',100,...
    'GrowthRate',0.2,...
    'FixIterNumberExploration',true,...
    'FixIterNumberConsolidation',true,...
    'UseAdaptiveGrowthRate',false);
[designBox,optimizationData] = sso_box_stochastic(designEvaluator,startingPoint,...
    designSpaceLowerBound,designSpaceUpperBound,options);


%% Plot Solution
figure;
goodTriangleX = [405000 627840 183000 405000];
goodTriangleY = [405000 627840 627840 405000];
plot(goodTriangleX,goodTriangleY,'g');
hold on;
grid minor;
plot_design_box_2d(gcf,designBox);
plot(startingPoint(1),startingPoint(2),'rx');
analyticalSolutionX = [291e3 516e3 516e3 291e3 291e3];
analyticalSolutionY = [516e3 516e3 628e3 628e3 516e3];
plot(analyticalSolutionX,analyticalSolutionY,'mo');
xlabel('F_1 [N]');
ylabel('F_2 [N]');
save_print_figure(gcf,[saveFolder,'SolutionBox'],'Size',figureSize);


%% Performance Metrics
algoData = postprocess_sso_box_stochastic(optimizationData);
plot_sso_box_stochastic_metrics(algoData,...
    'SaveFolder',saveFolder,...
    'SaveFigureOptions',{'Size',figureSize});


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

