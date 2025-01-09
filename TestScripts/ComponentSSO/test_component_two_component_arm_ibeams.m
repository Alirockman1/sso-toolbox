%TEST_COMPONENT_TWO_COMPONENT_ARM_IBEAMS 2-link Robot Arm with I-beams 
%   TEST_COMPONENT_TWO_COMPONENT_ARM_IBEAMS computes component solution spaces for the 
%   two-component arm problem where each component is an I-beam. The solution
%   regions refer to the stifnesses of each component.
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
saveFolder = save_diary_files(mfilename);
goldenRatio = (1+sqrt(5))/2;
figureSize = [goldenRatio 1]*8.5;


%% function call
systemFunction = @two_links_general_beam_stiffness;
systemParameter = {'IBeamVaryingCrossSection_Estimators.mat','SVM_PF','ANN_Mass',true};

% dvbox0 = [1781.396186, 1915.543356;...
%           53919221.423899, 62494590.746367;...
%           41705883.992516, 50083607.271770;...
%           561.355916, 656.869229;...
%           18875688.705177, 22625445.084016;...
%           14036861.194426, 17726738.777123]';
initialDesign = [1710.8 60426e3 42987e3 478.2648 23175e3 7253100];
%                        k1_11   k1_22   k1_44    k2_11    k2_22   k2_44
designSpaceLowerBound = [ 600.0   0.5e7     1e7      200   0.5e7   0.2e7];
designSpaceUpperBound = [3000.0     9e7     8e7     1600     7e7     8e7];
componentIndex = {[1,2,3],[4,5,6]};
%                         displ   m 
performanceLowerLimit = [ -inf -inf];
performanceUpperLimit = [    2  120];
requirementSpacesType = 'Omega2';


%% 
bottomUpMapping = BottomUpMappingFunction(systemFunction,systemParameter);
designEvaluator = DesignEvaluatorBottomUpMapping(bottomUpMapping,performanceLowerLimit,performanceUpperLimit);


%% Box Opt
tic
options = sso_stochastic_options('box',...
    'RequirementSpacesType',requirementSpacesType,...
    'NumberSamplesPerIterationExploration',300,...
    'NumberSamplesPerIterationConsolidation',300,...
    'FixIterNumberExploration',true,...
    'FixIterNumberConsolidation',true,...
    'MaxIterExploration',30,...
    'MaxIterConsolidation',30,...
    'UseAdaptiveGrowthRate',false,...
    'GrowthRate',0.02,...
    'ApplyLeanness','never',...
    'TrimmingOperationOptions',{'PassesCriterion','reduced'},...
    'TrimmingOrderOptions',{'OrderPreference','score'});

[solutionSpaceBox,problemDataBox,iterDataBox] = sso_box_stochastic(designEvaluator,...
    initialDesign,designSpaceLowerBound,designSpaceUpperBound,options);
toc


%% Component Opt
tic
options = sso_stochastic_options('component',...
    'RequirementSpacesType',requirementSpacesType,...
    'NumberSamplesPerIterationExploration',300,...
    'NumberSamplesPerIterationConsolidation',300,...
    'FixIterNumberExploration',true,...
    'FixIterNumberConsolidation',true,...
    'MaxIterExploration',30,...
    'MaxIterConsolidation',30,...
    'CandidateSpaceConstructorExploration',@CandidateSpaceConvexHull,...
    'CandidateSpaceConstructorConsolidation',@CandidateSpaceConvexHull,...
    'TrimmingMethodFunction',@component_trimming_method_planar_trimming,...
    'UseAdaptiveGrowthRate',false,...
    'GrowthRate',0.02,...
    'ApplyLeanness','never',...
    'UsePaddingSamplesInTrimming',true,...
    'UsePreviousEvaluatedSamplesConsolidation',false,...
    'UsePreviousPaddingSamplesConsolidation',false,...
    'TrimmingOperationOptions',{'PassesCriterion','reduced'},...
    'TrimmingOrderOptions',{'OrderPreference','score'});

[componentSolutionSpace,problemDataComponent,iterDataComponent] = sso_component_stochastic(designEvaluator,...
    solutionSpaceBox,designSpaceLowerBound,designSpaceUpperBound,componentIndex,options);
toc

%% Selective Design Space Projection
% desired plots by design variable index
plotTiles = [2,3]; % 2 rows, 3 columns of subplots
desiredPlots = [...
    1,2;...
    1,3;...
    2,3;...
    4,5;...
    4,6;...
    5,6];
axisLabel = {'k_{(1)11}','k_{(1)22}','k_{(1)44}','k_{(2)11}','k_{(2)22}','k_{(2)44}'};

[h,problemData,plotData] = plot_selective_design_space_projection(...
    designEvaluator,...
    solutionSpaceBox,...
    designSpaceLowerBound,...
    designSpaceUpperBound,...
    desiredPlots,...
    plotTiles,...
    'NumberSamplesPerPlot',10000,...
    'AxesLabels',axisLabel,...
    'MarkerColorsViolatedRequirements',{'r','b'},...
    'PlotOptionsBad',{'Marker','x'},...
    'PlotOptionsPhysicallyInfeasible',{'Marker','*'});
legend([h.GoodPerformance,h.PhysicallyInfeasible,h.BadPerformance],...
    {'Good Designs','Physically Infeasible Designs','Violate Displacement Requirement','Violate Mass Requirement'});
sgtitle('Two Component Arm - I-beams with Varying Cross Sections');
save_print_figure(h.MainFigure,[saveFolder,'SelectiveDesignSpaceProjection']);
clear h


%% Plot Comparison
volComponent1 = componentSolutionSpace(1).Measure;
volComponent2 = componentSolutionSpace(2).Measure;
volBox1 = prod(solutionSpaceBox(2,1:3)-solutionSpaceBox(1,1:3));
volBox2 = prod(solutionSpaceBox(2,4:6)-solutionSpaceBox(1,4:6));

% Original Solution Box + Convex Hull - Component 1
volume_string = sprintf(': \\mu(\\Omega)=%g ; Relative Increase \\approx %.2fx',volComponent1,volComponent1/volBox1);
figure;
plot_design_box_3d(gcf,solutionSpaceBox(:,1:3),'FaceAlpha',0.1,'FaceColor',[0 0 0]);
hold on; 
grid minor;
componentSolutionSpace(1).plot_candidate_space(gcf,'FaceColor','g','FaceAlpha',0.1,'EdgeColor','g');
% % Text
title(['Component 1',volume_string],'FontSize',14);
xlabel('k_{(1)11}','FontSize',16);
ylabel('k_{(1)22}','FontSize',16);
zlabel('k_{(1)44}','FontSize',16);
legend({'Original Solution Box','Component Solution Space'},'FontSize',14);
save_print_figure(gcf,[saveFolder,'ConvexHullBoxComponent1']);
clear ans

% Original Solution Box + Convex Hull - Component 2
volume_string = sprintf(': \\mu(\\Omega)=%g ; Relative Increase \\approx %.2fx',volComponent2,volComponent2/volBox2);
figure;
plot_design_box_3d(gcf,solutionSpaceBox(:,4:6),'FaceAlpha',0.1,'FaceColor',[0 0 0]);
hold on; 
grid minor;
% % Text
title(['Component 2',volume_string],'FontSize',14);
xlabel('k_{(2)11}','FontSize',16);
ylabel('k_{(2)22}','FontSize',16);
zlabel('k_{(2)44}','FontSize',16);
componentSolutionSpace(2).plot_candidate_space(gcf,'FaceColor','g','FaceAlpha',0.1,'EdgeColor','g');
legend({'Original Solution Box','Component Solution Space'},'FontSize',14);
save_print_figure(gcf,[saveFolder,'ConvexHullBoxComponent2']);
clear ans


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

