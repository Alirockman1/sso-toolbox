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
rngState = rng;
saveFolder = save_diary_files(mfilename);
goldenRatio = (1+sqrt(5))/2;
figureSize = [goldenRatio 1]*8.5;


%%
trussAnalysisChoice = '16-DoF';
computeDisplacement = true;
computeMass = true;
computeDisplacementAndMass = true;


%% function call
systemFunction = @truss_generic_2d_moving_node;

switch trussAnalysisChoice
    case 'six-bar-one-moving'
        systemParameter.ElementCrossSectionArea = 10; % [mm^2]
        systemParameter.ElementYoungsModulus = 210e3; % [MPa]
        systemParameter.ElementDensity = 7850e-9; % [kg/mm^3]
        systemParameter.BaseNodePosition = [...
            0 0;  ...
            1 0; ...
            2 0.5; ...
            nan nan; ...
            0 1];
        systemParameter.FixedDegreesOfFreedom = [...
            true true; ...
            false false; ...
            false false; ...
            false false; ...
            true true];
        systemParameter.NodeForce = [...
            0 0; ...
            0 0; ...
            0 -1000; ...
            0 0; ...
            0 0];
        systemParameter.NodeElement = [...
            1 2; ...
            2 3; ...
            3 4; ...
            4 5; ...
            2 5; ...
            2 4];
        initialDesign = [1,1];
        designSpaceLowerBoundDisplacement = [0.9 0.5];
        designSpaceUpperBoundDisplacement = [2.0 1.5];
        designSpaceLowerBoundMass = [0 -0.5];
        designSpaceUpperBoundMass = [2  1.5];
    case 'six-bar-two-moving'
        systemParameter.ElementCrossSectionArea = 10; % [mm^2]
        systemParameter.ElementYoungsModulus = 210e3; % [MPa]
        systemParameter.ElementDensity = 7850e-9; % [kg/mm^3]
        systemParameter.BaseNodePosition = [...
            0 0;  ...
            nan nan; ...
            2 0.5; ...
            nan nan; ...
            0 1];
        systemParameter.FixedDegreesOfFreedom = [...
            true true; ...
            false false; ...
            false false; ...
            false false; ...
            true true];
        systemParameter.NodeForce = [...
            0 0; ...
            0 0; ...
            0 -1000; ...
            0 0; ...
            0 0];
        systemParameter.NodeElement = [...
            1 2; ...
            2 3; ...
            3 4; ...
            4 5; ...
            2 5; ...
            2 4];
        initialDesign = [1,0,1,1];
        designSpaceLowerBoundDisplacement = [0.1 -0.5 0.9 0.5];
        designSpaceUpperBoundDisplacement = [1.5  0.5 2.0 1.5];
        designSpaceLowerBoundMass = [0 -0.5 0 -0.5];
        designSpaceUpperBoundMass = [2  1.5 2  1.5];
    case '16-DoF'
        systemParameter.ElementCrossSectionArea = 10; % [mm^2]
        systemParameter.ElementYoungsModulus = 210e3; % [MPa]
        systemParameter.ElementDensity = 7850e-9; % [kg/mm^3]
        systemParameter.BaseNodePosition = [...
	        nan nan; ... % (1)
            nan nan; ... % (2)
            nan nan; ... % (3)
            nan nan; ... % (4)
            nan nan; ... % (5)
            nan nan; ... % (6)
            nan nan; ... % (7)
            nan nan; ... % (8)
              5 0.5; ... % (9)
              0   0; ... % (10)
              0   1];    % (11)
        systemParameter.FixedDegreesOfFreedom = [...
	        false false; ... % (1)
	        false false; ... % (2)
	        false false; ... % (3)
	        false false; ... % (4)
	        false false; ... % (5)
	        false false; ... % (6)
            false false; ... % (7)
	        false false; ... % (8)
	        false false; ... % (9)
	         true  true;...  % (10)
             true  true];    % (11)
        systemParameter.NodeForce = [...
	        0 0; ... % (1)
	        0 0; ... % (2)
            0 0; ... % (3)
	        0 0; ... % (4)
            0 0; ... % (5)
	        0 0; ... % (6)
            0 0; ... % (7)
	        0 0; ... % (8)
	        0 -1000; ... % (9)
	        0 0; ... % (10)
	        0 0]; % (11) assumed [N]
        systemParameter.NodeElement = [...
	        11 1; ... % (1)
            11 2; ... % (2)
            10 2; ... % (3)
            1 2; ... % (4)
            1 4; ... % (5)
            1 3; ... % (6)
            2 4; ... % (7)
            3 4; ... % (8)
            3 6; ... % (9)
            3 5; ... % (10)
            4 6; ... % (11)
            5 6; ... % (12)
            5 8; ... % (13)
            5 7; ... % (14)
            6 8; ... % (15)
            7 8; ... % (16)
            7 9; ... % (17)
            8 9]; % (18)
        initialDesign = [...
            1, 1, ...
            1, 0, ...
            2, 1, ...
            2, 0, ...
            3, 1, ...
            3, 0, ...
            4, 1, ...
            4, 0];
        designSpaceLowerBoundDisplacement = initialDesign-0.9;
        designSpaceUpperBoundDisplacement = initialDesign+0.9;
        designSpaceLowerBoundMass = repmat([0 -0.5],1,8);
        designSpaceUpperBoundMass = repmat([10 1.5],1,8);
    case '36-DoF'
        systemParameter.ElementCrossSectionArea = 10; % [mm^2]
        systemParameter.ElementYoungsModulus = 210e3; % [MPa]
        systemParameter.ElementDensity = 7850e-9; % [kg/mm^3]
        systemParameter.BaseNodePosition = [...
            nan nan; ... % (1)
            nan nan; ... % (2)
            nan nan; ... % (3)
            nan nan; ... % (4)
            nan nan; ... % (5)
            nan nan; ... % (6)
            nan nan; ... % (7)
            nan nan; ... % (8)
            nan nan; ... % (9)
            nan nan; ... % (10)
            nan nan; ... % (11)
            nan nan; ... % (12)
            nan nan; ... % (13)
            nan nan; ... % (14)
            nan nan; ... % (15)
            nan nan; ... % (16)
            nan nan; ... % (17)
            nan nan; ... % (18)
             10 0.5; ... % (19)
              0   0; ... % (20)
              0   1];    % (21)
        systemParameter.FixedDegreesOfFreedom = [...
            false false; ... % (1)
            false false; ... % (2)
            false false; ... % (3)
            false false; ... % (4)
            false false; ... % (5)
            false false; ... % (6)
            false false; ... % (7)
            false false; ... % (8)
            false false; ... % (9)
            false false; ... % (10)
            false false; ... % (11)
            false false; ... % (12)
            false false; ... % (13)
            false false; ... % (14)
            false false; ... % (15)
            false false; ... % (16)
            false false; ... % (17)
            false false; ... % (18)
            false false; ... % (19)
             true  true;...  % (20)
             true  true];    % (21)
        systemParameter.NodeForce = [...
            0 0; ... % (1)
            0 0; ... % (2)
            0 0; ... % (3)
            0 0; ... % (4)
            0 0; ... % (5)
            0 0; ... % (6)
            0 0; ... % (7)
            0 0; ... % (8)
            0 0; ... % (9)
            0 0; ... % (10)
            0 0; ... % (11)
            0 0; ... % (12)
            0 0; ... % (13)
            0 0; ... % (14)
            0 0; ... % (15)
            0 0; ... % (16)
            0 0; ... % (17)
            0 0; ... % (18)
            0 -1000; ... % (19)
            0 0; ... % (20)
            0 0]; % (21) [N]
        systemParameter.NodeElement = [...
            21 1; ... % (1) -
            21 2; ... % (2) \ 
            20 2; ... % (3) -
            17 19; ... % (4) \
            18 19; ... % (5) /
            1 2; ... % (6) |
            1 4; ... % (7) \
            1 3; ... % (8) -
            2 4; ... % (9) -
            3 4; ... % (10) | 
            3 6; ... % (11) \
            3 5; ... % (12) -
            4 6; ... % (13) -
            5 6; ... % (14) |
            5 8; ... % (15) \ 
            5 7; ... % (16) -
            6 8; ... % (17) -
            7 8; ... % (18) |
            7 10; ... % (19) \ 
            7 9; ... % (20) - 
            8 10; ... % (21) -
            9 10; ... % (22) |
            9 12; ... % (23) \ 
            9 11; ... % (24) -
            10 12; ... % (25) - 
            11 12; ... % (26) |
            11 14; ... % (27) \ 
            11 13; ... % (28) -
            12 14; ... % (29) - 
            13 14; ... % (30) |
            13 16; ... % (31) \
            13 15; ... % (32) -
            14 16; ... % (33) -
            15 16; ... % (34) |
            15 18; ... % (35) \ 
            15 17; ... % (36) -
            16 18; ... % (37) - 
            17 18; ... % (38) |
            ];
        %   x  y
        initialDesign = [...
            1, 1, ...
            1, 0, ...
            2, 1, ...
            2, 0, ...
            3, 1, ...
            3, 0, ...
            4, 1, ...
            4, 0, ...
            5, 1, ...
            5, 0, ...
            6, 1, ...
            6, 0, ...
            7, 1, ...
            7, 0, ...
            8, 1, ...
            8, 0, ...
            9, 1, ...
            9, 0];
        designSpaceLowerBoundDisplacement = initialDesign-0.9;
        designSpaceUpperBoundDisplacement = initialDesign+0.9;
        designSpaceLowerBoundMass = repmat([0 -0.5],1,18);
        designSpaceUpperBoundMass = repmat([10 1.5],1,18);
end


%% starting process - compute necessary parts, visualize initial truss
% component generation - assume all position variables of each node form a component
isDesignVariable = isnan(systemParameter.BaseNodePosition);
nDimension = size(systemParameter.BaseNodePosition,2);
nComponent = sum(isDesignVariable(:,1),1);
componentIndex = cell(1,nComponent);
for i=1:nComponent
    componentIndex{i} = (1 + nDimension*(i-1) + [0:(nDimension-1)])';
end

% performance limits
nElement = size(systemParameter.NodeElement,1);
performanceLowerLimit = -inf;
performanceUpperLimit = [nan nan repmat(inf,1,nElement)];

% initial truss
plot_results_truss_generic_2d_moving_node(systemParameter,...
    'NodePositionInitial',initialDesign);
save_print_figure(gcf,[saveFolder,'InitialTruss']);

% deformed trusses
nodePositionInitial = systemParameter.BaseNodePosition;
nodePositionInitial(isDesignVariable) = column_vector_to_row_major_matrix(initialDesign',nDimension);
nodeDisplacementInitial = ...
    truss_analysis(...
	    nodePositionInitial,...
	    systemParameter.FixedDegreesOfFreedom,...
	    systemParameter.NodeForce,...
	    systemParameter.NodeElement,...
	    systemParameter.ElementCrossSectionArea,...
	    systemParameter.ElementYoungsModulus);

plot_results_truss_generic_2d_moving_node(systemParameter,...
    'NodePositionInitial',initialDesign,...
    'DeformationInitial',nodeDisplacementInitial);
save_print_figure(gcf,[saveFolder,'InitialTrussDeformation']);


%% find optimum
bottomUpMapping = BottomUpMappingFunction(systemFunction,'SystemParameter',systemParameter);

if(computeDisplacement)
    warning('off');
    [nodePositionOptimalDisplacement,displacementOptimal] = design_optimize_quantities_of_interest(...
        bottomUpMapping,...
        initialDesign,...
        designSpaceLowerBoundDisplacement,...
        designSpaceUpperBoundDisplacement,...
        @(performanceMeasure)[performanceMeasure(1)],...
        ...'InequalityConstraintFunction',@(performanceMeasure)[performanceMeasure(3:end)-performanceUpperLimit(3:end)],...
        'InequalityConstraintFunction',@(performanceMeasure)[-performanceMeasure(1)],...
        'OptimizationMethodFunction',@optimization_ga_wrapper,...
        'OptimizationMethodOptions',{'Display','diagnose'});
    warning('on');
end

if(computeMass)
    warning('off');
    [nodePositionOptimalMass,massOptimal] = design_optimize_quantities_of_interest(...
        bottomUpMapping,...
        initialDesign,...
        designSpaceLowerBoundMass,...
        designSpaceUpperBoundMass,...
        @(performanceMeasure)[performanceMeasure(2)],...
        ...'InequalityConstraintFunction',@(performanceMeasure)[performanceMeasure(3:end)-performanceUpperLimit(3:end)],...
        'OptimizationMethodFunction',@optimization_ga_wrapper,...
        'OptimizationMethodOptions',{'Display','diagnose'});
    warning('on');
end

if(computeDisplacement)
    % initial truss + optimized truss
    plot_results_truss_generic_2d_moving_node(systemParameter,...
        'NodePositionInitial',initialDesign,...
        'NodePositionOptimalDisplacement',nodePositionOptimalDisplacement);
    save_print_figure(gcf,[saveFolder,'InitialOptimizedTruss']);
    
    if(computeMass)
        % initial truss + optimized displacment truss + optimized mass truss
        plot_results_truss_generic_2d_moving_node(systemParameter,...
            'NodePositionInitial',initialDesign,...
            'NodePositionOptimalDisplacement',nodePositionOptimalDisplacement,...
            'NodePositionOptimalMass',nodePositionOptimalMass);
        save_print_figure(gcf,[saveFolder,'InitialOptimizedTrussWithMass']);
    end
end

if(computeDisplacement)
    nodePositionOptimized = systemParameter.BaseNodePosition;
    nodePositionOptimized(isDesignVariable) = column_vector_to_row_major_matrix(nodePositionOptimalDisplacement',nDimension);
    nodeDisplacementOptimal = ...
        truss_analysis(...
	        nodePositionOptimized,...
	        systemParameter.FixedDegreesOfFreedom,...
	        systemParameter.NodeForce,...
	        systemParameter.NodeElement,...
	        systemParameter.ElementCrossSectionArea,...
	        systemParameter.ElementYoungsModulus);

    plot_results_truss_generic_2d_moving_node(systemParameter,...
        'NodePositionInitial',initialDesign,...
        'NodePositionOptimalDisplacement',nodePositionOptimalDisplacement,...
        'DeformationInitial',nodeDisplacementInitial,...
        'DeformationOptimalDisplacement',nodeDisplacementOptimal);
    save_print_figure(gcf,[saveFolder,'InitialOptimizedTrussDeformation']);
end


%% establish upper performance limits
% update uppwer limit based on either optimal value or initial value
performanceMeasureInitial = bottomUpMapping.response(initialDesign);
performanceUpperLimitDisplacement = performanceUpperLimit;
performanceUpperLimitMass = performanceUpperLimit;
performanceUpperLimitDisplacementAndMass = performanceUpperLimit;

% displacement
if(computeDisplacement)
    performanceUpperLimitDisplacement(1) = performanceMeasureInitial(1);
    designEvaluatorDisplacement = DesignEvaluatorBottomUpMapping(...
        bottomUpMapping,...
        performanceLowerLimit,...
        performanceUpperLimitDisplacement);
end

% mass
if(computeMass)
    performanceUpperLimitMass(2) = performanceMeasureInitial(2);
    designEvaluatorMass = DesignEvaluatorBottomUpMapping(...
        bottomUpMapping,...
        performanceLowerLimit,...
        performanceUpperLimitMass);
end

if(computeDisplacementAndMass)
    performanceUpperLimitDisplacementAndMass(1) = performanceMeasureInitial(1);
    performanceUpperLimitDisplacementAndMass(2) = performanceMeasureInitial(2);
    designEvaluatorDisplacementAndMass = DesignEvaluatorBottomUpMapping(...
        bottomUpMapping,...
        performanceLowerLimit,...
        performanceUpperLimitDisplacementAndMass);
end


%% Box Opt
if(computeDisplacement)
    timeElapsedBox = tic;
    optionsBox = sso_stochastic_options('box',...
        'NumberSamplesPerIterationExploration',300,...
        'NumberSamplesPerIterationConsolidation',300,...
        'FixIterNumberExploration',true,...
        'FixIterNumberConsolidation',true,...
        'MaxIterExploration',30,...
        'MaxIterConsolidation',30,...
        'UseAdaptiveGrowthRate',false,...
        'GrowthRate',0.2,...
        'ApplyLeanness','never',...
        'TrimmingOperationOptions',{'PassesCriterion','reduced'},...
        'TrimmingOrderOptions',{'OrderPreference','score'});
    
    rng(rngState);
    [solutionSpaceBoxDisplacement,problemDataBoxDisplacement,iterDataBoxDisplacement] = sso_box_stochastic(designEvaluatorDisplacement,...
        initialDesign,designSpaceLowerBoundDisplacement,designSpaceUpperBoundDisplacement,optionsBox);
    toc(timeElapsedBox)
end

if(computeMass)
    timeElapsedBox = tic;
    optionsBox = sso_stochastic_options('box',...
        'NumberSamplesPerIterationExploration',300,...
        'NumberSamplesPerIterationConsolidation',300,...
        'FixIterNumberExploration',true,...
        'FixIterNumberConsolidation',true,...
        'MaxIterExploration',30,...
        'MaxIterConsolidation',30,...
        'UseAdaptiveGrowthRate',false,...
        'GrowthRate',0.05,...
        'ApplyLeanness','never',...
        'TrimmingOperationOptions',{'PassesCriterion','reduced'},...
        'TrimmingOrderOptions',{'OrderPreference','score'});
    
    rng(rngState);
    [solutionSpaceBoxMass,problemDataBoxMass,iterDataBoxMass] = sso_box_stochastic(designEvaluatorMass,...
        initialDesign,designSpaceLowerBoundMass,designSpaceUpperBoundMass,optionsBox);
    toc(timeElapsedBox)
end

if(computeDisplacementAndMass)
    timeElapsedBox = tic;
    optionsBox = sso_stochastic_options('box',...
        'NumberSamplesPerIterationExploration',300,...
        'NumberSamplesPerIterationConsolidation',300,...
        'FixIterNumberExploration',true,...
        'FixIterNumberConsolidation',true,...
        'MaxIterExploration',30,...
        'MaxIterConsolidation',30,...
        'UseAdaptiveGrowthRate',false,...
        'GrowthRate',0.2,...
        'ApplyLeanness','never',...
        'TrimmingOperationOptions',{'PassesCriterion','reduced'},...
        'TrimmingOrderOptions',{'OrderPreference','score'});
    
    rng(rngState);
    [solutionSpaceBoxDisplacementAndMass,problemDataBoxDisplacementAndMass,iterDataBoxDisplacementAndMass] = ...
        sso_box_stochastic(designEvaluatorDisplacementAndMass,...
        initialDesign,designSpaceLowerBoundDisplacement,designSpaceUpperBoundDisplacement,optionsBox);
    toc(timeElapsedBox)
end


%% Component Opt
if(computeDisplacement)
    timeElapsedComponent = tic;
    optionsComponent = sso_stochastic_options('component',...
        'NumberSamplesPerIterationExploration',300,...
        'NumberSamplesPerIterationConsolidation',300,...
        'FixIterNumberExploration',true,...
        'FixIterNumberConsolidation',true,...
        'MaxIterExploration',30,...
        'MaxIterConsolidation',30,...
        'CandidateSpaceConstructorExploration',@CandidateSpaceConvexHull,...
        'CandidateSpaceConstructorConsolidation',@CandidateSpaceConvexHull,...
        'TrimmingMethodFunction',@component_trimming_method_planar_trimming,...
        ... 'TrimmingMethodOptions',{'ReferenceDesigns','boundary-center'},...
        'UseAdaptiveGrowthRate',false,...
        'GrowthRate',0.2,...
        'ApplyLeanness','never',...
        'UsePaddingSamplesInTrimming',true,...
        'UsePreviousEvaluatedSamplesConsolidation',false,...
        'UsePreviousPaddingSamplesConsolidation',false,...
        'TrimmingOperationOptions',{'PassesCriterion','reduced'},...
        'TrimmingOrderOptions',{'OrderPreference','score'});
    
    rng(rngState);
    [componentSolutionSpaceDisplacement,problemDataComponentDisplacement,iterDataComponentDisplacement] = sso_component_stochastic(designEvaluatorDisplacement,...
        initialDesign,designSpaceLowerBoundDisplacement,designSpaceUpperBoundDisplacement,componentIndex,optionsComponent);
    toc(timeElapsedComponent)
end

if(computeMass)
    timeElapsedComponent = tic;
    optionsComponent = sso_stochastic_options('component',...
        'NumberSamplesPerIterationExploration',300,...
        'NumberSamplesPerIterationConsolidation',300,...
        'FixIterNumberExploration',true,...
        'FixIterNumberConsolidation',true,...
        'MaxIterExploration',30,...
        'MaxIterConsolidation',30,...
        'CandidateSpaceConstructorExploration',@CandidateSpaceConvexHull,...
        'CandidateSpaceConstructorConsolidation',@CandidateSpaceConvexHull,...
        'TrimmingMethodFunction',@component_trimming_method_planar_trimming,...
        ... 'TrimmingMethodOptions',{'ReferenceDesigns','boundary-center'},...
        'UseAdaptiveGrowthRate',false,...
        'GrowthRate',0.04,...
        'ApplyLeanness','never',...
        'UsePaddingSamplesInTrimming',true,...
        'UsePreviousEvaluatedSamplesConsolidation',false,...
        'UsePreviousPaddingSamplesConsolidation',false,...
        'TrimmingOperationOptions',{'PassesCriterion','reduced'},...
        'TrimmingOrderOptions',{'OrderPreference','score'});
    
    rng(rngState);
    [componentSolutionSpaceMass,problemDataComponentMass,iterDataComponentMass] = sso_component_stochastic(designEvaluatorMass,...
        initialDesign,designSpaceLowerBoundMass,designSpaceUpperBoundMass,componentIndex,optionsComponent);
    toc(timeElapsedComponent)
end

if(computeDisplacementAndMass)
    timeElapsedComponent = tic;
    optionsComponent = sso_stochastic_options('component',...
        'NumberSamplesPerIterationExploration',300,...
        'NumberSamplesPerIterationConsolidation',300,...
        'FixIterNumberExploration',true,...
        'FixIterNumberConsolidation',true,...
        'MaxIterExploration',30,...
        'MaxIterConsolidation',30,...
        'CandidateSpaceConstructorExploration',@CandidateSpaceConvexHull,...
        'CandidateSpaceConstructorConsolidation',@CandidateSpaceConvexHull,...
        'TrimmingMethodFunction',@component_trimming_method_planar_trimming,...
        ... 'TrimmingMethodOptions',{'ReferenceDesigns','boundary-center'},...
        'UseAdaptiveGrowthRate',false,...
        'GrowthRate',0.2,...
        'ApplyLeanness','never',...
        'UsePaddingSamplesInTrimming',true,...
        'UsePreviousEvaluatedSamplesConsolidation',false,...
        'UsePreviousPaddingSamplesConsolidation',false,...
        'TrimmingOperationOptions',{'PassesCriterion','reduced'},...
        'TrimmingOrderOptions',{'OrderPreference','score'});
    
    rng(rngState);
    [componentSolutionSpaceDisplacementAndMass,problemDataComponentDisplacementAndMass,iterDataComponentDisplacementAndMass] = ...
        sso_component_stochastic(designEvaluatorDisplacementAndMass,...
        initialDesign,designSpaceLowerBoundDisplacement,designSpaceUpperBoundDisplacement,componentIndex,optionsComponent);
    toc(timeElapsedComponent)
end


%% measure comparison



%% Plot Visualization
if(computeDisplacement)
    % initial truss + component solution spaces
    plot_results_truss_generic_2d_moving_node(systemParameter,...
        'NodePositionInitial',initialDesign,...
        'ComponentSolutionSpaceDisplacement',componentSolutionSpaceDisplacement);
    save_print_figure(gcf,[saveFolder,'InitialTrussComponent']);
    
    % initial truss + optimized truss + component solution spaces
    plot_results_truss_generic_2d_moving_node(systemParameter,...
        'NodePositionInitial',initialDesign,...
        'NodePositionOptimalDisplacement',nodePositionOptimalDisplacement,...
        'ComponentSolutionSpaceDisplacement',componentSolutionSpaceDisplacement);
    save_print_figure(gcf,[saveFolder,'InitialOptimizedTrussComponent']);
    
    % initial truss + optimized truss + box solution space + component solution spaces
    plot_results_truss_generic_2d_moving_node(systemParameter,...
        'NodePositionInitial',initialDesign,...
        'NodePositionOptimalDisplacement',nodePositionOptimalDisplacement,...
        'BoxSolutionSpaceDisplacement',solutionSpaceBoxDisplacement,...
        'ComponentSolutionSpaceDisplacement',componentSolutionSpaceDisplacement);
    save_print_figure(gcf,[saveFolder,'InitialOptimizedTrussBoxComponent']);
    
    % initial truss + box solution space + component solution spaces
    plot_results_truss_generic_2d_moving_node(systemParameter,...
        'NodePositionInitial',initialDesign,...
        'BoxSolutionSpaceDisplacement',solutionSpaceBoxDisplacement,...
        'ComponentSolutionSpaceDisplacement',componentSolutionSpaceDisplacement);
    save_print_figure(gcf,[saveFolder,'InitialTrussBoxComponent']);

    % initial truss + sample trusses + component solution space
    nRandomTruss = 2;
    randomTrussMovingNode = candidate_space_sampling_individual_feasible(componentSolutionSpaceDisplacement,componentIndex,nRandomTruss);
    
    plot_results_truss_generic_2d_moving_node(systemParameter,...
        'NodePositionInitial',initialDesign,...
        'NodePositionRandomDisplacement',randomTrussMovingNode,...
        'ComponentSolutionSpaceDisplacement',componentSolutionSpaceDisplacement);
    save_print_figure(gcf,[saveFolder,'InitialRandomTrussComponent']);
end

if(computeDisplacement & computeMass)
    % initial truss + component solution spaces + mass
    plot_results_truss_generic_2d_moving_node(systemParameter,...
        'NodePositionInitial',initialDesign,...
        'ComponentSolutionSpaceDisplacement',componentSolutionSpaceDisplacement,...
        'ComponentSolutionSpaceMass',componentSolutionSpaceMass);
    save_print_figure(gcf,[saveFolder,'InitialTrussComponentWithMass']);
    
    % initial truss + optimized truss + component solution spaces + mass
    plot_results_truss_generic_2d_moving_node(systemParameter,...
        'NodePositionInitial',initialDesign,...
        'NodePositionOptimalDisplacement',nodePositionOptimalDisplacement,...
        'NodePositionOptimalMass',nodePositionOptimalMass,...
        'ComponentSolutionSpaceDisplacement',componentSolutionSpaceDisplacement,...
        'ComponentSolutionSpaceMass',componentSolutionSpaceMass);
    save_print_figure(gcf,[saveFolder,'InitialOptimizedTrussComponentWithMass']);
    
    % initial truss + optimized truss + box solution space + component solution spaces + mass
    plot_results_truss_generic_2d_moving_node(systemParameter,...
        'NodePositionInitial',initialDesign,...
        'NodePositionOptimalDisplacement',nodePositionOptimalDisplacement,...
        'NodePositionOptimalMass',nodePositionOptimalMass,...
        'BoxSolutionSpaceDisplacement',solutionSpaceBoxDisplacement,...
        'BoxSolutionSpaceMass',solutionSpaceBoxMass,...
        'ComponentSolutionSpaceDisplacement',componentSolutionSpaceDisplacement,...
        'ComponentSolutionSpaceMass',componentSolutionSpaceMass);
    save_print_figure(gcf,[saveFolder,'InitialOptimizedTrussBoxComponentWithMass']);
    
    % initial truss + box solution space + component solution spaces + mass
    plot_results_truss_generic_2d_moving_node(systemParameter,...
        'NodePositionInitial',initialDesign,...
        'BoxSolutionSpaceDisplacement',solutionSpaceBoxDisplacement,...
        'ComponentSolutionSpaceDisplacement',componentSolutionSpaceDisplacement,...
        'BoxSolutionSpaceMass',solutionSpaceBoxMass,...
        'ComponentSolutionSpaceMass',componentSolutionSpaceMass);
    save_print_figure(gcf,[saveFolder,'InitialTrussBoxComponentWithMass']);
end

if(computeDisplacementAndMass)
    % initial truss + component solution spaces
    plot_results_truss_generic_2d_moving_node(systemParameter,...
        'NodePositionInitial',initialDesign,...
        'ComponentSolutionSpaceDisplacementAndMass',componentSolutionSpaceDisplacementAndMass);
    save_print_figure(gcf,[saveFolder,'InitialTrussComponentAndMass']);
    
    % initial truss + box solution space + component solution spaces
    plot_results_truss_generic_2d_moving_node(systemParameter,...
        'NodePositionInitial',initialDesign,...
        'BoxSolutionSpaceDisplacement',solutionSpaceBoxDisplacementAndMass,...
        'ComponentSolutionSpaceDisplacement',componentSolutionSpaceDisplacementAndMass);
    save_print_figure(gcf,[saveFolder,'InitialTrussBoxComponentAndMass']);

    % initial truss + sample trusses + component solution space
    nRandomTruss = 2;
    randomTrussMovingNode = candidate_space_sampling_individual_feasible(componentSolutionSpaceDisplacementAndMass,componentIndex,nRandomTruss);
    
    plot_results_truss_generic_2d_moving_node(systemParameter,...
        'NodePositionInitial',initialDesign,...
        'NodePositionRandomDisplacementAndMass',randomTrussMovingNode,...
        'ComponentSolutionSpaceDisplacementAndMass',componentSolutionSpaceDisplacementAndMass);
    save_print_figure(gcf,[saveFolder,'InitialRandomTrussComponentAndMass']);
end


%% plot candidate spaces samples



%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;


%% subfunction: plot truss
% plot additional random truss
% nRandomTruss = 1;
% randomTrussMovingNode = candidate_space_sampling_individual_feasible(componentSolutionSpaceDisplacement,Components,nRandomTruss);
% 
% plot_results_truss_generic_2d_moving_node(systemParameter,...
%     'NodePositionInitial',initialDesign,...
%     'NodePositionRandomDisplacement',randomTrussMovingNode,...
%     'ComponentSolutionSpaceDisplacement',componentSolutionSpaceDisplacement);

function figureHandle = plot_results_truss_generic_2d_moving_node(systemParameter,varargin)
    parser = inputParser;
    parser.addParameter('NodePositionInitial',[]);
    parser.addParameter('NodePositionOptimalDisplacement',[]);
    parser.addParameter('NodePositionOptimalMass',[]);
    parser.addParameter('NodePositionRandomDisplacement',[]);
    parser.addParameter('NodePositionRandomMass',[]);
    parser.addParameter('NodePositionRandomDisplacementAndMass',[]);
    parser.addParameter('DeformationInitial',[]);
    parser.addParameter('DeformationOptimalDisplacement',[]);
    parser.addParameter('DeformationOptimalMass',[]);
    parser.addParameter('BoxSolutionSpaceDisplacement',[]);
    parser.addParameter('ComponentSolutionSpaceDisplacement',[]);
    parser.addParameter('BoxSolutionSpaceMass',[]);
    parser.addParameter('ComponentSolutionSpaceMass',[]);
    parser.addParameter('BoxSolutionSpaceDisplacementAndMass',[]);
    parser.addParameter('ComponentSolutionSpaceDisplacementAndMass',[]);
    parser.addParameter('IncludeWall',true);
    parser.addParameter('IncludeAppliedForce',false);
    parser.addParameter('IncludeAxesInformation',false);
    parser.addParameter('IncludeLegend',false);
    parser.addParameter('WallOptions',{});
    parser.addParameter('AppliedForceOptions',{});
    parser.addParameter('InitialTrussOptions',{});
    parser.addParameter('OptimalTrussDisplacementOptions',{});
    parser.addParameter('OptimalTrussMassOptions',{});
    parser.addParameter('RandomTrussDisplacementOptions',{});
    parser.addParameter('RandomTrussMassOptions',{});
    parser.addParameter('RandomTrussDisplacementAndMassOptions',{});
    parser.addParameter('BoxSolutionSpaceDisplacementOptions',{});
    parser.addParameter('ComponentSolutionSpaceDisplacementOptions',{});
    parser.addParameter('BoxSolutionSpaceMassOptions',{});
    parser.addParameter('ComponentSolutionSpaceMassOptions',{});
    parser.addParameter('BoxSolutionSpaceDisplacementAndMassOptions',{});
    parser.addParameter('ComponentSolutionSpaceDisplacementAndMassOptions',{});
    parser.addParameter('LegendOptions',{});
    parser.parse(varargin{:});
    options = parser.Results;

    wallX = [0 0];
    wallY = [-0.1 1.1];

    defaultInitialTrussOptions = {'ColorUndeformed',[0.8 0.8 0.8],'ColorDeformed','c','MaximumLinewidth',3.0,'DisplacementScaleFactor',100};
    [~,initialTrussOptions] = merge_name_value_pair_argument(defaultInitialTrussOptions,options.InitialTrussOptions);

    defaultOptimalTrussDisplacementOptions = {'ColorUndeformed','b','ColorDeformed','m','MaximumLinewidth',3.0,'DisplacementScaleFactor',100};
    [~,optimalTrussDisplacementOptions] = merge_name_value_pair_argument(defaultOptimalTrussDisplacementOptions,options.OptimalTrussDisplacementOptions);

    defaultOptimalTrussMassOptions = {'ColorUndeformed','y','ColorDeformed','k','MaximumLinewidth',3.0,'DisplacementScaleFactor',100};
    [~,optimalTrussMassOptions] = merge_name_value_pair_argument(defaultOptimalTrussMassOptions,options.OptimalTrussMassOptions);

    defaultRandomTrussDisplacementOptions = {'MaximumLinewidth',3.0};
    [~,randomTrussDisplacementOptions] = merge_name_value_pair_argument(defaultRandomTrussDisplacementOptions,options.RandomTrussDisplacementOptions);

    defaultRandomTrussMassOptions = {'MaximumLinewidth',3.0};
    [~,randomTrussMassOptions] = merge_name_value_pair_argument(defaultRandomTrussMassOptions,options.RandomTrussMassOptions);

    defaultRandomTrussDisplacementAndMassOptions = {'MaximumLinewidth',3.0};
    [~,randomTrussDisplacementAndMassOptions] = merge_name_value_pair_argument(defaultRandomTrussDisplacementAndMassOptions,options.RandomTrussDisplacementAndMassOptions);

    defaultWallOptions = {'Linewidth',8.0,'Color',[0.7 0.7 0.7]};
    [~,wallOptions] = merge_name_value_pair_argument(defaultWallOptions,options.WallOptions);

    defaultAppliedForceOptions = {'Color','r','LineWidth',3.0};
    [~,appliedForceOptions] = merge_name_value_pair_argument(defaultAppliedForceOptions,options.AppliedForceOptions);

    defaultBoxSolutionDisplacementOptions = {'EdgeColor','c','Linewidth',2.0};
    [~,boxSolutionDisplacementOptions] = merge_name_value_pair_argument(defaultBoxSolutionDisplacementOptions,options.BoxSolutionSpaceDisplacementOptions);

    defaultBoxSolutionMassOptions = {'EdgeColor',[0.9290 0.6940 0.1250],'Linewidth',2.0};
    [~,boxSolutionMassOptions] = merge_name_value_pair_argument(defaultBoxSolutionMassOptions,options.BoxSolutionSpaceMassOptions);

    defaultBoxSolutionDisplacementAndMassOptions = {'EdgeColor','c','Linewidth',2.0};
    [~,boxSolutionDisplacementAndMassOptions] = merge_name_value_pair_argument(defaultBoxSolutionDisplacementAndMassOptions,options.BoxSolutionSpaceDisplacementAndMassOptions);

    defaultComponentSolutionDisplacementOptions = {'EdgeColor','g','FaceColor','none','FaceAlpha',0.5,'Linewidth',2.0};
    [~,componentSolutionDisplacementOptions] = merge_name_value_pair_argument(defaultComponentSolutionDisplacementOptions,options.ComponentSolutionSpaceDisplacementOptions);

    defaultComponentSolutionMassOptions = {'EdgeColor','k','FaceColor','none','FaceAlpha',0.5,'Linewidth',2.0};
    [~,componentSolutionMassOptions] = merge_name_value_pair_argument(defaultComponentSolutionMassOptions,options.ComponentSolutionSpaceMassOptions);

    defaultComponentSolutionDisplacementAndMassOptions = {'EdgeColor','g','FaceColor','none','FaceAlpha',0.5,'Linewidth',2.0};
    [~,componentSolutionDisplacementAndMassOptions] = merge_name_value_pair_argument(defaultComponentSolutionDisplacementAndMassOptions,options.ComponentSolutionSpaceDisplacementAndMassOptions);

    defaultLegendOptions = {'location','west'};
    [~,legendOptions] = merge_name_value_pair_argument(defaultLegendOptions,options.LegendOptions);

    isDesignVariable = isnan(systemParameter.BaseNodePosition);
    nDimension = size(systemParameter.BaseNodePosition,2);
    nComponent = sum(isDesignVariable(:,1),1);
    isTrussTip = (systemParameter.NodeForce~=0);

    figureHandle = figure;
    hold all;
    
    % initial truss
    handleInitial = [];
    handleInitialDeformed = [];
    if(~isempty(options.NodePositionInitial))
        nodePositionInitial = systemParameter.BaseNodePosition;
        nodePositionInitial(isDesignVariable) = column_vector_to_row_major_matrix(options.NodePositionInitial',nDimension);

        [handleInitial,handleInitialDeformed] = plot_truss_deformation(gcf,nodePositionInitial,systemParameter.NodeElement,options.DeformationInitial,initialTrussOptions{:});
    end

    % optimal truss  - displacement 
    handleOptimalDisplacement = [];
    handleOptimizedDisplacementDeformed = [];
    if(~isempty(options.NodePositionOptimalDisplacement))
        nodePositionOptimized = systemParameter.BaseNodePosition;
        nodePositionOptimized(isDesignVariable) = column_vector_to_row_major_matrix(options.NodePositionOptimalDisplacement',nDimension);

        [handleOptimalDisplacement,handleOptimizedDisplacementDeformed] = plot_truss_deformation(gcf,nodePositionOptimized,systemParameter.NodeElement,options.DeformationOptimalDisplacement,optimalTrussDisplacementOptions{:});
    end

    % optimal truss - mass 
    handleOptimalMass = [];
    handleOptimizedMassDeformed = [];
    if(~isempty(options.NodePositionOptimalMass))
        nodePositionOptimized = systemParameter.BaseNodePosition;
        nodePositionOptimized(isDesignVariable) = column_vector_to_row_major_matrix(options.NodePositionOptimalMass',nDimension);

        [handleOptimalMass,handleOptimizedMassDeformed] = plot_truss_deformation(gcf,nodePositionOptimized,systemParameter.NodeElement,options.DeformationOptimalMass,optimalTrussMassOptions{:});
    end

    % box-shaped solution space  - displacement 
    handleToleranceNodeDisplacementBox = [];
    if(~isempty(options.BoxSolutionSpaceDisplacement))
        for i=1:nComponent-1
            currentIndex = 1 + nDimension*(i-1) + [0:(nDimension-1)];
            plot_design_box_2d(gcf,options.BoxSolutionSpaceDisplacement(:,currentIndex),boxSolutionDisplacementOptions{:});
        end
        currentIndex = 1 + nDimension*(nComponent-1) + [0:(nDimension-1)];
        handleToleranceNodeDisplacementBox = plot_design_box_2d(gcf,options.BoxSolutionSpaceDisplacement(:,currentIndex),boxSolutionDisplacementOptions{:});
    end

    % box-shaped solution space - mass
    handleToleranceNodeMassBox = [];
    if(~isempty(options.BoxSolutionSpaceMass))
        for i=1:nComponent-1
            currentIndex = 1 + nDimension*(i-1) + [0:(nDimension-1)];
            plot_design_box_2d(gcf,options.BoxSolutionSpaceMass(:,currentIndex),boxSolutionMassOptions{:});
        end
        currentIndex = 1 + nDimension*(nComponent-1) + [0:(nDimension-1)];
        handleToleranceNodeMassBox = plot_design_box_2d(gcf,options.BoxSolutionSpaceMass(:,currentIndex),boxSolutionMassOptions{:});
    end

    % box-shaped solution space  - displacement and mass
    handleToleranceNodeDisplacementAndMassBox = [];
    if(~isempty(options.BoxSolutionSpaceDisplacementAndMass))
        for i=1:nComponent-1
            currentIndex = 1 + nDimension*(i-1) + [0:(nDimension-1)];
            plot_design_box_2d(gcf,options.BoxSolutionSpaceDisplacementAndMass(:,currentIndex),boxSolutionDisplacementAndMassOptions{:});
        end
        currentIndex = 1 + nDimension*(nComponent-1) + [0:(nDimension-1)];
        handleToleranceNodeDisplacementAndMassBox = plot_design_box_2d(gcf,options.BoxSolutionSpaceDisplacementAndMass(:,currentIndex),boxSolutionDisplacementAndMassOptions{:});
    end

    % component solution space - displacement
    handleToleranceNodeDisplacementComponent = [];
    if(~isempty(options.ComponentSolutionSpaceDisplacement))
        for i=1:nComponent-1
            options.ComponentSolutionSpaceDisplacement(i).plot_candidate_space(gcf,componentSolutionDisplacementOptions{:});
        end
        handleToleranceNodeDisplacementComponent = options.ComponentSolutionSpaceDisplacement(nComponent).plot_candidate_space(gcf,componentSolutionDisplacementOptions{:});
    end

    % component solution space - mass
    handleToleranceNodeMassComponent = [];
    if(~isempty(options.ComponentSolutionSpaceMass))
        for i=1:nComponent-1
            options.ComponentSolutionSpaceMass(i).plot_candidate_space(gcf,componentSolutionMassOptions{:});
        end
        handleToleranceNodeMassComponent = options.ComponentSolutionSpaceMass(nComponent).plot_candidate_space(gcf,componentSolutionMassOptions{:});
    end

    % component solution space - displacement and mass
    handleToleranceNodeDisplacementAndMassComponent = [];
    if(~isempty(options.ComponentSolutionSpaceDisplacementAndMass))
        for i=1:nComponent-1
            options.ComponentSolutionSpaceDisplacementAndMass(i).plot_candidate_space(gcf,componentSolutionDisplacementAndMassOptions{:});
        end
        handleToleranceNodeDisplacementAndMassComponent = options.ComponentSolutionSpaceDisplacementAndMass(nComponent).plot_candidate_space(gcf,componentSolutionDisplacementAndMassOptions{:});
    end

    % random trusses - displacement
    handleRandomDisplacement = [];
    if(~isempty(options.NodePositionRandomDisplacement))
        nRandomTruss = size(options.NodePositionRandomDisplacement,1);
        trussColor = rand(nRandomTruss,3);
        for i=1:nRandomTruss
            randomNodePosition = systemParameter.BaseNodePosition;
            randomNodePosition(isDesignVariable) = column_vector_to_row_major_matrix(options.NodePositionRandomDisplacement(i,:)',nDimension);

            handleRandomDisplacement(i) = plot_truss_deformation(gcf,randomNodePosition,systemParameter.NodeElement,'ColorUndeformed',trussColor(i,:),randomTrussDisplacementOptions{:});
        end
    end

    % random trusses - mass
    handleRandomMass = [];
    if(~isempty(options.NodePositionRandomMass))
        nRandomTruss = size(options.NodePositionRandomMass,1);
        trussColor = rand(nRandomTruss,3);
        for i=1:nRandomTruss
            randomNodePosition = systemParameter.BaseNodePosition;
            randomNodePosition(isDesignVariable) = column_vector_to_row_major_matrix(options.NodePositionRandomMass(i,:)',nDimension);

            handleRandomMass(i) = plot_truss_deformation(gcf,randomNodePosition,systemParameter.NodeElement,'ColorUndeformed',trussColor(i,:),randomTrussMassOptions{:});
        end
    end

    % random trusses - displacement and mass
    handleRandomDisplacementAndMass = [];
    if(~isempty(options.NodePositionRandomDisplacementAndMass))
        nRandomTruss = size(options.NodePositionRandomDisplacementAndMass,1);
        trussColor = rand(nRandomTruss,3);
        for i=1:nRandomTruss
            randomNodePosition = systemParameter.BaseNodePosition;
            randomNodePosition(isDesignVariable) = column_vector_to_row_major_matrix(options.NodePositionRandomDisplacementAndMass(i,:)',nDimension);

            handleRandomDisplacementAndMass(i) = plot_truss_deformation(gcf,randomNodePosition,systemParameter.NodeElement,'ColorUndeformed',trussColor(i,:),randomTrussDisplacementAndMassOptions{:});
        end
    end

    % wall
    handleWall = [];
    if(options.IncludeWall)
        handleWall = plot(wallX,wallY,wallOptions{:});
    end

    % applied force
    handleForce = [];
    if(options.IncludeAppliedForce)
        positionTip = systemParameter.BaseNodePosition(isTrussTip);
        handleForce = quiver(positionTip(1),positionTip(2),0,-0.5,appliedForceOptions{:});
    end

    % axis
    grid minor;
    if(~options.IncludeAxesInformation)
        set(gca,'XColor', 'none','YColor','none');
    end

    % legend
    if(options.IncludeLegend)
        handleObjectAll = {...
            handleInitial,...
            handleInitialDeformed,...
            handleOptimalDisplacement,...
            handleOptimizedDisplacementDeformed,...
            handleOptimalMass,...
            handleOptimizedMassDeformed,...
            handleWall,...
            handleForce,...
            handleToleranceNodeDisplacementBox,...
            handleToleranceNodeDisplacementComponent,...
            handleToleranceNodeMassBox,...
            handleToleranceNodeMassComponent,...
            handleToleranceNodeDisplacementAndMassBox,...
            handleToleranceNodeDisplacementAndMassComponent};

        legendTextAll = {...
            'Initial Truss',...
            'Initial Truss (Deformed)',...
            'Optimized Truss - Displacement',...
            'Optimized Truss (Deformed) - Displacement',...
            'Optimized Truss - Mass',...
            'Optimized Truss (Deformed) - Mass',...
            'Wall',...
            'Applied Force',...
            'Tolerance Region for the Node (Box) - Displacement',...
            'Tolerance Region for the Node (Component) - Displacement',...
            'Tolerance Region for the Node (Box) - Mass',...
            'Tolerance Region for the Node (Component) - Mass',...
            'Tolerance Region for the Node (Box) - Displacement + Mass',...
            'Tolerance Region for the Node (Component) - Displacement + Mass'};
        
        handleObject = [handleObjectAll{:}];
        legentText = {legendTextAll{~cellfun(@isempty,handleObjectAll)}};

        legend(handleObject,legentText,legendOptions{:});
    end

    if(nargout<1)
        clear figureHandle;
    end
end

