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
trussAnalysisChoice = '16-DoF-2D';
computeDisplacement = false;
computeMass = false;
computeDisplacementAndMass = true;
computeDelaunayComponent = false;
useBoxResultForComponent = false;


%% function call
systemFunction = @truss_generic_moving_node;

switch trussAnalysisChoice
    case '2-DoF-2D'
        maxIter = 30;
        nSample = 100;
        growthRateDisplacement = 0.2;
        growthRateMass = 0.2;
        trimmingPasses = 'full';
        nRandomTruss = 5;
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
    case '4-DoF-2D'
        maxIter = 30;
        nSample = 100;
        growthRateDisplacement = 0.2;
        growthRateMass = 0.2;
        trimmingPasses = 'full';
        nRandomTruss = 1;
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
    case '16-DoF-2D'
        maxIter = 150;
        nSample = 100;
        growthRateDisplacement = 0.2;
        growthRateMass = 0.2;
        trimmingPasses = 'reduced';
        nRandomTruss = 1;
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
    case '36-DoF-2D'
        maxIter = 30;
        nSample = 100;
        growthRateDisplacement = 0.1;
        growthRateMass = 0.1;
        trimmingPasses = 'reduced';
        nRandomTruss = 1;
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
    case '9-DoF-3D'
        maxIter = 50;
        nSample = 100;
        growthRateDisplacement = 0.2;
        growthRateMass = 0.2;
        trimmingPasses = 'reduced';
        nRandomTruss = 1;
        systemParameter.ElementCrossSectionArea = 10; % [mm^2]
        systemParameter.ElementYoungsModulus = 210e3; % [MPa]
        systemParameter.ElementDensity = 7850e-9; % [kg/mm^3]
        systemParameter.BaseNodePosition = [...
              0   0   0; % (1)
              0 0.5   1; % (2)
              0   1   0; % (3)
            nan nan nan; % (4)
            nan nan nan; % (5) 
            nan nan nan; % (6)
              2 0.5 0.5]; % (7)
        systemParameter.FixedDegreesOfFreedom = [...
            ...
            true true true; % (1) 
            true true true; % (2)
            true true true; % (3)
            false false false; % (4)
            false false false; % (5)
            false false false; % (6)
            false false false]; % (7)
        systemParameter.NodeForce = [...
            0 0     0; % (1)
            0 0     0; % (2)
            0 0     0; % (3)
            0 0     0; % (4)
            0 0     0; % (5)
            0 0     0; % (6)
            0 0 -1000]; % (7)
        systemParameter.NodeElement = [...
            1 4; % (1)
            2 5; % (2)
            3 6; % (3)
            1 5; % (4)
            3 4; % (5)
            3 5; % CONFIRM WITH ZM
            4 5; % (6)
            5 6; % (7)
            4 6; % (8)
            4 7; % (9)
            5 7; % (10)
            6 7]; % (11)
        %                x1 y1 z1 x2  y2 z2 x3 y3 z3
        initialDesign = [ 1  0  0  1 0.5  1  1  1  0];
        %                                     x1   y1   z1  x2   y2  z2  x3  y3   z3
        designSpaceLowerBoundDisplacement = [0.1 -0.5 -1.0 0.1  0.0 0.5 0.1 0.0 -1.0];
        designSpaceUpperBoundDisplacement = [2.0  0.5  0.5 2.0  1.0 1.5 2.0 2.0  1.0];
        designSpaceLowerBoundMass = [0 -0.5 -0.5 0 -0.5 -0.5 0 -0.5 -0.5];
        designSpaceUpperBoundMass = [2  1.5  1.5 2  1.5  1.5 2  1.5  1.5];
    case '36-DoF-3D'
        maxIter = 250;
        nSample = 100;
        growthRateDisplacement = 0.1;
        growthRateMass = 0.1;
        trimmingPasses = 'single';
        nRandomTruss = 1;
        systemParameter.ElementCrossSectionArea = 10; % [mm^2]
        systemParameter.ElementYoungsModulus = 210e3; % [MPa]
        systemParameter.ElementDensity = 7850e-9; % [kg/mm^3]
        systemParameter.BaseNodePosition = [...
              0   0   0; % (1)
              0 0.5   1; % (2)
              0   1   0; % (3)
            nan nan nan; % (4)
            nan nan nan; % (5) 
            nan nan nan; % (6)
            nan nan nan; % (7)
            nan nan nan; % (8) 
            nan nan nan; % (9)
            nan nan nan; % (10)
            nan nan nan; % (11) 
            nan nan nan; % (12)
            nan nan nan; % (13)
            nan nan nan; % (14) 
            nan nan nan; % (15)
              5 0.5 0.5]; % (16)
        systemParameter.FixedDegreesOfFreedom = [...
            true true true; % (1) 
            true true true; % (2)
            true true true; % (3)
            false false false; % (4)
            false false false; % (5)
            false false false; % (6)
            false false false; % (7)
            false false false; % (8)
            false false false; % (9)
            false false false; % (10)
            false false false; % (11)
            false false false; % (12)
            false false false; % (13)
            false false false; % (14)
            false false false; % (15)
            false false false]; % (16)
        systemParameter.NodeForce = [...
            0 0     0; % (1)
            0 0     0; % (2)
            0 0     0; % (3)
            0 0     0; % (4)
            0 0     0; % (5)
            0 0     0; % (6)
            0 0     0; % (7)
            0 0     0; % (8)
            0 0     0; % (9)
            0 0     0; % (10)
            0 0     0; % (11)
            0 0     0; % (12)
            0 0     0; % (13)
            0 0     0; % (14)
            0 0     0; % (15)
            0 0 -1000]; % (16)
        systemParameter.NodeElement = [...
            ... % 
            1 4; % (1) -
            2 5; % (2) -
            3 6; % (3) -
            1 5; % (4) /
            3 4; % (5) \
            3 5; % (6) /
            4 5; % (7) |
            5 6; % (8) \
            4 6; % (9) /
            ... %
            4 7;  % (10) -
            5 8;  % (11) -
            6 9;  % (12) -
            4 8;  % (13) /
            6 7;  % (14) \
            6 8;  % (15) /
            7 8;  % (16) |
            8 9;  % (17) \
            7 9;  % (18) /
            ... %
            7 10; % (19) -
            8 11; % (20) -
            9 12; % (21) -
            7 11; % (22) /
            9 10; % (23) \
            9 11; % (24) /
            10 11; % (25) |
            11 12; % (26) \
            10 12; % (27) /
            ... %
            10 13; % (28) -
            11 14; % (29) -
            12 15; % (30) -
            10 14; % (31) /
            12 13; % (32) \
            12 14; % (33) /
            13 14; % (34) |
            14 15; % (35) \
            13 15; % (36) /
            ... %
            13 16; % (37)
            14 16; % (38)
            15 16]; % (39)
        initialDesign = [...
            1    0  0 ... % (1)
            1  0.5  1 ... % (2)
            1    1  0 ... % (3)
            2    0  0 ... % (4)
            2  0.5  1 ... % (5)
            2    1  0 ... % (6)
            3    0  0 ... % (7)
            3  0.5  1 ... % (8)
            3    1  0 ... % (9)
            4    0  0 ... % (10)
            4  0.5  1 ... % (11)
            4    1  0]; % (12)
        designSpaceLowerBoundDisplacement = initialDesign - repmat([0.9 0.5 0.5],1,12);
        designSpaceUpperBoundDisplacement = initialDesign + repmat([0.9 0.5 0.5],1,12);
        designSpaceLowerBoundMass = repmat([0 -0.5 -0.5],1,12);
        designSpaceUpperBoundMass = repmat([5  1.5  1.5],1,12);
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
is3dPlot = (nDimension==3);

% performance limits
nElement = size(systemParameter.NodeElement,1);
performanceLowerLimit = -inf;
performanceUpperLimit = [nan nan repmat(inf,1,nElement)];

% initial truss
plot_results_truss_generic_moving_node(is3dPlot,systemParameter,...
    'NodePositionInitial',initialDesign);
save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialTruss']);
save_print_figure(gcf,[saveFolder,'InitialTruss'],'PrintFormat',{'png','pdf'});

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

plot_results_truss_generic_moving_node(is3dPlot,systemParameter,...
    'NodePositionInitial',initialDesign,...
    'DeformationInitial',nodeDisplacementInitial);
save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialTrussDeformation']);
save_print_figure(gcf,[saveFolder,'InitialTrussDeformation'],'PrintFormat',{'png','pdf'});


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
    plot_results_truss_generic_moving_node(is3dPlot,systemParameter,...
        'NodePositionInitial',initialDesign,...
        'NodePositionOptimalDisplacement',nodePositionOptimalDisplacement);
    save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialOptimizedTruss']);
    save_print_figure(gcf,[saveFolder,'InitialOptimizedTruss'],'PrintFormat',{'png','pdf'});
    
    if(computeMass)
        % initial truss + optimized displacment truss + optimized mass truss
        plot_results_truss_generic_moving_node(is3dPlot,systemParameter,...
            'NodePositionInitial',initialDesign,...
            'NodePositionOptimalDisplacement',nodePositionOptimalDisplacement,...
            'NodePositionOptimalMass',nodePositionOptimalMass);
        save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialOptimizedTrussWithMass']);
        save_print_figure(gcf,[saveFolder,'InitialOptimizedTrussWithMass'],'PrintFormat',{'png','pdf'});
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

    plot_results_truss_generic_moving_node(is3dPlot,systemParameter,...
        'NodePositionInitial',initialDesign,...
        'NodePositionOptimalDisplacement',nodePositionOptimalDisplacement,...
        'DeformationInitial',nodeDisplacementInitial,...
        'DeformationOptimalDisplacement',nodeDisplacementOptimal);
    save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialOptimizedTrussDeformation']);
    save_print_figure(gcf,[saveFolder,'InitialOptimizedTrussDeformation'],'PrintFormat',{'png','pdf'});
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
        'NumberSamplesPerIterationExploration',nSample,...
        'NumberSamplesPerIterationConsolidation',nSample,...
        'FixIterNumberExploration',true,...
        'FixIterNumberConsolidation',true,...
        'MaxIterExploration',maxIter,...
        'MaxIterConsolidation',maxIter,...
        'UseAdaptiveGrowthRate',true,...
        'GrowthRate',growthRateDisplacement,...
        'ApplyLeanness','never',...
        'TrimmingOperationOptions',{'PassesCriterion',trimmingPasses},...
        'TrimmingOrderOptions',{'OrderPreference','score'});
    
    rng(rngState);
    [solutionSpaceBoxDisplacement,problemDataBoxDisplacement,iterDataBoxDisplacement] = sso_box_stochastic(designEvaluatorDisplacement,...
        initialDesign,designSpaceLowerBoundDisplacement,designSpaceUpperBoundDisplacement,optionsBox);
    toc(timeElapsedBox)

    resultsFolder = [saveFolder,'BoxDisplacement/'];
    mkdir(resultsFolder);
    algoDataBoxDisplacement = postprocess_sso_box_stochastic(problemDataBoxDisplacement,iterDataBoxDisplacement);
    plot_sso_box_stochastic_metrics(algoDataBoxDisplacement,...
        'SaveFolder',resultsFolder,...
        'CloseFigureAfterSaving',true,...
        'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});
end

if(computeMass)
    timeElapsedBox = tic;
    optionsBox = sso_stochastic_options('box',...
        'NumberSamplesPerIterationExploration',nSample,...
        'NumberSamplesPerIterationConsolidation',nSample,...
        'FixIterNumberExploration',true,...
        'FixIterNumberConsolidation',true,...
        'MaxIterExploration',maxIter,...
        'MaxIterConsolidation',maxIter,...
        'UseAdaptiveGrowthRate',true,...
        'GrowthRate',growthRateMass,...
        'ApplyLeanness','never',...
        'TrimmingOperationOptions',{'PassesCriterion',trimmingPasses},...
        'TrimmingOrderOptions',{'OrderPreference','score'});
    
    rng(rngState);
    [solutionSpaceBoxMass,problemDataBoxMass,iterDataBoxMass] = sso_box_stochastic(designEvaluatorMass,...
        initialDesign,designSpaceLowerBoundMass,designSpaceUpperBoundMass,optionsBox);
    toc(timeElapsedBox)

    resultsFolder = [saveFolder,'BoxMass/'];
    mkdir(resultsFolder);
    algoDataBoxMass = postprocess_sso_box_stochastic(problemDataBoxMass,iterDataBoxMass);
    plot_sso_box_stochastic_metrics(algoDataBoxMass,...
        'SaveFolder',resultsFolder,...
        'CloseFigureAfterSaving',true,...
        'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});
end

if(computeDisplacementAndMass)
    timeElapsedBox = tic;
    optionsBox = sso_stochastic_options('box',...
        'NumberSamplesPerIterationExploration',nSample,...
        'NumberSamplesPerIterationConsolidation',nSample,...
        'FixIterNumberExploration',true,...
        'FixIterNumberConsolidation',true,...
        'MaxIterExploration',maxIter,...
        'MaxIterConsolidation',maxIter,...
        'UseAdaptiveGrowthRate',true,...
        'GrowthRate',growthRateDisplacement,...
        'ApplyLeanness','never',...
        'TrimmingOperationOptions',{'PassesCriterion',trimmingPasses},...
        'TrimmingOrderOptions',{'OrderPreference','score'});
    
    rng(rngState);
    [solutionSpaceBoxDisplacementAndMass,problemDataBoxDisplacementAndMass,iterDataBoxDisplacementAndMass] = ...
        sso_box_stochastic(designEvaluatorDisplacementAndMass,...
        initialDesign,designSpaceLowerBoundDisplacement,designSpaceUpperBoundDisplacement,optionsBox);
    toc(timeElapsedBox)

    resultsFolder = [saveFolder,'BoxDisplacementAndMass/'];
    mkdir(resultsFolder);
    algoDataBoxDisplacementAndMass = postprocess_sso_box_stochastic(problemDataBoxDisplacementAndMass,iterDataBoxDisplacementAndMass);
    plot_sso_box_stochastic_metrics(algoDataBoxDisplacementAndMass,...
        'SaveFolder',resultsFolder,...
        'CloseFigureAfterSaving',true,...
        'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});
end


%% Component Opt
if(computeDisplacement)
    % convex 
    timeElapsedComponent = tic;
    optionsComponent = sso_stochastic_options('component',...
        'NumberSamplesPerIterationExploration',nSample,...
        'NumberSamplesPerIterationConsolidation',nSample,...
        'FixIterNumberExploration',true,...
        'FixIterNumberConsolidation',true,...
        'MaxIterExploration',maxIter,...
        'MaxIterConsolidation',maxIter,...
        'CandidateSpaceConstructorExploration',@CandidateSpaceConvexHull,...
        'CandidateSpaceConstructorConsolidation',@CandidateSpaceConvexHull,...
        'TrimmingMethodFunction',@component_trimming_method_planar_trimming,...
        ... 'TrimmingMethodOptions',{'ReferenceDesigns','boundary-center'},...
        'UseAdaptiveGrowthRate',true,...
        'GrowthRate',growthRateDisplacement,...
        'ApplyLeanness','never',...
        'UsePaddingSamplesInTrimming',true,...
        'UsePreviousEvaluatedSamplesConsolidation',false,...
        'UsePreviousPaddingSamplesConsolidation',false,...
        'TrimmingOperationOptions',{'PassesCriterion',trimmingPasses},...
        'TrimmingOrderOptions',{'OrderPreference','score'});


    if(useBoxResultForComponent)
        initialDesignComponent = solutionSpaceBoxDisplacement;
    else
        initialDesignComponent = initialDesign;
    end
    
    rng(rngState);
    [componentSolutionSpaceConvexDisplacement,problemDataComponentConvexDisplacement,iterDataComponentConvexDisplacement] = sso_component_stochastic(designEvaluatorDisplacement,...
        initialDesignComponent,designSpaceLowerBoundDisplacement,designSpaceUpperBoundDisplacement,componentIndex,optionsComponent);
    toc(timeElapsedComponent)

    resultsFolder = [saveFolder,'ComponentConvexDisplacement/'];
    mkdir(resultsFolder);
    algoDataComponentConvexDisplacement = postprocess_sso_component_stochastic(problemDataComponentConvexDisplacement,iterDataComponentConvexDisplacement);
    plot_sso_component_stochastic_metrics(algoDataComponentConvexDisplacement,...
        'SaveFolder',resultsFolder,...
        'CloseFigureAfterSaving',true,...
        'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});

    % non-convex
    componentSolutionSpaceDelaunayDisplacement = [];
    if(computeDelaunayComponent)
        timeElapsedComponent = tic;
        optionsComponent = sso_stochastic_options('component',...
            'NumberSamplesPerIterationExploration',nSample,...
            'NumberSamplesPerIterationConsolidation',nSample,...
            'FixIterNumberExploration',true,...
            'FixIterNumberConsolidation',true,...
            'MaxIterExploration',maxIter,...
            'MaxIterConsolidation',maxIter,...
            'CandidateSpaceConstructorExploration',@CandidateSpaceDelaunay,...
            'CandidateSpaceConstructorConsolidation',@CandidateSpaceDelaunay,...
            'TrimmingMethodFunction',@component_trimming_method_corner_box_removal,...
            'UseAdaptiveGrowthRate',true,...
            'GrowthRate',growthRateDisplacement,...
            'ApplyLeanness','never',...
            'UsePaddingSamplesInTrimming',true,...
            'UsePreviousEvaluatedSamplesConsolidation',false,...
            'UsePreviousPaddingSamplesConsolidation',false,...
            'TrimmingOperationOptions',{'PassesCriterion',trimmingPasses},...
            'TrimmingOrderOptions',{'OrderPreference','score'});
        
        rng(rngState);
        [componentSolutionSpaceDelaunayDisplacement,problemDataComponentDelaunayDisplacement,iterDataComponentDelaunayDisplacement] = sso_component_stochastic(designEvaluatorDisplacement,...
            initialDesignComponent,designSpaceLowerBoundDisplacement,designSpaceUpperBoundDisplacement,componentIndex,optionsComponent);
        toc(timeElapsedComponent)

        resultsFolder = [saveFolder,'ComponentDelaunayDisplacement/'];
        mkdir(resultsFolder);
        algoDataComponentDelaunayDisplacement = postprocess_sso_component_stochastic(problemDataComponentDelaunayDisplacement,iterDataComponentDelaunayDisplacement);
        plot_sso_component_stochastic_metrics(algoDataComponentDelaunayDisplacement,...
            'SaveFolder',resultsFolder,...
            'CloseFigureAfterSaving',true,...
            'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});

        comparisonComponent = {algoDataComponentConvexDisplacement,algoDataComponentDelaunayDisplacement};
        componentLabel = {'Planar Trimming','Corner Box Removal'};
    else
        comparisonComponent = {algoDataComponentConvexDisplacement};
        componentLabel = {};
    end

    % comparison
    resultsFolder = [saveFolder,'ComparisonDisplacement/'];
    mkdir(resultsFolder);
    plot_sso_comparison_box_component_stochastic_metrics(...
        {algoDataBoxDisplacement},...
        comparisonComponent,...
        'ComponentLabel',componentLabel,...
        'SaveFolder',resultsFolder,...
        'CloseFigureAfterSaving',true,...
        'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});
end

if(computeMass)
    % convex
    timeElapsedComponent = tic;
    optionsComponent = sso_stochastic_options('component',...
        'NumberSamplesPerIterationExploration',nSample,...
        'NumberSamplesPerIterationConsolidation',nSample,...
        'FixIterNumberExploration',true,...
        'FixIterNumberConsolidation',true,...
        'MaxIterExploration',maxIter,...
        'MaxIterConsolidation',maxIter,...
        'CandidateSpaceConstructorExploration',@CandidateSpaceConvexHull,...
        'CandidateSpaceConstructorConsolidation',@CandidateSpaceConvexHull,...
        'TrimmingMethodFunction',@component_trimming_method_planar_trimming,...
        ... 'TrimmingMethodOptions',{'ReferenceDesigns','boundary-center'},...
        'UseAdaptiveGrowthRate',true,...
        'GrowthRate',growthRateMass,...
        'ApplyLeanness','never',...
        'UsePaddingSamplesInTrimming',true,...
        'UsePreviousEvaluatedSamplesConsolidation',false,...
        'UsePreviousPaddingSamplesConsolidation',false,...
        'TrimmingOperationOptions',{'PassesCriterion',trimmingPasses},...
        'TrimmingOrderOptions',{'OrderPreference','score'});

    if(useBoxResultForComponent)
        initialDesignComponent = solutionSpaceBoxMass;
    else
        initialDesignComponent = initialDesign;
    end
    
    rng(rngState);
    [componentSolutionSpaceConvexMass,problemDataComponentConvexMass,iterDataComponentConvexMass] = sso_component_stochastic(designEvaluatorMass,...
        initialDesignComponent,designSpaceLowerBoundMass,designSpaceUpperBoundMass,componentIndex,optionsComponent);
    toc(timeElapsedComponent)

    resultsFolder = [saveFolder,'ComponentConvexMass/'];
    mkdir(resultsFolder);
    algoDataComponentConvexMass = postprocess_sso_component_stochastic(problemDataComponentConvexMass,iterDataComponentConvexMass);
    plot_sso_component_stochastic_metrics(algoDataComponentConvexMass,...
        'SaveFolder',resultsFolder,...
        'CloseFigureAfterSaving',true,...
        'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});

    % delaunay
    componentSolutionSpaceDelaunayMass = [];
    if(computeDelaunayComponent)
        timeElapsedComponent = tic;
        optionsComponent = sso_stochastic_options('component',...
            'NumberSamplesPerIterationExploration',nSample,...
            'NumberSamplesPerIterationConsolidation',nSample,...
            'FixIterNumberExploration',true,...
            'FixIterNumberConsolidation',true,...
            'MaxIterExploration',maxIter,...
            'MaxIterConsolidation',maxIter,...
            'CandidateSpaceConstructorExploration',@CandidateSpaceDelaunay,...
            'CandidateSpaceConstructorConsolidation',@CandidateSpaceDelaunay,...
            'TrimmingMethodFunction',@component_trimming_method_corner_box_removal,...
            'UseAdaptiveGrowthRate',true,...
            'GrowthRate',growthRateMass,...
            'ApplyLeanness','never',...
            'UsePaddingSamplesInTrimming',true,...
            'UsePreviousEvaluatedSamplesConsolidation',false,...
            'UsePreviousPaddingSamplesConsolidation',false,...
            'TrimmingOperationOptions',{'PassesCriterion',trimmingPasses},...
            'TrimmingOrderOptions',{'OrderPreference','score'});
        
        rng(rngState);
        [componentSolutionSpaceDelaunayMass,problemDataComponentDelaunayMass,iterDataComponentDelaunayMass] = sso_component_stochastic(designEvaluatorMass,...
            initialDesignComponent,designSpaceLowerBoundMass,designSpaceUpperBoundMass,componentIndex,optionsComponent);
        toc(timeElapsedComponent)

        resultsFolder = [saveFolder,'ComponentDelaunayMass/'];
        mkdir(resultsFolder);
        algoDataComponentDelaunayMass = postprocess_sso_component_stochastic(problemDataComponentDelaunayMass,iterDataComponentDelaunayMass);
        plot_sso_component_stochastic_metrics(algoDataComponentDelaunayMass,...
            'SaveFolder',resultsFolder,...
            'CloseFigureAfterSaving',true,...
            'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});


        comparisonComponent = {algoDataComponentConvexMass,algoDataComponentDelaunayMass};
        componentLabel = {'Planar Trimming','Corner Box Removal'};
    else
        comparisonComponent = {algoDataComponentConvexMass};
        componentLabel = {};
    end

    % comparison
    resultsFolder = [saveFolder,'ComparisonMass/'];
    mkdir(resultsFolder);
    plot_sso_comparison_box_component_stochastic_metrics(...
        {algoDataBoxMass},...
        comparisonComponent,...
        'ComponentLabel',componentLabel,...
        'SaveFolder',resultsFolder,...
        'CloseFigureAfterSaving',true,...
        'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});
end

if(computeDisplacementAndMass)
    % convex
    timeElapsedComponent = tic;
    optionsComponent = sso_stochastic_options('component',...
        'NumberSamplesPerIterationExploration',nSample,...
        'NumberSamplesPerIterationConsolidation',nSample,...
        'FixIterNumberExploration',true,...
        'FixIterNumberConsolidation',true,...
        'MaxIterExploration',maxIter,...
        'MaxIterConsolidation',maxIter,...
        'CandidateSpaceConstructorExploration',@CandidateSpaceConvexHull,...
        'CandidateSpaceConstructorConsolidation',@CandidateSpaceConvexHull,...
        'TrimmingMethodFunction',@component_trimming_method_planar_trimming,...
        ... 'TrimmingMethodOptions',{'ReferenceDesigns','boundary-center'},...
        'UseAdaptiveGrowthRate',true,...
        'GrowthRate',growthRateDisplacement,...
        'ApplyLeanness','never',...
        'UsePaddingSamplesInTrimming',true,...
        'UsePreviousEvaluatedSamplesConsolidation',false,...
        'UsePreviousPaddingSamplesConsolidation',false,...
        'TrimmingOperationOptions',{'PassesCriterion',trimmingPasses},...
        'TrimmingOrderOptions',{'OrderPreference','score'});

    if(useBoxResultForComponent)
        initialDesignComponent = solutionSpaceBoxDisplacementAndMass;
    else
        initialDesignComponent = initialDesign;
    end
    
    rng(rngState);
    [componentSolutionSpaceConvexDisplacementAndMass,problemDataComponentConvexDisplacementAndMass,iterDataComponentConvexDisplacementAndMass] = ...
        sso_component_stochastic(designEvaluatorDisplacementAndMass,...
        initialDesignComponent,designSpaceLowerBoundDisplacement,designSpaceUpperBoundDisplacement,componentIndex,optionsComponent);
    toc(timeElapsedComponent)

    resultsFolder = [saveFolder,'ComponentConvexDisplacementAndMass/'];
    mkdir(resultsFolder);
    algoDataComponentConvexDisplacementAndMass = postprocess_sso_component_stochastic(problemDataComponentConvexDisplacementAndMass,iterDataComponentConvexDisplacementAndMass);
    plot_sso_component_stochastic_metrics(algoDataComponentConvexDisplacementAndMass,...
        'SaveFolder',resultsFolder,...
        'CloseFigureAfterSaving',true,...
        'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});

    % delaunay
    componentSolutionSpaceDelaunayDisplacementAndMass = [];
    if(computeDelaunayComponent)
        timeElapsedComponent = tic;
        optionsComponent = sso_stochastic_options('component',...
            'NumberSamplesPerIterationExploration',nSample,...
            'NumberSamplesPerIterationConsolidation',nSample,...
            'FixIterNumberExploration',true,...
            'FixIterNumberConsolidation',true,...
            'MaxIterExploration',maxIter,...
            'MaxIterConsolidation',maxIter,...
            'CandidateSpaceConstructorExploration',@CandidateSpaceDelaunay,...
            'CandidateSpaceConstructorConsolidation',@CandidateSpaceDelaunay,...
            'TrimmingMethodFunction',@component_trimming_method_corner_box_removal,...
            'UseAdaptiveGrowthRate',true,...
            'GrowthRate',growthRateDisplacement,...
            'ApplyLeanness','never',...
            'UsePaddingSamplesInTrimming',true,...
            'UsePreviousEvaluatedSamplesConsolidation',false,...
            'UsePreviousPaddingSamplesConsolidation',false,...
            'TrimmingOperationOptions',{'PassesCriterion',trimmingPasses},...
            'TrimmingOrderOptions',{'OrderPreference','score'});
        
        rng(rngState);
        [componentSolutionSpaceDelaunayDisplacementAndMass,problemDataComponentDelaunayDisplacementAndMass,iterDataComponentDelaunayDisplacementAndMass] = ...
            sso_component_stochastic(designEvaluatorDisplacementAndMass,...
            initialDesignComponent,designSpaceLowerBoundDisplacement,designSpaceUpperBoundDisplacement,componentIndex,optionsComponent);
        toc(timeElapsedComponent)

        resultsFolder = [saveFolder,'ComponentDelaunayDisplacementAndMass/'];
        mkdir(resultsFolder);
        algoDataComponentDelaunayDisplacementAndMass = postprocess_sso_component_stochastic(problemDataComponentDelaunayDisplacementAndMass,iterDataComponentDelaunayDisplacementAndMass);
        plot_sso_component_stochastic_metrics(algoDataComponentDelaunayDisplacementAndMass,...
            'SaveFolder',resultsFolder,...
            'CloseFigureAfterSaving',true,...
            'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});

        comparisonComponent = {algoDataComponentConvexDisplacementAndMass,algoDataComponentDelaunayDisplacementAndMass};
        componentLabel = {'Planar Trimming','Corner Box Removal'};
    else
        comparisonComponent = {algoDataComponentConvexDisplacementAndMass};
        componentLabel = {};
    end

    % comparison
    resultsFolder = [saveFolder,'ComparisonDisplacementAndMass/'];
    mkdir(resultsFolder);
    plot_sso_comparison_box_component_stochastic_metrics(...
        {algoDataBoxDisplacementAndMass},...
        comparisonComponent,...
        'ComponentLabel',componentLabel,...
        'SaveFolder',resultsFolder,...
        'CloseFigureAfterSaving',true,...
        'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});
end


%% measure comparison



%% Plot Visualization
if(computeDisplacement)
    % initial truss + component solution spaces
    plot_results_truss_generic_moving_node(is3dPlot,systemParameter,...
        'NodePositionInitial',initialDesign,...
        'ComponentSolutionSpaceConvexDisplacement',componentSolutionSpaceConvexDisplacement,...
        'ComponentSolutionSpaceDelaunayDisplacement',componentSolutionSpaceDelaunayDisplacement);
    save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialTrussComponent']);
    save_print_figure(gcf,[saveFolder,'InitialTrussComponent'],'PrintFormat',{'png','pdf'});
    
    % initial truss + optimized truss + component solution spaces
    plot_results_truss_generic_moving_node(is3dPlot,systemParameter,...
        'NodePositionInitial',initialDesign,...
        'NodePositionOptimalDisplacement',nodePositionOptimalDisplacement,...
        'ComponentSolutionSpaceConvexDisplacement',componentSolutionSpaceConvexDisplacement,...
        'ComponentSolutionSpaceDelaunayDisplacement',componentSolutionSpaceDelaunayDisplacement);
    save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialOptimizedTrussComponent']);
    save_print_figure(gcf,[saveFolder,'InitialOptimizedTrussComponent'],'PrintFormat',{'png','pdf'});
    
    % initial truss + optimized truss + box solution space + component solution spaces
    plot_results_truss_generic_moving_node(is3dPlot,systemParameter,...
        'NodePositionInitial',initialDesign,...
        'NodePositionOptimalDisplacement',nodePositionOptimalDisplacement,...
        'BoxSolutionSpaceDisplacement',solutionSpaceBoxDisplacement,...
        'ComponentSolutionSpaceConvexDisplacement',componentSolutionSpaceConvexDisplacement,...
        'ComponentSolutionSpaceDelaunayDisplacement',componentSolutionSpaceDelaunayDisplacement);
    save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialOptimizedTrussBoxComponent']);
    save_print_figure(gcf,[saveFolder,'InitialOptimizedTrussBoxComponent'],'PrintFormat',{'png','pdf'});
    
    % initial truss + box solution space + component solution spaces
    plot_results_truss_generic_moving_node(is3dPlot,systemParameter,...
        'NodePositionInitial',initialDesign,...
        'BoxSolutionSpaceDisplacement',solutionSpaceBoxDisplacement,...
        'ComponentSolutionSpaceConvexDisplacement',componentSolutionSpaceConvexDisplacement,...
        'ComponentSolutionSpaceDelaunayDisplacement',componentSolutionSpaceDelaunayDisplacement);
    save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialTrussBoxComponent']);
    save_print_figure(gcf,[saveFolder,'InitialTrussBoxComponent'],'PrintFormat',{'png','pdf'});

    % initial truss + sample trusses + component solution space
    randomTrussMovingNode = candidate_space_sampling_individual_feasible(componentSolutionSpaceConvexDisplacement,componentIndex,nRandomTruss);
    
    plot_results_truss_generic_moving_node(is3dPlot,systemParameter,...
        'NodePositionInitial',initialDesign,...
        'NodePositionRandomDisplacement',randomTrussMovingNode,...
        'ComponentSolutionSpaceConvexDisplacement',componentSolutionSpaceConvexDisplacement,...
        'ComponentSolutionSpaceDelaunayDisplacement',componentSolutionSpaceDelaunayDisplacement);
    save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialRandomTrussComponent']);
    save_print_figure(gcf,[saveFolder,'InitialRandomTrussComponent'],'PrintFormat',{'png','pdf'});
end

if(computeDisplacement & computeMass)
    % initial truss + component solution spaces + mass
    plot_results_truss_generic_moving_node(is3dPlot,systemParameter,...
        'NodePositionInitial',initialDesign,...
        'ComponentSolutionSpaceConvexDisplacement',componentSolutionSpaceConvexDisplacement,...
        'ComponentSolutionSpaceConvexMass',componentSolutionSpaceConvexMass,...
        'ComponentSolutionSpaceDelaunayDisplacement',componentSolutionSpaceDelaunayDisplacement,...
        'ComponentSolutionSpaceDelaunayMass',componentSolutionSpaceDelaunayMass);
    save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialTrussComponentWithMass']);
    save_print_figure(gcf,[saveFolder,'InitialTrussComponentWithMass'],'PrintFormat',{'png','pdf'});
    
    % initial truss + optimized truss + component solution spaces + mass
    plot_results_truss_generic_moving_node(is3dPlot,systemParameter,...
        'NodePositionInitial',initialDesign,...
        'NodePositionOptimalDisplacement',nodePositionOptimalDisplacement,...
        'NodePositionOptimalMass',nodePositionOptimalMass,...
        'ComponentSolutionSpaceConvexDisplacement',componentSolutionSpaceConvexDisplacement,...
        'ComponentSolutionSpaceConvexMass',componentSolutionSpaceConvexMass,...
        'ComponentSolutionSpaceDelaunayDisplacement',componentSolutionSpaceDelaunayDisplacement,...
        'ComponentSolutionSpaceDelaunayMass',componentSolutionSpaceDelaunayMass);
    save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialOptimizedTrussComponentWithMass']);
    save_print_figure(gcf,[saveFolder,'InitialOptimizedTrussComponentWithMass'],'PrintFormat',{'png','pdf'});
    
    % initial truss + optimized truss + box solution space + component solution spaces + mass
    plot_results_truss_generic_moving_node(is3dPlot,systemParameter,...
        'NodePositionInitial',initialDesign,...
        'NodePositionOptimalDisplacement',nodePositionOptimalDisplacement,...
        'NodePositionOptimalMass',nodePositionOptimalMass,...
        'BoxSolutionSpaceDisplacement',solutionSpaceBoxDisplacement,...
        'BoxSolutionSpaceMass',solutionSpaceBoxMass,...
        'ComponentSolutionSpaceConvexDisplacement',componentSolutionSpaceConvexDisplacement,...
        'ComponentSolutionSpaceConvexMass',componentSolutionSpaceConvexMass,...
        'ComponentSolutionSpaceDelaunayDisplacement',componentSolutionSpaceDelaunayDisplacement,...
        'ComponentSolutionSpaceDelaunayMass',componentSolutionSpaceDelaunayMass);
    save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialOptimizedTrussBoxComponentWithMass']);
    save_print_figure(gcf,[saveFolder,'InitialOptimizedTrussBoxComponentWithMass'],'PrintFormat',{'png','pdf'});
    
    % initial truss + box solution space + component solution spaces + mass
    plot_results_truss_generic_moving_node(is3dPlot,systemParameter,...
        'NodePositionInitial',initialDesign,...
        'BoxSolutionSpaceDisplacement',solutionSpaceBoxDisplacement,...
        'ComponentSolutionSpaceConvexDisplacement',componentSolutionSpaceConvexDisplacement,...
        'ComponentSolutionSpaceDelaunayDisplacement',componentSolutionSpaceDelaunayDisplacement,...
        'BoxSolutionSpaceMass',solutionSpaceBoxMass,...
        'ComponentSolutionSpaceConvexMass',componentSolutionSpaceConvexMass,...
        'ComponentSolutionSpaceDelaunayMass',componentSolutionSpaceDelaunayMass);
    save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialTrussBoxComponentWithMass']);
    save_print_figure(gcf,[saveFolder,'InitialTrussBoxComponentWithMass'],'PrintFormat',{'png','pdf'});
end

if(computeDisplacementAndMass)
    % initial truss + component solution spaces
    plot_results_truss_generic_moving_node(is3dPlot,systemParameter,...
        'NodePositionInitial',initialDesign,...
        'ComponentSolutionSpaceConvexDisplacementAndMass',componentSolutionSpaceConvexDisplacementAndMass,...
        'ComponentSolutionSpaceDelaunayDisplacementAndMass',componentSolutionSpaceDelaunayDisplacementAndMass);
    save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialTrussComponentAndMass']);
    save_print_figure(gcf,[saveFolder,'InitialTrussComponentAndMass'],'PrintFormat',{'png','pdf'});
    
    % initial truss + box solution space + component solution spaces
    plot_results_truss_generic_moving_node(is3dPlot,systemParameter,...
        'NodePositionInitial',initialDesign,...
        'BoxSolutionSpaceDisplacement',solutionSpaceBoxDisplacementAndMass,...
        'ComponentSolutionSpaceConvexDisplacement',componentSolutionSpaceConvexDisplacementAndMass,...
        'ComponentSolutionSpaceDelaunayDisplacement',componentSolutionSpaceDelaunayDisplacementAndMass);
    save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialTrussBoxComponentAndMass']);
    save_print_figure(gcf,[saveFolder,'InitialTrussBoxComponentAndMass'],'PrintFormat',{'png','pdf'});

    % initial truss + sample trusses + component solution space
    randomTrussMovingNode = candidate_space_sampling_individual_feasible(componentSolutionSpaceConvexDisplacementAndMass,componentIndex,nRandomTruss);
    
    plot_results_truss_generic_moving_node(is3dPlot,systemParameter,...
        'NodePositionInitial',initialDesign,...
        'NodePositionRandomDisplacementAndMass',randomTrussMovingNode,...
        'ComponentSolutionSpaceConvexDisplacementAndMass',componentSolutionSpaceConvexDisplacementAndMass,...
        'ComponentSolutionSpaceDelaunayDisplacementAndMass',componentSolutionSpaceDelaunayDisplacementAndMass);
    save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialRandomTrussComponentAndMass']);
    save_print_figure(gcf,[saveFolder,'InitialRandomTrussComponentAndMass'],'PrintFormat',{'png','pdf'});
end


%% plot candidate spaces samples



%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;


%% subfunction: plot truss
% plot additional random truss with deformation
nRandomTruss = 1;
randomTrussMovingNode = candidate_space_sampling_individual_feasible(componentSolutionSpaceConvexDisplacementAndMass,componentIndex,nRandomTruss);

% deformed trusses
nodePositionInitial = systemParameter.BaseNodePosition;
nodePositionInitial(isDesignVariable) = column_vector_to_row_major_matrix(randomTrussMovingNode',nDimension);
nodeDisplacementRandom = ...
    truss_analysis(...
	    nodePositionInitial,...
	    systemParameter.FixedDegreesOfFreedom,...
	    systemParameter.NodeForce,...
	    systemParameter.NodeElement,...
	    systemParameter.ElementCrossSectionArea,...
	    systemParameter.ElementYoungsModulus);

plot_results_truss_generic_2d_moving_node(systemParameter,...
    'NodePositionInitial',initialDesign,...
    'NodePositionOptimalDisplacement',randomTrussMovingNode,...
    'DeformationInitial',nodeDisplacementInitial,...
    'DeformationOptimalDisplacement',nodeDisplacementRandom, ...
    'ComponentSolutionSpaceConvexDisplacement',componentSolutionSpaceConvexDisplacementAndMass,...
    'OptimalTrussDisplacementOptions',{'ColorUndeformed',rand(1,3),'DisplacementScaleFactor',20},...
    'InitialTrussOptions',{'DisplacementScaleFactor',20});
save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialRandomTrussDeformationComponent3']);
save_print_figure(gcf,[saveFolder,'InitialRandomTrussDeformationComponent3'],'PrintFormat',{'png','pdf'});


function figureHandle = plot_results_truss_generic_moving_node(is3dPlot,systemParameter,varargin)
    if(~is3dPlot)
        figureHandle = plot_results_truss_generic_2d_moving_node(systemParameter,varargin{:});
    else
        figureHandle = plot_results_truss_generic_3d_moving_node(systemParameter,varargin{:});
        axis('equal','vis3d');
        camproj('perspective');
        cameratoolbar; % better adjust angle/perspective
    end

    if(nargout<1)
        clear figureHandle;
    end
end

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
    parser.addParameter('ComponentSolutionSpaceConvexDisplacement',[]);
    parser.addParameter('ComponentSolutionSpaceDelaunayDisplacement',[]);
    parser.addParameter('BoxSolutionSpaceMass',[]);
    parser.addParameter('ComponentSolutionSpaceConvexMass',[]);
    parser.addParameter('ComponentSolutionSpaceDelaunayMass',[]);
    parser.addParameter('BoxSolutionSpaceDisplacementAndMass',[]);
    parser.addParameter('ComponentSolutionSpaceConvexDisplacementAndMass',[]);
    parser.addParameter('ComponentSolutionSpaceDelaunayDisplacementAndMass',[]);
    parser.addParameter('IncludeWall',true);
    parser.addParameter('IncludeAppliedForce',true);
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
    parser.addParameter('ComponentSolutionSpaceConvexDisplacementOptions',{});
    parser.addParameter('ComponentSolutionSpaceDelaunayDisplacementOptions',{});
    parser.addParameter('BoxSolutionSpaceMassOptions',{});
    parser.addParameter('ComponentSolutionSpaceConvexMassOptions',{});
    parser.addParameter('ComponentSolutionSpaceDelaunayMassOptions',{});
    parser.addParameter('BoxSolutionSpaceDisplacementAndMassOptions',{});
    parser.addParameter('ComponentSolutionSpaceConvexDisplacementAndMassOptions',{});
    parser.addParameter('ComponentSolutionSpaceDelaunayDisplacementAndMassOptions',{});
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

    defaultComponentSolutionConvexDisplacementOptions = {'EdgeColor','g','FaceColor','none','FaceAlpha',0.5,'Linewidth',2.0};
    [~,componentSolutionConvexDisplacementOptions] = merge_name_value_pair_argument(defaultComponentSolutionConvexDisplacementOptions,options.ComponentSolutionSpaceConvexDisplacementOptions);

    defaultComponentSolutionDelaunayDisplacementOptions = {'EdgeColor','m','FaceColor','m','FaceAlpha',0.1,'Linewidth',2.0};
    [~,componentSolutionDelaunayDisplacementOptions] = merge_name_value_pair_argument(defaultComponentSolutionDelaunayDisplacementOptions,options.ComponentSolutionSpaceDelaunayDisplacementOptions);

    defaultComponentSolutionConvexMassOptions = {'EdgeColor','k','FaceColor','none','FaceAlpha',0.5,'Linewidth',2.0};
    [~,componentSolutionConvexMassOptions] = merge_name_value_pair_argument(defaultComponentSolutionConvexMassOptions,options.ComponentSolutionSpaceConvexMassOptions);

    defaultComponentSolutionDelaunayMassOptions = {'EdgeColor','y','FaceColor','y','FaceAlpha',0.1,'Linewidth',2.0};
    [~,componentSolutionDelaunayMassOptions] = merge_name_value_pair_argument(defaultComponentSolutionDelaunayMassOptions,options.ComponentSolutionSpaceDelaunayMassOptions);

    defaultComponentSolutionConvexDisplacementAndMassOptions = {'EdgeColor','g','FaceColor','none','FaceAlpha',0.5,'Linewidth',2.0};
    [~,componentSolutionConvexDisplacementAndMassOptions] = merge_name_value_pair_argument(defaultComponentSolutionConvexDisplacementAndMassOptions,options.ComponentSolutionSpaceConvexDisplacementAndMassOptions);

    defaultComponentSolutionDelaunayDisplacementAndMassOptions = {'EdgeColor','m','FaceColor','m','FaceAlpha',0.1,'Linewidth',2.0};
    [~,componentSolutionDelaunayDisplacementAndMassOptions] = merge_name_value_pair_argument(defaultComponentSolutionDelaunayDisplacementAndMassOptions,options.ComponentSolutionSpaceDelaunayDisplacementAndMassOptions);

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

    % component solution space - convex - displacement
    handleToleranceNodeDisplacementComponentConvex = [];
    if(~isempty(options.ComponentSolutionSpaceConvexDisplacement))
        for i=1:nComponent-1
            options.ComponentSolutionSpaceConvexDisplacement(i).plot_candidate_space(gcf,componentSolutionConvexDisplacementOptions{:});
        end
        handleToleranceNodeDisplacementComponentConvex = options.ComponentSolutionSpaceConvexDisplacement(nComponent).plot_candidate_space(gcf,componentSolutionConvexDisplacementOptions{:});
    end

    % component solution space - convex - mass
    handleToleranceNodeMassComponentConvex = [];
    if(~isempty(options.ComponentSolutionSpaceConvexMass))
        for i=1:nComponent-1
            options.ComponentSolutionSpaceConvexMass(i).plot_candidate_space(gcf,componentSolutionConvexMassOptions{:});
        end
        handleToleranceNodeMassComponentConvex = options.ComponentSolutionSpaceConvexMass(nComponent).plot_candidate_space(gcf,componentSolutionConvexMassOptions{:});
    end

    % component solution space - convex - displacement and mass
    handleToleranceNodeDisplacementAndMassComponentConvex = [];
    if(~isempty(options.ComponentSolutionSpaceConvexDisplacementAndMass))
        for i=1:nComponent-1
            options.ComponentSolutionSpaceConvexDisplacementAndMass(i).plot_candidate_space(gcf,componentSolutionConvexDisplacementAndMassOptions{:});
        end
        handleToleranceNodeDisplacementAndMassComponentConvex = options.ComponentSolutionSpaceConvexDisplacementAndMass(nComponent).plot_candidate_space(gcf,componentSolutionConvexDisplacementAndMassOptions{:});
    end

    % component solution space - delaunay - displacement
    handleToleranceNodeDisplacementComponentDelaunay = [];
    if(~isempty(options.ComponentSolutionSpaceDelaunayDisplacement))
        for i=1:nComponent-1
            options.ComponentSolutionSpaceDelaunayDisplacement(i).plot_candidate_space(gcf,componentSolutionDelaunayDisplacementOptions{:});
        end
        handleToleranceNodeDisplacementComponentDelaunay = options.ComponentSolutionSpaceDelaunayDisplacement(nComponent).plot_candidate_space(gcf,componentSolutionDelaunayDisplacementOptions{:});
    end

    % component solution space - delaunay - mass
    handleToleranceNodeMassComponentDelaunay = [];
    if(~isempty(options.ComponentSolutionSpaceDelaunayMass))
        for i=1:nComponent-1
            options.ComponentSolutionSpaceDelaunayMass(i).plot_candidate_space(gcf,componentSolutionConvexMassOptions{:});
        end
        handleToleranceNodeMassComponentDelaunay = options.ComponentSolutionSpaceDelaunayMass(nComponent).plot_candidate_space(gcf,componentSolutionDelaunayMassOptions{:});
    end

    % component solution space - delaunay - displacement and mass
    handleToleranceNodeDisplacementAndMassComponentDelaunay = [];
    if(~isempty(options.ComponentSolutionSpaceDelaunayDisplacementAndMass))
        for i=1:nComponent-1
            options.ComponentSolutionSpaceDelaunayDisplacementAndMass(i).plot_candidate_space(gcf,componentSolutionDelaunayDisplacementAndMassOptions{:});
        end
        handleToleranceNodeDisplacementAndMassComponentDelaunay = options.ComponentSolutionSpaceDelaunayDisplacementAndMass(nComponent).plot_candidate_space(gcf,componentSolutionDelaunayDisplacementAndMassOptions{:});
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
        positionTip = systemParameter.BaseNodePosition(any(isTrussTip,2),:);
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
            handleToleranceNodeDisplacementComponentConvex,...
            handleToleranceNodeDisplacementComponentDelaunay,...
            handleToleranceNodeMassBox,...
            handleToleranceNodeMassComponentConvex,...
            handleToleranceNodeMassComponentDelaunay,...
            handleToleranceNodeDisplacementAndMassBox,...
            handleToleranceNodeDisplacementAndMassComponentConvex,...
            handleToleranceNodeDisplacementAndMassComponentDelaunay};

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
            'Tolerance Region for the Node (Component) - Displacement',...
            'Tolerance Region for the Node (Box) - Mass',...
            'Tolerance Region for the Node (Component) - Mass',...
            'Tolerance Region for the Node (Component) - Mass',...
            'Tolerance Region for the Node (Box) - Displacement + Mass',...
            'Tolerance Region for the Node (Component) - Displacement + Mass',...
            'Tolerance Region for the Node (Component) - Displacement + Mass'};
        
        handleObject = [handleObjectAll{:}];
        legentText = {legendTextAll{~cellfun(@isempty,handleObjectAll)}};

        legend(handleObject,legentText,legendOptions{:});
    end

    if(nargout<1)
        clear figureHandle;
    end
end

function figureHandle = plot_results_truss_generic_3d_moving_node(systemParameter,varargin)
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
    parser.addParameter('ComponentSolutionSpaceConvexDisplacement',[]);
    parser.addParameter('ComponentSolutionSpaceDelaunayDisplacement',[]);
    parser.addParameter('BoxSolutionSpaceMass',[]);
    parser.addParameter('ComponentSolutionSpaceConvexMass',[]);
    parser.addParameter('ComponentSolutionSpaceDelaunayMass',[]);
    parser.addParameter('BoxSolutionSpaceDisplacementAndMass',[]);
    parser.addParameter('ComponentSolutionSpaceConvexDisplacementAndMass',[]);
    parser.addParameter('ComponentSolutionSpaceDelaunayDisplacementAndMass',[]);
    parser.addParameter('IncludeWall',true);
    parser.addParameter('IncludeAppliedForce',true);
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
    parser.addParameter('ComponentSolutionSpaceConvexDisplacementOptions',{});
    parser.addParameter('ComponentSolutionSpaceDelaunayDisplacementOptions',{});
    parser.addParameter('BoxSolutionSpaceMassOptions',{});
    parser.addParameter('ComponentSolutionSpaceConvexMassOptions',{});
    parser.addParameter('ComponentSolutionSpaceDelaunayMassOptions',{});
    parser.addParameter('BoxSolutionSpaceDisplacementAndMassOptions',{});
    parser.addParameter('ComponentSolutionSpaceConvexDisplacementAndMassOptions',{});
    parser.addParameter('ComponentSolutionSpaceDelaunayDisplacementAndMassOptions',{});
    parser.addParameter('LegendOptions',{});
    parser.parse(varargin{:});
    options = parser.Results;

    [wallY,wallZ] = meshgrid(-0.5:1.5);
    wallX = zeros(size(wallY));

    defaultInitialTrussOptions = {'ColorUndeformed',[0.8 0.8 0.8],'ColorDeformed','c','MaximumLinewidth',4.0,'DisplacementScaleFactor',100};
    [~,initialTrussOptions] = merge_name_value_pair_argument(defaultInitialTrussOptions,options.InitialTrussOptions);

    defaultOptimalTrussDisplacementOptions = {'ColorUndeformed','b','ColorDeformed','m','MaximumLinewidth',4.0,'DisplacementScaleFactor',100};
    [~,optimalTrussDisplacementOptions] = merge_name_value_pair_argument(defaultOptimalTrussDisplacementOptions,options.OptimalTrussDisplacementOptions);

    defaultOptimalTrussMassOptions = {'ColorUndeformed','y','ColorDeformed','k','MaximumLinewidth',4.0,'DisplacementScaleFactor',100};
    [~,optimalTrussMassOptions] = merge_name_value_pair_argument(defaultOptimalTrussMassOptions,options.OptimalTrussMassOptions);

    defaultRandomTrussDisplacementOptions = {'MaximumLinewidth',4.0};
    [~,randomTrussDisplacementOptions] = merge_name_value_pair_argument(defaultRandomTrussDisplacementOptions,options.RandomTrussDisplacementOptions);

    defaultRandomTrussMassOptions = {'MaximumLinewidth',4.0};
    [~,randomTrussMassOptions] = merge_name_value_pair_argument(defaultRandomTrussMassOptions,options.RandomTrussMassOptions);

    defaultRandomTrussDisplacementAndMassOptions = {'MaximumLinewidth',4.0};
    [~,randomTrussDisplacementAndMassOptions] = merge_name_value_pair_argument(defaultRandomTrussDisplacementAndMassOptions,options.RandomTrussDisplacementAndMassOptions);

    defaultWallOptions = {'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.9};
    [~,wallOptions] = merge_name_value_pair_argument(defaultWallOptions,options.WallOptions);

    defaultAppliedForceOptions = {'Color','r','LineWidth',3.0};
    [~,appliedForceOptions] = merge_name_value_pair_argument(defaultAppliedForceOptions,options.AppliedForceOptions);

    defaultBoxSolutionDisplacementOptions = {'FaceColor','c','FaceAlpha',0.2};
    [~,boxSolutionDisplacementOptions] = merge_name_value_pair_argument(defaultBoxSolutionDisplacementOptions,options.BoxSolutionSpaceDisplacementOptions);

    defaultBoxSolutionMassOptions = {'FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',0.2};
    [~,boxSolutionMassOptions] = merge_name_value_pair_argument(defaultBoxSolutionMassOptions,options.BoxSolutionSpaceMassOptions);

    defaultBoxSolutionDisplacementAndMassOptions = {'FaceColor','c','FaceAlpha',0.2};
    [~,boxSolutionDisplacementAndMassOptions] = merge_name_value_pair_argument(defaultBoxSolutionDisplacementAndMassOptions,options.BoxSolutionSpaceDisplacementAndMassOptions);

    defaultComponentSolutionConvexDisplacementOptions = {'FaceColor','g','FaceAlpha',0.2};
    [~,componentSolutionConvexDisplacementOptions] = merge_name_value_pair_argument(defaultComponentSolutionConvexDisplacementOptions,options.ComponentSolutionSpaceConvexDisplacementOptions);

    defaultComponentSolutionConvexMassOptions = {'FaceColor','k','FaceAlpha',0.2};
    [~,componentSolutionConvexMassOptions] = merge_name_value_pair_argument(defaultComponentSolutionConvexMassOptions,options.ComponentSolutionSpaceConvexMassOptions);

    defaultComponentSolutionConvexDisplacementAndMassOptions = {'FaceColor','g','FaceAlpha',0.2};
    [~,componentSolutionConvexDisplacementAndMassOptions] = merge_name_value_pair_argument(defaultComponentSolutionConvexDisplacementAndMassOptions,options.ComponentSolutionSpaceConvexDisplacementAndMassOptions);

    defaultComponentSolutionDelaunayDisplacementOptions = {'FaceColor','m','FaceAlpha',0.2};
    [~,componentSolutionDelaunayDisplacementOptions] = merge_name_value_pair_argument(defaultComponentSolutionDelaunayDisplacementOptions,options.ComponentSolutionSpaceDelaunayDisplacementOptions);

    defaultComponentSolutionDelaunayMassOptions = {'FaceColor','m','FaceAlpha',0.2};
    [~,componentSolutionDelaunayMassOptions] = merge_name_value_pair_argument(defaultComponentSolutionDelaunayMassOptions,options.ComponentSolutionSpaceDelaunayMassOptions);

    defaultComponentSolutionDelaunayDisplacementAndMassOptions = {'FaceColor','m','FaceAlpha',0.2};
    [~,componentSolutionDelaunayDisplacementAndMassOptions] = merge_name_value_pair_argument(defaultComponentSolutionDelaunayDisplacementAndMassOptions,options.ComponentSolutionSpaceDelaunayDisplacementAndMassOptions);

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
            plot_design_box_3d(gcf,options.BoxSolutionSpaceDisplacement(:,currentIndex),boxSolutionDisplacementOptions{:});
        end
        currentIndex = 1 + nDimension*(nComponent-1) + [0:(nDimension-1)];
        handleToleranceNodeDisplacementBox = plot_design_box_3d(gcf,options.BoxSolutionSpaceDisplacement(:,currentIndex),boxSolutionDisplacementOptions{:});
    end

    % box-shaped solution space - mass
    handleToleranceNodeMassBox = [];
    if(~isempty(options.BoxSolutionSpaceMass))
        for i=1:nComponent-1
            currentIndex = 1 + nDimension*(i-1) + [0:(nDimension-1)];
            plot_design_box_3d(gcf,options.BoxSolutionSpaceMass(:,currentIndex),boxSolutionMassOptions{:});
        end
        currentIndex = 1 + nDimension*(nComponent-1) + [0:(nDimension-1)];
        handleToleranceNodeMassBox = plot_design_box_3d(gcf,options.BoxSolutionSpaceMass(:,currentIndex),boxSolutionMassOptions{:});
    end

    % box-shaped solution space  - displacement and mass
    handleToleranceNodeDisplacementAndMassBox = [];
    if(~isempty(options.BoxSolutionSpaceDisplacementAndMass))
        for i=1:nComponent-1
            currentIndex = 1 + nDimension*(i-1) + [0:(nDimension-1)];
            plot_design_box_3d(gcf,options.BoxSolutionSpaceDisplacementAndMass(:,currentIndex),boxSolutionDisplacementAndMassOptions{:});
        end
        currentIndex = 1 + nDimension*(nComponent-1) + [0:(nDimension-1)];
        handleToleranceNodeDisplacementAndMassBox = plot_design_box_3d(gcf,options.BoxSolutionSpaceDisplacementAndMass(:,currentIndex),boxSolutionDisplacementAndMassOptions{:});
    end

    % component solution space - convex - displacement
    handleToleranceNodeDisplacementComponentConvex = [];
    if(~isempty(options.ComponentSolutionSpaceConvexDisplacement))
        for i=1:nComponent-1
            options.ComponentSolutionSpaceConvexDisplacement(i).plot_candidate_space(gcf,componentSolutionConvexDisplacementOptions{:});
        end
        handleToleranceNodeDisplacementComponentConvex = options.ComponentSolutionSpaceConvexDisplacement(nComponent).plot_candidate_space(gcf,componentSolutionConvexDisplacementOptions{:});
    end

    % component solution space - convex - mass
    handleToleranceNodeMassComponentConvex = [];
    if(~isempty(options.ComponentSolutionSpaceConvexMass))
        for i=1:nComponent-1
            options.ComponentSolutionSpaceConvexMass(i).plot_candidate_space(gcf,componentSolutionConvexMassOptions{:});
        end
        handleToleranceNodeMassComponentConvex = options.ComponentSolutionSpaceConvexMass(nComponent).plot_candidate_space(gcf,componentSolutionConvexMassOptions{:});
    end

    % component solution space - convex - displacement and mass
    handleToleranceNodeDisplacementAndMassComponentConvex = [];
    if(~isempty(options.ComponentSolutionSpaceConvexDisplacementAndMass))
        for i=1:nComponent-1
            options.ComponentSolutionSpaceConvexDisplacementAndMass(i).plot_candidate_space(gcf,componentSolutionConvexDisplacementAndMassOptions{:});
        end
        handleToleranceNodeDisplacementAndMassComponentConvex = options.ComponentSolutionSpaceConvexDisplacementAndMass(nComponent).plot_candidate_space(gcf,componentSolutionConvexDisplacementAndMassOptions{:});
    end

    % component solution space - delaunay - displacement
    handleToleranceNodeDisplacementComponentDelaunay = [];
    if(~isempty(options.ComponentSolutionSpaceDelaunayDisplacement))
        for i=1:nComponent-1
            options.ComponentSolutionSpaceDelaunayDisplacement(i).plot_candidate_space(gcf,componentSolutionDelaunayDisplacementOptions{:});
        end
        handleToleranceNodeDisplacementComponentDelaunay = options.ComponentSolutionSpaceDelaunayDisplacement(nComponent).plot_candidate_space(gcf,componentSolutionDelaunayDisplacementOptions{:});
    end

    % component solution space - delaunay - mass
    handleToleranceNodeMassComponentDelaunay = [];
    if(~isempty(options.ComponentSolutionSpaceDelaunayMass))
        for i=1:nComponent-1
            options.ComponentSolutionSpaceDelaunayMass(i).plot_candidate_space(gcf,componentSolutionDelaunayMassOptions{:});
        end
        handleToleranceNodeMassComponentDelaunay = options.ComponentSolutionSpaceDelaunayMass(nComponent).plot_candidate_space(gcf,componentSolutionDelaunayMassOptions{:});
    end

    % component solution space - delaunay - displacement and mass
    handleToleranceNodeDisplacementAndMassComponentDelaunay = [];
    if(~isempty(options.ComponentSolutionSpaceDelaunayDisplacementAndMass))
        for i=1:nComponent-1
            options.ComponentSolutionSpaceDelaunayDisplacementAndMass(i).plot_candidate_space(gcf,componentSolutionDelaunayDisplacementAndMassOptions{:});
        end
        handleToleranceNodeDisplacementAndMassComponentDelaunay = options.ComponentSolutionSpaceDelaunayDisplacementAndMass(nComponent).plot_candidate_space(gcf,componentSolutionDelaunayDisplacementAndMassOptions{:});
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
        handleWall = surf(wallX,wallY,wallZ,wallOptions{:});
    end

    % applied force
    handleForce = [];
    if(options.IncludeAppliedForce)
        positionTip = systemParameter.BaseNodePosition(any(isTrussTip,2),:);
        handleForce = quiver3(positionTip(1),positionTip(2),positionTip(3),0,0,-0.5,appliedForceOptions{:});
    end

    % axis
    grid minor;
    if(~options.IncludeAxesInformation)
        set(gca,'XColor', 'none','YColor','none','ZColor','none');
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
            handleToleranceNodeDisplacementComponentConvex,...
            handleToleranceNodeDisplacementComponentDelaunay,...
            handleToleranceNodeMassBox,...
            handleToleranceNodeMassComponentConvex,...
            handleToleranceNodeMassComponentDelaunay,...
            handleToleranceNodeDisplacementAndMassBox,...
            handleToleranceNodeDisplacementAndMassComponentConvex,...
            handleToleranceNodeDisplacementAndMassComponentDelaunay};

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
            'Tolerance Region for the Node (Component) - Displacement',...
            'Tolerance Region for the Node (Box) - Mass',...
            'Tolerance Region for the Node (Component) - Mass',...
            'Tolerance Region for the Node (Component) - Mass',...
            'Tolerance Region for the Node (Box) - Displacement + Mass',...
            'Tolerance Region for the Node (Component) - Displacement + Mass',...
            'Tolerance Region for the Node (Component) - Displacement + Mass'};
        
        handleObject = [handleObjectAll{:}];
        legentText = {legendTextAll{~cellfun(@isempty,handleObjectAll)}};

        legend(handleObject,legentText,legendOptions{:});
    end

    if(nargout<1)
        clear figureHandle;
    end
end

function save_3d_rotating_video_gif(is3dPlot,figureHandle,filename,varargin)
    if(~is3dPlot)
        return;
    end

    [defaultAzimuth,defaultElevation] = view(3);
    % rotate so truss can be better seen at the start
    rotatedAzimuth = defaultAzimuth + 60;
    view(rotatedAzimuth,defaultElevation);
    rotating_video(figureHandle,filename,varargin{:});
    rotating_gif(figureHandle,filename,varargin{:});
end

function rotating_video(figureHandle,filename,varargin)
    parser = inputParser;
    parser.addParameter('Duration',5);
    parser.addParameter('FramesPerSecond',60);
    parser.addParameter('VideoWriterProfile','Motion JPEG 2000');
    parser.addParameter('VideoWriterOptions',{'LosslessCompression',true});
    parser.parse(varargin{:});
    options = parser.Results;

    defaultVideoWriterOptions = {'FrameRate',options.FramesPerSecond};
    [~,videoWriterOptions] = merge_name_value_pair_argument(defaultVideoWriterOptions,options.VideoWriterOptions);

    % activate figure
    figure(figureHandle);
    [defaultAzimuth,defaultElevation] = view;

    % prepare angles for each frame
    totalFrame = options.FramesPerSecond*options.Duration;
    azimuthVideo = mod(defaultAzimuth + linspace(0,360,totalFrame),360);

    % initialize video object + process options
    videoHandle = VideoWriter(filename,options.VideoWriterProfile);
    for i=1:2:length(videoWriterOptions)
        videoHandle.(videoWriterOptions{i}) = videoWriterOptions{i+1};
    end

    % make video with each frame
    open(videoHandle);
    for i=1:totalFrame
        view(azimuthVideo(i),defaultElevation);
        currentFrame = getframe(gcf);
        writeVideo(videoHandle,currentFrame);
    end
    close(videoHandle);
end

function rotating_gif(figureHandle,filename,varargin)
    parser = inputParser;
    parser.addParameter('Duration',5);
    parser.addParameter('FramesPerSecond',30);
    parser.addParameter('ImwriteOptions',{})
    parser.parse(varargin{:});
    options = parser.Results;

    defaultImwriteOptions = {};
    [~,imwriteOptions] = merge_name_value_pair_argument(defaultImwriteOptions,options.ImwriteOptions);

    % activate figure
    figure(figureHandle);
    set(figureHandle, 'MenuBar', 'none');
    set(figureHandle, 'ToolBar', 'none');
    [defaultAzimuth,defaultElevation] = view;

    % prepare angles for each frame
    totalFrame = options.FramesPerSecond*options.Duration;
    azimuthVideo = mod(defaultAzimuth + linspace(0,360,totalFrame),360);
    frameDelay = 1/options.FramesPerSecond;

    % make gif with each frame
    for i=1:totalFrame
        view(azimuthVideo(i),defaultElevation);
        [A,map] = rgb2ind(frame2im(getframe(gcf)),256);

        if(i==1)
            imwrite(A,map,...
                [filename,'.gif'],'gif',...
                'LoopCount',inf,...
                'DelayTime',frameDelay,...
                imwriteOptions{:});
        else
            imwrite(A,map,...
                [filename,'.gif'],'gif',...
                'WriteMode','append', ...
                'DelayTime',frameDelay,...
                imwriteOptions{:});
        end
    end

    % reset menu/tool bars
    set(figureHandle, 'MenuBar', 'figure');
    set(figureHandle, 'ToolBar', 'auto');
end