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
useBoxResultForComponent = false;

computeDisplacement = true;
computeMass = false;
computeDisplacementAndMass = false;

computePlanarTrimmingComponent = true;
computeCornerBoxRemovalComponent = true;
computeHolePunchingComponent = false;


%% function call
systemFunction = @truss_generic_moving_node;

switch trussAnalysisChoice
    case '2-DoF-2D'
        nSample = 100;
        maxIterDisplacement = 50;
        maxIterMass = 20;
        maxIterDisplacementAndMass = 50;
        growthRateDisplacement = 0.07;
        growthRateMass = 0.07;
        trimmingPasses = 'full';
        nRandomTruss = 1;
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
        cameraPositionFigureSave = [];
    case '4-DoF-2D'
        nSample = 100;
        maxIterDisplacement = 30;
        maxIterMass = 30;
        maxIterDisplacementAndMass = 35;
        growthRateDisplacement = 0.06;
        growthRateMass = 0.06;
        trimmingPasses = 'full';
        nRandomTruss = 3;
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
        cameraPositionFigureSave = [];
    case '16-DoF-2D'
        nSample = 100;
        maxIterDisplacement = 250;
        maxIterMass = 250;
        maxIterDisplacementAndMass = 250;
        growthRateDisplacement = 0.01;
        growthRateMass = 0.01;
        trimmingPasses = 'reduced';
        nRandomTruss = 0;
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
        cameraPositionFigureSave = [];
    case '36-DoF-2D'
        nSample = 100;
        maxIterDisplacement = 300;
        maxIterMass = 300;
        maxIterDisplacementAndMass = 300;
        growthRateDisplacement = 0.007;
        growthRateMass = 0.007;
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
        cameraPositionFigureSave = [];
    case '9-DoF-3D'
        nSample = 100;
        maxIterDisplacement = 50;
        maxIterMass = 100;
        maxIterDisplacementAndMass = 50;
        growthRateDisplacement = 0.04;
        growthRateMass = 0.04;
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
        cameraPositionFigureSave = [4.7185  -11.2817    7.6322];
    case '36-DoF-3D'
        nSample = 100;
        maxIterDisplacement = 300;
        maxIterMass = 300;
        maxIterDisplacementAndMass = 300;
        growthRateDisplacement = 0.004;
        growthRateMass = 0.007;
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
        cameraPositionFigureSave = [7.3696  -14.9327    9.8425];
end


%% starting process - compute necessary parts, visualize initial truss
% component generation - assume all position variables of each node form a component
nDimension = size(systemParameter.BaseNodePosition,2);
isDesignVariable = isnan(systemParameter.BaseNodePosition);
nComponent = sum(isDesignVariable(:,1),1);
componentIndex = cell(1,nComponent);
for i=1:nComponent
    componentIndex{i} = (1 + nDimension*(i-1) + [0:(nDimension-1)])';
end
is3dPlot = (nDimension==3);


%% Analyse initial truss
plot_results_truss_generic_moving_node(systemParameter,...
    'NodePositionInitial',initialDesign);
save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialTruss'],cameraPositionFigureSave);
save_print_figure(gcf,[saveFolder,'InitialTruss'],'PrintFormat',{'png','pdf'},'Size',figureSize);

% deformed nitial truss
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

plot_results_truss_generic_moving_node(systemParameter,...
    'NodePositionInitial',initialDesign,...
    'DeformationInitial',nodeDisplacementInitial);
save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialTrussDeformation'],cameraPositionFigureSave);
save_print_figure(gcf,[saveFolder,'InitialTrussDeformation'],'PrintFormat',{'png','pdf'},'Size',figureSize);


%% establish upper performance limits
nElement = size(systemParameter.NodeElement,1);
performanceLowerLimit = [  0   0 repmat(-inf,1,2*nElement)];
performanceUpperLimit = [nan nan repmat(inf,1,2*nElement)];

% update uppwer limit based on either optimal value or initial value
bottomUpMapping = BottomUpMappingFunction(systemFunction,'SystemParameter',systemParameter);
performanceMeasureInitial = bottomUpMapping.response(initialDesign);


%% find optimum
% compute optimal displacement
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

    % initial truss + optimized truss
    plot_results_truss_generic_moving_node(systemParameter,...
        'NodePositionInitial',initialDesign,...
        'NodePositionOptimalDisplacement',nodePositionOptimalDisplacement);
    save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialOptimizedTrussDisplacement'],cameraPositionFigureSave);
    save_print_figure(gcf,[saveFolder,'InitialOptimizedTrussDisplacement'],'PrintFormat',{'png','pdf'},'Size',figureSize);
end

% compute optimal mass
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

    % initial truss + optimized displacment truss + optimized mass truss
    plot_results_truss_generic_moving_node(systemParameter,...
        'NodePositionInitial',initialDesign,...
        'NodePositionOptimalMass',nodePositionOptimalMass);
    save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialOptimizedTrussMass'],cameraPositionFigureSave);
    save_print_figure(gcf,[saveFolder,'InitialOptimizedTrussMass'],'PrintFormat',{'png','pdf'},'Size',figureSize);
end

% compute optimal displacement and mass
if(computeDisplacementAndMass)
    performanceUpperLimitDisplacementAndMass = performanceUpperLimit;
    performanceUpperLimitDisplacementAndMass(1) = performanceMeasureInitial(1);
    performanceUpperLimitDisplacementAndMass(2) = performanceMeasureInitial(2);
    designEvaluatorDisplacementAndMass = DesignEvaluatorBottomUpMapping(...
        bottomUpMapping,...
        performanceLowerLimit,...
        performanceUpperLimitDisplacementAndMass);

    warning('off');
    [nodePositionOptimalDisplacementAndMass,massOptimal] = design_optimize_performance_score(...
        designEvaluatorDisplacementAndMass,...
        initialDesign,...
        designSpaceLowerBoundDisplacement,...
        designSpaceUpperBoundDisplacement,...
        'OptimizationMethodFunction',@optimization_ga_wrapper,...
        'OptimizationMethodOptions',{'Display','diagnose'});
    warning('on');

    plot_results_truss_generic_moving_node(systemParameter,...
        'NodePositionInitial',initialDesign,...
        'NodePositionOptimalDisplacement',nodePositionOptimalDisplacementAndMass);
    save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialOptimizedTrussDisplacementAndMass'],cameraPositionFigureSave);
    save_print_figure(gcf,[saveFolder,'InitialOptimizedTrussDisplacementAndMass'],'PrintFormat',{'png','pdf'},'Size',figureSize);
end

% display deformed optimal displacement
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

    plot_results_truss_generic_moving_node(systemParameter,...
        'NodePositionInitial',initialDesign,...
        'NodePositionOptimalDisplacement',nodePositionOptimalDisplacement,...
        'DeformationInitial',nodeDisplacementInitial,...
        'DeformationOptimalDisplacement',nodeDisplacementOptimal);
    save_3d_rotating_video_gif(is3dPlot,gcf,[saveFolder,'InitialOptimizedTrussDeformation'],cameraPositionFigureSave);
    save_print_figure(gcf,[saveFolder,'InitialOptimizedTrussDeformation'],'PrintFormat',{'png','pdf'},'Size',figureSize);
end


%% Solve solution spaces problems
% displacement
if(computeDisplacement)
    performanceUpperLimitDisplacement = performanceUpperLimit;
    performanceUpperLimitDisplacement(1) = performanceMeasureInitial(1);
    designEvaluatorDisplacement = DesignEvaluatorBottomUpMapping(...
        bottomUpMapping,...
        performanceLowerLimit,...
        performanceUpperLimitDisplacement);

    [solutionSpaceBoxDisplacement,componentSolutionSpacePlanarTrimmingDisplacement,componentSolutionSpaceCornerBoxRemovalDisplacement,componentSolutionSpaceHolePunchingDisplacement] = ...
    compute_truss_solution_spaces('Displacement',designEvaluatorDisplacement,initialDesign,designSpaceLowerBoundDisplacement,...
        designSpaceUpperBoundDisplacement,componentIndex,nSample,maxIterDisplacement,growthRateDisplacement,trimmingPasses,...
        useBoxResultForComponent,computePlanarTrimmingComponent,computeCornerBoxRemovalComponent,computeHolePunchingComponent,...
        rngState,saveFolder,figureSize);

    plot_relevant_results_truss_moving_node('Displacement',systemParameter,initialDesign,nodePositionOptimalDisplacement,...
        solutionSpaceBoxDisplacement,componentSolutionSpacePlanarTrimmingDisplacement,componentSolutionSpaceCornerBoxRemovalDisplacement,componentSolutionSpaceHolePunchingDisplacement,...
        componentIndex,nRandomTruss,saveFolder,figureSize,cameraPositionFigureSave);
end

% mass
if(computeMass)
    performanceUpperLimitMass = performanceUpperLimit;
    performanceUpperLimitMass(2) = performanceMeasureInitial(2);
    designEvaluatorMass = DesignEvaluatorBottomUpMapping(...
        bottomUpMapping,...
        performanceLowerLimit,...
        performanceUpperLimitMass);

    [solutionSpaceBoxMass,componentSolutionSpacePlanarTrimmingMass,componentSolutionSpaceCornerBoxRemovalMass,componentSolutionSpaceHolePunchingMass] = ...
        compute_truss_solution_spaces('Mass',designEvaluatorMass,initialDesign,designSpaceLowerBoundMass,designSpaceUpperBoundMass,...
            componentIndex,nSample,maxIterMass,growthRateMass,trimmingPasses,useBoxResultForComponent,...
            computePlanarTrimmingComponent,computeCornerBoxRemovalComponent,computeHolePunchingComponent,...
            rngState,saveFolder,figureSize);

    plot_relevant_results_truss_moving_node('Mass',systemParameter,initialDesign,nodePositionOptimalMass,...
        solutionSpaceBoxMass,componentSolutionSpacePlanarTrimmingMass,componentSolutionSpaceCornerBoxRemovalMass,componentSolutionSpaceHolePunchingMass,...
        componentIndex,nRandomTruss,saveFolder,figureSize,cameraPositionFigureSave);
end

% displacement + mass
if(computeDisplacementAndMass)
    [solutionSpaceBoxDisplacementAndMass,componentSolutionSpacePlanarTrimmingDisplacementAndMass,componentSolutionSpaceCornerBoxRemovalDisplacementAndMass,componentSolutionSpaceHolePunchingDisplacementAndMass] = ...
    compute_truss_solution_spaces('DisplacementAndMass',designEvaluatorDisplacementAndMass,initialDesign,...
        designSpaceLowerBoundDisplacement,designSpaceUpperBoundDisplacement,componentIndex,nSample,maxIterDisplacementAndMass,...
        growthRateDisplacement,trimmingPasses,useBoxResultForComponent,...
        computePlanarTrimmingComponent,computeCornerBoxRemovalComponent,computeHolePunchingComponent,...
        rngState,saveFolder,figureSize);

    plot_relevant_results_truss_moving_node('DisplacementAndMass',systemParameter,initialDesign,nodePositionOptimalDisplacementAndMass,...
        solutionSpaceBoxDisplacementAndMass,componentSolutionSpacePlanarTrimmingDisplacementAndMass,componentSolutionSpaceCornerBoxRemovalDisplacementAndMass,...
        componentSolutionSpaceHolePunchingDisplacementAndMass,componentIndex,nRandomTruss,saveFolder,figureSize,cameraPositionFigureSave);
end


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;


%% Subfunctions
function [solutionSpaceBox,componentSolutionSpacePlanarTrimming,componentSolutionSpaceCornerBoxRemoval,componentSolutionSpaceHolePunching] = ...
    compute_truss_solution_spaces(typeName,designEvaluator,initialDesign,designSpaceLowerBound,designSpaceUpperBound,componentIndex,...
        nSample,maxIter,growthRate,trimmingPasses,useBoxResultForComponent,...
        computePlanarTrimmingComponent,computeCornerBoxRemovalComponent,computeHolePunchingComponent,...
        rngState,saveFolder,figureSize)
    
    %% box
    timeElapsedBox = tic;
    optionsBox = sso_stochastic_options('box',...
        'NumberSamplesPerIterationExploration',nSample,...
        'NumberSamplesPerIterationConsolidation',nSample,...
        'FixIterNumberExploration',true,...
        'FixIterNumberConsolidation',true,...
        'MaxIterExploration',maxIter,...
        'MaxIterConsolidation',maxIter,...
        'UseAdaptiveGrowthRate',true,...
        'GrowthRate',growthRate,...
        'ApplyLeanness','never',...
        'TrimmingOperationOptions',{'PassesCriterion',trimmingPasses},...
        'TrimmingOrderOptions',{'OrderPreference','score'});
    
    rng(rngState);
    [solutionSpaceBox,problemDataBox,iterDataBox] = sso_box_stochastic(designEvaluator,...
        initialDesign,designSpaceLowerBound,designSpaceUpperBound,optionsBox);
    toc(timeElapsedBox)

    resultsFolder = [saveFolder,sprintf('PerformanceBox%s/',typeName)];
    mkdir(resultsFolder);
    algoDataBox = postprocess_sso_box_stochastic(problemDataBox,iterDataBox);
    plot_sso_box_stochastic_metrics(algoDataBox,...
        'SaveFolder',resultsFolder,...
        'CloseFigureAfterSaving',true,...
        'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});

    if(useBoxResultForComponent)
        initialDesignComponent = solutionSpaceBox;
    else
        initialDesignComponent = initialDesign;
    end
    comparisonComponent = {};
    componentLabel = {};

    %% component planar trimming
    componentSolutionSpacePlanarTrimming = [];
    if(computePlanarTrimmingComponent)
        timeElapsedComponent = tic;
        optionsComponent = sso_stochastic_options('component',...
            'NumberSamplesPerIterationExploration',nSample,...
            'NumberSamplesPerIterationConsolidation',nSample,...
            'FixIterNumberExploration',true,...
            'FixIterNumberConsolidation',true,...
            'MaxIterExploration',maxIter,...
            'MaxIterConsolidation',maxIter,...
            'CandidateSpaceConstructor',@CandidateSpaceConvexHull,...
            'TrimmingMethodFunction',@component_trimming_method_planar_trimming,...
            ... 'TrimmingMethodOptions',{'ReferenceDesigns','boundary'},...
            'UseAdaptiveGrowthRate',true,...
            'GrowthRate',growthRate,...
            'ApplyLeanness','never',...
            'UsePaddingSamplesInTrimming',true,...
            'UsePreviousEvaluatedSamplesConsolidation',false,...
            'TrimmingOperationOptions',{'PassesCriterion',trimmingPasses},...
            'TrimmingOrderOptions',{'OrderPreference','score'},...
            'ShapeSamplesUsefulConsolidation',true);
        
        rng(rngState);
        [componentSolutionSpacePlanarTrimming,problemDataComponentPlanarTrimming,iterDataComponentPlanarTrimming] = sso_component_stochastic(designEvaluator,...
            initialDesignComponent,designSpaceLowerBound,designSpaceUpperBound,componentIndex,optionsComponent);
        toc(timeElapsedComponent)

        resultsFolder = [saveFolder,sprintf('PerformancePlanarTrimming%s/',typeName)];
        mkdir(resultsFolder);
        algoDataComponentPlanarTrimming = postprocess_sso_component_stochastic(problemDataComponentPlanarTrimming,iterDataComponentPlanarTrimming);
        plot_sso_component_stochastic_metrics(algoDataComponentPlanarTrimming,...
            'SaveFolder',resultsFolder,...
            'CloseFigureAfterSaving',true,...
            'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});

        comparisonComponent{end+1} = algoDataComponentPlanarTrimming;
        componentLabel{end+1} = 'Planar Trimming';
    end

    %% component corner box removal
    componentSolutionSpaceCornerBoxRemoval = [];
    if(computeCornerBoxRemovalComponent)
        timeElapsedComponent = tic;
        optionsComponent = sso_stochastic_options('component',...
            'NumberSamplesPerIterationExploration',nSample,...
            'NumberSamplesPerIterationConsolidation',nSample,...
            'FixIterNumberExploration',true,...
            'FixIterNumberConsolidation',true,...
            'MaxIterExploration',maxIter,...
            'MaxIterConsolidation',maxIter,...
            'CandidateSpaceConstructor',@CandidateSpaceCornerBoxRemoval,...
            'TrimmingMethodFunction',@component_trimming_method_corner_box_removal,...
            ...'TrimmingMethodOptions',{'CornersToTest','away'},...
            'UseAdaptiveGrowthRate',true,...
            'GrowthRate',growthRate,...
            'ApplyLeanness','never',...
            'UsePaddingSamplesInTrimming',true,...
            'UsePreviousEvaluatedSamplesConsolidation',true,...
            'TrimmingOperationOptions',{'PassesCriterion',trimmingPasses},...
            'TrimmingOrderOptions',{'OrderPreference','score'},...
            'ShapeSamplesUsefulConsolidation',false,...
            'MaximumGrowthAdaptationFactor',1.5);
        
        rng(rngState);
        [componentSolutionSpaceCornerBoxRemoval,problemDataComponentCornerBoxRemoval,iterDataComponentCornerBoxRemoval] = sso_component_stochastic(designEvaluator,...
            initialDesignComponent,designSpaceLowerBound,designSpaceUpperBound,componentIndex,optionsComponent);
        toc(timeElapsedComponent)

        resultsFolder = [saveFolder,sprintf('PerformanceCornerBoxRemoval%s/',typeName)];
        mkdir(resultsFolder);
        algoDataComponentCornerBoxRemoval = postprocess_sso_component_stochastic(problemDataComponentCornerBoxRemoval,iterDataComponentCornerBoxRemoval);
        plot_sso_component_stochastic_metrics(algoDataComponentCornerBoxRemoval,...
            'SaveFolder',resultsFolder,...
            'CloseFigureAfterSaving',true,...
            'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});

        comparisonComponent{end+1} = algoDataComponentCornerBoxRemoval;
        componentLabel{end+1} = 'Corner Box Removal';
    end

    %% component hole punching
    componentSolutionSpaceHolePunching = [];
    if(computeHolePunchingComponent)
        timeElapsedComponent = tic;
        optionsComponent = sso_stochastic_options('component',...
            'NumberSamplesPerIterationExploration',nSample,...
            'NumberSamplesPerIterationConsolidation',nSample,...
            'FixIterNumberExploration',true,...
            'FixIterNumberConsolidation',true,...
            'MaxIterExploration',maxIter,...
            'MaxIterConsolidation',maxIter,...
            'CandidateSpaceConstructor',@CandidateSpaceSvm,...
            'TrimmingMethodFunction',@component_trimming_method_hole_punching,...
            'UseAdaptiveGrowthRate',true,...
            'GrowthRate',growthRate,...
            'ApplyLeanness','never',...
            'UsePaddingSamplesInTrimming',true,...
            'UsePreviousEvaluatedSamplesConsolidation',true,...
            'TrimmingOperationOptions',{'PassesCriterion',trimmingPasses},...
            'TrimmingOrderOptions',{'OrderPreference','score'},...
            'ShapeSamplesUsefulConsolidation',false,...
            'NumberPaddingSamples',1000);
        
        rng(rngState);
        [componentSolutionSpaceHolePunching,problemDataComponentHolePunching,iterDataComponentHolePunching] = sso_component_stochastic(designEvaluator,...
            initialDesignComponent,designSpaceLowerBound,designSpaceUpperBound,componentIndex,optionsComponent);
        toc(timeElapsedComponent)

        resultsFolder = [saveFolder,sprintf('PerformanceHolePunching%s/',typeName)];
        mkdir(resultsFolder);
        algoDataComponentHolePunching = postprocess_sso_component_stochastic(problemDataComponentHolePunching,iterDataComponentHolePunching);
        plot_sso_component_stochastic_metrics(algoDataComponentHolePunching,...
            'SaveFolder',resultsFolder,...
            'CloseFigureAfterSaving',true,...
            'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});

        comparisonComponent{end+1} = algoDataComponentHolePunching;
        componentLabel{end+1} = 'Hole Punching';
    end

    % comparison
    resultsFolder = [saveFolder,sprintf('PerformanceComparison%s/',typeName)];
    mkdir(resultsFolder);
    plot_sso_comparison_box_component_stochastic_metrics(...
        {algoDataBox},...
        comparisonComponent,...
        'ComponentLabel',componentLabel,...
        'BoxColor','k',...
        'ComponentColor',color_palette_tol({'cyan','purple','yellow'}),...
        'SaveFolder',resultsFolder,...
        'CloseFigureAfterSaving',true,...
        'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});
end

function plot_relevant_results_truss_moving_node(typeName,systemParameter,initialDesign,nodePositionOptimal,...
    solutionSpaceBox,componentSolutionSpacePlanarTrimming,componentSolutionSpaceCornerBoxRemoval,componentSolutionSpaceHolePunching,...
    componentIndex,nRandomTruss,saveFolder,figureSize,cameraPositionFigureSave)

    resultsFolder = [saveFolder,sprintf('TrussResult%s/',typeName)];
    mkdir(resultsFolder);

    is3dPlot = (size(systemParameter.BaseNodePosition,2)==3);

    % % initial truss + component solution spaces
    % plot_results_truss_generic_moving_node(systemParameter,...
    %     'NodePositionInitial',initialDesign,...
    %     sprintf('PlanarTrimmingSolutionSpace%s',typeName),componentSolutionSpacePlanarTrimming,...
    %     sprintf('CornerBoxRemovalSolutionSpace%s',typeName),componentSolutionSpaceCornerBoxRemoval);
    % save_3d_rotating_video_gif(is3dPlot,gcf,[resultsFolder,'InitialTrussComponent'],cameraPositionFigureSave);
    % save_print_figure(gcf,[resultsFolder,'InitialTrussComponent'],'PrintFormat',{'png','pdf'},'Size',figureSize);
    % 
    % % initial truss + optimized truss + component solution spaces
    % plot_results_truss_generic_moving_node(systemParameter,...
    %     'NodePositionInitial',initialDesign,...
    %     sprintf('NodePositionOptimal%s',typeName),nodePositionOptimal,...
    %     sprintf('PlanarTrimmingSolutionSpace%s',typeName),componentSolutionSpacePlanarTrimming,...
    %     sprintf('CornerBoxRemovalSolutionSpace%s',typeName),componentSolutionSpaceCornerBoxRemoval);
    % save_3d_rotating_video_gif(is3dPlot,gcf,[resultsFolder,'InitialOptimizedTrussComponent'],cameraPositionFigureSave);
    % save_print_figure(gcf,[resultsFolder,'InitialOptimizedTrussComponent'],'PrintFormat',{'png','pdf'},'Size',figureSize);
    % 
    % % initial truss + optimized truss + box solution space + component solution spaces
    % plot_results_truss_generic_moving_node(systemParameter,...
    %     'NodePositionInitial',initialDesign,...
    %     sprintf('NodePositionOptimal%s',typeName),nodePositionOptimal,...
    %     sprintf('BoxSolutionSpace%s',typeName),solutionSpaceBox,...
    %     sprintf('PlanarTrimmingSolutionSpace%s',typeName),componentSolutionSpacePlanarTrimming,...
    %     sprintf('CornerBoxRemovalSolutionSpace%s',typeName),componentSolutionSpaceCornerBoxRemoval);
    % save_3d_rotating_video_gif(is3dPlot,gcf,[resultsFolder,'InitialOptimizedTrussBoxComponent'],cameraPositionFigureSave);
    % save_print_figure(gcf,[resultsFolder,'InitialOptimizedTrussBoxComponent'],'PrintFormat',{'png','pdf'},'Size',figureSize);
    
    % initial truss + box solution space + component solution spaces
    plot_results_truss_generic_moving_node(systemParameter,...
        'NodePositionInitial',initialDesign,...
        sprintf('BoxSolutionSpace%s',typeName),solutionSpaceBox,...
        sprintf('PlanarTrimmingSolutionSpace%s',typeName),componentSolutionSpacePlanarTrimming,...
        sprintf('CornerBoxRemovalSolutionSpace%s',typeName),componentSolutionSpaceCornerBoxRemoval,...
        sprintf('HolePunchingSolutionSpace%s',typeName),componentSolutionSpaceHolePunching);
    save_3d_rotating_video_gif(is3dPlot,gcf,[resultsFolder,'InitialTrussBoxComponent'],cameraPositionFigureSave);
    save_print_figure(gcf,[resultsFolder,'InitialTrussBoxComponent'],'PrintFormat',{'png','pdf'},'Size',figureSize);

    % initial truss + sample trusses + component solution space
    randomCandidateSpace = componentSolutionSpacePlanarTrimming;
    if(isempty(randomCandidateSpace))
        randomCandidateSpace = componentSolutionSpaceCornerBoxRemoval;
    end
    if(isempty(randomCandidateSpace))
        randomCandidateSpace = componentSolutionSpaceHolePunching;
    end
    nDimension = size(systemParameter.BaseNodePosition,2);
    isDesignVariable = isnan(systemParameter.BaseNodePosition);
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
    randomTrussMovingNode = candidate_space_sampling_individual_feasible(randomCandidateSpace,componentIndex,nRandomTruss);
    for i=1:nRandomTruss
        nodePositionRandom = systemParameter.BaseNodePosition;
        nodePositionRandom(isDesignVariable) = column_vector_to_row_major_matrix(randomTrussMovingNode(i,:)',nDimension);
        nodeDisplacementRandom = ...
            truss_analysis(...
                nodePositionInitial,...
                systemParameter.FixedDegreesOfFreedom,...
                systemParameter.NodeForce,...
                systemParameter.NodeElement,...
                systemParameter.ElementCrossSectionArea,...
                systemParameter.ElementYoungsModulus);
        currentName = sprintf('InitialRandomTrussComponent%0*d',get_number_digits_integer(nRandomTruss),i);

        plot_results_truss_generic_moving_node(systemParameter,...
            'NodePositionInitial',initialDesign,...
            'DeformationInitial',nodeDisplacementInitial,...
            sprintf('NodePositionRandom%s',typeName),randomTrussMovingNode(i,:),...
            sprintf('DeformationRandom%s',typeName),nodeDisplacementRandom,...
            sprintf('PlanarTrimmingSolutionSpace%s',typeName),randomCandidateSpace);
        save_3d_rotating_video_gif(is3dPlot,gcf,[resultsFolder,currentName],cameraPositionFigureSave);
        save_print_figure(gcf,[resultsFolder,currentName],'PrintFormat',{'png','pdf'},'Size',figureSize);
    end
end

function figureHandle = plot_results_truss_generic_moving_node(systemParameter,varargin)
    parser = inputParser;
    % (1) initial truss
    parser.addParameter('NodePositionInitial',[]);
    parser.addParameter('DeformationInitial',[]);
    parser.addParameter('InitialTrussOptions',{});
    % (2) optimal truss - displacement
    parser.addParameter('NodePositionOptimalDisplacement',[]);
    parser.addParameter('DeformationOptimalDisplacement',[]);
    parser.addParameter('OptimalTrussDisplacementOptions',{});
    % (3) optimal truss - mass
    parser.addParameter('NodePositionOptimalMass',[]);
    parser.addParameter('DeformationOptimalMass',[]);
    parser.addParameter('OptimalTrussMassOptions',{});
    % (4) optimal truss - displacement and mass
    parser.addParameter('NodePositionOptimalDisplacementAndMass',[]);
    parser.addParameter('DeformationOptimalDisplacementAndMass',[]);
    parser.addParameter('OptimalTrussDisplacementAndMassOptions',{});
    % (5) box - displacement 
    parser.addParameter('BoxSolutionSpaceDisplacement',[]);
    parser.addParameter('BoxSolutionSpaceDisplacementOptions',{});
    % (6) box - mass 
    parser.addParameter('BoxSolutionSpaceMass',[]);
    parser.addParameter('BoxSolutionSpaceMassOptions',{});
    % (7) box - displacement and mass
    parser.addParameter('BoxSolutionSpaceDisplacementAndMass',[]);
    parser.addParameter('BoxSolutionSpaceDisplacementAndMassOptions',{});
    % (8) planar trimming - displacement
    parser.addParameter('PlanarTrimmingSolutionSpaceDisplacement',[]);
    parser.addParameter('PlanarTrimmingSolutionSpaceDisplacementOptions',{});
    % (9) planar trimming - mass
    parser.addParameter('PlanarTrimmingSolutionSpaceMass',[]);
    parser.addParameter('PlanarTrimmingSolutionSpaceMassOptions',{});
    % (10) planar trimming - displacement and mass
    parser.addParameter('PlanarTrimmingSolutionSpaceDisplacementAndMass',[]);
    parser.addParameter('PlanarTrimmingSolutionSpaceDisplacementAndMassOptions',{});
    % (11) corner box removal - displacement
    parser.addParameter('CornerBoxRemovalSolutionSpaceDisplacement',[]);
    parser.addParameter('CornerBoxRemovalSolutionSpaceDisplacementOptions',{});
    % (12) corner box removal - mass
    parser.addParameter('CornerBoxRemovalSolutionSpaceMass',[]);
    parser.addParameter('CornerBoxRemovalSolutionSpaceMassOptions',{});
    % (13) corner box removal - displacement and mass
    parser.addParameter('CornerBoxRemovalSolutionSpaceDisplacementAndMass',[]);
    parser.addParameter('CornerBoxRemovalSolutionSpaceDisplacementAndMassOptions',{});
    % (14) hole punching - displacement
    parser.addParameter('HolePunchingSolutionSpaceDisplacement',[]);
    parser.addParameter('HolePunchingSolutionSpaceDisplacementOptions',{});
    % (15) hole punching - mass
    parser.addParameter('HolePunchingSolutionSpaceMass',[]);
    parser.addParameter('HolePunchingSolutionSpaceMassOptions',{});
    % (16) hole punching - displacement and mass
    parser.addParameter('HolePunchingSolutionSpaceDisplacementAndMass',[]);
    parser.addParameter('HolePunchingSolutionSpaceDisplacementAndMassOptions',{});
    % (17) random trusses - displacement
    parser.addParameter('NodePositionRandomDisplacement',[]);
    parser.addParameter('DeformationRandomDisplacement',[]);
    parser.addParameter('RandomTrussDisplacementOptions',{});
    % (18) random trusses - mass
    parser.addParameter('NodePositionRandomMass',[]);
    parser.addParameter('DeformationRandomMass',[]);
    parser.addParameter('RandomTrussMassOptions',{});
    % (19) random trusses - displacement and mass
    parser.addParameter('NodePositionRandomDisplacementAndMass',[]);
    parser.addParameter('DeformationRandomDisplacementAndMass',[]);
    parser.addParameter('RandomTrussDisplacementAndMassOptions',{});
    % (20) wall
    parser.addParameter('IncludeWall',true);
    parser.addParameter('WallOptions',{});
    % (21) applied force
    parser.addParameter('IncludeAppliedForce',true);
    parser.addParameter('AppliedForceOptions',{});
    % (22) axes
    parser.addParameter('IncludeAxesInformation',false);
    % (23) legend
    parser.addParameter('IncludeLegend',false);
    parser.addParameter('LegendOptions',{});
    
    parser.parse(varargin{:});
    options = parser.Results;

    isDesignVariable = isnan(systemParameter.BaseNodePosition);
    nDimension = size(systemParameter.BaseNodePosition,2);
    nComponent = sum(isDesignVariable(:,1),1);
    isTrussTip = (systemParameter.NodeForce~=0);
    is3dPlot = (size(systemParameter.BaseNodePosition,2)==3);

    %% common options
    % (1) initial truss
    defaultCommonInitialTrussOptions = {'TrussPlotOptions',{'Color',color_palette_tol('grey')},'MaximumLinewidth',3.0,'DisplacementScaleFactor',25};
    % (2) optimal truss - displacement
    defaultCommonOptimalTrussDisplacementOptions = {'TrussPlotOptions',{'Color',color_palette_tol('blue')},'MaximumLinewidth',3.0,'DisplacementScaleFactor',25};
    % (3) optimal truss - mass
    defaultCommonOptimalTrussMassOptions = {'TrussPlotOptions',{'Color',color_palette_tol('blue')},'MaximumLinewidth',3.0,'DisplacementScaleFactor',25};
    % (4) optimal truss - displacement and mass
    defaultCommonOptimalTrussDisplacementAndMassOptions = {'TrussPlotOptions',{'Color',color_palette_tol('blue')},'MaximumLinewidth',3.0,'DisplacementScaleFactor',25};
    % (5) box - displacement 
    defaultCommonBoxSolutionDisplacementOptions = {'EdgeColor','k','FaceAlpha',0.0};
    % (6) box - mass 
    defaultCommonBoxSolutionMassOptions = {'EdgeColor','k','FaceAlpha',0.0};
    % (7) box - displacement and mass
    defaultCommonBoxSolutionDisplacementAndMassOptions = {'EdgeColor','k','FaceAlpha',0.0};
    % (8) planar trimming - displacement
    defaultCommonPlanarTrimmingDisplacementOptions = {'FaceColor',color_palette_tol('cyan'),'EdgeColor',color_palette_tol('cyan')};
    % (9) planar trimming - mass
    defaultCommonPlanarTrimmingMassOptions = {'FaceColor',color_palette_tol('cyan'),'EdgeColor',color_palette_tol('cyan')};
    % (10) planar trimming - displacement and mass
    defaultCommonPlanarTrimmingDisplacementAndMassOptions = {'FaceColor',color_palette_tol('cyan'),'EdgeColor',color_palette_tol('cyan')};
    % (11) corner box removal - displacement
    defaultCommonCornerBoxRemovalDisplacementOptions = {'FaceColor',color_palette_tol('purple'),'EdgeColor',color_palette_tol('purple')};
    % (12) corner box removal - mass
    defaultCommonCornerBoxRemovalMassOptions = {'FaceColor',color_palette_tol('purple'),'EdgeColor',color_palette_tol('purple')};
    % (13) corner box removal - displacement and mass
    defaultCommonCornerBoxRemovalDisplacementAndMassOptions = {'FaceColor',color_palette_tol('purple'),'EdgeColor',color_palette_tol('purple')};
    % (14) hole punching - displacement
    defaultCommonHolePunchingDisplacementOptions = {'FaceColor',color_palette_tol('yellow'),'EdgeColor',color_palette_tol('yellow')};
    % (15) hole punching - mass
    defaultCommonHolePunchingMassOptions = {'FaceColor',color_palette_tol('yellow'),'EdgeColor',color_palette_tol('yellow')};
    % (16) hole punching - displacement and mass
    defaultCommonHolePunchingDisplacementAndMassOptions = {'FaceColor',color_palette_tol('yellow'),'EdgeColor',color_palette_tol('yellow')};
    % (17) random trusses - displacement
    defaultCommonRandomTrussDisplacementOptions = {'MaximumLinewidth',3.0,'DisplacementScaleFactor',25};
    % (18) random trusses - mass
    defaultCommonRandomTrussMassOptions = {'MaximumLinewidth',3.0,'DisplacementScaleFactor',25};
    % (19) random trusses - displacement and mass
    defaultCommonRandomTrussDisplacementAndMassOptions = {'MaximumLinewidth',3.0,'DisplacementScaleFactor',25};
    % (20) wall
    defaultCommonWallOptions = {};
    % (21) applied force
    defaultCommonAppliedForceOptions = {'Color',color_palette_tol('red'),'LineWidth',3.0};
    % (22) axes
    % (23) legend
    defaultCommonLegendOptions = {'location','west'};

    %% dimension-specific options
    if(is3dPlot)
        plotDesignBox = @plot_design_box_3d;
        [wallY,wallZ] = meshgrid([-0.5,1.5]);
        wallX = zeros(size(wallY));

        % (1) initial truss
        defaultDimensionInitialTrussOptions = {};
        % (2) optimal truss - displacement
        defaultDimensionOptimalTrussDisplacementOptions = {};
        % (3) optimal truss - mass
        defaultDimensionOptimalTrussMassOptions = {};
        % (4) optimal truss - displacement and mass
        defaultDimensionOptimalTrussDisplacementAndMassOptions = {};
        % (5) box - displacement 
        defaultDimensionBoxSolutionDisplacementOptions = {};
        % (6) box - mass 
        defaultDimensionBoxSolutionMassOptions = {};
        % (7) box - displacement and mass
        defaultDimensionBoxSolutionDisplacementAndMassOptions = {};
        % (8) planar trimming - displacement
        defaultDimensionPlanarTrimmingDisplacementOptions = {'FaceAlpha',0.2};
        % (9) planar trimming - mass
        defaultDimensionPlanarTrimmingMassOptions = {'FaceAlpha',0.2};
        % (10) planar trimming - displacement and mass
        defaultDimensionPlanarTrimmingDisplacementAndMassOptions = {'FaceAlpha',0.2};
        % (11) corner box removal - displacement
        defaultDimensionCornerBoxRemovalDisplacementOptions = {'FaceAlpha',0.2,'EdgeColor','none'};
        % (12) corner box removal - mass
        defaultDimensionCornerBoxRemovalMassOptions = {'FaceAlpha',0.2,'EdgeColor','none'};
        % (13) corner box removal - displacement and mass
        defaultDimensionCornerBoxRemovalDisplacementAndMassOptions = {'FaceAlpha',0.2,'EdgeColor','none'};
        % (14) hole punching - displacement
        defaultDimensionHolePunchingDisplacementOptions = {'FaceAlpha',0.2};
        % (15) hole punching - mass
        defaultDimensionHolePunchingMassOptions = {'FaceAlpha',0.2};
        % (16) hole punching - displacement and mass
        defaultDimensionHolePunchingDisplacementAndMassOptions = {'FaceAlpha',0.2};
        % (17) random trusses - displacement
        defaultDimensionRandomTrussDisplacementOptions = {};
        % (18) random trusses - mass
        defaultDimensionRandomTrussMassOptions = {};
        % (19) random trusses - displacement and mass
        defaultDimensionRandomTrussDisplacementAndMassOptions = {};
        % (20) wall
        defaultDimensionWallOptions = {'FaceColor','w','EdgeColor','k','FaceAlpha',0.7};
        % (21) applied force
        defaultDimensionAppliedForceOptions = {};
        % (22) axes
        % (23) legend
        defaultDimensionLegendOptions = {};
    else
        plotDesignBox = @plot_design_box_2d;
        wallX = [0 0];
        wallY = [-0.1 1.1];

        % (1) initial truss
        defaultDimensionInitialTrussOptions = {};
        % (2) optimal truss - displacement
        defaultDimensionOptimalTrussDisplacementOptions = {};
        % (3) optimal truss - mass
        defaultDimensionOptimalTrussMassOptions = {};
        % (4) optimal truss - displacement and mass
        defaultDimensionOptimalTrussDisplacementAndMassOptions = {};
        % (5) box - displacement 
        defaultDimensionBoxSolutionDisplacementOptions = {'Linewidth',2.0};
        % (6) box - mass 
        defaultDimensionBoxSolutionMassOptions = {'Linewidth',2.0};
        % (7) box - displacement and mass
        defaultDimensionBoxSolutionDisplacementAndMassOptions = {'Linewidth',2.0};
        % (8) planar trimming - displacement
        defaultDimensionPlanarTrimmingDisplacementOptions = {'FaceAlpha',0.0,'Linewidth',2.0};
        % (9) planar trimming - mass
        defaultDimensionPlanarTrimmingMassOptions = {'FaceAlpha',0.0,'Linewidth',2.0};
        % (10) planar trimming - displacement and mass
        defaultDimensionPlanarTrimmingDisplacementAndMassOptions = {'FaceAlpha',0.0,'Linewidth',2.0};
        % (11) corner box removal - displacement
        defaultDimensionCornerBoxRemovalDisplacementOptions = {'FaceAlpha',0.0,'Linewidth',2.0};
        % (12) corner box removal - mass
        defaultDimensionCornerBoxRemovalMassOptions = {'FaceAlpha',0.0,'Linewidth',2.0};
        % (13) corner box removal - displacement and mass
        defaultDimensionCornerBoxRemovalDisplacementAndMassOptions = {'FaceAlpha',0.0,'Linewidth',2.0};
        % (14) hole punching - displacement
        defaultDimensionHolePunchingDisplacementOptions = {'FaceAlpha',0.0,'Linewidth',2.0};
        % (15) hole punching - mass
        defaultDimensionHolePunchingMassOptions = {'FaceAlpha',0.0,'Linewidth',2.0};
        % (16) hole punching - displacement and mass
        defaultDimensionHolePunchingDisplacementAndMassOptions = {'FaceAlpha',0.0,'Linewidth',2.0};
        % (17) random trusses - displacement
        defaultDimensionRandomTrussDisplacementOptions = {};
        % (18) random trusses - mass
        defaultDimensionRandomTrussMassOptions = {};
        % (19) random trusses - displacement and mass
        defaultDimensionRandomTrussDisplacementAndMassOptions = {};
        % (20) wall
        defaultDimensionWallOptions = {'Linewidth',8.0,'Color','k'};
        % (21) applied force
        defaultDimensionAppliedForceOptions = {};
        % (22) axes
        % (23) legend 
        defaultDimensionLegendOptions = {};
    end

    %% merge options
    % (1) initial truss
    [~,initialTrussOptions] = merge_name_value_pair_argument(...
        defaultCommonInitialTrussOptions,...
        defaultDimensionInitialTrussOptions,...
        options.InitialTrussOptions);
    % (2) optimal truss - displacement
    [~,optimalTrussDisplacementOptions] = merge_name_value_pair_argument(...
        defaultCommonOptimalTrussDisplacementOptions,...
        defaultDimensionOptimalTrussDisplacementOptions,...
        options.OptimalTrussDisplacementOptions);
    % (3) optimal truss - mass
    [~,optimalTrussMassOptions] = merge_name_value_pair_argument(...
        defaultCommonOptimalTrussMassOptions,...
        defaultDimensionOptimalTrussMassOptions,...
        options.OptimalTrussMassOptions);
    % (4) optimal truss - displacement and mass
    [~,optimalTrussDisplacementAndMassOptions] = merge_name_value_pair_argument(...
        defaultCommonOptimalTrussDisplacementAndMassOptions,...
        defaultDimensionOptimalTrussDisplacementAndMassOptions,...
        options.OptimalTrussDisplacementAndMassOptions);
    % (5) box - displacement 
    [~,boxDisplacementOptions] = merge_name_value_pair_argument(...
        defaultCommonBoxSolutionDisplacementOptions,...
        defaultDimensionBoxSolutionDisplacementOptions,...
        options.BoxSolutionSpaceDisplacementOptions);
    % (6) box - mass 
    [~,boxMassOptions] = merge_name_value_pair_argument(...
        defaultCommonBoxSolutionMassOptions,...
        defaultDimensionBoxSolutionMassOptions,...
        options.BoxSolutionSpaceMassOptions);
    % (7) box - displacement and mass
    [~,boxDisplacementAndMassOptions] = merge_name_value_pair_argument(...
        defaultCommonBoxSolutionDisplacementAndMassOptions,...
        defaultDimensionBoxSolutionDisplacementAndMassOptions,...
        options.BoxSolutionSpaceDisplacementAndMassOptions);
    % (8) planar trimming - displacement
    [~,planarTrimmingDisplacementOptions] = merge_name_value_pair_argument(...
        defaultCommonPlanarTrimmingDisplacementOptions,...
        defaultDimensionPlanarTrimmingDisplacementOptions,...
        options.PlanarTrimmingSolutionSpaceDisplacementOptions);
    % (9) planar trimming - mass
    [~,planarTrimmingMassOptions] = merge_name_value_pair_argument(...
        defaultCommonPlanarTrimmingMassOptions,...
        defaultDimensionPlanarTrimmingMassOptions,...
        options.PlanarTrimmingSolutionSpaceMassOptions);
    % (10) planar trimming - displacement and mass
    [~,planarTrimmingDisplacementAndMassOptions] = merge_name_value_pair_argument(...
        defaultCommonPlanarTrimmingDisplacementAndMassOptions,...
        defaultDimensionPlanarTrimmingDisplacementAndMassOptions,...
        options.PlanarTrimmingSolutionSpaceDisplacementAndMassOptions);
    % (11) corner box removal - displacement
    [~,cornerBoxRemovalDisplacementOptions] = merge_name_value_pair_argument(...
        defaultCommonCornerBoxRemovalDisplacementOptions,...
        defaultDimensionCornerBoxRemovalDisplacementOptions,...
        options.CornerBoxRemovalSolutionSpaceDisplacementOptions);
    % (12) corner box removal - mass
    [~,cornerBoxRemovalMassOptions] = merge_name_value_pair_argument(...
        defaultCommonCornerBoxRemovalMassOptions,...
        defaultDimensionCornerBoxRemovalMassOptions,...
        options.CornerBoxRemovalSolutionSpaceMassOptions);
    % (13) corner box removal - displacement and mass
    [~,cornerBoxRemovalDisplacementAndMassOptions] = merge_name_value_pair_argument(...
        defaultCommonCornerBoxRemovalDisplacementAndMassOptions,...
        defaultDimensionCornerBoxRemovalDisplacementAndMassOptions,...
        options.CornerBoxRemovalSolutionSpaceDisplacementAndMassOptions);
    % (14) hole punching - displacement
    [~,holePunchingDisplacementOptions] = merge_name_value_pair_argument(...
        defaultCommonHolePunchingDisplacementOptions,...
        defaultDimensionHolePunchingDisplacementOptions,...
        options.HolePunchingSolutionSpaceDisplacementOptions);
    % (15) hole punching - mass
    [~,holePunchingMassOptions] = merge_name_value_pair_argument(...
        defaultCommonHolePunchingMassOptions,...
        defaultDimensionHolePunchingMassOptions,...
        options.HolePunchingSolutionSpaceMassOptions);
    % (16) hole punching - displacement and mass
    [~,holePunchingDisplacementAndMassOptions] = merge_name_value_pair_argument(...
        defaultCommonHolePunchingDisplacementAndMassOptions,...
        defaultDimensionHolePunchingDisplacementAndMassOptions,...
        options.HolePunchingSolutionSpaceDisplacementAndMassOptions);
    % (17) random trusses - displacement
    [~,randomTrussDisplacementOptions] = merge_name_value_pair_argument(...
        defaultCommonRandomTrussDisplacementOptions,...
        defaultDimensionRandomTrussDisplacementOptions,...
        options.RandomTrussDisplacementOptions);
    % (18) random trusses - mass
    [~,randomTrussMassOptions] = merge_name_value_pair_argument(...
        defaultCommonRandomTrussMassOptions,...
        defaultDimensionRandomTrussMassOptions,...
        options.RandomTrussMassOptions);
    % (19) random trusses - displacement and mass
    [~,randomTrussDisplacementAndMassOptions] = merge_name_value_pair_argument(...
        defaultCommonRandomTrussDisplacementAndMassOptions,...
        defaultDimensionRandomTrussDisplacementAndMassOptions,...
        options.RandomTrussDisplacementAndMassOptions);
    % (20) wall
    [~,wallOptions] = merge_name_value_pair_argument(...
        defaultCommonWallOptions,...
        defaultDimensionWallOptions,...
        options.WallOptions);
    % (21) applied force
    [~,appliedForceOptions] = merge_name_value_pair_argument(...
        defaultCommonAppliedForceOptions,...
        defaultDimensionAppliedForceOptions,...
        options.AppliedForceOptions);
    % (22) axes 
    % (23) legend
    [~,legendOptions] = merge_name_value_pair_argument(...
        defaultCommonLegendOptions,...
        defaultDimensionLegendOptions,...
        options.LegendOptions);


    %% plot each element (where applicable)
    figureHandle = figure;
    hold all;
    
    % (1) initial truss
    handleInitial = [];
    handleInitialDeformed = [];
    if(~isempty(options.NodePositionInitial))
        nodePositionInitial = systemParameter.BaseNodePosition;
        nodePositionInitial(isDesignVariable) = column_vector_to_row_major_matrix(options.NodePositionInitial',nDimension);

        [handleInitial,handleInitialDeformed] = plot_truss_deformation(gcf,nodePositionInitial,systemParameter.NodeElement,options.DeformationInitial,initialTrussOptions{:});
    end

    % (2) optimal truss - displacement
    handleOptimalDisplacement = [];
    handleOptimizedDisplacementDeformed = [];
    if(~isempty(options.NodePositionOptimalDisplacement))
        nodePositionOptimized = systemParameter.BaseNodePosition;
        nodePositionOptimized(isDesignVariable) = column_vector_to_row_major_matrix(options.NodePositionOptimalDisplacement',nDimension);

        [handleOptimalDisplacement,handleOptimizedDisplacementDeformed] = plot_truss_deformation(gcf,nodePositionOptimized,systemParameter.NodeElement,options.DeformationOptimalDisplacement,optimalTrussDisplacementOptions{:});
    end

    % (3) optimal truss - mass
    handleOptimalMass = [];
    handleOptimizedMassDeformed = [];
    if(~isempty(options.NodePositionOptimalMass))
        nodePositionOptimized = systemParameter.BaseNodePosition;
        nodePositionOptimized(isDesignVariable) = column_vector_to_row_major_matrix(options.NodePositionOptimalMass',nDimension);

        [handleOptimalMass,handleOptimizedMassDeformed] = plot_truss_deformation(gcf,nodePositionOptimized,systemParameter.NodeElement,options.DeformationOptimalMass,optimalTrussMassOptions{:});
    end

    % (4) optimal truss - displacement and mass
    handleOptimalDisplacementAndMass = [];
    handleOptimizedDisplacementAndMassDeformed = [];
    if(~isempty(options.NodePositionOptimalDisplacementAndMass))
        nodePositionOptimized = systemParameter.BaseNodePosition;
        nodePositionOptimized(isDesignVariable) = column_vector_to_row_major_matrix(options.NodePositionOptimalDisplacementAndMass',nDimension);

        [handleOptimalDisplacementAndMass,handleOptimizedDisplacementAndMassDeformed] = plot_truss_deformation(gcf,nodePositionOptimized,systemParameter.NodeElement,options.DeformationOptimalDisplacementAndMass,optimalTrussDisplacementAndMassOptions{:});
    end

    % (5) box - displacement 
    handleToleranceNodeDisplacementBox = [];
    if(~isempty(options.BoxSolutionSpaceDisplacement))
        for i=1:nComponent-1
            currentIndex = 1 + nDimension*(i-1) + [0:(nDimension-1)];
            plotDesignBox(gcf,options.BoxSolutionSpaceDisplacement(:,currentIndex),boxDisplacementOptions{:});

        end
        currentIndex = 1 + nDimension*(nComponent-1) + [0:(nDimension-1)];
        handleToleranceNodeDisplacementBox = plotDesignBox(gcf,options.BoxSolutionSpaceDisplacement(:,currentIndex),boxDisplacementOptions{:});
    end

    % (6) box - mass 
    handleToleranceNodeMassBox = [];
    if(~isempty(options.BoxSolutionSpaceMass))
        for i=1:nComponent-1
            currentIndex = 1 + nDimension*(i-1) + [0:(nDimension-1)];
            plotDesignBox(gcf,options.BoxSolutionSpaceMass(:,currentIndex),boxMassOptions{:});
        end
        currentIndex = 1 + nDimension*(nComponent-1) + [0:(nDimension-1)];
        handleToleranceNodeMassBox = plotDesignBox(gcf,options.BoxSolutionSpaceMass(:,currentIndex),boxMassOptions{:});
    end

    % (7) box - displacement and mass
    handleToleranceNodeDisplacementAndMassBox = [];
    if(~isempty(options.BoxSolutionSpaceDisplacementAndMass))
        for i=1:nComponent-1
            currentIndex = 1 + nDimension*(i-1) + [0:(nDimension-1)];
            plotDesignBox(gcf,options.BoxSolutionSpaceDisplacementAndMass(:,currentIndex),boxDisplacementAndMassOptions{:});
        end
        currentIndex = 1 + nDimension*(nComponent-1) + [0:(nDimension-1)];
        handleToleranceNodeDisplacementAndMassBox = plotDesignBox(gcf,options.BoxSolutionSpaceDisplacementAndMass(:,currentIndex),boxDisplacementAndMassOptions{:});
    end

    % (8) planar trimming - displacement
    handleToleranceNodeDisplacementComponentPlanarTrimming = [];
    if(~isempty(options.PlanarTrimmingSolutionSpaceDisplacement))
        for i=1:nComponent-1
            options.PlanarTrimmingSolutionSpaceDisplacement(i).plot_candidate_space(gcf,planarTrimmingDisplacementOptions{:});
        end
        handleToleranceNodeDisplacementComponentPlanarTrimming = options.PlanarTrimmingSolutionSpaceDisplacement(nComponent).plot_candidate_space(gcf,planarTrimmingDisplacementOptions{:});
    end

    % (9) planar trimming - mass
    handleToleranceNodeMassComponentPlanarTrimming = [];
    if(~isempty(options.PlanarTrimmingSolutionSpaceMass))
        for i=1:nComponent-1
            options.PlanarTrimmingSolutionSpaceMass(i).plot_candidate_space(gcf,planarTrimmingMassOptions{:});
        end
        handleToleranceNodeMassComponentPlanarTrimming = options.PlanarTrimmingSolutionSpaceMass(nComponent).plot_candidate_space(gcf,planarTrimmingMassOptions{:});
    end

    % (10) planar trimming - displacement and mass
    handleToleranceNodeDisplacementAndMassComponentPlanarTrimming = [];
    if(~isempty(options.PlanarTrimmingSolutionSpaceDisplacementAndMass))
        for i=1:nComponent-1
            options.PlanarTrimmingSolutionSpaceDisplacementAndMass(i).plot_candidate_space(gcf,planarTrimmingDisplacementAndMassOptions{:});
        end
        handleToleranceNodeDisplacementAndMassComponentPlanarTrimming = options.PlanarTrimmingSolutionSpaceDisplacementAndMass(nComponent).plot_candidate_space(gcf,planarTrimmingDisplacementAndMassOptions{:});
    end

    % (11) corner box removal - displacement
    handleToleranceNodeDisplacementComponentCornerBoxRemoval = [];
    if(~isempty(options.CornerBoxRemovalSolutionSpaceDisplacement))
        for i=1:nComponent-1
            options.CornerBoxRemovalSolutionSpaceDisplacement(i).plot_candidate_space(gcf,cornerBoxRemovalDisplacementOptions{:});
        end
        handleToleranceNodeDisplacementComponentCornerBoxRemoval = options.CornerBoxRemovalSolutionSpaceDisplacement(nComponent).plot_candidate_space(gcf,cornerBoxRemovalDisplacementOptions{:});
    end

    % (12) corner box removal - mass
    handleToleranceNodeMassComponentCornerBoxRemoval = [];
    if(~isempty(options.CornerBoxRemovalSolutionSpaceMass))
        for i=1:nComponent-1
            options.CornerBoxRemovalSolutionSpaceMass(i).plot_candidate_space(gcf,cornerBoxRemovalMassOptions{:});
        end
        handleToleranceNodeMassComponentCornerBoxRemoval = options.CornerBoxRemovalSolutionSpaceMass(nComponent).plot_candidate_space(gcf,cornerBoxRemovalMassOptions{:});
    end

    % (13) corner box removal - displacement and mass
    handleToleranceNodeDisplacementAndMassComponentCornerBoxRemoval = [];
    if(~isempty(options.CornerBoxRemovalSolutionSpaceDisplacementAndMass))
        for i=1:nComponent-1
            options.CornerBoxRemovalSolutionSpaceDisplacementAndMass(i).plot_candidate_space(gcf,cornerBoxRemovalDisplacementAndMassOptions{:});
        end
        handleToleranceNodeDisplacementAndMassComponentCornerBoxRemoval = options.CornerBoxRemovalSolutionSpaceDisplacementAndMass(nComponent).plot_candidate_space(gcf,cornerBoxRemovalDisplacementAndMassOptions{:});
    end

    % (14) hole punching - displacement
    handleToleranceNodeDisplacementComponentHolePunching = [];
    if(~isempty(options.HolePunchingSolutionSpaceDisplacement))
        for i=1:nComponent-1
            options.HolePunchingSolutionSpaceDisplacement(i).plot_candidate_space(gcf,holePunchingDisplacementOptions{:});
        end
        handleToleranceNodeDisplacementComponentHolePunching = options.HolePunchingSolutionSpaceDisplacement(nComponent).plot_candidate_space(gcf,holePunchingDisplacementOptions{:});
    end

    % (15) hole punching - mass
    handleToleranceNodeMassComponentHolePunching = [];
    if(~isempty(options.HolePunchingSolutionSpaceMass))
        for i=1:nComponent-1
            options.HolePunchingSolutionSpaceMass(i).plot_candidate_space(gcf,holePunchingMassOptions{:});
        end
        handleToleranceNodeMassComponentHolePunching = options.HolePunchingSolutionSpaceMass(nComponent).plot_candidate_space(gcf,holePunchingMassOptions{:});
    end

    % (16) hole punching - displacement and mass
    handleToleranceNodeDisplacementAndMassComponentHolePunching = [];
    if(~isempty(options.HolePunchingSolutionSpaceDisplacementAndMass))
        for i=1:nComponent-1
            options.HolePunchingSolutionSpaceDisplacementAndMass(i).plot_candidate_space(gcf,holePunchingDisplacementAndMassOptions{:});
        end
        handleToleranceNodeDisplacementAndMassComponentHolePunching = options.HolePunchingSolutionSpaceDisplacementAndMass(nComponent).plot_candidate_space(gcf,holePunchingDisplacementAndMassOptions{:});
    end

    % (17) random trusses - displacement
    handleRandomDisplacement = [];
    handleRandomDisplacementDeformed = [];
    if(~isempty(options.NodePositionRandomDisplacement))
        nRandomTruss = size(options.NodePositionRandomDisplacement,1);
        trussColor = rand(nRandomTruss,3);
        for i=1:nRandomTruss
            randomNodePosition = systemParameter.BaseNodePosition;
            randomNodePosition(isDesignVariable) = column_vector_to_row_major_matrix(options.NodePositionRandomDisplacement(i,:)',nDimension);

            [handleRandomDisplacement(i),handleRandomDisplacementDeformed(i)] = plot_truss_deformation(gcf,randomNodePosition,systemParameter.NodeElement,options.DeformationRandomDisplacement,'TrussPlotOptions',{'Color',trussColor(i,:)},randomTrussDisplacementOptions{:});
        end
    end

    % (18) random trusses - mass
    handleRandomMass = [];
    handleRandomMassDeformed = [];
    if(~isempty(options.NodePositionRandomMass))
        nRandomTruss = size(options.NodePositionRandomMass,1);
        trussColor = rand(nRandomTruss,3);
        for i=1:nRandomTruss
            randomNodePosition = systemParameter.BaseNodePosition;
            randomNodePosition(isDesignVariable) = column_vector_to_row_major_matrix(options.NodePositionRandomMass(i,:)',nDimension);

            [handleRandomMass(i),handleRandomMassDeformed(i)] = plot_truss_deformation(gcf,randomNodePosition,systemParameter.NodeElement,options.DeformationRandomMass,'TrussPlotOptions',{'Color',trussColor(i,:)},randomTrussMassOptions{:});
        end
    end

    % (19) random trusses - displacement and mass
    handleRandomDisplacementAndMass = [];
    handleRandomDisplacementAndMassDeformed = [];
    if(~isempty(options.NodePositionRandomDisplacementAndMass))
        nRandomTruss = size(options.NodePositionRandomDisplacementAndMass,1);
        trussColor = rand(nRandomTruss,3);
        for i=1:nRandomTruss
            randomNodePosition = systemParameter.BaseNodePosition;
            randomNodePosition(isDesignVariable) = column_vector_to_row_major_matrix(options.NodePositionRandomDisplacementAndMass(i,:)',nDimension);

            [handleRandomDisplacementAndMass(i),handleRandomDisplacementAndMassDeformed(i)] = plot_truss_deformation(gcf,randomNodePosition,systemParameter.NodeElement,options.DeformationRandomDisplacementAndMass,'TrussPlotOptions',{'Color',trussColor(i,:)},randomTrussDisplacementAndMassOptions{:});
        end
    end

    % (20) wall
    handleWall = [];
    if(options.IncludeWall)
        if(is3dPlot)
            handleWall = surf(wallX,wallY,wallZ,wallOptions{:});
        else
            handleWall = plot(wallX,wallY,wallOptions{:});
        end
    end

    % (21) applied force
    handleForce = [];
    if(options.IncludeAppliedForce)
        positionTip = systemParameter.BaseNodePosition(any(isTrussTip,2),:);
        if(is3dPlot)
            handleForce = quiver3(positionTip(1),positionTip(2),positionTip(3),0,0,-0.5,appliedForceOptions{:});
        else
            handleForce = quiver(positionTip(1),positionTip(2),0,-0.5,appliedForceOptions{:});
        end
    end

    % (22) axes
    grid('off');
    if(~options.IncludeAxesInformation)
        if(is3dPlot)
            set(gca,'XColor', 'none','YColor','none','ZColor','none');
        else
            set(gca,'XColor', 'none','YColor','none');
        end
    end

    % (23) legend
    if(options.IncludeLegend)
        handleObjectAll = {...
            handleInitial,...
            handleInitialDeformed,...
            handleOptimalDisplacement,...
            handleOptimizedDisplacementDeformed,...
            handleOptimalMass,...
            handleOptimizedMassDeformed,...
            handleOptimalDisplacementAndMass,...
            handleOptimizedDisplacementAndMassDeformed,...
            handleWall,...
            handleForce,...
            handleToleranceNodeDisplacementBox,...
            handleToleranceNodeDisplacementComponentPlanarTrimming,...
            handleToleranceNodeDisplacementComponentCornerBoxRemoval,...
            handleToleranceNodeDisplacementComponentHolePunching,...
            handleToleranceNodeMassBox,...
            handleToleranceNodeMassComponentPlanarTrimming,...
            handleToleranceNodeMassComponentCornerBoxRemoval,...
            handleToleranceNodeMassComponentHolePunching,...
            handleToleranceNodeDisplacementAndMassBox,...
            handleToleranceNodeDisplacementAndMassComponentPlanarTrimming,...
            handleToleranceNodeDisplacementAndMassComponentCornerBoxRemoval,...
            handleToleranceNodeDisplacementAndMassComponentHolePunching};

        legendTextAll = {...
            'Initial Truss',...
            'Initial Truss (Deformed)',...
            'Optimized Truss - Displacement',...
            'Optimized Truss (Deformed) - Displacement',...
            'Optimized Truss - Mass',...
            'Optimized Truss (Deformed) - Mass',...
            'Optimized Truss - Displacement and Mass',...
            'Optimized Truss (Deformed) - Displacement and Mass',...
            'Wall',...
            'Applied Force',...
            'Tolerance Region for the Node (Box) - Displacement',...
            'Tolerance Region for the Node (Planar Trimming) - Displacement',...
            'Tolerance Region for the Node (Corner Box Removal) - Displacement',...
            'Tolerance Region for the Node (Hole Punching) - Displacement',...
            'Tolerance Region for the Node (Box) - Mass',...
            'Tolerance Region for the Node (Planar Trimming) - Mass',...
            'Tolerance Region for the Node (Corner Box Removal) - Mass',...
            'Tolerance Region for the Node (Hole Punching) - Mass',...
            'Tolerance Region for the Node (Box) - Displacement + Mass',...
            'Tolerance Region for the Node (Planar Trimming) - Displacement + Mass',...
            'Tolerance Region for the Node (Corner Box Removal) - Displacement + Mass',...
            'Tolerance Region for the Node (Hole Punching) - Displacement + Mass'};
        
        handleObject = [handleObjectAll{:}];
        legentText = {legendTextAll{~cellfun(@isempty,handleObjectAll)}};

        legend(handleObject,legentText,legendOptions{:});
    end

    if(is3dPlot)
        axis('tight','equal','vis3d');
        camproj('perspective');
        cameratoolbar; % better adjust angle/perspective
    end

    if(nargout<1)
        clear figureHandle;
    end
end

function save_3d_rotating_video_gif(is3dPlot,figureHandle,filename,cameraPositionFigureSave,varargin)
    if(~is3dPlot)
        return;
    end

    [defaultAzimuth,defaultElevation] = view(3);
    % rotate so truss can be better seen at the start
    rotatedAzimuth = defaultAzimuth + 55;
    view(rotatedAzimuth,defaultElevation);
    rotating_video(figureHandle,filename,varargin{:});
    rotating_gif(figureHandle,filename,varargin{:});

    if(~isempty(cameraPositionFigureSave))
        campos(cameraPositionFigureSave);
    end
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

