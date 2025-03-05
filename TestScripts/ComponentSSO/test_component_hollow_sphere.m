%TEST_COMPONENT_HOLLOW_SPHERE Hollow sphere design problem
%   TEST_COMPONENT_HOLLOW_SPHERE uses a hollow sphere geometry to test the 
%   computation of component solution spaces. The problem involves finding
%   points that lie within a hollow spherical shell defined by inner and
%   outer radii. The solution space is computed using both box-shaped and
%   component-based approaches, and the results are compared.
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
% Define sphere parameters
sphereInnerRadius = 2;
sphereOuterRadius = 5;  % radius of the sphere
sphereCenter = [0, 0, 0];  % center coordinates

systemFunction = @distance_to_center;
systemParameter = sphereCenter;  % center coordinates [x,y,z]

% Design space bounds
designSpaceLowerBound = [-6 -6 -6];  % Lower bounds for x,y,z
designSpaceUpperBound = [6 6 6];     % Upper bounds for x,y,z

% Initial design at sphere center
initialDesign = [3 0 0];


%% Perform Solution Space Computation
bottomUpMapping = BottomUpMappingFunction(systemFunction,systemParameter);

performanceLowerLimit = sphereInnerRadius;
performanceUpperLimit = sphereOuterRadius;  % Points within sphere radius are good designs

optionsBox = sso_stochastic_options('box',...
    'UseAdaptiveGrowthRate',true,...
    'NumberSamplesPerIteration',100,...
    'GrowthRate',0.1,...
    'FixIterNumberExploration',true,...
    'FixIterNumberConsolidation',true,...
    'MaxIterExploration',20,...
    'MaxIterConsolidation',20,...
    'TrimmingOperationOptions',{'PassesCriterion','full'},...
    'TrimmingOrderOptions',{'OrderPreference','score'},...
    'LoggingLevel','all');

designEvaluator = DesignEvaluatorBottomUpMapping(bottomUpMapping,performanceLowerLimit,performanceUpperLimit);
[designBoxOptimal,optimizationDataBox] = sso_box_stochastic(designEvaluator,...
    initialDesign,designSpaceLowerBound,designSpaceUpperBound,optionsBox);


%% Component solution spaces setup - planar trimming
options = sso_stochastic_options('component',...
    'UseAdaptiveGrowthRate',true,...
    'TargetAcceptedRatioExploration',0.7,...
    'GrowthRate',0.1,...
    'FixIterNumberExploration',true,...
    'FixIterNumberConsolidation',true,...
    'MaxIterExploration',20,...
    'MaxIterConsolidation',20,...
    'NumberSamplesPerIteration',100,...
    'CandidateSpaceConstructor',@CandidateSpaceConvexHull,...
    'TrimmingMethodFunction',@component_trimming_method_planar_trimming,...
    'UsePaddingSamplesInTrimming',true,...
    'UsePreviousEvaluatedSamplesConsolidation',false,...
    'ApplyLeanness','never',...
    'LoggingLevel','all',...
    'TrimmingOperationOptions',{'PassesCriterion','full'},...
    'TrimmingOrderOptions',{'OrderPreference','score'});

% Define components - split into xy-plane and z coordinate
componentIndex = {[1,2],[3]};
        
[planarTrimmingSolutionSpace,optimizationDataPlanarTrimming] = sso_component_stochastic(designEvaluator,...
    initialDesign,designSpaceLowerBound,designSpaceUpperBound,componentIndex,options);


%% Component solution spaces - corner box removal
options = sso_stochastic_options('component',...
    'UseAdaptiveGrowthRate',true,...
    'TargetAcceptedRatioExploration',0.7,...
    'GrowthRate',0.1,...
    'FixIterNumberExploration',true,...
    'FixIterNumberConsolidation',true,...
    'MaxIterExploration',20,...
    'MaxIterConsolidation',20,...
    'NumberSamplesPerIteration',100,...
    'CandidateSpaceConstructor',@CandidateSpaceCornerBoxRemoval,...
    'TrimmingMethodFunction',@component_trimming_method_corner_box_removal,...
    'UsePaddingSamplesInTrimming',true,...
    'UsePreviousEvaluatedSamplesConsolidation',false,...
    'ApplyLeanness','never',...
    'LoggingLevel','all',...
    'TrimmingOperationOptions',{'PassesCriterion','full'},...
    'TrimmingOrderOptions',{'OrderPreference','score'});

% Define components - split into xy-plane and z coordinate
componentIndex = {[1,2],[3]};
        
[cornerBoxRemovalSolutionSpace,optimizationDataCornerBoxRemoval] = sso_component_stochastic(designEvaluator,...
    initialDesign,designSpaceLowerBound,designSpaceUpperBound,componentIndex,options);


%% visualization
% Create sphere surface for visualization
[sphereX, sphereY, sphereZ] = sphere(50);
outerSphereX = sphereX * sphereOuterRadius;
outerSphereY = sphereY * sphereOuterRadius;
outerSphereZ = sphereZ * sphereOuterRadius;
innerSphereX = sphereX * sphereInnerRadius;
innerSphereY = sphereY * sphereInnerRadius;
innerSphereZ = sphereZ * sphereInnerRadius;

% Plot inner and outer spheres in 3D
figure;
hold on;
surf(outerSphereX, outerSphereY, outerSphereZ, 'FaceColor', 'red', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
surf(innerSphereX, innerSphereY, innerSphereZ, 'FaceColor', 'red', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
axis equal;
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
title('Inner and Outer Spheres');
view(3);
save_print_figure(gcf,[saveFolder,'SpheresView']);


% Plot Component 1 (xy-plane)
figure;
hold on;
% Plot box solution space
plot_design_box_2d(gcf,designBoxOptimal(:,componentIndex{1}),'FaceAlpha',0.5,'FaceColor',[0.8 0.8 0.8]);
% Plot component solution spaces
planarTrimmingSolutionSpace(1).plot_candidate_space(gcf,'FaceColor','green','EdgeColor','green','FaceAlpha',0.1);
cornerBoxRemovalSolutionSpace(1).plot_candidate_space(gcf,'FaceColor','blue','EdgeColor','blue','FaceAlpha',0.1);
xlabel('x');
ylabel('y');
axis equal;
grid minor;
title('Component 1 (xy-plane)');
legend({'Box Solution Space', 'Planar Trimming', 'Corner Box Removal'});
save_print_figure(gcf,[saveFolder,'Component1View']);

% Plot Component 2 (z coordinate)
figure;
hold on;
% Plot box solution space
plot_design_box_1d(gcf,designBoxOptimal(:,componentIndex{2}),'Color',[0.8 0.8 0.8],'YValue',-1);
% Plot component solution spaces
planarTrimmingSolutionSpace(2).plot_candidate_space(gcf,'Color','green','YValue',0);
cornerBoxRemovalSolutionSpace(2).plot_candidate_space(gcf,'Color','blue','YValue',1);
axis equal;
grid minor;
title('Component 2 (z coordinate)');
legend({'Box Solution Space', 'Planar Trimming', 'Corner Box Removal'});
save_print_figure(gcf,[saveFolder,'Component2View']);


% Plot 3D view of all solution spaces
figure;
hold on;

% Plot spheres for reference
surf(outerSphereX, outerSphereY, outerSphereZ, 'FaceColor', 'red', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
surf(innerSphereX, innerSphereY, innerSphereZ, 'FaceColor', 'red', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% Plot box solution space
plot_design_box_3d(gcf, designBoxOptimal, 'FaceAlpha', 0.3, 'FaceColor', [0.8 0.8 0.8]);


% Plot component solution spaces by extruding 2D spaces

% Planar Trimming
% Create vertices for each point in xy combined with each z value
xyPointsPlanar = planarTrimmingSolutionSpace(1).ActiveDesign;
zPointsPlanar = [min(planarTrimmingSolutionSpace(2).ActiveDesign) ; max(planarTrimmingSolutionSpace(2).ActiveDesign)];
numXYPoints = size(xyPointsPlanar, 1);
numZPoints = length(zPointsPlanar);
verticesPlanar = zeros(numXYPoints * numZPoints, 3);
idx = 1;
for i = 1:numXYPoints
    for j = 1:numZPoints
        verticesPlanar(idx,:) = [xyPointsPlanar(i,:), zPointsPlanar(j)];
        idx = idx + 1;
    end
end
% Create a boundary around these points
k = boundary(verticesPlanar, 0.5); % Use boundary function to create surface
trisurf(k, verticesPlanar(:,1), verticesPlanar(:,2), verticesPlanar(:,3), ...
    'FaceColor', 'green', 'EdgeColor', 'green', 'FaceAlpha', 0.2);

% Corner Box Removal
% Create vertices for each point in xy combined with each z value
xyPointsCorner = cornerBoxRemovalSolutionSpace(1).ActiveDesign;
zPointsCorner = [min(cornerBoxRemovalSolutionSpace(2).ActiveDesign) ; max(cornerBoxRemovalSolutionSpace(2).ActiveDesign)];
numXYPoints = size(xyPointsCorner, 1);
numZPoints = length(zPointsCorner);
verticesCorner = zeros(numXYPoints * numZPoints, 3);
idx = 1;
for i = 1:numXYPoints
    for j = 1:numZPoints
        verticesCorner(idx,:) = [xyPointsCorner(i,:), zPointsCorner(j)];
        idx = idx + 1;
    end
end
% Create a boundary around these points
k = boundary(verticesCorner, 0.5);
trisurf(k, verticesCorner(:,1), verticesCorner(:,2), verticesCorner(:,3), ...
    'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.2);

axis equal;
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
title('3D View of Solution Spaces');
legend({'Outer Sphere', 'Inner Sphere', 'Box Solution Space', ...
        'Planar Trimming', 'Corner Box Removal'});
view(3);
save_print_figure(gcf,[saveFolder,'SolutionSpaces3DView']);



%% performance metrics
resultsFolder = [saveFolder,sprintf('ResultsBox/')];
mkdir(resultsFolder);
algoDataBox = postprocess_sso_box_stochastic(optimizationDataBox);
plot_sso_box_stochastic_metrics(algoDataBox,...
    'SaveFolder',resultsFolder,...
    'CloseFigureAfterSaving',true,...
    'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});

resultsFolder = [saveFolder,sprintf('ResultsPlanarTrimming/')];
mkdir(resultsFolder);
algoDataPlanarTrimming = postprocess_sso_component_stochastic(optimizationDataPlanarTrimming);
plot_sso_component_stochastic_metrics(algoDataPlanarTrimming,...
    'SaveFolder',resultsFolder,...
    'CloseFigureAfterSaving',true,...
    'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});

resultsFolder = [saveFolder,sprintf('ResultsCornerBoxRemoval/')];
mkdir(resultsFolder);
algoDataCornerBoxRemoval = postprocess_sso_component_stochastic(optimizationDataCornerBoxRemoval);
plot_sso_component_stochastic_metrics(algoDataCornerBoxRemoval,...
    'SaveFolder',resultsFolder,...
    'CloseFigureAfterSaving',true,...
    'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});

resultsFolder = [saveFolder,sprintf('ResultsComparison/')];
mkdir(resultsFolder);
plot_sso_comparison_box_component_stochastic_metrics(...
    {algoDataBox},...
    {algoDataPlanarTrimming,algoDataCornerBoxRemoval},...
    'ComponentLabel',{'Planar Trimming','Corner Box Removal'},...
    'BoxColor','k',...
    'ComponentColor',color_palette_tol({'green','blue'}),...
    'SaveFolder',resultsFolder,...
    'CloseFigureAfterSaving',false,...
    'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

