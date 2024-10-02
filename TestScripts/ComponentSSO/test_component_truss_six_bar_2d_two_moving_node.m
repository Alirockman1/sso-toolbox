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


%% function call
%
systemFunction = @truss_six_bar_2d_two_moving_node;
systemParameter = [10,210e3,7850e-9]; % [mm^2],[MPa],[kg/mm^3]
%                         x1   y1  x2  y2
designSpaceLowerBound = [0.1 -0.5 0.9 0.5];
designSpaceUpperBound = [1.5  0.5 2.0 1.5];
Components = {[1,2]',[3,4]'};
%
performanceLowerLimit = -inf;
performanceUpperLimit = [nan nan repmat(250,1,6)];
%
initialDesign = [1,0,1,1];


%% find optimum
bottomUpMapping = BottomUpMappingFunction(systemFunction,systemParameter);

[nodePositionOptimal,displacementOptimal] = design_optimize_quantities_of_interest(...
    bottomUpMapping,...
    initialDesign,...
    designSpaceLowerBound,...
    designSpaceUpperBound,...
    @(performanceMeasure)[performanceMeasure(1)],...
    'InequalityConstraintFunction',@(performanceMeasure)[performanceMeasure(3:end)-performanceUpperLimit(3:end)],...
    'OptimizationMethodOptions',{'Display','iter-detailed'});


%% plot optimum and deformations
baseNode = [...
    0 0;  ...
    nan nan; ...
    2 0.5; ...
    nan nan; ...
    0 1];
fixedDegreesOfFreedom = [...
    true true; ...
    false false; ...
    false false; ...
    false false; ...
    true true];
nodeForce = [...
    0 0; ...
    0 0; ...
    0 -1000; ...
    0 0; ...
    0 0];
nodeElement = [...
    1 2; ...
    2 3; ...
    3 4; ...
    4 5; ...
    2 5; ...
    2 4];
elementCrossSectionArea = systemParameter(1); % [mm^2]
elementYoungsModulus = systemParameter(2); % [MPa]

% initial truss
plot_results_truss_six_bar_2d_two_moving_node(baseNode,nodeElement,...
    'NodePositionInitial',initialDesign);
save_print_figure(gcf,[saveFolder,'InitialTruss']);

% initial truss + optimized truss
plot_results_truss_six_bar_2d_two_moving_node(baseNode,nodeElement,...
    'NodePositionInitial',initialDesign,...
    'NodePositionOptimal',nodePositionOptimal);
save_print_figure(gcf,[saveFolder,'InitialOptimizedTruss']);


% deformed trusses
nodePositionInitial = baseNode;
nodePositionInitial(2,:) = initialDesign(:,[1,2]);
nodePositionInitial(4,:) = initialDesign(:,[3,4]);
nodeDisplacementInitial = ...
	truss_analysis(...
		nodePositionInitial,...
		fixedDegreesOfFreedom,...
		nodeForce,...
		nodeElement,...
		elementCrossSectionArea,...
		elementYoungsModulus);

nodePositionOptimized = baseNode;
nodePositionOptimized(2,:) = nodePositionOptimal(:,[1,2]);
nodePositionOptimized(4,:) = nodePositionOptimal(:,[3,4]);
nodeDisplacementOptimal = ...
    truss_analysis(...
	    nodePositionOptimized,...
	    fixedDegreesOfFreedom,...
	    nodeForce,...
	    nodeElement,...
	    elementCrossSectionArea,...
	    elementYoungsModulus);

plot_results_truss_six_bar_2d_two_moving_node(baseNode,nodeElement,...
    'NodePositionInitial',initialDesign,...
    'NodePositionOptimal',nodePositionOptimal,...
    'NodeDisplacementInitial',nodeDisplacementInitial,...
    'NodeDisplacementOptimal',nodeDisplacementOptimal);
save_print_figure(gcf,[saveFolder,'InitialOptimizedTrussDeformation']);


%% establish upper performance limits
% update uppwer limit based on either optimal value or initial value
performanceMeasureInitial = bottomUpMapping.response(initialDesign);

% displacement
performanceUpperLimit(1) = performanceMeasureInitial(1);
...performanceUpperLimit(1) = displacementOptimal*1.1;

% mass
...performanceUpperLimit(2) = performanceMeasureInitial(2);

designEvaluator = DesignEvaluatorBottomUpMapping(...
        bottomUpMapping,...
        performanceLowerLimit,...
        performanceUpperLimit);


%% Box Opt
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
    'TrimmingOperationOptions',{'PassesCriterion','full'},...
    'TrimmingOrderOptions',{'OrderPreference','score'});

rng(rngState);
[solutionSpaceBox,problemDataBox,iterDataBox] = sso_box_stochastic(designEvaluator,...
    initialDesign,designSpaceLowerBound,designSpaceUpperBound,optionsBox);
toc(timeElapsedBox)


%% Component Opt
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
[componentSolutionSpace,problemDataComponent,iterDataComponent] = sso_component_stochastic(designEvaluator,...
    initialDesign,designSpaceLowerBound,designSpaceUpperBound,Components,optionsComponent);
toc(timeElapsedComponent)


%% measure comparison
measureComponent1 = componentSolutionSpace(1).Measure;
measureComponent2 = componentSolutionSpace(2).Measure;
measureBox1 = prod(solutionSpaceBox(2,[1,2])-solutionSpaceBox(1,[1,2]));
measureBox2 = prod(solutionSpaceBox(2,[3,4])-solutionSpaceBox(1,[3,4]));
fprintf('\nVolume Increases - Node 1: %.3gx ; Node 2: %.gx\n',measureComponent1/measureBox1,measureComponent2/measureBox2);


%% Plot Visualization
% initial truss + component solution spaces
plot_results_truss_six_bar_2d_two_moving_node(baseNode,nodeElement,...
    'NodePositionInitial',initialDesign,...
    'ComponentSolutionSpace',componentSolutionSpace);
save_print_figure(gcf,[saveFolder,'InitialTrussComponent']);

% initial truss + optimized truss + box solution space + component solution spaces
plot_results_truss_six_bar_2d_two_moving_node(baseNode,nodeElement,...
    'NodePositionInitial',initialDesign,...
    'NodePositionOptimal',nodePositionOptimal,...
    'BoxSolutionSpace',solutionSpaceBox,...
    'ComponentSolutionSpace',componentSolutionSpace);
save_print_figure(gcf,[saveFolder,'InitialOptimizedTrussBoxComponent']);

% initial truss + box solution space + component solution spaces
plot_results_truss_six_bar_2d_two_moving_node(baseNode,nodeElement,...
    'NodePositionInitial',initialDesign,...
    'BoxSolutionSpace',solutionSpaceBox,...
    'ComponentSolutionSpace',componentSolutionSpace);
save_print_figure(gcf,[saveFolder,'InitialTrussBoxComponent']);

% initial truss + sample trusses + component solution space
nRandomTruss = 5;
randomTrussMovingNode = candidate_space_sampling_individual_feasible(componentSolutionSpace,Components,nRandomTruss);

plot_results_truss_six_bar_2d_two_moving_node(baseNode,nodeElement,...
    'NodePositionInitial',initialDesign,...
    'NodePositionRandom',randomTrussMovingNode,...
    'ComponentSolutionSpace',componentSolutionSpace);
save_print_figure(gcf,[saveFolder,'InitialRandomTrussComponent']);


%% plot candidate spaces samples
% candidate space 1
designSample = componentSolutionSpace(1).DesignSampleDefinition;
labelSample = componentSolutionSpace(1).IsInsideDefinition;
figure;
componentSolutionSpace(1).plot_candidate_space(gcf,'EdgeColor','none','FaceColor','g','FaceAlpha',0.5,'Linewidth',2.0);
hold all;
plot(designSample(labelSample,1),designSample(labelSample,2),'g.');
plot(designSample(~labelSample,1),designSample(~labelSample,2),'r.');
grid minor;
save_print_figure(gcf,[saveFolder,'ComponentSpace1TrimmingPlot']);

% candidate space 2
designSample = componentSolutionSpace(2).DesignSampleDefinition;
labelSample = componentSolutionSpace(2).IsInsideDefinition;
figure;
componentSolutionSpace(2).plot_candidate_space(gcf,'EdgeColor','none','FaceColor','g','FaceAlpha',0.5,'Linewidth',2.0);
hold all;
plot(designSample(labelSample,1),designSample(labelSample,2),'g.');
plot(designSample(~labelSample,1),designSample(~labelSample,2),'r.');
grid minor;
save_print_figure(gcf,[saveFolder,'ComponentSpace2TrimmingPlot']);


%% 
algoDataBox = postprocess_sso_box_stochastic(problemDataBox,iterDataBox);
plot_sso_box_stochastic_metrics(algoDataBox,'SaveFolder',saveFolder,'CloseFigureAfterSaving',true);

algoDataComponent = postprocess_sso_component_stochastic(problemDataComponent,iterDataComponent);
plot_sso_component_stochastic_metrics(algoDataComponent,'SaveFolder',saveFolder,'CloseFigureAfterSaving',true);



%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;


%% subfunction: plot truss
function figureHandle = plot_results_truss_six_bar_2d_two_moving_node(baseNode,nodeElement,varargin)
    parser = inputParser;
    parser.addParameter('NodePositionInitial',[]);
    parser.addParameter('NodePositionOptimal',[]);
    parser.addParameter('NodePositionRandom',[]);
    parser.addParameter('NodeDisplacementInitial',[]);
    parser.addParameter('NodeDisplacementOptimal',[]);
    parser.addParameter('BoxSolutionSpace',[]);
    parser.addParameter('ComponentSolutionSpace',[]);
    parser.addParameter('IncludeWall',true);
    parser.addParameter('IncludeAppliedForce',false);
    parser.addParameter('IncludeAxesInformation',false);
    parser.addParameter('IncludeLegend',false);
    parser.addParameter('WallOptions',{});
    parser.addParameter('AppliedForceOptions',{});
    parser.addParameter('InitialTrussOptions',{});
    parser.addParameter('OptimalTrussOptions',{});
    parser.addParameter('RandomTrussOptions',{});
    parser.addParameter('BoxSolutionSpaceOptions',{});
    parser.addParameter('ComponentSolutionSpaceOptions',{});
    parser.addParameter('LegendOptions',{});
    parser.parse(varargin{:});
    options = parser.Results;

    wallX = [0 0];
    wallY = [-0.1 1.1];

    defaultInitialTrussOptions = {'ColorUndeformed',[0.8 0.8 0.8],'ColorDeformed','c','MaximumLinewidth',3.0,'DisplacementScaleFactor',500};
    [~,initialTrussOptions] = merge_name_value_pair_argument(defaultInitialTrussOptions,options.InitialTrussOptions);

    defaultOptimalTrussOptions = {'ColorUndeformed','b','ColorDeformed','m','MaximumLinewidth',3.0,'DisplacementScaleFactor',500};
    [~,optimalTrussOptions] = merge_name_value_pair_argument(defaultOptimalTrussOptions,options.OptimalTrussOptions);

    defaultRandomTrussOptions = {'MaximumLinewidth',3.0};
    [~,randomTrussOptions] = merge_name_value_pair_argument(defaultRandomTrussOptions,options.RandomTrussOptions);

    defaultWallOptions = {'Linewidth',8.0,'Color',[0.7 0.7 0.7]};
    [~,wallOptions] = merge_name_value_pair_argument(defaultWallOptions,options.WallOptions);

    defaultAppliedForceOptions = {'Color','r','LineWidth',3.0};
    [~,appliedForceOptions] = merge_name_value_pair_argument(defaultAppliedForceOptions,options.AppliedForceOptions);

    defaultBoxSolutionOptions = {'EdgeColor','c','Linewidth',2.0};
    [~,boxSolutionOptions] = merge_name_value_pair_argument(defaultBoxSolutionOptions,options.BoxSolutionSpaceOptions);

    defaultComponentSolutionOptions = {'EdgeColor','g','FaceColor','none','FaceAlpha',0.5,'Linewidth',2.0};
    [~,componentSolutionOptions] = merge_name_value_pair_argument(defaultComponentSolutionOptions,options.ComponentSolutionSpaceOptions);

    defaultLegendOptions = {'location','west'};
    [~,legendOptions] = merge_name_value_pair_argument(defaultLegendOptions,options.LegendOptions);

    figureHandle = figure;
    hold all;
    
    % initial truss
    handleInitial = [];
    handleInitialDeformed = [];
    if(~isempty(options.NodePositionInitial))
        nodePositionInitial = baseNode;
        nodePositionInitial(2,:) = options.NodePositionInitial(:,[1,2]);
        nodePositionInitial(4,:) = options.NodePositionInitial(:,[3,4]);

        [handleInitial,handleInitialDeformed] = plot_truss_deformation(gcf,nodePositionInitial,nodeElement,options.NodeDisplacementInitial,initialTrussOptions{:});
    end

    % optimal truss
    handleOptimal = [];
    handleOptimizedDeformed = [];
    if(~isempty(options.NodePositionOptimal))
        nodePositionOptimized = baseNode;
        nodePositionOptimized(2,:) = options.NodePositionOptimal(:,[1,2]);
        nodePositionOptimized(4,:) = options.NodePositionOptimal(:,[3,4]);

        [handleOptimal,handleOptimizedDeformed] = plot_truss_deformation(gcf,nodePositionOptimized,nodeElement,options.NodeDisplacementOptimal,optimalTrussOptions{:});
    end

    % box-shaped solution space
    handleToleranceNodeBox = [];
    if(~isempty(options.BoxSolutionSpace))
        plot_design_box_2d(gcf,options.BoxSolutionSpace(:,[1,2]),boxSolutionOptions{:});
        handleToleranceNodeBox = plot_design_box_2d(gcf,options.BoxSolutionSpace(:,[3,4]),boxSolutionOptions{:});
    end

    % component solution space
    handleToleranceNodeComponent = [];
    if(~isempty(options.ComponentSolutionSpace))
        options.ComponentSolutionSpace(1).plot_candidate_space(gcf,componentSolutionOptions{:});
        handleToleranceNodeComponent = options.ComponentSolutionSpace(2).plot_candidate_space(gcf,componentSolutionOptions{:});
    end

    % random trusses
    handleRandom = [];
    if(~isempty(options.NodePositionRandom))
        nRandomTruss = size(options.NodePositionRandom,1);
        trussColor = rand(nRandomTruss,3);
        for i=1:nRandomTruss
            randomNodePosition = baseNode;
            randomNodePosition(2,:) = options.NodePositionRandom(i,[1,2]);
            randomNodePosition(4,:) = options.NodePositionRandom(i,[3,4]);
            handleRandom(i) = plot_truss_deformation(gcf,randomNodePosition,nodeElement,'ColorUndeformed',trussColor(i,:),randomTrussOptions{:});
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
        handleForce = quiver(2,0.5,0,-0.5,appliedForceOptions{:});
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
            handleOptimal,...
            handleOptimizedDeformed,...
            handleWall,...
            handleForce,...
            handleToleranceNodeBox,...
            handleToleranceNodeComponent};

        legendTextAll = {...
            'Initial Truss',...
            'Initial Truss (Deformed)',...
            'Optimized Truss',...
            'Optimized Truss (Deformed)',...
            'Wall',...
            'Applied Force',...
            'Tolerance Region for the Node (Box)',...
            'Tolerance Region for the Node (Component)'};
        
        handleObject = [handleObjectAll{:}];
        legentText = {legendTextAll{~cellfun(@isempty,handleObjectAll)}};

        legend(handleObject,legentText,legendOptions{:});
    end

    if(nargout<1)
        clear figureHandle;
    end
end

