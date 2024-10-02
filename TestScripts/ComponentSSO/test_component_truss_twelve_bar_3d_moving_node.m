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
systemFunction = @truss_twelve_bar_3d_moving_node;
systemParameter = [100,210e3,795,7850e-9]; % [mm^2],[MPa],[mm^4],[kg/mm^3]
%                         x1   y1   z1  x2   y2  z2  x3  y3   z3
designSpaceLowerBound = [0.5 -0.5 -1.0 0.1 -0.5 0.5 0.1 0.0 -1.0];
designSpaceUpperBound = [2.0  0.5  0.5 2.0  0.5 1.5 2.0 2.0  1.0];
Components = {[1,2,3]',[4,5,6]',[7,8,9]'};
%
performanceLowerLimit = -inf;
performanceUpperLimit = [nan nan repmat(250,1,12) ones(1,12)];
%                x1 y1 z1 x2  y2 z2 x3 y3 z3
initialDesign = [ 1  0  0  1   0  1  1  1  0];


%% find optimum
bottomUpMapping = BottomUpMappingFunction(systemFunction,systemParameter);

[nodePositionOptimal,displacementOptimal] = design_optimize_quantities_of_interest(...
    bottomUpMapping,...
    initialDesign,...
    designSpaceLowerBound,...
    designSpaceUpperBound,...
    @(performanceMeasure)[performanceMeasure(1)],...
    'InequalityConstraintFunction',@(performanceMeasure)[performanceMeasure(3:end)-performanceUpperLimit(3:end)],...
    'OptimizationMethodFunction',@optimization_ga_wrapper,...
    'OptimizationMethodOptions',{'Display','iter'});


%% plot optimum and deformations
baseNode = [...
      0   0   0; % (1)
      0   0   1; % (2)
      0   1   0; % (3)
    nan nan nan; % (4)
    nan nan nan; % (5) 
    nan nan nan; % (6)
      2 0.5 0.5]; % (7)
fixedDegreesOfFreedom = [...
    true true true; % (1) 
    true true true; % (2)
    true true true; % (3)
    false false false; % (4)
    false false false; % (5)
    false false false; % (6)
    false false false]; % (7)
nodeForce = [...
    0 0     0; % (1)
    0 0     0; % (2)
    0 0     0; % (3)
    0 0     0; % (4)
    0 0     0; % (5)
    0 0     0; % (6)
    0 0 -1000]; % (7)
nodeElement = [...
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
elementCrossSectionArea = systemParameter(1); % [mm^2]
elementYoungsModulus = systemParameter(2); % [MPa]

% initial truss
plot_results_truss_twelve_bar_3d_moving_node(baseNode,nodeElement,...
    'NodePositionInitial',initialDesign);
save_print_figure(gcf,[saveFolder,'InitialTruss']);
rotating_video(gcf,[saveFolder,'InitialTruss']);
rotating_gif(gcf,[saveFolder,'InitialTruss']);

% initial truss + optimized truss
plot_results_truss_twelve_bar_3d_moving_node(baseNode,nodeElement,...
    'NodePositionInitial',initialDesign,...
    'NodePositionOptimal',nodePositionOptimal);
save_print_figure(gcf,[saveFolder,'InitialOptimizedTruss']);
rotating_video(gcf,[saveFolder,'InitialOptimizedTruss']);
rotating_gif(gcf,[saveFolder,'InitialOptimizedTruss']);

% deformed trusses
nodePositionInitial = baseNode;
nodePositionInitial(4,:) = initialDesign(1:3);
nodePositionInitial(5,:) = initialDesign(4:6);
nodePositionInitial(6,:) = initialDesign(7:9);
nodeDisplacementInitial = ...
	truss_analysis(...
		nodePositionInitial,...
		fixedDegreesOfFreedom,...
		nodeForce,...
		nodeElement,...
		elementCrossSectionArea,...
		elementYoungsModulus);

nodePositionOptimized = baseNode;
nodePositionOptimized(4,:) = nodePositionOptimal(1:3);
nodePositionOptimized(5,:) = nodePositionOptimal(4:6);
nodePositionOptimized(6,:) = nodePositionOptimal(7:9);
nodeDisplacementOptimal = ...
    truss_analysis(...
	    nodePositionOptimized,...
	    fixedDegreesOfFreedom,...
	    nodeForce,...
	    nodeElement,...
	    elementCrossSectionArea,...
	    elementYoungsModulus);

plot_results_truss_twelve_bar_3d_moving_node(baseNode,nodeElement,...
    'NodePositionInitial',initialDesign,...
    'NodePositionOptimal',nodePositionOptimal,...
    'NodeDisplacementInitial',nodeDisplacementInitial,...
    'NodeDisplacementOptimal',nodeDisplacementOptimal);
save_print_figure(gcf,[saveFolder,'InitialOptimizedTrussDeformation']);
rotating_video(gcf,[saveFolder,'InitialOptimizedTrussDeformation']);
rotating_gif(gcf,[saveFolder,'InitialOptimizedTrussDeformation']);


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
    'GrowthRate',0.1,...
    'ApplyLeanness','never',...
    'TrimmingOperationOptions',{'PassesCriterion','reduced'},...
    'TrimmingOrderOptions',{'OrderPreference','score'});

rng(rngState);
[solutionSpaceBox,problemDataBox,iterDataBox] = sso_box_stochastic(designEvaluator,...
    initialDesign,designSpaceLowerBound,designSpaceUpperBound,optionsBox);
toc(timeElapsedBox)


%% component opt
timeElapsedComponent = tic;
options = sso_stochastic_options('component',...
    'NumberSamplesPerIterationExploration',300,...
    'NumberSamplesPerIterationConsolidation',300,...
    'FixIterNumberExploration',true,...
    'FixIterNumberConsolidation',true,...
    'MaxIterExploration',30,...
    'MaxIterConsolidation',30,...
    ... 'CandidateSpaceConstructorExploration',@CandidateSpaceDelaunay,...
    ... 'CandidateSpaceConstructorConsolidation',@CandidateSpaceDelaunay,...
    ... 'TrimmingMethodFunction',@component_trimming_method_corner_box_removal,...
    'CandidateSpaceConstructorExploration',@CandidateSpaceConvexHull,...
    'CandidateSpaceConstructorConsolidation',@CandidateSpaceConvexHull,...
    'TrimmingMethodFunction',@component_trimming_method_planar_trimming,...
    ... 'TrimmingMethodOptions',{'ReferenceDesigns','boundary-center'},...
    'UseAdaptiveGrowthRate',false,...
    'GrowthRate',0.1,...
    'ApplyLeanness','never',...
    'UsePaddingSamplesInTrimming',true,...
    'UsePreviousEvaluatedSamplesConsolidation',false,...
    'UsePreviousPaddingSamplesConsolidation',false,...
    'TrimmingOperationOptions',{'PassesCriterion','reduced'},...
    'TrimmingOrderOptions',{'OrderPreference','score'});

rng(rngState);
[componentSolutionSpace,problemData,iterData] = sso_component_stochastic(designEvaluator,...
    initialDesign,designSpaceLowerBound,designSpaceUpperBound,Components,options);
toc(timeElapsedComponent)


%% measure comparison
measureComponent1 = componentSolutionSpace(1).Measure;
measureComponent2 = componentSolutionSpace(2).Measure;
measureComponent3 = componentSolutionSpace(3).Measure;
measureBox1 = prod(solutionSpaceBox(2,[1,2,3])-solutionSpaceBox(1,[1,2,3]));
measureBox2 = prod(solutionSpaceBox(2,[4,5,6])-solutionSpaceBox(1,[4,5,6]));
measureBox3 = prod(solutionSpaceBox(2,[7,8,9])-solutionSpaceBox(1,[7,8,9]));
fprintf('\nVolume Increases - Node 1: %.3gx ; Node 2: %.gx ; Node 3: %.gx\n',...
    measureComponent1/measureBox1,measureComponent2/measureBox2,measureComponent3/measureBox3);


%% Plot Visualization
% initial truss + component solution spaces
plot_results_truss_twelve_bar_3d_moving_node(baseNode,nodeElement,...
    'NodePositionInitial',initialDesign,...
    'ComponentSolutionSpace',componentSolutionSpace);
save_print_figure(gcf,[saveFolder,'InitialTrussComponent']);
rotating_video(gcf,[saveFolder,'InitialTrussComponent']);
rotating_gif(gcf,[saveFolder,'InitialTrussComponent']);

% initial truss + optimized truss + box solution space + component solution spaces
plot_results_truss_twelve_bar_3d_moving_node(baseNode,nodeElement,...
    'NodePositionInitial',initialDesign,...
    'NodePositionOptimal',nodePositionOptimal,...
    'BoxSolutionSpace',solutionSpaceBox,...
    'ComponentSolutionSpace',componentSolutionSpace);
save_print_figure(gcf,[saveFolder,'InitialOptimizedTrussBoxComponent']);
rotating_video(gcf,[saveFolder,'InitialOptimizedTrussBoxComponent']);
rotating_gif(gcf,[saveFolder,'InitialOptimizedTrussBoxComponent']);

% initial truss + box solution space + component solution spaces
plot_results_truss_twelve_bar_3d_moving_node(baseNode,nodeElement,...
    'NodePositionInitial',initialDesign,...
    'BoxSolutionSpace',solutionSpaceBox,...
    'ComponentSolutionSpace',componentSolutionSpace);
save_print_figure(gcf,[saveFolder,'InitialTrussBoxComponent']);
rotating_video(gcf,[saveFolder,'InitialTrussBoxComponent']);
rotating_gif(gcf,[saveFolder,'InitialTrussBoxComponent']);

% initial truss + sample trusses + component solution space
nRandomTruss = 5;
randomTrussMovingNode = candidate_space_sampling_individual_feasible(componentSolutionSpace,Components,nRandomTruss);

plot_results_truss_twelve_bar_3d_moving_node(baseNode,nodeElement,...
    'NodePositionInitial',initialDesign,...
    'NodePositionRandom',randomTrussMovingNode,...
    'ComponentSolutionSpace',componentSolutionSpace);
save_print_figure(gcf,[saveFolder,'InitialRandomTrussComponent']);
rotating_video(gcf,[saveFolder,'InitialRandomTrussComponent']);
rotating_gif(gcf,[saveFolder,'InitialRandomTrussComponent']);


%% plot candidate spaces samples
% candidate space 1
designSample = componentSolutionSpace(1).DesignSampleDefinition;
labelSample = componentSolutionSpace(1).IsInsideDefinition;
figure;
componentSolutionSpace(1).plot_candidate_space(gcf,'EdgeColor','none','FaceColor','g','FaceAlpha',0.5,'Linewidth',2.0);
hold all;
plot3(designSample(labelSample,1),designSample(labelSample,2),designSample(labelSample,3),'g.');
plot3(designSample(~labelSample,1),designSample(~labelSample,2),designSample(~labelSample,3),'r.');
grid minor;
xlabel('x');
ylabel('y');
zlabel('z');
axis('equal');
axis('vis3d');
view(3);
save_print_figure(gcf,[saveFolder,'ComponentSpace1TrimmingPlot']);

% candidate space 2
designSample = componentSolutionSpace(2).DesignSampleDefinition;
labelSample = componentSolutionSpace(2).IsInsideDefinition;
figure;
componentSolutionSpace(2).plot_candidate_space(gcf,'EdgeColor','none','FaceColor','g','FaceAlpha',0.5,'Linewidth',2.0);
hold all;
plot3(designSample(labelSample,1),designSample(labelSample,2),designSample(labelSample,3),'g.');
plot3(designSample(~labelSample,1),designSample(~labelSample,2),designSample(~labelSample,3),'r.');
grid minor;
xlabel('x');
ylabel('y');
zlabel('z');
axis('equal');
axis('vis3d');
view(3);
save_print_figure(gcf,[saveFolder,'ComponentSpace2TrimmingPlot']);

% candidate space 3
designSample = componentSolutionSpace(3).DesignSampleDefinition;
labelSample = componentSolutionSpace(3).IsInsideDefinition;
figure;
componentSolutionSpace(3).plot_candidate_space(gcf,'EdgeColor','none','FaceColor','g','FaceAlpha',0.5,'Linewidth',2.0);
hold all;
plot3(designSample(labelSample,1),designSample(labelSample,2),designSample(labelSample,3),'g.');
plot3(designSample(~labelSample,1),designSample(~labelSample,2),designSample(~labelSample,3),'r.');
grid minor;
xlabel('x');
ylabel('y');
zlabel('z');
axis('equal');
axis('vis3d');
view(3);
save_print_figure(gcf,[saveFolder,'ComponentSpace3TrimmingPlot']);


%% 
algoDataBox = postprocess_sso_box_stochastic(problemDataBox,iterDataBox);
plot_sso_box_stochastic_metrics(algoDataBox,'SaveFolder',saveFolder,'CloseFigureAfterSaving',true);

algoData = postprocess_sso_component_stochastic(problemData,iterData);
plot_sso_component_stochastic_metrics(algoData,'SaveFolder',saveFolder,'CloseFigureAfterSaving',true);


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;


%% subfunction: plot truss
function figureHandle = plot_results_truss_twelve_bar_3d_moving_node(baseNode,nodeElement,varargin)
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

    [wallY,wallZ] = meshgrid(-0.5:1.5);
    wallX = zeros(size(wallY));

    defaultInitialTrussOptions = {'ColorUndeformed',[0.8 0.8 0.8],'ColorDeformed','c','MaximumLinewidth',4.0,'DisplacementScaleFactor',500};
    [~,initialTrussOptions] = merge_name_value_pair_argument(defaultInitialTrussOptions,options.InitialTrussOptions);

    defaultOptimalTrussOptions = {'ColorUndeformed','b','ColorDeformed','m','MaximumLinewidth',4.0,'DisplacementScaleFactor',500};
    [~,optimalTrussOptions] = merge_name_value_pair_argument(defaultOptimalTrussOptions,options.OptimalTrussOptions);

    defaultRandomTrussOptions = {'MaximumLinewidth',4.0};
    [~,randomTrussOptions] = merge_name_value_pair_argument(defaultRandomTrussOptions,options.RandomTrussOptions);

    defaultWallOptions = {'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.9};
    [~,wallOptions] = merge_name_value_pair_argument(defaultWallOptions,options.WallOptions);

    defaultAppliedForceOptions = {'Color','r','LineWidth',3.0};
    [~,appliedForceOptions] = merge_name_value_pair_argument(defaultAppliedForceOptions,options.AppliedForceOptions);

    defaultBoxSolutionOptions = {'FaceAlpha',0.2,'FaceColor','c'};
    [~,boxSolutionOptions] = merge_name_value_pair_argument(defaultBoxSolutionOptions,options.BoxSolutionSpaceOptions);

    defaultComponentSolutionOptions = {'FaceAlpha',0.2,'FaceColor','g'};
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
        nodePositionInitial(4,:) = options.NodePositionInitial(:,[1,2,3]);
        nodePositionInitial(5,:) = options.NodePositionInitial(:,[4,5,6]);
        nodePositionInitial(6,:) = options.NodePositionInitial(:,[7,8,9]);

        [handleInitial,handleInitialDeformed] = plot_truss_deformation(gcf,nodePositionInitial,nodeElement,options.NodeDisplacementInitial,initialTrussOptions{:});
    end

    % optimal truss
    handleOptimal = [];
    handleOptimizedDeformed = [];
    if(~isempty(options.NodePositionOptimal))
        nodePositionOptimized = baseNode;
        nodePositionOptimized(4,:) = options.NodePositionOptimal(:,[1,2,3]);
        nodePositionOptimized(5,:) = options.NodePositionOptimal(:,[4,5,6]);
        nodePositionOptimized(6,:) = options.NodePositionOptimal(:,[7,8,9]);

        [handleOptimal,handleOptimizedDeformed] = plot_truss_deformation(gcf,nodePositionOptimized,nodeElement,options.NodeDisplacementOptimal,optimalTrussOptions{:});
    end

    % box-shaped solution space
    handleToleranceNodeBox = [];
    if(~isempty(options.BoxSolutionSpace))
        plot_design_box_3d(gcf,options.BoxSolutionSpace(:,[1,2,3]),boxSolutionOptions{:});
        plot_design_box_3d(gcf,options.BoxSolutionSpace(:,[4,5,6]),boxSolutionOptions{:});
        handleToleranceNodeBox = plot_design_box_3d(gcf,options.BoxSolutionSpace(:,[7,8,9]),boxSolutionOptions{:});
    end

    % component solution space
    handleToleranceNodeComponent = [];
    if(~isempty(options.ComponentSolutionSpace))
        options.ComponentSolutionSpace(1).plot_candidate_space(gcf,componentSolutionOptions{:});
        options.ComponentSolutionSpace(2).plot_candidate_space(gcf,componentSolutionOptions{:});
        handleToleranceNodeComponent = options.ComponentSolutionSpace(3).plot_candidate_space(gcf,componentSolutionOptions{:});
    end

    % random trusses
    handleRandom = [];
    if(~isempty(options.NodePositionRandom))
        nRandomTruss = size(options.NodePositionRandom,1);
        trussColor = rand(nRandomTruss,3);
        for i=1:nRandomTruss
            randomNodePosition = baseNode;
            randomNodePosition(4,:) = options.NodePositionRandom(i,[1,2,3]);
            randomNodePosition(5,:) = options.NodePositionRandom(i,[4,5,6]);
            randomNodePosition(6,:) = options.NodePositionRandom(i,[7,8,9]);
            handleRandom(i) = plot_truss_deformation(gcf,randomNodePosition,nodeElement,'ColorUndeformed',trussColor(i,:),randomTrussOptions{:});
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
        handleForce = quiver3(2,0.5,0.5,0,0,-0.5,appliedForceOptions{:});
    end

    % axis
    grid minor;
    if(~options.IncludeAxesInformation)
        set(gca,'XColor', 'none','YColor','none','ZColor','none');
    else
        xlabel('x');
        ylabel('y');
        zlabel('z');
    end

    % adjust perspective
    axis('equal');
    axis('vis3d');
    camproj('perspective');
    cameratoolbar; % better adjust angle/perspective

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

function rotating_video(figureHandle,filename,varargin)
    parser = inputParser;
    parser.addParameter('Duration',5);
    parser.addParameter('FramesPerSecond',60);
    parser.addParameter('VideoWriterProfile','Motion JPEG 2000');
    parser.addParameter('VideoWriterOptions',{});
    parser.parse(varargin{:});
    options = parser.Results;

    defaultVideoWriterOptions = {'LosslessCompression',true,'FrameRate',options.FramesPerSecond};
    [~,videoWriterOptions] = merge_name_value_pair_argument(defaultVideoWriterOptions,options.VideoWriterOptions);

    % activate figure
    figure(figureHandle);
    [defaultAzimuth,defaultElevation] = view(3);

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
    parser.addParameter('FramesPerSecond',60);
    parser.addParameter('ExportgraphicsOptions',{})
    parser.parse(varargin{:});
    options = parser.Results;

    defaultExportgraphicsOptions = {};
    [~,exportgraphicsOptions] = merge_name_value_pair_argument(defaultExportgraphicsOptions,options.ExportgraphicsOptions);

    % activate figure
    figure(figureHandle);
    set(figureHandle, 'MenuBar', 'none');
    set(figureHandle, 'ToolBar', 'none');
    [defaultAzimuth,defaultElevation] = view(3);

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
                'DelayTime',frameDelay);
        else
            imwrite(A,map,...
                [filename,'.gif'],'gif',...
                'WriteMode','append', ...
                'DelayTime',frameDelay);
        end
    end

    % reset menu/tool bars
    set(figureHandle, 'MenuBar', 'figure');
    set(figureHandle, 'ToolBar', 'auto');
end