%TUTORIAL_10_COMPONENT_SOLUTION_SPACE_ADVANCED 3D Truss Component SSO
%   TUTORIAL_10_COMPONENT_SOLUTION_SPACE_ADVANCED finds tolerance regions for
%   node positions of a 3D truss. This is a non-linear problem and with 
%   components that are three dimensional - a new development.


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
systemParameter = [1000,210e3]; % [mm^2],[MPa]

%                         x1   y1   z1  x2   y2  z2  x3  y3   z3
designSpaceLowerBound = [0.5 -0.5 -0.5 0.5 -0.5 0.5 0.5 0.5 -0.5];
designSpaceUpperBound = [1.5  0.5  0.5 1.5  0.5 1.5 1.5 1.5  0.5];
componentIndex = {[1,2,3]',[4,5,6]',[7,8,9]'};

%
performanceLowerLimit = 0;
performanceUpperLimit = [nan repmat(250,1,12)];

%                x1 y1 z1 x2 y2 z2 x3 y3 z3
initialDesign = [ 1  0  0  1  0  1  1  1  0];

RequirementSpacesType = 'Omega1';


%% find optimum
bottomUpMapping = BottomUpMappingFunction(systemFunction,systemParameter);

[nodePositionOptimal,displacementOptimal] = design_optimize_quantities_of_interest(...
    bottomUpMapping,...
    initialDesign,...
    designSpaceLowerBound,...
    designSpaceUpperBound,...
    @(performanceMeasure)[performanceMeasure(1)],...
    'InequalityConstraintFunction',@(performanceMeasure)[performanceMeasure(2:end)-250],...
    'OptimizationMethodOptions',{'Display','iter-detailed'});


%% Component Opt - Function
% update uppwer limit based on optimal value
performanceUpperLimit(1) = displacementOptimal*1.1;
designEvaluator = DesignEvaluatorBottomUpMapping(...
        bottomUpMapping,...
        performanceLowerLimit,...
        performanceUpperLimit);

tic
optionsBox = sso_stochastic_options('box',...
    'NumberSamplesPerIteration',300,...
    'MaxIterExploration',30,...
    'MaxIterConsolidation',30,...
    'GrowthRate',0.1,...
    'TrimmingOperationOptions',{'PassesCriterion','reduced'});

[solutionSpaceBox,problemDataBox,iterDataBox] = sso_box_stochastic(designEvaluator,...
    nodePositionOptimal,designSpaceLowerBound,designSpaceUpperBound,optionsBox);
toc

tic
options = sso_stochastic_options('component',...
    'NumberSamplesPerIteration',300,...
    'MaxIterExploration',30,...
    'MaxIterConsolidation',30,...
    'GrowthRate',0.1,...
    'TrimmingOperationOptions',{'PassesCriterion','reduced'});

[componentSolutionSpace,problemData,iterData] = sso_component_stochastic(designEvaluator,...
    nodePositionOptimal,designSpaceLowerBound,designSpaceUpperBound,componentIndex,options);
toc


%% Plot Visualization
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

nodePositionInitial = [...
             0   0   0; % (1)
             0   0   1; % (2)
             0   1   0; % (3)
    initialDesign(1:3); % (4)
    initialDesign(4:6); % (5) 
    initialDesign(7:9); % (6)
             2 0.5 0.5]; % (7)

nodePositionOptimized = [...
                   0   0   0; % (1)
                   0   0   1; % (2)
                   0   1   0; % (3)
    nodePositionOptimal(1:3); % (4)
    nodePositionOptimal(4:6); % (5) 
    nodePositionOptimal(7:9); % (6)
                   2 0.5 0.5]; % (7)

figure;
hold all;
handleInitial = plot_truss_deformation(gcf,nodePositionInitial,nodeElement);
handleOptimal = plot_truss_deformation(gcf,nodePositionOptimized,nodeElement,'ColorUndeformed','r');
handleToleranceNode1Box = plot_design_box_3d(gcf,solutionSpaceBox(:,[1,2,3]),'FaceAlpha',0.2,'FaceColor','c');
handleToleranceNode2Box = plot_design_box_3d(gcf,solutionSpaceBox(:,[4,5,6]),'FaceAlpha',0.2,'FaceColor','c');
handleToleranceNode3Box = plot_design_box_3d(gcf,solutionSpaceBox(:,[7,8,9]),'FaceAlpha',0.2,'FaceColor','c');
handleToleranceNode1Component = componentSolutionSpace(1).plot_candidate_space(gcf,'FaceAlpha',0.2,'FaceColor','g');
handleToleranceNode2Component = componentSolutionSpace(2).plot_candidate_space(gcf,'FaceAlpha',0.2,'FaceColor','g');
handleToleranceNode3Component = componentSolutionSpace(3).plot_candidate_space(gcf,'FaceAlpha',0.2,'FaceColor','g');
grid minor;
legend([handleInitial,handleOptimal,...
    handleToleranceNode1Box,handleToleranceNode2Box,handleToleranceNode3Box,...
    handleToleranceNode1Component,handleToleranceNode2Component,handleToleranceNode3Component],...
    {'Initial Truss','Optimized Truss',...
    'Tolerance Region for Node 1 (Box)','Tolerance Region for Node 2 (Box)','Tolerance Region for Node 3 (Box)',...
    'Tolerance Region for Node 1 (Component)','Tolerance Region for Node 2 (Component)','Tolerance Region for Node 3 (Component)'});
save_print_figure(gcf,[saveFolder,'TrussVisualization']);


%% 
algoDataBox = postprocess_sso_box_stochastic(problemDataBox,iterDataBox);
plot_sso_box_stochastic_metrics(algoDataBox,'SaveFolder',saveFolder);

algoData = postprocess_sso_component_stochastic(problemData,iterData);
plot_sso_component_stochastic_metrics(algoData,'SaveFolder',saveFolder);


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

