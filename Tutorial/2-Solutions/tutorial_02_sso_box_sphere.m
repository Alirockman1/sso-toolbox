%TUTORIAL_02_SSO_BOX_SPHERE Find optimal box inside of sphere
%   TUTORIAL_02_SSO_BOX_SPHERE shows how to use the functions to obtain a box-
%	shaped solution space and be able to visualize it.

%% Cleanup
fclose all;
close all;
clear all;
clc;
more off;
diary off;


%% debugging
rng(7);


%% Documentation / Archive
rngState = rng;
saveFolder = save_diary_files(mfilename);
goldenRatio = (1+sqrt(5))/2;
figureSize = [goldenRatio 1]*8.5;


%% function call
% system information
systemFunction = @tutorial_01_euclidean_distance_3d;
systemParameter = [0,0,0];

bottomUpMapping = BottomUpMappingFunction(systemFunction,systemParameter);

%
performanceLowerLimit = nan;
performanceUpperLimit = 5;

designEvaluator = DesignEvaluatorBottomUpMapping(bottomUpMapping,...
    performanceLowerLimit,performanceUpperLimit);


% design space
designSpaceLowerBound = [-6 -6 -6];
designSpaceUpperBound = [ 6  6  6];
initialDesign = [3,0,0];


%% optimization
options = sso_stochastic_options('box','SamplingMethodFunction',@sampling_random);

%%REMOVE
[designBox,problemData,iterData] = sso_box_stochastic(...
    designEvaluator,...
    initialDesign,...
    designSpaceLowerBound,...
    designSpaceUpperBound,...
    options);


%% Plot Results
% Solution with normal Solution Spaces
boxSideAnalytical = performanceUpperLimit*sqrt(4/3);
boxAnalyticalSolution = [...
    -boxSideAnalytical/2 ,  -boxSideAnalytical/2 ,  -boxSideAnalytical/2;...
    boxSideAnalytical/2 , boxSideAnalytical/2 , boxSideAnalytical/2];

% Plot
solutionSpaceHandle = [];
figure;
hold on;
goodRegionHandle = plot_ellipsoid_3d(gcf,systemParameter,performanceUpperLimit,...
    'NumberPoints',20,'SurfOptions',{'FaceAlpha',0.3,'FaceColor','g','EdgeColor','none'});

solutionSpaceHandle = plot_design_box_3d(gcf,designBox,'FaceAlpha',0.5,'FaceColor','b','EdgeColor','b');

analyticalHandle = plot_design_box_3d(gcf,boxAnalyticalSolution,'FaceAlpha',0.5,'FaceColor','r','EdgeColor','r','Linestyle','--');
initialDesignHandle = plot3(initialDesign(1),initialDesign(2),initialDesign(3),'rx');
xlabel('x_1');
ylabel('x_2');
ylabel('x_3');
legend(...
   [goodRegionHandle,solutionSpaceHandle,analyticalHandle,initialDesignHandle],...
    {'Limit of Good Designs','Solution Box','Analytical Solution','Initial Design'},...
    'Location','eastoutside');
grid minor;
axis('equal');
save_print_figure(gcf,[saveFolder,'SolutionBox'],'Size',figureSize*1.5);


%% Selective Design Space Projection
% desired plots by design variable index
plotTiles = [1,3]; % 2 rows, 3 columns of subplots
desiredPlots = [...
    1,2;...
    1,3;...
    2,3];
axisLabel = {'x_1','x_2','x_3',};

[h,PlotData] = plot_selective_design_space_projection(...
    designEvaluator,...
    designBox,...
    designSpaceLowerBound,...
    designSpaceUpperBound,...
    desiredPlots,...
    plotTiles,...
    'NumberSamplesPerPlot',10000,...
    'AxesLabels',axisLabel,...
    'MarkerColorsViolatedRequirements','r',...
    'PlotOptionsBad',{'Marker','x'});

legend([h.GoodPerformance,h.BadPerformance],...
    {'Good Designs','Violate Distance Requirement'});
sgtitle('3D Sphere - Box Decomposition');
save_print_figure(h.MainFigure,[saveFolder,'SelectiveDesignSpaceProjection']);
clear h


%% Performance Metrics
algoData = postprocess_sso_box_stochastic(problemData,iterData);
plot_sso_box_stochastic_metrics(algoData,...
    'SaveFolder',saveFolder,...
    'SaveFigureOptions',{'Size',figureSize});


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

