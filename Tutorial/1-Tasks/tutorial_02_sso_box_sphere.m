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


%% bottom-up mapping definition
%TODO: construct bottom-up  mapping from the system response function
% definition: bottomUpMapping

%% design evaluator
%TODO: construct design evaluator from the bottom-up mapping
% definition: performanceLowerLimit,performanceUpperLimit,designEvaluator

%% design space
%TODO: define design space
% definition: designSpaceLowerBound,designSpaceUpperBound,initialDesign


%% box-shaped SSO
%TODO: solve the SSO problem
% definition: designBox,problemData,iterData


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

%TODO: visualize box in 3D
solutionSpaceHandle = [];

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
plotTiles = [1,3]; % 2 rows, 3 columns of subplots
desiredPlots = [...
    1,2;...
    1,3;...
    2,3];
axisLabel = {'x_1','x_2','x_3',};

%TODO: visualize with selective design space projection
% definition: figurePlotHandle


legend([figurePlotHandle.GoodPerformance,figurePlotHandle.BadPerformance],...
    {'Good Designs','Violate Distance Requirement'});
sgtitle('3D Sphere - Box Decomposition');
save_print_figure(figurePlotHandle.MainFigure,[saveFolder,'SelectiveDesignSpaceProjection']);
clear figurePlotHandle


%% Performance Metrics
%TODO: visualize evolution of performance metrics


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

