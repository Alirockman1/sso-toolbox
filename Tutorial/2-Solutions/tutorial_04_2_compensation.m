%TUTORIAL_04_2_COMPENSATION Solution-compensation spaces for a sphere problem 
%   TUTORIAL_04_2_COMPENSATION computes a solution box with solution-
%   compensation spaces for a sphere problem.

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

designEvaluatorBase = DesignEvaluatorBottomUpMapping(bottomUpMapping,...
    performanceLowerLimit,performanceUpperLimit);


% design space
designSpaceLowerBound = [-6 -6 -6];
designSpaceUpperBound = [ 6  6  6];
initialDesign = [3,0,0];

% compensation 
compensationAspaceIndex = [true,true,false];

[designEvaluator,aspaceLowerBound,aspaceUpperBound,aspaceInitialDesign] = DesignEvaluatorCompensation(...
    designEvaluatorBase,compensationAspaceIndex,designSpaceLowerBound,designSpaceUpperBound,initialDesign);


%% optimization
options = sso_stochastic_options('box',...
    'FixIterNumberExploration',true,...
    'FixIterNumberConsolidation',true,...
    'NumberSamplesPerIteration',100);

[designBoxNormal,problemDataNormal,iterDataNormal] = sso_box_stochastic(...
    designEvaluatorBase,...
    initialDesign,...
    designSpaceLowerBound,...
    designSpaceUpperBound,...
    options);

[designBoxCompensation,problemDataCompensation,iterDataCompensation] = sso_box_stochastic(...
    designEvaluator,...
    aspaceInitialDesign,...
    aspaceLowerBound,...
    aspaceUpperBound,...
    options);


%% Plot Results
% Analytical Solution
boxSideCompensation = performanceUpperLimit*sqrt(2);
compensationSolution = [-boxSideCompensation/2 ,  -boxSideCompensation/2 ; boxSideCompensation/2 , boxSideCompensation/2];
% Solution with normal Solution Spaces
boxSideNormal = performanceUpperLimit*sqrt(4/3);
normalSolution = [-boxSideNormal/2 ,  -boxSideNormal/2 ; boxSideNormal/2 , boxSideNormal/2];

% Plot
figure;
plot_design_box_2d(gcf,designBoxCompensation,'EdgeColor','k');
hold on;
grid minor;
plot_design_box_2d(gcf,compensationSolution,'EdgeColor','m','Linestyle','--');
plot_design_box_2d(gcf,designBoxNormal(:,compensationAspaceIndex),'EdgeColor','b');
plot_design_box_2d(gcf,normalSolution,'EdgeColor','r','Linestyle','--');
plot(initialDesign(1),initialDesign(2),'rx');
xlabel('x_1');
ylabel('x_2');
legend({'Solution Box (Compensation)','Analytical Solution (Compensation)',...
        'Solution Box (Normal)','Analytical Solution (Normal)',...
        'Initial Design'},'Location','eastoutside');
axis('equal');
save_print_figure(gcf,[saveFolder,'SolutionBox'],'Size',figureSize*1.5);


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

