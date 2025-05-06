%TUTORIAL_05_POINT_BASED_OPTIMAL_DESIGN Find optimal design for car crash 2D
%   TUTORIAL_05_POINT_BASED_OPTIMAL_DESIGN uses special optimization wrappers
%   to find the optimal design (based on pure performance and score) for the
%   2D car crash problem.


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


%% problem setup
systemFunction = @car_crash_2d;
%               m d1c d2c   v0
systemParameter = [2000 0.3 0.3 15.6];
%                       E_rem    a_max order
performanceLowerLimit = [-inf        0  -inf];
performanceUpperLimit = [   0  32*9.81     0];

bottomUpMapping = BottomUpMappingFunction(systemFunction,systemParameter);
designEvaluator = DesignEvaluatorBottomUpMapping(bottomUpMapping,...
    performanceLowerLimit,...
    performanceUpperLimit,...
    ...'PerformanceNormalizationFactor',[5.3e3 3.1392 1.7e4]...
    'PerformanceNormalizationFactor',[5.3e5 313.92 1.7e6]...
    );


%% design space
%                            F1     F2   
designSpaceLowerBound = [   1e5    1e5];
designSpaceUpperBound = [   1e6    1e6];
initialDesign         = [3.55e5 3.55e5];


%% optimization methods
optimizationFunction = @optimization_fmincon_wrapper;
optimizationOptions = {'Display','iter-detailed'};

% optimizationFunction = @optimization_sampling;
% optimizationOptions = {};

% optimizationFunction = @optimization_ga_wrapper;
% optimizationOptions = {'Display','diagnose'};

% optimizationFunction = @optimization_paretosearch_wrapper;
% optimizationOptions = {'Display','diagnose'};

% optimizationFunction = @optimization_patternsearch_wrapper;
% optimizationOptions = {'Display','diagnose'};


%% start point-based optimization
%TODO: do point-based optimization
% definition: designOptimalQuantityOfInterest, designOptimalScore


%% Plot Solution
figure;
goodTriangleX = [405000 627840 183000 405000];
goodTriangleY = [405000 627840 627840 405000];
plot(goodTriangleX,goodTriangleY,'g');
hold on;
grid minor;
plot(initialDesign(1),initialDesign(2),'rx');
plot(designOptimalQuantityOfInterest(1),designOptimalQuantityOfInterest(2),'bo');
plot(designOptimalScore(1),designOptimalScore(2),'mo');
xlabel('F_1 [N]');
ylabel('F_2 [N]');
legend('Good Region','Initial Design','Optimal (Quantities of Interest)','Optimal (Score)',...
    'location','eastoutside')
save_print_figure(gcf,[saveFolder,'OptimalSolutionByScore'],'Size',figureSize*1.2);


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

