%TUTORIAL_06_PYTHON_BOTTOM_UP_MAPPING Solve same problem with Python response
%   TUTORIAL_06_PYTHON_BOTTOM_UP_MAPPING shows how to use a Python system
%   response function instead of MATLAB. 

%% Cleanup
close all;
fclose all;
clear all;
clc;
more off;
diary off;


%% debugging
rng(6);


%% Documentation / Archive
rngState = rng;
saveFolder = save_diary_files(mfilename);
goldenRatio = (1+sqrt(5))/2;
figureSize = [goldenRatio 1]*8.5;

%                     m d1c d2c   v0
systemParameter = [2000 0.3 0.3 15.6];

%% problem setup
% systemFunction = @car_crash_2d;
% bottomUpMapping = BottomUpMappingFunction(systemFunction,systemParameter);
systemFunction = 'car_crash_2d_python.py';
bottomUpMapping = BottomUpMappingPython(systemFunction,systemParameter);

%                            F1     F2 
designSpaceLowerBound = [   1e5    1e5];
designSpaceUpperBound = [   1e6    1e6];
initialDesign         = [4.05e5 4.05e5];

%                       E_rem    a_max order
performanceLowerLimit = [-inf        0  -inf];
performanceUpperLimit = [   0  32*9.81     0];


%% optimization call
designEvaluator = DesignEvaluatorBottomUpMapping(bottomUpMapping,performanceLowerLimit,performanceUpperLimit);

options = sso_stochastic_options('box',...
    'NumberSamplesPerIteration',100,...
    'FixIterNumberExploration',true,...
    'FixIterNumberConsolidation',true,...
    'UseAdaptiveGrowthRate',false);
[designBox,problemData,iterData] = sso_box_stochastic(...
    designEvaluator,initialDesign,designSpaceLowerBound,designSpaceUpperBound,options);


%% Plot Solution
figure;
goodRegionTriangleX = [405000 627840 183000 405000];
goodRegionTriangleY = [405000 627840 627840 405000];
plot(goodRegionTriangleX,goodRegionTriangleY,'g');
hold on;
grid minor;
plot_design_box_2d(gcf,designBox,'EdgeColor','k');
plot(initialDesign(1),initialDesign(2),'rx');
analyticalSolutionX = [291e3 516e3 516e3 291e3 291e3];
analyticalSolutionY = [516e3 516e3 628e3 628e3 516e3];
plot(analyticalSolutionX,analyticalSolutionY,'mo');
xlabel('F_1 [N]');
ylabel('F_2 [N]');
save_print_figure(gcf,[saveFolder,'SolutionBox'],'Size',figureSize);


%% Performance Metrics
algoData = postprocess_sso_box_stochastic(problemData,iterData);
plot_sso_box_stochastic_metrics(algoData,...
    'SaveFolder',saveFolder,...
    'SaveFigureOptions',{'Size',figureSize});


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

