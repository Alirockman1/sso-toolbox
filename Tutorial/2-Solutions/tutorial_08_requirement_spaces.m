%TUTORIAL_08_REQUIREMENT_SPACES Requirement spaces computation with two ellipses
%   TUTORIAL_08_REQUIREMENT_SPACES uses a bottom-up mapping with two defined  
%   ellipses one which describes a region where performance requirements are 
%   satisfied and another which describes where designs are physically feasible,
%   to test the computation of requirement spaces with the stochastic algorithm.

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


%% function call
systemFunction = @two_ellipses_requirement_space;

designSpaceLowerBound = [-2 -2];
designSpaceUpperBound = [14 12];
initialPoint = [10 6];

performanceLowerLimit = -inf;
performanceUpperLimit = 1;
physicalFeasibilityLowerLimit = -inf;
physicalFeasibilityUpperLimit = 1;

bottomUpMapping = BottomUpMappingFunction(systemFunction);
designEvaluator = DesignEvaluatorBottomUpMapping(...
    bottomUpMapping,...
    performanceLowerLimit,...
    performanceUpperLimit,...
    'PhysicalFeasibilityLowerLimit',physicalFeasibilityLowerLimit,...
    'PhysicalFeasibilityUpperLimit',physicalFeasibilityUpperLimit);

requirementSpacesType = 'Omega2';
options = sso_stochastic_options('box',...
    'RequirementSpacesType',requirementSpacesType,...
    'ApplyLeanness','end-only',...
    'NumberSamplesPerIteration',200,...
    'FixIterNumberExploration',true,...
    'FixIterNumberConsolidation',true);

[designBox,problemData,iterData] = sso_box_stochastic(designEvaluator,initialPoint,...
    designSpaceLowerBound,designSpaceUpperBound,options);


%% Plot Solution
plot_ellipse_2d(figure,[11,5],[7,2],1.8,'PatchOptions',{'FaceColor','none','EdgeColor','g'});
hold all;
plot_ellipse_2d(gcf,[5,6],[7,2],0.2,'PatchOptions',{'FaceColor','none','EdgeColor','b'});
grid minor;
hold on;
plot_design_box_2d(gcf,designBox);
plot(initialPoint(1),initialPoint(2),'rx');
xlabel('x_1');
ylabel('x_2');
save_print_figure(gcf,[saveFolder,'BasicSolutionBox'],'Size',figureSize);


%% Performance Metrics
algoData = postprocess_sso_box_stochastic(problemData,iterData);
plot_sso_box_stochastic_metrics(algoData,...
    'SaveFolder',saveFolder,...
    'SaveFigureOptions',{'Size',figureSize});


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

