%TUTORIAL_09_COMPONENT_SOLUTION_SPACE_SIMPLE Component solution spaces sphere 
%   TUTORIAL_09_COMPONENT_SOLUTION_SPACE_SIMPLE computes a component solution 
%   spaces for a sphere problem.

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
RNGstate = rng;
savefolder = save_diary_files(mfilename);
goldenratio = (1+sqrt(5))/2;
figure_size = [goldenratio 1]*8.5;


%% function call
%
systemFunction = @tutorial_01_euclidean_distance_3d;
systemParameter = [0,0,0];
%                        x1 x2 x3
designSpaceLowerBound = [-6 -6 -6];
designSpaceUpperBound = [ 6  6  6];
%
performanceLowerLimit = -inf;
performanceUpperLimit = 5;
%
initialDesign = [0,0,0];

bottomUpMapping = BottomUpMappingFunction(systemFunction,systemParameter);
designEvaluator = DesignEvaluatorBottomUpMapping(...
    bottomUpMapping,...
    performanceLowerLimit,...
    performanceUpperLimit);


%% Component Opt - Function
%TODO: define component index, set the right options, and call the component SSO method


%% Plot Comparison
% analytical solution
radiusAnalytical = performanceUpperLimit*sqrt(2/3);
heightAnalytical = performanceUpperLimit*sqrt(1/3);

% component 1
designSampleDefinition = componentSolutionSpace(1).DesignSampleDefinition;
isInsideDefinition = componentSolutionSpace(1).IsInsideDefinition;
figure;
plot(designSampleDefinition(isInsideDefinition,1),designSampleDefinition(isInsideDefinition,2),'g.');
hold on;
plot(designSampleDefinition(~isInsideDefinition,1),designSampleDefinition(~isInsideDefinition,2),'.','color',[227 114 34]/255);
plot_ellipse_2d(gcf,systemParameter(componentIndex{1}),radiusAnalytical,'PatchOptions',{'FaceColor','none','EdgeColor','m','linewidth',2.0});
componentSolutionSpace(1).plot_candidate_space(gcf);
xlabel('x_1');
ylabel('x_3');
axis([-5 +5 -5 +5]);
grid minor;
legend({'Inside Component Space','Outside Component Space','Analytical Solution','Decision Boundary'});
save_print_figure(gcf,[savefolder,'Component1TrimmingPlot']);

% component 2
designSampleDefinition = componentSolutionSpace(2).DesignSampleDefinition;
isInsideDefinition = componentSolutionSpace(2).IsInsideDefinition;
figure;
plot(designSampleDefinition(isInsideDefinition,1),zeros(size(designSampleDefinition(isInsideDefinition,1))),'g.');
hold on;
plot(designSampleDefinition(~isInsideDefinition,1),zeros(size(designSampleDefinition(~isInsideDefinition,1))),'.','color',[227 114 34]/255);
plot([-heightAnalytical heightAnalytical],[0 0],'mx','linewidth',2.0);
xlabel('x_2');
axis([-5 +5 -1 +1]);
grid minor;
legend({'Inside Component Space','Outside Component Space','Analytical Solution'});
save_print_figure(gcf,[savefolder,'Component2TrimmingPlot']);


%% Performance metrics
%TODO: plot performance metrics for component SSO procedure


%% Save and Stop Transcripting
save([savefolder,'Data.mat']);
diary off;

