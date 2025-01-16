%TUTORIAL_04_1_MULTIFIDELITY Solving the problem with a multi-fidelity approach
%   TUTORIAL_04_1_MULTIFIDELITY solves the SSO box problem with now a multi-
%   fidelity approach, where we first train a surrogate model for the problem,
%   and then combine it with the analytical evaluator.
%
%   Copyright 2025 Eduardo Rodrigues Della Noce
%   SPDX-License-Identifier: Apache-2.0

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


%% bottom-up mapping definition
%TODO: construct bottom-up  mapping from the system response function
% definition: bottomUpMapping

%% performance limits
%TODO: implement performance limits
% definition: performanceLowerLimit,performanceUpperLimit

%% design space
%TODO: define design space
% definition: designSpaceLowerBound,designSpaceUpperBound,initialDesign

%% Active Learning Options
%TODO: choose options for the active learning training 
% definition: options

%% Active learning training
%TODO: train a model using active learning
% definition: designFastForward,problemData,iterData

%% approximate the threshold value of score where errors start happening
%TODO: find an estimate as to the score threshold when errors start occuring
% definition: scoreThreshold,iterDataUncertainty


%% make multi-fidelity evaluator
%TODO: construct multiple evaluators: high-, low- and multi-fidelity
% definition: designEvaluatorHighFidelity, designEvaluatorLowFidelity, designEvaluatorMultiFidelity


%% solve sphere problem with three evaluators to compare results
options = sso_stochastic_options('box');

initialRngState = rng;
[designBoxHighFidelity,problemDataHighFidelity,iterDataHighFidelity] = ...
    sso_box_stochastic(...
    designEvaluatorHighFidelity,...
    initialDesign,...
    designSpaceLowerBoundary,...
    designSpaceUpperBoundary,...
    options);

rng(initialRngState);
[designBoxLowFidelity,problemDataLowFidelity,iterDataLowFidelity] = ...
    sso_box_stochastic(...
    designEvaluatorLowFidelity,...
    initialDesign,...
    designSpaceLowerBoundary,...
    designSpaceUpperBoundary,...
    options);

rng(initialRngState);
[designBoxMultiFidelity,problemDataMultiFidelity,iterDataMultiFidelity] = ...
    sso_box_stochastic(...
    designEvaluatorMultiFidelity,...
    initialDesign,...
    designSpaceLowerBoundary,...
    designSpaceUpperBoundary,...
    options);


%% known analytical solution
boxSideAnalytical = performanceUpperLimit*sqrt(4/3);
solutionAnalytical = [-boxSideAnalytical/2 ,  -boxSideAnalytical/2 ,  -boxSideAnalytical/2;...
    boxSideAnalytical/2 , boxSideAnalytical/2 , boxSideAnalytical/2];


%% check total number of function evaluations
functionEvaluationHighFidelity = size(vertcat(iterDataHighFidelity(:).EvaluatedDesignSamples),1);

functionEvaluationLowFidelity = size(vertcat(iterDataTraining(:).EvaluatedSamples),1);

functionEvaluationUncertain = size(vertcat(iterDataUncertainty(:).EvaluatedSamples),1);

functionEvaluationMultiFidelity = sum(vertcat(vertcat(iterDataMultiFidelity(:).EvaluationOutput).IsUncertain));


%% report in text
fprintf(['\n',repelem('=',80),'\n']);
fprintf(['%14s - High-Fidelity Function Evaluations = %g\n',...
         '%14s - High-Fidelity Function Evaluations = %g\n',...
         '%14s - High-Fidelity Function Evaluations = %g + %g + %g = %g\n'],...
         'High Fidelity',functionEvaluationHighFidelity,...
         'Low Fidelity',functionEvaluationLowFidelity,...
         'Multi-Fidelity',...
            functionEvaluationLowFidelity,functionEvaluationUncertain,functionEvaluationMultiFidelity,...
            functionEvaluationLowFidelity + functionEvaluationUncertain + functionEvaluationMultiFidelity);


%% see results
% sphere points
[sphereX,sphereY,sphereZ] = sphere;
sphereX = systemParameter(1)*sphereX;
sphereY = systemParameter(1)*sphereY;
sphereZ = systemParameter(1)*sphereZ;

figure;
hold all;
plot_ellipsoid_3d(gcf,systemParameter,performanceUpperLimit,'NumberPoints',30,'SurfOptions',{'FaceAlpha',0.1,'FaceColor','black','EdgeColor','none'});
plot_design_box_3d(gcf,designBoxHighFidelity,'FaceAlpha',0.1,'FaceColor','blue');
plot_design_box_3d(gcf,designBoxLowFidelity,'FaceAlpha',0.1,'FaceColor','yellow');
plot_design_box_3d(gcf,designBoxMultiFidelity,'FaceAlpha',0.1,'FaceColor','cyan');
plot_design_box_3d(gcf,solutionAnalytical,'FaceAlpha',0.1,'FaceColor','green');
xlabel('x_1');
ylabel('x_2');
zlabel('x_3');
grid minor;
legend({'Sphere','High-Fidelity','Low-Fidelity','Multi-Fidelity','Analytical'});

figure;
hold all;
subplot(2,2,1);
plot_design_box_2d(gcf,designBoxHighFidelity(:,[1,2]),'EdgeColor','blue','Linewidth',2.0);
plot_design_box_2d(gcf,designBoxLowFidelity(:,[1,2]),'EdgeColor','yellow','Linewidth',2.0);
plot_design_box_2d(gcf,designBoxMultiFidelity(:,[1,2]),'EdgeColor','cyan','Linewidth',2.0);
plot_design_box_2d(gcf,solutionAnalytical(:,[1,2]),'EdgeColor','green','Linewidth',2.0);
axis([-4 4 -4 4]);
xlabel('x_1');
ylabel('x_2');
grid minor;
subplot(2,2,2);
plot_design_box_2d(gcf,designBoxHighFidelity(:,[1,3]),'EdgeColor','blue','Linewidth',2.0);
plot_design_box_2d(gcf,designBoxLowFidelity(:,[1,3]),'EdgeColor','yellow','Linewidth',2.0);
plot_design_box_2d(gcf,designBoxMultiFidelity(:,[1,3]),'EdgeColor','cyan','Linewidth',2.0);
plot_design_box_2d(gcf,solutionAnalytical(:,[1,3]),'EdgeColor','green','Linewidth',2.0);
axis([-4 4 -4 4]);
xlabel('x_1');
ylabel('x_3');
grid minor;
subplot(2,2,4);
plot_design_box_2d(gcf,designBoxHighFidelity(:,[2,3]),'EdgeColor','blue','Linewidth',2.0);
plot_design_box_2d(gcf,designBoxLowFidelity(:,[2,3]),'EdgeColor','yellow','Linewidth',2.0);
plot_design_box_2d(gcf,designBoxMultiFidelity(:,[2,3]),'EdgeColor','cyan','Linewidth',2.0);
plot_design_box_2d(gcf,solutionAnalytical(:,[2,3]),'EdgeColor','green','Linewidth',2.0);
axis([-4 4 -4 4]);
xlabel('x_2');
ylabel('x_3');
grid minor;
legend({'High-Fidelity','Low-Fidelity','Multi-Fidelity','Analytical'});


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

