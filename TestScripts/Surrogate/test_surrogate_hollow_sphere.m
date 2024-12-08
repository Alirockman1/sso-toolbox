%test_surrogate_hollow_sphere Surrogate model of the hollow sphere problem
%   test_surrogate_hollow_sphere creates a surrogate model for the 
%   hollow sphere problem. This is done using an active set strategy.
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
systemFunction = @distance_to_center;
systemParameter = [0,0,0];

%        x1  x2  x3
designSpaceLowerBoundary = [-6 -6 -6];
designSpaceUpperBoundary = [ 6  6  6];

performanceLowerLimit = 2;
performanceUpperLimit = 5;


%% Active Learning Options
options = active_learning_model_training_options(...
    'MaxIter',20,...
    'NumberSamplesEvaluation',20,...
    'FixIter',true,...
    'LoggingLevel','all',...
    'ActiveLearningSamplingFunction',@active_learning_sampling_line_search,...
    'ExplorationToBoundaryRatio',0.0,...
    'FastForwardModelInitial',@DesignFastForwardAnn);


%%
bottomUpMapping = BottomUpMappingFunction(systemFunction,systemParameter);

[fastForwardModel,problemData,iterData] = active_learning_model_training(...
    bottomUpMapping,...
    designSpaceLowerBoundary,...
    designSpaceUpperBoundary,...
    performanceLowerLimit,...
    performanceUpperLimit,...
    options);

methodDescription = sprintf('Training: Active Learning %d Samples',...
    size(fastForwardModel.DesignSampleTrain,1));


%% Model Information
% analytical good designs boundary
[sphereX,sphereY,sphereZ] = sphere;
innerSphereX = performanceLowerLimit*sphereX;
innerSphereY = performanceLowerLimit*sphereY;
innerSphereZ = performanceLowerLimit*sphereZ;
outerSphereX = performanceUpperLimit*sphereX;
outerSphereY = performanceUpperLimit*sphereY;
outerSphereZ = performanceUpperLimit*sphereZ;


%% Plot Data used for Training
designSampleTrain = fastForwardModel.DesignSampleTrain;
labelTrain = fastForwardModel.LabelTrain;
figure;
plot3(designSampleTrain(labelTrain,1),designSampleTrain(labelTrain,2),designSampleTrain(labelTrain,3),'g.','MarkerSize',10);
hold on;
grid minor;
plot3(designSampleTrain(~labelTrain,1),designSampleTrain(~labelTrain,2),designSampleTrain(~labelTrain,3),'r.','MarkerSize',10);
surf(innerSphereX,innerSphereY,innerSphereZ,'FaceAlpha',0.1,'FaceColor','black','EdgeColor','none');
surf(outerSphereX,outerSphereY,outerSphereZ,'FaceAlpha',0.1,'FaceColor','black','EdgeColor','none');
xlabel('x_1');
ylabel('x_2');
zlabel('x_3');
title(['Surrogate Modeling - ',methodDescription,'; Training Data']);
axis('equal','vis3d','square','tight');
save_print_figure(gcf,[saveFolder,'TrainingSamples']);


%% Test Classifiers
% Test Samples
testSample = sampling_random([designSpaceLowerBoundary;designSpaceUpperBoundary],100000);

% True Labels / Quantities of Interest
performanceMeasureTrue = bottomUpMapping.response(testSample);
performanceDeficitTrue = design_measure_to_deficit(performanceMeasureTrue,performanceLowerLimit,performanceUpperLimit);
[labelTrue,scoreTrue] = design_deficit_to_label_score(performanceDeficitTrue);

% Predicted Label - SVM
performanceDeficitPredicted = fastForwardModel.predict_designs(testSample);
[labelPredicted,scorePredicted] = design_deficit_to_label_score(performanceDeficitPredicted);

% Find Errors - False Positive and False Negatives
[falsePositive,falseNegative,~,~] = classification_confusion_matrix(labelTrue,labelPredicted);


%% approximate the threshold value of score where errors start happening
[scoreThreshold,iterDataUncertainty] = design_fast_forward_find_uncertainty_score(bottomUpMapping,fastForwardModel,...
    'NumberCandidateSamples',1000,'NumberNeighborEvaluations',3);
isUncertain = (abs(scorePredicted)<=scoreThreshold);


%% Plot False-Predictions
figure;
surf(innerSphereX,innerSphereY,innerSphereZ,'FaceAlpha',0.1,'FaceColor','black','EdgeColor','none','DisplayName','Boundary (Inner)');
hold on;
grid minor;
surf(outerSphereX,outerSphereY,outerSphereZ,'FaceAlpha',0.1,'FaceColor','black','EdgeColor','none','DisplayName','Boundary (Outer)');
plot3(testSample(falsePositive,1),testSample(falsePositive,2),testSample(falsePositive,3),'c.','DisplayName','False Positive');
plot3(testSample(falseNegative,1),testSample(falseNegative,2),testSample(falseNegative,3),'m.','DisplayName','False Negative');
plot3(testSample(isUncertain,1),testSample(isUncertain,2),testSample(isUncertain,3),'ko','DisplayName','Detected Uncertainty');
xlabel('x_1');
ylabel('x_2');
zlabel('x_3');
falsePredictionTitle = sprintf(...
    ['Surrogate Modeling - ',methodDescription,...
    '\nPrediction Error Rate: %g%%',...
    '\nUncertainty Rate: %g%%'],...
    100*sum(falsePositive|falseNegative)/size(falsePositive,1),...
    100*sum(isUncertain)/size(isUncertain,1));
title(falsePredictionTitle);
legend('location','southeast');
axis('equal','vis3d','square','tight');
save_print_figure(gcf,[saveFolder,'FalsePredictions']);


%% Metrics
algoData = postprocess_active_learning_model_training(problemData,iterData);
plot_active_learning_model_training_metrics(algoData,'SaveFolder',saveFolder);


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

