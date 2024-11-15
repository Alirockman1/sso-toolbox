%TUTORIAL_03_SURROGATE_MODELING Surrogate model of the hollow sphere problem
%   TUTORIAL_03_SURROGATE_MODELING creates a surrogate model for the 
%   hollow sphere problem. This is done using an active set strategy.

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
% definition: designSpaceLowerBound,designSpaceUpperBound

%% Active Learning 
%TODO: train a model using active learning
% definition: fastForwardModel,problemData,iterData


%% Model Information
methodDescription = sprintf('Training: Active Learning %d Samples',...
    size(fastForwardModel.DesignSampleTrain,1));

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

% Plot False-Predictions
figure;
surf(innerSphereX,innerSphereY,innerSphereZ,'FaceAlpha',0.1,'FaceColor','black','EdgeColor','none','DisplayName','Boundary (Inner)');
hold on;
grid minor;
surf(outerSphereX,outerSphereY,outerSphereZ,'FaceAlpha',0.1,'FaceColor','black','EdgeColor','none','DisplayName','Boundary (Outer)');
plot3(testSample(falsePositive,1),testSample(falsePositive,2),testSample(falsePositive,3),'c.','DisplayName','False Positive');
plot3(testSample(falseNegative,1),testSample(falseNegative,2),testSample(falseNegative,3),'m.','DisplayName','False Negative');
xlabel('x_1');
ylabel('x_2');
zlabel('x_3');
falsePredictionTitle = sprintf(...
    ['Surrogate Modeling - ',methodDescription,...
    '\nPrediction Error Rate: %g%%'],...
    100*sum(falsePositive|falseNegative)/size(falsePositive,1));
title(falsePredictionTitle);
legend('location','southeast');
save_print_figure(gcf,[saveFolder,'FalsePredictions']);


%% Performance Metrics
%TODO : visualize evolution of performance metrics


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

