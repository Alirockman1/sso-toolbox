%% Cleanup
close all;
fclose all;
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


%% generate sample - 2D  
designSampleGood =  rand(40,2);
designSampleBad = rand(40,2);
designSample = [designSampleGood;designSampleBad];
shrinkFactor = 1.0;

[boundaryIndex,measure] = boundary(designSampleGood,shrinkFactor);

[inArea,onEdge] = inpolygon(designSample(:,1),designSample(:,2),designSample(boundaryIndex,1),designSample(boundaryIndex,2));
inside = inArea | onEdge;

figure;
plot(designSample(:,1),designSample(:,2),'.');
hold on;
grid minor;
plot(designSample(boundaryIndex,1),designSample(boundaryIndex,2));
plot(designSample(inside,1),designSample(inside,2),'m*');


%% generate sample - 3D  
designSampleGood =  rand(40,3);
designSampleBad = rand(40,3);
designSample = [designSampleGood;designSampleBad];
shrinkFactor = 1.0;

[boundaryIndex,measure] = boundary(designSampleGood,shrinkFactor);

figure;
plot3(designSample(:,1),designSample(:,2),designSample(:,3),'.');
hold on;
grid minor;
trisurf(boundaryIndex,designSample(:,1),designSample(:,2),designSample(:,3),'FaceColor','g','FaceAlpha',0.2);

[boundaryIndex,measure] = boundary(designSampleGood,0);

figure;
plot3(designSample(:,1),designSample(:,2),designSample(:,3),'.');
hold on;
grid minor;
trisurf(boundaryIndex,designSample(:,1),designSample(:,2),designSample(:,3),'FaceColor','g','FaceAlpha',0.2);

shape = alphaShape(designSampleGood);
alphaRadius = criticalAlpha(shape,'all-points');
possibleAlpha = alphaSpectrum(shape);
boundaryIndex = shape.boundaryFacets;

figure;
plot3(designSample(:,1),designSample(:,2),designSample(:,3),'.');
hold on;
grid minor;
trisurf(boundaryIndex,designSample(:,1),designSample(:,2),designSample(:,3),'FaceColor','g','FaceAlpha',0.2);

shape = alphaShape(designSampleGood,0.5);
boundaryIndex = shape.boundaryFacets;

figure;
plot3(designSample(:,1),designSample(:,2),designSample(:,3),'.');
hold on;
grid minor;
trisurf(boundaryIndex,designSample(:,1),designSample(:,2),designSample(:,3),'FaceColor','g','FaceAlpha',0.2);


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

