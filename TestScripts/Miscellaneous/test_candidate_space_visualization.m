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


%%
designSpaceLowerBoundSample = [0 0];
designSpaceUpperBoundSample = [1 1];
designSpaceLowerBoundSpace = [-1 -1];
designSpaceUpperBoundSpace = [2 2];


%% visualize - convex hull
nTest = 20;
designSample = sampling_random([designSpaceLowerBoundSample; designSpaceUpperBoundSample],nTest);

nSampleInside = 8;
candidateSpace = CandidateSpaceConvexHull(designSpaceLowerBoundSpace,designSpaceUpperBoundSpace);
candidateSpace = candidateSpace.define_candidate_space(designSample(1:nSampleInside,:));
inCandidateSpace = candidateSpace.is_in_candidate_space(designSample);
candidateSpace = candidateSpace.define_candidate_space(designSample(inCandidateSpace,:));

figure;
hold all;
plot(designSample(inCandidateSpace,1),designSample(inCandidateSpace,2),'o','Color',color_palette_tol('green'),'MarkerSize',10,'linewidth',4);
plot(designSample(~inCandidateSpace,1),designSample(~inCandidateSpace,2),'x','Color',color_palette_tol('red'),'MarkerSize',10,'linewidth',4);
candidateSpace.plot_candidate_space(gcf,'EdgeColor',color_palette_tol('cyan'),'linewidth',4);
legend({'Inside Candidate Space','Outside Candidate Space','Convex Hull Boundary'},'location','south');
grid minor;
set(gca,'XColor', 'none','YColor','none');
h = get(gca,'Children');
set(gca,'Children',[h(3) h(2) h(1)]);
save_print_figure(gcf,[saveFolder,'CandidateSpaceConvexHullExample'],'Size',figureSize*1.25,'PrintFormat',{'png','pdf'});

% growth
growthRate = 0.1;

% find growth direction
insideSample = designSample(inCandidateSpace,:);
convexHullCenter = mean(insideSample,1);
distances = insideSample - convexHullCenter;
directionGrowth = distances./vecnorm(distances,2,2);
designSpaceFactor = designSpaceUpperBoundSpace - designSpaceLowerBoundSpace;
growthVector = growthRate.*designSpaceFactor.*directionGrowth;

insideSampleGrown = insideSample + growthVector;
candidateSpaceGrown = candidateSpace.grow_candidate_space(growthRate);
figure;
hold all;
plot(convexHullCenter(1),convexHullCenter(2),'b.','MarkerSize',20);
candidateSpace.plot_candidate_space(gcf,'EdgeColor','k','linewidth',4);
plot(insideSample(:,1),insideSample(:,2),'o','color',color_palette_tol('green'),'MarkerSize',10,'linewidth',4);
quiver(insideSample(:,1),insideSample(:,2),growthVector(:,1),growthVector(:,2),'AutoScale','off');
candidateSpaceGrown.plot_candidate_space(gcf,'EdgeColor',color_palette_tol('cyan'),'linewidth',4);
plot(insideSampleGrown(:,1),insideSampleGrown(:,2),'o','color',color_palette_tol('yellow'),'MarkerSize',10,'linewidth',4);
legend({'Convex Hull Center','Original Convex Hull','Inside Points','Growth Vector','Grown Convex Hull','New Inside Points'},'location','northwest')
grid minor;
set(gca,'XColor', 'none','YColor','none');
save_print_figure(gcf,[saveFolder,'CandidateSpaceConvexHullGrowth'],'Size',figureSize*1.25,'PrintFormat',{'png','pdf'});

% isShapeDefinition = candidateSpace.IsShapeDefinition;
% boundarySample = candidateSpace.DesignSampleDefinition(isShapeDefinition,:);
% boundaryGrowth = growthVector(isShapeDefinition,:);
% newBoundarySample = boundarySample + boundaryGrowth;
% candidateSpaceGrown = CandidateSpaceConvexHull(designSpaceLowerBoundSpace,designSpaceUpperBoundSpace);
% candidateSpaceGrown = candidateSpaceGrown.define_candidate_space(newBoundarySample);
% figure;
% hold all;
% plot(convexHullCenter(1),convexHullCenter(2),'b.','MarkerSize',20);
% candidateSpace.plot_candidate_space(gcf,'EdgeColor','k');
% plot(boundarySample(:,1),boundarySample(:,2),'g.','MarkerSize',20);
% quiver(boundarySample(:,1),boundarySample(:,2),boundaryGrowth(:,1),boundaryGrowth(:,2),'AutoScale','off');
% candidateSpaceGrown.plot_candidate_space(gcf,'EdgeColor','c');
% plot(newBoundarySample(:,1),newBoundarySample(:,2),'y.','MarkerSize',20);
% legend({'Convex Hull Center','Original Convex Hull','Convex Hull Boundary Points','Growth Vector','Grown Convex Hull','New Boundary Points'},'location','northwest')
% grid minor;
% set(gca,'XColor', 'none','YColor','none');
% save_print_figure(gcf,[saveFolder,'CandidateSpaceConvexHullGrowth'],'Size',figureSize*1.25,'PrintFormat',{'png','pdf'});


%% visualize - delaunay
nSampleInside = 11;

isInside = [true(nSampleInside,1);false(nTest-nSampleInside,1)];
candidateSpace = CandidateSpaceDelaunay(designSpaceLowerBoundSpace,designSpaceUpperBoundSpace);
candidateSpace = candidateSpace.define_candidate_space(designSample,isInside);

% find the simplices that were removed to plot those too for visualization
insideSample = designSample(1:nSampleInside,:);
delaunayIndex = delaunayn(insideSample);

% remove all simpleces with 'outside' designs inside them
outsideDelaunayIndex = [];
if(any(~isInside))
    % find which simplices have 'outside' designs inside
    outsideInsideSimplex = tsearchn(insideSample,delaunayIndex,designSample(~isInside,:));
    isInsideWrong = unique(outsideInsideSimplex(~isnan(outsideInsideSimplex)));

    % delete simplices where bad designs are inside
    outsideDelaunayIndex = delaunayIndex(isInsideWrong,:);
end

figure;
hold all;
plot(designSample(isInside,1),designSample(isInside,2),'g.','MarkerSize',20);
plot(designSample(~isInside,1),designSample(~isInside,2),'r.','MarkerSize',20);
candidateSpace.plot_candidate_space(gcf,'FaceColor','k','FaceAlpha',0.2);
patch('Faces',outsideDelaunayIndex,'Vertices',insideSample,'LineStyle','--','EdgeColor','m','FaceColor','m','FaceAlpha',0.2);
legend({'Inside Candidate Space','Outside Candidate Space','Simplices Inside Area','Removed Simplices'},'location','south');
grid minor;
set(gca,'XColor', 'none','YColor','none');
h = get(gca,'Children');
set(gca,'Children',[h(4) h(3) h(2) h(1)]);
save_print_figure(gcf,[saveFolder,'CandidateSpaceDelaunayExample'],'Size',figureSize*1.25,'PrintFormat',{'png','pdf'});

% growth 
growthRate = 0.1;

% for each simplex, get the growth vector
originalVertices = [];
growthVector = [];
for i=1:size(candidateSpace.DelaunaySimplex,2)
    simplexVertices = candidateSpace.DelaunaySimplex(i).Vertices;
    simplexCenter = mean(simplexVertices,1);
    distanceVertice = simplexVertices - simplexCenter;
    directionGrowth = distanceVertice./vecnorm(distanceVertice,2,2);

    originalVertices = [originalVertices;simplexVertices];
    growthVector = [growthVector;growthRate.*designSpaceFactor.*directionGrowth];
end

candidateSpaceGrown = candidateSpace.grow_candidate_space(growthRate);

figure;
hold all;
candidateSpace.plot_candidate_space(gcf,'FaceColor','k','FaceAlpha',0.2);
candidateSpaceGrown.plot_candidate_space(gcf,'FaceColor','c','FaceAlpha',0.2);
quiver(originalVertices(:,1),originalVertices(:,2),growthVector(:,1),growthVector(:,2),'AutoScale','off','Color',[0.9290 0.6940 0.1250]);
legend({'Original Simplices','Grown Simplices','Vertices Growth'},'location','northwest');
grid minor;
set(gca,'XColor', 'none','YColor','none');
save_print_figure(gcf,[saveFolder,'CandidateSpaceDelaunayGrowth'],'Size',figureSize*1.25,'PrintFormat',{'png','pdf'});


%% visualize - svm
nTest = 200;

% create separable region
designSample = designSpaceLowerBoundSample + rand(nTest,2).*(designSpaceUpperBoundSample-designSpaceLowerBoundSample);

systemFunction = @distance_to_center;
systemParameter = [0.5 0.5];
bottomUpMapping = BottomUpMappingFunction(systemFunction,systemParameter);

performanceLowerLimit = 0.15;
performanceUpperLimit = 0.4;
designEvaluator = DesignEvaluatorBottomUpMapping(bottomUpMapping,performanceLowerLimit,performanceUpperLimit);
labelSample = design_deficit_to_label_score(designEvaluator.evaluate(designSample));

% train candidate space
candidateSpace = CandidateSpaceSvm(designSpaceLowerBoundSpace,designSpaceUpperBoundSpace);
candidateSpace = candidateSpace.define_candidate_space(designSample,labelSample);
isShapeDefinition = candidateSpace.IsShapeDefinition;

figure;
plot(designSample(labelSample,1),designSample(labelSample,2),'g.','MarkerSize',20);
hold all;
plot(designSample(~labelSample,1),designSample(~labelSample,2),'r.','MarkerSize',20);
candidateSpace.plot_candidate_space(gcf,'FaceColor','k','FaceAlpha',0.2,'EdgeColor','none');
legend({'Inside Candidate Space','Outside Candidate Space','Positive Region of Classifier'});
grid minor;
set(gca,'XColor', 'none','YColor','none');
h = get(gca,'Children');
set(gca,'Children',[h(3) h(2) h(1)]);
save_print_figure(gcf,[saveFolder,'CandidateSpaceSvmExample'],'Size',figureSize*1.25,'PrintFormat',{'png','pdf'});

% growth
