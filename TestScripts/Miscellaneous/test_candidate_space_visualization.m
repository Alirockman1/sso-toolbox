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


%% color choice
colorCandidateSpaceInside = color_palette_tol('cyan');
colorCandidateSpaceOutside = color_palette_tol('purple');
colorPointInside = color_palette_tol('green');
colorPointOutside = color_palette_tol('red');
colorGrowthPoint = color_palette_tol('yellow');
colorVector = color_palette_tol('yellow');

optionsPointInside = {'linestyle','none','Marker','*','color',colorPointInside,'MarkerSize',10,'linewidth',2.0};
optionsPointOutside = {'linestyle','none','Marker','+','color',colorPointOutside,'MarkerSize',10,'linewidth',2.0};
optionsVector = {'AutoScale','off','Color',colorVector,'linewidth',1.5};


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
plot(designSample(inCandidateSpace,1),designSample(inCandidateSpace,2),optionsPointInside{:});
plot(designSample(~inCandidateSpace,1),designSample(~inCandidateSpace,2),optionsPointOutside{:});
candidateSpace.plot_candidate_space(gcf,'EdgeColor',colorCandidateSpaceInside,'linewidth',4);
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
candidateSpace.plot_candidate_space(gcf,'EdgeColor',colorCandidateSpaceInside,'LineStyle','--','linewidth',4);
plot(insideSample(:,1),insideSample(:,2),optionsPointInside{:});
quiver(insideSample(:,1),insideSample(:,2),growthVector(:,1),growthVector(:,2),optionsVector{:});
candidateSpaceGrown.plot_candidate_space(gcf,'EdgeColor',colorCandidateSpaceInside,'linewidth',4);
plot(insideSampleGrown(:,1),insideSampleGrown(:,2),optionsPointInside{:},'color',colorGrowthPoint);
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
nTest = 10;
designSample = sampling_random([designSpaceLowerBoundSample; designSpaceUpperBoundSample],nTest);

nSampleInside = 6;

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
plot(designSample(isInside,1),designSample(isInside,2),optionsPointInside{:});
plot(designSample(~isInside,1),designSample(~isInside,2),optionsPointOutside{:});
candidateSpace.plot_candidate_space(gcf,'FaceColor',colorCandidateSpaceInside,'FaceAlpha',0.5,'EdgeColor',colorCandidateSpaceInside,'FreeFacetOnly',false);
patch('Faces',outsideDelaunayIndex,'Vertices',insideSample,'LineStyle','--','EdgeColor','k','FaceColor','none','linewidth',1.5);
legend({'Inside Candidate Space','Outside Candidate Space','Simplices Inside Area','Removed Simplices'},'location','northwest');
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
newVertices = [];
for i=1:size(candidateSpace.DelaunaySimplex,2)
    simplexVertices = candidateSpace.DelaunaySimplex(i).Vertices;
    simplexCenter = mean(simplexVertices,1);
    distanceVertice = simplexVertices - simplexCenter;
    directionGrowth = distanceVertice./vecnorm(distanceVertice,2,2);

    originalVertices = [originalVertices;simplexVertices];
    growthVector = [growthVector;growthRate.*designSpaceFactor.*directionGrowth];

    insideSampleGrown = simplexVertices + directionGrowth.*growthRate.*designSpaceFactor;
    newVertices = [newVertices;insideSampleGrown];
end

candidateSpaceGrown = candidateSpace.grow_candidate_space(growthRate);


figure;
hold all;
plot(originalVertices(:,1),originalVertices(:,2),optionsPointInside{:});
plot(newVertices(:,1),newVertices(:,2),optionsPointInside{:},'color',colorGrowthPoint);
candidateSpace.plot_candidate_space(gcf,'FaceColor','none','EdgeColor',colorCandidateSpaceInside,'LineStyle','--','linewidth',2,'FreeFacetOnly',false);
candidateSpaceGrown.plot_candidate_space(gcf,'FaceColor','none','EdgeColor',colorCandidateSpaceInside,'linewidth',2);
quiver(originalVertices(:,1),originalVertices(:,2),growthVector(:,1),growthVector(:,2),optionsVector{:});
legend({'Original Vertices','Grown Vertices','Original Simplices','Grown Simplices','Vertices Growth'},'location','northwest');
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
plot(designSample(labelSample,1),designSample(labelSample,2),optionsPointInside{:});
hold all;
plot(designSample(~labelSample,1),designSample(~labelSample,2),optionsPointOutside{:});
candidateSpace.plot_candidate_space(gcf,'FaceColor',colorCandidateSpaceInside,'FaceAlpha',0.5,'EdgeColor','none');
legend({'Inside Candidate Space','Outside Candidate Space','Positive Region of Classifier'},'Location','southwest');
grid minor;
set(gca,'XColor', 'none','YColor','none');
h = get(gca,'Children');
set(gca,'Children',[h(3) h(2) h(1)]);
save_print_figure(gcf,[saveFolder,'CandidateSpaceSvmExample'],'Size',figureSize*1.25,'PrintFormat',{'png','pdf'});

% growth
