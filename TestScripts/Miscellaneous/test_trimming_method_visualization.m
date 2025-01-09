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
designSpaceLowerBound1d = 0;
designSpaceUpperBound1d = 1;
designSpace1d = [designSpaceLowerBound1d;designSpaceUpperBound1d];

designSpaceLowerBound2d = [0 0];
designSpaceUpperBound2d = [1 1];
designSpace2d = [designSpaceLowerBound2d;designSpaceUpperBound2d];

designSpaceLowerBound3d = [0 0 0];
designSpaceUpperBound3d = [1 1 1];
designSpace3d = [designSpaceLowerBound3d;designSpaceUpperBound3d];

anchorPosition1d = 0.3;
anchorPosition2d = [0.3 0.3];
anchorPosition3d = [0.3 0.3 0.3];


%% color choice
colorCandidateSpaceInside = color_palette_tol('cyan');
colorCandidateSpaceOutside = color_palette_tol('grey');
colorPointInside = color_palette_tol('green');
colorPointOutside = color_palette_tol('red');
colorGrowthPoint = color_palette_tol('yellow');
colorVector = color_palette_tol('purple');

optionsCandidateSpaceInside1d = {'-','color',colorCandidateSpaceInside,'linewidth',3.0};
optionsCandidateSpaceOutside1d = {'-','color',colorCandidateSpaceOutside,'linewidth',3.0};
optionsCandidateSpaceInside2d3d = {'FaceColor',colorCandidateSpaceInside,'FaceAlpha',0.4,'LineStyle','none'};
optionsCandidateSpaceOutside2d3d = {'FaceColor',colorCandidateSpaceOutside,'FaceAlpha',0.6,'LineStyle','none'};
optionsPointOutside = {'x','color',colorPointOutside,'MarkerSize',10,'linewidth',2.0};
optionsPointInside = {'o','color',colorPointInside,'MarkerSize',10,'linewidth',2.0};


%% planar trimming - "pure"
% 1d
figure;
hold all;
plot([anchorPosition1d designSpaceUpperBound1d],[-1 -1],optionsCandidateSpaceOutside1d{:});
plot([designSpaceLowerBound1d anchorPosition1d],[-1 -1],optionsCandidateSpaceInside1d{:});
plot([anchorPosition1d],[-1],optionsPointOutside{:});
plot([designSpaceLowerBound1d anchorPosition1d],[1 1],optionsCandidateSpaceOutside1d{:},'HandleVisibility','off');
plot([anchorPosition1d designSpaceUpperBound1d],[1 1],optionsCandidateSpaceInside1d{:},'HandleVisibility','off');
quiver(...
    [anchorPosition1d anchorPosition1d], [-1 1], ...
    [-0.15 0.15], [0 0], ...
    'AutoScale','off','LineWidth',2.0,'Color',colorVector);
plot([anchorPosition1d],[1],optionsPointOutside{:},'HandleVisibility','off');
set(gca,'XColor', 'none','YColor','none');
grid off;
ylim([-2 2]);
legend({'Removed Region','Kept Region','Normal Vector','Removed Design'},'location','east');
save_print_figure(gcf,[saveFolder,'PlanarTrimming1D'],'Size',[figureSize(1) figureSize(2)/2],'PrintFormat',{'png','pdf'});


% 2d
trimmingNormalsInside = [1 1; -3 3; 0.5 -0.7; 0 1; 1 0; -1 -1];

trimmingNormalsInside = trimmingNormalsInside./vecnorm(trimmingNormalsInside,2,2);

cornerCombination = logical(round(fullfact(2*ones(2,1)) - 1));
cornerCoordinate = repmat(designSpaceLowerBound2d,size(cornerCombination,1),1);
cornerCoordinateUpper = repmat(designSpaceUpperBound2d,size(cornerCombination,1),1);
cornerCoordinate(cornerCombination) = cornerCoordinateUpper(cornerCombination);
distanceCorner = anchorPosition2d - cornerCoordinate;

figure;
for i=1:size(trimmingNormalsInside,1)
    subplot(3,2,i);
    hold all;

    planeVector = [trimmingNormalsInside(i,2), -trimmingNormalsInside(i,1)];

    designSpaceLimit = region_limit_line_search([],[anchorPosition2d;anchorPosition2d],[planeVector;-planeVector],designSpace2d);
    designSpaceIntersection = anchorPosition2d + designSpaceLimit.*[planeVector;-planeVector];

    dotProduct = sum(trimmingNormalsInside(i,:).*distanceCorner,2);

    pointInside = [anchorPosition2d;designSpaceIntersection;cornerCoordinate(dotProduct<=0,:)];
    pointOutside = [anchorPosition2d;designSpaceIntersection;cornerCoordinate(dotProduct>0,:)];

    convexHullInside = CandidateSpaceConvexHull(designSpaceLowerBound2d,designSpaceUpperBound2d);
    convexHullInside = convexHullInside.generate_candidate_space(pointInside);

    convexHullOutside = CandidateSpaceConvexHull(designSpaceLowerBound2d,designSpaceUpperBound2d);
    convexHullOutside = convexHullOutside.generate_candidate_space(pointOutside);

    convexHullOutside.plot_candidate_space(gcf,optionsCandidateSpaceOutside2d3d{:});
    convexHullInside.plot_candidate_space(gcf,optionsCandidateSpaceInside2d3d{:});
    quiver(...
        anchorPosition2d(1),anchorPosition2d(2),...
        trimmingNormalsInside(i,1),trimmingNormalsInside(i,2),...
        'AutoScale','on','AutoScaleFactor',0.3,'LineWidth',2.0,'Color',colorVector);
    plot(anchorPosition2d(1),anchorPosition2d(2),optionsPointOutside{:});
    set(gca,'XColor', 'none','YColor','none');
    grid minor;
end
legend({'Removed Region','Kept Region','Normal Vector','Removed Design'});
save_print_figure(gcf,[saveFolder,'PlanarTrimming2D'],'Size',[figureSize(1) 2*figureSize(2)],'PrintFormat',{'png','pdf'});

% 3d
trimmingNormalsInside = [1 1 1;  -1 -1 -1; -3 3 -3; 0.5 -0.7 -.2; 0 1 0; 1 0 0];

trimmingNormalsInside = trimmingNormalsInside./vecnorm(trimmingNormalsInside,2,2);

cornerCombination = logical(round(fullfact(2*ones(3,1)) - 1));
cornerCoordinate = repmat(designSpaceLowerBound3d,size(cornerCombination,1),1);
cornerCoordinateUpper = repmat(designSpaceUpperBound3d,size(cornerCombination,1),1);
cornerCoordinate(cornerCombination) = cornerCoordinateUpper(cornerCombination);
distanceCorner = anchorPosition3d - cornerCoordinate;

figure;
for i=1:size(trimmingNormalsInside,1)
    subplot(3,2,i);
    hold all;

    % plane vectors
    planeVectorNull = null(trimmingNormalsInside(i,:))';
    intersectionInterpolation = linspace(0,1,100)';
    planeVector = [...
        intersectionInterpolation.*planeVectorNull(1,:) + (1-intersectionInterpolation).*planeVectorNull(2,:);...
        intersectionInterpolation.*(-planeVectorNull(1,:)) + (1-intersectionInterpolation).*planeVectorNull(2,:);...
        intersectionInterpolation.*(-planeVectorNull(1,:)) + (1-intersectionInterpolation).*(-planeVectorNull(2,:));...
        intersectionInterpolation.*planeVectorNull(1,:) + (1-intersectionInterpolation).*(-planeVectorNull(2,:))];

    designSpaceLimit = region_limit_line_search([],...
        repmat(anchorPosition3d,size(planeVector,1),1),...
        planeVector,...
        designSpace3d);
    designSpaceIntersection = anchorPosition3d + designSpaceLimit.*planeVector;

    dotProduct = sum(trimmingNormalsInside(i,:).*distanceCorner,2);

    pointInside = [anchorPosition3d;designSpaceIntersection;cornerCoordinate(dotProduct<=0,:)];
    pointOutside = [anchorPosition3d;designSpaceIntersection;cornerCoordinate(dotProduct>0,:)];

    tolerance = 1e-5;

    convexHullInside = CandidateSpaceConvexHull(designSpaceLowerBound3d-tolerance,designSpaceUpperBound3d+tolerance);
    convexHullInside = convexHullInside.generate_candidate_space(pointInside+tolerance*trimmingNormalsInside(i,:));

    convexHullOutside = CandidateSpaceConvexHull(designSpaceLowerBound3d-tolerance,designSpaceUpperBound3d+tolerance);
    convexHullOutside = convexHullOutside.generate_candidate_space(pointOutside-tolerance*trimmingNormalsInside(i,:));

    convexHullOutside.plot_candidate_space(gcf,optionsCandidateSpaceOutside2d3d{:},'linestyle','none');
    convexHullInside.plot_candidate_space(gcf,optionsCandidateSpaceInside2d3d{:});
    quiver3(...
        anchorPosition3d(1),anchorPosition3d(2),anchorPosition3d(3),...
        trimmingNormalsInside(i,1),trimmingNormalsInside(i,2),trimmingNormalsInside(i,3),...
        'AutoScale','on','AutoScaleFactor',0.3,'LineWidth',2.0,'Color',colorVector);
    plot3(anchorPosition3d(1),anchorPosition3d(2),anchorPosition3d(3),optionsPointOutside{:});
    set(gca,'XColor','none','YColor','none','ZColor','none');
    grid minor;
    view(3);
end
legend({'Removed Region','Kept Region','Normal Vector','Removed Design'});
save_print_figure(gcf,[saveFolder,'PlanarTrimming3D'],'Size',[figureSize(1) 2*figureSize(2)],'PrintFormat',{'png','pdf'});


%% planar trimming - with designs / preferential direction
designSpaceLowerBound = 0;
designSpaceUpperBound = +1;

designSample = [0.12 0.22 0.35 0.44 0.52 0.56 0.78 0.91 0.99]';
isGood = [true false false true true true false true true]';
trimmingOrder = find(~isGood);
candidateSpace = CandidateSpacePlanarTrimming(designSpaceLowerBound,designSpaceUpperBound);
trimmedCandidateSpace = component_trimming_operation(designSample,isGood,trimmingOrder,{[1]},candidateSpace,'TrimmingMethodOptions',{'TrimmingSlack',1});

anchorPoint = trimmedCandidateSpace.AnchorPoint(1,:);
planeOrientation = trimmedCandidateSpace.PlaneOrientationAnchor(1,:);

figure;
hold all;
if(planeOrientation<0)
    plot([anchorPoint designSpaceUpperBound],[0 0],optionsCandidateSpaceOutside1d{:});
    plot([designSpaceLowerBound anchorPoint],[0 0],optionsCandidateSpaceInside1d{:});
else
    plot([designSpaceLowerBound anchorPoint],[0 0],optionsCandidateSpaceOutside1d{:});
    plot([anchorPoint designSpaceUpperBound],[0 0],optionsCandidateSpaceInside1d{:},'HandleVisibility','off');
end
plot(designSample(isGood),zeros(sum(isGood),1),optionsPointInside{:});
plot(designSample(~isGood),zeros(sum(~isGood),1),optionsPointOutside{:});
plot(anchorPoint,0,optionsPointOutside{:},'MarkerSize',12,'linewidth',3.0);
daspect([13,3,1]);
pbaspect([13,3,1]);
arrow3([anchorPoint,0], [anchorPoint+planeOrientation*0.1 0], 'f', 2, 2);
set(gca,'XColor', 'none','YColor','none');
grid off;
% ylim([-1 1]);
save_print_figure(gcf,[saveFolder,'PlanarTrimming1D-Preferential'],'Size',[13 3],'PrintFormat',{'png','pdf'});

clear designSpaceLowerBound designSpaceUpperBound designSample isGood trimmingOrder
clear candidateSpace trimmedCandidateSpace anchorPoint planeOrientation

% 2d
designSpaceLowerBound = [0,0];
designSpaceUpperBound = [+1,+1];
designSpace = [designSpaceLowerBound;designSpaceUpperBound];

designSample = [0.45,0.87; 0.91,0.54; 0.36,0.69; 0.31,0.42; 0.44,0.22; 0.62,0.31; 0.81,0.23; 0.12,0.61; 0.84,0.77];
isGood = [true true false true true true false true false]';
trimmingOrder = find(~isGood);
candidateSpace = CandidateSpacePlanarTrimming(designSpaceLowerBound,designSpaceUpperBound);
trimmedCandidateSpace = component_trimming_operation(designSample,isGood,trimmingOrder,{[1,2]},candidateSpace,'TrimmingMethodOptions',{'TrimmingSlack',1});

anchorPoint = trimmedCandidateSpace.AnchorPoint(1,:);
planeOrientation = trimmedCandidateSpace.PlaneOrientationAnchor(1,:);

cornerCombination = logical(round(fullfact(2*ones(2,1)) - 1));
cornerCoordinate = repmat(designSpaceLowerBound,size(cornerCombination,1),1);
cornerCoordinateUpper = repmat(designSpaceUpperBound,size(cornerCombination,1),1);
cornerCoordinate(cornerCombination) = cornerCoordinateUpper(cornerCombination);
distanceCorner = anchorPoint - cornerCoordinate;

planeVector = [planeOrientation(2), -planeOrientation(1)];

designSpaceLimit = region_limit_line_search([],[anchorPoint;anchorPoint],[planeVector;-planeVector],designSpace);
designSpaceIntersection = anchorPoint + designSpaceLimit.*[planeVector;-planeVector];

dotProduct = sum(planeOrientation.*distanceCorner,2);

pointInside = [anchorPoint;designSpaceIntersection;cornerCoordinate(dotProduct<=0,:)];
pointOutside = [anchorPoint;designSpaceIntersection;cornerCoordinate(dotProduct>0,:)];

convexHullInside = CandidateSpaceConvexHull(designSpaceLowerBound,designSpaceUpperBound);
convexHullInside = convexHullInside.generate_candidate_space(pointInside);

convexHullOutside = CandidateSpaceConvexHull(designSpaceLowerBound,designSpaceUpperBound);
convexHullOutside = convexHullOutside.generate_candidate_space(pointOutside);

figure;
convexHullOutside.plot_candidate_space(gcf,optionsCandidateSpaceOutside2d3d{:});
convexHullInside.plot_candidate_space(gcf,optionsCandidateSpaceInside2d3d{:});
plot(designSample(isGood,1),designSample(isGood,2),optionsPointInside{:});
plot(designSample(~isGood,1),designSample(~isGood,2),optionsPointOutside{:});
plot(anchorPoint(1),anchorPoint(2),optionsPointOutside{:},'MarkerSize',12,'linewidth',3.0);
daspect([1,1,1]);
pbaspect([1,1,1]);
arrow3([anchorPoint(1),anchorPoint(2)], [anchorPoint(1)+0.2*planeOrientation(1) anchorPoint(2)+0.2*planeOrientation(2)], 'f', 2, 2);
set(gca,'XColor', 'none','YColor','none');
grid off;
% ylim([-1 1]);
save_print_figure(gcf,[saveFolder,'PlanarTrimming2D-Preferential'],'Size',[10 10],'PrintFormat',{'png','pdf'});

clear designSpaceLowerBound designSpaceUpperBound designSample isGood trimmingOrder
clear candidateSpace trimmedCandidateSpace anchorPoint planeOrientation
clear cornerCombination cornerCoordinate cornerCoordinateUpper distanceCorner
clear planeVector designSpaceLimit designSpaceIntersection dotProduct
clear pointInside pointOutside convexHullInside convexHullOutside


% 3d
designSpaceLowerBound = [0,0,0];
designSpaceUpperBound = [+1,+1,+1];
designSpace = [designSpaceLowerBound;designSpaceUpperBound];

designSample = sampling_latin_hypercube(designSpace,9);
isGood = [true true false true true true false true false]';
trimmingOrder = find(~isGood);
candidateSpace = CandidateSpacePlanarTrimming(designSpaceLowerBound,designSpaceUpperBound);
trimmedCandidateSpace = component_trimming_operation(designSample,isGood,trimmingOrder,{[1,2,3]},candidateSpace,'TrimmingMethodOptions',{'TrimmingSlack',1});

anchorPoint = trimmedCandidateSpace.AnchorPoint(1,:);
planeOrientation = trimmedCandidateSpace.PlaneOrientationAnchor(1,:);

cornerCombination = logical(round(fullfact(2*ones(3,1)) - 1));
cornerCoordinate = repmat(designSpaceLowerBound,size(cornerCombination,1),1);
cornerCoordinateUpper = repmat(designSpaceUpperBound,size(cornerCombination,1),1);
cornerCoordinate(cornerCombination) = cornerCoordinateUpper(cornerCombination);
distanceCorner = anchorPoint - cornerCoordinate;

% plane vectors
planeVectorNull = null(planeOrientation)';
intersectionInterpolation = linspace(0,1,100)';
planeVector = [...
    intersectionInterpolation.*planeVectorNull(1,:) + (1-intersectionInterpolation).*planeVectorNull(2,:);...
    intersectionInterpolation.*(-planeVectorNull(1,:)) + (1-intersectionInterpolation).*planeVectorNull(2,:);...
    intersectionInterpolation.*(-planeVectorNull(1,:)) + (1-intersectionInterpolation).*(-planeVectorNull(2,:));...
    intersectionInterpolation.*planeVectorNull(1,:) + (1-intersectionInterpolation).*(-planeVectorNull(2,:))];

designSpaceLimit = region_limit_line_search([],...
    repmat(anchorPoint,size(planeVector,1),1),...
    planeVector,...
    designSpace);
designSpaceIntersection = anchorPoint + designSpaceLimit.*planeVector;

dotProduct = sum(planeOrientation.*distanceCorner,2);

pointInside = [anchorPoint;designSpaceIntersection;cornerCoordinate(dotProduct<=0,:)];
pointOutside = [anchorPoint;designSpaceIntersection;cornerCoordinate(dotProduct>0,:)];

tolerance = 1e-5;

convexHullInside = CandidateSpaceConvexHull(designSpaceLowerBound-tolerance,designSpaceUpperBound+tolerance);
convexHullInside = convexHullInside.generate_candidate_space(pointInside+tolerance*planeOrientation);

convexHullOutside = CandidateSpaceConvexHull(designSpaceLowerBound-tolerance,designSpaceUpperBound+tolerance);
convexHullOutside = convexHullOutside.generate_candidate_space(pointOutside-tolerance*planeOrientation);

figure;
convexHullOutside.plot_candidate_space(gcf,optionsCandidateSpaceOutside2d3d{:},'linestyle','none');
convexHullInside.plot_candidate_space(gcf,optionsCandidateSpaceInside2d3d{:});
plot3(designSample(isGood,1),designSample(isGood,2),designSample(isGood,3),optionsPointInside{:});
plot3(designSample(~isGood,1),designSample(~isGood,2),designSample(~isGood,3),optionsPointOutside{:});
plot3(anchorPoint(1),anchorPoint(2),anchorPoint(3),optionsPointOutside{:},'MarkerSize',12,'linewidth',3.0);
daspect([1,1,1]);
pbaspect([1,1,1]);
arrow3([anchorPoint(1),anchorPoint(2),anchorPoint(3)], [anchorPoint(1)+0.3*planeOrientation(1),anchorPoint(2)+0.3*planeOrientation(2),anchorPoint(3)+0.3*planeOrientation(3)], 'f', 2, 2);
set(gca,'XColor', 'none','YColor','none','ZColor','none');
grid off;
view(3);
save_print_figure(gcf,[saveFolder,'PlanarTrimming3D-Preferential'],'Size',[10 10],'PrintFormat',{'png','pdf'});
 
clear designSpaceLowerBound designSpaceUpperBound designSample isGood trimmingOrder
clear candidateSpace trimmedCandidateSpace anchorPoint planeOrientation
clear cornerCombination cornerCoordinate cornerCoordinateUpper distanceCorner
clear planeVector designSpaceLimit designSpaceIntersection dotProduct
clear pointInside pointOutside convexHullInside convexHullOutside


%% corner box removal
% 1d
figure;
hold all;
plot([designSpaceLowerBound1d anchorPosition1d],[-1 -1],optionsCandidateSpaceInside1d{:});
plot([anchorPosition1d designSpaceUpperBound1d],[-1 -1],optionsCandidateSpaceOutside1d{:});
plot([anchorPosition1d],[-1],optionsPointOutside{:});
plot([designSpaceLowerBound1d anchorPosition1d],[1 1],optionsCandidateSpaceOutside1d{:});
plot([anchorPosition1d designSpaceUpperBound1d],[1 1],optionsCandidateSpaceInside1d{:});
plot([anchorPosition1d],[1],optionsPointOutside{:});
set(gca,'XColor', 'none','YColor','none');
grid off;
ylim([-2 2])
legend({'Kept Region','Removed Region','Removed Design'},'location','east');
save_print_figure(gcf,[saveFolder,'CornerBoxRemoval1D'],'Size',[figureSize(1) figureSize(2)/2],'PrintFormat',{'png','pdf'});

% 2d
cornerCombination = [...
    false true;... % left top corner
    false false;... % left bottom corner
    true true;... % right top corner
    true false]; % right bottom corner
figure;
for i=1:size(cornerCombination,1)
    subplot(2,2,i);
    hold all;

    for j=1:size(cornerCombination,1)
        cornerBox = designSpace2d;
        
        if(cornerCombination(j,1))
            cornerBox(2,2) = anchorPosition2d(1);
        else
            cornerBox(1,2) = anchorPosition2d(1);
        end

        if(cornerCombination(j,2))
            cornerBox(2,1) = anchorPosition2d(2);
        else
            cornerBox(1,1) = anchorPosition2d(2);
        end

        if(j==i)
            handleRemoved = plot_design_box_2d(gcf,cornerBox,optionsCandidateSpaceOutside2d3d{:});
        else
            handleKept = plot_design_box_2d(gcf,cornerBox,optionsCandidateSpaceInside2d3d{:});
        end
    end
    handleAnchor = plot(anchorPosition2d(1),anchorPosition2d(2),optionsPointOutside{:});
    set(gca,'XColor', 'none','YColor','none');
    grid minor;
end
legend([handleRemoved,handleKept,handleAnchor],{'Removed Region','Kept Region','Removed Design'});
save_print_figure(gcf,[saveFolder,'CornerBoxRemoval2D'],'Size',[figureSize(1) figureSize(2)],'PrintFormat',{'png','pdf'});

% 3d - version 1
cornerCombination = logical(round(fullfact(2*ones(3,1)) - 1));
tileMapping = [1 3 5 7 2 4 6 8];
figure;
for i=1:size(cornerCombination,1)
    subplot(4,2,tileMapping(i));
    hold all;

    for j=1:size(cornerCombination,1)
        cornerBox = designSpace3d;
        
        if(cornerCombination(j,1))
            cornerBox(2,1) = anchorPosition3d(1);
        else
            cornerBox(1,1) = anchorPosition3d(1);
        end

        if(cornerCombination(j,2))
            cornerBox(2,2) = anchorPosition3d(2);
        else
            cornerBox(1,2) = anchorPosition3d(2);
        end

        if(cornerCombination(j,3))
            cornerBox(2,3) = anchorPosition3d(3);
        else
            cornerBox(1,3) = anchorPosition3d(3);
        end

        if(j==i)
            handleRemoved = plot_design_box_3d(gcf,cornerBox,optionsCandidateSpaceOutside2d3d{:});
        else
            handleKept = plot_design_box_3d(gcf,cornerBox,optionsCandidateSpaceInside2d3d{:});
        end
    end
    handleAnchor = plot3(anchorPosition3d(1),anchorPosition3d(2),anchorPosition3d(3),optionsPointOutside{:});
    set(gca,'XColor', 'none','YColor','none','ZColor','none');
    grid minor;
    view(3);
end
legend([handleRemoved,handleKept,handleAnchor],{'Removed Region','Kept Region','Removed Design'});
save_print_figure(gcf,[saveFolder,'CornerBoxRemoval3D-1'],'Size',[figureSize(1) 2*figureSize(2)],'PrintFormat',{'png','pdf'});

% 3d - version 2
cornerCombination = logical(round(fullfact(2*ones(3,1)) - 1));
tileMapping = [1 3 5 7 2 4 6 8];
figure;
for i=1:size(cornerCombination,1)
    subplot(4,2,tileMapping(i));
    hold all;
    tolerance = -1e-3;
    goodRegion = designSpace3d + [+tolerance;-tolerance];
    handleKept = plot_design_box_3d(gcf,goodRegion,optionsCandidateSpaceInside2d3d{:});

    cornerBox = designSpace3d;
    if(cornerCombination(i,1))
        cornerBox(2,1) = anchorPosition3d(1);
    else
        cornerBox(1,1) = anchorPosition3d(1);
    end
    if(cornerCombination(i,2))
        cornerBox(2,2) = anchorPosition3d(2);
    else
        cornerBox(1,2) = anchorPosition3d(2);
    end
    if(cornerCombination(i,3))
        cornerBox(2,3) = anchorPosition3d(3);
    else
        cornerBox(1,3) = anchorPosition3d(3);
    end
    handleRemoved = plot_design_box_3d(gcf,cornerBox,optionsCandidateSpaceOutside2d3d{:});

    handleAnchor = plot3(anchorPosition3d(1),anchorPosition3d(2),anchorPosition3d(3),optionsPointOutside{:});
    set(gca,'XColor', 'none','YColor','none','ZColor','none');
    grid minor;
    view(3);
end
legend([handleRemoved,handleKept,handleAnchor],{'Removed Region','Kept Region','Removed Design'});
save_print_figure(gcf,[saveFolder,'CornerBoxRemoval3D-2'],'Size',[figureSize(1) 2*figureSize(2)],'PrintFormat',{'png','pdf'});


%% hole punching
holeSize = 0.1;

% 1d
figure;
hold all;
plot([designSpaceLowerBound1d anchorPosition1d-holeSize],[0 0],optionsCandidateSpaceInside1d{:});
plot([anchorPosition1d+holeSize designSpaceUpperBound1d],[0 0],optionsCandidateSpaceInside1d{:},'HandleVisibility','off');
plot([anchorPosition1d-holeSize anchorPosition1d+holeSize],[0 0],optionsCandidateSpaceOutside1d{:});
plot(anchorPosition1d,0,optionsPointOutside{:});
set(gca,'XColor', 'none','YColor','none');
grid off;
legend({'Kept Region','Removed Region','Removed Design'},'location','northeast');
save_print_figure(gcf,[saveFolder,'HolePunching1D'],'Size',[figureSize(1) figureSize(2)/2],'PrintFormat',{'png','pdf'});

% 2d
figure;
hold all;
plot_design_box_2d(gcf,designSpace2d,optionsCandidateSpaceInside2d3d{:});
plot_design_box_2d(gcf,anchorPosition2d+[-holeSize -holeSize ;+holeSize +holeSize],optionsCandidateSpaceOutside2d3d{:});
plot(anchorPosition2d(1),anchorPosition2d(2),optionsPointOutside{:});
set(gca,'XColor', 'none','YColor','none');
grid minor;
legend({'Kept Region','Removed Region','Removed Design'},'location','east');
save_print_figure(gcf,[saveFolder,'HolePunching2D'],'Size',[figureSize(1) figureSize(2)],'PrintFormat',{'png','pdf'});

% 3d
figure;
hold all;
plot_design_box_3d(gcf,designSpace3d,optionsCandidateSpaceInside2d3d{:});
plot_design_box_3d(gcf,anchorPosition3d+[-holeSize -holeSize -holeSize;+holeSize +holeSize +holeSize],optionsCandidateSpaceOutside2d3d{:});
plot3(anchorPosition3d(1),anchorPosition3d(2),anchorPosition3d(3),optionsPointOutside{:});
set(gca,'XColor','none','YColor','none','ZColor','none');
grid minor;
view(3)
legend({'Kept Region','Removed Region','Removed Design'},'location','east');
save_print_figure(gcf,[saveFolder,'HolePunching3D'],'Size',[figureSize(1) figureSize(2)],'PrintFormat',{'png','pdf'});


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

