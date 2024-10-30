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


%% planar trimming
% 1d
figure;
hold all;
plot([designSpaceLowerBound1d anchorPosition1d],[-1 -1],'g-','linewidth',3.0);
plot([anchorPosition1d designSpaceUpperBound1d],[-1 -1],'r-','linewidth',3.0);
plot([anchorPosition1d],[-1],'r.','MarkerSize',20);
plot([designSpaceLowerBound1d anchorPosition1d],[1 1],'r-','linewidth',3.0,'HandleVisibility','off');
plot([anchorPosition1d designSpaceUpperBound1d],[1 1],'g-','linewidth',3.0,'HandleVisibility','off');
plot([anchorPosition1d],[1],'r.','MarkerSize',20,'HandleVisibility','off');
quiver(...
    [anchorPosition1d anchorPosition1d], [-1 1], ...
    [-0.1 0.1], [0 0], ...
    'AutoScale','off','LineWidth',2.0,'Color',[0.9290 0.6940 0.1250]);
set(gca,'XColor', 'none','YColor','none');
grid minor;
ylim([-2 2])
legend({'Kept Region','Removed Region','Removed Design','Normal Vector'},'location','east');
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
    convexHullInside = convexHullInside.define_candidate_space(pointInside);

    convexHullOutside = CandidateSpaceConvexHull(designSpaceLowerBound2d,designSpaceUpperBound2d);
    convexHullOutside = convexHullOutside.define_candidate_space(pointOutside);

    convexHullOutside.plot_candidate_space(gcf,'FaceColor','r','FaceAlpha',0.2,'LineStyle','none');
    convexHullInside.plot_candidate_space(gcf,'FaceColor','g','FaceAlpha',0.2,'LineStyle','none');
    quiver(...
        anchorPosition2d(1),anchorPosition2d(2),...
        trimmingNormalsInside(i,1),trimmingNormalsInside(i,2),...
        'AutoScale','on','AutoScaleFactor',0.3,'LineWidth',2.0,'Color',[0.9290 0.6940 0.1250]);
    plot(anchorPosition2d(1),anchorPosition2d(2),'r.','MarkerSize',20);
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

    convexHullInside = CandidateSpaceConvexHull(designSpaceLowerBound3d,designSpaceUpperBound3d);
    convexHullInside = convexHullInside.define_candidate_space(pointInside);

    convexHullOutside = CandidateSpaceConvexHull(designSpaceLowerBound3d,designSpaceUpperBound3d);
    convexHullOutside = convexHullOutside.define_candidate_space(pointOutside);

    convexHullOutside.plot_candidate_space(gcf,'FaceColor','r','FaceAlpha',0.2,'LineStyle','none');
    convexHullInside.plot_candidate_space(gcf,'FaceColor','g','FaceAlpha',0.2,'LineStyle','none');
    quiver3(...
        anchorPosition3d(1),anchorPosition3d(2),anchorPosition3d(3),...
        trimmingNormalsInside(i,1),trimmingNormalsInside(i,2),trimmingNormalsInside(i,3),...
        'AutoScale','on','AutoScaleFactor',0.3,'LineWidth',2.0,'Color',[0.9290 0.6940 0.1250]);
    plot3(anchorPosition3d(1),anchorPosition3d(2),anchorPosition3d(3),'r.','MarkerSize',20);
    set(gca,'XColor','none','YColor','none','ZColor','none');
    grid minor;
    view(3);
end
legend({'Removed Region','Kept Region','Normal Vector','Removed Design'});
save_print_figure(gcf,[saveFolder,'PlanarTrimming3D'],'Size',[figureSize(1) 2*figureSize(2)],'PrintFormat',{'png','pdf'});



%% corner box removal
% 1d
figure;
hold all;
plot([designSpaceLowerBound1d anchorPosition1d],[-1 -1],'g-','linewidth',3.0);
plot([anchorPosition1d designSpaceUpperBound1d],[-1 -1],'r-','linewidth',3.0);
plot([anchorPosition1d],[-1],'r.','MarkerSize',20);
plot([designSpaceLowerBound1d anchorPosition1d],[1 1],'r-','linewidth',3.0);
plot([anchorPosition1d designSpaceUpperBound1d],[1 1],'g-','linewidth',3.0);
plot([anchorPosition1d],[1],'r.','MarkerSize',20);
set(gca,'XColor', 'none','YColor','none');
grid minor;
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
            handleRemoved = plot_design_box_2d(gcf,cornerBox,'FaceColor','r','FaceAlpha',0.2,'LineStyle','none');
        else
            handleKept = plot_design_box_2d(gcf,cornerBox,'FaceColor','g','FaceAlpha',0.2,'LineStyle','none');
        end
    end
    handleAnchor = plot(anchorPosition2d(1),anchorPosition2d(2),'r.','MarkerSize',20);
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
            handleRemoved = plot_design_box_3d(gcf,cornerBox,'FaceColor','r','FaceAlpha',0.2,'LineStyle','none');
        else
            handleKept = plot_design_box_3d(gcf,cornerBox,'FaceColor','g','FaceAlpha',0.2,'LineStyle','none');
        end
    end
    handleAnchor = plot3(anchorPosition3d(1),anchorPosition3d(2),anchorPosition3d(3),'r.','MarkerSize',20);
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
    handleKept = plot_design_box_3d(gcf,goodRegion,'FaceColor','g','FaceAlpha',0.2,'LineStyle','none');

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
    handleRemoved = plot_design_box_3d(gcf,cornerBox,'FaceColor','r','FaceAlpha',0.2,'LineStyle','none');

    handleAnchor = plot3(anchorPosition3d(1),anchorPosition3d(2),anchorPosition3d(3),'r.','MarkerSize',20);
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
plot([designSpaceLowerBound1d designSpaceUpperBound1d],[0 0],'g-','linewidth',3.0);
plot([anchorPosition1d-holeSize anchorPosition1d+holeSize],[0 0],'r-','linewidth',3.0);
plot(anchorPosition1d,0,'r.','MarkerSize',20);
set(gca,'XColor', 'none','YColor','none');
grid minor;
legend({'Kept Region','Removed Region','Removed Design'},'location','northeast');
save_print_figure(gcf,[saveFolder,'HolePunching1D'],'Size',[figureSize(1) figureSize(2)/2],'PrintFormat',{'png','pdf'});

% 2d
figure;
hold all;
plot_design_box_2d(gcf,designSpace2d,'FaceColor','g','FaceAlpha',0.2,'LineStyle','none');
plot_design_box_2d(gcf,anchorPosition2d+[-holeSize -holeSize ;+holeSize +holeSize],'FaceColor','r','FaceAlpha',0.2,'LineStyle','none');
plot(anchorPosition2d(1),anchorPosition2d(2),'r.','MarkerSize',20);
set(gca,'XColor', 'none','YColor','none');
grid minor;
legend({'Kept Region','Removed Region','Removed Design'},'location','east');
save_print_figure(gcf,[saveFolder,'HolePunching2D'],'Size',[figureSize(1) figureSize(2)],'PrintFormat',{'png','pdf'});

% 3d
figure;
hold all;
plot_design_box_3d(gcf,designSpace3d,'FaceColor','g','FaceAlpha',0.2,'LineStyle','none');
plot_design_box_3d(gcf,anchorPosition3d+[-holeSize -holeSize -holeSize;+holeSize +holeSize +holeSize],'FaceColor','r','FaceAlpha',0.2,'LineStyle','none');
plot3(anchorPosition3d(1),anchorPosition3d(2),anchorPosition3d(3),'r.','MarkerSize',20);
set(gca,'XColor','none','YColor','none','ZColor','none');
grid minor;
view(3)
legend({'Kept Region','Removed Region','Removed Design'},'location','east');
save_print_figure(gcf,[saveFolder,'HolePunching3D'],'Size',[figureSize(1) figureSize(2)],'PrintFormat',{'png','pdf'});

