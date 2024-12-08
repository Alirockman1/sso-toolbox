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
colorVector = color_palette_tol('yellow');

optionsCandidateSpaceInside1d = {'-','color',colorCandidateSpaceInside,'linewidth',3.0};
optionsCandidateSpaceOutside1d = {'-','color',colorCandidateSpaceOutside,'linewidth',3.0};
optionsCandidateSpaceInside2d3d = {'FaceColor',colorCandidateSpaceInside,'FaceAlpha',0.4,'LineStyle','none'};
optionsCandidateSpaceOutside2d3d = {'FaceColor',colorCandidateSpaceOutside,'FaceAlpha',0.6,'LineStyle','none'};
optionsPointOutside = {'x','color',colorPointOutside,'MarkerSize',10,'linewidth',2.0};


%% planar trimming
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
    convexHullInside = convexHullInside.define_candidate_space(pointInside);

    convexHullOutside = CandidateSpaceConvexHull(designSpaceLowerBound2d,designSpaceUpperBound2d);
    convexHullOutside = convexHullOutside.define_candidate_space(pointOutside);

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
    convexHullInside = convexHullInside.define_candidate_space(pointInside+tolerance*trimmingNormalsInside(i,:));

    convexHullOutside = CandidateSpaceConvexHull(designSpaceLowerBound3d-tolerance,designSpaceUpperBound3d+tolerance);
    convexHullOutside = convexHullOutside.define_candidate_space(pointOutside-tolerance*trimmingNormalsInside(i,:));

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

