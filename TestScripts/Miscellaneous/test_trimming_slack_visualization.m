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
colorPointInside = color_palette_tol('green');
colorPointOutside = color_palette_tol('red');
optionsPointOutside = {'x','color',colorPointOutside,'MarkerSize',10,'linewidth',2.0};
optionsPointInside = {'o','color',colorPointInside,'MarkerSize',10,'linewidth',2.0};


%% CHOSEN SAMPLES
designSpaceLowerBound = [0 0];
designSpaceUpperBound = [1 1];
anchorPoint = [0.21,0.32];
planeOrientation = [1,1.5];
designSample = [0.31,0.61;0.47,0.35;0.18,0.75;0.66,0.33;0.62,0.64];
designPadding = sampling_latin_hypercube([designSpaceLowerBound;designSpaceUpperBound],100);


%% PLANAR TRIMMING
% normalize plane orientation
planeOrientation = planeOrientation./vecnorm(planeOrientation,2,2);

% find maximum slack
dotProduct = sum(planeOrientation.*(anchorPoint - designSample),2);
allowedSlack = min(-dotProduct);
anchorSlack = allowedSlack*planeOrientation;

anchorPointHalf = anchorPoint + (1-0.5)*anchorSlack;
anchorPointNull = anchorPoint + (1-0)*anchorSlack;

designSpaceIntersectionOne = find_design_space_intersection_plane(designSpaceLowerBound,designSpaceUpperBound,anchorPoint,planeOrientation);
designSpaceIntersectionHalf = find_design_space_intersection_plane(designSpaceLowerBound,designSpaceUpperBound,anchorPointHalf,planeOrientation);
designSpaceIntersectionNull = find_design_space_intersection_plane(designSpaceLowerBound,designSpaceUpperBound,anchorPointNull,planeOrientation);

figure;
hold all;
plot(designSpaceIntersectionOne(:,1),designSpaceIntersectionOne(:,2),'Color',color_palette_tol('grey'),'linewidth',1.5);
plot(designSpaceIntersectionHalf(:,1),designSpaceIntersectionHalf(:,2),'Color',color_palette_tol('purple'),'linewidth',1.5);
plot(designSpaceIntersectionNull(:,1),designSpaceIntersectionNull(:,2),'Color',color_palette_tol('yellow'),'linewidth',1.5);
plot(anchorPoint(1),anchorPoint(2),optionsPointOutside{:});
% plot(anchorPointHalf(1),anchorPointHalf(2),optionsPointOutside{:},'linewidth',1.0,'HandleVisibility','off');
% plot(anchorPointNull(1),anchorPointNull(2),optionsPointOutside{:},'linewidth',1.0,'HandleVisibility','off');
plot(designSample(:,1),designSample(:,2),optionsPointInside{:});
plot(designPadding(:,1),designPadding(:,2),'.','Color',color_palette_tol('blue'),'MarkerSize',10);
axis([0.1 0.7 0.1 0.7])
daspect([1,1,1]);
pbaspect([1,1,1]);
arrow3([anchorPoint(1),anchorPoint(2)], [anchorPoint(1)+0.05*planeOrientation(1) anchorPoint(2)+0.05*planeOrientation(2)], 'f', 0.5, 0.5);
grid off;
set(gca,'XColor', 'none','YColor','none');
legend({'Plane at Slack=1','Plane at Slack=0.5','Plane at Slack=0','Removed Design','Remaining Good Designs','Padding Points','Plane Orientation'});
save_print_figure(gcf,[saveFolder,'PlanarTrimming-Slack'],'PrintFormat',{'png','pdf'});


%% CORNER BOX REMOVAL
combination = [false,false];

% find maximum slack
lowerBound = min(designSample,[],1);
upperBound = max(designSample,[],1);
normalizationFactor = upperBound-lowerBound;

distanceToAnchorBase = designSample-anchorPoint;
distanceToAnchor = distanceToAnchorBase;
distanceToAnchor(:,combination) = -distanceToAnchor(:,combination);
insideToAnchorDistance = distanceToAnchor./normalizationFactor;

nDesignVariable = size(designSample,2);
maximumSlack = nan(1,nDesignVariable);
[maximumSlackInside,iDimensionMaximumSlack] = max(insideToAnchorDistance,[],2);
for j=1:nDesignVariable
    allowedSlack = maximumSlackInside(iDimensionMaximumSlack==j).*normalizationFactor(j);
    
    if(isempty(allowedSlack))
        if(~combination(j))
            maximumSlack(j) = upperBound(j) - anchorPoint(j);
        else
            maximumSlack(j) = anchorPoint(j) - lowerBound(j);
        end
    else
        maximumSlack(j) = min(allowedSlack);
    end
end
anchorSlack(~combination) = maximumSlack(~combination);
anchorSlack(combination) = -maximumSlack(combination);

anchorPointFull = anchorPoint + (1-1)*anchorSlack;
anchorPointHalf = anchorPoint + (1-0.5)*anchorSlack;
anchorPointNull = anchorPoint + (1-0)*anchorSlack;

boxRemovalFull = find_design_space_intersection_corner(designSpaceLowerBound,designSpaceUpperBound,anchorPointFull,combination);
boxRemovalHalf = find_design_space_intersection_corner(designSpaceLowerBound,designSpaceUpperBound,anchorPointHalf,combination);
boxRemovalNull = find_design_space_intersection_corner(designSpaceLowerBound,designSpaceUpperBound,anchorPointNull,combination);

figure;
hold all;
plot_design_box_2d(gcf,boxRemovalFull,'EdgeColor',color_palette_tol('grey'),'FaceColor','none','linewidth',1.5);
plot_design_box_2d(gcf,boxRemovalHalf,'EdgeColor',color_palette_tol('purple'),'FaceColor','none','linewidth',1.5);
plot_design_box_2d(gcf,boxRemovalNull,'EdgeColor',color_palette_tol('yellow'),'FaceColor','none','linewidth',1.5);
plot(anchorPoint(1),anchorPoint(2),optionsPointOutside{:});
% plot(anchorPointHalf(1),anchorPointHalf(2),optionsPointOutside{:},'linewidth',1.0,'HandleVisibility','off');
% plot(anchorPointNull(1),anchorPointNull(2),optionsPointOutside{:},'linewidth',1.0,'HandleVisibility','off');
plot(designSample(:,1),designSample(:,2),optionsPointInside{:});
axis([0.1 0.8 0.1 0.8]);
% daspect([1,1,1]);
% pbaspect([1,1,1]);
grid off;
set(gca,'XColor', 'none','YColor','none');
legend({'Corner Removed at Slack=1','Corner Removed at Slack=0.5','Corner Removed at Slack=0','Removed Design','Remaining Good Designs'},'location','southeast');
save_print_figure(gcf,[saveFolder,'CornerBoxRemoval-Slack'],'PrintFormat',{'png','pdf'});


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;


%%
function designSpaceIntersection = find_design_space_intersection_plane(designSpaceLowerBound,designSpaceUpperBound,anchorPoint,planeOrientation)
    planeVector = [planeOrientation(2), -planeOrientation(1)];
    
    designSpace = [designSpaceLowerBound;designSpaceUpperBound];
    designSpaceLimit = region_limit_line_search([],[anchorPoint;anchorPoint],[planeVector;-planeVector],designSpace);
    designSpaceIntersection = anchorPoint + designSpaceLimit.*[planeVector;-planeVector];
end

function boxRemoval = find_design_space_intersection_corner(designSpaceLowerBound,designSpaceUpperBound,anchorPoint,combination)
    boxRemoval = [designSpaceLowerBound;designSpaceUpperBound];
    
    nDimension = size(designSpaceLowerBound,2);
    for i=1:nDimension
        if(combination(i))
            boxRemoval(1,i) = anchorPoint(i);
        else
            boxRemoval(2,i) = anchorPoint(i);
        end
    end
end