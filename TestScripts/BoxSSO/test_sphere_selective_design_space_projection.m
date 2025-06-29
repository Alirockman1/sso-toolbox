%% cleanup
close all;
fclose all;
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
defaultFigureSize = [goldenRatio 1]*8.5;


%% 
bottomUpMapping = BottomUpMappingFunction(@distance_to_center,[0,0,0]);
designEvaluator = DesignEvaluatorBottomUpMapping(bottomUpMapping,-inf,5);

designSpaceLowerBound = [-6 -6 -6];
designSpaceUpperBound = [+6 +6 +6];
designBox = [-4 -3 -2; 2 3 4];

desiredPairs = [1 2; 1 3; 2 3];
plotGrid = [1,3];

plot_selective_design_space_projection(designEvaluator,designBox,designSpaceLowerBound,designSpaceUpperBound,desiredPairs,plotGrid,'PlotIntervals',true);


%% stop transcript 
save([saveFolder,'Data.mat']);
diary off;

