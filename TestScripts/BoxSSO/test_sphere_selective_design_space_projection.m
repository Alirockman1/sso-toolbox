%TEST_SPHERE_SELECTIVE_DESIGN_SPACE_PROJECTION solution space visualization
%   TEST_SPHERE_SELECTIVE_DESIGN_SPACE_PROJECTION uses a design problem where
%	all good designs are inside a ball, and generates selective design space
%	projection plots for this problem. This is used to illustrate the concept of
%	selective design space projection with an example.
%
%   Copyright 2025 Eduardo Rodrigues Della Noce
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
bottomUpMapping = BottomUpMappingFunction(@distance_to_center,[0,0,0;0,0,0;0,0,0]);
designEvaluator = DesignEvaluatorBottomUpMapping(bottomUpMapping,[-inf,-inf,3],[5,inf,inf]);

designSpaceLowerBound = [-6 -6 -6];
designSpaceUpperBound = [+6 +6 +6];
designBox = [-4 -3 -2; 2 3 4];

desiredPairs = [1 2; 1 3; 2 3];
plotGrid = [1,3];

h = plot_selective_design_space_projection(...
	designEvaluator,...
	designBox,...
	designSpaceLowerBound,...
	designSpaceUpperBound,...
	desiredPairs,...
	plotGrid,...
	'PlotIntervals',true,...
    'AxesLabels',{'x1','x2','x3'});
legend([h.GoodPerformance,h.BadPerformance(1),h.BadPerformance(3)],{'Good Designs',...
    'Outside Outer Radius','Inside Inner Radius'});
clear h

%% stop transcript 
save([saveFolder,'Data.mat']);
diary off;

