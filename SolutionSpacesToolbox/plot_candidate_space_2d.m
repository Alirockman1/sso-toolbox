function plotHandle = plot_candidate_space_2d(graphicsHandle,candidateSpace,varargin)
%PLOT_CANDIDATE_SPACE_2D Visualize the boundary of a candidate space
%   PLOT_CANDIDATE_SPACE_2D plots the boundary of a candidate space in the
%   given figure. It estimates where the boundary is via the scores applied to
%   a large sample, and using contour to plot where that score is closest to 0.
%
%   PLOT_CANDIDATE_SPACE_2D(GRAPHICSHANDLE,CANDIDATESPACE) plots in 
%   GRAPHICSHANDLE the boundary of the candidate space CANDIDATESPACE.
%
%   PLOT_CANDIDATE_SPACE_2D(...,NAME,VALUE,...) allows for setting
%   additional options for the plot operation; said options should refer to the
%   'contour' function.
%
%   PLOTHANDLE = PLOT_CANDIDATE_SPACE_2D(...) returns the handle of the
%   generated plot PLOTHANDLE as a result of 'contour', which can be used later 
%   for 'legend', for example.
%
%   Inputs:
%       - GRAPHICSHANDLE : Figure
%       - CANDIDATESPACE : CandidateSpaceBase
%       - Name-value pair arguments: passed directly to 'contour'.
%
%   Output:
%       - PLOTHANDLE : Line
%
%   See also contour, ClassificationSVM, legend.
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

    parser = inputParser;
    parser.KeepUnmatched = true;
    parser.addParameter('StepSize',0.001);
    parser.addParameter('FixedVariables',[]);
    parser.parse(varargin{:});

    stepSize = parser.Results.StepSize;
    fixedVariables = parser.Results.FixedVariables;
    inputPlotOptions = namedargs2cell(parser.Unmatched);

    activate_graphics_object(graphicsHandle);
    hold on;

    % make the intervals for each variable
    [~,positiveRegionBox] = design_bounding_box(candidateSpace.DesignSampleDefinition,candidateSpace.IsInsideDefinition);
    positiveRegionBox(1,:) = max(positiveRegionBox(1,:),candidateSpace.DesignSpaceLowerBound);
    positiveRegionBox(2,:) = min(positiveRegionBox(2,:),candidateSpace.DesignSpaceUpperBound);

    if(~isempty(fixedVariables))
        nDimension = size(positiveRegionBox,2);
        isFixed = ~isnan(fixedVariables);
        positiveRegionBox = positiveRegionBox(:,~isFixed);
    end

    xInterval = positiveRegionBox(1,1) + (0:stepSize:1)*(positiveRegionBox(2,1)-positiveRegionBox(1,1));
    yInterval = positiveRegionBox(1,2) + (0:stepSize:1)*(positiveRegionBox(2,2)-positiveRegionBox(1,2));

    % generate grid for predictions at finer sample rate
    [xGrid, yGrid] = meshgrid(xInterval,yInterval);
    fullGrid = [xGrid(:),yGrid(:)];

    if(~isempty(fixedVariables))
        nPoints = size(fullGrid,1);
        fullGridComplete = nan(nPoints,nDimension);
        fullGridComplete(:,~isFixed) = fullGrid;
        fullGridComplete(:,isFixed) = repmat(fixedVariables(isFixed),nPoints,1);
        fullGrid = fullGridComplete;
        clear fullGridComplete;
    end

    % get scores
    [~,score] = candidateSpace.is_in_candidate_space(fullGrid);

    % reshape to same grid size as the input
    score = reshape(score(:,1), size(xGrid));

    % plot options
    defaultPlotOptions = {};
    [~,plotOptions] = merge_name_value_pair_argument(defaultPlotOptions,inputPlotOptions);

    % find decision boundary surface
    contourMatrix = contourc(xInterval,yInterval,score,[0 0]);

    % process data
    countourXY = [];
    iEntry = 0;
    startCurve = [];
    endCurve = [];
    while(~isempty(contourMatrix))
        iEntry = iEntry + 1;

        index = contourMatrix(2)+1;
        countourXY{iEntry} = [contourMatrix(:,2:index)'];
        contourMatrix(:,1:index) = [];
        
        startCurve = [startCurve;countourXY{iEntry}(1,:)];
        endCurve = [endCurve;countourXY{iEntry}(end,:)];
    end

    % determine the order of entries
    orderedCurve = countourXY{1};
    availableQuery = true(iEntry,1);
    availableQuery(1) = false;
    for i=1:iEntry-1
        % find start of next curve closest to end of current curve
        [iClosestAvailableForward,distanceForward] = knnsearch(startCurve(availableQuery,:),orderedCurve(end,:));
        [iClosestAvailableBackward,distanceBackward] = knnsearch(endCurve(availableQuery,:),orderedCurve(end,:));
        
        if(distanceForward<distanceBackward)
            iClosest = convert_index_base(availableQuery,iClosestAvailableForward,'backward');
            xyCoordinate = countourXY{iClosest};
        else
            iClosest = convert_index_base(availableQuery,iClosestAvailableBackward,'backward');
            xyCoordinate = flip(countourXY{iClosest},1);
        end
        orderedCurve = [orderedCurve;xyCoordinate];
        availableQuery(iClosest) = false;
    end

    % plot options
    defaultPlotOptions = {};
    [~,plotOptions] = merge_name_value_pair_argument(defaultPlotOptions,inputPlotOptions);

    activate_graphics_object(graphicsHandle);
    hold on;
    plotHandle = patch('XData',orderedCurve(:,1),'YData',orderedCurve(:,2),plotOptions{:});

    if(nargout<1)
        clear plotHandle
    end
 end