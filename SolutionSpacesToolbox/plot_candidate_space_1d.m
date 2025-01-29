function plotHandle = plot_candidate_space_1d(figureHandle, candidateSpace, varargin)
%PLOT_CANDIDATE_SPACE_1D Visualize the boundary of a 1D candidate space
%   PLOT_CANDIDATE_SPACE_1D plots the intervals of a 1D candidate space in the
%   given figure. It estimates where the boundary is via the scores applied to
%   a large sample to identify potentially disconnected intervals.
%
%   PLOT_CANDIDATE_SPACE_1D(FIGUREHANDLE,CANDIDATESPACE) plots in figure 
%   FIGUREHANDLE the intervals of the candidate space CANDIDATESPACE.
%
%   PLOT_CANDIDATE_SPACE_1D(...,NAME,VALUE,...) allows for setting
%   additional options for the plot operation; said options should refer to the
%   'line' function. The default properties are:
%       - 'Color': black ([0 0 0])
%       - 'LineWidth': 2
%       - 'YValue': 0 (vertical position of the line)
%   These can be changed, as well as any other 'line' option.
%
%   PLOTHANDLE = PLOT_CANDIDATE_SPACE_1D(...) returns the handle of the
%   generated plot PLOTHANDLE as a result of 'line', which can be used later 
%   for 'legend', for example.
%
%   Inputs:
%       - FIGUREHANDLE : Figure
%       - CANDIDATESPACE : CandidateSpaceBase
%       - Name-value pair arguments: passed directly to 'line'.
%
%   Output:
%       - PLOTHANDLE : Line
%
%   See also line, ClassificationSVM, legend.
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

    inputPlotOptions = parser_variable_input_to_structure(varargin{:});

    % set step size for fine sampling
    stepSize = 0.001;

    % make the intervals for each variable
    [~,positiveRegionBox] = design_bounding_box(candidateSpace.DesignSampleDefinition,candidateSpace.IsInsideDefinition);
    positiveRegionBox(1) = max(positiveRegionBox(1),candidateSpace.DesignSpaceLowerBound);
    positiveRegionBox(2) = min(positiveRegionBox(2),candidateSpace.DesignSpaceUpperBound);
    xInterval = positiveRegionBox(1) + (0:stepSize:1)*(positiveRegionBox(2)-positiveRegionBox(1));

    % get scores
    [~,score] = candidateSpace.is_in_candidate_space(xInterval');

    % find zero crossings to identify interval boundaries
    signChanges = diff(sign(score));
    boundaryPoints = xInterval(abs(signChanges) > 0);
    
    % Add endpoints if they are inside the space
    if score(1) < 0
        boundaryPoints = [xInterval(1), boundaryPoints];
    end
    if score(end) < 0
        boundaryPoints = [boundaryPoints, xInterval(end)];
    end

    % Reshape into pairs of points defining intervals
    nIntervals = length(boundaryPoints)/2;
    intervals = reshape(boundaryPoints, [2, nIntervals])';

    % plot options
    defaultPlotOptions = {'Color',[0 0 0],'LineWidth',2,'YValue',0};
    [~,plotOptions] = merge_name_value_pair_argument(defaultPlotOptions,inputPlotOptions);

    % Extract and remove YValue from options since it's not a line property
    yValueIndex = find(strcmpi(plotOptions, 'YValue'));
    if ~isempty(yValueIndex)
        yValue = plotOptions{yValueIndex + 1};
        plotOptions(yValueIndex:yValueIndex+1) = [];
    else
        yValue = 0;
    end

    % Plot each interval
    figure(figureHandle);
    hold on;
    plotHandle = gobjects(nIntervals,1);
    for i = 1:nIntervals
        xData = intervals(i,:);
        yData = [yValue yValue];
        plotHandle(i) = line(xData, yData, plotOptions{:});
    end

    if(nargout<1)
        clear plotHandle
    end
end 