function plotHandle = plot_convex_hull_1d(figureHandle, hullSample, convexHullIndex, varargin)
%PLOT_CONVEX_HULL_1D Visualize a 1-dimensional convex hull
%   PLOT_CONVEX_HULL_1D plots in a figure the interval of a convex hull.
%   
%   PLOT_CONVEX_HULL_1D(FIGUREHANDLE,HULLSAMPLE,CONVEXHULLINDEX) plots in 
%   FIGUREHANDLE the convex hull defined by the samples used in its computation
%   HULLSAMPLE and the resulting connections CONVEXHULLINDEX. The hull is 
%   represented by a black line segment at y=0.
%
%   PLOT_CONVEX_HULL_1D(...,NAME,VALUE,...) also allows one to choose the 
%   options of the plot done for the convex hull; these options refer to the 
%   'line' function. The default properties are:
%       - 'Color': black ([0 0 0])
%       - 'LineWidth': 2
%       - 'YValue': 0 (vertical position of the line)
%   These can be changed, as well as any other 'line' option.
%
%   PLOTHANDLE = PLOT_CONVEX_HULL_1D(...) returns the handle of the resulting
%   plot, which can be used later for 'legend', for example.
%
%   Input:
%       - FIGUREHANDLE : Figure
%       - HULLSAMPLE : (nSample,1) double
%       - CONVEXHULLINDEX : (2,1) integer
%       - Name-value pair arguments : passed directly to 'line'.
%   Output:
%       - PLOTHANDLE : Line
%
%   See also line, compute_convex_hull, legend, plot_convex_hull_2d.
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

    % Get interval endpoints
    xData = [min(hullSample), max(hullSample)];
    yData = [yValue, yValue];

    figure(figureHandle);
    hold on;
    plotHandle = line(xData, yData, plotOptions{:});

    if(nargout<1)
        clear plotHandle
    end
end 