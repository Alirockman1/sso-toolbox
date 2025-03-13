function plotHandle = plot_design_box_1d(graphicsHandle, designBox, varargin)
%PLOT_DESIGN_BOX_1D Visualize a 1-dimensional design box
%   PLOT_DESIGN_BOX_1D plots in a figure the interval of a design box.
%
%   PLOT_DESIGN_BOX_1D(FIGUREHANDLE,DESIGNBOX) plots in figure FIGUREHANDLE
%   the design box DESIGNBOX. Said box is shown as a black line segment at y=0.
%
%   PLOT_DESIGN_BOX_1D(...NAME,VALUE,...) also allows for the choice of options 
%   regarding how the design box is displayed; said options refer to the 'line'
%   function. The default properties are:
%       - 'Color': black ([0 0 0])
%       - 'LineWidth': 2
%       - 'YValue': 0 (vertical position of the line)
%   These can be changed, as well as any other 'line' option.
%   
%   PLOTHANDLE = PLOT_DESIGN_BOX_1D(...) returns the handle of the plot created
%   by 'line', which can be used later for 'legend', for example.
%
%   Input:
%       - FIGUREHANDLE : Figure
%       - DESIGNBOX : (2,1) double
%           -- Row (1) : lower boundary of the design box
%           -- Row (2) : upper boundary of the design box
%       - Name-value pair arguments : passed directly to 'line'.
%
%   Output:
%       - PLOTHANDLE : Line
%
%   See also line, legend, plot_design_box_2d, plot_design_box_3d.
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

    % set options
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

    % make plot
    activate_graphics_object(graphicsHandle);
    hold on;
    plotHandle = line([designBox(1) designBox(2)], [yValue yValue], plotOptions{:});

    if(nargout<1)
        clear plotHandle
    end
end 