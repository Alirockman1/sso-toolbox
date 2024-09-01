function plotHandle = plot_design_box_2d(figureHandle,designBox,varargin)
%PLOT_DESIGN_BOX_2D Visualize a 2-dimensional design box
%	PLOT_DESIGN_BOX_2D plots in a figure the edges of a design box.
%
%	PLOT_DESIGN_BOX_2D(FIGUREHANDLE,DESIGNBOX) plots in figure FIGUREHANDLE
%	the design box DESIGNBOX. Said box is shown with its edges in black.
%
%	PLOT_DESIGN_BOX_2D(...NAME,VALUE,...) also allows for the choice of options 
%	regarding how the design box is displayed; said options refer to the 'patch'
%	function. Per default, the only option changed is 'FaceAlpha', which assumes
%	a value of 0 (transparent). This and all other 'patch' options can be 
%	changed as desired.
%	
%	PLOTHANDLE = PLOT_DESIGN_BOX_2D(...) returns the handle of the plot created
%	by 'patch', which can be used later for 'legend', for example.
%
%   Input:
%       - FIGUREHANDLE : Figure
%		- DESIGNBOX : (2,2) double
%           -- Row (1) : lower boundary of the design box
%           -- Row (2) : upper boundary of the design box
%		- Name-value pair arguments : passed directly to 'patch'.
%
%	Output:
%       - PLOTHANDLE : Line
%
%   See also patch, legend, plot_design_box_3d.
%
%   Copyright 2024 Eduardo Rodrigues Della Noce
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

	% separate box into an anchor corner and its edges
	boxCorner = min(designBox,[],1);
	boxEdges = (designBox(2,:)-designBox(1,:));

	% find box vertices
	%               x y
	verticesBase = [0 0;
					1 0;
					1 1;
					0 1];
	faces = [1 2 3 4];
	vertices = boxCorner + boxEdges.*verticesBase;

	% set options
	defaultPlotOptions = {'FaceAlpha',0};
	[~,plotOptions] = merge_name_value_pair_argument(defaultPlotOptions,inputPlotOptions);
	
	% make plot
	figure(figureHandle);
	plotHandle = patch('Vertices', vertices,'Faces', faces, plotOptions{:});

	if(nargout<1)
        clear plotHandle
    end
end