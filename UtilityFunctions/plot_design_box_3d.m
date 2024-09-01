function h = plot_design_box_3d(figureHandle,designBox,varargin)
%PLOT_DESIGN_BOX_3D Visualize a 3-dimensional design box
%   PLOT_DESIGN_BOX_3D plots in a figure the edges of a design box.
%
%   PLOT_DESIGN_BOX_3D(FIGUREHANDLE,DESIGNBOX) plots in figure FIGUREHANDLE
%   the design box DESIGNBOX. Said box is shown with its edges in black.
%
%   PLOT_DESIGN_BOX_3D(...NAME,VALUE,...) also allows for the choice of options 
%   regarding how the design box is displayed; said options refer to the 'patch'
%   function.
%   
%   PLOTHANDLE = PLOT_DESIGN_BOX_3D(...) returns the handle of the plot created
%   by 'patch', which can be used later for 'legend', for example.
%
%   For more information on how 'patch' works, see:
%   https://mathworks.com/help/matlab/visualize/multifaceted-patches.html
%
%   Input:
%       - FIGUREHANDLE : Figure
%       - DESIGNBOX : (2,3) double
%           -- Row (1) : lower boundary of the design box
%           -- Row (2) : upper boundary of the design box
%       - OPTIONS : 'patch' options, structure OR name-value pair cell, optional
%   Output:
%       - PLOTHANDLE : Line
%
%   See also patch, legend, plot_design_box_2d.
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

    % separate box into an anchor corner and its edges
    boxCorner = min(designBox,[],1);
    boxEdges = (designBox(2,:)-designBox(1,:));

    % find box vertices
    %               x y z
    verticesBase = [0 0 0;
                    1 0 0;
                    1 1 0;
                    0 1 0;
                    0 0 1;
                    1 0 1;
                    1 1 1;
                    0 1 1];            
   faces = [1 2 6 5;
            2 3 7 6;
            3 4 8 7;
            4 1 5 8;
            1 2 3 4;
            5 6 7 8];
   vertices = boxCorner + boxEdges.*verticesBase;
    
    % set options
    inputPlotOptions = parser_variable_input_to_structure(varargin{:});
    defaultPlotOptions = {};
    [~,plotOptions] = merge_name_value_pair_argument(defaultPlotOptions,inputPlotOptions);

    % make plot
    figure(figureHandle);
    h = patch('Vertices', vertices, 'Faces', faces, plotOptions{:});
    view(3);

    if(nargout<1)
        clear h
    end
end