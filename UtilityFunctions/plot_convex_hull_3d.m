function plotHandle = plot_convex_hull_3d(figureHandle,hullSample,convexHullIndex,varargin)
%PLOT_CONVEX_HULL_3D Visualize a 3-dimensional convex hull
%   PLOT_CONVEX_HULL_3D plots in a figure the surfaces of a convex hull.
%   
%   PLOT_CONVEX_HULL_3D(FIGUREHANDLE,HULLSAMPLE,CONVEXHULLINDEX) plots in 
%   FIGUREHANDLE the convex hull defined by the samples used in its computation
%   HULLSAMPLE and the resulting connections CONVEXHULLINDEX. The hull is 
%   represented by black faces.
%
%   PLOT_CONVEX_HULL_3D(...,NAME,VALUE,...) also allows one to choose the 
%   options of the plot done for the convex hull; these options refer to the 
%   'patch' function. The one property with a default value is 'FaceColor', 
%   which is set to black ([0 0 0]); this can be changed, as well as any other 
%   'patch' option.
%
%   PLOTHANDLE = PLOT_CONVEX_HULL_3D(...) returns the handle of the resulting
%   plot from 'patch', which can be used later for 'legend', for example.
%
%   Input:
%       - FIGUREHANDLE : Figure
%       - HULLSAMPLE : (nSample,2) double
%       - CONVEXHULLINDEX : (nHull,2) integer
%       - Name-value pair arguments: passed directly to 'patch'.
%   Output:
%       - PLOTHANDLE : Line
%
%   See also patch, compute_convex_hull, legend, plot_convex_hull_2d.
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
    
    defaultPlotOptions = {'FaceColor',[0 0 0]};
    [~,plotOptions] = merge_name_value_pair_argument(defaultPlotOptions,inputPlotOptions);

    figure(figureHandle);
    hold on;
	plotHandle = patch('Faces',convexHullIndex,'Vertices',hullSample,...
        plotOptions{:});

    if(nargout<1)
        clear plotHandle
    end
end