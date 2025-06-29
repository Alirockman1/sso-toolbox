function plotHandle = plot_candidate_space_3d(fighan,candidateSpace,varargin)
%PLOT_CANDIDATE_SPACE_3D Visualize the boundary of a candidate space
%   PLOT_CANDIDATE_SPACE_3D plots the boundary of a candidate space in the
%   given figure. It estimates where the boundary is via the scores applied to
%   a large sample, and using isosurface and patch to plot where that score is 
%   closest to 0.
%
%   PLOT_CANDIDATE_SPACE_3D(FIGUREHANDLE,CANDIDATESPACE) plots in figure 
%   FIGUREHANDLE the boundary of the candidate space CANDIDATESPACE.
%
%   PLOT_CANDIDATE_SPACE_3D(...,NAME,VALUE,...) allows for setting
%   additional options for the plot operation; said options should refer to the
%   'patch' function.
%
%   PLOTHANDLE = PLOT_CANDIDATE_SPACE_3D(...) returns the handle of the
%   generated plot PLOTHANDLE as a result of 'patch', which can be used later 
%   for 'legend', for example.
%
%   Inputs:
%       - FIGUREHANDLE : Figure
%       - CANDIDATESPACE : CandidateSpaceBase
%       - Name-value pair arguments: passed directly to 'patch'.
%
%   Output:
%       - PLOTHANDLE : Line
%
%   See also patch, ClassificationSVM, legend, isosurface.
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

    figure(fighan);
    hold on

    % set step size for fine sampling
    stepSize = 0.01;

    % make the intervals for each variable
    definitionSample = candidateSpace.DesignSampleDefinition;
    definitionLabel = candidateSpace.IsInsideDefinition;
    [~,positiveRegionBox] = design_bounding_box(definitionSample,definitionLabel);
    xInterval = positiveRegionBox(1,1) + (0:stepSize:1)*(positiveRegionBox(2,1)-positiveRegionBox(1,1));
    yInterval = positiveRegionBox(1,2) + (0:stepSize:1)*(positiveRegionBox(2,2)-positiveRegionBox(1,2));
    zInterval = positiveRegionBox(1,3) + (0:stepSize:1)*(positiveRegionBox(2,3)-positiveRegionBox(1,3));

    % generate grid for predictions at finer sample rate
    [xGrid, yGrid, zGrid] = meshgrid(xInterval,yInterval,zInterval);
    fullGrid = [xGrid(:),yGrid(:),zGrid(:)];

    % get scores
    [~,score] = candidateSpace.is_in_candidate_space(fullGrid);

    % reshape to same grid size as the input
    score = reshape(score(:,1), size(xGrid));

    % plot options
    defaultPlotOptions = {};
    [~,plotOptions] = merge_name_value_pair_argument(defaultPlotOptions,inputPlotOptions);

    % plot decision boundary surface
    [faces,vertices] = isosurface(xGrid, yGrid, zGrid, score, 0);
    plotHandle = patch('Vertices', vertices, 'Faces', faces, plotOptions{:});
    box on
    view(3)

    if(nargout<1)
        clear plotHandle
    end
 end