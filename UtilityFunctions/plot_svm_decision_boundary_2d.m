function plotHandle = plot_svm_decision_boundary_2d(graphicsHandle,svm,varargin)
%PLOT_SVM_DECISION_BOUNDARY_2D Visualize the decision boundary of a SVM
%   PLOT_SVM_DECISION_BOUNDARY_2D plots a visualization of a Support Vector 
%   Machine decision boundary in 2D. It estimates where the boundary is via the
%   scores applied to a large sample, using 'contourc' to find where that score 
%   is closest to 0 and 'patch' to plot it.
%
%   PLOT_SVM_DECISION_BOUNDARY_2D(FIGUREHANDLE,SVM) plots in figure FIGUREHANDLE
%   the decision boundary of Support Vector Machine SVM.
%
%   PLOT_SVM_DECISION_BOUNDARY_2D(...,NAME,VALUE,...) allows for setting
%   additional options for the plot operation; said options should refer to the
%   'patch' function.
%
%   PLOTHANDLE = PLOT_SVM_DECISION_BOUNDARY_2D(...) returns the handle of the
%   generated plot PLOTHANDLE as a result of 'contour', which can be used later 
%   for 'legend', for example.
%
%   Inputs:
%       - FIGUREHANDLE : Figure
%       - SVM : ClassificationSVM
%       - Name-value pair arguments: passed directly to 'contour'.
%
%   Output:
%       - PLOTHANDLE : patch object
%
%   See also patch, contourc, ClassificationSVM, legend.
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
    xInterval = min(svm.X(:,1)) + (0:stepSize:1)*(max(svm.X(:,1))-min(svm.X(:,1)));
    yInterval = min(svm.X(:,2)) + (0:stepSize:1)*(max(svm.X(:,2))-min(svm.X(:,2)));

    % generate grid for predictions at finer sample rate
    [xGrid, yGrid] = meshgrid(xInterval,yInterval);
    fullGrid = [xGrid(:),yGrid(:)];

    % get scores
    [~,score] = predict(svm,fullGrid);

    % reshape to same grid size as the input
    score = reshape(score(:,1), size(xGrid));

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