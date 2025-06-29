function plotHandle = plot_svm_training_data_3d(figureHandle,svm,options)
%PLOT_SVM_TRAINING_DATA_3D Visualize the training data of a SVM
%   PLOT_SVM_TRAINING_DATA_3D plots all the samples used in the training of a
%	Support Vector Machine in a color-coded way; it also highlights the points
%	used as support vectors. 
%
%   PLOT_SVM_TRAINING_DATA_3D(FIGUREHANDLE,SVM) plots in figure FIGUREHANDLE
%   the training data of Support Vector Machine SVM. Designs labeled as 
%	'true' are shown as green, 'false' are shown as red, and the support vector
%	have an additional blue circle around them.
%
%   PLOT_SVM_TRAINING_DATA_3D(FIGUREHANDLE,SVM,OPTIONS) allows for setting
%   additional options for the plot operation; said options should refer to the
%   'plot' function. OPTIONS should have four fields:
%		- 'OverallOptions' : 'plot3' options which apply to all plots
%		- 'TrueLabelOptions' : 'plot3' options which apply to the plot of 'true'
%		label points
%		- 'FalseLabelOptions' : 'plot3' options which apply to the plot of 
%		'false' label points
%		- 'SupportVectorOptions' : 'plot3' options which apply to the plot of 
%		the support vectors
%	In terms of priority, 'OverallOptions' has the lowest priority, the default
%	options for each plot has medium priority, and user-definied specific 
%	options have the highest priority. 
%	By default:
%		- 'OverallOptions' has the 'Marker' set to '.'
%		- 'TrueLabelOptions' has the 'Color' set to green 'g'
%		- 'FalseLabelOptions' has the 'Color' set to red 'r'
%		- 'SupportVectorOptions' has the 'Marker' set to 'o' and 'Color' to 'b'
%
%   PLOTHANDLE = PLOT_SVM_TRAINING_DATA_3D(...) returns the handles of the
%   generated plots PLOTHANDLE as a result of 'plot', which can be used later 
%   for 'legend', for example. These are the plots of: (1) the 'true' points,
%	(2) the 'false' points, and (3) the support vectors.
%
%   Inputs:
%       - FIGUREHANDLE : Figure
%       - SVM : ClassificationSVM
%		- OPTIONS : (1,4) structure OR name-value pair cell, optional
%			-- 'OverallOptions' : 'plot' options, name-value pair
%			-- 'TrueLabelOptions' : 'plot' options, name-value pair
%			-- 'FalseLabelOptions' : 'plot' options, name-value pair
%			-- 'SupportVectorOptions' : 'plot' options, name-value pair
%
%   Output:
%       - PLOTHANDLE : (1,3) Line
%
%   See also plot3, ClassificationSVM, legend, plot_svm_training_data_2d.
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

    % Assume class labels are 1 and 0 and convert to logical
    label = logical(svm.Y);

    % process options
    defaultOptions = {'TrueLabelOptions',[],'FalseLabelOptions',[],...
    	'SupportVectorOptions',[],'OverallOptions',[]};
    options = merge_name_value_pair_argument(defaultOptions,options);

    defaultOverallOptions = {'Marker','.'};
    [~,overallOptions] = merge_name_value_pair_argument(defaultOverallOptions,options.OverallOptions);

    defaultTrueOptions = {'Color','g'};
    [~,trueOptions] = merge_name_value_pair_argument(overallOptions,defaultTrueOptions,...
    	options.TrueLabelOptions);

    defaultFalseOptions = {'Color','r'};
    [~,falseOptions] = merge_name_value_pair_argument(overallOptions,defaultFalseOptions,...
    	options.FalseLabelOptions);

    defaultSupportVectorOptions = {'Marker','o','Color','b'};
    [~,supportVectorOptions] = merge_name_value_pair_argument(overallOptions,defaultSupportVectorOptions,...
    	options.SupportVectorOptions);

    % Plot data points, color by class label
    plotHandle(1) = plot3(svm.X(label,1), svm.X(label,2), svm.X(label,3), 'LineStyle','none', trueOptions);
    plotHandle(2) = plot3(svm.X(~label,1), svm.X(~label,2), svm.X(~label,3), 'LineStyle','none', falseOptions);

    % Re-scale the support vectors and plot them
    sv = svm.SupportVectors.*svm.Sigma + svm.Mu;
    plotHandle(3) = plot3(sv(:, 1), sv(:, 2), sv(:, 3), 'LineStyle','none', supportVectorOptions);

    if(nargout<1)
        clear plotHandle
    end
end