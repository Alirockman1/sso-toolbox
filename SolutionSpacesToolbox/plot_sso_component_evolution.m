function figureHandles = plot_sso_component_evolution(optimizationData, varargin)
%PLOT_SSO_COMPONENT_EVOLUTION Plot evolution of components across iterations
%   PLOT_SSO_COMPONENT_EVOLUTION visualizes how components evolve during the
%   stochastic solution space optimization process. For each iteration and
%   component, it plots the initial component (same as trimmed component from
%   previous iteration), grown component (where applicable), and trimmed
%   component, all in the same figure with fixed axes to the design space.
%
%   FIGUREHANDLES = PLOT_SSO_COMPONENT_EVOLUTION(OPTIMIZATIONDATA) creates
%   figures showing the evolution of each component across iterations using
%   the data in OPTIMIZATIONDATA from sso_component_stochastic.
%
%   FIGUREHANDLES = PLOT_SSO_COMPONENT_EVOLUTION(OPTIMIZATIONDATA, NAME, VALUE)
%   allows specification of additional options:
%       - 'ComponentIndices': Indices of components to plot. Default: all.
%       - 'IterationIndices': Indices of iterations to plot. Default: all.
%       - 'PlotGrid': Whether to show grid. Default: true.
%       - 'PlotMinorGrid': Whether to show minor grid. Default: false.
%       - 'PlotTitle': Whether to show title. Default: true.
%       - 'PlotLegend': Whether to show legend. Default: true.
%       - 'VariableNames': Cell array of variable names for axis labels. Default: [] (uses LaTeX x₁, x₂, etc.).
%       - 'PlotEvaluatedPoints': Whether to plot evaluated design points. Default: false.
%       - 'AcceptablePointsColor': Color for acceptable points. Default: 'g'.
%       - 'UnacceptablePointsColor': Color for unacceptable points. Default: 'r'.
%       - 'PointSize': Size of plotted points. Default: 20.
%       - 'PointMarker': Marker style for points. Default: '.'.
%       - 'PointAlpha': Transparency of points (0-1). Default: 0.5.
%       - 'InitialColor': Color for initial component. Default: 'b'.
%       - 'GrownColor': Color for grown component. Default: 'r'.
%       - 'TrimmedColor': Color for trimmed component. Default: 'g'.
%       - 'InitialLineStyle': Line style for initial component. Default: '-'.
%       - 'GrownLineStyle': Line style for grown component. Default: '--'.
%       - 'TrimmedLineStyle': Line style for trimmed component. Default: '-'.
%       - 'InitialLineWidth': Line width for initial component. Default: 1.
%       - 'GrownLineWidth': Line width for grown component. Default: 1.
%       - 'TrimmedLineWidth': Line width for trimmed component. Default: 2.
%       - 'SaveFolder': Folder to save figures. Default: [].
%       - 'CloseFigureAfterSaving': Whether to close figures after saving. Default: false.
%       - 'SaveFigureOptions': Options for saving figures. Default: {}.
%       - 'InitialPlotOptions': Additional options for plotting initial component. Default: {}.
%       - 'GrownPlotOptions': Additional options for plotting grown component. Default: {}.
%       - 'TrimmedPlotOptions': Additional options for plotting trimmed component. Default: {}.
%
%   Inputs:
%       - OPTIMIZATIONDATA: struct from sso_component_stochastic
%       - Various NAME, VALUE pairs as described above
%
%   Outputs:
%       - FIGUREHANDLES: Array of figure handles
%
%   See also sso_component_stochastic, CandidateSpaceBase.plot_candidate_space.
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

    % Parse input arguments
    parser = inputParser;
    parser.addRequired('optimizationData', @isstruct);
    parser.addParameter('IterationIndices', [], @(x) isempty(x) || (isnumeric(x) && all(x > 0)));
    parser.addParameter('ComponentIndices', [], @(x) isempty(x) || (isnumeric(x) && all(x > 0)));
    parser.addParameter('GridStyle', 'minor');
    parser.addParameter('PlotTitle', true, @islogical);
    parser.addParameter('PlotLegend', true, @islogical);
    parser.addParameter('VariableNames', [], @(x) isempty(x) || iscell(x));
    parser.addParameter('PlotEvaluatedPoints', true, @islogical);
    parser.addParameter('AcceptablePointsColor', 'g', @(x) ischar(x) || isstring(x) || (isnumeric(x) && length(x) == 3));
    parser.addParameter('UnacceptablePointsColor', 'r', @(x) ischar(x) || isstring(x) || (isnumeric(x) && length(x) == 3));
    parser.addParameter('ScatterPointSize', 40, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parser.addParameter('ScatterOptions', {}, @iscell);
    parser.addParameter('InitialPlotOptions', {}, @iscell);
    parser.addParameter('GrownPlotOptions', {}, @iscell);
    parser.addParameter('TrimmedPlotOptions', {}, @iscell);
    parser.addParameter('LabelOptions', {'interpreter','latex'}, @iscell);
    parser.addParameter('SaveFolder', [], @(x) isempty(x) || ischar(x) || isstring(x));
    parser.addParameter('CloseFigureAfterSaving', false, @islogical);
    parser.addParameter('SaveFigureOptions', {}, @iscell);
    parser.parse(optimizationData, varargin{:});
    options = parser.Results;

    if(isempty(options.IterationIndices) && isempty(options.ComponentIndices))
        warning('No iteration or component indices provided, returning.');
        figureHandles = [];
        return;
    end
    
    % Get iteration indices
    nIterations = length(optimizationData.IterationData);
    nComponents = length(optimizationData.ComponentIndex);

    if(isempty(options.IterationIndices))
        options.IterationIndices = 1:nIterations;
    end 
    if(isempty(options.ComponentIndices))
        options.ComponentIndices = 1:nComponents;
    end

    % Initialize figure handles array
    figureHandles = gobjects(nComponents, nIterations);
    
    % Loop through components and iterations
    for iComponent = options.ComponentIndices        
        % Create component folder if saving is enabled
        componentFolder = [];
        if ~isempty(options.SaveFolder)
            componentFolder = [options.SaveFolder, sprintf('Component%d/', iComponent)];
            if ~exist(componentFolder, 'dir')
                mkdir(componentFolder);
            end
        end

        % Create variable names for axis labels
        currentComponentIndex = optimizationData.ComponentIndex{iComponent};
        nDimension = length(currentComponentIndex);
        componentVariableNames = cell(1, nDimension);
        if(~isempty(options.VariableNames) && length(options.VariableNames) >= max(currentComponentIndex))
            % Use provided variable names
            for i = 1:nDimension
                componentVariableNames{i} = options.VariableNames{currentComponentIndex(i)};
            end
        else
            % Use default LaTeX formatted variable names (x₁, x₂, etc.)
            for i = 1:nDimension
                componentVariableNames{i} = ['$$x_{',num2str(currentComponentIndex(i)),'}$$'];
            end
        end
        
        for iIteration = options.IterationIndices
            % Create figure
            figureHandles(iComponent, iIteration) = figure;
            
            % Set title if requested
            if options.PlotTitle
                if optimizationData.IterationData(iIteration).Phase == 1
                    phaseType = 'Exploration';
                else
                    phaseType = 'Consolidation';
                end
                title(sprintf('Component %d - %s Phase - Total Iteration %d', iComponent, phaseType, iIteration));
            end

            % Plot components
            hold on;
            
            % Plot initial component (previous iteration's trimmed component)
            initialPlot = [];
            if(iIteration > 1)
                initialComponent = optimizationData.IterationData(iIteration-1).CandidateSpacesAfterTrim(iComponent);
                initialPlot = initialComponent.plot_candidate_space(gcf, options.InitialPlotOptions{:});
            end
            
            % Plot grown component (only in exploration phase)
            grownPlot = [];
            if(~isempty(optimizationData.IterationData(iIteration).CandidateSpacesBeforeTrim) && optimizationData.IterationData(iIteration).Phase == 1)
                grownComponent = optimizationData.IterationData(iIteration).CandidateSpacesBeforeTrim(iComponent);
                grownPlot = grownComponent.plot_candidate_space(gcf, options.GrownPlotOptions{:});
            end
            
            % Plot trimmed component
            trimmedPlot = [];
            if(~isempty(optimizationData.IterationData(iIteration).CandidateSpacesAfterTrim))
                trimmedComponent = optimizationData.IterationData(iIteration).CandidateSpacesAfterTrim(iComponent);
                trimmedPlot = trimmedComponent.plot_candidate_space(gcf, options.TrimmedPlotOptions{:});
            end
            
            % Plot evaluated points if requested
            acceptablePointsPlot = [];
            unacceptablePointsPlot = [];
            if options.PlotEvaluatedPoints && ~isempty(optimizationData.IterationData(iIteration).EvaluatedDesignSamples)
                % Get the evaluated design samples and their acceptability status
                designSamples = optimizationData.IterationData(iIteration).EvaluatedDesignSamples;
                isAcceptable = optimizationData.IterationData(iIteration).IsAcceptable;
                
                % Plot points based on dimensionality
                if nDimension == 1
                    % 1D case
                    acceptablePoints = designSamples(isAcceptable, currentComponentIndex);
                    unacceptablePoints = designSamples(~isAcceptable, currentComponentIndex);
                    
                    if ~isempty(acceptablePoints)
                        acceptablePointsPlot = scatter(acceptablePoints, zeros(size(acceptablePoints)), ...
                            options.ScatterPointSize, options.AcceptablePointsColor, options.ScatterOptions{:});
                    end
                    
                    if ~isempty(unacceptablePoints)
                        unacceptablePointsPlot = scatter(unacceptablePoints, zeros(size(unacceptablePoints)), ...
                            options.ScatterPointSize, options.UnacceptablePointsColor, options.ScatterOptions{:});
                    end
                elseif nDimension == 2
                    % 2D case
                    acceptablePoints = designSamples(isAcceptable, currentComponentIndex);
                    unacceptablePoints = designSamples(~isAcceptable, currentComponentIndex);
                    
                    if ~isempty(acceptablePoints)
                        acceptablePointsPlot = scatter(acceptablePoints(:,1), acceptablePoints(:,2), ...
                            options.ScatterPointSize, options.AcceptablePointsColor, options.ScatterOptions{:});
                    end
                    
                    if ~isempty(unacceptablePoints)
                        unacceptablePointsPlot = scatter(unacceptablePoints(:,1), unacceptablePoints(:,2), ...
                            options.ScatterPointSize, options.UnacceptablePointsColor, options.ScatterOptions{:});
                    end
                elseif nDimension == 3
                    % 3D case
                    acceptablePoints = designSamples(isAcceptable, currentComponentIndex);
                    unacceptablePoints = designSamples(~isAcceptable, currentComponentIndex);
                    
                    if ~isempty(acceptablePoints)
                        acceptablePointsPlot = scatter3(acceptablePoints(:,1), acceptablePoints(:,2), acceptablePoints(:,3), options.ScatterPointSize, options.AcceptablePointsColor, options.ScatterOptions{:});
                    end
                    
                    if ~isempty(unacceptablePoints)
                        unacceptablePointsPlot = scatter3(unacceptablePoints(:,1), unacceptablePoints(:,2), unacceptablePoints(:,3), options.ScatterPointSize, options.UnacceptablePointsColor, options.ScatterOptions{:});
                    end
                end
            end
            
            % Set axis limits based on component's design space boundaries
            % Use the most recent component to get the design space boundaries
            componentReference = [];
            if(~isempty(initialPlot))
                componentReference = initialComponent;
            elseif(~isempty(grownPlot))
                componentReference = grownComponent;
            elseif(~isempty(trimmedPlot))
                componentReference = trimmedComponent;
            end
            
            if(~isempty(componentReference))
                if(nDimension == 1)
                    xlim([componentReference.DesignSpaceLowerBound, componentReference.DesignSpaceUpperBound]);
                    xlabel(componentVariableNames{1}, options.LabelOptions{:});
                elseif(nDimension == 2)
                    xlim([componentReference.DesignSpaceLowerBound(1), componentReference.DesignSpaceUpperBound(1)]);
                    ylim([componentReference.DesignSpaceLowerBound(2), componentReference.DesignSpaceUpperBound(2)]);
                    xlabel(componentVariableNames{1}, options.LabelOptions{:});
                    ylabel(componentVariableNames{2}, options.LabelOptions{:});
                elseif(nDimension == 3)
                    xlim([componentReference.DesignSpaceLowerBound(1), componentReference.DesignSpaceUpperBound(1)]);
                    ylim([componentReference.DesignSpaceLowerBound(2), componentReference.DesignSpaceUpperBound(2)]);
                    zlim([componentReference.DesignSpaceLowerBound(3), componentReference.DesignSpaceUpperBound(3)]);
                    xlabel(componentVariableNames{1}, options.LabelOptions{:});
                    ylabel(componentVariableNames{2}, options.LabelOptions{:});
                    zlabel(componentVariableNames{3}, options.LabelOptions{:});
                    view(3); % Set 3D view
                end
            end
            
            % Add grid if requested
            grid(options.GridStyle);
            
            % Add legend if requested
            if options.PlotLegend
                legendEntries = {};
                legendHandles = {};
                
                if(~isempty(initialPlot))
                    legendHandles = [legendHandles, initialPlot];
                    legendEntries{end+1} = 'Initial Component';
                end
                
                if(~isempty(grownPlot))
                    legendHandles = [legendHandles, grownPlot];
                    legendEntries{end+1} = 'Grown Component';
                end
                
                if(~isempty(trimmedPlot))
                    legendHandles = [legendHandles, trimmedPlot];
                    legendEntries{end+1} = 'Trimmed Component';
                end
                
                if(~isempty(acceptablePointsPlot))
                    legendHandles = [legendHandles, acceptablePointsPlot];
                    legendEntries{end+1} = 'Acceptable Points';
                end
                
                if(~isempty(unacceptablePointsPlot))
                    legendHandles = [legendHandles, unacceptablePointsPlot];
                    legendEntries{end+1} = 'Unacceptable Points';
                end
                
                if(~isempty(legendEntries))
                    legend(legendHandles, legendEntries, 'Location', 'best');
                end
            end
            
            % Save figure if requested
            if(~isempty(options.SaveFolder))
                figureName = [componentFolder, sprintf('Iteration%d.png', iIteration)];
                save_print_figure(figureHandles(iComponent, iIteration), figureName, options.SaveFigureOptions{:});
                
                if(options.CloseFigureAfterSaving)
                    close(figureHandles(iComponent, iIteration));
                end
            end
        end
    end
    
    % Clear output if not requested
    if(nargout < 1)
        clear figureHandles;
    end
end 