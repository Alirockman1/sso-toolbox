function [figureElementHandle,outputData] = plot_selective_design_space_projection(varargin)
%PLOT_SELECTIVE_DESIGN_SPACE_PROJECTION Pair-view of box-shaped solution space
%   PLOT_SELECTIVE_DESIGN_SPACE_PROJECTION projects the solution space design
%   box into 2D-planes and allows for the visualization of them together with
%   designs around it. 
%
%   PLOT_SELECTIVE_DESIGN_SPACE_PROJECTION(DESIGNEVALUATOR,DESIGNBOX,
%   DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND,DESIREDPAIRS,PLOTGRID) uses the
%   design evaluator DESIGNEVALUATOR to evaluate design sample points inside / 
%   outside the solution space box DESIGNBOX in the design space defined by
%   DESIGNSPACELOWERBOUND and DESIGNSPACEUPPERBOUND, and plots those sample
%   points with appropriate colors depending on the result. The plots produced
%   are pair-wise plots as defined by DESIREDPAIRS, and are on a grid of 
%   subplots as specified by PLOTGRID.
%
%   PLOT_SELECTIVE_DESIGN_SPACE_PROJECTION(DESIGNEVALUATOR,DESIGNBOX,
%   DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND,DESIREDPAIRS,PLOTGRID,
%   FIGUREHANDLE) additionally allows one to specify an existing figure to 
%   use for the plots, instead of creating a new one.
%
%   PLOT_SELECTIVE_DESIGN_SPACE_PROJECTION(...NAME,VALUE,...) allows the 
%   specification of additional options. These are:
%       - 'NumberSamplesPerPlot' : how many sample points are used in each 
%       projection plot. Default: 3000.
%       - 'AxesLabels' : labels for each design variable. If empty, nothing is 
%       shown. Default is empty.
%       - 'MarkerColorsViolatedRequirements' : color of the points for each
%       performance requirement which is violated. If only one color is given,
%       it is used for all performance requirements. Default: 'r'.
%       - 'MarkerColorsCriterion' : criterion to select a color when multiple
%       performance requirements are violated. Two options are available: 
%       'random', where colors are picked randomly from the violated
%       requirements; and 'worst', where the color of the biggest (absolute) 
%       deficit is used. Default is ' worst'.
%       - 'SamplingMethod' : sampling method used to generate the samples for
%       the plots. Default: @sampling_latin_hypercube.
%       - 'SamplingOptions' : extra options for the sampling method. Default
%       is empty.
%       - 'PlotIntervals' : true/false flag indicating whether the dashed 
%       intervals for the design variables should be printed or not. Default:
%       true.
%       - 'PlotOptionsGood' : options for the plotting of design points with 
%       good performance. Default: {'Linestyle','none','Marker','.','Color','g'}
%       - 'PlotOptionsBad' : options for the plotting of design points with 
%       bad performance. Default: {'Linestyle','none','Marker','.','Color','r'}.
%       - 'PlotOptionsPhysicallyInfeasible' : options for the plotting of design 
%       points that are physically infeasible. Default: {'Linestyle','none',
%       'Marker','.','Color',[0.8 0.8 0.8]}.
%       - 'PlotOptionsIntervals' : options for the plotting of the intervals of
%       the solution space design box. Default: {'LineStyle','--','Color','k',
%       'Linewidth',1.5}.
%       - 'PlotOptionsBox' : options for the plotting of the solution box. Uses
%       the 'plot_design_box_2d' function. Default: {'Linestyle','-',
%       'EdgeColor','k','Linewidth',2.0}.
%
%   FIGUREELEMENTHANDLE = PLOT_SELECTIVE_DESIGN_SPACE_PROJECTION(...) returns 
%   the handles for the main elements of the figure as a whole. These are:
%       - MainFigure : handle for the main figure.
%       - IndividualPlots : handle for each subplot.
%       - PhysicallyInfeasible : handle for physically infeasible points.
%       - GoodPerformance : handle for good performance points.
%       - BadPerformance : handles for each type of bad performance points.
%
%   [FIGUREELEMENTHANDLE,PROBLEMDATA] = 
%   PLOT_SELECTIVE_DESIGN_SPACE_PROJECTION(...) additionally returns the fixed 
%   problem data in PROBLEMDATA.
%
%   [FIGUREELEMENTHANDLE,PROBLEMDATA,PLOTDATA] = 
%   PLOT_SELECTIVE_DESIGN_SPACE_PROJECTION(...) additionally returns all the 
%   data generated for each plot.
%       - DesignSample : design sample points evaluated.
%       - PerformanceDeficit : performance deficit for each design sample point
%       - PhysicalFeasibilityDeficit : physical feasibility deficit for each 
%       design sample point.
%       - EvaluatorOutput : extra output from the evaluator; depends on the 
%       class.
%
%   Inputs:
%       - DESIGNEVALUATOR : DesignEvaluatorBase
%       - DESIGNBOX : (2,nDesignVariable) double
%       - DESIGNSPACELOWERBOUND : (1,nDesignVariable) double
%       - DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double
%       - DESIREDPAIRS : (nPlot,2) integer
%       - PLOTGRID : (1,2) integer
%       - 'NumberSamplesPerPlot' : integer
%       - 'AxesLabels' : (1,nDesignVariable) cell
%       - 'MarkerColorsViolatedRequirements' : (1,nPerformance) cell
%       - 'MarkerColorsCriterion' : char OR string
%       - 'SamplingMethod' : function_handle
%       - 'SamplingOptions' : (1,nOptionSampling) cell
%       - 'PlotOptionsGood' : (1,nOptionPlotGood) cell
%       - 'PlotOptionsBad' : (1,nOptionPlotBad) cell
%       - 'PlotOptionsPhysicallyInfeasible' : (1,nOptionPlotPhysicalFeasibility)
%       cell
%       - 'PlotOptionsIntervals' : (1,nOptionPlotInterval) cell
%       - 'PlotOptionsBox' : (1,nOptionPlotBox) cell
%
%   Outputs:
%       - FIGUREELEMENTHANDLE : Figure
%           -- MainFigure : Figure
%           -- IndividualPlots : (1,nPlot) Figure
%           -- PhysicallyInfeasible : Line
%           -- GoodPerformance : Line
%           -- BadPerformance : (1,nPerformance) Line
%       - PLOTDATA : (1,nPlot) structure 
%           -- DesignSample : (nSample,nDesignVariable) double
%           -- PerformanceDeficit : (nSample,nPerformance) double
%           -- PhysicalFeasibilityDeficit : (nSample,nPhysicalFeasibility) 
%           double
%           -- EvaluatorOutput : class-defined
%
%   See also sso_box_stochastic.
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

    parser = inputParser;
    parser.KeepUnmatched = true;
    parser.addRequired('DesignEvaluator',@(x)isa(x,'DesignEvaluatorBase'));
    parser.addRequired('DesignBox',@(x)isnumeric(x)&&size(x,1)==2);
    parser.addRequired('DesignSpaceLowerBound',@(x)isnumeric(x)&&size(x,1)==1);
    parser.addRequired('DesignSpaceUpperBound',@(x)isnumeric(x)&&size(x,1)==1);
    parser.addRequired('DesiredPairs',@(x)(isnumeric(x)&&size(x,2)==2));
    parser.addRequired('PlotGrid',@(x)(isnumeric(x)&&length(x)==2));
<<<<<<< HEAD
    parser.addParameter('NumberSamplesPerPlot',3000,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    parser.addParameter('AxesLabels',{},@(x)iscell(x));
    parser.addParameter('MarkerColorsViolatedRequirements','r');
    parser.addParameter('MarkerColorsCriterion','worst',@(x)any(stcpmi(x,{'worst','random'})));
    parser.addParameter('SamplingMethod',@sampling_latin_hypercube,@(x)isa(x,'function_handle'));
    parser.addParameter('SamplingOptions',{},@(x)iscell(x));
    parser.addParameter('PlotOptionsGood',{},@(x)iscell(x));
    parser.addParameter('PlotOptionsBad',{},@(x)iscell(x));
    parser.addParameter('PlotOptionsPhysicallyInfeasible',{},@(x)iscell(x));
    parser.addParameter('PlotOptionsIntervals',{},@(x)iscell);
    parser.addParameter('PlotOptionsBox',{},@(x)iscell(x));
    parser.parse(designEvaluator, designBox, designSpaceLowerBound, designSpaceUpperBound, desiredPairs, plotGrid, varargin{:});
    options = parser.Results;

    defaultPlotOptionsGood = {'Linestyle','none','Marker','.','Color','g'};
    [~,plotOptionsGood] = merge_name_value_pair_argument(defaultPlotOptionsGood,options.PlotOptionsGood);

    defaultPlotOptionsBad = {'Linestyle','none','Marker','.','Color','r'};
    [~,plotOptionsBad] = merge_name_value_pair_argument(defaultPlotOptionsBad,options.PlotOptionsBad);

    defaultPlotOptionsPhysicallyInfeasible = {'Linestyle','none','Marker','.','Color',[0.8 0.8 0.8]};
    [~,plotOptionsPhysicallyInfeasible] = merge_name_value_pair_argument(defaultPlotOptionsPhysicallyInfeasible,...
        options.PlotOptionsPhysicallyInfeasible);

    defaultPlotOptionsIntervals = {'LineStyle','--','Color','k','Linewidth',1.5};
    [~,plotOptionsIntervals] = merge_name_value_pair_argument(defaultPlotOptionsIntervals,options.PlotOptionsIntervals);

    defaultPlotOptionsBox = {'Linestyle','-','Color','k','Linewidth',2.0};
    [~,plotOptionBox] = merge_name_value_pair_argument(defaultPlotOptionsBox,options.PlotOptionsBox);
=======
    parser.addOptional('CurrentFigure',[],@(x)isa(x,'matlab.ui.Figure')||isempty(x));
    parser.parse(varargin{:});
    
    designEvaluator = parser.Results.DesignEvaluator;
    designBox = parser.Results.DesignBox;
    designSpaceLowerBound = parser.Results.DesignSpaceLowerBound;
    designSpaceUpperBound = parser.Results.DesignSpaceUpperBound;
    desiredPairs = parser.Results.DesiredPairs;
    plotGrid = parser.Results.PlotGrid;
    currentFigure = parser.Results.CurrentFigure;
    options = parser.Unmatched;
>>>>>>> f6a30e3627663a9ce0c7714b18eba9603ade2e9e
    
    
    %% Pre-allocate important arrays
    isOutputData = (nargout>=2);
    if(isOutputData)
        outputData = struct(...
            'DesignEvaluator',designEvaluator,...
            'DesignBox',designBox,...
            'DesignSpaceLowerBound',designSpaceLowerBound,...
            'DesignSpaceUpperBound',designSpaceUpperBound,...
            'DesiredPairwisePlots',desiredPairs,...
            'PlotGrid',plotGrid,...
            'Options',options,...
            'InitialRNGState',rng);
    end

    evaluationFields = {'SamplingMethod','SamplingOptions','NumberSamplesPerPlot','MarkerColorsCriterion'};
    evaluationOptions = namedargs2cell(structure_extract_fields(options,evaluationFields));
    evaluationInput = [{designEvaluator,designBox,designSpaceLowerBound,designSpaceUpperBound,desiredPairs},evaluationOptions];
    if(isOutputData)
        [plotData,outputData.EvaluationData] = evaluate_selective_design_space_projection(evaluationInput{:});
    else
        plotData = evaluate_selective_design_space_projection(evaluationInput{:});
    end

    
    %% Start Figure
    if(isempty(currentFigure))
        figureElementHandle.MainFigure = figure;
    else
        figureElementHandle.MainFigure = currentFigure;
    end
    figure(figureElementHandle.MainFigure);
    hold all; 
    drawnow;

    plotFields = {'MarkerColorsViolatedRequirements','PlotIntervals','PlotOptionsGood','PlotOptionsBad','PlotOptionsPhysicallyInfeasible','PlotOptionsIntervals','PlotOptionsBox'};
    plotOptionsBase = namedargs2cell(structure_extract_fields(options,plotFields));
    for j=1:size(desiredPairs,1)
        %% Random Sampling for Relevant Scatter Plot
        currentPair = [desiredPairs(j,1),desiredPairs(j,2)];
        if(isfield(options,'AxesLabels') && ~isempty(options.AxesLabels))
            currentPairLabels = options.AxesLabels(currentPair);
            plotOptions = [plotOptionsBase,{'AxesLabels',currentPairLabels}];
        else
            plotOptions = plotOptionsBase;
        end

        %% Plot Solution
        figureElementHandle.IndividualPlots(j) = subplot(plotGrid(1),plotGrid(2),j);
        plotTag = sprintf('Plot%d', j);
        set(gca, 'Tag', plotTag); % Tag for identification
        hold all;
        designTypeHandle = plot_design_space_projection_axis(plotData(j),plotOptions{:});

        % merge handles
        if(j==1)
            figureElementHandle.GoodPerformance = designTypeHandle.GoodPerformance;
            figureElementHandle.PhysicallyInfeasible = designTypeHandle.PhysicallyInfeasible;
            figureElementHandle.BadPerformance = designTypeHandle.BadPerformance;
        else
            if(isa(figureElementHandle.GoodPerformance,'matlab.graphics.GraphicsPlaceholder') && ...
               ~isa(designTypeHandle.GoodPerformance,'matlab.graphics.GraphicsPlaceholder'))
                figureElementHandle.GoodPerformance = designTypeHandle.GoodPerformance;
            end
            
            if(isa(figureElementHandle.PhysicallyInfeasible,'matlab.graphics.GraphicsPlaceholder') && ...
               ~isa(designTypeHandle.PhysicallyInfeasible,'matlab.graphics.GraphicsPlaceholder'))
                figureElementHandle.PhysicallyInfeasible = designTypeHandle.PhysicallyInfeasible;
            end
            
            for k=1:length(designTypeHandle.BadPerformance)
                if(k>length(figureElementHandle.BadPerformance) || ...
                    (isa(figureElementHandle.BadPerformance(k),'matlab.graphics.GraphicsPlaceholder') && ...
                    ~isa(figureElementHandle.BadPerformance(k),'matlab.graphics.GraphicsPlaceholder')))
                    figureElementHandle.BadPerformance(k) = designTypeHandle.BadPerformance(k);
                end
            end
        end
    end


    if(nargout<1)
        clear figureElementHandle
    end
end


