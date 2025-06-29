function designTypeHandle = plot_design_space_projection_axis(plotData,varargin)
% PLOT_DESIGN_SPACE_PROJECTION_AXIS Visualizes 2D projection of design space samples.
%
%   designTypeHandle = PLOT_DESIGN_SPACE_PROJECTION_AXIS(plotData, varargin) 
%   plots design samples in a 2D axis projection, categorizing points based 
%   on performance and feasibility criteria. It supports extensive plotting 
%   customization via name-value pair arguments.
%
%   INPUTS:
%       plotData : Struct with fields:
%           - DesignSample            : Nx2 array of design variable samples.
%           - IsPhysicallyFeasible    : Nx1 logical vector, true if design physically feasible.
%           - RequirementViolated     : Nx1 integer vector indicating violated requirement index (0 if none).
%           - DesignSpaceLowerBound   : 2x1 vector, lower bounds of design space.
%           - DesignSpaceUpperBound   : 2x1 vector, upper bounds of design space.
%           - DesignBox               : 2x2 matrix, bounds of design box [xmin xmax; ymin ymax].
%
%   OPTIONAL NAME-VALUE PAIR ARGUMENTS:
%       'DataManager'                  : SolutionSpaceData instance or empty (default: new SolutionSpaceData).
%       'CurrentAxis'                  : MATLAB axes handle to plot on (default: []).
%       'AxesLabels'                  : Numeric or string array for axis labels (default: {}).
%       'MarkerColorsViolatedRequirements' : Cell array or color specification for violation markers (default: red).
%       'PlotIntervals'               : Logical flag to plot dashed interval lines (default: true).
%       'PlotOptionsGood'             : Cell array of line/marker properties for good performance points (default: green circles).
%       'PlotOptionsBad'              : Cell array of line/marker properties for bad performance points (default: red crosses).
%       'PlotOptionsPhysicallyInfeasible' : Cell array of line/marker properties for infeasible points (default: grey dots).
%       'PlotOptionsIntervals'        : Cell array of line properties for interval lines (default: black dashed).
%       'PlotOptionsBox'              : Cell array of line properties for design box contour (default: black solid).
%
%   OUTPUT:
%       designTypeHandle : Struct containing plot handles for categories:
%           - PhysicallyInfeasible  : Plot handle for infeasible designs.
%           - GoodPerformance      : Plot handle for good designs.
%           - BadPerformance       : Array of plot handles, one per violated requirement.
%
%   DESCRIPTION:
%       This function plots design points color-coded and styled according 
%       to their performance and feasibility status. It plots dashed lines 
%       to denote intervals of the design box and can label axes. It supports 
%       interactive callbacks for data selection on the plot.
%
%   DEPENDENCIES:
%       - color_palette_tol: Function returning colors by name.
%       - merge_name_value_pair_argument: Helper for merging default and user plot options.
%       - select_data_from_plot: Callback to handle plot point selection.
%       - SolutionSpaceData: Class or struct for design space data management.
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor, Main Author)
%   Copyright 2025 Ali Abbas Kapadia (Contributor)
%   SPDX-License-Identifier: Apache-2.0

    parser = inputParser;
    parser.addParameter('DataManager', SolutionSpaceData(), @(x) isa(x, 'SolutionSpaceData') || isempty(x));
    parser.addParameter('CurrentAxis',[],@(x)isa(x,'matlab.graphics.axis.Axes'));
    parser.addParameter('AxesLabels',{},@(x) isnumeric(x) || isstring(x));
    parser.addParameter('MarkerColorsViolatedRequirements',color_palette_tol('red'));
    parser.addParameter('PlotIntervals',true,@(x)islogical(x));
    parser.addParameter('PlotOptionsGood',{},@(x)iscell(x));
    parser.addParameter('PlotOptionsBad',{},@(x)iscell(x));
    parser.addParameter('PlotOptionsPhysicallyInfeasible',{},@(x)iscell(x));
    parser.addParameter('PlotOptionsIntervals',{},@(x)iscell);
    parser.addParameter('PlotOptionsBox',{},@(x)iscell(x));
    parser.parse(varargin{:});
    options = parser.Results;

    defaultPlotOptionsGood = {'Linestyle','none','Marker','o','Color',color_palette_tol('green')};
    [~,plotOptionsGood] = merge_name_value_pair_argument(defaultPlotOptionsGood,options.PlotOptionsGood);

    defaultPlotOptionsBad = {'Linestyle','none','Marker','x'};
    [~,plotOptionsBad] = merge_name_value_pair_argument(defaultPlotOptionsBad,options.PlotOptionsBad);

    defaultPlotOptionsPhysicallyInfeasible = {'Linestyle','none','Marker','.','Color',color_palette_tol('grey')};
    [~,plotOptionsPhysicallyInfeasible] = merge_name_value_pair_argument(defaultPlotOptionsPhysicallyInfeasible,...
        options.PlotOptionsPhysicallyInfeasible);

    defaultPlotOptionsIntervals = {'LineStyle','--','Color','k','Linewidth',1.5};
    [~,plotOptionsIntervals] = merge_name_value_pair_argument(defaultPlotOptionsIntervals,options.PlotOptionsIntervals);

    defaultPlotOptionsBox = {'Linestyle','-','EdgeColor','k','Linewidth',2.0};
    [~,plotOptionBox] = merge_name_value_pair_argument(defaultPlotOptionsBox,options.PlotOptionsBox);


    %% Extract Data
    designSample = plotData.DesignSample;
    isPhysicallyFeasible = plotData.IsPhysicallyFeasible;
    iRequirementViolated = plotData.RequirementViolated;
    designSpaceLowerBound = plotData.DesignSpaceLowerBound;
    designSpaceUpperBound = plotData.DesignSpaceUpperBound;
    designBox = plotData.DesignBox;
    

    isGoodPerformance = (iRequirementViolated==0);
    nRequirement = max(max(iRequirementViolated),length(options.MarkerColorsViolatedRequirements));
    

    %% Plotting
    currentAxis = parser.Results.CurrentAxis;

    % physically infeasible designs
    if(any(~isPhysicallyFeasible))
        designTypeHandle.PhysicallyInfeasible = plot(currentAxis,...
            designSample(~isPhysicallyFeasible,1),...
            designSample(~isPhysicallyFeasible,2),...
            plotOptionsPhysicallyInfeasible{:},...
            'ButtonDownFcn', @(src,event) select_data_from_plot(src,event,options.DataManager), ...
            'PickableParts', 'all', ...
            'UserData', find(~isPhysicallyFeasible));
    else
        designTypeHandle.PhysicallyInfeasible = matlab.graphics.GraphicsPlaceholder;
    end
    
    hold(currentAxis, 'on');
    

    % good designs
    if(any(isGoodPerformance & isPhysicallyFeasible))
        designTypeHandle.GoodPerformance = plot(currentAxis,...
            designSample(isGoodPerformance & isPhysicallyFeasible,1),...
            designSample(isGoodPerformance & isPhysicallyFeasible,2),...
            plotOptionsGood{:},...
            'ButtonDownFcn', @(src,event) select_data_from_plot(src,event,options.DataManager), ...
            'PickableParts', 'all', ...
            'UserData', find(isGoodPerformance & isPhysicallyFeasible));
    else
        designTypeHandle.GoodPerformance = matlab.graphics.GraphicsPlaceholder;
    end

   hold(currentAxis, 'on');
    
    % bad designs
    for i=1:nRequirement
        if(iscell(options.MarkerColorsViolatedRequirements))
            violateColor = options.MarkerColorsViolatedRequirements{i};
        else
            violateColor = options.MarkerColorsViolatedRequirements;
        end
        
        iViolated = find(iRequirementViolated==i);
        
        if(~isempty(iViolated))
            designTypeHandle.BadPerformance(i) = plot(currentAxis,...
                designSample(iViolated,1),...
                designSample(iViolated,2),...
                'Color',violateColor,...
                plotOptionsBad{:},...
                'ButtonDownFcn', @(src,event) select_data_from_plot(src,event,options.DataManager), ...
                'PickableParts', 'all', ...
                'UserData', find(iRequirementViolated == i));
        else
            designTypeHandle.BadPerformance(i) = matlab.graphics.GraphicsPlaceholder;
        end
    end

    hold(currentAxis, 'on');
    
    % dashed box limits
    if(options.PlotIntervals)
        plot(currentAxis,...
            [designBox(1,1) designBox(1,1)],...
            [designSpaceLowerBound(2) designSpaceUpperBound(2)],...
            plotOptionsIntervals{:}); % vertical left line

        hold(currentAxis, 'on');


        plot(currentAxis,...
            [designBox(2,1) designBox(2,1)],...
            [designSpaceLowerBound(2) designSpaceUpperBound(2)],...
            plotOptionsIntervals{:}); % vertical right line

        hold(currentAxis, 'on');

        plot(currentAxis,...
            [designSpaceLowerBound(1) designSpaceUpperBound(1)],...
            [designBox(1,2) designBox(1,2)],...
            plotOptionsIntervals{:}); % horizontal lower line
        
        hold(currentAxis, 'on');

        plot(currentAxis,...
            [designSpaceLowerBound(1) designSpaceUpperBound(1)],...
            [designBox(2,2) designBox(2,2)],...
            plotOptionsIntervals{:}); % horizontal upper line
    end

    hold(currentAxis, 'on');

    % box contour
    % plot_design_box_2d(currentAxis,designBox,plotOptionBox{:});
    
    % axis lengths and labels
    axis(currentAxis,[...
        designSpaceLowerBound(1) ...
        designSpaceUpperBound(1) ...
        designSpaceLowerBound(2) ...
        designSpaceUpperBound(2)]);
    if(~isempty(options.AxesLabels))
        xlabel(currentAxis, sprintf('x_%d', options.AxesLabels(1)));
        ylabel(currentAxis, sprintf('x_%d', options.AxesLabels(2)));
    end

    drawnow;
end
