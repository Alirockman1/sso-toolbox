function designTypeHandle = plot_design_space_projection_axis(plotData,varargin)
    parser = inputParser;
    parser.addOptional('CurrentAxis',gca,@(x)isa(x,'matlab.graphics.axis.Axes'));
    parser.addParameter('AxesLabels',{},@(x)iscell(x));
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
    if isa(currentAxis, 'matlab.ui.Figure') || isa(currentAxis, 'matlab.ui.container.Figure')
		figure(currentAxis);
	elseif isa(currentAxis, 'matlab.graphics.axis.Axes') || isa(currentAxis, 'matlab.ui.control.UIAxes')
		axes(currentAxis);
	end
    hold all;

    % physically infeasible designs
    if(any(~isPhysicallyFeasible))
        designTypeHandle.PhysicallyInfeasible = plot(currentAxis,...
            designSample(~isPhysicallyFeasible,1),...
            designSample(~isPhysicallyFeasible,2),...
            plotOptionsPhysicallyInfeasible{:});
    else
        designTypeHandle.PhysicallyInfeasible = matlab.graphics.GraphicsPlaceholder;
    end
    

    % good designs
    if(any(isGoodPerformance & isPhysicallyFeasible))
        designTypeHandle.GoodPerformance = plot(currentAxis,...
            designSample(isGoodPerformance & isPhysicallyFeasible,1),...
            designSample(isGoodPerformance & isPhysicallyFeasible,2),...
            plotOptionsGood{:});
    else
        designTypeHandle.GoodPerformance = matlab.graphics.GraphicsPlaceholder;
    end
    
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
                plotOptionsBad{:});
        else
            designTypeHandle.BadPerformance(i) = matlab.graphics.GraphicsPlaceholder;
        end
    end
    
    % dashed box limits
    if(options.PlotIntervals)
        plot(currentAxis,...
            [designBox(1,1) designBox(1,1)],...
            [designSpaceLowerBound(2) designSpaceUpperBound(2)],...
            plotOptionsIntervals{:}); % vertical left line


        plot(currentAxis,...
            [designBox(2,1) designBox(2,1)],...
            [designSpaceLowerBound(2) designSpaceUpperBound(2)],...
            plotOptionsIntervals{:}); % vertical right line



        plot(currentAxis,...
            [designSpaceLowerBound(1) designSpaceUpperBound(1)],...
            [designBox(1,2) designBox(1,2)],...
            plotOptionsIntervals{:}); % horizontal lower line
        

        plot(currentAxis,...
            [designSpaceLowerBound(1) designSpaceUpperBound(1)],...
            [designBox(2,2) designBox(2,2)],...
            plotOptionsIntervals{:}); % horizontal upper line
    end

    % box contour
    % plot_design_box_2d(currentAxis,designBox,plotOptionBox{:});
    
    % axis lengths and labels
    axis(currentAxis,[...
        designSpaceLowerBound(1) ...
        designSpaceUpperBound(1) ...
        designSpaceLowerBound(2) ...
        designSpaceUpperBound(2)]);
    if(~isempty(options.AxesLabels))
        xlabel(currentAxis,options.AxesLabels{1});
        ylabel(currentAxis,options.AxesLabels{2});
    end
    grid('off');
    drawnow;
end
