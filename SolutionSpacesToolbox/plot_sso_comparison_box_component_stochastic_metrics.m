function plot_sso_comparison_box_component_stochastic_metrics(algoDataBox,algoDataComponent,varargin)
	parser = inputParser;
    parser.addParameter('SaveFolder',[]);
    parser.addParameter('CloseFigureAfterSaving',false);
    parser.addParameter('SaveFigureOptions',{});
    parser.addParameter('BoxLabel',{});
    parser.addParameter('ComponentLabel',{});
    parser.addParameter('BoxColor',[]);
    parser.addParameter('ComponentColor',[]);
    parser.addParameter('GeneralPlotOptions',{});
    parser.addParameter('BoxPlotOptions',{});
    parser.addParameter('ComponentPlotOptions',{});
    parser.parse(varargin{:});
    options = parser.Results;

    % create different labels for boxes
    nBox = length(algoDataBox);
    if(isempty(options.BoxLabel))
        % create numerical labels if none were given and there's more than one entry
        if(nBox>1)
            options.BoxLabel = cellfun(@num2str,num2cell(1:nBox));
        else
            options.BoxLabel = {{}};
        end
    end
    for i=1:nBox
        % add space at the end if not already included
        if(~isempty(options.BoxLabel{i}) && ~strcmpi(options.BoxLabel{i}(end),' '))
            options.BoxLabel{i} = [options.BoxLabel{i},' '];
        end
    end
    if(isempty(options.BoxColor))
        options.BoxColor = color_palette_tol(1:nBox)
    end

    % create different labels for components
    nComponent = length(algoDataComponent);
    if(isempty(options.ComponentLabel))
        % create numerical labels if none were given and there's more than one entry
        if(nComponent>1)
            options.ComponentLabel = cellfun(@num2str,num2cell(1:nBox));
        else
            options.ComponentLabel = {{}};
        end
    end
    for i=1:nComponent
        % add space at the end if not already included
        if(~isempty(options.ComponentLabel{i}) && ~strcmpi(options.ComponentLabel{i}(end),' '))
            options.ComponentLabel{i} = [options.ComponentLabel{i},' '];
        end
    end
    if(isempty(options.ComponentColor))
        options.BoxColor = color_palette_tol((nBox+1):(nBox+1+nComponent))
    end

    defaultGeneralPlotOptions = {'linewidth',1.1,'Marker','.','MarkerSize',8};
    [~,generalPlotOptions] = merge_name_value_pair_argument(defaultGeneralPlotOptions,options.GeneralPlotOptions);

    defaultBoxPlotOptions = {};
    [~,boxPlotOptions] = merge_name_value_pair_argument(defaultBoxPlotOptions,...
        generalPlotOptions,options.BoxPlotOptions);

    defaultComponentPlotOptions = {};
    [~,componentPlotOptions] = merge_name_value_pair_argument(defaultComponentPlotOptions,...
        generalPlotOptions,options.ComponentPlotOptions);

    %% Plots
    % Volume of Box
    figureHandle(1) = figure;
    hold all;
    legendEntry = {};
    xAxisMinLimit = inf;
    xAxisMaxLimit = -inf;
    explorationEnd = [];
    for i=1:nBox
        plot([algoDataBox{i}.MeasureBeforeTrim],'Color',options.BoxColor(i,:),...
            boxPlotOptions{:});
        plot([algoDataBox{i}.MeasureAfterTrim],'Color',options.BoxColor(i,:),...
            'linestyle','--',boxPlotOptions{:});

        [legendEntry{end+1},legendEntry{end+2}] = deal(...
            ['Box ',char(options.BoxLabel{i}),'- Before Trimming Operation'],...
            ['Box ',char(options.BoxLabel{i}),'- After Trimming Operation']);

        xAxisMinLimit = min(xAxisMinLimit,algoDataBox{i}.IndexExplorationStart);
        xAxisMaxLimit = max(xAxisMaxLimit,algoDataBox{i}.IndexConsolidationEnd);
        explorationEnd = [explorationEnd,algoDataBox{i}.IndexExplorationEnd];
    end
    for i=1:nComponent
        plot(algoDataComponent{i}.TotalMeasureBeforeTrim,'Color',options.ComponentColor(i,:),...
            componentPlotOptions{:});
        plot(algoDataComponent{i}.TotalMeasureAfterTrim,'Color',options.ComponentColor(i,:),...
            'linestyle','--',componentPlotOptions{:});

        [legendEntry{end+1},legendEntry{end+2}] = deal(...
            ['Component ',char(options.ComponentLabel{i}),'- Before Trimming Operation'],...
            ['Component ',char(options.ComponentLabel{i}),'- After Trimming Operation']);

        xAxisMinLimit = min(xAxisMinLimit,algoDataComponent{i}.IndexExplorationStart);
        xAxisMaxLimit = max(xAxisMaxLimit,algoDataComponent{i}.IndexConsolidationEnd);
        explorationEnd = [explorationEnd,algoDataComponent{i}.IndexExplorationEnd];
    end
    lim = axis;
    for i=1:size(explorationEnd,2)
        plot([explorationEnd(i) explorationEnd(i)],[lim(3) lim(4)],'k-.','HandleVisibility','off');
    end
    axis([xAxisMinLimit xAxisMaxLimit lim(3) lim(4)]);
    grid minor;
    xlabel('Iteration Step');
    ylabel('Measure');
    legend(legendEntry);

    % Volume of Box (normalized by Design Space Volume)
    figureHandle(2) = figure;
    hold all;
    legendEntry = {};
    xAxisMinLimit = inf;
    xAxisMaxLimit = -inf;
    explorationEnd = [];
    for i=1:nBox
        plot([algoDataBox{i}.MeasureBeforeTrimNormalized],'Color',options.BoxColor(i,:),...
            boxPlotOptions{:});
        plot([algoDataBox{i}.MeasureAfterTrimNormalized],'Color',options.BoxColor(i,:),...
            'linestyle','--',boxPlotOptions{:});

        [legendEntry{end+1},legendEntry{end+2}] = deal(...
            ['Box ',char(options.BoxLabel{i}),'- Before Trimming Operation'],...
            ['Box ',char(options.BoxLabel{i}),'- After Trimming Operation']);

        xAxisMinLimit = min(xAxisMinLimit,algoDataBox{i}.IndexExplorationStart);
        xAxisMaxLimit = max(xAxisMaxLimit,algoDataBox{i}.IndexConsolidationEnd);
        explorationEnd = [explorationEnd,algoDataBox{i}.IndexExplorationEnd];
    end
    for i=1:nComponent
        plot(algoDataComponent{i}.TotalMeasureBeforeTrimNormalized,'Color',options.ComponentColor(i,:),...
            componentPlotOptions{:});
        plot(algoDataComponent{i}.TotalMeasureAfterTrimNormalized,'Color',options.ComponentColor(i,:),...
            'linestyle','--',componentPlotOptions{:});

        [legendEntry{end+1},legendEntry{end+2}] = deal(...
            ['Component ',char(options.ComponentLabel{i}),'- Before Trimming Operation'],...
            ['Component ',char(options.ComponentLabel{i}),'- After Trimming Operation']);

        xAxisMinLimit = min(xAxisMinLimit,algoDataComponent{i}.IndexExplorationStart);
        xAxisMaxLimit = max(xAxisMaxLimit,algoDataComponent{i}.IndexConsolidationEnd);
        explorationEnd = [explorationEnd,algoDataComponent{i}.IndexExplorationEnd];
    end
    lim = axis;
    for i=1:size(explorationEnd,2)
        plot([explorationEnd(i) explorationEnd(i)],[lim(3) lim(4)],'k-.','HandleVisibility','off');
    end
    axis([xAxisMinLimit xAxisMaxLimit lim(3) lim(4)]);
    grid minor;
    xlabel('Iteration Step');
    ylabel('Normalized Measure (V/V_{ds})');
    legend(legendEntry);

    % growth rate
    figureHandle(3) = figure;
    hold all;
    legendEntry = {};
    xAxisMinLimit = inf;
    xAxisMaxLimit = -inf;
    explorationEnd = [];
    for i=1:nBox
        plot([algoDataBox{i}.GrowthRate],'Color',options.BoxColor(i,:),boxPlotOptions{:});
        legendEntry{end+1} = ['Box ',char(options.BoxLabel{i})];

        xAxisMinLimit = min(xAxisMinLimit,algoDataBox{i}.IndexExplorationStart);
        xAxisMaxLimit = max(xAxisMaxLimit,algoDataBox{i}.IndexConsolidationEnd);
        explorationEnd = [explorationEnd,algoDataBox{i}.IndexExplorationEnd];
    end
    for i=1:nComponent
        plot(algoDataComponent{i}.GrowthRate,'Color',options.ComponentColor(i,:),componentPlotOptions{:});
        legendEntry{end+1} = ['Component ',char(options.ComponentLabel{i})];

        xAxisMinLimit = min(xAxisMinLimit,algoDataComponent{i}.IndexExplorationStart);
        xAxisMaxLimit = max(xAxisMaxLimit,algoDataComponent{i}.IndexConsolidationEnd);
        explorationEnd = [explorationEnd,algoDataComponent{i}.IndexExplorationEnd];
    end
    lim = axis;
    for i=1:size(explorationEnd,2)
        plot([explorationEnd(i) explorationEnd(i)],[lim(3) lim(4)],'k-.','HandleVisibility','off');
    end
    axis([xAxisMinLimit xAxisMaxLimit lim(3) lim(4)]);
    grid minor;
    xlabel('Iteration Step');
    ylabel('Growth Rate');
    legend(legendEntry);
    
    % Sampled Points
    figureHandle(4) = figure;
    hold all;
    legendEntry = {};
    xAxisMinLimit = inf;
    xAxisMaxLimit = -inf;
    explorationEnd = [];
    for i=1:nBox
        if(algoDataBox{i}.IsUsingRequirementSpaces)
            plot(algoDataBox{i}.NumberAcceptableAndUsefulDesigns,boxPlotOptions{:});
            plot(algoDataBox{i}.NumberAcceptableDesigns,boxPlotOptions{:});
            plot(algoDataBox{i}.NumberUsefulDesigns,boxPlotOptions{:});
            [legendEntry{end+1},legendEntry{end+2},legendEntry{end+3}] = deal(...
                ['Box ',char(options.BoxLabel{i}),'- Number of Acceptable & Useful Samples'],...
                ['Box ',char(options.BoxLabel{i}),'- Number of Acceptable Samples'],...
                ['Box ',char(options.BoxLabel{i}),'- Number of Useful Samples']);
        else
            plot([algoDataBox{i}.NumberAcceptableAndUsefulDesigns],'Color',options.BoxColor(i,:),...
                boxPlotOptions{:});
            legendEntry{end+1} = ['Box ',char(options.BoxLabel{i}),' - Number of Good Samples'];
        end

        plot([algoDataBox{i}.NumberEvaluatedSamples]);
        legendEntry{end+1} = ['Box ',char(options.BoxLabel{i}),' - Sample Size'];

        xAxisMinLimit = min(xAxisMinLimit,algoDataBox{i}.IndexExplorationStart);
        xAxisMaxLimit = max(xAxisMaxLimit,algoDataBox{i}.IndexConsolidationEnd);
        explorationEnd = [explorationEnd,algoDataBox{i}.IndexExplorationEnd];
    end
    for i=1:nComponent
        if(algoDataComponent{i}.IsUsingRequirementSpaces)
            plot(algoDataComponent{i}.NumberAcceptableAndUsefulDesigns,componentPlotOptions{:});
            plot(algoDataComponent{i}.NumberAcceptableDesigns,componentPlotOptions{:});
            plot(algoDataComponent{i}.NumberUsefulDesigns,componentPlotOptions{:});
            [legendEntry{end+1},legendEntry{end+2},legendEntry{end+3}] = deal(...
                ['Component ',char(options.ComponentLabel{i}),'- Number of Acceptable & Useful Samples'],...
                ['Component ',char(options.ComponentLabel{i}),'- Number of Acceptable Samples'],...
                ['Component ',char(options.ComponentLabel{i}),'- Number of Useful Samples']);
        else
            plot([algoDataComponent{i}.NumberAcceptableAndUsefulDesigns],'Color',options.ComponentColor(i,:),...
                componentPlotOptions{:});
            legendEntry{end+1} = ['Component ',char(options.ComponentLabel{i}),' - Number of Good Samples'];
        end

        plot([algoDataComponent{i}.NumberEvaluatedSamples]);
        legendEntry{end+1} = ['Component ',char(options.ComponentLabel{i}),' - Sample Size'];

        xAxisMinLimit = min(xAxisMinLimit,algoDataComponent{i}.IndexExplorationStart);
        xAxisMaxLimit = max(xAxisMaxLimit,algoDataComponent{i}.IndexConsolidationEnd);
        explorationEnd = [explorationEnd,algoDataComponent{i}.IndexExplorationEnd];
    end
    lim = axis;
    for i=1:size(explorationEnd,2)
        plot([explorationEnd(i) explorationEnd(i)],[lim(3) lim(4)],'k-.','HandleVisibility','off');
    end
    axis([xAxisMinLimit xAxisMaxLimit lim(3) lim(4)]);
    grid minor;
    xlabel('Iteration Step');
    ylabel('Number of Sample Points');
    legend(legendEntry);

    % Relative Increases
    figureHandle(5) = figure;
    hold all;
    legendEntry = {};
    xAxisMinLimit = inf;
    xAxisMaxLimit = -inf;
    explorationEnd = [];
    for i=1:nBox
        plot(algoDataBox{i}.RatioMeasureChangeBeforeTrim,'Color',options.BoxColor(i,:),boxPlotOptions{:});
        plot(algoDataBox{i}.SamplePurity,'Color',options.BoxColor(i,:),'linestyle','--',boxPlotOptions{:});

        [legendEntry{end+1},legendEntry{end+2}] = deal(...
            ['Box ',char(options.BoxLabel{i}),'- V_i/V_{i-1}'],...
            ['Box ',char(options.BoxLabel{i}),'- Sample Purity']);

        xAxisMinLimit = min(xAxisMinLimit,algoDataBox{i}.IndexExplorationStart);
        xAxisMaxLimit = max(xAxisMaxLimit,algoDataBox{i}.IndexConsolidationEnd);
        explorationEnd = [explorationEnd,algoDataBox{i}.IndexExplorationEnd];
    end
    for i=1:nComponent
        plot(algoDataComponent{i}.RatioTotalMeasureChangeBeforeTrim,'Color',options.ComponentColor(i,:),...
            componentPlotOptions{:});
        plot(algoDataComponent{i}.SamplePurity,'Color',options.ComponentColor(i,:),...
            'linestyle','--',componentPlotOptions{:});

        [legendEntry{end+1},legendEntry{end+2}] = deal(...
            ['Component ',char(options.ComponentLabel{i}),'- V_i/V_{i-1}'],...
            ['Component ',char(options.ComponentLabel{i}),'- Sample Purity']);

        xAxisMinLimit = min(xAxisMinLimit,algoDataComponent{i}.IndexExplorationStart);
        xAxisMaxLimit = max(xAxisMaxLimit,algoDataComponent{i}.IndexConsolidationEnd);
        explorationEnd = [explorationEnd,algoDataComponent{i}.IndexExplorationEnd];
    end
    lim = axis;
    for i=1:size(explorationEnd,2)
        plot([explorationEnd(i) explorationEnd(i)],[lim(3) lim(4)],'k-.','HandleVisibility','off');
    end
    plot([xAxisMinLimit xAxisMaxLimit],[1 1],'k--');
    axis([xAxisMinLimit xAxisMaxLimit max(0,lim(3)) min(2,lim(4))]);
    grid minor;
    xlabel('Iteration Step');
    legend(legendEntry);

    % Normalized Volume over m/N
    figureHandle(6) = figure;
    hold all;
    legendEntry = {};
    for i=1:nBox
        plot(algoDataBox{i}.SamplePurity,algoDataBox{i}.MeasureBeforeTrimNormalized,'-',...
            'Color',options.BoxColor(i,:),'Marker','.',boxPlotOptions{:});
        legendEntry{end+1} = ['Box ',char(options.BoxLabel{i})];
    end
    for i=1:nComponent
        plot(algoDataComponent{i}.SamplePurity,algoDataComponent{i}.TotalMeasureBeforeTrimNormalized,'-',...
            'Color',options.ComponentColor(i,:),'Marker','.',componentPlotOptions{:});
        legendEntry{end+1} = ['Component ',char(options.ComponentLabel{i})];
    end
    grid minor;
    xlabel('Sample Purity');
    ylabel('Normalized Volume (V/V_{ds})');
    legend(legendEntry);
    
    % 3D plot with purity/size/cost
    figureHandle(7) = figure;
    hold all;
    legendEntry = {};
    for i=1:nBox
        plot3(...
            algoDataBox{i}.TotalFunctionEvaluations,...
            algoDataBox{i}.SamplePurity,...
            algoDataBox{i}.MeasureBeforeTrimNormalized,...
            '-','Color',options.BoxColor(i,:),'Marker','.',boxPlotOptions{:});
        legendEntry{end+1} = ['Box ',char(options.BoxLabel{i})];
    end
    for i=1:nComponent
        plot3(...
            algoDataComponent{i}.TotalFunctionEvaluations,...
            algoDataComponent{i}.SamplePurity,...
            algoDataComponent{i}.TotalMeasureBeforeTrimNormalized,...
            '-','Color',options.ComponentColor(i,:),'Marker','.',componentPlotOptions{:});
        legendEntry{end+1} = ['Component ',char(options.ComponentLabel{i})];
    end
    grid minor;
    xlabel('Total Number of Function Evaluations');
    ylabel('Sample Purity');
    zlabel('Normalized Volume (V/V_{ds})');
    legend(legendEntry);

    % save if required
    if(~isempty(options.SaveFolder))
        save_print_figure(figureHandle(1),[options.SaveFolder,'Comparison-Metrics-TotalMeasure'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(2),[options.SaveFolder,'Comparison-Metrics-TotalMeasureNormalized'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(3),[options.SaveFolder,'Comparison-Metrics-GrowthRate'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(4),[options.SaveFolder,'Comparison-Metrics-NumberLabelSamplePoints'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(5),[options.SaveFolder,'Comparison-Metrics-TotalMeasureRelativeRatio'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(6),[options.SaveFolder,'Comparison-Metrics-NormalizedMeasureSample'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(7),[options.SaveFolder,'Comparison-Metrics-TotalEvaluationsPurityNormalizedMeasure'],options.SaveFigureOptions{:});

        if(options.CloseFigureAfterSaving)
            close(figureHandle);
        end
    end

    if(nargout<1)
        clear figureHandle
    end
end