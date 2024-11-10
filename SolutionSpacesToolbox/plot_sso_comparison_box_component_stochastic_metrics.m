function plot_sso_comparison_box_component_stochastic_metrics(algoDataBox,algoDataComponent,varargin)
	parser = inputParser;
    parser.addParameter('SaveFolder',[]);
    parser.addParameter('CloseFigureAfterSaving',false);
    parser.addParameter('SaveFigureOptions',{});
    parser.addParameter('BoxLabel',{});
    parser.addParameter('ComponentLabel',{});
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

    %% Plots
    % Volume of Box
    figureHandle(1) = figure;
    hold all;
    legendEntry = {};
    xAxisMinLimit = inf;
    xAxisMaxLimit = -inf;
    explorationEnd = [];
    for i=1:nBox
        plot([algoDataBox{i}.MeasureBeforeTrim]);
        plot([algoDataBox{i}.MeasureAfterTrim]);

        [legendEntry{end+1},legendEntry{end+2}] = deal(...
            ['Box ',char(options.BoxLabel{i}),'- Before Trimming Operation'],...
            ['Box ',char(options.BoxLabel{i}),'- After Trimming Operation']);

        xAxisMinLimit = min(xAxisMinLimit,algoDataBox{i}.IndexExplorationStart);
        xAxisMaxLimit = max(xAxisMaxLimit,algoDataBox{i}.IndexConsolidationEnd);
        explorationEnd = [explorationEnd,algoDataBox{i}.IndexExplorationEnd];
    end
    for i=1:nComponent
        plot(algoDataComponent{i}.TotalMeasureBeforeTrim);
        plot(algoDataComponent{i}.TotalMeasureAfterTrim);

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
        plot([algoDataBox{i}.MeasureBeforeTrimNormalized]);
        plot([algoDataBox{i}.MeasureAfterTrimNormalized]);

        [legendEntry{end+1},legendEntry{end+2}] = deal(...
            ['Box ',char(options.BoxLabel{i}),'- Before Trimming Operation'],...
            ['Box ',char(options.BoxLabel{i}),'- After Trimming Operation']);

        xAxisMinLimit = min(xAxisMinLimit,algoDataBox{i}.IndexExplorationStart);
        xAxisMaxLimit = max(xAxisMaxLimit,algoDataBox{i}.IndexConsolidationEnd);
        explorationEnd = [explorationEnd,algoDataBox{i}.IndexExplorationEnd];
    end
    for i=1:nComponent
        plot(algoDataComponent{i}.TotalMeasureBeforeTrimNormalized);
        plot(algoDataComponent{i}.TotalMeasureAfterTrimNormalized);

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
        plot([algoDataBox{i}.GrowthRate]);
        legendEntry{end+1} = ['Box ',char(options.BoxLabel{i})];

        xAxisMinLimit = min(xAxisMinLimit,algoDataBox{i}.IndexExplorationStart);
        xAxisMaxLimit = max(xAxisMaxLimit,algoDataBox{i}.IndexConsolidationEnd);
        explorationEnd = [explorationEnd,algoDataBox{i}.IndexExplorationEnd];
    end
    for i=1:nComponent
        plot(algoDataComponent{i}.GrowthRate);
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
            plot(algoDataBox{i}.NumberAcceptableAndUsefulDesigns);
            plot(algoDataBox{i}.NumberAcceptableDesigns);
            plot(algoDataBox{i}.NumberUsefulDesigns);
            [legendEntry{end+1},legendEntry{end+2},legendEntry{end+3}] = deal(...
                ['Box ',char(options.BoxLabel{i}),'- Number of Acceptable & Useful Samples'],...
                ['Box ',char(options.BoxLabel{i}),'- Number of Acceptable Samples'],...
                ['Box ',char(options.BoxLabel{i}),'- Number of Useful Samples']);
        else
            plot([algoDataBox{i}.NumberAcceptableAndUsefulDesigns]);
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
            plot(algoDataComponent{i}.NumberAcceptableAndUsefulDesigns);
            plot(algoDataComponent{i}.NumberAcceptableDesigns);
            plot(algoDataComponent{i}.NumberUsefulDesigns);
            [legendEntry{end+1},legendEntry{end+2},legendEntry{end+3}] = deal(...
                ['Component ',char(options.ComponentLabel{i}),'- Number of Acceptable & Useful Samples'],...
                ['Component ',char(options.ComponentLabel{i}),'- Number of Acceptable Samples'],...
                ['Component ',char(options.ComponentLabel{i}),'- Number of Useful Samples']);
        else
            plot([algoDataComponent{i}.NumberAcceptableAndUsefulDesigns]);
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
        plot(algoDataBox{i}.RatioMeasureChangeBeforeTrim);
        plot(algoDataBox{i}.SamplePurity);

        [legendEntry{end+1},legendEntry{end+2}] = deal(...
            ['Box ',char(options.BoxLabel{i}),'- V_i/V_{i-1}'],...
            ['Box ',char(options.BoxLabel{i}),'- Sample Purity']);

        xAxisMinLimit = min(xAxisMinLimit,algoDataBox{i}.IndexExplorationStart);
        xAxisMaxLimit = max(xAxisMaxLimit,algoDataBox{i}.IndexConsolidationEnd);
        explorationEnd = [explorationEnd,algoDataBox{i}.IndexExplorationEnd];
    end
    for i=1:nComponent
        plot(algoDataComponent{i}.RatioTotalMeasureChangeBeforeTrim);
        plot(algoDataComponent{i}.SamplePurity);

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
        plot(algoDataBox{i}.SamplePurity,algoDataBox{i}.MeasureBeforeTrimNormalized,'.-');
        legendEntry{end+1} = ['Box ',char(options.BoxLabel{i})];
    end
    for i=1:nComponent
        plot(algoDataComponent{i}.SamplePurity,algoDataComponent{i}.TotalMeasureBeforeTrimNormalized,'.-');
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
            '.-');
        legendEntry{end+1} = ['Box ',char(options.BoxLabel{i})];
    end
    for i=1:nComponent
        plot3(...
            algoDataComponent{i}.TotalFunctionEvaluations,...
            algoDataComponent{i}.SamplePurity,...
            algoDataComponent{i}.TotalMeasureBeforeTrimNormalized,...
            '.-');
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