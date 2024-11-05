function plot_sso_comparison_box_component_stochastic_metrics(algoDataBox,algoDataComponent,varargin)
	parser = inputParser;
    parser.addParameter('SaveFolder',[]);
    parser.addParameter('CloseFigureAfterSaving',false);
    parser.addParameter('SaveFigureOptions',{});
    parser.parse(varargin{:});
    options = parser.Results;

    %% Plots
    % Volume of Box
    figureHandle(1) = figure;
    plot([algoDataBox.MeasureBeforeTrim]);
    hold on;
    plot(algoDataComponent.TotalMeasureBeforeTrim);
    plot([algoDataBox.MeasureAfterTrim]);
    plot(algoDataComponent.TotalMeasureAfterTrim);
    lim = axis;
    plot([algoDataBox.IndexExplorationEnd algoDataBox.IndexExplorationEnd],[lim(3) lim(4)],'k-.');
    plot([algoDataComponent.IndexExplorationEnd algoDataComponent.IndexExplorationEnd],[lim(3) lim(4)],'k-.');
    axis([min(algoDataBox.IndexExplorationStart,algoDataComponent.IndexExplorationStart) max(algoDataBox.IndexConsolidationEnd,algoDataComponent.IndexConsolidationEnd) lim(3) lim(4)]);
    grid minor;
    xlabel('Iteration Step');
    ylabel('Measure');
    legend({'Box - Before Trimming Operation','Component - Before Trimming Operation',...
    	'Box - After Trimming Operation','Component - After Trimming Operation'});

    % Volume of Box (normalized by Design Space Volume)
    figureHandle(2) = figure;
    plot(algoDataBox.MeasureBeforeTrimNormalized);
    hold on;
    plot(algoDataComponent.TotalMeasureBeforeTrimNormalized);
    plot(algoDataBox.MeasureAfterTrimNormalized);
    plot(algoDataComponent.TotalMeasureAfterTrimNormalized);
    lim = axis;
    plot([algoDataBox.IndexExplorationEnd algoDataBox.IndexExplorationEnd],[lim(3) lim(4)],'k-.');
    plot([algoDataComponent.IndexExplorationEnd algoDataComponent.IndexExplorationEnd],[lim(3) lim(4)],'k-.');
    axis([min(algoDataBox.IndexExplorationStart,algoDataComponent.IndexExplorationStart) max(algoDataBox.IndexConsolidationEnd,algoDataComponent.IndexConsolidationEnd) lim(3) lim(4)]);
    grid minor;
    xlabel('Iteration Step');
    ylabel('Normalized Measure (V/V_{ds})');
    legend({'Box - Before Trimming Operation','Component - Before Trimming Operation',...
    	'Box - After Trimming Operation','Component - After Trimming Operation'});

    % growth rate
    figureHandle(3) = figure;
    plot(algoDataBox.GrowthRate);
    hold on;
    plot(algoDataComponent.GrowthRate);
    lim = axis;
    plot([algoDataBox.IndexExplorationEnd algoDataBox.IndexExplorationEnd],[lim(3) lim(4)],'k-.');
    plot([algoDataComponent.IndexExplorationEnd algoDataComponent.IndexExplorationEnd],[lim(3) lim(4)],'k-.');
    axis([min(algoDataBox.IndexExplorationStart,algoDataComponent.IndexExplorationStart) max(algoDataBox.IndexConsolidationEnd,algoDataComponent.IndexConsolidationEnd) lim(3) lim(4)]);
    grid minor;
    xlabel('Iteration Step');
    ylabel('Growth Rate');
    legend({'Box','Component'});
    
    % Sampled Points
    figureHandle(4) = figure;
    plot(algoDataBox.NumberAcceptableAndUsefulDesigns);
    hold on;
    if(algoDataBox.IsUsingRequirementSpaces)
        plot(algoDataBox.NumberAcceptableDesigns);
        plot(algoDataBox.NumberUsefulDesigns);
    end
    plot(algoDataBox.NumberEvaluatedSamples);
    plot(algoDataComponent.NumberAcceptableAndUsefulDesigns);
    if(algoDataComponent.IsUsingRequirementSpaces)
        plot(algoDataComponent.NumberAcceptableDesigns);
        plot(algoDataComponent.NumberUsefulDesigns);
    end
    plot(algoDataComponent.NumberEvaluatedSamples);
    lim = axis;
    plot([algoDataBox.IndexExplorationEnd algoDataBox.IndexExplorationEnd],[lim(3) lim(4)],'k-.');
    plot([algoDataComponent.IndexExplorationEnd algoDataComponent.IndexExplorationEnd],[lim(3) lim(4)],'k-.');
    axis([min(algoDataBox.IndexExplorationStart,algoDataComponent.IndexExplorationStart) max(algoDataBox.IndexConsolidationEnd,algoDataComponent.IndexConsolidationEnd) lim(3) lim(4)]);
    grid minor;
    xlabel('Iteration Step');
    ylabel('Number of Sample Points');
    if(algoDataBox.IsUsingRequirementSpaces)
        legend({...
        	'Box - Number of Acceptable & Useful Samples',...
        	'Box - Number of Acceptable Samples',...
        	'Box - Number of Useful Samples',...
        	'Box - Sample Size',...
        	'Component - Number of Acceptable & Useful Samples',...
        	'Component - Number of Acceptable Samples',...
        	'Component - Number of Useful Samples',...
        	'Component - Sample Size'});
    else
        legend({...
        	'Box - Number of Good Samples',...
        	'Box - Sample Size',...
        	'Component - Number of Good Samples',...
        	'Component - Sample Size'});
    end

    % Relative Increases
    figureHandle(5) = figure;
    plot(algoDataBox.RatioMeasureChangeBeforeTrim);
    hold on;
    plot(algoDataBox.SamplePurity);
    plot(algoDataComponent.RatioTotalMeasureChangeBeforeTrim);
    plot(algoDataComponent.SamplePurity);
    lim = axis;
    plot([algoDataBox.IndexExplorationEnd algoDataBox.IndexExplorationEnd],[lim(3) lim(4)],'k-.');
    plot([algoDataComponent.IndexExplorationEnd algoDataComponent.IndexExplorationEnd],[lim(3) lim(4)],'k-.');
    plot([algoDataBox.IndexExplorationStart algoDataBox.IndexConsolidationEnd],[1 1],'k--');
    axis([...
        min(algoDataBox.IndexExplorationStart,algoDataComponent.IndexExplorationStart) ... % x_min
        max(algoDataBox.IndexConsolidationEnd,algoDataComponent.IndexConsolidationEnd) ... % x_max
        max(0,lim(3)) ... % y_min
        min(2,lim(4))]); % y_max
    grid minor;
    xlabel('Iteration Step');
    legend({'Box - V_i/V_{i-1}','Box - Sample Purity','Component - V_i/V_{i-1}','Component - Sample Purity'})

    % Normalized Volume over m/N
    figureHandle(6) = figure;
    plot(algoDataBox.SamplePurity,...
        algoDataBox.MeasureBeforeTrimNormalized,...
        '.-');
    hold on;
    plot(algoDataComponent.SamplePurity,algoDataComponent.TotalMeasureBeforeTrimNormalized,'.-');
    grid minor;
    xlabel('Sample Purity');
    ylabel('Normalized Volume (V/V_{ds})');
    legend({'Box','Component'});
    
    % 3D plot with purity/size/cost
    figureHandle(7) = figure;
    plot3(algoDataBox.TotalFunctionEvaluations,...
        algoDataBox.SamplePurity,...
        algoDataBox.MeasureBeforeTrimNormalized,...
        '-','Marker','.');
    hold on;
    plot3(algoDataComponent.TotalFunctionEvaluations,...
        algoDataComponent.SamplePurity,...
        algoDataComponent.TotalMeasureBeforeTrimNormalized,...
        '-','Marker','.');
    grid minor;
    xlabel('Total Number of Function Evaluations');
    ylabel('Sample Purity');
    zlabel('Normalized Volume (V/V_{ds})');
    legend({'Box','Component'});

    % save if required
    if(~isempty(options.SaveFolder))
        save_print_figure(figureHandle(1),[options.SaveFolder,'Comparsion-BoxComponent-Metrics-TotalMeasure'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(2),[options.SaveFolder,'Comparsion-BoxComponent-Metrics-TotalMeasureNormalized'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(3),[options.SaveFolder,'Comparsion-BoxComponent-Metrics-GrowthRate'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(4),[options.SaveFolder,'Comparsion-BoxComponent-Metrics-NumberLabelSamplePoints'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(5),[options.SaveFolder,'Comparsion-BoxComponent-Metrics-TotalMeasureRelativeRatio'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(6),[options.SaveFolder,'Comparsion-BoxComponent-Metrics-NormalizedMeasureSample'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(7),[options.SaveFolder,'Comparsion-BoxComponent-Metrics-TotalEvaluationsPurityNormalizedMeasure'],options.SaveFigureOptions{:});

        if(options.CloseFigureAfterSaving)
            close(figureHandle);
        end
    end

    if(nargout<1)
        clear figureHandle
    end
end