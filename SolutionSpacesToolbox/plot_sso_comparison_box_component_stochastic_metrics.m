function plot_sso_comparison_box_component_stochastic_metrics(algoDataBox,algoDataComponent,varargin)
%PLOT_SSO_COMPARISON_BOX_COMPONENT_STOCHASTIC_METRICS Compare stochastic results
%   PLOT_SSO_COMPARISON_BOX_COMPONENT_STOCHASTIC_METRICS automatically generates
%   a series of plots to analyze the performance of a stochastic algorithm that  
%   finds (1) box-based solutions and (2) component-based solutions. The 
%   function creates seven figures that compare both solution approaches along 
%   iteration steps.
%
%   These plots are:
%       (1) Iteration Step vs. Measure
%       (2) Iteration Step vs. Normalized Measure (V/V_{ds})
%       (3) Iteration Step vs. Growth Rate
%       (4) Iteration Step vs. Number of Sampled Points
%       (5) Iteration Step vs. Measure Ratio (V_i/V_{i-1}) and Sample Purity
%       (6) Sample Purity vs. Normalized Measure (V/V_{ds})
%       (7) Total Function Evaluations vs. Sample Purity vs. Normalized Measure
%
%   PLOT_SSO_COMPARISON_BOX_COMPONENT_STOCHASTIC_METRICS(ALGODATABOX,
%   ALGODATACOMPONENT) expects two structures:
%       - ALGODATABOX : (1,nBox) cell array of structs, each struct holding 
%       iterative metrics for a particular box-based approach.
%       - ALGODATACOMPONENT: (1,nComponent) cell array of structs, each struct 
%       holding iterative metrics for a particular component-based approach.
%
%   PLOT_SSO_COMPARISON_BOX_COMPONENT_STOCHASTIC_METRICS(...,'NAME',VALUE,...) 
%   allows additional name-value pairs to customize the output. Key options are:
%       - 'SaveFolder' : if specified, all figures are automatically saved in 
%       this folder using 'save_print_figure'. Default is [] (no saving).
%       - 'CloseFigureAfterSaving' : if true each figure is closed after saving. 
%       Default is false.
%       - 'SaveFigureOptions' : Cell array of options passed to 
%       'save_print_figure'.
%       - 'BoxLabel','ComponentLabel' : cell arrays of labels for each 
%       box/component.
%       - 'BoxColor','ComponentColor' : arrays of colors for plotting each 
%       box/component.
%       - 'GeneralPlotOptions','BoxPlotOptions','ComponentPlotOptions' : cell  
%       arrays of additional name-value pairs for customizing plot appearance.
%
%   FIGUREHANDLE = PLOT_SSO_COMPARISON_BOX_COMPONENT_STOCHASTIC_METRICS(...)
%   returns a (1,7) array of figure handles in FIGUREHANDLE. If no output 
%   variable is specified, no figure handles are returned.
%
%   Inputs:
%       - ALGODATABOX : (1,nBox) cell array of structs
%       - ALGODATACOMPONENT : (1,nComponent) cell array of structs
%       - 'SaveFolder' : char or string
%       - 'CloseFigureAfterSaving' : logical 
%       - 'SaveFigureOptions' : cell
%       - 'BoxLabel','ComponentLabel' : cell of char/string
%       - 'BoxColor','ComponentColor' : (nBox,3) and (nComponent,3) double 
%       (RGB color arrays)
%       - 'GeneralPlotOptions','BoxPlotOptions','ComponentPlotOptions' : cell 
%       arrays of plot options
%
%   Outputs:
%       - FIGUREHANDLE : (1,7) Figure
%
%   See also: sso_box_stochastic, postprocess_sso_box_stochastic, 
%   sso_component_stochastic, postprocess_component_stochastic, 
%   save_print_figure.
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
        save_print_figure(figureHandle(1),[options.SaveFolder,'TotalMeasure'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(2),[options.SaveFolder,'TotalMeasureNormalized'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(3),[options.SaveFolder,'GrowthRate'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(4),[options.SaveFolder,'NumberLabelSamplePoints'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(5),[options.SaveFolder,'TotalMeasureRelativeRatio'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(6),[options.SaveFolder,'NormalizedMeasureSample'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(7),[options.SaveFolder,'TotalEvaluationsPurityNormalizedMeasure'],options.SaveFigureOptions{:});

        if(options.CloseFigureAfterSaving)
            close(figureHandle);
        end
    end

    if(nargout<1)
        clear figureHandle
    end
end