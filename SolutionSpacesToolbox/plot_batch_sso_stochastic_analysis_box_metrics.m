function plot_batch_sso_stochastic_analysis_box_metrics(batchOptions,algoDataBox,varargin)
%PLOT_BATCH_SSO_STOCHASTIC_ANALYSIS_BOX_METRICS Visualize batch results
%   PLOT_BATCH_SSO_STOCHASTIC_ANALYSIS_BOX_METRICS produces plots to allow
%   one to visualize the results of a batch analysis for component SSO.
%   The following plots are produced:
%       - Individual plots of each test separately for their metrics 
%       ('plot_sso_component_stochastic_metrics').
%       - Collated metric plots with all the components and the three most 
%       important metrics: total number of function evaluations, purity and
%       measure.
%       - Sensitivity plots of the evaluations/purity/measure for each option, 
%       including all tests.
%       - Sensitivity plots of the evaluations/purity/measure for each option,
%       but exclusively including tests were that option was the only value 
%       changed from the reference.
%
%   PLOT_BATCH_SSO_STOCHASTIC_ANALYSIS_BOX_METRICS(BATCHOPTIONS,
%   ALGODATACOMPONENT) produces the plots as described above.
%
%   PLOT_BATCH_SSO_STOCHASTIC_ANALYSIS_BOX_METRICS(...NAME,VALUE,...)
%   allows the specification of additional options. These are:
%       - 'SaveFolder' : folder were pictures can be saved (if empty, figures
%       are not saved). Default is empty.
%       - 'CloseFigureAfterSaving' : flag as to whether pictures should be 
%       closed after they are saved. Default: true.
%       - 'SaveFigureOptions' : options for 'save_print_figure'. Default is 
%       empty.
%       - 'PlotOptionsIndividualComponentMetrics' : options for 
%       'plot_component_stochastic_metrics'. Default is empty, but the previous
%       options ('SaveFolder','CloseFigureAfterSaving','SaveFigureOptions') are
%       included.
%       - 'PlotOptionsCollatedComponentMetrics' : options for 'plot'/'plot3' of
%       collated plots. Default: {'Linestyle','--','Marker','.'}.
%       - 'PlotOptionsScatterMetrics' : options for 'plot' when checking
%       the final values of the metrics only. Default: {'Linestyle','none',
%       'Marker','.','MarkerSize',20}.
%       - 'PlotOptionsSensitivityAll' : options for 'plot' when checking
%       the sensitivity of the metrics with all tests. Default: {'Linestyle',
%       'none','Marker','.','MarkerSize',20}.
%       - 'PlotOptionsSensitivityIndividual' : options for 'plot' when checking
%       the sensitivity of the metrics with only tests where each option is
%       the only parameter changed. Default: {'Linestyle','--','Marker','.'};
%       - 'BoxMeasureComparison' : measure of a reference box; if given, all
%       measures will be normalized relative to this.
%       - 'BoxAlgoData' : algorithm data for the box stochastic algorithm; if 
%       given, this information will be plotted together with the component
%       algorithm data for comparisons. Additionally, the measures of the
%       component SSO tests will be normalized with the final measure from the
%       box.
%       
%
%   Input:
%       - BATCHOPTIONS : structure
%           -- ID : char OR string
%           -- ReferenceID : char OR string
%           -- Options : structure
%       - ALGODATACOMPONENT : (nTest,1) structure
%       - 'SaveFolder' : char OR string
%       - 'CloseFigureAfterSaving' : logical
%       - 'SaveFigureOptions' : (1,nOptionFigure) cell
%       - 'PlotOptionsIndividualComponentMetrics' : (1,nOptionPlot1) cell
%       - 'PlotOptionsCollatedComponentMetrics' : (1,nOptionPlot2) cell
%       - 'PlotOptionsScatterMetrics' : (1,nOptionPlot3) cell
%       - 'PlotOptionsSensitivityAll' : (1,nOptionPlot4) cell
%       - 'PlotOptionsSensitivityIndividual' : (1,nOptionPlot5) cell
%       - 'BoxMeasureComparison' : double
%       - 'BoxAlgoData' : structure
%
%   See also batch_sso_stochastic_analysis.
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

	parser = inputParser;
    parser.addParameter('SaveFolder',[]);
    parser.addParameter('CloseFigureAfterSaving',true);
    parser.addParameter('SaveFigureOptions',{});
	parser.addParameter('PlotOptionsIndividualComponentMetrics',{});
    parser.addParameter('PlotOptionsCollatedComponentMetrics',{});
    parser.addParameter('PlotOptionsScatterMetrics',{});
    parser.addParameter('PlotOptionsSensitivityAll',{});
    parser.addParameter('PlotOptionsSensitivityIndividual',{});
	parser.parse(varargin{:});
	options = parser.Results;


    defaultPlotOptionsIndividualComponentMetrics = {};
    [~,plotOptionsComponentMetrics] = merge_name_value_pair_argument(...
        defaultPlotOptionsIndividualComponentMetrics,options.PlotOptionsIndividualComponentMetrics);

    defaultPlotOptionsCollatedComponentMetrics = {'Linestyle','-','Marker','.'};
    [~,plotOptionsCollatedComponentMetrics] = merge_name_value_pair_argument(...
        defaultPlotOptionsCollatedComponentMetrics,options.PlotOptionsCollatedComponentMetrics);

    defaultPlotOptionsScatterMetrics = {'Linestyle','none','Marker','.','MarkerSize',20};
    [~,plotOptionsScatterMetrics] = merge_name_value_pair_argument(...
        defaultPlotOptionsScatterMetrics,options.PlotOptionsScatterMetrics);

    defaultPlotOptionsSensitivityAll = {'Linestyle','none','Marker','.','MarkerSize',20};
    [~,plotOptionsSensitivityAll] = merge_name_value_pair_argument(...
        defaultPlotOptionsSensitivityAll,options.PlotOptionsSensitivityAll);

    defaultPlotOptionsSensitivityIndividual = {'Linestyle','--','Marker','.','MarkerSize',20};
    [~,plotOptionsSensitivityIndividual] = merge_name_value_pair_argument(...
        defaultPlotOptionsSensitivityIndividual,options.PlotOptionsSensitivityIndividual);


	%% create performance metric plots for each result independently
    nTest = length(algoDataBox);
    for i=1:nTest
    	% create folder to save data
    	testFolder = options.SaveFolder;
    	if(~isempty(testFolder))
	        testFolder = [testFolder,batchOptions(i).ID,'/'];
	        mkdir(testFolder);
	    end
        
        plot_sso_box_stochastic_metrics(...
        	algoDataBox(i),...
        	'SaveFolder',testFolder,...
        	'CloseFigureAfterSaving',options.CloseFigureAfterSaving,...
            'SaveFigureOptions',options.SaveFigureOptions,...
        	plotOptionsComponentMetrics{:});
    end

	%% prepare legend
    for j=1:nTest
        legendEntry{j} = sprintf('Box Solution %s',batchOptions(j).ID);
    end


    %%
    figure;
    hold on;
    for j=1:nTest
        plot(algoDataBox(j).SamplePurity,...
            algoDataBox(j).MeasureBeforeTrim,...
            plotOptionsCollatedComponentMetrics{:});
    end
    grid minor;
    xlabel('Purity - Acceptable/Total Samples (m/N)');
    ylabel('Measure (V)');
    legend(legendEntry);
    if(~isempty(options.SaveFolder))
        save_print_figure(gcf,[options.SaveFolder,'PerformanceComparison-Purity-Measure'],options.SaveFigureOptions{:});
        if(options.CloseFigureAfterSaving)
            close(gcf);
        end
    end
    
    
    %%
    figure;
    hold on;
    for j=1:nTest
        plot(algoDataBox(j).TotalFunctionEvaluations,...
            algoDataBox(j).MeasureBeforeTrim,...
            plotOptionsCollatedComponentMetrics{:});
    end
    grid minor;
    xlabel('Total Number of Function Evaluations');
    ylabel('Measure (V)');
    legend(legendEntry);
    if(~isempty(options.SaveFolder))
        save_print_figure(gcf,[options.SaveFolder,'PerformanceComparison-TotalEvaluations-Measure'],options.SaveFigureOptions{:});
        if(options.CloseFigureAfterSaving)
            close(gcf);
        end
    end
    
    
    %%
    figure;
    hold on;
    for j=1:nTest
        plot(algoDataBox(j).SamplePurity(end),...
            algoDataBox(j).MeasureBeforeTrim(end),...
            plotOptionsScatterMetrics{:});
    end
    grid minor;
    xlabel('Purity');
    ylabel('Measure (V)');
    legend(legendEntry,'location','southwest');
    if(~isempty(options.SaveFolder))
        save_print_figure(gcf,[options.SaveFolder,'PerformanceComparisonFinal-Purity-Measure'],options.SaveFigureOptions{:});
        if(options.CloseFigureAfterSaving)
            close(gcf);
        end
    end
    
    
    %%
    figure;
    hold on;
    for j=1:nTest
        plot3(algoDataBox(j).TotalFunctionEvaluations,...
            algoDataBox(j).SamplePurity,...
            algoDataBox(j).MeasureBeforeTrim,...
            plotOptionsCollatedComponentMetrics{:});
    end
    grid minor;
    xlabel('Total Number of Function Evaluations');
    ylabel('Purity');
    zlabel('Measure (V)');
    legend(legendEntry);
    if(~isempty(options.SaveFolder))
        save_print_figure(gcf,[options.SaveFolder,'PerformanceComparison-TotalEvaluations-Purity-Volume'],options.SaveFigureOptions{:});
        if(options.CloseFigureAfterSaving)
            close(gcf);
        end
    end


    %% make plots of changes of individual numerical options
    idEntry = {batchOptions.ID}';
    referenceIdEntry = {batchOptions.ReferenceID}'; 
    nameOption = fieldnames(batchOptions(1).Options);
    optionsTable = struct2table([batchOptions.Options]);
    optionsIndex = (1:length(nameOption))';

    allReference = unique(referenceIdEntry);
    for i=1:length(allReference)
        referenceEntry = find(strcmpi(idEntry,allReference{i}));
        derivedEntry = find(strcmpi(referenceIdEntry,allReference{i}));

        if(isempty(referenceEntry))
            continue;
        end

        for j=1:length(optionsIndex)
            if(~isnumeric(optionsTable{referenceEntry,j}))
                continue;
            end

            optionSaveFolder = [options.SaveFolder,'OptionSensitivity/',nameOption{j},'/',allReference{i},'/'];
            mkdir(optionSaveFolder);

            %
            allRelatedToReference = [referenceEntry;derivedEntry];
            optionValue = [optionsTable{allRelatedToReference,j}];
            [optionValue,sortedIndex] = sort(optionValue);
            allRelatedToReference = allRelatedToReference(sortedIndex);
            nRelated = length(allRelatedToReference);

            %
            measureValues = nan(nRelated,1);
            for k=1:nRelated
                measureValues(k) = algoDataBox(allRelatedToReference(k)).MeasureBeforeTrim(end);
            end
            figure;
            hold on;
            plot(optionValue,measureValues,plotOptionsSensitivityAll{:});
            grid minor;
            xlabel(nameOption{j});
            ylabel('Measure (V)');
            title(idEntry{referenceEntry});
            if(~isempty(options.SaveFolder))
                save_print_figure(gcf,[optionSaveFolder,'SensitivityAllMeasure'],options.SaveFigureOptions{:});
                if(options.CloseFigureAfterSaving)
                    close(gcf);
                end
            end
            
            %
            totalFunctionEvaluationValues = nan(nRelated,1);
            for k=1:nRelated
                totalFunctionEvaluationValues(k) = algoDataBox(allRelatedToReference(k)).TotalFunctionEvaluations(end);
            end
            figure;
            hold on;
            plot(optionValue,totalFunctionEvaluationValues,plotOptionsSensitivityAll{:});
            grid minor;
            xlabel(nameOption{j});
            ylabel('Total Number of Function Evaluations');
            title(idEntry{referenceEntry});
            if(~isempty(options.SaveFolder))
                save_print_figure(gcf,[optionSaveFolder,'SensitivityAllTotalEvaluations'],options.SaveFigureOptions{:});
                if(options.CloseFigureAfterSaving)
                    close(gcf);
                end
            end
            
            %
            purityValues = nan(nRelated,1);
            for k=1:nRelated
                purityValues(k) = algoDataBox(allRelatedToReference(k)).SamplePurity(end);
            end
            figure;
            hold on;
            plot(optionValue,purityValues,plotOptionsSensitivityAll{:});
            grid minor;
            xlabel(nameOption{j});
            ylabel('Purity');
            title(idEntry{referenceEntry});
            if(~isempty(options.SaveFolder))
                save_print_figure(gcf,[optionSaveFolder,'SensitivityAllPurity'],options.SaveFigureOptions{:});
                if(options.CloseFigureAfterSaving)
                    close(gcf);
                end
            end

            otherOptions = optionsIndex;
            otherOptions(j) = [];
            
            % from the derived entries, find all the ones that have only the 
            % value of this particular option different
            nDerivedEntry = length(derivedEntry);
            selectForPlot = false(nDerivedEntry,1);
            for k=1:nDerivedEntry
                differentValueThisOption = ~isequal(optionsTable(referenceEntry,j),...
                    optionsTable(derivedEntry(k),j));
                sameValueOtherOption = isequal(optionsTable(referenceEntry,otherOptions),...
                    optionsTable(derivedEntry(k),otherOptions));
                selectForPlot(k) = differentValueThisOption && sameValueOtherOption;
            end

            % get relevant values for this option
            allSelectForPlotIndex = [referenceEntry;derivedEntry(selectForPlot)];
            optionValue = [optionsTable{allSelectForPlotIndex,j}];
            [optionValue,sortedIndex] = sort(optionValue);
            allSelectForPlotIndex = allSelectForPlotIndex(sortedIndex);
            nPlotPoint = length(allSelectForPlotIndex);

            %
            measureValues = nan(nPlotPoint,1);
            for k=1:nPlotPoint
                measureValues(k) = algoDataBox(allSelectForPlotIndex(k)).MeasureBeforeTrim(end);
            end
            figure;
            hold on;
            plot(optionValue,measureValues,plotOptionsSensitivityIndividual{:});
            grid minor;
            xlabel(nameOption{j});
            ylabel('Measure (V)');
            title(idEntry{referenceEntry});
            if(~isempty(options.SaveFolder))
                save_print_figure(gcf,[optionSaveFolder,'SensitivityIndividualMeasure'],options.SaveFigureOptions{:});
                if(options.CloseFigureAfterSaving)
                    close(gcf);
                end
            end
            
            %
            totalFunctionEvaluationValues = nan(nPlotPoint,1);
            for k=1:nPlotPoint
                totalFunctionEvaluationValues(k) = algoDataBox(allSelectForPlotIndex(k)).TotalFunctionEvaluations(end);
            end
            figure;
            hold on;
            plot(optionValue,totalFunctionEvaluationValues,plotOptionsSensitivityIndividual{:});
            grid minor;
            xlabel(nameOption{j});
            ylabel('Total Number of Function Evaluations');
            title(idEntry{referenceEntry});
            if(~isempty(options.SaveFolder))
                save_print_figure(gcf,[optionSaveFolder,'SensitivityIndividualTotalEvaluations'],options.SaveFigureOptions{:});
                if(options.CloseFigureAfterSaving)
                    close(gcf);
                end
            end
            
            %
            purityValues = nan(nPlotPoint,1);
            for k=1:nPlotPoint
                purityValues(k) = algoDataBox(allSelectForPlotIndex(k)).SamplePurity(end);
            end
            figure;
            hold on;
            plot(optionValue,purityValues,plotOptionsSensitivityIndividual{:});
            grid minor;
            xlabel(nameOption{j});
            ylabel('Purity');
            title(idEntry{referenceEntry});
            if(~isempty(options.SaveFolder))
                save_print_figure(gcf,[optionSaveFolder,'SensitivityIndividualPurity'],options.SaveFigureOptions{:});
                if(options.CloseFigureAfterSaving)
                    close(gcf);
                end
            end
        end
    end
end