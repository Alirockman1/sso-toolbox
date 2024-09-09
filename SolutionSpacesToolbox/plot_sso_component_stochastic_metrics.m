function figureHandle = plot_sso_component_stochastic_metrics(algorithmData,varargin)
%PLOT_SSO_COMPONENT_STOCHASTIC_METRICS Plot performance of stochastic component SSO 
%   PLOT_SSO_COMPONENT_STOCHASTIC_METRICS is used to automatically create the most  
%   common necessary plots for analyzing the performance of the stochastic 
%   algorithm used to find the optimal solution boxes. These plots are:
%       - (1) : Iterations x Measure
%       - (2) : Iterations x Normalized Measure (V/V_ds)
%       - (3) : Iterations x Component Normalized Measure (V_i/V_ds,i)
%       - (4) : Iterations x Candidate Space Growth Rate
%       - (5) : Iterations x Number of Acceptable Points
%       - (6) : Iterations x Measure Ratio (current/previous), Sample Purity
%       - (7) : Sample Purity x Normalized Measure (V/V_ds)
%       - (8) : Total Function Evaluations x Sample Purity x Normalized Measure
%       - (9) : Iterations x Ratio of Padding Samples
%   If requirement spaces were used and there are useful/useless designs,
%   those are plotted together.
%
%   PLOT_SSO_COMPONENT_STOCHASTIC_METRICS(ALGORITHMDATA) takes the algorithm 
%   performance data in ALGORITHMDATA and makes figures with the plots specified
%   above.
%
%   PLOT_SSO_COMPONENT_STOCHASTIC_METRICS(...NAME,VALUE,...) also allows the 
%   specification of additional options. These are:
%       - 'SaveFolder' : if it is desired for the pictures to be automatically
%       saved, this must be set to the path of the folder of where to save
%       these plots. If empty, figures are not saved. Default is empty.
%       - 'CloseFigureAfterSaving' : if set to true, the figures will be closed
%       automatically after being saved. Default: false.
%       - 'SaveFigureOptions' : any options for saving the figures; see 
%       'save_print_figure'.
%
%   FIGUREHANDLE = PLOT_SSO_COMPONENT_STOCHASTIC_METRICS(...) also  
%   returns an array with the figure handles in FIGUREHANDLE.
%
%   Inputs:
%       - ALGORITHMDATA : structure
%       - 'SaveFolder' : char OR string
%       - 'CloseFigureAfterSaving' : logical
%       - 'SaveFigureOptions' : cell
%
%   Outputs:
%       - FIGUREHANDLE : (1,9) Figure
%
%   See also sso_component_stochastic, postprocess_component_stochastic, 
%   save_print_figure.
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
    parser.addParameter('CloseFigureAfterSaving',false);
    parser.addParameter('SaveFigureOptions',{});
    parser.parse(varargin{:});
    options = parser.Results;
    
    %% Plots
    % Volume of Candidate Space
    figureHandle(1) = figure;
    plot(algorithmData.TotalMeasureBeforeTrim);
    hold on;
    plot(algorithmData.TotalMeasureAfterTrim);
    lim = axis;
    plot([algorithmData.IndexExplorationEnd algorithmData.IndexExplorationEnd],[lim(3) lim(4)],'k-.');
    axis([algorithmData.IndexExplorationStart algorithmData.IndexConsolidationEnd lim(3) lim(4)]);
    grid minor;
    xlabel('Iteration Step');
    ylabel('Total Measure');
    legend({'Before Trimming Operation','After Trimming Operation'});

    % Volume of Candidate Space (normalized by Design Space Volume)
    figureHandle(2) = figure;
    plot(algorithmData.TotalMeasureBeforeTrimNormalized);
    hold on;
    plot(algorithmData.TotalMeasureAfterTrimNormalized);
    lim = axis;
    plot([algorithmData.IndexExplorationEnd algorithmData.IndexExplorationEnd],[lim(3) lim(4)],'k-.');
    axis([algorithmData.IndexExplorationStart algorithmData.IndexConsolidationEnd lim(3) lim(4)]);
    grid minor;
    xlabel('Iteration Step');
    ylabel('Normalized Total Measure (V/V_{ds})');
    legend({'Before Trimming Operation','After Trimming Operation'});
    
    % Volume of Individual Candidate Space (normalized by Design Space Volume)
    legTextBefore = {};
    legTextAfter = {};
    for i=1:size(algorithmData.ComponentMeasureBeforeTrim,2)
        legTextBefore{end+1} = sprintf('Component %d Before Trimming',i);
        legTextAfter{end+1} = sprintf('Component %d After Trimming',i);
    end
    figureHandle(3) = figure;
    plot(algorithmData.ComponentMeasureBeforeTrimNormalized);
    hold on;
    plot(algorithmData.ComponentMeasureAfterTrimNormalized);
    lim = axis;
    plot([algorithmData.IndexExplorationEnd algorithmData.IndexExplorationEnd],[lim(3) lim(4)],'k-.');
    axis([algorithmData.IndexExplorationStart algorithmData.IndexConsolidationEnd lim(3) lim(4)]);
    grid minor;
    xlabel('Iteration Step');
    ylabel('Normalized Component Measures (V_i/V_{i,ds})');
    legend([legTextBefore,legTextAfter]);

    % growth rate
    figureHandle(4) = figure;
    plot(algorithmData.GrowthRate);
    hold on;
    lim = axis;
    plot([algorithmData.IndexExplorationEnd algorithmData.IndexExplorationEnd],[lim(3) lim(4)],'k-.');
    axis([algorithmData.IndexExplorationStart algorithmData.IndexConsolidationEnd lim(3) lim(4)]);
    grid minor;
    xlabel('Iteration Step');
    ylabel('Growth Rate');
    
    % Sampled Points
    figureHandle(5) = figure;
    plot(algorithmData.NumberAcceptableAndUsefulDesigns);
    hold on;
    if(algorithmData.IsUsingRequirementSpaces)
        plot(algorithmData.NumberAcceptableDesigns);
        plot(algorithmData.NumberUsefulDesigns);
    end
    plot(algorithmData.NumberEvaluatedSamples);
    lim = axis;
    plot([algorithmData.IndexExplorationEnd algorithmData.IndexExplorationEnd],[lim(3) lim(4)],'k-.');
    axis([algorithmData.IndexExplorationStart algorithmData.IndexConsolidationEnd lim(3) lim(4)]);
    grid minor;
    xlabel('Iteration Step');
    ylabel('Number of Sample Points');
    if(algorithmData.IsUsingRequirementSpaces)
        legend({'Number of Acceptable & Useful Samples','Number of Acceptable Samples','Number of Useful Samples','Sample Size'});
    else
        legend({'Number of Good Samples','Sample Size'});
    end

    % Relative Increases
    figureHandle(6) = figure;
    plot(algorithmData.RatioTotalMeasureChangeBeforeTrim);
    hold on;
    plot(algorithmData.SamplePurity);
    lim = axis;
    plot([algorithmData.IndexExplorationEnd algorithmData.IndexExplorationEnd],[lim(3) lim(4)],'k-.');
    plot([algorithmData.IndexExplorationStart algorithmData.IndexConsolidationEnd],[1 1],'k--');
    axis([algorithmData.IndexExplorationStart algorithmData.IndexConsolidationEnd 0 2]);
    grid minor;
    xlabel('Iteration Step');
    legend({'V_i/V_{i-1}','Sample Purity'})

    % Normalized Volume over m/N
    figureHandle(7) = figure;
    plot(algorithmData.SamplePurity,algorithmData.TotalMeasureBeforeTrimNormalized,'.-');
    hold on;
    grid minor;
    xlabel('Sample Purity');
    ylabel('Normalized Volume (V/V_{ds})');

    % 3D plot with purity/size/cost
    figureHandle(8) = figure;
    plot3(algorithmData.TotalFunctionEvaluations,...
        algorithmData.SamplePurity,...
        algorithmData.TotalMeasureBeforeTrimNormalized,...
        '-','Marker','.');
    hold on;
    grid minor;
    xlabel('Total Number of Function Evaluations');
    ylabel('Sample Purity');
    zlabel('Normalized Total Component Volume (V/V_{ds})');
    
    % padding samples
    figureHandle(9) = figure;
    plot(algorithmData.RatioPadding);
    hold on;
    lim = axis;
    plot([algorithmData.IndexExplorationEnd algorithmData.IndexExplorationEnd],[lim(3) lim(4)],'k-.');
    axis([algorithmData.IndexExplorationStart algorithmData.IndexConsolidationEnd lim(3) lim(4)]);
    xlabel('Iteration Step')
    ylabel('Ratio of Padding Samples')
    grid minor;

    % save if required
    if(~isempty(options.SaveFolder))
        save_print_figure(figureHandle(1),[options.SaveFolder,'Component-Metrics-TotalMeasure'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(2),[options.SaveFolder,'Component-Metrics-TotalMeasureNormalized'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(3),[options.SaveFolder,'Component-Metrics-ComponentMeasures'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(4),[options.SaveFolder,'Component-Metrics-GrowthRate'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(5),[options.SaveFolder,'Component-Metrics-NumberLabelSamplePoints'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(6),[options.SaveFolder,'Component-Metrics-TotalMeasureRelativeRatio'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(7),[options.SaveFolder,'Component-Metrics-NormalizedMeasureSample'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(8),[options.SaveFolder,'Component-Metrics-TotalEvaluationsPurityNormalizedMeasure'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(9),[options.SaveFolder,'Component-Metrics-PaddingRatio'],options.SaveFigureOptions{:});

        if(options.CloseFigureAfterSaving)
            close(figureHandle);
        end
    end

    if(nargout<1)
        clear figureHandle
    end
end

