function figureHandle = plot_sso_box_stochastic_metrics(algorithmData,varargin)
%PLOT_SSO_BOX_STOCHASTIC_METRICS Plot performance metrics of stochastic box SSO 
%   PLOT_SSO_BOX_STOCHASTIC_METRICS is used to automatically create the most  
%   common necessary plots for analyzing the performance of the stochastic 
%   algorithm used to find the optimal solution boxes. These plots are:
%       - (1) : Iterations x Measure
%       - (2) : Iterations x Normalized Measure (V/V_ds)
%       - (3) : Iterations x Design Variable Normalized Measure (V_i/V_ds,i)
%       - (4) : Iterations x Candidate Box Growth Rate
%       - (5) : Iterations x Number of Acceptable Points
%       - (6) : Iterations x Measure Ratio (current/previous), Sample Purity
%       - (7) : Sample Purity x Normalized Measure (V/V_ds)
%       - (8) : Total Function Evaluations x Sample Purity x Normalized Measure
%   If requirement spaces were used and there are useful/useless designs,
%   those are plotted together.
%
%   PLOT_SSO_BOX_STOCHASTIC_METRICS(ALGORITHMDATA) takes the algorithm 
%   performance data in ALGORITHMDATA and makes figures with the plots specified
%   above.
%
%   PLOT_SSO_BOX_STOCHASTIC_METRICS(...NAME,VALUE,...) also allows the 
%   specification of additional options. These are:
%       - 'SaveFolder' : if it is desired for the pictures to be automatically
%       saved, this must be set to the path of the folder of where to save
%       these plots. If empty, figures are not saved. Default is empty.
%       - 'CloseFigureAfterSaving' : if set to true, the figures will be closed
%       automatically after being saved. Default: false.
%       - 'SaveFigureOptions' : any options for saving the figures; see 
%       'save_print_figure'.
%
%   FIGUREHANDLE = PLOT_SSO_BOX_STOCHASTIC_METRICS(ALGORITHMDATA) also returns 
%   an array with the figure handles in FIGUREHANDLE.
%
%   Inputs:
%       - ALGORITHMDATA : structure
%       - 'SaveFolder' : char OR string
%       - 'CloseFigureAfterSaving' : logical
%       - 'SaveFigureOptions' : cell
%
%   Outputs:
%       - FIGUREHANDLE : (1,8) Figure
%
%   See also sso_box_stochastic, postprocess_sso_box_stochastic.
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
    % Volume of Box
    figureHandle(1) = figure;
    plot([algorithmData.MeasureBeforeTrim]);
    hold on;
    plot([algorithmData.MeasureAfterTrim]);
    lim = axis;
    plot([algorithmData.IndexExplorationEnd algorithmData.IndexExplorationEnd],[lim(3) lim(4)],'k-.');
    axis([algorithmData.IndexExplorationStart algorithmData.IndexConsolidationEnd lim(3) lim(4)]);
    grid minor;
    xlabel('Iteration Step');
    ylabel('Box Measure');
    legend({'Before Trimming Operation','After Trimming Operation'});

    % Volume of Box (normalized by Design Space Volume)
    figureHandle(2) = figure;
    plot(algorithmData.MeasureBeforeTrimNormalized);
    hold on;
    plot(algorithmData.MeasureAfterTrimNormalized);
    lim = axis;
    plot([algorithmData.IndexExplorationEnd algorithmData.IndexExplorationEnd],[lim(3) lim(4)],'k-.');
    axis([algorithmData.IndexExplorationStart algorithmData.IndexConsolidationEnd lim(3) lim(4)]);
    grid minor;
    xlabel('Iteration Step');
    ylabel('Normalized Box Measure (V/V_{ds})');
    legend({'Before Trimming Operation','After Trimming Operation'});

    % Interval Sizes of Individual Design Variables (normalized by Design Space Sizes)
    legTextBefore = {};
    legTextAfter = {};
    for i=1:size(algorithmData.DesignVariableIntervalSizeBeforeTrim,2)
        legTextBefore{end+1} = sprintf('Design Variable %d Before Trimming',i);
        legTextAfter{end+1} = sprintf('Design Variable %d After Trimming',i);
    end
    figureHandle(3) = figure;
    plot(algorithmData.DesignVariableIntervalSizeBeforeTrimNormalized);
    hold on;
    plot(algorithmData.DesignVariableIntervalSizeAfterTrimNormalized);
    lim = axis;
    plot([algorithmData.IndexExplorationEnd algorithmData.IndexExplorationEnd],[lim(3) lim(4)],'k-.');
    axis([algorithmData.IndexExplorationStart algorithmData.IndexConsolidationEnd lim(3) lim(4)]);
    grid minor;
    xlabel('Iteration Step');
    ylabel('Normalized Design Variable Interval Sizes (I_i/I_{i,ds})');
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
    plot(algorithmData.RatioMeasureChangeBeforeTrim);
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
    plot(algorithmData.SamplePurity,...
        algorithmData.MeasureBeforeTrimNormalized,...
        '.-');
    hold on;
    grid minor;
    xlabel('Sample Purity');
    ylabel('Normalized Box Volume (V/V_{ds})');
    
    % 3D plot with purity/size/cost
    figureHandle(8) = figure;
    plot3(algorithmData.TotalFunctionEvaluations,...
        algorithmData.SamplePurity,...
        algorithmData.MeasureBeforeTrimNormalized,...
        '-','Marker','.');
    hold on;
    grid minor;
    xlabel('Total Number of Function Evaluations');
    ylabel('Sample Purity');
    zlabel('Normalized Box Volume (V/V_{ds})');

    % save if required
    if(~isempty(options.SaveFolder))
        save_print_figure(figureHandle(1),[options.SaveFolder,'Box-Metrics-TotalMeasure'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(2),[options.SaveFolder,'Box-Metrics-TotalMeasureNormalized'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(3),[options.SaveFolder,'Box-Metrics-IntervalSizeNormalized'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(4),[options.SaveFolder,'Box-Metrics-GrowthRate'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(5),[options.SaveFolder,'Box-Metrics-NumberLabelSamplePoints'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(6),[options.SaveFolder,'Box-Metrics-TotalMeasureRelativeRatio'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(7),[options.SaveFolder,'Box-Metrics-NormalizedMeasureSample'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(8),[options.SaveFolder,'Box-Metrics-TotalEvaluationsPurityNormalizedMeasure'],options.SaveFigureOptions{:});

        if(options.CloseFigureAfterSaving)
            close(figureHandle);
        end
    end

    if(nargout<1)
        clear figureHandle
    end
end

