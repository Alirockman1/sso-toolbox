function figureHandle = plot_active_learning_model_training_metrics(algoData,varargin)
%PLOT_ACTIVE_LEARNING_MODEL_TRAINING_METRICS Visualization of performance
%   PLOT_ACTIVE_LEARNING_MODEL_TRAINING_METRICS is used to automatically create 
%   the most common necessary plots for analyzing the performance of the active 
%   learning algorithm used to find fast models of the complete spaces. These 
%   plots are:
%       - Iterations x Number of Designs used in Training
%       - Iterations x Ratio (Good/Bad) Designs used in Training
%       - Iterations x Number of Classification Errors Commited
%       - Iterations x Evaluation Metrics (Exploratory)
%       - Iterations x Evaluation Metrics (Boundary)
%   
%   PLOT_ACTIVE_LEARNING_MODEL_TRAINING_METRICS(ALGORITHMDATA) takes the  
%   algorithm performance data in ALGORITHMDATA and makes figures with the plots 
%   specified above.
%
%   PLOT_ACTIVE_LEARNING_MODEL_TRAINING_METRICS(...NAME,VALUE,...) also allows  
%   the specification of additional options. These are:
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
%       - ALGODATA : struct
%       - 'SaveFolder' : char OR string
%       - 'CloseFigureAfterSaving' : logical
%       - 'SaveFigureOptions' : cell
%
%   Outputs:
%       - FIGUREHANDLE : (1,5) Figure 
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
    figureHandle(1) = figure;
    plot(algoData.NumberPositiveSamplesInTrainingDataset);
    hold on;
    plot(algoData.NumberNegativeSamplesInTrainingDataset);
    lim = axis;
    axis([algoData.IndexStart algoData.IndexEnd lim(3) lim(4)]);
    grid minor;
    xlabel('Iteration Step');
    ylabel('Number of Designs Present in Training Dataset');
    legend({'Good Designs','Bad Designs'});
    
    figureHandle(2) = figure;
    plot(algoData.TrainingDataPositiveRatio);
    hold on;
    plot([algoData.IndexStart algoData.IndexEnd], [0.5 0.5], 'k--');
    lim = axis;
    axis([algoData.IndexStart algoData.IndexEnd lim(3) lim(4)]);
    grid minor;
    xlabel('Iteration Step');
    ylabel('Ratio (Positive/Total) of Designs used in Training');
    legend({'Realized','Target'});
    
    figureHandle(3) = figure;
    plot(algoData.NumberFalsePositiveExploration);
    hold on;
    plot(algoData.NumberFalseNegativeExploration);
    plot(algoData.NumberFalsePositiveBoundary);
    plot(algoData.NumberFalseNegativeBoundary);
    lim = axis;
    axis([algoData.IndexStart algoData.IndexEnd lim(3) lim(4)]);
    grid minor;
    xlabel('Iteration Step');
    ylabel('Number of Classification Errors Commited');
    legend({'False Positive (Exploratory)','False Negative (Exploratory)',...
        'False Positive (Boundary)','False Negative (Boundary)'});
    
    figureHandle(4) = figure;
    plot(algoData.AccuracyExploration);
    hold on;
    plot(algoData.PrecisionExploration);
    plot(algoData.RecallExploration);
    plot(algoData.SpecificityExploration);
    plot(algoData.F1scoreExploration);
    lim = axis;
    axis([algoData.IndexStart algoData.IndexEnd lim(3) lim(4)]);
    grid minor;
    xlabel('Iteration Step');
    ylabel('Evaluation Metrics (Exploratory)');
    legend({'Accuracy','Precision','Recall','Specificity','F_1-score'});
    
    figureHandle(5) = figure;
    plot(algoData.AccuracyBoundary);
    hold on;
    plot(algoData.PrecisionBoundary);
    plot(algoData.RecallBoundary);
    plot(algoData.SpecificityBoundary);
    plot(algoData.F1scoreBoundary);
    lim = axis;
    axis([algoData.IndexStart algoData.IndexEnd lim(3) lim(4)]);
    grid minor;
    xlabel('Iteration Step');
    ylabel('Evaluation Metrics (Boundary)');
    legend({'Accuracy','Precision','Recall','Specificity','F_1-score'});

    % save if required
    if(~isempty(options.SaveFolder))
        save_print_figure(figureHandle(1),[options.SaveFolder,'ActiveLearning-Metrics-NumberTraining'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(2),[options.SaveFolder,'ActiveLearning-Metrics-RatioTraining'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(3),[options.SaveFolder,'ActiveLearning-Metrics-TypeIandIIErrors'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(4),[options.SaveFolder,'ActiveLearning-Metrics-EvaluationExploratory'],options.SaveFigureOptions{:});
        save_print_figure(figureHandle(5),[options.SaveFolder,'ActiveLearning-Metrics-EvaluationBoundary'],options.SaveFigureOptions{:});

        if(options.CloseFigureAfterSaving)
            close(figureHandle);
        end
    end

    if(nargout<1)
        clear figureHandle
    end
end