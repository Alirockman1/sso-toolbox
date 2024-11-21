%   Some references:
%       - https://distill.pub/2020/bayesian-optimization/ 
%       - https://machinelearningmastery.com/what-is-bayesian-optimization/
%       - https://krasserm.github.io/2018/03/21/bayesian-optimization/

%% Cleanup
fclose all;
close all;
clear all;
clc;
more off;
diary off;


%% debugging
rng(4);


%% Documentation / Archive
rngState = rng;
saveFolder = save_diary_files(mfilename);
goldenRatio = (1+sqrt(5))/2;
figureSize = [goldenRatio 1]*8.5;


%%
designSpaceLowerBound = [0];
designSpaceUpperBound = [1];
objectiveFunction = @(x) -x.^2 .* sin(5*pi.*x).^6.0;

options = {'MaxIter',10,'InitialExplorationFactor',0.05,'ExplorationFactorUpdateFunction',[]};

[designOptimal,objectiveOptimal,optimizationOutput] = optimization_bayesian(...
	objectiveFunction,[],designSpaceLowerBound,designSpaceUpperBound,[],...
	options{:});


%%
visualize_bayesian_optimization_1d(optimizationOutput,'SaveFolder',saveFolder,'SaveFigureOptions',{'Size',figureSize});


% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;


%
function visualize_bayesian_optimization_1d(optimizationOutput,varargin)
	parser = inputParser;
    parser.addParameter('SaveFolder',[]);
    parser.addParameter('CloseFigureAfterSaving',false);
    parser.addParameter('SaveFigureOptions',{});
    parser.addParameter('NumberEvaluationPoint',1000);
    parser.parse(varargin{:});
    options = parser.Results;

    designSample = linspace(optimizationOutput.ProblemData.DesignSpaceLowerBound,...
    	optimizationOutput.ProblemData.DesignSpaceUpperBound,options.NumberEvaluationPoint)';
    objectiveValue = optimizationOutput.ProblemData.ObjectiveFunction(designSample);

    [predictedObjectiveValue,predictedStandardDeviation] = ...
    	optimizationOutput.InitialData.ObjectiveRegressionModel.predict(designSample);
    figure;
    hold all;
    grid minor;
    plot(designSample,objectiveValue,'-','color',color_palette_tol('blue'),'linewidth',2.0);
    plot(designSample,predictedObjectiveValue,'-','color',color_palette_tol('green'),'linewidth',2.0);
    plot(designSample,predictedObjectiveValue+predictedStandardDeviation,':','color',color_palette_tol('red'),'linewidth',2.0,'HandleVisibility','off');
    plot(designSample,predictedObjectiveValue-predictedStandardDeviation,':','color',color_palette_tol('red'),'linewidth',2.0);
    plot(optimizationOutput.InitialData.EvaluatedSample,optimizationOutput.InitialData.EvaluatedObjective,'.','color',color_palette_tol('cyan'),'MarkerSize',20);
    legend({'Objective Function','Model Prediction','Standard Deviation','Evaluated Points'},'location','southwest');
    xlim([min(designSample) max(designSample)]);
    title('Initial Data');

    if(~isempty(options.SaveFolder))
        save_print_figure(gcf,[options.SaveFolder,'InitialData'],options.SaveFigureOptions{:});
        if(options.CloseFigureAfterSaving)
            close(figureHandle);
        end
    end

    evaluatedSample = optimizationOutput.InitialData.EvaluatedSample;
    evaluatedObjective = optimizationOutput.InitialData.EvaluatedObjective;

    nIter = length(optimizationOutput.IterationData);
    for i=1:nIter
        [predictedObjectiveValue,predictedStandardDeviation] = ...
            optimizationOutput.IterationData(i).ObjectiveRegressionModelNew.predict(designSample);

        if(i==1)
            expectedImprovement = bayesian_acquisition_gaussian_expected_improvement(...
                designSample,optimizationOutput.InitialData.ObjectiveRegressionModel,...
                optimizationOutput.InitialData.OptimalObjective,optimizationOutput.IterationData(i).ExplorationFactor);
        else
            expectedImprovement = bayesian_acquisition_gaussian_expected_improvement(...
                designSample,optimizationOutput.IterationData(i-1).ObjectiveRegressionModelNew,...
                optimizationOutput.IterationData(i-1).OptimalObjective,optimizationOutput.IterationData(i).ExplorationFactor);
        end

        figure;
        hold all;
        grid minor;
        plot(designSample,objectiveValue,'-','color',color_palette_tol('blue'),'linewidth',2.0);
        plot(designSample,predictedObjectiveValue,'-','color',color_palette_tol('green'),'linewidth',2.0);
        plot(designSample,predictedObjectiveValue+predictedStandardDeviation,':','color',color_palette_tol('red'),'linewidth',2.0,'HandleVisibility','off');
        plot(designSample,predictedObjectiveValue-predictedStandardDeviation,':','color',color_palette_tol('red'),'linewidth',2.0);
        plot(evaluatedSample,evaluatedObjective,'k.','MarkerSize',20);
        plot(optimizationOutput.IterationData(i).EvaluatedSample,optimizationOutput.IterationData(i).EvaluatedObjective,'.','color',color_palette_tol('cyan'),'MarkerSize',20);
        plot(designSample,expectedImprovement,'-','color',color_palette_tol('purple'),'linewidth',2.0);
        legend({'Objective Function','Model Prediction','Standard Deviation',...
            'Previous Points','Newly Evaluated Points','Previous Expected Improvement'},'location','southwest');
        xlim([min(designSample) max(designSample)]);
        title(sprintf('Iteration %d',i));

        if(~isempty(options.SaveFolder))
            filename = sprintf('Iteration-%0*d',get_number_digits_integer(nIter),i);
            save_print_figure(gcf,[options.SaveFolder,filename],options.SaveFigureOptions{:});
            if(options.CloseFigureAfterSaving)
                close(figureHandle);
            end
        end

        evaluatedSample = [evaluatedSample;optimizationOutput.IterationData(i).EvaluatedSample];
        evaluatedObjective = [evaluatedObjective;optimizationOutput.IterationData(i).EvaluatedObjective];
    end
end

function visualize_bayesian_optimization_2d(optimizationOutput)
end