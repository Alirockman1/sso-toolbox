function [fastForwardModel,problemData,iterData] = active_learning_model_training(bottomUpMapping,designSpaceLowerBound,designSpaceUpperBound,varargin)
%ACTIVE_LEARNING_MODEL_TRAINING Design Fast-Forward model with active learning 
%   ACTIVE_LEARNING_MODEL_TRAINING uses an active learning strategy to train a 
%   given design fast-forward model. This is accomplished with an uncertainty 
%   sampling query strategy mixed with stochatsic sampling for general 
%   exploration. In the context of top-down design, these regions of uncertainty
%   are the boundaries between good and bad / physically feasible and physically
%   infeasible designs.
%
%   FASTFORWARDMODEL = ACTIVE_LEARNING_MODEL_TRAINING(BOTTOMUPMAPPING,
%   DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND) creates a design fast-forward 
%   model FASTFORWARDMODEL of the bottom-up mapping function BOTTOMUPMAPPING  
%   in the region defined by the lower and upper boundaries 
%   DESIGNSPACELOWERBOUND and DESIGNSPACEUPPERBOUND. It is assumed there is no 
%   lower critical limit value for the measures, while the upper critical limit 
%   value is assumed to be 0. "Positive" designs are those whose measures 
%   satisfy said critical limits, and "negative" designs are those that don't. 
%   For training, default options are used.
%
%   FASTFORWARDMODEL = ACTIVE_LEARNING_MODEL_TRAINING(BOTTOMUPMAPPING,
%   DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND,MEASURELOWERLIMIT) allows 
%   one to specify the lower limit of the measures MEASURELOWERLIMIT. This can
%   be given as a single value thta applies to all measures or individual values
%   for each measure. When not given or empty('[]'), it's assumed there's no 
%   lower limit. One can use 'nan' entries for '-inf'.
%
%   FASTFORWARDMODEL = ACTIVE_LEARNING_MODEL_TRAINING(BOTTOMUPMAPPING,
%   DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND,MEASURELOWERLIMIT,
%   MEASUREUPPERLIMIT) allows one to specify the upper limit of the measures 
%   MEASURELOWERLIMIT. This can be given as a single value that applies to all 
%   measures or individual values for each measure. When not given or 
%   empty('[]'), 0 is used. One can use 'nan' entries for '+inf'.
%
%   FASTFORWARDMODEL = ACTIVE_LEARNING_MODEL_TRAINING(BOTTOMUPMAPPING,
%   DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND,MEASURELOWERLIMIT,
%   MEASUREUPPERLIMIT,OPTIONS) allows one to set additional options for the 
%   training. See more on 'active_learning_model_training_options'. Don't give 
%   this argument or set it to '[]' to use default options.
%
%   [FASTFORWARDMODEL,PROBLEMDATA] = ACTIVE_LEARNING_MODEL_TRAINING(...) returns
%   a structure with the processed parameters used during training. This can
%   be later used for reproducibility, to check for any issues in the input,
%   and/or for plotting the performance of the algorithm. Its fields are:
%       - BottomUpMapping : bottom-up mapping as specified in the input
%       - DesignSpaceLowerBound : lower boundary of the design space
%       - DesignSpaceUpperBound : upper boundary of the design space
%       - AlgorithmOptions : options as used by the algorithm
%       - MeasureLowerCriticalLimit : lower critical limit for measures
%       - MeasureUpperCriticalLimit : upper critical limit for measures
%       - UsePerformance : flag to indicate whether performance measures are 
%       being used for training, or the physical feasibility ones
%       - InitialRNGState : 'rng' information as taken at the start of the call
%
%   [FASTFORWARDMODEL,PROBLEMDATA,ITERDATA] = 
%   ACTIVE_LEARNING_MODEL_TRAINING(...) also returns the data gathered during
%   each iteration, which can be used for checking the performance/behavior of  
%   the algorithm. Its fields are:
%       - NumberGeneratedCandidates : number of candidates generated in sampling
%       - FastForwardModel : fast-forward model before training
%       - EvaluatedSamples : new samples evaluated for training
%       - PerformanceMeasure : performance measures of said samples
%       - PhysicalFeasibilityMeasure : physical feasibility measures of samples
%       - SelectedBoundary : boolean label indicating for each sample if it 
%       was selected for being supposedly close to the decision boundary (true)
%       or not (false).
%       - ScorePredicted : score as predicted by the fast-forward model for each
%       new sample
%       - LabelPredicted : label as predicted by the fast-forward model for each
%       new sample
%       - ScoreEvaluated : score as evaluated with the bottom-up mapping results
%       - LabelEvaluated : label as evaluated with the bottom-up mapping results
%
%   Inputs:
%       - BOTTOMUPMAPPING : 'run_bottom_up_mapping'-compatible input
%       - PARAMETER : user-defined
%       - DESIGNSPACELOWERBOUND : (1,nDesignVariable) double
%       - DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double
%       - MEASURELOWERLIMIT : (1,nMeasure) double, optional
%       - MEASUREUPPERLIMIT : (1,nMeasure) double, optional
%       - OPTIONS : structure OR name-value pair arguments, optional
%
% Outputs:
%       - FASTFORWARDMODEL : DesignFastForward
%       - PROBLEMDATA : structure
%           -- BottomUpMapping : function handle
%           -- ConstantParameters : user-defined
%           -- DesignSpaceLowerBound : (1,nDesignVariable) double
%           -- DesignSpaceUpperBound : (1,nDesignVariable) double
%           -- AlgorithmOptions : struct
%           -- MeasureLowerCriticalLimit : (1,nMeasure) double
%           -- MeasureUpperCriticalLimit : (1,nMeasure) double
%           -- UsePerformance : logical
%           -- InitialRNGState : struct
%       - ITERDATA : (1,nIter) structure
%           -- NumberGeneratedCandidates : double
%           -- FastForwardModel : DesignFastForward
%           -- EvaluatedSamples : (nSample,nDesignVariable) double
%           -- PerformanceMeasure : (nSample,nPerformance) double
%           -- PhysicalFeasibilityMeasure : (nSample,nComponent) double
%           -- SelectedBoundary : (nSample,1) logical
%           -- ScorePredicted : (nSample,1) double
%           -- LabelPredicted : (nSample,1) logical
%           -- ScoreEvaluated : (nSample,1) double
%           -- LabelEvaluated : (nSample,1) logical
%
%   See also run_bottom_up_mapping, active_learning_model_training_options.
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

    %% Input Processing/Validation
    % Constant system parameters
    % Quantities of Interest - Critical Lower Limit
    if(size(varargin,2)<1 || isempty(varargin{1}))
        measureLowerLimit = -inf;
    else
        measureLowerLimit = varargin{1};
    end
    measureLowerLimit(isnan(measureLowerLimit)) = -inf;
    
    % Quantities of Interest - Critical Upper Limit
    if(size(varargin,2)<2 || isempty(varargin{2}))
        measureUpperLimit = 0;
    else
        measureUpperLimit = varargin{2};
    end
    measureUpperLimit(isnan(measureUpperLimit)) = +inf;
    
    % Options
    defaultOptions = active_learning_model_training_options;
    if(size(varargin,2)<3 || isempty(varargin{3}))
        options = defaultOptions;
    else
        inputOptions = parser_variable_input_to_structure(varargin{3:end});
        options = merge_name_value_pair_argument(defaultOptions,inputOptions);
    end

    % set up console logging
    console = ConsoleLogging(options.LoggingLevel);

    % decide if performance or physical feasibility will be used for training
    usePhysicalFeasibility = false(2,1);
    usePhysicalFeasibility(1) = (isnumeric(options.MeasureChoice) && options.MeasureChoice==2);
    usePhysicalFeasibility(2) = ((ischar(options.MeasureChoice) || isstring(options.MeasureChoice)) ...
        && strcmpi(erase(options.MeasureChoice,{'-','_'}),'physicalfeasibility'));
    usePerformance = ~any(usePhysicalFeasibility);

    designSpace = [designSpaceLowerBound;designSpaceUpperBound];
    
    % Initial Model
    if(isa(options.FastForwardModelInitial,'function_handle'))
        % create initial object
        fastForwardModel = options.FastForwardModelInitial(measureLowerLimit,measureUpperLimit,options.FastForwardModelOptions{:});

        % perform initial sampling to create new object
        initialRegion = ternary(isempty(options.NewModelInitialRegion),designSpace,options.NewModelInitialRegion);
        grownRegion = design_box_grow_fixed(initialRegion,designSpaceLowerBound,designSpaceUpperBound,options.NewModelInitialRegionGrowth);
        initialSample = sampling_random(grownRegion,options.NumberSamplesEvaluation(1));
        initialSample = unique(initialSample,'rows');
        [performanceMeasure,physicalFeasibilityMeasure] = bottomUpMapping.response(initialSample);
        measure = ternary(usePerformance,performanceMeasure,physicalFeasibilityMeasure);
        fastForwardModel = fastForwardModel.train_model(initialSample,measure);
    elseif(isa(options.FastForwardModel,'DesignFastForward'))
        fastForwardModel = options.FastForwardModel;
    else
        console.error('ActiveLearning:FFModel:Unrecognized','Fast forward model not recognized.');
    end

    if(nargout>=2)
        problemData = struct(...
            'BottomUpMapping',bottomUpMapping,...
            'DesignSpaceLowerBound',designSpaceLowerBound,...
            'DesignSpaceUpperBound',designSpaceUpperBound,...
            'MeasureLowerCriticalLimit',measureLowerLimit,...
            'MeasureUpperCriticalLimit',measureUpperLimit,...
            'AlgorithmOptions',options,...
            'UsePerformance',usePerformance,...
            'InitialRNGState',rng);
    end
    
    % Initialize Iteration Data Log
    isOutputIterData = (nargout>=3);
    if(isOutputIterData)
        iterData = struct(...
            'NumberGeneratedCandidates',[],...
            'FastForwardModel',[],...
            'EvaluatedSamples',[],...
            'PerformanceMeasure',[],...
            'PhysicalFeasibilityMeasure',[],...
            'BottomUpMappingOutput',[],...
            'SelectedBoundary',[],...
            'ScorePredicted',[],...
            'LabelPredicted',[],...
            'ScoreEvaluated',[],...
            'LabelEvaluated',[]);
    end
    
    
    %% Active Learning Strategy
    nCandidate = options.NumberCandidatesInitial;
    iteration = 1;
    hasConverged = false;
    cycleCounter = 0;
    while(~hasConverged && (iteration<=options.MaxIter))
        nPositiveTrain = sum(fastForwardModel.LabelTrain);
        nNegativeTrain = sum(~fastForwardModel.LabelTrain);

        console.info([repelem('=',120),'\n']);
        console.info('Current Ratio of Training Data: %g (%d Positive / %d Negative; %d Total)\n',...
            nPositiveTrain/nNegativeTrain,nPositiveTrain,nNegativeTrain,nNegativeTrain+nPositiveTrain);

        console.info('Determining number of samples to be evaluated and candidates to be generated... ');
        tic
        nSample = get_current_array_entry(options.NumberSamplesEvaluation,iteration);
        nCandidate = options.FunctionUpdateNumberCandidates(options.NumberCandidatesInitial,nCandidate,iteration);
        console.info('Elapsed time is %g seconds.\n',toc);

        console.info('Deciding on ratio between expected positive and negative designs...');
        tic
        nSampleExploratory = round(nSample*options.ExplorationToBoundaryRatio);
        if(nSampleExploratory<=nNegativeTrain-nPositiveTrain) % too many negative designs already in training set
            nNegativeEvaluate = 0;
            nPositiveEvaluate = nSampleExploratory;
        elseif(nSampleExploratory<=nPositiveTrain-nNegativeTrain) % too many positive designs already in training set
            nPositiveEvaluate = 0;
            nNegativeEvaluate = nSampleExploratory;
        else % close balance, try to match numbers
            nPositiveEvaluate = round((nSampleExploratory - (nPositiveTrain-nNegativeTrain))/2);
            nNegativeEvaluate = nSampleExploratory - nPositiveEvaluate;
        end
        nBoundaryEvaluate = max(nSample - nSampleExploratory,0);
        console.info('Elapsed time is %g seconds.\n',toc);

        % Sampling region
        console.info('Generating initial candidate starting points and directions... ');
        tic
        nDesiredSample = struct(...
            'Positive',nPositiveEvaluate,...
            'Negative',nNegativeEvaluate,...
            'Boundary',nBoundaryEvaluate);
        [designSample,scorePrediction,labelBoundary] = options.ActiveLearningSamplingFunction(fastForwardModel,designSpace,nDesiredSample,nCandidate,options.ActiveLearningSamplingOptions{:});
        labelPrediction = (scorePrediction<=0);
        console.info('Elapsed time is %g seconds.\n',toc);

        % evaluation
        console.info('Evaluating new points... ');
        tic
        [performanceMeasure,physicalFeasibilityMeasure,bottomUpMappingOutput] = bottomUpMapping.response(designSample);
        measure = ternary(usePerformance,performanceMeasure,physicalFeasibilityMeasure);
        console.info('Elapsed time is %g seconds.\n',toc);
        [labelEvaluated,scoreEvaluated] = fastForwardModel.classify_measure(measure);
        console.debug('- Number of true positive boundary designs: %d\n',sum(labelBoundary & labelPrediction & labelEvaluated));
        console.debug('- Number of true negative boundary designs: %d\n',sum(labelBoundary & ~labelPrediction & ~labelEvaluated));
        console.debug('- Number of true positive exploratory designs: %d\n',sum(~labelBoundary & labelPrediction & labelEvaluated));
        console.debug('- Number of true negative exploratory designs: %d\n',sum(~labelBoundary & ~labelPrediction & ~labelEvaluated));

        % Training
        console.info('Training model with new information... ');
        tic
        fastForwardModelOld = fastForwardModel;
        fastForwardModel = fastForwardModel.update_model(designSample,measure);
        console.info('Elapsed time is %g seconds.\n',toc);
        
        % Check overall performance
        predictionsExploration = labelPrediction(~labelBoundary);
        evaluatedExploration = labelEvaluated(~labelBoundary);

        % exploratory confusion matrix
        [falsePositiveExploratory,falseNegativeExploratory,~,~] = ...
            classification_confusion_matrix(evaluatedExploration,predictionsExploration);
        
        % Check convergence
        % Note: mu here is not a good representation of the actual volume,
        % but as it converges, it means the region is not growing anymore
        console.info('Checking convergence... ');
        tic
        type1ErrorRatio = sum(falsePositiveExploratory)/size(designSample,1);
        type2ErrorRatio = sum(falseNegativeExploratory)/size(designSample,1);
        convergenceCriteria = false(3,1);
        convergenceCriteria(1) = (type1ErrorRatio<=options.ExplorationErrorTolerance);
        convergenceCriteria(2) = (type2ErrorRatio<=options.ExplorationErrorTolerance);
        convergenceCriteria(3) = ~options.FixIter;
        if(all(convergenceCriteria))
            cycleCounter = cycleCounter + 1;
        else
            cycleCounter = max([0,cycleCounter-1]);
        end
        if(cycleCounter>=options.CyclesToConverge)
            hasConverged = true;
        end
        console.info('Elapsed time is %g seconds.\n',toc);
        console.debug('- Number of False Positive Samples (Exploratory): %d\n',sum(falsePositiveExploratory));
        console.debug('- Number of False Negative Samples (Exploratory): %d\n',sum(falseNegativeExploratory));
        
        % Save Data
        if(isOutputIterData)
            console.info('Logging relevant information... ');
            tic
            iterData(iteration) = struct(...
                'NumberGeneratedCandidates',nCandidate,...
                'FastForwardModel',fastForwardModelOld,...
                'EvaluatedSamples',designSample,...
                'PerformanceMeasure',performanceMeasure,...
                'PhysicalFeasibilityMeasure',physicalFeasibilityMeasure,...
                'BottomUpMappingOutput',bottomUpMappingOutput,...
                'SelectedBoundary',labelBoundary,...
                'ScorePredicted',scorePrediction,...
                'LabelPredicted',labelPrediction,...
                'ScoreEvaluated',scoreEvaluated,...
                'LabelEvaluated',labelEvaluated);
            console.info('Elapsed time is %g seconds.\n',toc);
        end

        % Update Loop Information
        console.info('Done with iteration %d!\n\n',iteration);
        iteration = iteration + 1;
    end
    
    console.info('Done in %d iterations!\n',iteration-1);
end