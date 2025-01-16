function [designOptimal,objectiveOptimal,optimizationOutput] = optimization_bayesian(objectiveFunction,initialDesign,designSpaceLowerBound,designSpaceUpperBound,constraintFunction,varargin)
%OPTIMIZATION_BAYESIAN using Gaussian Process Regression (Minimization Variant)
%	OPTIMIZATION_BAYESIAN minimizes the given objective function using
%	bayesian optimization with a gaussian process regression surrogate model. 
%	NOTE: contrary to most literature, this assumes a bayesian optimization
%	being used to minimize the objective function instead of maximize.
%
%	DESIGNOPTIMAL = OPTIMIZATION_BAYESIAN(OBJECTIVEFUNCTION,
%	INITIALDESIGN,DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND,
%	CONSTRAINTFUNCTION) 
%
%	DESIGNOPTIMAL = OPTIMIZATION_BAYESIAN(...NAME,VALUE,...) allows for
%	the specification of additional options for this algorithm. Available 
%	options are:
%		- 'ObjectiveFunctionThresholdValue' : if a design is evaluated and has
%		an objective value which is lower than this threshold and constraints 
%		are satisfied, the optimization procedure stops and that design is 
%		returned. Default: -inf.
%		- 'MaxIter' : maximum number of iterations allowed. Default value: 30.
%		- 'NumberInitialEvaluation' : number of points to be evaluated at
%		the start to create the inital surrogate model. Default: 10.
%		- 'BatchTrialSize' : number of points to be randomly created and then
%		used as initial points to find the optimum of the acquisition function.
%		Default: 100.
%		- 'EvaluateFromBatch' : number of points from the batch to evaluate,
%		going from those that have the largest value of improvement to the 
%		smallest. Default: 1.
%		- 'SamplingFunction' : function used to create the batch sample points.
%		Default: @sampling_latin_hypercube.
%		- 'SamplingOptions' : options for the sampling function. Default is 
%		empty.
%		- 'InternalOptimizationFunction' : function used to try to maximze the
%		acquisition function. Default: optimization_fmincon_wrapper.
%		- 'InternalOptimizationOptions' : options for the function used to 
%		maximize the acquisition function. Default: {'Display','none'}.
%		- 'GaussianRegressionModelTrainingOptions' : options for 'fitrgp'.
%		By default, the hyperparameters are optimized without showing any
%		figures or printing to console.
%		- 'InitialExplorationFactor' : exploration factor to use at the first
%		iteration. Default: 0.1.
%		- 'ExplorationFactorUpdateFunction' : function to update the exploration
%		factor at each subsequent iteration. Must have the structure: 
%		'explorationFactorNew = f(explorationFactor,iteration,options)'. By
%		default, the exploration factor is always cut in half.
%		- 'ExplorationFactorUpdateOptions' : options for the update function
%		of the exploration factor. Default is empty.
%		- 'AcquisitionFunction' : function handle to be used for acquisition.
%		This function will be maximized to find the most promising candidates.
%		Default: bayesian_acquisition_gaussian_expected_improvement.
%		- 'UseSurrogateConstraint' : flag on how to deal with constraints. If
%		set to false, constraints are directly computed for each design being
%		considered. If true, a surrogate model is also created for the 
%		constraints (together with the objective), and that model is used 
%		instead when finding the best candidates for evaluation with the
%		acquisition function. Default: false.
%
%	[DESIGNOPTIMAL,OBJECTIVEOPTIMAL] = OPTIMIZATION_BAYESIAN(...) also returns
%	the value of the objective function OBJECTIVEOPTIMAL for the optimal design.
%
%	[DESIGNOPTIMAL,OBJECTIVEOPTIMAL,OPTIMIZATIONOUTPUT] = 
%	OPTIMIZATION_BAYESIAN(...) returns additional information regarding
%	the optimization procedure.
%
%   Input:
%		- OBJECTIVEFUNCTION : function_handle
%		- INITIALDESIGN : (1,nDesignVariable) double
%		- DESIGNSPACELOWERBOUND : (1,nDesignVariable) double
%		- DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double
%		- CONSTRAINTFUNCTION : function_handle
%		- Name-value pair arguments : 
%			-- 'ObjectiveFunctionThresholdValue' : double
%			-- 'MaxIter' : integer
%			-- 'NumberInitialEvaluation' : integer
%			-- 'BatchTrialSize' : integer
%			-- 'EvaluateFromBatch' : integer
%			-- 'SamplingFunction' : function_handle
%			-- 'SamplingOptions' : cell
%			-- 'InternalOptimizationFunction' : function_handle
%			-- 'InternalOptimizationOptions' : cell
%			-- 'GaussianRegressionModelTrainingOptions' : cell
%			-- 'InitialExplorationFactor' : double
%			-- 'ExplorationFactorUpdateFunction' : function_handle
%			-- 'ExplorationFactorUpdateOptions' : cell
%			-- 'AcquisitionFunction' : function_handle
%			-- 'UseSurrogateConstraint' : logical
%
%   Output:
%		- DESIGNOPTIMAL : (1,nDesignVariable) double
%		- OBJECTIVEOPTIMAL : double
%		- OPTIMIZATIONOUTPUT : structure
%			-- ProblemData : structure
%				--- DesignSpaceLowerBound : (1,nDimension) double
%				--- DesignSpaceUpperBound : (1,nDimension) double
%				--- InitialDesign : (nInitial,nDimension) double
%				--- ObjectiveFunction : function_handle
%				--- ConstraintFunction : function_handle
%				--- Options : structure
%			-- InitialData : structure
%				--- EvaluatedSample : (nInitial,nDimension) double
%				--- EvaluatedObjective : (nInitial,1) double
%				--- EvaluatedInequalityConstraint : (nInitial,nInequality) 
%				double
%				--- EvaluatedEqualityConstraint : (nInitial,nEquality) double
%				--- ObjectiveRegressionModel : RegressionGP
%				--- InequalityConstraintRegressionModel : (1,nInequality) cell
%				--- EqualityConstraintRegressionModel : (1,nEquality) cell
%				--- OptimumCandidate : (1,nDimension) double
%				--- OptimalObjective : double
%			-- IterationData : (1,nIter) structure
%				--- ExplorationFactor : double
%				--- BatchSample : (nBatch,nDimension) double
%				--- CandidateSample : (nBatch,nDimension) double
%				--- PredictedCandidateImprovement : (nBatch,1) double
%				--- InternalOptimizationOutput : cell
%				--- EvaluatedSample : (nInitial,nDimension) double
%				--- EvaluatedObjective : (nInitial,1) double
%				--- EvaluatedInequalityConstraint : (nInitial,nInequality) 
%				double
%				--- EvaluatedEqualityConstraint : (nInitial,nEquality) double
%				--- ObjectiveRegressionModelNew : RegressionGP
%				--- InequalityConstraintRegressionModelNew : (1,nInequality) 
%				cell
%				--- EqualityConstraintRegressionModelNew : (1,nEquality) cell
%				--- OptimumCandidate : (1,nDimension) double
%				--- OptimalObjective : double
%
%   See also design_optimize_quantities_of_interest,
%	design_optimize_performance_score.
%	bayesian_acquisition_gaussian_expected_improvement, 
%	bayesian_acquisition_gaussian_probability_of_improvement.
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
	parser.addParameter('ObjectiveFunctionThresholdValue',-inf);
	parser.addParameter('MaxIter',30);
    parser.addParameter('NumberInitialEvaluation',10);
	parser.addParameter('BatchTrialSize',100);
	parser.addParameter('EvaluateFromBatch',1);
	parser.addParameter('SamplingFunction',@sampling_latin_hypercube);
	parser.addParameter('SamplingOptions',{});
	parser.addParameter('InternalOptimizationFunction',@optimization_fmincon_wrapper);
	parser.addParameter('InternalOptimizationOptions',{'Display','none'});
	parser.addParameter('GaussianRegressionModelTrainingOptions',...
        {'Standardize',true,'OptimizeHyperparameters',{'BasisFunction','KernelFunction','KernelScale','Sigma'},...
        'HyperparameterOptimizationOptions',struct('ShowPlots',false,'Verbose',0)});
	parser.addParameter('InitialExplorationFactor',0.2);
	parser.addParameter('ExplorationFactorUpdateFunction',@(explorationFactor,iteration)explorationFactor/2);
	parser.addParameter('ExplorationFactorUpdateOptions',{});
    parser.addParameter('AcquisitionFunction',@bayesian_acquisition_gaussian_expected_improvement);
    parser.addParameter('UseSurrogateConstraint',false);
	parser.parse(varargin{:});
    options = parser.Results;

    rngState = rng;

	% do initial sampling and evaluation
	nDimension = size(designSpaceLowerBound,2);
	designSpace = [designSpaceLowerBound;designSpaceUpperBound];
	nInitialBatch = options.NumberInitialEvaluation - size(initialDesign,1);
	initialBatch = [];
	if(nInitialBatch>0)
		initialBatch = options.SamplingFunction(designSpace,nInitialBatch,options.SamplingOptions{:});
	end
	designSampleTrain = [initialDesign;initialBatch];

	objectiveValueTrain = objectiveFunction(designSampleTrain);

	if(isempty(constraintFunction))
		inequalityConstraintValueTrain = [];
		equalityConstraintValueTrain = [];
	elseif(nargout(constraintFunction)==1)
		inequalityConstraintValueTrain = constraintFunction(designSampleTrain);
		equalityConstraintValueTrain = [];
	else
		[inequalityConstraintValueTrain,equalityConstraintValueTrain] = constraintFunction(designSampleTrain);
	end

	% train initial models
	gaussianRegressionModelObjective = fitrgp(designSampleTrain,objectiveValueTrain,options.GaussianRegressionModelTrainingOptions{:});
	
    gaussianRegressionModelInequalityConstraint = {};
	gaussianRegressionModelEqualityConstraint = {};
	if(options.UseSurrogateConstraint)
		nInequalityConstraint = size(inequalityConstraintValueTrain,2);
	    for i=1:nInequalityConstraint
			gaussianRegressionModelInequalityConstraint{i} = fitrgp(designSampleTrain,inequalityConstraintValueTrain(:,i),options.GaussianRegressionModelTrainingOptions{:});
		end

		nEqualityConstraint = size(equalityConstraintValueTrain,2);
	    for i=1:nEqualityConstraint
			gaussianRegressionModelEqualityConstraint{i} = fitrgp(designSampleTrain,equalityConstraintValueTrain(:,i),options.GaussianRegressionModelTrainingOptions{:});
        end
    end
    gaussianRegressionModelInequalityConstraintNew = {};
	gaussianRegressionModelEqualityConstraintNew = {};

	% check which design is currently optimal
	satisfiesConstraint = all(inequalityConstraintValueTrain<=0,2);
	if(~any(satisfiesConstraint) || isempty(inequalityConstraintValueTrain))
		% if nothing satisfied constraint, keep all for choosing best
		isValidSample = true(size(designSampleTrain,1),1);
	else
		% keep only those that satisfy constraints
		isValidSample = satisfiesConstraint;
	end
	[objectiveOptimal,iOptimalValid] = min(objectiveValueTrain(isValidSample));
	iOptimal = convert_index_base(isValidSample,iOptimalValid,'backward');
	designOptimal = designSampleTrain(iOptimal,:);

	isOutputOptimizationInformation = (nargout>2);
	if(isOutputOptimizationInformation)
		optimizationOutput.ProblemData = struct(...
			'RNGState',rngState,...
			'DesignSpaceLowerBound',designSpaceLowerBound,...
			'DesignSpaceUpperBound',designSpaceUpperBound,...
			'InitialDesign',initialDesign,...
			'ObjectiveFunction',objectiveFunction,...
			'ConstraintFunction',constraintFunction,...
			'Options',options);

		optimizationOutput.InitialData = struct(...
			'EvaluatedSample',designSampleTrain,...
			'EvaluatedObjective',objectiveValueTrain,...
			'EvaluatedInequalityConstraint',inequalityConstraintValueTrain,...
			'EvaluatedEqualityConstraint',equalityConstraintValueTrain,...
			'ObjectiveRegressionModel',gaussianRegressionModelObjective,...
			'InequalityConstraintRegressionModel',{gaussianRegressionModelInequalityConstraint},...
			'EqualityConstraintRegressionModel',{gaussianRegressionModelEqualityConstraint},...
			'OptimumCandidate',designOptimal,...
			'OptimalObjective',objectiveOptimal);
	end

	% loop
	iteration = 1;
	foundOptimal = false;
	explorationFactor = options.InitialExplorationFactor;
	while(~foundOptimal && iteration<=options.MaxIter)
		% update exploration factor 
		if(iteration>1 && ~isempty(options.ExplorationFactorUpdateFunction))
			explorationFactor = options.ExplorationFactorUpdateFunction(explorationFactor,iteration,options.ExplorationFactorUpdateOptions{:});
		end

		% find next points to evaluate
		designSampleNew = options.SamplingFunction(designSpace,options.BatchTrialSize,options.SamplingOptions{:});
		internalObjectiveFunction = @(x)-options.AcquisitionFunction(x,gaussianRegressionModelObjective,objectiveOptimal,explorationFactor);

		if(options.UseSurrogateConstraint)
			internalConstraintFunction = @(x)gaussian_regression_constraint(gaussianRegressionModelInequalityConstraint,gaussianRegressionModelEqualityConstraint,x);
		else
			internalConstraintFunction = constraintFunction;
		end

		predictedDesignOptimal = nan(options.BatchTrialSize,nDimension);
		predictedCandidateImprovement = nan(options.BatchTrialSize,1);
		predictedOptimizationOutput = cell(options.BatchTrialSize,1);
		for i=1:options.BatchTrialSize
			[predictedDesignOptimal(i,:),predictedCandidateImprovement(i),predictedOptimizationOutput{i}] = ...
				options.InternalOptimizationFunction(internalObjectiveFunction,designSampleNew(i,:),...
					designSpaceLowerBound,designSpaceUpperBound,internalConstraintFunction,...
					options.InternalOptimizationOptions{:});
		end

		% evaluate chosen samples
		[~,iSortedBestPrediction] = sort(predictedCandidateImprovement,'ascend'); % smallest negative -> largest improvement
		designEvaluate = predictedDesignOptimal(iSortedBestPrediction(1:options.EvaluateFromBatch),:);
		designEvaluate = unique(designEvaluate,'rows');
		
		objectiveEvaluate = objectiveFunction(designEvaluate);
		
		if(isempty(constraintFunction))
			inequalityConstraintEvaluate = [];
			equalityConstraintEvaluate = [];
		elseif(nargout(constraintFunction)==1)
			inequalityConstraintEvaluate = constraintFunction(designEvaluate);
			equalityConstraintEvaluate = [];
		else
			[inequalityConstraintEvaluate,equalityConstraintEvaluate] = constraintFunction(designEvaluate);
		end

		% re-train model with new data
		designSampleTrainNew = [designSampleTrain;designEvaluate];
		objectiveValueTrainNew = [objectiveValueTrain;objectiveEvaluate];
		gaussianRegressionModelObjectiveNew = fitrgp(designSampleTrainNew,objectiveValueTrainNew,...
			options.GaussianRegressionModelTrainingOptions{:});

        inequalityConstraintValueTrainNew = [inequalityConstraintValueTrain;inequalityConstraintEvaluate];
        equalityConstraintValueTrainNew = [equalityConstraintValueTrain;equalityConstraintEvaluate];
		if(options.UseSurrogateConstraint)
			for i=1:nInequalityConstraint
				gaussianRegressionModelInequalityConstraintNew{i} = fitrgp(designSampleTrainNew,inequalityConstraintValueTrainNew(:,i),...
					options.GaussianRegressionModelTrainingOptions{:});
            end
			for i=1:nEqualityConstraint
				gaussianRegressionModelEqualityConstraintNew{i} = fitrgp(designSampleTrainNew,equalityConstraintValueTrainNew(:,i),...
					options.GaussianRegressionModelTrainingOptions{:});
			end
		end

		% check for convergence 

		% log
		if(isOutputOptimizationInformation)
			optimizationOutput.IterationData(iteration) = struct(...
				'ExplorationFactor',explorationFactor,...
				'BatchSample',designSampleNew,...
				'CandidateSample',predictedDesignOptimal,...
				'PredictedCandidateImprovement',predictedCandidateImprovement,...
				'InternalOptimizationOutput',{predictedOptimizationOutput},...
				'EvaluatedSample',designEvaluate,...
				'EvaluatedObjective',objectiveEvaluate,...
				'EvaluatedInequalityConstraint',inequalityConstraintEvaluate,...
				'EvaluatedEqualityConstraint',equalityConstraintEvaluate,...
				'ObjectiveRegressionModelNew',gaussianRegressionModelObjectiveNew,...
				'InequalityConstraintRegressionModelNew',{gaussianRegressionModelInequalityConstraintNew},...
				'EqualityConstraintRegressionModelNew',{gaussianRegressionModelEqualityConstraintNew},...
				'OptimumCandidate',designOptimal,...
				'OptimalObjective',objectiveOptimal);
		end

		% update values
		designSampleTrain = designSampleTrainNew;
		objectiveValueTrain = objectiveValueTrainNew;
		gaussianRegressionModelObjective = gaussianRegressionModelObjectiveNew;
		inequalityConstraintValueTrain = inequalityConstraintValueTrainNew;
		equalityConstraintValueTrain = equalityConstraintValueTrainNew;
		gaussianRegressionModelInequalityConstraint = gaussianRegressionModelInequalityConstraintNew;
		gaussianRegressionModelEqualityConstraint = gaussianRegressionModelEqualityConstraintNew;
		iteration = iteration + 1;

		% check which design is currently optimal
		satisfiesConstraint = all(inequalityConstraintValueTrain<=0,2);
		if(~any(satisfiesConstraint) || isempty(inequalityConstraintValueTrain))
			% if nothing satisfied constraint, keep all for choosing best
			isValidSample = true(size(objectiveValueTrain,1),1);
		else
			% keep only those that satisfy constraints
			isValidSample = satisfiesConstraint;
		end
		[objectiveOptimal,iOptimalValid] = min(objectiveValueTrain(isValidSample));
		iOptimal = convert_index_base(isValidSample,iOptimalValid,'backward');
		designOptimal = designSampleTrain(iOptimal,:);
    end
end

function [predictedInequalityConstraint,predictedEqualityConstraint] = gaussian_regression_constraint(gaussianRegressionModelInequalityConstraint,gaussianRegressionModelEqualityConstraint,designSample)
	nInequalityConstraint = length(gaussianRegressionModelInequalityConstraint);
	nEqualityConstraint = length(gaussianRegressionModelEqualityConstraint);
	nSample = size(designSample,1);

	predictedInequalityConstraint = nan(nSample,nInequalityConstraint);
	for i=1:nInequalityConstraint
		predictedInequalityConstraint(:,i) = gaussianRegressionModelInequalityConstraint{i}.predict(designSample);
    end

    predictedEqualityConstraint = nan(nSample,nEqualityConstraint);
	for i=1:nEqualityConstraint
		predictedEqualityConstraint(:,i) = gaussianRegressionModelEqualityConstraint{i}.predict(designSample);
    end
end