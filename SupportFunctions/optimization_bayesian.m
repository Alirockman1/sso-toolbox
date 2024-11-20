function [designOptimal,objectiveOptimal,optimizationOutput] = optimization_bayesian(objectiveFunction,initialDesign,designSpaceLowerBound,designSpaceUpperBound,constraintFunction,varargin)
	parser = inputParser;
	parser.addParameter('ObjectiveFunctionThresholdValue',-inf);
	parser.addParameter('MaxIter',10);
	parser.addParameter('BatchTrialSize',10);
	parser.addParameter('EvaluateFromBatch',1);
	parser.addParameter('SamplingFunction',@sampling_latin_hypercube);
	parser.addParameter('SamplingOptions',{});
	parser.addParameter('InternalOptimizationFunction',@optimization_sampling);
	parser.addParameter('InternalOptimizationOptions',{});
	parser.addParameter('GaussianRegressionModelTrainingOptions',{});
	parser.addParameter('InitialExplorationFactor',0.2);
	parser.addParameter('ExplorationFactorUpdateFunction',@(explorationFactor,iteration)explorationFactor/2);
	parser.addParameter('ExplorationFactorUpdateOptions',{});
	parser.parse(varargin{:});
    options = parser.Results;

    rngState = rng;

	% do initial sampling and evaluation
	nDimension = size(designSpaceLowerBound,2);
	designSpace = [designSpaceLowerBound;designSpaceUpperBound];
	nInitialBatch = options.BatchTrialSize - size(initialDesign,1);
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
	
	nInequalityConstraint = size(inequalityConstraintValueTrain,2);
	gaussianRegressionModelInequalityConstraint = {};
    for i=1:nInequalityConstraint
		gaussianRegressionModelInequalityConstraint{i} = fitrgp(designSampleTrain,inequalityConstraintValueTrain(:,i),options.GaussianRegressionModelTrainingOptions{:});
	end

	nEqualityConstraint = size(equalityConstraintValueTrain,2);
	gaussianRegressionModelEqualityConstraint = {};
    for i=1:nEqualityConstraint
		gaussianRegressionModelEqualityConstraint{i} = fitrgp(designSampleTrain,equalityConstraintValueTrain(:,i),options.GaussianRegressionModelTrainingOptions{:});
	end

	% check which design is currently optimal
	satisfiesConstraint = all(inequalityConstraintValueTrain<=0,2);
	if(~any(satisfiesConstraint) || isempty(inequalityConstraintValueTrain))
		% if nothing satisfied constraint, keep all for choosing best
		isValidSample = true(size(designSampleTrain,1),1);
	else
		% keep only those that satisfy constraints
		isValidSample = satisfiesConstraint;
	end
	[objectiveOptimal,iOptimalValid] = max(objectiveValueTrain(isValidSample));
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
		internalObjectiveFunction = @(x)-bayesian_acquisition_gaussian_expected_improvement(gaussianRegressionModelObjective,x,objectiveOptimal,explorationFactor);
		internalConstraintFunction = @(x)gaussian_regression_constraint(gaussianRegressionModelInequalityConstraint,gaussianRegressionModelEqualityConstraint,x);

		predictedDesignOptimal = nan(options.BatchTrialSize,nDimension);
		predictedObjectiveOptimal = nan(options.BatchTrialSize,1);
		predictedOptimizationOutput = cell(options.BatchTrialSize,1);
		for i=1:options.BatchTrialSize
			[predictedDesignOptimal(i,:),predictedObjectiveOptimal(i),predictedOptimizationOutput{i}] = ...
				options.InternalOptimizationFunction(internalObjectiveFunction,designSampleNew(i,:),...
					designSpaceLowerBound,designSpaceUpperBound,internalConstraintFunction,...
					options.InternalOptimizationOptions{:});
		end

		% evaluate chosen samples
		[~,iSortedBestPrediction] = sort(predictedObjectiveOptimal,'ascend');
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
		gaussianRegressionModelInequalityConstraintNew = {};
		for i=1:nInequalityConstraint
			gaussianRegressionModelInequalityConstraintNew{i} = fitrgp(designSampleTrainNew,inequalityConstraintValueTrainNew(:,i),...
				options.GaussianRegressionModelTrainingOptions{:});
		end

		equalityConstraintValueTrainNew = [equalityConstraintValueTrain;equalityConstraintEvaluate];
		gaussianRegressionModelEqualityConstraintNew = {};
		for i=1:nEqualityConstraint
			gaussianRegressionModelEqualityConstraintNew{i} = fitrgp(designSampleTrainNew,equalityConstraintValueTrainNew(:,i),...
				options.GaussianRegressionModelTrainingOptions{:});
		end

		% check for convergence 

		% log
		if(isOutputOptimizationInformation)
			optimizationOutput.IterationData(iteration) = struct(...
				'ExplorationFactor',explorationFactor,...
				'BatchSample',designSampleNew,...
				'CandidateSample',predictedDesignOptimal,...
				'PredictedObjective',predictedObjectiveOptimal,...
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
		[objectiveOptimal,iOptimalValid] = max(objectiveValueTrain(isValidSample));
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