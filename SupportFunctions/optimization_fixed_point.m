function [designOptimal,objectiveOptimal,optimizationOutput] = optimization_fixed_point(objectiveFunction,initialDesign,designSpaceLowerBound,designSpaceUpperBound,constraintFunction,varargin)
	parser = inputParser;
	parser.addParameter('ObjectiveFunctionThresholdValue',-inf);
	parser.addParameter('GroupedEvaluations',[]);
	parser.addParameter('PointsToEvaluate',[]);

	parser.parse(varargin{:});
    options = parser.Results;

   	if(isempty(options.PointsToEvaluate))
   		options.PointsToEvaluate = [initialDesign;designSpaceLowerBound;designSpaceUpperBound];
   	end

   	designSample = options.PointsToEvaluate;
   	nSample = size(designSample,2);
   	if(isempty(options.GroupedEvaluations))
    	options.GroupedEvaluations = ceil(sqrt(nSample));
   	end
   	
   	if(isinf(options.ObjectiveFunctionThresholdValue) || ...
		isempty(options.ObjectiveFunctionThresholdValue) || ...
		isnan(options.ObjectiveFunctionThresholdValue))
		
		% evaluate all maximum samples together
		objectiveValue = objectiveFunction(designSample);
		constraintValue = constraintFunction(designSample);

		satisfiesConstraint = all(constraintValue<=0,2);
		if(~any(satisfiesConstraint))
			validSample = true(nSample,1);
		else
			validSample = satisfiesConstraint;
		end
		[objectiveOptimal,indexOptimalValid] = min(objectiveValue(validSample));
		indexOptimal = convert_index_base(validSample,indexOptimalValid,'backward');
		designOptimal = designSample(indexOptimal,:);
	else
		objectiveValue = nan(nSample,1);
		constraintValue = [];

		iStart = 1;
		iEnd = min([iStart+options.GroupedEvaluations-1,nSample]);
		foundOptimal = false;
		while(~foundOptimal && iStart<=nSample)
			iEvaluate = iStart:iEnd;

			objectiveValue(iEvaluate) = objectiveFunction(designSample(iEvaluate,:));
			if(ismember(1,iEvaluate))
				constraintValue1 = constraintFunction(designSample(iEvaluate,:));

			    constraintValue = nan(nSample,size(constraintValue1,2));
			    constraintValue(iEvaluate,:) = constraintValue1;
			    clear constraintValue1
			else
			    constraintValue(iEvaluate,:) = constraintFunction(designSample(iEvaluate,:));
			end

			satisfiesConstraint = all(constraintValue<=0,2);
			if(~any(satisfiesConstraint))
				% if nothing satisfied constraint, keep all for choosing best
				validSample = true(nSample,1); 
			else
				% keep only those that satisfy constraints
				validSample = satisfiesConstraint; 
			end
			[objectiveOptimal,indexOptimalValid] = min(objectiveValue(validSample));
			indexOptimal = convert_index_base(validSample,indexOptimalValid,'backward');
			designOptimal = designSample(indexOptimal,:);
			if(objectiveOptimal<=options.ObjectiveFunctionThresholdValue && satisfiesConstraint(indexOptimal))
				foundOptimal = true;
				designSample(iEvaluate(end)+1:end,:) = [];
				objectiveValue(iEvaluate(end)+1:end) = [];
				constraintValue(iEvaluate(end)+1:end,:) = [];
			else
				iStart = iEnd+1;
				iEnd = min([iStart+options.GroupedEvaluations-1,nSample]);
			end
		end
	end

	optimizationOutput.EvaluatedSamples = designSample;
	optimizationOutput.ObjectiveValue = objectiveValue;
	optimizationOutput.ConstraintValue = constraintValue;
end
