function [designOptimal,objectiveOptimal,optimizationOutput] = optimization_sampling(objectiveFunction,initialDesign,designSpaceLowerBound,designSpaceUpperBound,constraintFunction,varargin)
%OPTIMIZATION_SAMPLING for bottom-up mapping / evaluator optimization
%	OPTIMIZATION_SAMPLING uses a brute-force optimization approach where
%	it samples the design space and returns the best value from the tested
%	options.
%
%	DESIGNOPTIMAL = OPTIMIZATION_SAMPLING(OBJECTIVEFUNCTION,INITIALDESIGN,
%	DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND,CONSTRAINTFUNCTION) minimizes 
%	OBJECTIVEFUNCTION starting on INITIALDESIGN, on the design space defined by 
%	DESIGNSPACELOWERBOUND and DESIGNSPACEUPPERBOUND, subject to the constraint 
%	functions CONSTRAINTFUNCTION being less than or equal to zero, returning the
%	optimal design DESIGNOPTIMAL. 
%
%	DESIGNOPTIMAL = OPTIMIZATION_SAMPLING(...NAME,VALUE,...) allows for the
%	specification of additional options. These include:
%		- 'ObjectiveFunctionThresholdValue' : if a design is evaluated and has
%		an objective value which is lower than this threshold and constraints 
%		are satisfied, the optimization procedure stops and that design is 
%		returned. Default: -inf.
%		- 'MaxFunctionEvaluations' : maximum number of function evaluations. 
%		Default value: 1000.
%		- 'GroupedEvaluations' : if a threshold value for the objective has been
%		set, then for each iteration, only 'GroupedEvaluations' samples are
%		evaluated at once. This may allow for better parallelization, while
%		still being able to stop the optimization if a suitable design is found.
%		For example, if the objective function or constraint functions have
%		an expensive process that can be parallelized for each CPU core, then 
%		setting this number to your number of cores may make sense.
%		Default is the square root of the maximum number of function 
%		evaluations (which will lead to also the same number of iterations, at 
%		most).
%		- 'SamplingFunction' : function handle on how to perform the sampling.
%		Default: @sampling_latin_hypercube.
%		- 'SamplingOptions' : options to be used with the sampling function.
%		Default is empty.
%
%	[DESIGNOPTIMAL,OBJECTIVEOPTIMAL] = OPTIMIZATION_SAMPLING(...) also
%	returns the value of the objective function OBJECTIVEOPTIMAL for the optimal
%	design.
%
%	[DESIGNOPTIMAL,OBJECTIVEOPTIMAL,OPTIMIZATIONOUTPUT] = 
%	OPTIMIZATION_SAMPLING(...) returns additional information regarding
%	the optimization procedure. In particular:
%		- EvaluatedSamples : design sample points evaluated
%		- ObjectiveValue : objective value for each design sample point
%		- ConstraintValue : constraint values for each design sample point
%
%   Input:
%		- OBJECTIVEFUNCTION : function_handle
%		- INITIALDESIGN : (1,nDesignVariable) double
%		- DESIGNSPACELOWERBOUND : (1,nDesignVariable) double
%		- DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double
%		- CONSTRAINTFUNCTION : function_handle
%		- 'ObjectiveFunctionThresholdValue' : double
%		- 'MaxFunctionEvaluations' : integer
%		- 'GroupedEvaluations' : integer
%		- 'SamplingFunction' : function_handle
%		- 'SamplingOptions' : (1,nOptions) cell
%
%   Output:
%		- DESIGNOPTIMAL : (1,nDesignVariable) double
%		- OBJECTIVEOPTIMAL : double
%		- OPTIMIZATIONOUTPUT : structure
%			-- EvaluatedSamples : (nSampleEvaluated,nDesignVariable) double
%			-- ObjectiveValue : (nSampleEvaluated,1) double
%			-- ConstraintValue : (nSampleEvaluated,nConstraint) double
%
%   See also design_optimize_quantities_of_interest, 
%	design_optimize_performance_score.
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
	parser.addParameter('ObjectiveFunctionThresholdValue',-inf);
	parser.addParameter('MaxFunctionEvaluations',1000);
	parser.addParameter('GroupedEvaluations',[]);
	parser.addParameter('SamplingFunction',@sampling_latin_hypercube);
	parser.addParameter('SamplingOptions',{});

	parser.parse(varargin{:});
    options = parser.Results;

    if(isempty(options.GroupedEvaluations))
    	options.GroupedEvaluations = ceil(sqrt(options.MaxFunctionEvaluations));
   	end

   	nSample = options.MaxFunctionEvaluations;
	designSample = options.SamplingFunction([designSpaceLowerBound;designSpaceUpperBound],...
		nSample-1,options.SamplingOptions{:});
	designSample = [initialDesign;designSample];

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
		iEnd = min([iStart+options.GroupedEvaluations-1,options.MaxFunctionEvaluations]);
		foundOptimal = false;
		while(~foundOptimal && iStart<=options.MaxFunctionEvaluations)
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
				validSample = true(nSample,1);
			else
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
				iEnd = min([iStart+options.GroupedEvaluations-1,options.MaxFunctionEvaluations]);
			end
		end
	end

	optimizationOutput.EvaluatedSamples = designSample;
	optimizationOutput.ObjectiveValue = objectiveValue;
	optimizationOutput.ConstraintValue = constraintValue;
end
