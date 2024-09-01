function [scoreThreshold,iterData] = design_fast_forward_find_uncertainty_score(bottomUpMapping,designFastForward,varargin)
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
	
	parser = inputParser;
	parser.addParameter('MeasureChoice','solution',@(x)ischar(x)||isstring(x));
	parser.addParameter('SamplingFunction',@sampling_random,@(x)isa(x,'function_handle'));
    parser.addParameter('SamplingOptions',{},@(x)iscell(x));
    parser.addParameter('NumberCandidateSamples',100000);
    parser.addParameter('NumberNeighborEvaluations',3);

	parser.parse(varargin{:});
	options = parser.Results;

	% decide if performance or physical feasibility will be checked
	usePhysicalFeasibility = false(2,1);
    usePhysicalFeasibility(1) = (isnumeric(options.MeasureChoice) && options.MeasureChoice==2);
    usePhysicalFeasibility(2) = ((ischar(options.MeasureChoice) || isstring(options.MeasureChoice)) ...
        && strcmpi(erase(options.MeasureChoice,{'-','_'}),'physicalfeasibility'));
    usePerformance = ~any(usePhysicalFeasibility);

    isOutputIterData = (nargout>=2);
    if(isOutputIterData)
        iterData = struct(...
            'EvaluatedSamples',[],...
            'PerformanceMeasure',[],...
            'PhysicalFeasibilityMeasure',[],...
            'BottomUpMappingOutput',[],...
            'DeficitPredicted',[],...
            'ScorePredicted',[],...
            'LabelPredicted',[],...
            'ScoreEvaluated',[],...
            'LabelEvaluated',[]);
        iteration = 1;
    end

	% generate large sample
	trainingSample = designFastForward.DesignSampleTrain;
	labelSample = designFastForward.LabelTrain;
	[~,samplingBox] = design_bounding_box(trainingSample,labelSample);

	designSample = options.SamplingFunction(samplingBox,options.NumberCandidateSamples,options.SamplingOptions{:});

	% predict samples according to fast forward
	deficitPredicted = designFastForward.predict_designs(designSample);
	[labelPredicted,scorePredicted] = design_deficit_to_label_score(deficitPredicted);

	% order test data from lowest to highest absolute score
	[~,iOrder] = sort(abs(scorePredicted));
	designSample = designSample(iOrder,:);
	deficitPredicted = deficitPredicted(iOrder,:);
	scorePredicted = scorePredicted(iOrder);
	labelPredicted = labelPredicted(iOrder);

	% allocate vector to save results already evaluated
	labelEvaluated = nan(options.NumberCandidateSamples,1);
	scoreEvaluated = nan(options.NumberCandidateSamples,1);

	% perform binary search
	iStart = 1;
	iEnd = size(scorePredicted,1);
	converged = false;
	while (~converged)
        % update current reference 
		iMiddle = iStart + floor((iEnd-iStart)/2);

        % determine if any of these points have been evaluated before; if not, don't evaluate again
		iCandidateEvaluate = ((iMiddle-options.NumberNeighborEvaluations):(iMiddle+options.NumberNeighborEvaluations))';
		validEvaluation = (iCandidateEvaluate>=1) & (iCandidateEvaluate<=options.NumberCandidateSamples);
        iCheckError = iCandidateEvaluate(validEvaluation);
		needsEvaluation = (isnan(scoreEvaluated(iCheckError)));
        iEvaluate = iCheckError(needsEvaluation);
        
		% evaluate points
		if(any(needsEvaluation))
			[performanceMeasure,physicalFeasibilityMeasure,bottomUpMappingOutput] = ...
				bottomUpMapping.response(designSample(iEvaluate,:));
	        measure = ternary(usePerformance,performanceMeasure,physicalFeasibilityMeasure);
	        [labelEvaluated(iEvaluate),scoreEvaluated(iEvaluate)] = ...
	        	designFastForward.classify_measure(measure);
        else
            performanceMeasure = [];
            physicalFeasibilityMeasure = [];
            bottomUpMappingOutput = [];
	    end

        % see if errors occured
        isError = (labelPredicted(iCheckError) ~= labelEvaluated(iCheckError));

        % save evaluations
        if(isOutputIterData)
	        iterData(iteration) = struct(...
	            'EvaluatedSamples',designSample(iEvaluate,:),...
	            'PerformanceMeasure',performanceMeasure,...
	            'PhysicalFeasibilityMeasure',physicalFeasibilityMeasure,...
	            'BottomUpMappingOutput',bottomUpMappingOutput,...
                'DeficitPredicted',deficitPredicted(iEvaluate),...
	            'ScorePredicted',scorePredicted(iEvaluate),...
	            'LabelPredicted',labelPredicted(iEvaluate),...
	            'ScoreEvaluated',scoreEvaluated(iEvaluate),...
	            'LabelEvaluated',labelEvaluated(iEvaluate));
	        iteration = iteration + 1;
	    end

        % update search region
	    if(any(isError))
	    	iStart = iMiddle;
	    else
	    	iEnd = iMiddle;
	    end
	    converged = (iEnd-iStart<=1);
    end
    scoreThreshold = abs(scorePredicted(iEnd));
end