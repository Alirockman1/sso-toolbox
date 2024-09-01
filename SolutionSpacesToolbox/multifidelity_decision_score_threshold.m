function isUncertain =  multifidelity_decision_score_threshold(designSample,performanceDeficit,physicalFeasibilityDeficit,varargin);
%MULTIFIDELITY_DECISION_SCORE_THRESHOLD Decide uncertainty based on score value
%	MULTIFIDELITY_DECISION_SCORE_THRESHOLD decides whether or not the result of
%	a design evaluator is uncertain in nature (correctly classifies a design
%	as being good/bad or physically feasible/infeasible) based on its score.
%
%	ISUNCERTAIN =  MULTIFIDELITY_DECISION_SCORE_THRESHOLD(DESIGNSAMPLE,
%	PERFORMANCEDEFICIT,PHYSICALFEASIBILITYDEFICIT) receives the design sample 
%	points in DESIGNSAMPLE, and the results of the design evaluator in terms of 
%	performance PERFORMANCEDEFICIT and physical feasibility 
%	PHYSICALFEASIBILITYDEFICIT, and uses that to decide for each sample point
%	if the result is uncertain or not, returning the logical flag ISUNCERTAIN.
%	By default, nothing will be considered uncertain.
%
%	ISUNCERTAIN =  MULTIFIDELITY_DECISION_SCORE_THRESHOLD(...NAME,VALUE,...)
%	allows the specification of additional parameters as name-value pair 
%	arguments. These are:
%		- 'PerformanceDeficitWeight' : weight when computing the score from the
%		deficits for performance. Default is empty.
%		- 'PerformanceScoreThreshold' : either: an absolute value of score, 
%		where any design that has a score in absolute value lower than that is 
%		considered uncertain; or two values, the first being a lower threshold
%		and the second an upper threshold, where any design that has a score
%		between those two is considered uncertain. Default value is empty.
%		- 'PerformanceScoreLowerThreshold' : lower threshold of score for 
%		performance. Default value is empty.
%		- 'PerformanceScoreUpperThreshold' : upper threshold of score for 
%		performance. Default value is empty.
%		- 'PhysicalFeasibilityDeficitWeight' : weight when computing the score 
%		from the deficits for physical feasibility. Default is empty.
%		- 'PhysicalFeasibilityScoreThreshold' : either: an absolute value of 
%		score, where any design that has a score in absolute value lower than
%		that is considered uncertain; or two values, the first being a lower 
%		threshold and the second an upper threshold, where any design that
%		has a score between those two is considered uncertain. Default value is 
%		empty.
%		- 'PhysicalFeasibilityScoreLowerThreshold' : lower threshold of score  
%		for physical feasibility. Default value is empty.
%		- 'PhysicalFeasibilityScoreUpperThreshold' : upper threshold of score  
%		for physical feasibility. Default value is empty.
%
%   Input:
%		- DESIGNSAMPLE : (nSample,nDesignVariable) double
%		- PERFORMANCEDEFICIT : (nSample,nPerformance) double
%		- PHYSICALFEASIBILITYDEFICIT : (nSample,nPhysicalFeasibility)
%		double
%		- EVALUATOROUTPUTLOWFIDELITY : class-dependent
%		- 'PerformanceDeficitWeight' : (1,nPerformance) double
%		- 'PerformanceScoreThreshold' : double OR (1,2) double
%		- 'PerformanceScoreLowerThreshold' : double
%		- 'PerformanceScoreUpperThreshold' : double
%		- 'PhysicalFeasibilityDeficitWeight' : (1,nPhysicalFeasibility) double
%		- 'PhysicalFeasibilityScoreThreshold' : double OR (1,2) double
%		- 'PhysicalFeasibilityScoreLowerThreshold' : double
%		- 'PhysicalFeasibilityScoreUpperThreshold' : double
%
%   Output:
%		- ISUNCERTAIN : (nSample,1) logical
%
%   See also DesignEvaluatorMultiFidelity, 
%	design_fast_forward_find_uncertainty_score.
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
	parser.StructExpand = false;
	parser.addOptional('EvaluatorOutputLowFidelity',[]);
	parser.addParameter('PerformanceDeficitWeight',[]);
	parser.addParameter('PerformanceScoreThreshold',[]);
	parser.addParameter('PerformanceScoreLowerThreshold',[]);
	parser.addParameter('PerformanceScoreUpperThreshold',[]);
	parser.addParameter('PhysicalFeasibilityDeficitWeight',[]);
	parser.addParameter('PhysicalFeasibilityScoreThreshold',[]);
	parser.addParameter('PhysicalFeasibilityScoreLowerThreshold',[]);
	parser.addParameter('PhysicalFeasibilityScoreUpperThreshold',[]);
	parser.parse(varargin{:});
	options = parser.Results;

	nSample = size(designSample,1);


	% performance uncertainty
    if(~isempty(options.PerformanceScoreThreshold))
	    if(length(options.PerformanceScoreThreshold) == 1)
		    options.PerformanceScoreThreshold(2) = options.PerformanceScoreThreshold(1);
		    options.PerformanceScoreThreshold(1) = -options.PerformanceScoreThreshold(1);
	    end
	    if(isempty(options.PerformanceScoreLowerThreshold))
		    options.PerformanceScoreLowerThreshold = options.PerformanceScoreThreshold(1);
	    end
	    if(isempty(options.PerformanceScoreUpperThreshold))
		    options.PerformanceScoreUpperThreshold = options.PerformanceScoreThreshold(2);
        end
    end

	if(isempty(options.PerformanceScoreLowerThreshold))
		options.PerformanceScoreLowerThreshold = -inf;
	end
	if(isempty(options.PerformanceScoreUpperThreshold))
		options.PerformanceScoreUpperThreshold = +inf;
	end

	if(isinf(options.PerformanceScoreLowerThreshold) && isinf(options.PerformanceScoreUpperThreshold))
		isUncertainPerformance = false(nSample,1);
	else
		[~,scorePerformance] = design_deficit_to_label_score(...
			performanceDeficit,options.PerformanceDeficitWeight);
		isUncertainPerformance = (scorePerformance>=options.PerformanceScoreLowerThreshold) & ...
			(scorePerformance<=options.PerformanceScoreUpperThreshold);
	end


	% physical feasibility uncertainty
    if(~isempty(options.PhysicalFeasibilityScoreThreshold))
	    if(length(options.PhysicalFeasibilityScoreThreshold) == 1)
		    options.PhysicalFeasibilityScoreThreshold(2) = options.PhysicalFeasibilityScoreThreshold(1);
		    options.PhysicalFeasibilityScoreThreshold(1) = -options.PhysicalFeasibilityScoreThreshold(1);
	    end
	    if(isempty(options.PhysicalFeasibilityScoreLowerThreshold))
		    options.PhysicalFeasibilityScoreLowerThreshold = options.PhysicalFeasibilityScoreThreshold(1);
	    end
	    if(isempty(options.PhysicalFeasibilityScoreUpperThreshold))
		    options.PhysicalFeasibilityScoreUpperThreshold = options.PhysicalFeasibilityScoreThreshold(2);
        end
    end

	if(isempty(options.PhysicalFeasibilityScoreLowerThreshold))
		options.PhysicalFeasibilityScoreLowerThreshold = -inf;
	end
	if(isempty(options.PhysicalFeasibilityScoreUpperThreshold))
		options.PhysicalFeasibilityScoreUpperThreshold = +inf;
	end

	if(isinf(options.PhysicalFeasibilityScoreLowerThreshold) && isinf(options.PhysicalFeasibilityScoreUpperThreshold))
		isUncertainPhysicalFeasibility = false(nSample,1);
	else
		[~,scorePhysicalFeasibility] = design_deficit_to_label_score(...
			physicalFeasibilityDeficit,options.PhysicalFeasibilityDeficitWeight);
		isUncertainPhysicalFeasibility = (scorePhysicalFeasibility>=options.PhysicalFeasibilityScoreLowerThreshold) & ...
			(scorePhysicalFeasibility<=options.PhysicalFeasibilityScoreUpperThreshold);
	end

	% aggregate
	isUncertain = isUncertainPerformance | isUncertainPhysicalFeasibility;
end