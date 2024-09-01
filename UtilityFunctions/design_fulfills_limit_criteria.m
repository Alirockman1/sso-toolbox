function [label,score] = design_fulfills_limit_criteria(measure,lowerLimit,upperLimit)
%DESIGN_FULFILLS_LIMIT_CRITERIA Check if measures of designs are within limit
%	DESIGN_FULFILLS_LIMIT_CRITERIA can return labels for given measures of 
%	designs to classify them as being either within the estabilished limits
%	or not. Additionally, it can also compute a score for said designs, which 
%	is a quantitative value of how much the limits are violated or not.
%
%	LABEL = DESIGN_FULFILLS_LIMIT_CRITERIA(MEASURE) takes the measures MEASURE 
%	resulting from the evaluation of a design, and returns labels LABEL 
%	indicating whether or not all given measures are within the design limits.
%	Labels of 'true' indicate the measures are within the limits, and 'false'
%	indicate at least one is not. In the absence of specification for the lower 
%	and upper limits, it is assumed there is no lower limit, and that the upper  
%	limit is 0.
%
%	LABEL = DESIGN_FULFILLS_LIMIT_CRITERIA(MEASURE,LOWERLIMIT) allows for the 
%	specification of the lower limits LOWERLIMIT for the measures.
%
%	LABEL = DESIGN_FULFILLS_LIMIT_CRITERIA(MEASURE,LOWERLIMIT,UPPERLIMIT) allows
%	for the specification of the upper limits UPPERLIMIT for the measures.
%
%	[LABEL,SCORE] = DESIGN_FULFILLS_LIMIT_CRITERIA(...) also returns the score
%	SCORE for each design. Negative values of score indicate all measures are 
%	within the set limits, and positive indicate the oppositve.
%
%	Input:
%		- MEASURE : (nSample,nRequirement) double
%		- LOWERLIMIT : (1,nRequirement) double, optional
%		- UPPERLIMIT : (1,nRequirement) double, optional
%	Output:
%		- LABEL : (nSample,1) logical
%		- SCORE : (nSample,1) double
%
%   See also design_measure_to_deficit, design_deficit_label_score.
%
%   Copyright 2024 Eduardo Rodrigues Della Noce
%   SPDX-License-Identifier: Apache-2.0

	if(nargin<3)
		upperLimit = [];
		if(nargin<2)
			lowerLimit = [];
		end
	end
	measureDeficit = design_measure_to_deficit(measure,lowerLimit,upperLimit);
	[label,score] = design_deficit_to_label_score(measureDeficit);
end