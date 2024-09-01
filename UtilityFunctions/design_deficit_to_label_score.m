function [label,score] = design_deficit_to_label_score(measureDeficit,deficitWeight)
%DESIGN_DEFICIT_TO_LABEL_SCORE Get the label and score for designs given deficit
%   DESIGN_DEFICIT_TO_LABEL_SCORE uses the given measure deficits of design 
%   samples to label them as either fulfilling the given requirements or not. It
%   can also produce a quantitative value to represent how much said measures 
%   are within their supposed limits.
%
%   LABEL = DESIGN_DEFICIT_TO_LABEL_SCORE(MEASUREDEFICIT) gets the deficits of 
%   the measures of design samples in MEASUREDEFICIT and labels them in LABEL as 
%   fulfilling the requirements ('true') if all deficits are negative, or not 
%   ('false') if at least one deficit is positive.
%
%   [LABEL,SCORE] = DESIGN_DEFICIT_TO_LABEL_SCORE(...) also returns a SCORE for 
%   each design. Said score is a negative value for designs that meet the 
%   requirements, and positive value if one requirement is violated. In 
%   particular:
%       - For designs that meet the requirements, the score is computed as a 
%       worst-case function (maximum value of deficit - 'weakest link').
%       - For designs that do not meet at least one requirement, the score is
%       computed as the root square sum of all negative deficits (so all 
%       deficits of measures that violated the requirements).
%
%   [LABEL,SCORE] = DESIGN_DEFICIT_TO_LABEL_SCORE(MEASUREDEFICIT,DEFICITWEIGHT)
%   allows one to choose weights for each deficit, which determines how much
%   they are going to influence the score. By default, all deficits are given 
%   weight equal to 1.
%
%   Input:
%       - MEASUREDEFICIT : (nSample,nRequirement) double
%   Output:
%       - LABEL : (nSample,1) logical
%       - SCORE : (nSample,1) double
%
%   See also DesignEvaluatorBottomUpMapping.
%
%   Copyright 2024 Eduardo Rodrigues Della Noce
%   SPDX-License-Identifier: Apache-2.0

    % set default weight
    nDeficit  = size(measureDeficit,2);
    if(nargin<2 || isempty(deficitWeight))
        deficitWeight = 1;
    end

    % calculate score based on worst-case objective meta function and 
    % label accordingly
    weightedDeficit = deficitWeight.*measureDeficit;
    score = max(weightedDeficit,[],2);
    label = (score<=0);

    % for bad designs, use the sum of violated limits instead of worst-case
    violatedLimit = max(weightedDeficit,0);
    score(~label) = sqrt(sum(violatedLimit(~label,:).^2,2));
end