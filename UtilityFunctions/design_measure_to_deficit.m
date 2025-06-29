function [deficit,normalizationFactor] = design_measure_to_deficit(measure,lowerLimit,upperLimit,normalizationFactor)
%DESIGN_MEASURE_TO_DEFICIT Deficit of design given measures and the requirements
%   DESIGN_MEASURE_TO_DEFICIT is used to compute the (performance/physical
%   feasibility) deficit of a design given its (performance/physical 
%   feasibility) measures previously computed and the limits imposed on the
%   system.
%   A normalization of the measures is done, and the deficit can is computed
%   in the following way:
%       - lower limit deficit: (lower limit - measure) / normalization factor
%       - upper limit deficit: (measure - upper limit) / normalization factor
%       - deficit : maximum between lower limit and upper limit deficit
%   In this way, measures that are within their given limits are given deficit 
%   values <0 and otherwise deficit values >0.
%
%   DEFICIT = DESIGN_MEASURE_TO_DEFICIT(MEASURE) receives the measures of design
%   sample points in MEASURE and returns the deficits of each measure in 
%   DEFICIT. All other inputs assume their default values, where lower limits 
%   are assumed to not apply, upper limits are assumed to be 0 and the 
%   normalization factor is computed on-thy-fly based on the given measures and
%   limits.
%
%   DEFICIT = DESIGN_MEASURE_TO_DEFICIT(MEASURE,LOWERLIMIT) allows for the 
%   specification of minimum values required for the measures in LOWERLIMIT.
%   Its default value is '-inf'. If any entry is set to 'nan', that entry will 
%   be assumed to be '-inf' as well.
%
%   DEFICIT = DESIGN_MEASURE_TO_DEFICIT(MEASURE,LOWERLIMIT,UPPERLIMIT) allows for  
%   the specification of maximum values required for the measures in UPPERLIMIT.
%   Its default value is '0'. If any entry is set to 'nan', that entry will 
%   be assumed to be '+inf'.
%
%   DEFICIT = DESIGN_MEASURE_TO_DEFICIT(MEASURE,LOWERLIMIT,UPPERLIMIT,
%   NORMALIZATIONFACTOR) also allows for the specification of the normalization
%   factors to be used in deficit computation in NORMALIZATIONFACTOR. Default
%   value is 'nan'. If any entry is set to 'nan', that entry will be computed
%   on-the-fly based on the given measures and limits.
%
%   Inputs:
%       - MEASURE :  (nSample,nMeasure) double
%       - LOWERLIMIT : (1,nMeasure) double, optional
%       - UPPERLIMIT : (1,nMeasure) double, optional
%       - NORMALIZATIONFACTOR : (1,nMeasure) double, optional
%
%   Output:
%       - DEFICIT : (nSample,nMeasure) double
%
%   See also design_measure_normalization_factor.
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

    if(nargin<4)
        normalizationFactor = [];
        if(nargin<3)
            upperLimit = [];
            if(nargin<2)
                lowerLimit = [];
            end
        end
    end
    normalizationFactor = conditional_default_value_assignment(normalizationFactor,nan(1,size(measure,2)));
    lowerLimit = conditional_default_value_assignment(lowerLimit,-inf);
    upperLimit = conditional_default_value_assignment(upperLimit,0);

    % change NaNs to respective no-limit equivalent
    lowerLimit(isnan(lowerLimit)) = -inf;
    upperLimit(isnan(upperLimit)) = +inf;

    nMeasure = size(measure,2);
    if(size(lowerLimit,2)==1)
        lowerLimit = lowerLimit*ones(1,nMeasure);
    end
    if(size(upperLimit,2)==1)
        upperLimit = upperLimit*ones(1,nMeasure);
    end

    % Find lower and upper scores, as well as normalization factor
    computeNormalization = (isinf(normalizationFactor) | isnan(normalizationFactor) | (normalizationFactor==0));
    if(any(computeNormalization))
        normalizationFactor(computeNormalization) = design_measure_normalization_factor(...
            measure(:,computeNormalization),...
            lowerLimit(computeNormalization),...
            upperLimit(computeNormalization));
    end

    % compute deficits for each limit
    lowerLimitDeficit = (lowerLimit - measure)./normalizationFactor;
    upperLimitDeficit = (measure - upperLimit)./normalizationFactor;

    % take worst-case for each 
    deficit = max(lowerLimitDeficit,upperLimitDeficit);
end
