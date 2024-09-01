function normalizationFactor = design_measure_normalization_factor(measure,lowerLimit,upperLimit)
%design_measure_normalization_factor Normalization factor for measures
%   design_measure_normalization_factor calculates an appropriate normalization  
%   factor which may be used for the measures of a design sample; it does so 
%   based on the given limits and, if necessary, on the measure acquired.
%   - For values with well-defined limits: upper limit - lower limit;
%   - If either is not suitable: absolute value of the other limit;
%   - If both are not suitable: maximum measure - minimum measure;
%   - If that result is not suitable: absolute value of the mean of measures; 
%   - Last case: 1.
%   
%   NORMALIZATIONFACTOR = design_measure_normalization_factor(MEASURE) receives
%   the measures of the design samples in MEASURE and computes the normalization
%   factor, returning it in NORMALIZATIONFACTOR.
%
%   NORMALIZATIONFACTOR = design_measure_normalization_factor(MEASURE,
%   LOWERLIMIT) allows for the specification of the desired lower limit for the
%   measures in LOWERLIMIT. Default value is 'inf'. If any entry is set to 
%   'nan', that entry will be assumed to be '-inf' as well.
%
%   NORMALIZATIONFACTOR = design_measure_normalization_factor(MEASURE,
%   LOWERLIMIT,UPPERLIMIT) allows for the specification of the desired upper 
%   limit for the measures in UPPERLIMIT. Default value is '0'. If any entry  
%   is set to 'nan', that entry will be assumed to be '+inf'.
%
%   Inputs:
%       - MEASURE : (nSample,nMeasure) double
%       - LOWERLIMIT : (1,nMeasure) double
%       - UPPERLIMIT : (1,nMeasure) double
%
%   Output:
%       - NORMALIZATIONFACTOR : (1,nMeasure) double
%
%   See also design_measure_to_deficit.
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

    if(nargin<3)
        upperLimit = [];
        if(nargin<2)
            lowerLimit = [];
        end
    end
    lowerLimit = ternary(isempty(lowerLimit),-inf,lowerLimit);
    upperLimit = ternary(isempty(upperLimit),0,upperLimit);

    % change NaNs to respective no-limit equivalent
    lowerLimit(isnan(lowerLimit)) = -inf;
    upperLimit(isnan(upperLimit)) = +inf;

    % Calculate Normalization Factor
    normalizationFactor = upperLimit - lowerLimit;
    
    % check if each entry is suitable; if not, try different types
    unsuitableFactor = @(entry) [isinf(entry) || isnan(entry) || (entry==0)];
    for i=1:length(normalizationFactor)
        attempt = 1;
        while(unsuitableFactor(normalizationFactor(i)))
            switch attempt
                case 1
                    normalizationFactor(i) = abs(upperLimit(i));
                case 2
                    normalizationFactor(i) = abs(lowerLimit(i));
                case 3
                    normalizationFactor(i) = max(measure(:,i)) - min(measure(:,i));
                case 4
                    normalizationFactor(i) = abs(mean(measure(:,i)));
                case 5
                    normalizationFactor(i) = 1;
            end
            attempt = attempt + 1;
        end
    end
end