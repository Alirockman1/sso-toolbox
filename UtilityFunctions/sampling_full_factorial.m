function designSample = sampling_full_factorial(designBox,nSample,varargin)
%SAMPLING_FULL_FACTORIAL Full-factorial sampling inside a design box
%   SAMPLING_FULL_FACTORIAL generates a full factorial design of experiments in 
%   the given desired boundaries.
%   The number of designs desired must be compatible with the number of 
%   dimensions of the problem. In other words, (number of samples)^(1/dimension)
%   must be an integer.
%   Optionally, one can also perform a condensed full factorial design of 
%   experiments; by 'condensed', it means the extremities are not used as part  
%   of the places to sample. For example, for a sample of a single variable with
%   two levels between 0 and 1, full factorial would produce the sample [0,1]; 
%   the condensed full factorial instead produces [1/3, 2/3]. For three levels, 
%   full factorial would produce [0, 0.5, 1], and condensed full factorial 
%   [0.25,0.50,0.75], and so on.
%
%   DESIGNSAMPLE = SAMPLING_FULL_FACTORIAL(DESIGNBOX,NSAMPLE) produces NSAMPLE
%   design samples inside design box DESIGNBOX, returning that as DESIGNSAMPLE.
%   Said samples are distributed equally by levels for each design variable, 
%   with the number of levels of the full factorial being computed as 
%   NSAMPLE^(1/nDesignVariable).
%
%   DESIGNSAMPLE = SAMPLING_FULL_FACTORIAL(...NAME,VALUE,...) allows
%   for the use of the condensed full factorial if the property 'Condensed' is
%   set to 'true'.
%
%   Inputs:
%       - DESIGNBOX : (2,nDesignVariable) double
%           -- (1) : lower boundary of the design box
%           -- (2) : upper boundary of the design box
%       - NSAMPLE : integer
%       - 'Condensed' : logical
%
%   Output:
%       - DESIGNSAMPLE : (nSample,nDesignVariable) double
%
%   See also fullfact, nthroot, sampling_full_factorial_condensed.
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

    options = parser_variable_input_to_structure(varargin{:});

    % options processing
    defaultOptions = {'Condensed',false};
    [~,sampleOptions] = merge_name_value_pair_argument(defaultOptions,options);

    % Problem Dimension
    nDesignVariable = size(designBox,2);
    
    % Determine if it is possible to generate nSample with problem of size dim
    % nSample = levels^d -> levels = nSample^(1/d) should be an integer
    approximatedLevel = nthroot(nSample,nDesignVariable);
    fullFactorialLevel = round(approximatedLevel);
    if(fullFactorialLevel^nDesignVariable ~= nSample)
        error('SFF:NotPossible',['Cannot make uniform full factorial distribution with given sample size\n',...
            'Requested number of samples: %d\nProblem dimension: %d\nApproximate number of levels: %g'],...
            nSample,nDesignVariable,approximatedLevel);
    end
    
    % Design of Experiments
    if(sampleOptions.Condensed)
        designOfExperiments = (fullfact(fullFactorialLevel*ones(nDesignVariable,1)))/(fullFactorialLevel+1);
    else
        designOfExperiments = (fullfact(fullFactorialLevel*ones(nDesignVariable,1))-1)/(fullFactorialLevel-1);
    end

    % Rescale to desired region
    designSample = designBox(1,:) + (designBox(2,:)-designBox(1,:)).*designOfExperiments;
end