function designSample = sampling_random(designBox,nSample,varargin)
%SAMPLING_RANDOM Full-factorial sampling inside a design box
%   SAMPLING_RANDOM generates a uniformly distributed random design of 
%   experiments in the given desired boundaries.
%
%   DESIGNSAMPLE = SAMPLING_RANDOM(DESIGNBOX,NSAMPLE) produces NSAMPLE
%   design samples inside design box DESIGNBOX, returning that as DESIGNSAMPLE.
%
%   DESIGNSAMPLE = SAMPLING_RANDOM(...NAME,VALUE,...) allows for the 
%   specification of options for the design of experiments using 'rand'.
%
%   Inputs:
%       - DESIGNBOX : (2,nDesignVariable) double
%           -- (1) : lower boundary of the design box
%           -- (2) : upper boundary of the design box
%       - NSAMPLE : integer
%       - Name-value pair options: passed directly to 'rand'.
%
%   Output:
%       - DESIGNSAMPLE : (nSample,nDesignVariable) double
%
%   See also rand, sampling_latin_hypercube.

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
    defaultOptions = {};
    [~,sampleOptions] = merge_name_value_pair_argument(defaultOptions,options);

    % Problem Dimension
    nDesignVariable = size(designBox,2);
    
    % Design of Experiments
    designOfExperiments = rand(nSample,nDesignVariable,sampleOptions{:});

    % Rescale to desired region
    designSample = designBox(1,:) + (designBox(2,:)-designBox(1,:)).*designOfExperiments;
end