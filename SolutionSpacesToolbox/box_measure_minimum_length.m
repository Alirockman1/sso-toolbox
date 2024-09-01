function measure = box_measure_minimum_length(designBox, varargin)
%BOX_MEASURE_MINIMUM_LENGTH gives the minimum side length as measure of a box
%   BOX_MEASURE_MINIMUM_LENGTH uses as measure for a candidate design box the 
%   minimum side length (upper boundary - lower boundary) of said box. In an 
%   box-optimization algorithm, this will tend to generate square (all sides 
%   with the same length) solutions.
%   
%   MEASURE = BOX_MEASURE_MINIMUM_LENGTH(DESIGNBOX) receives the candidate 
%   design box in DESIGNBOX and returns its smallest length in MEASURE.
%
%   MEASURE = BOX_MEASURE_MINIMUM_LENGTH(DESIGNBOX,FRACTIONUSEFUL) additionally 
%   receives how many designs inside the current box are useful. For this 
%   measure, this does not affect the result.
%
%   MEASURE = BOX_MEASURE_MINIMUM_LENGTH(...,NAME,VALUE,...) also
%   allows the specification of the additional name-value pair arguments.
%   These are:
%       - 'Weight' : weights for each design variable interval length, where
%       the length is divided by said weight when computing the measure. For 
%       higher weights, that dimension becomes a higher priority. If no value 
%       is given, 1 is assumed for all lengths. 
%
%   Input:
%       - DESIGNBOX : (2,nDesignVariable) double
%           -- (1) : lower boundary of the design box
%           -- (2) : upper boundary of the design box
%       - 'Weight' : (1,nDesignVariable) double
%
%   Output:
%       - MEASURE : double
%
%   See also box_measure_volume.
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
    parser.addOptional('FractionUseful',[]);
    parser.addParameter('Weight',[],@(x)isnumeric(x));
    parser.parse(varargin{:});

    weight = conditional_default_value_assignment(parser.Results.Weight,ones(1,size(designBox,2)));
    
    % edges length
    norm_lengths = (designBox(2,:) - designBox(1,:))./weight;
    measure = min(abs(norm_lengths));
end