function measure = box_measure_volume(designBox, varargin)
%BOX_MEASURE_VOLUME uses the design box area/volume/hypervolume as its measure 
%   BOX_MEASURE_VOLUME uses as measure for a candidate design box its (useful) 
%   volume. In an box-optimization algorithm, this will tend to generate 
%   the largest solutions in terms of the cartesian products of each design 
%   variable.
%
%   MEASURE = BOX_MEASURE_VOLUME(DESIGNBOX) receives the design box in DESIGNBOX
%   and returns its volume (calculated as the product of its lengths) in 
%   MEASURE.
%
%   MEASURE = BOX_MEASURE_VOLUME(DESIGNBOX,FRACTIONUSEFUL) additionally receives
%   the fraction of the volume which is useful in FRACTIONUSEFUL. The resulting
%   measure computation will be multiplied by this factor. By default, it 
%   assumes the value of 1.
%
%   MEASURE = BOX_MEASURE_MINIMUM_LENGTH(...,NAME,VALUE,...) also
%   allows the specification of the additional name-value pair arguments.
%   These are:
%       - 'Weight' : weights for each design variable interval length, where
%       each length is exponentiated by its weight. For higher weights, that 
%       dimension becomes a higher priority. If no value is given, 1 is assumed 
%       for all intervals. 
%
%   Input:
%       - DESIGNBOX : (2,nDesignVariable) double
%           -- (1) : lower boundary of the design box
%           -- (2) : upper boundary of the design box
%       - FRACTIONUSEFUL : double
%       - 'Weight' : (1,nDesignVariable) double
%
%   Output:
%       - MEASURE : double
%
%   See also sso_box_stochastic, box_measure_minimum_length.
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

    volumeFactor = conditional_default_value_assignment(parser.Results.FractionUseful,1);
    weight = conditional_default_value_assignment(parser.Results.Weight,ones(1,size(designBox,2)));
    
    % volume
    weightedLengths = (designBox(2,:) - designBox(1,:)).^weight;
    measure = volumeFactor * prod(weightedLengths);
end