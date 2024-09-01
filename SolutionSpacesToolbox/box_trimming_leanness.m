function candidateBoxLean = box_trimming_leanness(designSample,labelKeep,trimmingOrder,varargin)
%box_apply_leanness Application of leanness condition for design boxes
%   box_apply_leanness applies the leanness condition of requirement spaces to a
%   given design box, ensuring that there are useful designs in all edges.
%
%   CANDIDATEBOXLEAN = box_apply_leanness(DESIGNBOX,DESIGNSAMPLE,LABELUSE) receives
%   the current design box in DESIGNBOX, the design samples in DESIGNSAMPLE, and
%   their respective usefulness labels in LABELUSE, and returns the design box
%   after the leanness condition is applied in CANDIDATEBOXLEAN. For LABELUSE,
%   values of 'true' indicate that the respective sample is useful, and 'false'
%   indicates the opposite.
%
%   Inputs:
%       - DESIGNBOX : (2,nDesignVariable) double
%           -- (1) : lower boundary of the design box
%           -- (2) : upper boundary of the design box
%       - DESIGNSAMPLE : (nSample,nDesignVariable) double
%       - LABELUSE : (nSample,1) logical
%
%   Outputs:
%       - CANDIDATEBOXLEAN : (2,nDesignVariable) double
%           -- (1) : lower boundary of the design box
%           -- (2) : upper boundary of the design box
%
%   See also design_bounding_box, box_trimming_classic, sso_box_stochastic.
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

    % parse inputs
    parser = inputParser;
    parser.addRequired('designSample',@(x)isnumeric(x));
    parser.addRequired('labelKeep',@(x)islogical(x));
    parser.addRequired('trimmingOrder',@(x)(isnumeric(x)&&(size(x,2)==1))||isempty(x));
    parser.addOptional('CandidateBox',[],@(x)isnumeric(x)&&(size(x,1)==2 || isempty(x)));
    parser.addParameter('Slack',0.5,@(x)isnumeric(x)&&isscalar(x)&&(x>=0)&&(x<=1));
    
    parser.parse(designSample,labelKeep,trimmingOrder,varargin{:});
    options = parser.Results;

    % assign conditional default values
    candidateBoxLean = conditional_default_value_assignment(options.CandidateBox,design_bounding_box(designSample));
    isInsideBox = is_in_design_box(designSample,candidateBoxLean);

    % if all designs inside box are useful, no need to apply leanness
    if(all(labelKeep(isInsideBox)))
        return;
    else
        isInsideBoxKeep = isInsideBox & labelKeep;
    end

    % get sizes
    nDesignVariable = size(designSample,2);
    nExclude = size(trimmingOrder,1);
    
    % if not, trim to closest useful design
    for i=1:nExclude
        iExclude = trimmingOrder(i);
        if(~isInsideBox(iExclude))
            continue;
        end
        
        for j=1:nDesignVariable
            for k=[1,2]
                candidateBoxCurrent = candidateBoxLean;

                % move lower/upper bound
                candidateBoxCurrent(k,j) = designSample(iExclude,j);
                isInsideBoxCurrent = is_in_design_box(designSample,candidateBoxCurrent);
                isInsideBoxKeepCurrent = isInsideBoxCurrent & labelKeep;

                % only trim if no useful design that should be kept is trimmed
                if(any(isInsideBoxKeepCurrent~=isInsideBoxKeep))
                    continue;
                else
                    candidateBoxLean = candidateBoxCurrent;
                    isInsideBox = isInsideBoxCurrent;
                end
            end
        end
    end

    % use slack
    strictBox = design_bounding_box(designSample,isInsideBox & labelKeep);
    candidateBoxLean = candidateBoxLean*options.Slack + strictBox*(1-options.Slack);
end