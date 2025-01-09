function candidateBoxLean = box_trimming_leanness(designSample,isKeep,trimmingOrder,varargin)
%BOX_TRIMMING_LEANNESS Application of leanness condition for design boxes
%   BOX_TRIMMING_LEANNESS applies the leanness condition of requirement spaces 
%   to a given design box, ensuring that there are useful designs in all edges.
%   This works similarly to the usual trimming operation, but here, the designs
%   marked to be removed are only removed if no designs with the "keep" label
%   would be removed by said operation; otherwise, the trimming is not done.
%
%   CANDIDATEBOXLEAN = BOX_TRIMMING_LEANNESS(DESIGNSAMPLE,ISKEEP,
%   TRIMMINGORDER) receives the design sample points DESIGNSAMPLE, their label
%   in terms of whether it should be kept or not ISKEEP, and the order of
%   designs that may be removed TRIMMINGORDER, and returns the design box
%   after the leanness condition is applied in CANDIDATEBOXLEAN. For ISKEEP,
%   values of 'true' indicate that the respective sample must be kept, and 
%   'false' indicates the opposite. The assumed candidate box is a bounding box
%   around the designs that should be kept.
%
%   CANDIDATEBOXLEAN = BOX_TRIMMING_LEANNESS(DESIGNSAMPLE,ISKEEP,
%   TRIMMINGORDER,CANDIDATEBOX) allows direct specification of the initial 
%   candidate box CANDIDATEBOX.
%
%   ... = BOX_TRIMMING_LEANNESS(...,NAME,VALUE,...) allows the specification of
%   additional values. These include:
%       - 'Slack': a value between 0 and 1 which specifies where the boundary
%       should be relocated to when performing the trimming operation: for 0,
%       the boundary of the design box gets relocated to the closest acceptable
%       design; for 1, the boundary gets relocated to the unnaceptable design
%       being removed; and proportionally in-between for values between 0 and 1.
%
%   Inputs:
%       - DESIGNSAMPLE : (nSample,nDesignVariable) double
%       - ISKEEP : (nSample,1) logical
%       - TRIMMINGORDER : (nRemove,1) integer
%       - CANDIDATEBOX : (2,nDesignVariable) double
%           -- (1) : lower boundary of the design box
%           -- (2) : upper boundary of the design box
%       - 'Slack' : double
%
%   Outputs:
%       - CANDIDATEBOXLEAN : (2,nDesignVariable) double
%           -- (1) : lower boundary of the design box
%           -- (2) : upper boundary of the design box
%
%   See also design_bounding_box, box_trimming_operation, sso_box_stochastic.
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
    parser.addRequired('isKeep',@(x)islogical(x));
    parser.addRequired('trimmingOrder',@(x)(isnumeric(x)&&(size(x,2)==1))||isempty(x));
    parser.addOptional('CandidateBox',[],@(x)isnumeric(x)&&(size(x,1)==2 || isempty(x)));
    parser.addParameter('Slack',0.0,@(x)isnumeric(x)&&isscalar(x)&&(x>=0)&&(x<=1));
    parser.parse(designSample,isKeep,trimmingOrder,varargin{:});
    options = parser.Results;

    % assign conditional default values
    candidateBoxLean = options.CandidateBox;
    if(isempty(candidateBoxLean))
        [~,candidateBoxLean] = design_bounding_box(designSample,isKeep)
    end
    isInsideBox = is_in_design_box(designSample,candidateBoxLean);

    % if all designs inside box are useful, no need to apply leanness
    if(all(isKeep(isInsideBox)))
        return;
    else
        isInsideBoxKeep = isInsideBox & isKeep;
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

                closestKeep = [];
                if(options.Slack<1)
                    isInsideBoxCurrent = is_in_design_box(designSample,candidateBoxCurrent);
                    if(l==1) % lower boundary
                        remainRegion = (designSample(:,j)>=designSample(iExclude,j));
                        closestKeep = min(designSample(isInsideBoxCurrent&remainRegion&isKeep,j));
                    else % upper boundary
                        remainRegion = (designSample(:,j)<=designSample(iExclude,j));
                        closestKeep = max(designSample(isInsideBoxCurrent&remainRegion&isKeep,j));
                    end
                end

                if(isempty(closestKeep))
                    closestKeep = designSample(iExclude,k);
                end

                % move lower/upper bound
                candidateBoxCurrent(k,j) = designSample(iExclude,j)*options.Slack + closestKeep*(1-options.Slack);
                isInsideBoxCurrent = is_in_design_box(designSample,candidateBoxCurrent);
                isInsideBoxKeepCurrent = isInsideBoxCurrent & isKeep;

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
end