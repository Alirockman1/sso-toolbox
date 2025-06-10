function iChoice = component_trimming_cost(designSample,isViable,isExclude,isInsideComponent,isInsideAll,removalCandidate,ineligibleCandidate,varargin)
%COMPONENT_TRIMMING_COST Choose removal operation based on multiple criteria
%   COMPONENT_TRIMMING_COST determines which one of multiple candidate removal
%   operations yields the best outcome for a component, according to a 
%   prioritized list of selection criteria. Each criterion helps break ties 
%   among equally good candidates, in sequence.
%
%   ICHOICE = COMPONENT_TRIMMING_COST(DESIGNSAMPLE,ISVIABLE,ISEXCLUDE,
%   ISINSIDECOMPONENT,ISINSIDEALL,REMOVALCANDIDATE,INELIGIBLECANDIDATE) returns
%   the index ICHOICE corresponding to the best candidate removal operation,
%   considering the set of design points in DESIGNSAMPLE and various logical
%   masks:
%       - ISVIABLE : indicates viable points for the component.
%       - ISEXCLUDE : indicates points that are excluded.
%       - ISINSIDECOMPONENT : indicates points inside the component region.
%       - ISINSIDEALL : indicates points currently inside all components.
%   REMOVALCANDIDATE is an logical matrix whose columns represent different 
%   candidate operations. INELIGIBLECANDIDATE indicates which of these 
%   candidates are invalid or already excluded.
%
%   The function attempts to pick a single best candidate index by evaluating
%   the criteria listed in the 'SelectionCriteria' name-value pair (in order).
%   By default, 'SelectionCriteria' = 
%       {'NumberInsideViable','NumberInsideExclude','NumberExclude', ...
%        'NumberViable','NumberInsideNotExclude'}
%
%   ICHOICE = COMPONENT_TRIMMING_COST(...,'SelectionCriteria',CRITERIA) allows
%   specifying a custom sequence of criteria (CRITERIA can be a cell array of 
%   strings or function handles). For each criterion, the cost is computed for 
%   each candidate, and the function narrows down to whichever candidates 
%   minimize that cost. The process continues until only one candidate remains
%   or all criteria are exhausted.
%
%   Inputs:
%       - DESIGNSAMPLE : (nSample,nDesignVariable) double
%       - ISVIABLE : (nSample,1) logical
%       - ISEXCLUDE : (nSample,1) logical
%       - ISINSIDECOMPONENT : (nSample,1) logical
%       - ISINSIDEALL : (nSample,1) logical
%       - REMOVALCANDIDATE : (nSample,nRemovalCandidate) logical
%       - INELIGIBLECANDIDATE : (1,nRemovalCandidate) logical
%       - 'SelectionCriteria' : cell array of strings/function handles
%
%   Output:
%       - ICHOICE : integer
%
%   See also component_trimming_choice, component_trimming_operation.
%   
%   Copyright 2025 Eduardo Rodrigues Della Noce
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
    parser.addParameter('SelectionCriteria',{'NumberInsideViable','NumberInsideExclude','NumberInsideNotExclude','NumberExclude','NumberViable'});
    parser.parse(varargin{:});
    options = parser.Results;
    
    nCandidate = size(removalCandidate,2);
    nCriteria = length(options.SelectionCriteria);

    if(isempty(ineligibleCandidate) || all(ineligibleCandidate))
        isTie = true(1,nCandidate);
    else
        isTie = ~ineligibleCandidate;
    end
    nTie = sum(isTie);

    i = 1;
    while(nTie>1 && i<=nCriteria)
        criterionCurrent = options.SelectionCriteria{i};
        removalCandidateCurrent = removalCandidate(:,isTie);
        nCandidateCurrent = size(removalCandidateCurrent,2);

        removalCost = nan(1,nCandidateCurrent);
        if(isa(criterionCurrent,'function_handle'))
            removalCost = options.CostType(designSample,isViable,isExclude,isInsideComponent,isInsideAll,removalCandidateCurrent);
        elseif(strcmpi(criterionCurrent,'NumberInsideViable'))
            removalCost = -sum(isViable & isInsideAll & ~removalCandidateCurrent,1);
        elseif(strcmpi(criterionCurrent,'NumberViable'))
            removalCost = -sum(isViable & isInsideComponent & ~removalCandidateCurrent,1);
        elseif(strcmpi(criterionCurrent,'NumberInsideExclude'))
            removalCost = -sum(isExclude & isInsideAll & removalCandidateCurrent,1);
        elseif(strcmpi(criterionCurrent,'NumberExclude'))
            removalCost = -sum(isExclude & isInsideComponent & removalCandidateCurrent,1);
        elseif(strcmpi(criterionCurrent,'NumberInsideNotExclude'))
            removalCost = -sum(isInsideComponent & ~isExclude & ~removalCandidateCurrent,1);
        elseif(strcmpi(criterionCurrent,'VolumeInsideViable'))
            for j=1:nCandidateCurrent
                isInsideRemainViable = isViable & isInsideAll & ~removalCandidateCurrent(:,j);
                if(~any(isInsideRemainViable))
                    removalCost(j)=0;
                    continue;
                end
                
                boundingBox = design_bounding_box(designSample,isInsideRemainViable);
                isInBoundingBox = is_in_design_box(designSample,boundingBox);
                volumeFraction = sum(isInsideRemainViable)/sum(isInBoundingBox);
                removalCost(j) = -volumeFraction*prod(boundingBox(2,:)-boundingBox(1,:));
            end
        end

        isTie(isTie) = (removalCost==min(removalCost));
        nTie = sum(isTie);
        i = i+1;
    end
    iChoice = find(isTie,1);
end

