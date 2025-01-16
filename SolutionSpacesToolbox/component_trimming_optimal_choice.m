function iChoice = component_trimming_optimal_choice(designSample,componentIndex,isViable,isExclude,isInsideComponent,varargin)
%COMPONENT_TRIMMING_CHOICE Choice of which component trimming to perform
%	COMPONENT_TRIMMING_CHOICE compares the costs of the candidate trimming 
%	operation to be performed in every component and chooses the one with the
%	minimum weighted cost. 
%
%	ICHOICE = COMPONENT_TRIMMING_CHOICE(COST,COMPONENT) receives the cost of
%	trimming for each component in COST and the component definition itself in
%	COMPONENT, and return the index ICHOICE that has the minimum cost without
%	any specific weighting.
%
%	ICHOICE = COMPONENT_TRIMMING_CHOICE(...NAME,VALUE,...) allows the 
%	specification of name-value pair arguments. These can be:
%		- 'WeightedCostType' : how to weight the cost for each component. This 
%		can be one of the following:
%			-- 'SimpleCost' : no special weighting, costs are used as-is.
%			-- 'ComponentDimension' : costs are multiplied by a factor of 
%			2^(component dimension).
%		Alternatively, a function handle can be used to define custom weighting.
%		This function must have the form 'weightedCost = f(cost,component)'.
%
%   Input:
%		- COST : (1,nComponent) double
%		- COMPONENT : (1,nComponent) cell
%		- 'WeightedCostType' : char OR string OR function_handle
%
%   Output:
%		- ICHOICE : integer
%
%   See also component_trimming_cost, component_trimming_operation.
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
	parser.addParameter('SelectionCriteria',{'NumberInsideViable','NumberInsideNotExclude'});
	parser.parse(varargin{:});
	options = parser.Results;

    nCriteria = length(options.SelectionCriteria);
    nCandidate = length(isInsideComponent);

    nSample = size(designSample,1);
    isInsideAll = false(nSample,nCandidate);
    for i=1:nCandidate
        isInsideComponentCurrent = isInsideComponent{i};
        if(~isempty(isInsideComponentCurrent))
            isInsideAll(:,i) = all(isInsideComponentCurrent,2);
        end
    end

    isTie = any(isInsideAll,1);
    nTie = sum(isTie);

    i = 1;
    while(nTie>1 && i<=nCriteria)
        criterionCurrent = options.SelectionCriteria{i};
        isInsideComponentCurrent = isInsideComponent(isTie);
        isInsideAllCurrent = isInsideAll(:,isTie);
        nCandidateCurrent = size(isInsideComponentCurrent,2);

        removalCost = nan(1,nCandidateCurrent);
        if(isa(criterionCurrent,'function_handle'))
            removalCost = options.CostType(designSample,componentIndex,isViable,isExclude,isInsideComponentCurrent);
        elseif(strcmpi(criterionCurrent,'NumberInsideViable'))
            removalCost = -sum(isViable & isInsideAllCurrent,1);
        elseif(strcmpi(criterionCurrent,'NumberInsideNotExclude'))
            removalCost = -sum(isInsideAllCurrent & ~isExclude,1);
        elseif(strcmpi(criterionCurrent,'VolumeInsideViable'))
            for j=1:nCandidate
                isViableRemain = isViable & isInsideAllCurrent(:,j);
                if(~any(isViableRemain))
                    removalCost(j) = 0;
                    continue;
                end
                
                boundingBox = design_bounding_box(designSample,isViableRemain);
                isInBoundingBox = is_in_design_box(designSample,boundingBox);
                volumeFraction = sum(isViableRemain)/sum(isInBoundingBox);
                removalCost(j) = -volumeFraction*prod(boundingBox(2,:)-boundingBox(1,:));
            end
        end

        isTie(isTie) = (removalCost==min(removalCost));
        nTie = sum(isTie);
        i = i+1;
    end
    iChoice = find(isTie,1);
end