function iChoice = component_trimming_choice(designSample,componentIndex,isViable,isExclude,isInsideComponent,componentRemoval,ineligibleCandidate,varargin)
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
	parser.addParameter('SelectionCriteria',{'NumberInsideViable','NumberInsideExclude','NumberExclude','NumberViable','NumberInsideNotExclude'});
	parser.parse(varargin{:});
	options = parser.Results;

	nComponent = size(componentRemoval,2);
    removalCost = nan(1,nComponent);
    nCriteria = length(options.SelectionCriteria);
    isInsideAll = all(isInsideComponent,2);

    if(isempty(ineligibleCandidate) || all(ineligibleCandidate))
        isTie = true(1,nComponent);
    else
        isTie = ~ineligibleCandidate;
    end
    nTie = sum(isTie);

    i = 1;
    while(nTie>1 && i<=nCriteria)
        criterionCurrent = options.SelectionCriteria{i};
        removalCandidateCurrent = componentRemoval(:,isTie);
        componentIndexCurrent = componentIndex(isTie);
        isInsideComponentCurrent = isInsideComponent(:,isTie);
        componentRemovalCurrent = componentRemoval(:,isTie);

        if(isa(criterionCurrent,'function_handle'))
            removalCost = options.CostType(designSample,componentIndexCurrent,isViable,isExclude,isInsideComponentCurrent,componentRemovalCurrent);
        elseif(strcmpi(criterionCurrent,'NumberInsideViable'))
            removalCost = -sum(isViable & isInsideAll & ~componentRemovalCurrent,1);
        elseif(strcmpi(criterionCurrent,'NumberViable'))
            removalCost = -sum(isViable & isInsideComponentCurrent & ~componentRemovalCurrent,1);
        elseif(strcmpi(criterionCurrent,'NumberInsideExclude'))
            removalCost = -sum(isExclude & isInsideAll & componentRemovalCurrent,1);
        elseif(strcmpi(criterionCurrent,'NumberExclude'))
            removalCost = -sum(isExclude & isInsideComponentCurrent & componentRemovalCurrent,1);
        elseif(strcmpi(criterionCurrent,'NumberInsideNotExclude'))
            removalCost = -sum(isInsideComponentCurrent & ~isExclude & ~componentRemovalCurrent,1);
        end

        isTie(isTie) = (removalCost==min(removalCost));
        nTie = sum(isTie);
        i = i+1;
    end
    iChoice = find(isTie,1);
end