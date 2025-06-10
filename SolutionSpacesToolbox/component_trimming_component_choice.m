function iChoice = component_trimming_component_choice(designSample,componentIndex,isViable,isExclude,isInsideComponent,componentRemoval,ineligibleCandidate,varargin)
%COMPONENT_TRIMMING_COMPONENT_CHOICE Select best component-trimming operation
%   COMPONENT_TRIMMING_COMPONENT_CHOICE evaluates multiple candidate trimming
%   operations across different components, using a set of selection criteria to 
%   break ties. The function returns an index identifying which candidate yields 
%   the best outcome according to the listed criteria.
%
%   ICHOICE = COMPONENT_TRIMMING_COMPONENT_CHOICE(DESIGNSAMPLE,COMPONENTINDEX,
%   ISVIABLE,ISEXCLUDE,ISINSIDECOMPONENT,COMPONENTREMOVAL,INELIGIBLECANDIDATE)
%   processes the design-space samples DESIGNSAMPLE and the logical arrays:
%       - ISVIABLE : which points are viable.
%       - ISEXCLUDE : which points are excluded.
%       - ISINSIDECOMPONENT : which points lie inside each component.
%   COMPONENTREMOVAL is an logical matrix, whose columns correspond to different 
%   trimming operations. INELIGIBLECANDIDATE indicates which operations are 
%   invalid. The function identifies a single best candidate among all columns 
%   using tie-breaking logic.
%
%   ICHOICE = COMPONENT_TRIMMING_COMPONENT_CHOICE(...,'SelectionCriteria',
%   CRITERIA) allows the user to specify CRITERIA, a cell array listing how to 
%   compare and break ties among candidates. Each criterion is one of the 
%   recognized strings (e.g., 'NumberInsideViable', 'NumberExclude', etc.) or a 
%   function handle with the form:
%       cost = f(designSample,componentIndex,isViable,isExclude,
%       isInsideComponent,componentRemoval)
%   By default, CRITERIA is:
%       {'NumberInsideViable','NumberInsideExclude','NumberExclude', ...
%        'NumberViable','NumberInsideNotExclude'}.
%
%   Inputs:
%       - DESIGNSAMPLE : (nSample,nDesignVariable) double
%       - COMPONENTINDEX : (1,nComponent) or (nSample,1) identifying components
%       - ISVIABLE : (nSample,1) logical
%       - ISEXCLUDE : (nSample,1) logical
%       - ISINSIDECOMPONENT : (nSample,nComponent) logical
%       - COMPONENTREMOVAL : (nSample,nCandidate) logical
%       - INELIGIBLECANDIDATE : (1,nCandidate) logical
%       - 'SelectionCriteria' : cell array of strings/function handles
%
%   Output:
%       - ICHOICE : integer
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
	parser.addParameter('SelectionCriteria',{'NumberInsideViable','NumberInsideExclude','NumberInsideNotExclude','NumberExclude','NumberViable'});
	parser.parse(varargin{:});
	options = parser.Results;

	nComponent = size(componentRemoval,2);
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
        componentIndexCurrent = componentIndex(isTie);
        isInsideComponentCurrent = isInsideComponent(:,isTie);
        componentRemovalCurrent = componentRemoval(:,isTie);
        nCandidateCurrent = size(componentRemovalCurrent,2);

        removalCost = nan(1,nCandidateCurrent);
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
        elseif(strcmpi(criterionCurrent,'VolumeInsideViable'))
            for j=1:nCandidateCurrent
                isViableRemain = isViable & isInsideAll & ~componentRemovalCurrent(:,j);
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