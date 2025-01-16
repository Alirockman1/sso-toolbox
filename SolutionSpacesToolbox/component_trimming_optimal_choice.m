function iChoice = component_trimming_optimal_choice(designSample,componentIndex,isViable,isExclude,isInsideComponent,varargin)
%COMPONENT_TRIMMING_OPTIMAL_CHOICE Select a single best trimming option
%   COMPONENT_TRIMMING_OPTIMAL_CHOICE analyzes multiple full trimming passes 
%   along with the design space samples and logical masks to choose which 
%   single trimming pass is optimal based on user-provided selection criteria.
%
%   ICHOICE = COMPONENT_TRIMMING_OPTIMAL_CHOICE(DESIGNSAMPLE,COMPONENTINDEX,
%   ISVIABLE,ISEXCLUDE,ISINSIDECOMPONENT) uses the design points in 
%   DESIGNSAMPLE, the logical arrays ISVIABLE and ISEXCLUDE indicating which 
%   points are viable or excluded, and ISINSIDECOMPONENT, a cell array whose 
%   elements each represent which design points are inside the components
%   for each pass. This function computes the best component trimming 
%   option by iterating through specified criteria to break ties. 
%   The chosen index is returned as ICHOICE.
%
%   ICHOICE = COMPONENT_TRIMMING_OPTIMAL_CHOICE(...,'SelectionCriteria',
%   CRITERIA) 
%   allows specifying CRITERIA, a cell array of strings or function handles 
%   that define how to evaluate and break ties. Each criterion is computed 
%   across all potential component definitions. By default, CRITERIA is:
%       {'NumberInsideViable','NumberInsideNotExclude'}
%   Additional recognized strings might include 'VolumeInsideViable'. 
%   A custom function handle can also be supplied, having the signature:
%       cost = f(designSample,componentIndex,isViable,isExclude,
%       isInsideComponentSubset)
%
%   Inputs:
%       - DESIGNSAMPLE : (nSample,nDesignVariable) double
%       - COMPONENTINDEX : Identifies the component to which points belong
%       - ISVIABLE : (nSample,1) logical
%       - ISEXCLUDE : (nSample,1) logical
%       - ISINSIDECOMPONENT : (1,nCandidate) cell
%       - 'SelectionCriteria' : cell array of strings or function handles
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