function iChoice = component_trimming_choice(cost,designSample,componentIndex,isKeep,isRemainComponent,componentRemoval,varargin)
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
	parser.addParameter('WeightedCostType','SimpleCost');
	parser.parse(varargin{:});
	options = parser.Results;
	
	weightedCost = nan(size(cost));
	if(isa(options.WeightedCostType,'function_handle'))
		weightedCost = options.WeightedCostType(cost,componentIndex);
	elseif(strcmpi(options.WeightedCostType,'ComponentDimension'))
		for i=1:size(componentIndex,2)
			weightedCost(i) = cost(i)*2^length(componentIndex{i});
		end
	else%if(strcmpi(options.WeightedCostType,'SimpleCost'))
		weightedCost = cost;
	end

	% select index with minimum weighted cost
	[sortedCost,iCost] = sort(weightedCost);
    iTieBreaker = (sortedCost==sortedCost(1));
    iTieBreaker = iTieBreaker(iCost);

    nTieBreaker = sum(iTieBreaker);
    if(nTieBreaker>1)
        % tie-breaker: total amount of points eliminated
        %removalCostTieBreaker = sum(componentRemoval(:,iTieBreaker),1);
        
        % tie-breaker: volume that remains
        nComponent = length(componentIndex);
        componentMeasureBase = nan(1,nComponent);
        componentMeasureTrimmed = nan(1,nComponent);
        for i=1:nComponent
        	designSampleComponent = designSample(:,componentIndex{i});

        	% compute base measure -> no removal
        	boundingBox = design_bounding_box(designSampleComponent,isRemainComponent(:,i));
        	pointInsideBox = is_in_design_box(designSampleComponent,boundingBox);
        	volumeFraction = sum(pointInsideBox & isRemainComponent(:,i))/sum(pointInsideBox);
        	componentMeasureBase(i) = volumeFraction*prod(boundingBox(2,:) - boundingBox(1,:));

        	% compute trimmed measure -> with removal
        	boundingBox = design_bounding_box(designSampleComponent,isRemainComponent(:,i) & ~componentRemoval(:,i));
        	pointInsideBox = is_in_design_box(designSampleComponent,boundingBox);
        	volumeFraction = sum(pointInsideBox & isRemainComponent(:,i) & ~componentRemoval(:,i))/sum(pointInsideBox);
        	componentMeasureTrimmed(i) = volumeFraction*prod(boundingBox(2,:) - boundingBox(1,:));
        end
        removalCostTieBreaker = nan(1,nTieBreaker);
        for i=1:nTieBreaker
            iCurrent = convert_index_base(iTieBreaker',i,'backward');
            removalCostTieBreaker(i) = -prod(componentMeasureBase)*componentMeasureTrimmed(iCurrent)/componentMeasureBase(iCurrent);
        end

        [~,iCostTieBreaker] = sort(removalCostTieBreaker);
        iChoice = convert_index_base(iTieBreaker',iCostTieBreaker(1),'backward');
    else
        iChoice = iCost(1);
    end
end