function iChoice = component_trimming_choice(cost,component,varargin)
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
		weightedCost = options.WeightedCostType(cost,component);
	elseif(strcmpi(options.WeightedCostType,'ComponentDimension'))
		for i=1:size(component,2)
			weightedCost(i) = cost(i)*2^length(component{i});
		end
	else%if(strcmpi(options.WeightedCostType,'SimpleCost'))
		weightedCost = cost;
	end

	% select index with minimum weighted cost
	[~,iChoice] = min(weightedCost);
end