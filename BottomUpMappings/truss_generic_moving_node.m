function performanceMeasure = truss_generic_moving_node(designSample,systemParameter)
%TRUSS_GENERIC_MOVING_NODE Bottom-up Mapping (Node Positioning Problem)
%	TRUSS_GENERIC_MOVING_NODE calculates the relevant responses of a truss
%	with a force applied at its tip, where the positions of the nodes are the
%	design variables of the problem. 
%	o---0---0---0---0
%	   /|  /|  /|  /| \
%	  /	| /	| /	| /	|  o |
%	 /  |/  |/  |/  | /  V F
%	o---0---0---0---0
%	For the definition of the truss, the following information is necessary:
%		- The base position of the nodes (using 'nan' for moving nodes).
%		- The external force applied on each node.
%		- The connections of nodes which make up each element of the truss.
%		- The degrees of freedom which are fixed (either as logical, or with
%		numbers for fixed displacements and nan for free nodes).
%		- The cross-section area of each element.
%		- Young's Modulus of each element.
%		- Density of each element (optional, mass is considered -inf if not 
%		given).
%		- Moment of inertia of each element (optional, buckling ratio is 
%		considered -inf if not given).
%	As for what is computed from this:
%		- The displacement on the tip.
%		- The total mass of the structure.
%		- Absolute value of stresses in each element.
%		- Buckling ratio (current compression force / buckling factor limit) in
%		each element.
%
%	PERFORMANCEMEASURE = TRUSS_GENERIC_MOVING_NODE(DESIGNSAMPLE,
%	SYSTEMPARAMETER) receives the positions of the moving nodes in DESIGNSAMPLE
%	and the truss information in SYSTEMPARAMETER, returning the tip 
%	displacement, total mass of the structure, absolute stress on each element,
%	and buckling ratio for each element in PERFORMANCEMEASURE.
%
%	Input:
%		- DESIGNSAMPLE : (nSample,dimension*nMovingNode) double
%		- SYSTEMPARAMETER : struct
%			-- BaseNodePosition : (nNode,2) OR (nNode,3) double
%			-- NodeForce : (nNode,2) OR (nNode,3) double
%			-- NodeElement : (nElement,2) integer
%			-- FixedDegreesOfFreedom : (nNode,2) OR (nNode,3) logical OR double
%			-- ElementCrossSectionArea : double OR (nElement,1) double
%			-- ElementYoungsModulus : double OR (nElement,1) double
%			-- ElementDensity : double OR (nElement,1) double
%			-- ElementMomentOfInertia : double OR (nElement,1) double
%
%	Output:
%		- PERFORMANCEMEASURE : (nSample,2+2*nElement) double
%			-- (1) : displacement of tip in the direction of the applied force
%			-- (2) : total mass of the structure
%			-- (2+i) : absolute value of stress in element i
%			-- (2+nElement+i) : buckling ratio (F/F_lim) in element i
%
%	See also truss_analysis, truss_two_bar_hollow_circle.
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

	% check which entries are design variables based on nan entry
	isDesignVariable = isnan(systemParameter.BaseNodePosition);
	isTrussTip = (systemParameter.NodeForce~=0);
	nDimension = size(systemParameter.BaseNodePosition,2);
	nElement = size(systemParameter.NodeElement,1);

	nodePosition = systemParameter.BaseNodePosition;
	nSample = size(designSample,1);
	performanceMeasure = nan(nSample,2+2*nElement);
	for i=1:nSample
		% assign design variables to build current truss
		nodePosition(isDesignVariable) = column_vector_to_row_major_matrix(designSample(i,:)',nDimension);

		% solve truss system --> find displacements and forces
		[nodeDisplacement,~,elementAxialForce] = ...
			truss_analysis(...
				nodePosition,...
				systemParameter.FixedDegreesOfFreedom,...
				systemParameter.NodeForce,...
				systemParameter.NodeElement,...
				systemParameter.ElementCrossSectionArea,...
				systemParameter.ElementYoungsModulus);

		% find stresses
		elementStress = truss_deformed_stress(elementAxialForce,systemParameter.ElementCrossSectionArea);

		% compute element lengths for later calculations
		nodeDistance = nodePosition(systemParameter.NodeElement(:,2),:)-nodePosition(systemParameter.NodeElement(:,1),:);
		elementLength = vecnorm(nodeDistance,2,2);

		% find mass
		if(isfield(systemParameter,'ElementDensity'))
			totalMass = sum(systemParameter.ElementDensity.*elementLength.*systemParameter.ElementCrossSectionArea);
		else
			totalMass = -inf;
		end

		% find buckling ratio
		if(isfield(systemParameter,'ElementMomentOfInertia'))
			bucklingCriticalLoad = pi^2.*systemParameter.ElementYoungsModulus.*systemParameter.ElementMomentOfInertia...
				./(elementLength.^2);
	        bucklingRatio = elementAxialForce(:,1)./bucklingCriticalLoad;
	    else
	    	bucklingRatio = -inf(nElement,1);
	    end

        % deal with numerical errors
        displacement = -nodeDisplacement(isTrussTip);
        displacement(isnan(displacement)) = +inf;
        elementStress(isnan(elementStress)) = +inf;
        bucklingRatio(isnan(bucklingRatio)) = +inf;

		performanceMeasure(i,:) = [displacement,totalMass,abs(elementStress'),bucklingRatio'];
	end
end

