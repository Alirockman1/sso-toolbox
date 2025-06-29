function [nodeDisplacement,nodeReactionForce,elementAxialForce] = truss_analysis(nodePosition,nodeFixedDegreesOfFreedom,nodeForce,nodeElement,elementCrossSectionArea,elementYoungsModulus)
%TRUSS_ANALYSIS Solve truss problem using the direct stiffness method
%	TRUSS_ANALYSIS uses the direct stiffness method to find the displacements,
%	reaction forces and element axial forces of a given truss system.
%
%	NODEDISPLACEMENT = TRUSS_ANALYSIS(NODEPOSITION,NODEFIXEDDEGREESOFFREEDOM,
%	NODEFORCE,NODEELEMENT,ELEMENTCROSSSECTIONAREA,ELEMENTYOUNGSMODULUS) 
%	receives the node positions of the truss in NODEPOSITION, the fixed degrees
%	of freedom of those nodes in NODEFIXEDDEGREESOFFREEDOM, the external forces
%	applied to the nodes in NODEFORCE, the nodes that form each element in 
%	NODEELEMENT, the cross-section areas of the elements ELEMENTCROSSSECTIONAREA 
%	and the elasticity module of the elements ELEMENTYOUNGSMODULUS, returning
%	the displacement of each node NODEDISPLACEMENT after solving the system.
%	Of note regarding these inputs:
%		- NODEPOSITION is given as an array where each column is one dimension
%		and each row is a new node.
%		- NODEFIXEDDEGREESOFFREEDOM can be specified in one of two forms: it
%		can be given either as a logical array with the same size as 
%		NODEPOSITION, where 'true' indicates that degree of freedom is fixed
%		and the prescribed displacement is zero, and 'false' indicates it is
%		free; or as an array of numerical values, where those values indicate
%		the prescribed displacement, and nan entries indicate that degree of
%		freedom is actuallz free.
%		- NODEELEMENT is a two-column array, where each column has one node
%		ID entry (row number as specified in NODEPOSITION) and each row is one
%		element.
%		- ELEMENTCROSSSECTIONAREA and ELEMENTYOUNGSMODULUS can each be given
%		as either a single value (which is assumed to be the same for all 
%		elements) or as individual values for each element as column arrays.
%
%	[NODEDISPLACEMENT,NODEREACTIONFORCE] = TRUSS_ANALYSIS(...) also returns
%	the node reaction force for the fixed degrees of freedom NODEREACTIONFORCE.
%	The summation of all NODEREACTIONFORCE and NODEFORCE for each dimension 
%	should be zero.
%
%	[NODEDISPLACEMENT,NODEREACTIONFORCE,ELEMENTAXIALFORCE] = TRUSS_ANALYSIS(...)
%	also returns the forces applied at each element ELEMENTAXIALFORCE. These 
%	forces are specified in a two-columns array, one for each node that composes
%	the element. With the chosen local coordinate system, positive values on the
%	first column indicating compression force, and negative values indicating
%	tension force. The value on the second column should be the same but with 
%	an opposite sign, as locally the truss should be in equilibrium.
%
%	References:
%	- https://people.duke.edu/~hpgavin/cee421/truss-3d.pdf
%	- https://learnaboutstructures.com/Stiffness-Method-for-One-Dimensional-Truss-Elements 
%
%	Input:
%		- NODEPOSITION : (nNode,nDimension) double
%		- NODEFIXEDDEGREESOFFREEDOM : (nNode,nDimension) logical OR 
%		(nNode,nDimension) double
%		- NODEFORCE : (nNode,nDimension) double
%		- NODEELEMENT : (nElement,2) double
%		- ELEMENTCROSSSECTIONAREA : double OR (nElement,1) double
%		- ELEMENTYOUNGSMODULUS : double OR (nElement,1) double
%	
%	Output: 
%		- NODEDISPLACEMENT : (nNode,nDimension) double
%		- NODEREACTIONFORCE : (nNode,nDimension) double
%		- ELEMENTAXIALFORCE : (nElement,2) double
%
%	See also plot_truss_deformation, plot_truss_element_response.
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

	% get important problem dimensions
	nNode = size(nodePosition,1);
	nDimension = size(nodePosition,2);
	nElement = size(nodeElement,1);

	% use same area and Yong's Modulus for every element if only one value was given
	if(length(elementCrossSectionArea)==1)
		elementCrossSectionArea = elementCrossSectionArea * ones(nElement,1);
	end
	if(length(elementYoungsModulus)==1)
		elementYoungsModulus = elementYoungsModulus * ones(nElement,1);
    end

    % check for compatibility of inputs
    if(any(size(nodePosition)~=size(nodeFixedDegreesOfFreedom)))
        error('TrussAnalysis:WrongSizeInput:DOF','Size of array with fixed degrees of freedom incompatible with defined nodes.');
    end
    if(any(size(nodePosition)~=size(nodeForce)))
        error('TrussAnalysis:WrongSizeInput:Force','Size of array with external forces incompatible with defined nodes.');
    end
    if(size(nodeElement,2)~=2)
        error('TrussAnalysis:WrongSizeInput:Element','Elements can only connect 2 nodes to each other, check the definition.');
    end
    if(size(nodeElement,1)~=size(elementCrossSectionArea,1))
        error('TrussAnalysis:WrongSizeInput:Area','Size of array with cross-section areas incompatible with defined elements.');
    end
    if(size(nodeElement,1)~=size(elementYoungsModulus,1))
        error('TrussAnalysis:WrongSizeInput:Elasticity','Size of array with Young''s Modulus incompatible with defined elements.');
    end

    % convert input into separate information: whether it is fixed or not, and
    % the prescribed displacement
    if(all(islogical(nodeFixedDegreesOfFreedom)))
    	isFixed = nodeFixedDegreesOfFreedom;

    	prescribedDisplacement = nan(size(nodeFixedDegreesOfFreedom));
    	prescribedDisplacement(isFixed) = 0;
    end
    if(all(isnumeric(nodeFixedDegreesOfFreedom)))
    	isFixed = ~isnan(nodeFixedDegreesOfFreedom);

    	prescribedDisplacement = nodeFixedDegreesOfFreedom;
    end

    % find node cartesian distances for each element
	nodeDistance = nodePosition(nodeElement(:,2),:)-nodePosition(nodeElement(:,1),:);

	% start by building global stiffness matrix
	nSystemDegreesOfFreedom = nNode * nDimension;
	globalStiffness = zeros(nSystemDegreesOfFreedom,nSystemDegreesOfFreedom);

	% build local stiffnesses and add them to global matrix
	for i=1:nElement
		[elementLocalStiffness,localToGlobalSpatialTransformation] = ...
			truss_local_stiffness_global_transformation(...
				nodeDistance(i,:),...
				elementYoungsModulus(i),...
				elementCrossSectionArea(i));

		% get the element's contribution to global stiffness
		% Kg = T^T * kl * T
		elementGlobalStiffness = ...
			localToGlobalSpatialTransformation'*...
			elementLocalStiffness*...
			localToGlobalSpatialTransformation;

		% add that contribution to the global stiffness matrix
		iGlobalNode1 = [0:(nDimension-1)] + nDimension*(nodeElement(i,1)-1) + 1;
		iGlobalNode2 = [0:(nDimension-1)] + nDimension*(nodeElement(i,2)-1) + 1;
		iGlobalElement = [iGlobalNode1,iGlobalNode2];
		globalStiffness(iGlobalElement,iGlobalElement) = ...
			globalStiffness(iGlobalElement,iGlobalElement) + ...
			elementGlobalStiffness;
	end

	% get force and DoF in equivalent column-form
	globalForce = row_major_matrix_to_column_vector(nodeForce);
	globalIsFixed = row_major_matrix_to_column_vector(isFixed);
	globalPrescribedDisplacement = row_major_matrix_to_column_vector(prescribedDisplacement);

	% find displacements: u_free = K_free^-1 * (F_free - K_fixed*u_fixed)
	nodeDisplacement = globalPrescribedDisplacement;
	nodeDisplacement(~globalIsFixed) = globalStiffness(~globalIsFixed,~globalIsFixed) \ ...
		(globalForce(~globalIsFixed) - globalStiffness(~globalIsFixed,globalIsFixed)*globalPrescribedDisplacement(globalIsFixed));

	% find reaction forces: F = K*u
	nodeReactionForce = zeros(nSystemDegreesOfFreedom,1);
	nodeReactionForce(globalIsFixed) = globalStiffness(globalIsFixed,:)*nodeDisplacement;

	% convert back to normal notation where each row is a node and each column is a spatial dimension
	nodeDisplacement = column_vector_to_row_major_matrix(nodeDisplacement,nDimension);
	nodeReactionForce = column_vector_to_row_major_matrix(nodeReactionForce,nDimension);

	% find axial forces for each element
	elementAxialForce = nan(nElement,2);
	for i=1:nElement
		[elementLocalStiffness,localToGlobalSpatialTransformation] = ...
			truss_local_stiffness_global_transformation(...
				nodeDistance(i,:),...
				elementYoungsModulus(i),...
				elementCrossSectionArea(i));

		% local displacement
		globalNodeDisplacement = nodeDisplacement([nodeElement(i,1);nodeElement(i,2)],:);
		localNodeDisplacement = localToGlobalSpatialTransformation*row_major_matrix_to_column_vector(globalNodeDisplacement);

        % local force based on displacement
		elementAxialForce(i,:) = (elementLocalStiffness * localNodeDisplacement)';
	end
end

function [elementLocalStiffness,localToGlobalSpatialTransformation] = truss_local_stiffness_global_transformation(nodeDistance,elementYoungsModulus,elementCrossSectionArea)
	nDimension = size(nodeDistance,2);
	elementLength = vecnorm(nodeDistance,2,2);

	% truss element stiffness in local coordinates
	elementStiffnessCoefficient = elementYoungsModulus*elementCrossSectionArea/elementLength;
	elementLocalStiffness = elementStiffnessCoefficient * [1 -1;-1 1];

	% transformation matrix from local to global coordinates 
	localToGlobalSpatialTransformation = [...
		nodeDistance./elementLength , zeros(1,nDimension);...
		zeros(1,nDimension) , nodeDistance./elementLength];
end