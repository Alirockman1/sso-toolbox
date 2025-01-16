function elementStrain = truss_deformed_strain(nodePosition,nodeDisplacement,nodeElement)
%TRUSS_DEFORMED_STRAIN compute (engineering) strain from truss displacements
%	TRUSS_DEFORMED_STRAIN computes the strain for each element of a truss given
%	the original/undeformed position of its nodes and the displacement as a 
%	result of a force.
%
%	ELEMENTSTRAIN = TRUSS_DEFORMED_STRAIN(NODEPOSITION,NODEDISPLACEMENT,
%	NODEELEMENT) receives the truss undeformed node positions NODEPOSITION and
%	their respective displacements when deformed NODEDISPLACEMENT, along the
%	node-element relation NODEELEMENT, and returns the strain for each element
%	in ELEMENTSTRAIN.
%
%   Input:
%		- NODEPOSITION : (nNode,nDimension) double
%		- NODEDISPLACEMENT : (nNode,nDimension) double
%		- NODEELEMENT : (nElement,2) double
%
%   Output:
%		- ELEMENTSTRAIN : (nElement,1) double
%
%   See also truss_analysis, plot_truss_element_response.
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

    nodeDistanceUndeformed = nodePosition(nodeElement(:,2),:) - nodePosition(nodeElement(:,1),:);
	elementLengthUndeformed = vecnorm(nodeDistanceUndeformed,2,2);

	nodeDistanceDeformed = nodeDistanceUndeformed + ...
		nodeDisplacement(nodeElement(:,2),:) - nodeDisplacement(nodeElement(:,1),:);
	elementLengthDeformed = vecnorm(nodeDistanceDeformed,2,2);

	elementStrain = (elementLengthDeformed-elementLengthUndeformed)./elementLengthUndeformed;
end