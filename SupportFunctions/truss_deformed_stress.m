function elementStress = truss_deformed_stress(elementAxialForce,elementCrossSectionArea)
%TRUSS_DEFORMED_STRESS compute (engineering) stress from truss displacements
%   TRUSS_DEFORMED_STRESS computes the stress for each element of a truss given
%   the resulting axial forces and the cross-section areas.
%
%   ELEMENTSTRESS = TRUSS_DEFORMED_STRESS(ELEMENTAXIALFORCE,
%   ELEMENTCROSSSECTIONAREA) receives the forces in local coordinates for each
%   element ELEMENTAXIALFORCE and the cross-section area for each element
%   ELEMENTCROSSSECTIONAREA, returning the stress on that truss element in 
%   ELEMENTSTRESS. Positive values indicate tensile stress, and negative values
%   compressive stress.
%
%   Input:
%       - ELEMENTAXIALFORCE : (nElement,2) double
%       - ELEMENTCROSSSECTIONAREA : double OR (nElement,1) double
%
%   Output:
%       - ELEMENTSTRESS : (nElement,1) double
%
%   See also truss_analysis, plot_truss_element_response.
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

    elementStress = -elementAxialForce(:,1)./elementCrossSectionArea;
end