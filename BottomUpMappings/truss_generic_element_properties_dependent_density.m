function [performanceMeasure, physicalFeasibilityMeasure] = truss_generic_element_properties_dependent_density(designSample, systemParameter)
%TRUSS_GENERIC_ELEMENT_PROPERTIES_DEPENDENT_DENSITY Bottom-up Mapping (Density from Young's Modulus)
%   This variant calculates material density based on Young's modulus using a
%   provided estimation function. Design variables per element are:
%   - Radius (r)
%   - Thickness (t)
%   - Young's modulus (E)
%   (Density ρ is calculated from E via systemParameter.EstimateMassGivenYoungsModulus)
%
%   PERFORMANCEMEASURE = TRUSS_GENERIC_ELEMENT_PROPERTIES_DEPENDENT_DENSITY(DESIGNSAMPLE,
%   SYSTEMPARAMETER) takes:
%       - DESIGNSAMPLE : (nSample,3*nElement) matrix with [E, r, t] per element
%       - SYSTEMPARAMETER : struct containing:
%           -- NodePosition: (nNode, 2/3) node coordinates
%           -- NodeForce: (nNode, 2/3) external forces
%           -- NodeElement: (nElement, 2) element connectivity
%           -- FixedDegreesOfFreedom: (nNode, 2/3) boundary conditions
%           -- EstimateMassGivenYoungsModulus: Function handle for ρ = f(E)
%
%   Outputs match truss_generic_element_properties format with calculated density.
%
%	See also truss_generic_element_properties, truss_analysis, truss_two_bar_hollow_circle.
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

    nElements = size(designSample, 2) / 3;
    
    % Extract Young's modulus from every 3rd column starting at column 1
    youngsModulus = designSample(:, 3:3:3*nElements);
    
    % Calculate density using provided estimation function
    elementDensity = systemParameter.EstimateMassGivenYoungsModulus(youngsModulus);
    yieldStrength = systemParameter.EstimateYieldStrengthGivenYoungsModulus(youngsModulus);
    
    % Update system parameters with calculated density
    systemParameter.ElementDensity = elementDensity;
    systemParameter.ElementYieldStrength = yieldStrength;
    
    % Call fixed density version with updated parameters
    [performanceMeasure, physicalFeasibilityMeasure] = truss_generic_element_properties_fixed_density(designSample, systemParameter);
end