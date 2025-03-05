function [performanceMeasure, physicalFeasibilityMeasure] = truss_generic_element_properties_fixed_density(designSample, systemParameter)
%TRUSS_GENERIC_ELEMENT_PROPERTIES_FIXED_DENSITY Bottom-up Mapping (Fixed Density)
%   Variant with fixed material density from system parameters. Design variables
%   per element are:
%   - Radius (r)
%   - Thickness (t)
%   - Young's modulus (E)
%   (Density ρ comes from systemParameter.ElementDensity)
%
%   PERFORMANCEMEASURE = TRUSS_GENERIC_ELEMENT_PROPERTIES_FIXED_DENSITY(DESIGNSAMPLE,
%   SYSTEMPARAMETER) takes:
%       - DESIGNSAMPLE : (nSample,3*nElement) matrix with [E, r, t] per element
%       - SYSTEMPARAMETER : struct containing:
%           -- NodePosition: (nNode, 2/3) node coordinates
%           -- NodeForce: (nNode, 2/3) external forces
%           -- NodeElement: (nElement, 2) element connectivity
%           -- FixedDegreesOfFreedom: (nNode, 2/3) boundary conditions
%           -- ElementDensity: scalar or (nElement,1) density values
%
%   Outputs match truss_generic_element_properties format using provided density.
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

    % get important sizes
    nElements = size(designSample,2)/3;
    nSample = size(designSample,1);
    
    % get element density
    elementDensity = systemParameter.ElementDensity;
    if(isscalar(elementDensity))
        elementDensity = elementDensity * ones(1, nSample*nElements);
    else
        elementDensity = row_major_matrix_to_column_vector(elementDensity);
    end

    % get element yield strength
    elementYieldStrength = systemParameter.ElementYieldStrength;
    if(isscalar(elementYieldStrength))
        elementYieldStrength = elementYieldStrength * ones(1, nSample*nElements);
    else
        elementYieldStrength = row_major_matrix_to_column_vector(elementYieldStrength);
    end

    % Reshape the input array to separate elements and insert density
    reshapedSamplePerElement = reshape(designSample',3,[])';
    
    % Create density array matching dimensions and concatenate
    expandedSample = [reshapedSamplePerElement,elementDensity,elementYieldStrength];
    
    % Reshape back to original format with added density column
    designSampleWithDensity = reshape(expandedSample', 5*nElements, [])';

    % Call the original function with the expanded design sample
    [performanceMeasure, physicalFeasibilityMeasure] = truss_generic_element_properties(designSampleWithDensity, systemParameter);
end 