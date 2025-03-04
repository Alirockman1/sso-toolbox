function [performanceMeasure, physicalFeasibilityMeasure] = truss_generic_element_properties(designSample, systemParameter)
%TRUSS_GENERIC_ELEMENT_PROPERTIES Bottom-up Mapping (Element Properties)
%   This function calculates truss performance with element properties as
%   design variables. Each element's properties are defined by:
%   - Young's modulus (E)
%   - Radius (r)
%   - Thickness (t)
%   - Density (ρ)
%
%   Input:
%       - designSample: (nSample, 5*nElement) matrix where each group of 5 columns
%         represents [E, r, t, ρ, sigmaY] for an element
%       - systemParameter: struct containing:
%           -- NodePosition: (nNode, 2/3) node coordinates
%           -- NodeForce: (nNode, 2/3) external forces
%           -- NodeElement: (nElement, 2) element connectivity
%           -- FixedDegreesOfFreedom: (nNode, 2/3) boundary conditions
%
%   Output:
%       - performanceMeasure: (nSample, 2+2*nElement) matrix containing:
%           (1) Tip displacement, (2) Total mass, 
%           (3:2+nElement) Ratio of absolute stresses and yield strength, 
%           (3+nElement:end) Buckling ratios
%
%	See also truss_analysis, truss_two_bar_hollow_circle.
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

    nElement = size(systemParameter.NodeElement, 1);
    nSample = size(designSample, 1);
    performanceMeasure = nan(nSample, 2 + 2*nElement);
    
    % Precompute element lengths based on fixed node positions
    nodePosition = systemParameter.NodePosition;
    elementNodes = systemParameter.NodeElement;
    nodeDistance = nodePosition(elementNodes(:,2),:) - nodePosition(elementNodes(:,1),:);
    elementLength = vecnorm(nodeDistance, 2, 2);

    % physical feasibility 
    physicalFeasibilityMeasure = nan(nSample,nElement);
    includeMaterialPhysicalFeasibility = isfield(systemParameter, 'PhysicalFeasibilityEstimator') && ~isempty(systemParameter.PhysicalFeasibilityEstimator);
    if(includeMaterialPhysicalFeasibility)
        physicalFeasibilityMeasure = [physicalFeasibilityMeasure,nan(nSample,nElement)];
    end

    % Find tip node with applied force
    isTrussTip = (systemParameter.NodeForce~=0);
    for i = 1:nSample
        % Extract element properties from design sample
        elementProperties = reshape(designSample(i,:), 5, nElement)';
        youngsModulus = elementProperties(:,1);
        elementRadius = elementProperties(:,2);
        elementThickness = elementProperties(:,3);
        elementDensity = elementProperties(:,4);
        elementYieldStrength = elementProperties(:,5);
        
        % Calculate cross-sectional properties
        elementArea = pi*(elementRadius.^2 - (elementRadius - elementThickness).^2);  % Area of hollow circle
        elementMomentOfInertia = (pi/4)*(elementRadius.^4 - (elementRadius - elementThickness).^4);  % Moment of inertia
        
        % Perform truss analysis
        [nodeDisplacement, ~, axialForce] = truss_analysis(...
            nodePosition,...
            systemParameter.FixedDegreesOfFreedom,...
            systemParameter.NodeForce,...
            elementNodes,...
            elementArea,...
            youngsModulus);
        
        % Calculate performance metrics
        stress = truss_deformed_stress(axialForce,elementArea);
        totalMass = sum(elementDensity .* elementArea .* elementLength);
        
        % Buckling calculations
        bucklingCriticalLoad = (pi^2 * youngsModulus .* elementMomentOfInertia) ./ (elementLength.^2);

        stressRatio = abs(stress) ./ elementYieldStrength;
        bucklingRatio = axialForce(:,1) ./ bucklingCriticalLoad;
        tipDisplacement = -nodeDisplacement(isTrussTip);
        
        % Handle numerical issues
        tipDisplacement(isnan(tipDisplacement)) = inf;
        stressRatio(isnan(stressRatio)) = inf;
        bucklingRatio(isnan(bucklingRatio)) = inf;
        
        performanceMeasure(i,:) = [tipDisplacement, totalMass, stressRatio', bucklingRatio'];

        physicalFeasibilityMeasure(i,1:nElement) = (elementThickness - elementRadius)';
        if(includeMaterialPhysicalFeasibility)
            physicalFeasibilityMeasure(i,(nElement+1):end) = systemParameter.PhysicalFeasibilityEstimator(youngsModulus,elementDensity);
        end
    end
end 