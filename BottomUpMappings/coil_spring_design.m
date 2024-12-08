function performanceMeasure = coil_spring_design(designSample,systemParameter)
%COIL_SPRING_DESIGN Bottom-up Mapping (Coil Spring Design Problem)
%	COIL_SPRING_DESIGN analyzes the most important characteristics of a coil
%	spring (stiffness, mass, operation space, and limiting deformations) for the
%	given designs.
%
%	PERFORMANCEMEASURE = COIL_SPRING_DESIGN(DESIGNSAMPLE,SYSTEMPARAMETER)  
%   receives the diameters, number of windings, spacing of wires, and material 
%	properties of the spring in DESIGNSAMPLE, and returns the coil spring 
%	characteristics in PERFORMANCEMEASURE.
%
%	Input:
%       - DESIGNSAMPLE : (nSample,7) double
%           -- (1) : wire diameter (d)
%           -- (2) : coil diameter (D)
%           -- (3) : number of windings (n)
%           -- (4) : spacing of wires in relaxed state (s)
%           -- (5) : material density (rho)
%           -- (6) : material shear modulus (G)
%           -- (7) : material yield stres (sigmaY)
%
%   Output:
%       - PERFORMANCEMEASURE : (nSample,8) double
%           -- (1) : stiffness (K)
%           -- (2) : mass (m)
%           -- (3) : height of operation space (H)
%           -- (4) : outer diameter of operation space (Wo)
%           -- (5) : inner diameter of operation space (Wi)
%           -- (6) : deformation to compaction (uc)
%           -- (7) : deformation to plastic yield (uY)
%           -- (8) : deformation to end of linearity (uL)
%
%   See also coil_spring_design_given_material.
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

	%% unwrap design variables
	wireDiameter = designSample(:,1); % d
	coilDiameter = designSample(:,2); % D
	nWindings = designSample(:,3); % n
	spacingWiresRelaxed = designSample(:,4); % s
	density = designSample(:,5); % rho
	shearModulus = designSample(:,6); % G
	yieldStress = designSample(:,7); % sigmaY

	%% bottom-up mappings
	% k + G*d^4/(8*D^3*n)
	stiffness = (shearModulus.*wireDiameter.^4)./(8.*coilDiameter.^3.*nWindings);
	% m = pi^2/4*rho*d^2*D*n
	mass = (pi^2/4).*density.*(wireDiameter.^2).*coilDiameter.*nWindings;
	% Wo = D + d
	outerDiameterOperationSpace = coilDiameter + wireDiameter;
	% Wi = D - d
	innerDiameterOperationSpace = coilDiameter - wireDiameter;
	% uc = -(s-d)*n
	deformationCompaction = -(spacingWiresRelaxed - wireDiameter).*nWindings;
	% uY = sigmaY*(2*pi*d^3)/(K*D) (tau_max = d/2 * MT/Jp = MT/(pi*d^3) = K*uY*D/(2*pi*d^3) =~ sigma_eq^max)
	deformationPlasticYield = yieldStress.*(2*pi.*wireDiameter.^3)./(stiffness.*coilDiameter);
	% uL = 0.1*D*n (change of angle: delta alpha =~ u/(Dn) << 0.1)
	deformationEndLinearity = 0.1.*coilDiameter.*nWindings;
	% H + L0 + uop (uop + min(uY,uL))
	lengthRelaxedState = spacingWiresRelaxed.*nWindings;
	heightOperationSpace = lengthRelaxedState + min(deformationPlasticYield,deformationEndLinearity);


	%% wrap quantities of interest
	performanceMeasure(:,1) = stiffness; % K
	performanceMeasure(:,2) = mass; % m
	performanceMeasure(:,3) = heightOperationSpace; % H
	performanceMeasure(:,4) = outerDiameterOperationSpace; % Wo
	performanceMeasure(:,5) = innerDiameterOperationSpace; % Wi
	performanceMeasure(:,6) = deformationCompaction; % uc
	performanceMeasure(:,7) = deformationPlasticYield; % uY
	performanceMeasure(:,8) = deformationEndLinearity; % uL
end