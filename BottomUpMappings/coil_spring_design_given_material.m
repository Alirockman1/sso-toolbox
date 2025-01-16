function performanceMeasure = coil_spring_design_given_material(designSample,systemParameter)
%COIL_SPRING_DESIGN_GIVEN_MATERIAL Bottom-up Mapping (Coil Spring Problem)
%	COIL_SPRING_DESIGN_GIVEN_MATERIAL works as a wrapper for 
%	'coil_spring_design', where the material properties are fixed.
%
%	PERFORMANCEMEASURE = COIL_SPRING_DESIGN_GIVEN_MATERIAL(DESIGNSAMPLE,
%	SYSTEMPARAMETER) receives the diameters, number of windings and spacing of 
%	wires in DESIGNSAMPLE, and material properties of the spring in 
%	SYSTEMPARAMETER, and returns the coil spring characteristics in 
%	PERFORMANCEMEASURE.
%
%	Input:
%       - DESIGNSAMPLE : (nSample,4) double
%           -- (1) : wire diameter (d)
%           -- (2) : coil diameter (D)
%           -- (3) : number of windings (n)
%           -- (4) : spacing of wires in relaxed state (s)
%		- SYSTEMPARAMETER : (1,3) double
%           -- (1) : material density (rho)
%           -- (2) : material shear modulus (G)
%           -- (3) : material yield stres (sigmaY)
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
%   See also coil_spring_design.
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

	nSample = size(designSample,1);
	performanceMeasure = coil_spring_design([designSample,repmat(systemParameter,nSample,1)]);
end