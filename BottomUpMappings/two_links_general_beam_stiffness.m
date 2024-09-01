function [performanceMeasure,physicalFeasibilityMeasure] = two_links_general_beam_stiffness(designSample,systemParameter)
%TWO_LINKS_GENERAL_BEAM_STIFFNESS Bottom-Up Mapping (2-Links Tip Displacement)
%   TWO_LINKS_GENERAL_BEAM_STIFFNESS is the bottom-up mapping function for the 
%   problem of the tip displacement of a mechanical structure with two links.
%	As design variables, the main stiffness values for each component are given. 
%	The displacement is analytically calculated with the given stiffnesses, and 
%	the mass for each component is estimated with an Aritificial Neural Network.
%   More on how the mass estimation works given the desired stiffness can be 
%   found on Dr. Lukas Krischer's paper "Decomposition and optimization of 
%	linear structures using meta-models".
%
%            component 1                   component 2           |
%   |o----------------------------0----------------------------o V F
%
%   A force of 50N is applied at the tip. Total tip displacement and total mass   
%   are the quantities of interest.
%
%   PERFORMANCEMEASURE = TWO_LINKS_GENERAL_BEAM_STIFFNESS(DESIGNSAMPLE,
%	SYSTEMPARAMETER) receives the desired main stiffness values for each  
%	component specified in DESIGNSAMPLE and the files with the mass / physical 
%	feasibility estimator in SYSTEMPARAMETER, and uses those to compute the   
%	total displacement and structural mass, returning both in 
%	PERFORMANCEMEASURE. Additionally, SYSTEMPARAMETER can have a second entry
%	which allows one to choose how to handle physically infeasible designs. If 
%	set to 'true' and the design is physically infeasible, calculations are 
%	skipped and +inf displacement and mass are returned for that design. If set 
%	to 'false', calculations proceed as normal. Default value is 'true'.
%
%   [PERFORMANCEMEASURE,PHYSICALFEASIBILITYMEASURE] = 
%   TWO_LINKS_GENERAL_BEAM_STIFFNESS(...) also returns the physical  
%   feasibility measure for the two components PHYSICALFEASIBILITYMEASURE. Its 
%   value is negative when the component is physically feasible, and positive 
%   otherwise.
%
%   Input: 
%       - DESIGNSAMPLE : (nSample,6) double
%           -- (1) : Stiffness value k11 of component 1
%           -- (2) : Stiffness value k22 of component 1
%           -- (3) : Stiffness value k44 of component 1
%           -- (4) : Stiffness value k11 of component 2
%			-- (5) : Stiffness value k22 of component 2
%           -- (6) : Stiffness value k44 of component 2
%       - SYSTEMPARAMETER : char OR struct OR cell OR (1,2) cell
%			-- (1) : name of file with estimators OR structure with estimators
%			in fields named 'PhysicalFeasibilityEstimator' (ClassificationSVM 
%			or similar) and 'MassEstimator' (network).
%			-- (2) : logical flag to skip physically infeasible designs
%
%   Output:
%       - PERFORMANCEMEASURE : (nSample,2) double
%           -- (1) total tip displacement (u)
%           -- (2) total mass of the structure (m)
%       - PHYSICALFEASIBILITYMEASURE : (nSample,2) double
%           -- (1) physical feasibility measure of component 1
%           -- (2) physical feasibility measure of component 2
%
%   See also two_links_ibeams_constant_web_thickness.
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

    persistent estimatorFile

	skipPhysicallyInfeasible = true;
	if(iscell(systemParameter))
		estimatorFileName = systemParameter{1};
        physicalFeasibilityEstimatorName = systemParameter{2};
        massEstimatorName = systemParameter{3};

		if(length(systemParameter)==4)
			skipPhysicallyInfeasible = systemParameter{4};
		end
	else
		estimatorFileName = systemParameter;
	end
    
    
	% Unwrap inputs
	stiffnessK11Component1 = designSample(:,1); % [N/mm]
	stiffnessK22Component1 = designSample(:,2); % [N*mm/rad]
	stiffnessK44Component1 = designSample(:,3); % [N*mm/rad]
	stiffnessK11Component2 = designSample(:,4); % [N/mm]
	stiffnessK22Component2 = designSample(:,5); % [N*mm/rad]
	stiffnessK44Component2 = designSample(:,6); % [N*mm/rad]
	
    % Load physicalFeasibilityEstimator classifier which indicates physical feasibility
    if(isempty(estimatorFile))
        if(~isstruct(estimatorFileName))
		    estimatorFile = load(estimatorFileName);
        else
            estimatorFile = estimatorFileName;
        end
    end
	
    % Force applied at tip interface / Component length
	appliedForceTip = [0;0;-50;0]; %[N,N*mm,N,N*mm]
	length1 = 300; % [mm]
	length2 = 300; % [mm]
	
    % Check physical feasibility of each component
	physicalFeasibilityScore1 = estimate_physical_feasibility(stiffnessK11Component1,...
		stiffnessK22Component1,stiffnessK44Component1,estimatorFile.(physicalFeasibilityEstimatorName));
    physicalFeasibilityScore2 = estimate_physical_feasibility(stiffnessK11Component2,...
    	stiffnessK22Component2,stiffnessK44Component2,estimatorFile.(physicalFeasibilityEstimatorName));
    
    % Save physical feasibility vector
    nSample = size(designSample,1);
    physicalFeasibilityMeasure = nan(nSample,2);
    physicalFeasibilityMeasure(:,1) = physicalFeasibilityScore1(:,1);
    physicalFeasibilityMeasure(:,2) = physicalFeasibilityScore2(:,1);

    % Compute total estimated mass
	mass1 = estimate_mass(stiffnessK11Component1,stiffnessK22Component1,...
		stiffnessK44Component1,estimatorFile.(massEstimatorName));
	mass2 = estimate_mass(stiffnessK11Component2,stiffnessK22Component2,...
		stiffnessK44Component2,estimatorFile.(massEstimatorName));
	massTotal = mass1 + mass2;

    % Get remaining values of each stiffness matrix
    [stiffnessK21Component1,stiffnessK31Component1,stiffnessK41Component1,...
    	stiffnessK32Component1,stiffnessK42Component1,stiffnessK33Component1,stiffnessK43Component1] = ...
    	build_element_stiffness_matrix(length1,stiffnessK11Component1,stiffnessK22Component1,stiffnessK44Component1);
    [stiffnessK21Component2,stiffnessK31Component2,stiffnessK41Component2,...
    	stiffnessK32Component2,stiffnessK42Component2,stiffnessK33Component2,stiffnessK43Component2] = ...
    	build_element_stiffness_matrix(length2,stiffnessK11Component2,stiffnessK22Component2,stiffnessK44Component2);
    
    % Loop over all samples to get displacement for each
	tipDisplacement = nan(nSample,1);
	for i = 1:nSample
		if (~skipPhysicallyInfeasible || (physicalFeasibilityScore1(i)<=0 && physicalFeasibilityScore2(i)<=0))
			% build element matrices
			stiffnessMatrix1 = [...
				stiffnessK11Component1(i), stiffnessK21Component1(i), stiffnessK31Component1(i), stiffnessK41Component1(i);
				stiffnessK21Component1(i), stiffnessK22Component1(i), stiffnessK32Component1(i), stiffnessK42Component1(i);
				stiffnessK31Component1(i), stiffnessK32Component1(i), stiffnessK33Component1(i), stiffnessK43Component1(i);
				stiffnessK41Component1(i), stiffnessK42Component1(i), stiffnessK43Component1(i), stiffnessK44Component1(i)];
			  
			stiffnessMatrix2 = [...
				stiffnessK11Component2(i), stiffnessK21Component2(i), stiffnessK31Component2(i), stiffnessK41Component2(i);
				stiffnessK21Component2(i), stiffnessK22Component2(i), stiffnessK32Component2(i), stiffnessK42Component2(i);
				stiffnessK31Component2(i), stiffnessK32Component2(i), stiffnessK33Component2(i), stiffnessK43Component2(i);
				stiffnessK41Component2(i), stiffnessK42Component2(i), stiffnessK43Component2(i), stiffnessK44Component2(i)];
			
			% Get global matrix from the other two and solve the system
			globalStiffnessMatrix = build_global_stiffness_matrix(stiffnessMatrix1,stiffnessMatrix2);
			displacement = globalStiffnessMatrix\appliedForceTip;
			tipDisplacement(i) = -displacement(3); % d = d1 + theta1*l2 + d2
		else
			tipDisplacement(i) = inf; 
			massTotal(i)  = inf;
		end
	end
	
	% Wrap outputs
	performanceMeasure = nan(nSample,2);
	performanceMeasure(:,1) = tipDisplacement;
	performanceMeasure(:,2) = massTotal;
end


%% Function wrapper to evaluate physical feasibility
function physicallyFeasibilityScore = estimate_physical_feasibility(stiffnessK11,stiffnessK22,stiffnessK44,physicalFeasibilityEstimator)
	[~,physicallyFeasibilityScore] = predict(physicalFeasibilityEstimator,[stiffnessK11(:),stiffnessK22(:),stiffnessK44(:)]);
end


%% Function wrapper to estimate mass
function massComponent = estimate_mass(stiffnessK11,stiffnessK22,stiffnessK44,massEstimator)
	massComponent = massEstimator([stiffnessK11,stiffnessK22,stiffnessK44]')';
end


%% Function to build the element stiffness matrix from its given main values
function [stiffnessK21,stiffnessK31,stiffnessK41,stiffnessK32,stiffnessK42,stiffnessK33,stiffnessK43] = ...
	build_element_stiffness_matrix(lengthComponent,stiffnessK11,stiffnessK22,stiffnessK44)
	stiffnessK21(:) =  (stiffnessK11.*lengthComponent.^2 + stiffnessK22 - stiffnessK44)./(2.*lengthComponent);
	stiffnessK31(:) = -stiffnessK11;
	stiffnessK41(:) =  (stiffnessK11.*lengthComponent.^2 - stiffnessK22 + stiffnessK44)./(2.*lengthComponent);
	stiffnessK32(:) = -(stiffnessK11.*lengthComponent.^2 + stiffnessK22 - stiffnessK44)./(2.*lengthComponent);
	stiffnessK42(:) =  (stiffnessK11.*lengthComponent.^2 - stiffnessK22 - stiffnessK44)./(2);
	stiffnessK33(:) =  stiffnessK11;
	stiffnessK43(:) = -(stiffnessK11.*lengthComponent.^2 - stiffnessK22 + stiffnessK44)./(2.*lengthComponent);
end


%% Function to build the global stiffness matrix
function globalStiffnessMatrix = build_global_stiffness_matrix(K1,K2)
    % 2DoF per Interface * (4 Interfaces - 1 Common Interface) = 6DoF initially
	globalStiffnessMatrix = zeros(6,6);
    
    % Common interface : indexes [3,4], must sum both element matrices
	globalStiffnessMatrix(1:4,1:4) = globalStiffnessMatrix(1:4,1:4) + K1;
	globalStiffnessMatrix(3:6,3:6) = globalStiffnessMatrix(3:6,3:6) + K2;
    
    % DoF indexes [1,2] have known solution d=[0;0] since the component is
    % clamped in that interface
	globalStiffnessMatrix = globalStiffnessMatrix(3:6,3:6);
end

