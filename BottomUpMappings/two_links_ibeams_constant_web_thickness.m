function [performanceMeasure,physicalFeasibilityMeasure] = two_links_ibeams_constant_web_thickness(designSample,systemParameter)
%TWO_LINKS_IBEAMS_CONSTANT_WEB_THICKNESS Bottom-up Mapping (Tip Displacement)
%   TWO_LINKS_IBEAMS_CONSTANT_WEB_THICKNESS is the bottom-up mapping function   
%   for the problem of the tip displacement of a mechanical structure with two  
%   links, both with an I-beam profile. A force of 50N is applied at the tip. 
%   Both components are assumed to have a length of 300mm, are made of steel, 
%   have a total beam height of 40mm and have a flange thickness of 3mm. The 
%   design variables are web thickness and flange width for each component.
%
%                     <------- Flange Width (W) ------->
%                   ^ +--------------------------------+ ^
%                   | |            flange              | | Flange Thickness (t)
%                   | +--------------------------------+ V
%                   |  ^            |   |
%   Beam Height (H) |  | Web      ->| w |<- Web Thickness (w)
%                   |  | Height     | e |
%                   |  | (h)        | b |
%                   |  V            |   |
%                   | +--------------------------------+
%                   | |                                |
%                   V +--------------------------------+
%
%            component 1                   component 2           |
%   |o----------------------------0----------------------------o V F
%
%   Total displacement of the tip and mass are the relevant performance 
%   indicators.
%
%   PERFORMANCEMEASURE = TWO_LINKS_IBEAMS_CONSTANT_WEB_THICKNESS(DESIGNSAMPLE) 
%   takes the minor/major lengths of each I-beam profile (web thickness /
%   flange width) specified in DESIGNSAMPLE and computes the total displacement
%   and structural mass, returning both in PERFORMANCEMEASURE. 
%
%   PERFORMANCEMEASURE = TWO_LINKS_IBEAMS_CONSTANT_WEB_THICKNESS(DESIGNSAMPLE,
%   SYSTEMPARAMETER) also allows one to choose how to handle physically 
%   infeasible designs in SYSTEMPARAMETER. If set to 'true' and the design is 
%   physically infeasible, calculations are skipped and +inf displacement and 
%   mass are returned for that design. If set to 'false', calculations proceed 
%   as normal. Default value is 'true'.
%
%   [PERFORMANCEMEASURE,PHYSICALFEASIBILITYMEASURE] = 
%   TWO_LINKS_IBEAMS_CONSTANT_WEB_THICKNESS(...) also returns the physical  
%   feasibility measure for the two components PHYSICALFEASIBILITYMEASURE. Its 
%   value is negative when the component is physically feasible, and positive 
%   otherwise.
%
%   Input: 
%       - DESIGNSAMPLE : (nSample,4) double
%           -- (1) web thickness (minor length) of component 1 (w1)
%           -- (2) flange width (major length) of component 1 (W1)
%           -- (3) web thickness (minor length) of component 2 (w2)
%           -- (4) flange width (major length) of component 2 (W2)
%       - SYSTEMPARAMETER : logical 
%
%   Output:
%       - PERFORMANCEMEASURE : (nSample,2) double
%           -- (1) total tip displacement (u)
%           -- (2) total mass of the structure (m)
%       - PHYSICALFEASIBILITYMEASURE : (nSample,2) double
%           -- (1) physical feasibility measure of component 1
%           -- (2) physical feasibility measure of component 2
%
%   See also two_links_ibeam_hollow_circle.
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

    if(nargin<2)
        systemParameter = true;
    end
    
    % Unwrap Inputs - Design Variables
    webThickness1 = designSample(:,1); % w1 [mm]
    flangeWidth1 = designSample(:,2); % W1 [mm]
    webThickness2 = designSample(:,3); % w2 [mm]
    flangeWidth2 = designSample(:,4); % W2 [mm]
    
    % Force applied at tip interface
	forceAppliedTip = [0;0;-50;0]; % F/M [N,N*mm,N,N*mm]
    
    % Component Properties - Component 1
	length1 = 300; % l1 [mm]
    youngsModulus1 = 70e3; % E1 [MPa]
    beamHeight1 = 40; % H1 [mm]
    flangeThickness1 = 3; % t1 [mm]
    webHeight1 = beamHeight1 - 2*flangeThickness1; % h1 [mm]
    density1 = 2.7e-9; % rho1 [t/mm^3]
    
    % Component Properties - Component 2
    length2 = 300; % l2 [mm]
    youngsModulus2 = 70e3; % E2 [MPa]
    beamHeight2 = 40; % H2 [mm]
    flangeThickness2 = 3; % t2 [mm]
    webHeight2 = H2 - 2*th2; % h2 [mm]
    density2 = 2.7e-9; % rho2 [t/mm^3]
    
    % Stiffness Element Array (Beam)
    % k = E*[12 , 6*l , -12 ,  4*l^2, -6*l1 , 2*l^2]/(l^3);
    stiffnessEntry1 = youngsModulus1*[12 , 6*length1 , -12 ,  4*length1^2, -6*length1 , 2*length1^2]./(length1^3);
    stiffnessEntry2 = youngsModulus2*[12 , 6*length2 , -12 ,  4*length2^2, -6*length2 , 2*length2^2]./(length2^3);
    
    % Check Physical Feasbility
    % physically feasible if major length is larger than minor length (W>=w)
    feasibility1 = webThickness1 - flangeWidth1; 
    feasibility2 = webThickness2 - flangeWidth2;
    
    % Loop over all samples to get displacement and mass for each
    tipDisplacement = nan(size(designSample,1),1);
    massTotal = nan(size(designSample,1),1);
    for i = 1:size(designSample,1)
        if (~systemParameter || (feasibility1(i)<=0 && feasibility2(i)<=0))
			% build element matrices
            % I = 1/12 * (H^3*W - h^3*(W-w))
            % K = I * [k1 k2 k3 k2]
            %         [k2 k4 k5 k6]
            %         [k3 k5 k1 k5]
            %         [k2 k6 k5 k4]
            emptyWidth1 = flangeWidth1-webThickness1(i);
            momentOfInertia1 = 1/12*(beamHeight1^3*flangeWidth1(i) - webHeight1^3*emptyWidth1);
			stiffnessMatrix1 = momentOfInertia1*...
                  [stiffnessEntry1(1) stiffnessEntry1(2) stiffnessEntry1(3) stiffnessEntry1(2)
                   stiffnessEntry1(2) stiffnessEntry1(4) stiffnessEntry1(5) stiffnessEntry1(6)
                   stiffnessEntry1(3) stiffnessEntry1(5) stiffnessEntry1(1) stiffnessEntry1(5)
                   stiffnessEntry1(2) stiffnessEntry1(6) stiffnessEntry1(5) stiffnessEntry1(4)];
			
            emptyWidth2 = flangeWidth2-webThickness2(i);
            momentOfInertia2 = 1/12*(beamHeight2^3*flangeWidth2(i) - webHeight2^3*emptyWidth2);
			stiffnessMatrix2 = momentOfInertia2*...
                  [stiffnessEntry2(1) stiffnessEntry2(2) stiffnessEntry2(3) stiffnessEntry2(2)
                   stiffnessEntry2(2) stiffnessEntry2(4) stiffnessEntry2(5) stiffnessEntry2(6)
                   stiffnessEntry2(3) stiffnessEntry2(5) stiffnessEntry2(1) stiffnessEntry2(5)
                   stiffnessEntry2(2) stiffnessEntry2(6) stiffnessEntry2(5) stiffnessEntry2(4)];
			
			% Get global matrix from the other two and solve the system
			globalStiffnessMatrix = build_global_stiffness_matrix(stiffnessMatrix1,stiffnessMatrix2);
			displacement = globalStiffnessMatrix\forceAppliedTip;
			tipDisplacement(i) = -displacement(3);
			
			% Compute total estimated mass
            % m = rho*l*A = rho*l*(2*W*t + h*w)
			mass1 = density1*length1*(2*flangeWidth1(i)*flangeThickness1 + webHeight1*webThickness1(i));
			mass2 = density2*length2*(2*flangeWidth2(i)*flangeThickness2 + webHeight2*webThickness2(i));
			massTotal(i) = mass1 + mass2;
        else
			tipDisplacement(i) = inf;
			massTotal(i) = inf;
        end
    end
    
    % Wrap quantities of interest array
	performanceMeasure(:,1) = tipDisplacement;
	performanceMeasure(:,2) = massTotal;

    % Wrap physical feasibility indicator array
    physicalFeasibilityMeasure(:,1) = feasibility1;
    physicalFeasibilityMeasure(:,2) = feasibility2;
end


%% Function to build the global stiffness matrix
function globalStiffnessMatrix = build_global_stiffness_matrix(stiffnessMatrix1,stiffnessMatrix2)
    % 2DoF per Interface * (4 Interfaces - 1 Common Interface) = 6DoF initially
	globalStiffnessMatrix = zeros(6,6);
    
    % Common interface : indexes [3,4], must sum both element matrices
	globalStiffnessMatrix(1:4,1:4) = globalStiffnessMatrix(1:4,1:4) + stiffnessMatrix1;
	globalStiffnessMatrix(3:6,3:6) = globalStiffnessMatrix(3:6,3:6) + stiffnessMatrix2;
    
    % DoF indexes [1,2] have known solution d=[0;0] since the component is
    % clamped in that interface
	globalStiffnessMatrix = globalStiffnessMatrix(3:6,3:6);
end

