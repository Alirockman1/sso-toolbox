function [performanceMeasure,physicalFeasibilityMeasure] = two_links_ibeam_hollow_circle(designSample,systemParameter)
%TWO_LINKS_IBEAM_HOLLOW_CIRCLE Bottom-up Mapping (2-Links Tip Displacement)
%   TWO_LINKS_IBEAM_HOLLOW_CIRCLE is the bottom-up mapping function for the 
%   problem of the tip displacement of a mechanical structure with two links, 
%   one with a I-beam profile and the other with a hollow circular beam profile.
%   A force of 50N is applied at the tip. Both components are assumed to have
%   a length of 300mm and made of steel. The I-beam component has a total height
%   of 40mm and flange thickness of 3mm, and its design variables are web 
%   thickness and flange width.
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
%   The hollow circle's design variables are the inner and outer diameter of 
%   the cross-section of the bar.
%              I-beam                    hollow circle           |
%   |o----------------------------0----------------------------o V F
%
%   Total displacement of the tip and mass are the relevant performance 
%   indicators.
%
%   PERFORMANCEMEASURE = TWO_LINKS_IBEAM_HOLLOW_CIRCLE(DESIGNSAMPLE) takes the 
%   minor/major lengths of the I-beam profile (web thickness / flange width) and
%   internal/external radii of the hollow circle profile specified in 
%   DESIGNSAMPLE and computes the total displacement and structural mass, 
%   returning both in PERFORMANCEMEASURE. 
%
%   PERFORMANCEMEASURE = 
%   TWO_LINKS_IBEAM_HOLLOW_CIRCLE(DESIGNSAMPLE,SYSTEMPARAMETER) also allows one
%   to choose how to handle physically infeasible designs in SYSTEMPARAMETER.
%   If set to 'true' and the design is physically infeasible, calculations are 
%   skipped and +inf displacement and mass are returned for that design. If
%   set to 'false', calculations proceed as normal. Default value is 'true'.
%
%   [PERFORMANCEMEASURE,PHYSICALFEASIBILITYMEASURE] = 
%   TWO_LINKS_IBEAM_HOLLOW_CIRCLE(...) also returns the physical feasibility 
%   measure for the two components PHYSICALFEASIBILITYMEASURE. Its value is 
%   negative when the component is physically feasible, and positive otherwise.
%
%   Input: 
%       - DESIGNSAMPLE : (nSample,4) double
%           -- (1) web thickness (minor length) of I-beam profile (w)
%           -- (2) flange width (major length) of I-beam profile (W)
%           -- (3) internal radius of hollow circular profile (r)
%           -- (4) external radius of hollow circular profile (R)
%       - SYSTEMPARAMETER : logical 
%
%   Output:
%       - PERFORMANCEMEASURE : (nSample,2) double
%           -- (1) total tip displacement (u)
%           -- (2) total mass of the structure (m)
%       - PHYSICALFEASIBILITYMEASURE : (nSample,2) double
%           -- (1) physical feasibility measure of I-beam component
%           -- (2) physical feasibility measure of hollow circular profile
%
%   See also two_links_ibeams_fixed_thickness.
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
    
    % Unwrap Inputs
    webThickness1 = designSample(:,1); % w1 [mm]
    flangeWidth1 = designSample(:,2); % W1 [mm]
    innerRadius2 = designSample(:,3); % r1 [mm]
    outerRadius2 = designSample(:,4); % R2 [mm]
    
    % Force applied at tip interface
	forceAppliedTip = [0;0;-50;0]; % F/M [N,N*mm,N,N*mm]
    
    % Component Properties - Component 1 (I-beam)
	length1 = 300; % l1 [mm]
    youngsModulus1 = 70e3; % E1 [MPa]
    beamHeight1 = 40; % H1 [mm]
    flangeThickness1 = 3; % t1 [mm]
    webHeight1 = beamHeight1 - 2*flangeThickness1; % h1[mm]
    density1 = 2.7e-9; % rho1 [t/mm^3]
    
    % Component Properties - Component 2 (Circular)
    length2 = 300; % [mm]
    youngsModulus2 = 70e3; % [MPa]
    density2 = 2.7e-9; % [t/mm^3]
    
    % Stiffness Element Array (Beam)
    % k = E*[12 , 6*l , -12 ,  4*l^2, -6*l1 , 2*l^2]/(l^3);
    stiffnessEntry1 = youngsModulus1*[12 , 6*length1 , -12 ,  4*length1^2, -6*length1 , 2*length1^2]./(length1^3);
    stiffnessEntry2 = youngsModulus2*[12 , 6*length2 , -12 ,  4*length2^2, -6*length2 , 2*length2^2]./(length2^3);
    
    % Check Feasbility
    % I-beam: physically feasible if major length is larger than minor length (W>=w)
    feasibility1 = webThickness1 - flangeWidth1; 
    % Hollow circle: physically feasible if outer radius is larger than inner radius (R>=r)
    feasibility2 = innerRadius2 - outerRadius2; 
    
    % Loop over all samples to get displacement and mass for each
    tipDisplacement = nan(size(webThickness1,1),1);
    massTotal = nan(size(webThickness1,1),1);
    for i = 1:size(webThickness1,1)
        if (~systemParameter || (feasibility1(i)<=0 && feasibility2(i)<=0))
			% build element matrices
            % I-beam: I = 1/12 * (H^3*W - h^3*(W-w))
            emptyWidth1 = flangeWidth1-webThickness1(i);
            momentOfInertia1 = 1/12*(beamHeight1^3*flangeWidth1(i) - webHeight1^3*emptyWidth1); % I-beam

            % Hollow cirle: I = pi/4 * (R^4 - r^4)
            momentOfInertia2 = pi/4*(outerRadius2(i)^4 - innerRadius2(i)^4); % hollow circle

            % K = I * [k1 k2 k3 k2]
            %         [k2 k4 k5 k6]
            %         [k3 k5 k1 k5]
            %         [k2 k6 k5 k4]
			stiffnessMatrix1 = momentOfInertia1*...
                [stiffnessEntry1(1) stiffnessEntry1(2) stiffnessEntry1(3) stiffnessEntry1(2)
                 stiffnessEntry1(2) stiffnessEntry1(4) stiffnessEntry1(5) stiffnessEntry1(6)
                 stiffnessEntry1(3) stiffnessEntry1(5) stiffnessEntry1(1) stiffnessEntry1(5)
                 stiffnessEntry1(2) stiffnessEntry1(6) stiffnessEntry1(5) stiffnessEntry1(4)];
			
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
            % m = rho*l*A
            % I-beam: A = 2*W*t + h*w
			mass1 = density1*length1*(2*flangeWidth1(i)*flangeThickness1 + webHeight1*webThickness1(i)); % I-beam
            % Hollow cirle: A = pi*(R^2 - r^2)
			mass2 = density2*length2*pi*(outerRadius2(i)^2 - innerRadius2(i)^2); % circular
			massTotal(i) = mass1 + mass2;
        else
			tipDisplacement(i) = inf;
			massTotal(i) = inf;
        end
    end
    
    % Wrap quantities of interest array
	performanceMeasure(:,1) = tipDisplacement;
	performanceMeasure(:,2) = massTotal;
    
    % Wrap physical feasibility array
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

