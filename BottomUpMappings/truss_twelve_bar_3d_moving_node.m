function performanceMeasure = truss_twelve_bar_3d_moving_node(designSample,systemParameter)
%TRUSS_SIX_BAR_2D_MOVING_NODE Bottom-Up Mapping (2D Node Position Problem)
%	TRUSS_SIX_BAR_2D_MOVING_NODE uses a pre-defined truss and determines the
%	total vertical displacement of the tip and stresses in the bars for 
%	different positions of the two free nodes:
%	o---0
%	   /|\
%	  /	| o | 
%	 /  |/	V F	  
%	o---0
%	The tip node is at position (x=2,y=0.5) on this coordinate system, and
%	the force applied at the tip is 1000 downward.
%
%	PERFORMANCEMEASURE = TRUSS_SIX_BAR_2D_MOVING_NODE(DESIGNSAMPLE,
%	SYSTEMPARAMETER) receives the positions of the two nodes in DESIGNSAMPLE
%	and the cross-section area and Young's modulus of the bars in 
%	SYSTEMPARAMETER, returning the total vertical displacement of the tip and 
%	the stresses in each bar in PERFORMANCEMEASURE.
%
%	Input:
%		- DESIGNSAMPLE : (nSample,4) double
%			-- (1) : horizontal position of free node 1 (nx1)
%			-- (2) : vertical position of free node 1 (ny1)
%			-- (3) : horizontal position of free node 2 (nx2)
%			-- (4) : vertical position of free node 2 (ny2)
%		- SYSTEMPARAMETER : (1,2) double OR (6,2) double 
%			-- (1) : cross-section area of each bar
%			-- (2) : Young's modulus of each bar
%
%	Output:
%		- PERFORMANCEMEASURE : (nSample,7) double
%			-- (1) : vertical displacement of tip node (uy)
%			-- (2) : absolute value of stress on bar 1
%			-- (3) : absolute value of stress on bar 2
%			-- (4) : absolute value of stress on bar 3
%			-- (5) : absolute value of stress on bar 4
%			-- (6) : absolute value of stress on bar 5
%			-- (7) : absolute value of stress on bar 6
%
%	See also truss_analysis, truss_two_bar_hollow_circle.
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

	nodePosition = [...
          0   0   0; % (1)
          0   0   1; % (2)
          0   1   0; % (3)
        nan nan nan; % (4)
        nan nan nan; % (5) 
        nan nan nan; % (6)
          2 0.5 0.5]; % (7)
    fixedDegreesOfFreedom = [...
        true true true; % (1) 
        true true true; % (2)
        true true true; % (3)
        false false false; % (4)
        false false false; % (5)
        false false false; % (6)
        false false false]; % (7)
    nodeForce = [...
        0 0     0; % (1)
        0 0     0; % (2)
        0 0     0; % (3)
        0 0     0; % (4)
        0 0     0; % (5)
        0 0     0; % (6)
        0 0 -1000]; % (7)
    nodeElement = [...
        1 4; % (1)
        2 5; % (2)
        3 6; % (3)
        1 5; % (4)
        3 4; % (5)
        3 5; % CONFIRM WITH ZM
        4 5; % (6)
        5 6; % (7)
        4 6; % (8)
        4 7; % (9)
        5 7; % (10)
        6 7]; % (11)
	elementCrossSectionArea = systemParameter(:,1); % assumed [mm^2]
	elementYoungsModulus = systemParameter(:,2); % assumed [MPa]

	nSample = size(designSample,1);
	performanceMeasure = nan(nSample,13);
	for i=1:nSample
		nodePosition(4,:) = designSample(i,[1,2,3]);
		nodePosition(5,:) = designSample(i,[4,5,6]);
        nodePosition(6,:) = designSample(i,[7,8,9]);

		[nodeDisplacement,~,elementAxialForce] = ...
			truss_analysis(...
				nodePosition,...
				fixedDegreesOfFreedom,...
				nodeForce,...
				nodeElement,...
				elementCrossSectionArea,...
				elementYoungsModulus);
		elementStress = truss_deformed_stress(elementAxialForce,elementCrossSectionArea);

		performanceMeasure(i,:) = [-nodeDisplacement(7,3),abs(elementStress')];
	end
end