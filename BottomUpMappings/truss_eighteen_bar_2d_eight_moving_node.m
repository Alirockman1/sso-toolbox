function performanceMeasure = truss_eighteen_bar_2d_eight_moving_node(designSample,systemParameter)
%TRUSS_SIX_BAR_2D_TWO_MOVING_NODE Bottom-Up Mapping (2D Node Position Problem)
%	TRUSS_SIX_BAR_2D_TWO_MOVING_NODE uses a pre-defined truss and determines the
%	total vertical displacement of the tip, total mass of the structure and 
%	stresses in the bars for different positions of the two free nodes:
%	o---0---0---0---0
%	   /|  /|  /|  /| \
%	  /	| /	| /	| /	|  o |
%	 /  |/  |/  |/  | /  V F
%	o---0---0---0---0
%	The tip node is at position (x=2,y=0.5) on this coordinate system, and
%	the force applied at the tip is 1000 downward.
%
%	PERFORMANCEMEASURE = TRUSS_SIX_BAR_2D_TWO_MOVING_NODE(DESIGNSAMPLE,
%	SYSTEMPARAMETER) receives the positions of the two nodes in DESIGNSAMPLE
%	and the cross-section area, Young's modulus and density of the bars in 
%	SYSTEMPARAMETER, returning the total vertical displacement of the tip, total
%	structure mass and the stresses in each bar in PERFORMANCEMEASURE.
%
%	Input:
%		- DESIGNSAMPLE : (nSample,4) double
%			-- (1) : horizontal position of free node 1 (nx1)
%			-- (2) : vertical position of free node 1 (ny1)
%			-- (3) : horizontal position of free node 2 (nx2)
%			-- (4) : vertical position of free node 2 (ny2)
%		- SYSTEMPARAMETER : (1,3) double OR (6,3) double 
%			-- (1) : cross-section area of each bar
%			-- (2) : Young's modulus of each bar
%			-- (3) : density of each bar
%
%	Output:
%		- PERFORMANCEMEASURE : (nSample,8) double
%			-- (1) : vertical displacement of tip node (uy)
%			-- (2) : total mass of the structure
%			-- (3) : absolute value of stress on bar 1
%			-- (4) : absolute value of stress on bar 2
%			-- (5) : absolute value of stress on bar 3
%			-- (6) : absolute value of stress on bar 4
%			-- (7) : absolute value of stress on bar 5
%			-- (8) : absolute value of stress on bar 6
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
		nan nan; ... % (1)
        nan nan; ... % (2)
        nan nan; ... % (3)
        nan nan; ... % (4)
        nan nan; ... % (5)
        nan nan; ... % (6)
        nan nan; ... % (7)
        nan nan; ... % (8)
          5 0.5; ... % (9)
          0   0; ... % (10)
          0   1];    % (11)
          
	fixedDegreesOfFreedom = [...
		false false; ... % (1)
		false false; ... % (2)
		false false; ... % (3)
		false false; ... % (4)
		false false; ... % (5)
		false false; ... % (6)
        false false; ... % (7)
		false false; ... % (8)
		false false; ... % (9)
		 true  true;...  % (10)
         true  true];    % (11)
	nodeForce = [...
		0 0; ... % (1)
		0 0; ... % (2)
        0 0; ... % (3)
		0 0; ... % (4)
        0 0; ... % (5)
		0 0; ... % (6)
        0 0; ... % (7)
		0 0; ... % (8)
		0 -1000; ... % (9)
		0 0; ... % (10)
		0 0]; % (11) assumed [N]
	nodeElement = [...
		11 1; ... % (1)
        11 2; ... % (2)
        10 2; ... % (3)
        1 2; ... % (4)
        1 4; ... % (5)
        1 3; ... % (6)
        2 4; ... % (7)
        3 4; ... % (8)
        3 6; ... % (9)
        3 5; ... % (10)
        4 6; ... % (11)
        5 6; ... % (12)
        5 8; ... % (13)
        5 7; ... % (14)
        6 8; ... % (15)
        7 8; ... % (16)
        7 9; ... % (17)
        8 9]; % (18)

	elementCrossSectionArea = systemParameter(:,1); % assumed [mm^2]
	elementYoungsModulus = systemParameter(:,2); % assumed [MPa]
	elementDensity = systemParameter(:,3); % assumed [MPa]

	nSample = size(designSample,1);
	performanceMeasure = nan(nSample,20);
	for i=1:nSample
		nodePosition(1,:) = designSample(i,[1,2]);
		nodePosition(2,:) = designSample(i,[3,4]);
        nodePosition(3,:) = designSample(i,[5,6]);
		nodePosition(4,:) = designSample(i,[7,8]);
        nodePosition(5,:) = designSample(i,[9,10]);
		nodePosition(6,:) = designSample(i,[11,12]);
        nodePosition(7,:) = designSample(i,[13,14]);
		nodePosition(8,:) = designSample(i,[15,16]);

		[nodeDisplacement,~,elementAxialForce] = ...
			truss_analysis(...
				nodePosition,...
				fixedDegreesOfFreedom,...
				nodeForce,...
				nodeElement,...
				elementCrossSectionArea,...
				elementYoungsModulus);
		elementStress = truss_deformed_stress(elementAxialForce,elementCrossSectionArea);

		nodeDistance = nodePosition(nodeElement(:,2),:)-nodePosition(nodeElement(:,1),:);
	    elementLength = vecnorm(nodeDistance,2,2);
		totalMass = sum(elementDensity.*elementLength.*elementCrossSectionArea);

		performanceMeasure(i,:) = [-nodeDisplacement(9,2),totalMass,abs(elementStress')];
	end
end