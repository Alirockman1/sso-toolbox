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


%% Cleanup
fclose all;
close all;
clear all;
clc;
more off;
diary off;


%% debugging
rng(4);


%% Documentation / Archive
rngState = rng;
saveFolder = save_diary_files(mfilename);
goldenRatio = (1+sqrt(5))/2;
figureSize = [goldenRatio 1]*8.5;


%% two-bar truss 
% information
nodePosition = [0                0;...
				0     1000*sqrt(3);...% alpha = 60°
				1000            0];... [mm]
fixedDegreesOfFreedom = [true    true; 
						 true    true; 
						 false false];
nodeForce = [0      0; 
		     0      0; 
			 0 -1000]; % [N]
nodeElement = [1 3; ...
			   2 3];
elementCrossSectionArea = 100; % [mm^2]
elementYoungsModulus = 210e3; % [MPa]


% function call
fprintf([repelem('=',80),'\n']);
fprintf('TWO-BAR 2D TRUSS\n');
[nodeDisplacement,nodeReactionForce,elementAxialForce] = ...
	truss_analysis(...
		nodePosition,...
		fixedDegreesOfFreedom,...
		nodeForce,...
		nodeElement,...
		elementCrossSectionArea,...
		elementYoungsModulus) % display
elementStrain = truss_deformed_strain(nodePosition,nodeDisplacement,nodeElement) % display
elementStress = truss_deformed_stress(elementAxialForce,elementCrossSectionArea) % display


% plot results
figure;
[handleUndeformed,handleDeformed] = plot_truss_deformation(gcf,nodePosition,nodeElement,nodeDisplacement,...
    'DisplacementScaleFactor',1000);
grid minor;
legend([handleUndeformed,handleDeformed],{'Undeformed Truss','Deformed Truss'},'location','northeast');

plot_truss_element_response(nodePosition,nodeElement,elementStrain);
title('Strain');
grid minor;
xlabel('X');
ylabel('Y');

plot_truss_element_response(nodePosition,nodeElement,elementStress);
title('Stress');
grid minor;
xlabel('X');
ylabel('Y');


%% six-bar 2D truss 
% information
nodePosition = [0 0; 0 1; 1 0; 1 1; 2 0.5];
fixedDegreesOfFreedom = [true true; true true; false false; false false; false false];
nodeForce = [0 0; 0 0; 0 0; 0 0; 0 -1000];
nodeElement = [1 3; 2 4; 2 3; 3 4; 3 5; 4 5];


% function call
fprintf([repelem('=',80),'\n']);
fprintf('SIX-BAR 2D TRUSS\n');
[nodeDisplacement,nodeReactionForce,elementAxialForce] = ...
	truss_analysis(...
		nodePosition,...
		fixedDegreesOfFreedom,...
		nodeForce,...
		nodeElement,...
		elementCrossSectionArea,...
		elementYoungsModulus) % display
elementStrain = truss_deformed_strain(nodePosition,nodeDisplacement,nodeElement) % display
elementStress = truss_deformed_stress(elementAxialForce,elementCrossSectionArea) % display


% plot results
figure;
[handleUndeformed,handleDeformed] = plot_truss_deformation(gcf,nodePosition,nodeElement,nodeDisplacement,...
    'DisplacementScaleFactor',1000);
grid minor;
legend([handleUndeformed,handleDeformed],{'Undeformed Truss','Deformed Truss'},'location','northeast');

plot_truss_element_response(nodePosition,nodeElement,elementStrain);
title('Strain');
grid minor;
xlabel('X');
ylabel('Y');

plot_truss_element_response(nodePosition,nodeElement,elementStress);
title('Stress');
grid minor;
xlabel('X');
ylabel('Y');


%% three-bar 3D truss 
% information
nodePosition = [0 0 0; 0 0.5 1; 0 1 0; 0.5 0.5 0.5];
fixedDegreesOfFreedom = [true true true; true true true; true true true; false false false];
nodeForce = [0 0 0; 0 0 0; 0 0 0; 0 0 -1000];
nodeElement = [1 4; 2 4; 3 4];


% function call
fprintf([repelem('=',80),'\n']);
fprintf('THREE-BAR 3D TRUSS\n');
[nodeDisplacement,nodeReactionForce,elementAxialForce] = ...
	truss_analysis(...
		nodePosition,...
		fixedDegreesOfFreedom,...
		nodeForce,...
		nodeElement,...
		elementCrossSectionArea,...
		elementYoungsModulus) % display
elementStrain = truss_deformed_strain(nodePosition,nodeDisplacement,nodeElement) % display
elementStress = truss_deformed_stress(elementAxialForce,elementCrossSectionArea) % display


% plot results
figure;
[handleUndeformed,handleDeformed] = plot_truss_deformation(gcf,nodePosition,nodeElement,nodeDisplacement,...
    'DisplacementScaleFactor',1000);
grid minor;
legend([handleUndeformed,handleDeformed],{'Undeformed Truss','Deformed Truss'},'location','northeast');
xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);

plot_truss_element_response(nodePosition,nodeElement,elementStrain);
title('Strain');
grid minor;
xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);

plot_truss_element_response(nodePosition,nodeElement,elementStress);
title('Stress');
grid minor;
xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);


%% twelve-bar 3D truss 
% information
nodePosition = [0 0 0; % (1)
                0 0 1; % (2)
                0 1 0; % (3)
                1 0 0; % (4)
                1 0 1; % (5) 
                1 1 0; % (6)
                2 0.5 0.5]; % (7)
fixedDegreesOfFreedom = [true true true; % (1) 
                         true true true; % (2)
                         true true true; % (3)
                         false false false; % (4)
                         false false false; % (5)
                         false false false; % (6)
                         false false false]; % (7)
nodeForce = [0 0 0; % (1)
             0 0 0; % (2)
             0 0 0; % (3)
             0 0 0; % (4)
             0 0 0; % (5)
             0 0 0; % (6)
             0 0 -1000]; % (7)
nodeElement = [1 4; % (1)
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


% function call
fprintf([repelem('=',80),'\n']);
fprintf('TWELVE-BAR 3D TRUSS\n');
[nodeDisplacement,nodeReactionForce,elementAxialForce] = ...
	truss_analysis(...
		nodePosition,...
		fixedDegreesOfFreedom,...
		nodeForce,...
		nodeElement,...
		elementCrossSectionArea,...
		elementYoungsModulus) % display
elementStrain = truss_deformed_strain(nodePosition,nodeDisplacement,nodeElement) % display
elementStress = truss_deformed_stress(elementAxialForce,elementCrossSectionArea) % display


% plot results
figure;
[handleUndeformed,handleDeformed] = plot_truss_deformation(gcf,nodePosition,nodeElement,nodeDisplacement,...
    'DisplacementScaleFactor',1000);
grid minor;
legend([handleUndeformed,handleDeformed],{'Undeformed Truss','Deformed Truss'},'location','northeast');
xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);

plot_truss_element_response(nodePosition,nodeElement,elementStrain);
title('Strain');
grid minor;
xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);

plot_truss_element_response(nodePosition,nodeElement,elementStress);
title('Stress');
grid minor;
xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);

