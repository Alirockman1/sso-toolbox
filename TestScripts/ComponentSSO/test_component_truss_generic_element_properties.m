%TEST_COMPONENT_TRUSS_GENERIC_MOVABLE_NODE Truss design problem w/ component SSO 
%   TEST_COMPONENT_TRUSS_GENERIC_MOVABLE_NODE
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


%% Considered material properties
% Material properties arrays
youngsModulusMaterial = [200e9, 70e9, 116e9, 45e9]; % [Pa] Steel, Aluminum, Titanium, Magnesium
densityMaterial = [7850, 2700, 4500, 1740];         % [kg/m^3] Steel, Aluminum, Titanium, Magnesium
yieldStrengthMaterial = [250e6, 110e6, 1200e6, 200e6]; % [Pa] Steel, Aluminum, Titanium, Magnesium


[youngsModulusMaterial,materialOrder] = sort(youngsModulusMaterial);
densityMaterial = densityMaterial(materialOrder);
yieldStrengthMaterial = yieldStrengthMaterial(materialOrder);

figure;
plot(youngsModulusMaterial,densityMaterial,'-o');
xlabel('Young''s Modulus [Pa]');
ylabel('Density [kg/m^3]');
save_print_figure(gcf,[saveFolder,'MaterialRelationEtoRho'],'PrintFormat',{'png','pdf'},'Size',figureSize);

figure;
plot(youngsModulusMaterial,yieldStrengthMaterial,'-o');
xlabel('Young''s Modulus [Pa]');
ylabel('Yield Strength [Pa]');
save_print_figure(gcf,[saveFolder,'MaterialRelationEtoSigmaY'],'PrintFormat',{'png','pdf'},'Size',figureSize);


%%
trussAnalysisChoice = '4-Bar-2D';
fixRadius = true;
useBoxResultForComponent = false;

computeDisplacement = true;
computeMass = false;
computeDisplacementAndMass = false;

computeBoxSolution = true;
computePlanarTrimmingComponent = true;
computeCornerBoxRemovalComponent = true;
computeHolePunchingComponent = false;


%% function call
systemFunction = @truss_generic_element_properties_dependent_density;

switch trussAnalysisChoice
    case '2-Bar-2D'
        nSample = 500;
        maxIterDisplacement = 150;
        maxIterMass = 20;
        maxIterDisplacementAndMass = 100;
        growthRateDisplacement = 0.01;
        growthRateMass = 0.07;
        trimmingPasses = 'single';
        requirementSpacesType = 'Omega1';
        systemParameter.NodePosition = [...
            0   0;  ...
            0   1; ...
            1 0.5];
        systemParameter.FixedDegreesOfFreedom = [...
             true  true; ...
             true  true;
            false false];
        systemParameter.NodeForce = [...
            0 0; ...
            0 0; ...
            0 -1000];
        systemParameter.NodeElement = [...
            1 3; ...
            2 3];
    case '4-Bar-2D'
        nSample = 400;
        maxIterDisplacement = 150;
        maxIterMass = 20;
        maxIterDisplacementAndMass = 200;
        growthRateDisplacement = 0.07;
        growthRateMass = 0.07;
        trimmingPasses = 'single';
        requirementSpacesType = 'Omega1';
        systemParameter.NodePosition = [...
            0 0;  ...
            1 0; ...
            1 1; ...
            0 1];
        systemParameter.FixedDegreesOfFreedom = [...
            true true; ...
            false false; ...
            false false; ...
            true true];
        systemParameter.NodeForce = [...
            0 0; ...
            0 0; ...
            0 -1000; ...
            0 0];
        systemParameter.NodeElement = [...
            1 2; ...
            2 3; ...
            3 4; ...
            2 4];
    case '6-Bar-2D'
        nSample = 600;
        maxIterDisplacement = 200;
        maxIterMass = 20;
        maxIterDisplacementAndMass = 30;
        growthRateDisplacement = 0.07;
        growthRateMass = 0.07;
        trimmingPasses = 'single';
        requirementSpacesType = 'Omega1';
        systemParameter.NodePosition = [...
            0 0;  ...
            1 0; ...
            2 0.5; ...
            1 1; ...
            0 1];
        systemParameter.FixedDegreesOfFreedom = [...
            true true; ...
            false false; ...
            false false; ...
            false false; ...
            true true];
        systemParameter.NodeForce = [...
            0 0; ...
            0 0; ...
            0 -1000; ...
            0 0; ...
            0 0];
        systemParameter.NodeElement = [...
            1 2; ...
            2 3; ...
            3 4; ...
            4 5; ...
            2 5; ...
            2 4];
    case '18-Bar-2D'
        nSample = 300;
        maxIterDisplacement = 500;
        maxIterMass = 250;
        maxIterDisplacementAndMass = 250;
        growthRateDisplacement = 0.01;
        growthRateMass = 0.01;
        trimmingPasses = 'single';
        requirementSpacesType = 'Omega1';
        systemParameter.NodePosition = [...
	          1, 1; ... % (1)
              1, 0; ... % (2)
              2, 1; ... % (3)
              2, 0; ... % (4)
              3, 1; ... % (5)
              3, 0; ... % (6)
              4, 1; ... % (7)
              4, 0; ... % (8)
              5 0.5; ... % (9)
              0   0; ... % (10)
              0   1];    % (11)
        systemParameter.FixedDegreesOfFreedom = [...
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
        systemParameter.NodeForce = [...
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
        systemParameter.NodeElement = [...
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
    case '38-Bar-2D'
        nSample = 300;
        maxIterDisplacement = 500;
        maxIterMass = 300;
        maxIterDisplacementAndMass = 300;
        growthRateDisplacement = 0.007;
        growthRateMass = 0.007;
        trimmingPasses = 'reduced';
        requirementSpacesType = 'Omega1';
        systemParameter.NodePosition = [...
              1, 1; ... % (1)
              1, 0; ... % (2)
              2, 1; ... % (3)
              2, 0; ... % (4)
              3, 1; ... % (5)
              3, 0; ... % (6)
              4, 1; ... % (7)
              4, 0; ... % (8)
              5, 1; ... % (9)
              5, 0; ... % (10)
              6, 1; ... % (11)
              6, 0; ... % (12)
              7, 1; ... % (13)
              7, 0; ... % (14)
              8, 1; ... % (15)
              8, 0; ... % (16)
              9, 1; ... % (17)
              9, 0; ... % (18)
             10 0.5; ... % (19)
              0   0; ... % (20)
              0   1];    % (21)
        systemParameter.FixedDegreesOfFreedom = [...
            false false; ... % (1)
            false false; ... % (2)
            false false; ... % (3)
            false false; ... % (4)
            false false; ... % (5)
            false false; ... % (6)
            false false; ... % (7)
            false false; ... % (8)
            false false; ... % (9)
            false false; ... % (10)
            false false; ... % (11)
            false false; ... % (12)
            false false; ... % (13)
            false false; ... % (14)
            false false; ... % (15)
            false false; ... % (16)
            false false; ... % (17)
            false false; ... % (18)
            false false; ... % (19)
             true  true;...  % (20)
             true  true];    % (21)
        systemParameter.NodeForce = [...
            0 0; ... % (1)
            0 0; ... % (2)
            0 0; ... % (3)
            0 0; ... % (4)
            0 0; ... % (5)
            0 0; ... % (6)
            0 0; ... % (7)
            0 0; ... % (8)
            0 0; ... % (9)
            0 0; ... % (10)
            0 0; ... % (11)
            0 0; ... % (12)
            0 0; ... % (13)
            0 0; ... % (14)
            0 0; ... % (15)
            0 0; ... % (16)
            0 0; ... % (17)
            0 0; ... % (18)
            0 -1000; ... % (19)
            0 0; ... % (20)
            0 0]; % (21) [N]
        systemParameter.NodeElement = [...
            21 1; ... % (1) -
            21 2; ... % (2) \ 
            20 2; ... % (3) -
            17 19; ... % (4) \
            18 19; ... % (5) /
            1 2; ... % (6) |
            1 4; ... % (7) \
            1 3; ... % (8) -
            2 4; ... % (9) -
            3 4; ... % (10) | 
            3 6; ... % (11) \
            3 5; ... % (12) -
            4 6; ... % (13) -
            5 6; ... % (14) |
            5 8; ... % (15) \ 
            5 7; ... % (16) -
            6 8; ... % (17) -
            7 8; ... % (18) |
            7 10; ... % (19) \ 
            7 9; ... % (20) - 
            8 10; ... % (21) -
            9 10; ... % (22) |
            9 12; ... % (23) \ 
            9 11; ... % (24) -
            10 12; ... % (25) - 
            11 12; ... % (26) |
            11 14; ... % (27) \ 
            11 13; ... % (28) -
            12 14; ... % (29) - 
            13 14; ... % (30) |
            13 16; ... % (31) \
            13 15; ... % (32) -
            14 16; ... % (33) -
            15 16; ... % (34) |
            15 18; ... % (35) \ 
            15 17; ... % (36) -
            16 18; ... % (37) - 
            17 18; ... % (38) |
            ];
    case '11-Bar-3D'
        nSample = 300;
        maxIterDisplacement = 500;
        maxIterMass = 100;
        maxIterDisplacementAndMass = 100;
        growthRateDisplacement = 0.04;
        growthRateMass = 0.04;
        trimmingPasses = 'reduced';
        requirementSpacesType = 'Omega1';
        systemParameter.NodePosition = [...
              0   0   0; % (1)
              0 0.5   1; % (2)
              0   1   0; % (3)
              1 0 0; % (4)
              1 0.5 1; % (5) 
              1 1 1; % (6)
              2 0.5 0.5]; % (7)
        systemParameter.FixedDegreesOfFreedom = [...
            ...
            true true true; % (1) 
            true true true; % (2)
            true true true; % (3)
            false false false; % (4)
            false false false; % (5)
            false false false; % (6)
            false false false]; % (7)
        systemParameter.NodeForce = [...
            0 0     0; % (1)
            0 0     0; % (2)
            0 0     0; % (3)
            0 0     0; % (4)
            0 0     0; % (5)
            0 0     0; % (6)
            0 0 -1000]; % (7)
        systemParameter.NodeElement = [...
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
    case '21-Bar-3D'
        nSample = 300;
        maxIterDisplacement = 500;
        maxIterMass = 250;
        maxIterDisplacementAndMass = 300;
        growthRateDisplacement = 0.005;
        growthRateMass = 0.005;
        trimmingPasses = 'single';
        requirementSpacesType = 'Omega1';
        systemParameter.NodePosition = [...
              0   0   0; % (1)
              0 0.5   1; % (2)
              0   1   0; % (3)
              1    0  0; % (4)
              1  0.5  1; % (5) 
              1    1  0; % (6)
              2    0  0; % (7)
              2  0.5  1; % (8) 
              2    1  0; % (9)
              3 0.5 0.5]; % (10)
        systemParameter.FixedDegreesOfFreedom = [...
            true true true; % (1) 
            true true true; % (2)
            true true true; % (3)
            false false false; % (4)
            false false false; % (5)
            false false false; % (6)
            false false false; % (7)
            false false false; % (8)
            false false false; % (9)
            false false false]; % (10)
        systemParameter.NodeForce = [...
            0 0     0; % (1)
            0 0     0; % (2)
            0 0     0; % (3)
            0 0     0; % (4)
            0 0     0; % (5)
            0 0     0; % (6)
            0 0     0; % (7)
            0 0     0; % (8)
            0 0     0; % (9)
            0 0 -1000]; % (10)
        systemParameter.NodeElement = [...
            ... % 
            1 4; % (1) -
            2 5; % (2) -
            3 6; % (3) -
            1 5; % (4) /
            3 4; % (5) \
            3 5; % (6) /
            4 5; % (7) |
            5 6; % (8) \
            4 6; % (9) /
            ... %
            4 7;  % (10) -
            5 8;  % (11) -
            6 9;  % (12) -
            4 8;  % (13) /
            6 7;  % (14) \
            6 8;  % (15) /
            7 8;  % (16) |
            8 9;  % (17) \
            7 9;  % (18) /
            ... %
            7 10; % (19)
            8 10; % (20)
            9 10]; % (21)
    case '39-Bar-3D'
        nSample = 300;
        maxIterDisplacement = 500;
        maxIterMass = 300;
        maxIterDisplacementAndMass = 300;
        growthRateDisplacement = 0.004;
        growthRateMass = 0.007;
        trimmingPasses = 'single';
        requirementSpacesType = 'Omega1';
        systemParameter.NodePosition = [...
              0   0   0; % (1)
              0 0.5   1; % (2)
              0   1   0; % (3)
              1    0  0; % (4)
              1  0.5  1; % (5) 
              1    1  0; % (6)
              2    0  0; % (7)
              2  0.5  1; % (8) 
              2    1  0; % (9)
              3    0  0; % (10)
              3  0.5  1; % (11) 
              3    1  0; % (12)
              4    0  0; % (13)
              4  0.5  1; % (14) 
              4    1  0; % (15)
              5 0.5 0.5]; % (16)
        systemParameter.FixedDegreesOfFreedom = [...
            true true true; % (1) 
            true true true; % (2)
            true true true; % (3)
            false false false; % (4)
            false false false; % (5)
            false false false; % (6)
            false false false; % (7)
            false false false; % (8)
            false false false; % (9)
            false false false; % (10)
            false false false; % (11)
            false false false; % (12)
            false false false; % (13)
            false false false; % (14)
            false false false; % (15)
            false false false]; % (16)
        systemParameter.NodeForce = [...
            0 0     0; % (1)
            0 0     0; % (2)
            0 0     0; % (3)
            0 0     0; % (4)
            0 0     0; % (5)
            0 0     0; % (6)
            0 0     0; % (7)
            0 0     0; % (8)
            0 0     0; % (9)
            0 0     0; % (10)
            0 0     0; % (11)
            0 0     0; % (12)
            0 0     0; % (13)
            0 0     0; % (14)
            0 0     0; % (15)
            0 0 -1000]; % (16)
        systemParameter.NodeElement = [...
            ... % 
            1 4; % (1) -
            2 5; % (2) -
            3 6; % (3) -
            1 5; % (4) /
            3 4; % (5) \
            3 5; % (6) /
            4 5; % (7) |
            5 6; % (8) \
            4 6; % (9) /
            ... %
            4 7;  % (10) -
            5 8;  % (11) -
            6 9;  % (12) -
            4 8;  % (13) /
            6 7;  % (14) \
            6 8;  % (15) /
            7 8;  % (16) |
            8 9;  % (17) \
            7 9;  % (18) /
            ... %
            7 10; % (19) -
            8 11; % (20) -
            9 12; % (21) -
            7 11; % (22) /
            9 10; % (23) \
            9 11; % (24) /
            10 11; % (25) |
            11 12; % (26) \
            10 12; % (27) /
            ... %
            10 13; % (28) -
            11 14; % (29) -
            12 15; % (30) -
            10 14; % (31) /
            12 13; % (32) \
            12 14; % (33) /
            13 14; % (34) |
            14 15; % (35) \
            13 15; % (36) /
            ... %
            13 16; % (37)
            14 16; % (38)
            15 16]; % (39)
end


%% base parameters
nElement = size(systemParameter.NodeElement,1);

minimumRadius = 0.5;
minimumThickness = 0.001;
minimumYoungsModulus = min(youngsModulusMaterial);

maximumRadius = 1.5;
maximumThickness = 0.5;
maximumYoungsModulus = max(youngsModulusMaterial);

designSpaceLowerBound = repmat([minimumRadius,minimumThickness,minimumYoungsModulus],1,nElement);
designSpaceUpperBound = repmat([maximumRadius,maximumThickness,maximumYoungsModulus],1,nElement);
initialDesign = repmat([minimumRadius,minimumThickness,minimumYoungsModulus]+[maximumRadius,maximumThickness,maximumYoungsModulus],1,nElement)/2;

componentIndex = arrayfun(@(i) (3*(i-1)+1):3*i, 1:nElement, 'UniformOutput', false);

systemParameter.EstimateMassGivenYoungsModulus = @(E) interp1(youngsModulusMaterial,densityMaterial,E);
systemParameter.EstimateYieldStrengthGivenYoungsModulus = @(E) interp1(youngsModulusMaterial,yieldStrengthMaterial,E);


%% Analyse initial truss
% deformed initial truss
elementRadius = initialDesign(1:3:end)';
elementThickness = initialDesign(2:3:end)';
elementYoungsModulus = initialDesign(3:3:end)';
elementCrossSectionArea = pi*elementRadius.^2 - pi*(elementRadius-elementThickness).^2;

nodeDisplacementInitial = ...
    truss_analysis(...
	    systemParameter.NodePosition,...
	    systemParameter.FixedDegreesOfFreedom,...
	    systemParameter.NodeForce,...
	    systemParameter.NodeElement,...
	    elementCrossSectionArea,...
	    elementYoungsModulus);

is3dPlot = (size(systemParameter.NodePosition,2) == 3);

figure;

plot_truss_deformation(gcf,systemParameter.NodePosition,systemParameter.NodeElement,'TrussPlotOptions',{'Color',color_palette_tol('blue')},'MaximumLinewidth',3.0,'ShowBarNumber',true,'BarNumberOptions',{'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',16});

isTrussTip = (systemParameter.NodeForce~=0);
positionTip = systemParameter.NodePosition(any(isTrussTip,2),:);
appliedForceOptions = {'Color',color_palette_tol('red'),'LineWidth',3.0};
if(is3dPlot)
    [wallY,wallZ] = meshgrid([-0.5,1.5]);
    wallX = zeros(size(wallY));
    wallOptions = {'FaceColor','w','EdgeColor','k','FaceAlpha',0.7};
    handleWall = surf(wallX,wallY,wallZ,wallOptions{:});

    handleForce = quiver3(positionTip(1),positionTip(2),positionTip(3),0,0,-0.5,appliedForceOptions{:});

    set(gca,'XColor', 'none','YColor','none','ZColor','none');

    axis('tight','equal','vis3d');
    camproj('perspective');
    cameratoolbar; % better adjust angle/perspective
else
    wallX = [0 0];
    wallY = [-0.1 1.1];
    wallOptions = {'Linewidth',8.0,'Color','k'};
    handleWall = plot(wallX,wallY,wallOptions{:});

    handleForce = quiver(positionTip(1),positionTip(2),0,-0.5,appliedForceOptions{:});

    set(gca,'XColor', 'none','YColor','none');

    axis('tight','equal');
end

save_print_figure(gcf,[saveFolder,'InitialTrussDeformation'],'PrintFormat',{'png','pdf'},'Size',figureSize);


%% establish upper performance limits
nElement = size(systemParameter.NodeElement,1);
performanceLowerLimit = [  0   0 -inf(1,2*nElement)];
performanceUpperLimit = [nan nan ones(1,nElement) ones(1,nElement)];

% update uppwer limit based on either optimal value or initial value
bottomUpMapping = BottomUpMappingFunction(systemFunction,'SystemParameter',systemParameter);
if(fixRadius)
    [bottomUpMapping,designSpaceLowerBound,designSpaceUpperBound,initialDesign,componentIndex] = ...
        BottomUpMappingFixedVariables(bottomUpMapping,repmat([true false false],1,nElement),initialDesign(1:3:end),designSpaceLowerBound,designSpaceUpperBound,initialDesign,componentIndex);
end
performanceMeasureInitial = bottomUpMapping.response(initialDesign);


%% Solve solution spaces problems
% displacement
if(computeDisplacement)
    performanceUpperLimitDisplacement = performanceUpperLimit;
    performanceUpperLimitDisplacement(1) = performanceMeasureInitial(1)*1.1;

    performanceLowerLimitDisplacement = performanceLowerLimit;
    performanceLowerLimitDisplacement(1) = performanceMeasureInitial(1)*0.9;

    designEvaluatorDisplacement = DesignEvaluatorBottomUpMapping(...
        bottomUpMapping,...
        performanceLowerLimitDisplacement,...
        performanceUpperLimitDisplacement);

    [solutionSpaceBoxDisplacement,componentSolutionSpacePlanarTrimmingDisplacement,componentSolutionSpaceCornerBoxRemovalDisplacement,componentSolutionSpaceHolePunchingDisplacement] = ...
    compute_truss_solution_spaces('Displacement',designEvaluatorDisplacement,initialDesign,designSpaceLowerBound,designSpaceUpperBound,...
        componentIndex,nSample,maxIterDisplacement,growthRateDisplacement,trimmingPasses,requirementSpacesType,...
        useBoxResultForComponent,computeBoxSolution,computePlanarTrimmingComponent,computeCornerBoxRemovalComponent,computeHolePunchingComponent,...
        rngState,saveFolder,figureSize);

    plot_relevant_results_truss_moving_node('Displacement',...
        solutionSpaceBoxDisplacement,componentSolutionSpacePlanarTrimmingDisplacement,componentSolutionSpaceCornerBoxRemovalDisplacement,componentSolutionSpaceHolePunchingDisplacement,...
        saveFolder,figureSize);
end

% mass
if(computeMass)
    performanceUpperLimitMass = performanceUpperLimit;
    performanceUpperLimitMass(2) = performanceMeasureInitial(2);
    designEvaluatorMass = DesignEvaluatorBottomUpMapping(...
        bottomUpMapping,...
        performanceLowerLimit,...
        performanceUpperLimitMass);

    [solutionSpaceBoxMass,componentSolutionSpacePlanarTrimmingMass,componentSolutionSpaceCornerBoxRemovalMass,componentSolutionSpaceHolePunchingMass] = ...
        compute_truss_solution_spaces('Mass',designEvaluatorMass,initialDesign,designSpaceLowerBound,designSpaceUpperBound,...
            componentIndex,nSample,maxIterMass,growthRateMass,trimmingPasses,requirementSpacesType,useBoxResultForComponent,...
            computeBoxSolution,computePlanarTrimmingComponent,computeCornerBoxRemovalComponent,computeHolePunchingComponent,...
            rngState,saveFolder,figureSize);

    plot_relevant_results_truss_moving_node('Mass',...
        solutionSpaceBoxMass,componentSolutionSpacePlanarTrimmingMass,componentSolutionSpaceCornerBoxRemovalMass,componentSolutionSpaceHolePunchingMass,...
        saveFolder,figureSize);
end

% displacement + mass
if(computeDisplacementAndMass)
    performanceUpperLimitDisplacementAndMass = performanceUpperLimit;
    performanceUpperLimitDisplacementAndMass(1) = performanceMeasureInitial(1);
    performanceUpperLimitDisplacementAndMass(2) = performanceMeasureInitial(2);
    designEvaluatorDisplacementAndMass = DesignEvaluatorBottomUpMapping(...
        bottomUpMapping,...
        performanceLowerLimit,...
        performanceUpperLimitDisplacementAndMass);

    [solutionSpaceBoxDisplacementAndMass,componentSolutionSpacePlanarTrimmingDisplacementAndMass,componentSolutionSpaceCornerBoxRemovalDisplacementAndMass,componentSolutionSpaceHolePunchingDisplacementAndMass] = ...
    compute_truss_solution_spaces('DisplacementAndMass',designEvaluatorDisplacementAndMass,initialDesign,...
        designSpaceLowerBound,designSpaceUpperBound,componentIndex,nSample,maxIterDisplacementAndMass,...
        growthRateDisplacement,trimmingPasses,requirementSpacesType,useBoxResultForComponent,...
        computeBoxSolution,computePlanarTrimmingComponent,computeCornerBoxRemovalComponent,computeHolePunchingComponent,...
        rngState,saveFolder,figureSize);

    plot_relevant_results_truss_moving_node('DisplacementAndMass',...
        solutionSpaceBoxDisplacementAndMass,componentSolutionSpacePlanarTrimmingDisplacementAndMass,componentSolutionSpaceCornerBoxRemovalDisplacementAndMass,...
        componentSolutionSpaceHolePunchingDisplacementAndMass,saveFolder,figureSize);
end


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;


%% Subfunctions
function [solutionSpaceBox,componentSolutionSpacePlanarTrimming,componentSolutionSpaceCornerBoxRemoval,componentSolutionSpaceHolePunching] = ...
    compute_truss_solution_spaces(typeName,designEvaluator,initialDesign,designSpaceLowerBound,designSpaceUpperBound,componentIndex,...
        nSample,maxIter,growthRate,trimmingPasses,requirementSpacesType,useBoxResultForComponent,...
        computeBoxSolution,computePlanarTrimmingComponent,computeCornerBoxRemovalComponent,computeHolePunchingComponent,...
        rngState,saveFolder,figureSize)
    
    %% box
    solutionSpaceBox = [];
    if(computeBoxSolution)
        timeElapsedBox = tic;
        optionsBox = sso_stochastic_options('box',...
            'RequirementSpacesType',requirementSpacesType,...
            'NumberSamplesPerIterationExploration',nSample,...
            'NumberSamplesPerIterationConsolidation',nSample,...
            'FixIterNumberExploration',true,...
            'FixIterNumberConsolidation',true,...
            'MaxIterExploration',maxIter,...
            'MaxIterConsolidation',maxIter,...
            'UseAdaptiveGrowthRate',true,...
            'GrowthRate',growthRate,...
            'MaximumGrowthAdaptationFactor',1.5,...
            'ApplyLeanness','never',...
            'TrimmingOperationOptions',{'PassesCriterion',trimmingPasses},...
            'TrimmingOrderOptions',{'OrderPreference','score'},...
            'LoggingLevel','all');
        
        rng(rngState);
        [solutionSpaceBox,optimizationDataBox] = sso_box_stochastic(designEvaluator,...
            initialDesign,designSpaceLowerBound,designSpaceUpperBound,optionsBox);
        toc(timeElapsedBox)

        resultsFolder = [saveFolder,sprintf('PerformanceBox%s/',typeName)];
        mkdir(resultsFolder);
        algoDataBox = postprocess_sso_box_stochastic(optimizationDataBox);
        plot_sso_box_stochastic_metrics(algoDataBox,...
            'SaveFolder',resultsFolder,...
            'CloseFigureAfterSaving',true,...
            'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});
    end

    if(useBoxResultForComponent && ~isempty(solutionSpaceBox))
        initialDesignComponent = solutionSpaceBox;
    else
        if(useBoxResultForComponent && isempty(solutionSpaceBox))
            warning('Box solution requested for component initialization but box solution was not computed. Using initial design instead.');
        end
        initialDesignComponent = initialDesign;
    end
    comparisonComponent = {};
    componentLabel = {};

    %% component planar trimming
    componentSolutionSpacePlanarTrimming = [];
    if(computePlanarTrimmingComponent)
        timeElapsedComponent = tic;
        optionsComponent = sso_stochastic_options('component',...
            'RequirementSpacesType',requirementSpacesType,...
            'NumberSamplesPerIterationExploration',nSample,...
            'NumberSamplesPerIterationConsolidation',nSample,...
            'FixIterNumberExploration',true,...
            'FixIterNumberConsolidation',true,...
            'MaxIterExploration',maxIter,...
            'MaxIterConsolidation',maxIter,...
            'CandidateSpaceConstructor',@CandidateSpaceConvexHull,...
            'CandidateSpaceOptions',{'NormalizeVariables',true,'NormalizeGrowthDirection',true},...
            'TrimmingMethodFunction',@component_trimming_method_planar_trimming,...
            'TrimmingMethodOptions',{'NormalizeVariables',true},...
            'UseAdaptiveGrowthRate',true,...
            'GrowthRate',growthRate,...
            'MaximumGrowthAdaptationFactor',1.5,...
            'ApplyLeanness','never',...
            'MaximumNumberPaddingSamples',10000,...
            'UsePaddingSamplesInTrimming',true,...
            'UsePreviousEvaluatedSamplesConsolidation',false,...
            'TrimmingOperationOptions',{'PassesCriterion',trimmingPasses},...
            'TrimmingOrderOptions',{'OrderPreference','score'},...
            'ShapeSamplesUsefulConsolidation',true,...
            'LoggingLevel','all');
        
        rng(rngState);
        [componentSolutionSpacePlanarTrimming,optimizationDataComponentPlanarTrimming] = sso_component_stochastic(designEvaluator,...
            initialDesignComponent,designSpaceLowerBound,designSpaceUpperBound,componentIndex,optionsComponent);
        toc(timeElapsedComponent)

        resultsFolder = [saveFolder,sprintf('PerformancePlanarTrimming%s/',typeName)];
        mkdir(resultsFolder);
        algoDataComponentPlanarTrimming = postprocess_sso_component_stochastic(optimizationDataComponentPlanarTrimming);
        plot_sso_component_stochastic_metrics(algoDataComponentPlanarTrimming,...
            'SaveFolder',resultsFolder,...
            'CloseFigureAfterSaving',true,...
            'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});

        comparisonComponent{end+1} = algoDataComponentPlanarTrimming;
        componentLabel{end+1} = 'Planar Trimming';
    end

    %% component corner box removal
    componentSolutionSpaceCornerBoxRemoval = [];
    if(computeCornerBoxRemovalComponent)
        timeElapsedComponent = tic;
        optionsComponent = sso_stochastic_options('component',...
            'RequirementSpacesType',requirementSpacesType,...
            'NumberSamplesPerIterationExploration',nSample,...
            'NumberSamplesPerIterationConsolidation',nSample,...
            'FixIterNumberExploration',true,...
            'FixIterNumberConsolidation',true,...
            'MaxIterExploration',maxIter,...
            'MaxIterConsolidation',maxIter,...
            'CandidateSpaceConstructor',@CandidateSpaceCornerBoxRemoval,...
            'CandidateSpaceOptions',{'NormalizeGrowthDirection',true,'CheckRedundantTrimmingGrowth',false,'CheckRedundantTrimmingUpdate',true,'CheckDuplicatePointsGrowth',false,'CheckDuplicatePointsUpdate',true,'MeasureEstimationFactor',2},...
            'TrimmingMethodFunction',@component_trimming_method_corner_box_removal,...
            'TrimmingMethodOptions',{'NormalizeVariables',true},...
            'UseAdaptiveGrowthRate',true,...
            'GrowthRate',growthRate,...
            'MaximumGrowthAdaptationFactor',1.5,...
            'ApplyLeanness','never',...
            'MaximumNumberPaddingSamples',10000,...
            'UsePaddingSamplesInTrimming',true,...
            'UsePreviousEvaluatedSamplesConsolidation',false,...
            'TrimmingOperationOptions',{'PassesCriterion',trimmingPasses},...
            'TrimmingOrderOptions',{'OrderPreference','score'},...
            'ShapeSamplesUsefulConsolidation',true,...
            'LoggingLevel','all');
        
        rng(rngState);
        [componentSolutionSpaceCornerBoxRemoval,optimizationDataComponentCornerBoxRemoval] = sso_component_stochastic(designEvaluator,...
            initialDesignComponent,designSpaceLowerBound,designSpaceUpperBound,componentIndex,optionsComponent);
        toc(timeElapsedComponent)

        resultsFolder = [saveFolder,sprintf('PerformanceCornerBoxRemoval%s/',typeName)];
        mkdir(resultsFolder);
        algoDataComponentCornerBoxRemoval = postprocess_sso_component_stochastic(optimizationDataComponentCornerBoxRemoval);
        plot_sso_component_stochastic_metrics(algoDataComponentCornerBoxRemoval,...
            'SaveFolder',resultsFolder,...
            'CloseFigureAfterSaving',true,...
            'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});

        comparisonComponent{end+1} = algoDataComponentCornerBoxRemoval;
        componentLabel{end+1} = 'Corner Box Removal';
    end

    %% component hole punching
    componentSolutionSpaceHolePunching = [];
    if(computeHolePunchingComponent)
        timeElapsedComponent = tic;
        optionsComponent = sso_stochastic_options('component',...
            'RequirementSpacesType',requirementSpacesType,...
            'NumberSamplesPerIterationExploration',nSample,...
            'NumberSamplesPerIterationConsolidation',nSample,...
            'FixIterNumberExploration',true,...
            'FixIterNumberConsolidation',true,...
            'MaxIterExploration',maxIter,...
            'MaxIterConsolidation',maxIter,...
            'CandidateSpaceConstructor',@CandidateSpaceSvm,...
            'TrimmingMethodFunction',@component_trimming_method_hole_punching,...
            'UseAdaptiveGrowthRate',true,...
            'GrowthRate',growthRate,...
            'MaximumGrowthAdaptationFactor',1.5,...
            'ApplyLeanness','never',...
            'MaximumNumberPaddingSamples',10000,...
            'UsePaddingSamplesInTrimming',true,...
            'UsePreviousEvaluatedSamplesConsolidation',true,...
            'TrimmingOperationOptions',{'PassesCriterion',trimmingPasses},...
            'TrimmingOrderOptions',{'OrderPreference','score'},...
            'ShapeSamplesUsefulConsolidation',false,...
            'NumberPaddingSamples',1000,...
            'LoggingLevel','all');
        
        rng(rngState);
        [componentSolutionSpaceHolePunching,optimizationDataComponentHolePunching] = sso_component_stochastic(designEvaluator,...
            initialDesignComponent,designSpaceLowerBound,designSpaceUpperBound,componentIndex,optionsComponent);
        toc(timeElapsedComponent)

        resultsFolder = [saveFolder,sprintf('PerformanceHolePunching%s/',typeName)];
        mkdir(resultsFolder);
        algoDataComponentHolePunching = postprocess_sso_component_stochastic(optimizationDataComponentHolePunching);
        plot_sso_component_stochastic_metrics(algoDataComponentHolePunching,...
            'SaveFolder',resultsFolder,...
            'CloseFigureAfterSaving',true,...
            'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});

        comparisonComponent{end+1} = algoDataComponentHolePunching;
        componentLabel{end+1} = 'Hole Punching';
    end

    % comparison
    resultsFolder = [saveFolder,sprintf('PerformanceComparison%s/',typeName)];
    mkdir(resultsFolder);
    
    boxData = {};
    if(computeBoxSolution)
        boxData = {algoDataBox};
    end
    
    plot_sso_comparison_box_component_stochastic_metrics(...
        boxData,...
        comparisonComponent,...
        'ComponentLabel',componentLabel,...
        'BoxColor','k',...
        'ComponentColor',color_palette_tol({'purple','cyan','yellow'}),...
        'SaveFolder',resultsFolder,...
        'CloseFigureAfterSaving',true,...
        'SaveFigureOptions',{'Size',figureSize,'PrintFormat',{'png','pdf'}});
end

function plot_relevant_results_truss_moving_node(typeName,...
    solutionSpaceBox,componentSolutionSpacePlanarTrimming,componentSolutionSpaceCornerBoxRemoval,componentSolutionSpaceHolePunching,...
    saveFolder,figureSize)

    resultsFolder = [saveFolder,sprintf('TrussResult%s/',typeName)];
    mkdir(resultsFolder);
    
    % box solution space + component solution spaces
    figureHandle = plot_results_truss_generic_moving_node(...
        sprintf('BoxSolutionSpace%s',typeName),solutionSpaceBox,...
        sprintf('PlanarTrimmingSolutionSpace%s',typeName),componentSolutionSpacePlanarTrimming,...
        sprintf('CornerBoxRemovalSolutionSpace%s',typeName),componentSolutionSpaceCornerBoxRemoval,...
        sprintf('HolePunchingSolutionSpace%s',typeName),componentSolutionSpaceHolePunching);

    for i=1:length(figureHandle)
        figureName = sprintf('InitialTrussBoxComponent%0*d',get_number_digits_integer(length(figureHandle)),i);
        save_print_figure(figureHandle(i),[resultsFolder,figureName],'PrintFormat',{'png','pdf'},'Size',figureSize);
    end
end

function figureHandle = plot_results_truss_generic_moving_node(varargin)
    parser = inputParser;
    % (1) box - displacement 
    parser.addParameter('BoxSolutionSpaceDisplacement',[]);
    parser.addParameter('BoxSolutionSpaceDisplacementOptions',{});
    % (2) box - mass 
    parser.addParameter('BoxSolutionSpaceMass',[]);
    parser.addParameter('BoxSolutionSpaceMassOptions',{});
    % (3) box - displacement and mass
    parser.addParameter('BoxSolutionSpaceDisplacementAndMass',[]);
    parser.addParameter('BoxSolutionSpaceDisplacementAndMassOptions',{});
    % (4) planar trimming - displacement
    parser.addParameter('PlanarTrimmingSolutionSpaceDisplacement',[]);
    parser.addParameter('PlanarTrimmingSolutionSpaceDisplacementOptions',{});
    % (5) planar trimming - mass
    parser.addParameter('PlanarTrimmingSolutionSpaceMass',[]);
    parser.addParameter('PlanarTrimmingSolutionSpaceMassOptions',{});
    % (6) planar trimming - displacement and mass
    parser.addParameter('PlanarTrimmingSolutionSpaceDisplacementAndMass',[]);
    parser.addParameter('PlanarTrimmingSolutionSpaceDisplacementAndMassOptions',{});
    % (7) corner box removal - displacement
    parser.addParameter('CornerBoxRemovalSolutionSpaceDisplacement',[]);
    parser.addParameter('CornerBoxRemovalSolutionSpaceDisplacementOptions',{});
    % (8) corner box removal - mass
    parser.addParameter('CornerBoxRemovalSolutionSpaceMass',[]);
    parser.addParameter('CornerBoxRemovalSolutionSpaceMassOptions',{});
    % (9) corner box removal - displacement and mass
    parser.addParameter('CornerBoxRemovalSolutionSpaceDisplacementAndMass',[]);
    parser.addParameter('CornerBoxRemovalSolutionSpaceDisplacementAndMassOptions',{});
    % (10) hole punching - displacement
    parser.addParameter('HolePunchingSolutionSpaceDisplacement',[]);
    parser.addParameter('HolePunchingSolutionSpaceDisplacementOptions',{});
    % (11) hole punching - mass
    parser.addParameter('HolePunchingSolutionSpaceMass',[]);
    parser.addParameter('HolePunchingSolutionSpaceMassOptions',{});
    % (12) hole punching - displacement and mass
    parser.addParameter('HolePunchingSolutionSpaceDisplacementAndMass',[]);
    parser.addParameter('HolePunchingSolutionSpaceDisplacementAndMassOptions',{});
    % (13) axes
    parser.addParameter('IncludeAxesInformation',true);
    % (14) legend
    parser.addParameter('IncludeLegend',false);
    parser.addParameter('LegendOptions',{});
    % (15) title
    parser.addParameter('ShowTitle',true);
    parser.addParameter('TitleOptions',{});

    parser.parse(varargin{:});
    options = parser.Results;

    nElement(1) = length(options.PlanarTrimmingSolutionSpaceDisplacement);
    nElement(2) = length(options.PlanarTrimmingSolutionSpaceMass);
    nElement(3) = length(options.PlanarTrimmingSolutionSpaceDisplacementAndMass);
    nElement(4) = length(options.CornerBoxRemovalSolutionSpaceDisplacement);
    nElement(5) = length(options.CornerBoxRemovalSolutionSpaceMass);
    nElement(6) = length(options.CornerBoxRemovalSolutionSpaceDisplacementAndMass);
    nElement(7) = length(options.HolePunchingSolutionSpaceDisplacement);
    nElement(8) = length(options.HolePunchingSolutionSpaceMass);
    nElement(9) = length(options.HolePunchingSolutionSpaceDisplacementAndMass);
    [nElement,iEntry] = max(nElement);

    switch iEntry
        case 1
            nDimension = size(options.PlanarTrimmingSolutionSpaceDisplacement(1).DesignSpaceLowerBound,2);
        case 2
            nDimension = size(options.PlanarTrimmingSolutionSpaceMass(1).DesignSpaceLowerBound,2);
        case 3
            nDimension = size(options.PlanarTrimmingSolutionSpaceDisplacementAndMass(1).DesignSpaceLowerBound,2);
        case 4
            nDimension = size(options.CornerBoxRemovalSolutionSpaceDisplacement(1).DesignSpaceLowerBound,2);
        case 5
            nDimension = size(options.CornerBoxRemovalSolutionSpaceMass(1).DesignSpaceLowerBound,2);
        case 6
            nDimension = size(options.CornerBoxRemovalSolutionSpaceDisplacementAndMass(1).DesignSpaceLowerBound,2);
        case 7
            nDimension = size(options.HolePunchingSolutionSpaceDisplacement(1).DesignSpaceLowerBound,2);
        case 8
            nDimension = size(options.HolePunchingSolutionSpaceMass(1).DesignSpaceLowerBound,2);
        case 9
            nDimension = size(options.HolePunchingSolutionSpaceDisplacementAndMass(1).DesignSpaceLowerBound,2);
    end
    if(nDimension==2)
        boxPlotFunction = @plot_design_box_2d;
    elseif(nDimension==3)
        boxPlotFunction = @plot_design_box_3d;
    end

    %% common options
    % (1) box - displacement 
    defaultCommonBoxSolutionDisplacementOptions = {'EdgeColor','k','FaceAlpha',0.0};
    % (2) box - mass 
    defaultCommonBoxSolutionMassOptions = {'EdgeColor','k','FaceAlpha',0.0};
    % (3) box - displacement and mass
    defaultCommonBoxSolutionDisplacementAndMassOptions = {'EdgeColor','k','FaceAlpha',0.0};
    % (4) planar trimming - displacement
    defaultCommonPlanarTrimmingDisplacementOptions = {'FaceColor',color_palette_tol('purple'),'EdgeColor',color_palette_tol('purple'),'FaceAlpha',0.25};
    % (5) planar trimming - mass
    defaultCommonPlanarTrimmingMassOptions = {'FaceColor',color_palette_tol('purple'),'EdgeColor',color_palette_tol('purple'),'FaceAlpha',0.25};
    % (6) planar trimming - displacement and mass
    defaultCommonPlanarTrimmingDisplacementAndMassOptions = {'FaceColor',color_palette_tol('purple'),'EdgeColor',color_palette_tol('purple'),'FaceAlpha',0.25};
    % (7) corner box removal - displacement
    defaultCommonCornerBoxRemovalDisplacementOptions = {'FaceColor',color_palette_tol('cyan'),'EdgeColor',color_palette_tol('cyan'),'FaceAlpha',0.4,'EdgeColor','none'};
    % (8) corner box removal - mass
    defaultCommonCornerBoxRemovalMassOptions = {'FaceColor',color_palette_tol('cyan'),'EdgeColor',color_palette_tol('cyan'),'FaceAlpha',0.4,'EdgeColor','none'};
    % (9) corner box removal - displacement and mass
    defaultCommonCornerBoxRemovalDisplacementAndMassOptions = {'FaceColor',color_palette_tol('cyan'),'EdgeColor',color_palette_tol('cyan'),'FaceAlpha',0.4,'EdgeColor','none'};
    % (10) hole punching - displacement
    defaultCommonHolePunchingDisplacementOptions = {'FaceColor',color_palette_tol('yellow'),'EdgeColor',color_palette_tol('yellow'),'FaceAlpha',0.3};
    % (11) hole punching - mass
    defaultCommonHolePunchingMassOptions = {'FaceColor',color_palette_tol('yellow'),'EdgeColor',color_palette_tol('yellow'),'FaceAlpha',0.3};
    % (12) hole punching - displacement and mass
    defaultCommonHolePunchingDisplacementAndMassOptions = {'FaceColor',color_palette_tol('yellow'),'EdgeColor',color_palette_tol('yellow'),'FaceAlpha',0.3};
    % (13) axes
    % (14) legend
    defaultCommonLegendOptions = {'location','west'};
    % (15) title
    defaultCommonTitleOptions = {'FontSize',14};

    %% merge options
    % (1) box - displacement 
    [~,boxDisplacementOptions] = merge_name_value_pair_argument(...
        defaultCommonBoxSolutionDisplacementOptions,...
        options.BoxSolutionSpaceDisplacementOptions);
    % (2) box - mass 
    [~,boxMassOptions] = merge_name_value_pair_argument(...
        defaultCommonBoxSolutionMassOptions,...
        options.BoxSolutionSpaceMassOptions);
    % (3) box - displacement and mass
    [~,boxDisplacementAndMassOptions] = merge_name_value_pair_argument(...
        defaultCommonBoxSolutionDisplacementAndMassOptions,...
        options.BoxSolutionSpaceDisplacementAndMassOptions);
    % (4) planar trimming - displacement
    [~,planarTrimmingDisplacementOptions] = merge_name_value_pair_argument(...
        defaultCommonPlanarTrimmingDisplacementOptions,...
        options.PlanarTrimmingSolutionSpaceDisplacementOptions);
    % (5) planar trimming - mass
    [~,planarTrimmingMassOptions] = merge_name_value_pair_argument(...
        defaultCommonPlanarTrimmingMassOptions,...
        options.PlanarTrimmingSolutionSpaceMassOptions);
    % (6) planar trimming - displacement and mass
    [~,planarTrimmingDisplacementAndMassOptions] = merge_name_value_pair_argument(...
        defaultCommonPlanarTrimmingDisplacementAndMassOptions,...
        options.PlanarTrimmingSolutionSpaceDisplacementAndMassOptions);
    % (7) corner box removal - displacement
    [~,cornerBoxRemovalDisplacementOptions] = merge_name_value_pair_argument(...
        defaultCommonCornerBoxRemovalDisplacementOptions,...
        options.CornerBoxRemovalSolutionSpaceDisplacementOptions);
    % (8) corner box removal - mass
    [~,cornerBoxRemovalMassOptions] = merge_name_value_pair_argument(...
        defaultCommonCornerBoxRemovalMassOptions,...
        options.CornerBoxRemovalSolutionSpaceMassOptions);
    % (9) corner box removal - displacement and mass
    [~,cornerBoxRemovalDisplacementAndMassOptions] = merge_name_value_pair_argument(...
        defaultCommonCornerBoxRemovalDisplacementAndMassOptions,...
        options.CornerBoxRemovalSolutionSpaceDisplacementAndMassOptions);
    % (10) hole punching - displacement
    [~,holePunchingDisplacementOptions] = merge_name_value_pair_argument(...
        defaultCommonHolePunchingDisplacementOptions,...
        options.HolePunchingSolutionSpaceDisplacementOptions);
    % (11) hole punching - mass
    [~,holePunchingMassOptions] = merge_name_value_pair_argument(...
        defaultCommonHolePunchingMassOptions,...
        options.HolePunchingSolutionSpaceMassOptions);
    % (12) hole punching - displacement and mass
    [~,holePunchingDisplacementAndMassOptions] = merge_name_value_pair_argument(...
        defaultCommonHolePunchingDisplacementAndMassOptions,...
        options.HolePunchingSolutionSpaceDisplacementAndMassOptions);
    % (13) axes 
    % (14) legend
    [~,legendOptions] = merge_name_value_pair_argument(...
        defaultCommonLegendOptions,...
        options.LegendOptions);
    % (15) title
    [~,titleOptions] = merge_name_value_pair_argument(...
        defaultCommonTitleOptions,...
        options.TitleOptions);


    %% plot each element (where applicable)
    for i=1:nElement
        figureHandle(i) = figure;
        hold all;
    end

    % (1) box - displacement 
    handleDisplacementBox = [];
    if(~isempty(options.BoxSolutionSpaceDisplacement))
        for i=1:nElement
            figure(figureHandle(i));
            currentIndex = 1 + nDimension*(i-1) + [0:(nDimension-1)];
            if i < nElement
                boxPlotFunction(gcf,options.BoxSolutionSpaceDisplacement(:,currentIndex),boxDisplacementOptions{:});
            else
                handleDisplacementBox = boxPlotFunction(gcf,options.BoxSolutionSpaceDisplacement(:,currentIndex),boxDisplacementOptions{:});
            end
        end
    end

    % (2) box - mass 
    handleMassBox = [];
    if(~isempty(options.BoxSolutionSpaceMass))
        for i=1:nElement
            figure(figureHandle(i));
            currentIndex = 1 + nDimension*(i-1) + [0:(nDimension-1)];
            if i < nElement
                boxPlotFunction(gcf,options.BoxSolutionSpaceMass(:,currentIndex),boxMassOptions{:});
            else
                handleMassBox = boxPlotFunction(gcf,options.BoxSolutionSpaceMass(:,currentIndex),boxMassOptions{:});
            end
        end
    end

    % (3) box - displacement and mass
    handleDisplacementAndMassBox = [];
    if(~isempty(options.BoxSolutionSpaceDisplacementAndMass))
        for i=1:nElement
            figure(figureHandle(i));
            currentIndex = 1 + nDimension*(i-1) + [0:(nDimension-1)];
            if i < nElement
                boxPlotFunction(gcf,options.BoxSolutionSpaceDisplacementAndMass(:,currentIndex),boxDisplacementAndMassOptions{:});
            else
                handleDisplacementAndMassBox = boxPlotFunction(gcf,options.BoxSolutionSpaceDisplacementAndMass(:,currentIndex),boxDisplacementAndMassOptions{:});
            end
        end
    end

    % (4) planar trimming - displacement
    handleDisplacementComponentPlanarTrimming = [];
    if(~isempty(options.PlanarTrimmingSolutionSpaceDisplacement))
        for i=1:nElement
            figure(figureHandle(i));
            if i < nElement
                options.PlanarTrimmingSolutionSpaceDisplacement(i).plot_candidate_space(gcf,planarTrimmingDisplacementOptions{:});
            else
                handleDisplacementComponentPlanarTrimming = options.PlanarTrimmingSolutionSpaceDisplacement(i).plot_candidate_space(gcf,planarTrimmingDisplacementOptions{:});
            end
        end
    end

    % (5) planar trimming - mass
    handleMassComponentPlanarTrimming = [];
    if(~isempty(options.PlanarTrimmingSolutionSpaceMass))
        for i=1:nElement
            figure(figureHandle(i));
            if i < nElement
                options.PlanarTrimmingSolutionSpaceMass(i).plot_candidate_space(gcf,planarTrimmingMassOptions{:});
            else
                handleMassComponentPlanarTrimming = options.PlanarTrimmingSolutionSpaceMass(i).plot_candidate_space(gcf,planarTrimmingMassOptions{:});
            end
        end
    end

    % (6) planar trimming - displacement and mass
    handleDisplacementAndMassComponentPlanarTrimming = [];
    if(~isempty(options.PlanarTrimmingSolutionSpaceDisplacementAndMass))
        for i=1:nElement
            figure(figureHandle(i));
            if i < nElement
                options.PlanarTrimmingSolutionSpaceDisplacementAndMass(i).plot_candidate_space(gcf,planarTrimmingDisplacementAndMassOptions{:});
            else
                handleDisplacementAndMassComponentPlanarTrimming = options.PlanarTrimmingSolutionSpaceDisplacementAndMass(i).plot_candidate_space(gcf,planarTrimmingDisplacementAndMassOptions{:});
            end
        end
    end

    % (7) corner box removal - displacement
    handleDisplacementComponentCornerBoxRemoval = [];
    if(~isempty(options.CornerBoxRemovalSolutionSpaceDisplacement))
        for i=1:nElement
            figure(figureHandle(i));
            if i < nElement
                options.CornerBoxRemovalSolutionSpaceDisplacement(i).plot_candidate_space(gcf,cornerBoxRemovalDisplacementOptions{:});
            else
                handleDisplacementComponentCornerBoxRemoval = options.CornerBoxRemovalSolutionSpaceDisplacement(i).plot_candidate_space(gcf,cornerBoxRemovalDisplacementOptions{:});
            end
        end
    end

    % (8) corner box removal - mass
    handleMassComponentCornerBoxRemoval = [];
    if(~isempty(options.CornerBoxRemovalSolutionSpaceMass))
        for i=1:nElement
            figure(figureHandle(i));
            if i < nElement
                options.CornerBoxRemovalSolutionSpaceMass(i).plot_candidate_space(gcf,cornerBoxRemovalMassOptions{:});
            else
                handleMassComponentCornerBoxRemoval = options.CornerBoxRemovalSolutionSpaceMass(i).plot_candidate_space(gcf,cornerBoxRemovalMassOptions{:});
            end
        end
    end

    % (9) corner box removal - displacement and mass
    handleDisplacementAndMassComponentCornerBoxRemoval = [];
    if(~isempty(options.CornerBoxRemovalSolutionSpaceDisplacementAndMass))
        for i=1:nElement
            figure(figureHandle(i));
            if i < nElement
                options.CornerBoxRemovalSolutionSpaceDisplacementAndMass(i).plot_candidate_space(gcf,cornerBoxRemovalDisplacementAndMassOptions{:});
            else
                handleDisplacementAndMassComponentCornerBoxRemoval = options.CornerBoxRemovalSolutionSpaceDisplacementAndMass(i).plot_candidate_space(gcf,cornerBoxRemovalDisplacementAndMassOptions{:});
            end
        end
    end

    % (10) hole punching - displacement
    handleDisplacementComponentHolePunching = [];
    if(~isempty(options.HolePunchingSolutionSpaceDisplacement))
        for i=1:nElement
            figure(figureHandle(i));
            if i < nElement
                options.HolePunchingSolutionSpaceDisplacement(i).plot_candidate_space(gcf,holePunchingDisplacementOptions{:});
            else
                handleDisplacementComponentHolePunching = options.HolePunchingSolutionSpaceDisplacement(i).plot_candidate_space(gcf,holePunchingDisplacementOptions{:});
            end
        end
    end

    % (11) hole punching - mass
    handleMassComponentHolePunching = [];
    if(~isempty(options.HolePunchingSolutionSpaceMass))
        for i=1:nElement
            figure(figureHandle(i));
            if i < nElement
                options.HolePunchingSolutionSpaceMass(i).plot_candidate_space(gcf,holePunchingMassOptions{:});
            else
                handleMassComponentHolePunching = options.HolePunchingSolutionSpaceMass(i).plot_candidate_space(gcf,holePunchingMassOptions{:});
            end
        end
    end

    % (12) hole punching - displacement and mass
    handleDisplacementAndMassComponentHolePunching = [];
    if(~isempty(options.HolePunchingSolutionSpaceDisplacementAndMass))
        for i=1:nElement
            figure(figureHandle(i));
            if i < nElement
                options.HolePunchingSolutionSpaceDisplacementAndMass(i).plot_candidate_space(gcf,holePunchingDisplacementAndMassOptions{:});
            else
                handleDisplacementAndMassComponentHolePunching = options.HolePunchingSolutionSpaceDisplacementAndMass(i).plot_candidate_space(gcf,holePunchingDisplacementAndMassOptions{:});
            end
        end
    end

    % (13) axes
    grid('off');
    if(~options.IncludeAxesInformation)
        if(nDimension==2)
            set(gca,'XColor', 'none','YColor','none');
        elseif(nDimension==3)
            set(gca,'XColor', 'none','YColor','none','ZColor','none');
        end
    end
    for i=1:nElement
        figure(figureHandle(i));
        if(nDimension==2)
            xlabel('Element Thickness [m]');
            ylabel('Young''s Modulus [Pa]');
        elseif(nDimension==3)
            xlabel('Element Radius [m]');
            ylabel('Element Thickness [m]');
            zlabel('Young''s Modulus [Pa]');
            axis('tight','vis3d');
            camproj('perspective');
        end
    end
    

    % (14) legend
    if(options.IncludeLegend)
        handleObjectAll = {...
            handleDisplacementBox,...
            handleDisplacementComponentPlanarTrimming,...
            handleDisplacementComponentCornerBoxRemoval,...
            handleDisplacementComponentHolePunching,...
            handleMassBox,...
            handleMassComponentPlanarTrimming,...
            handleMassComponentCornerBoxRemoval,...
            handleMassComponentHolePunching,...
            handleDisplacementAndMassBox,...
            handleDisplacementAndMassComponentPlanarTrimming,...
            handleDisplacementAndMassComponentCornerBoxRemoval,...
            handleDisplacementAndMassComponentHolePunching};

        legendTextAll = {...
            'Box - Displacement',...
            'Planar Trimming - Displacement',...
            'Corner Box Removal - Displacement',...
            'Hole Punching - Displacement',...
            'Box - Mass',...
            'Planar Trimming - Mass',...
            'Corner Box Removal - Mass',...
            'Hole Punching - Mass',...
            'Box - Displacement + Mass',...
            'Planar Trimming - Displacement + Mass',...
            'Corner Box Removal - Displacement + Mass',...
            'Hole Punching - Displacement + Mass'};
        
        handleObject = [handleObjectAll{:}];
        legentText = {legendTextAll{~cellfun(@isempty,handleObjectAll)}};

        legend(handleObject,legentText,legendOptions{:});
    end

    if(options.ShowTitle)
        for i=1:nElement
            figure(figureHandle(i));
            title(['Element ',num2str(i)],titleOptions{:});
        end
    end

    if(nargout<1)
        clear figureHandle;
    end
end