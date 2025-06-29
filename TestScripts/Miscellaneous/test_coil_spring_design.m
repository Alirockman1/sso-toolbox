%test_coil_spring_design
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

%                      rho    G sigmaY
systemParameter = [7.85e-6 79e3    250];
%                          d   D   n    s
designSpaceLowerBound = [0.2  15   3 0.25];
designSpaceUpperBound = [ 12 320  15 12.6];

performanceLowerLimit = [5.1 -inf -inf -inf   26 -inf    5   50];
performanceUpperLimit = [5.1 +inf   95   38 +inf  -20 +inf +inf];

bottomUpMapping = BottomUpMappingFunction(@coil_spring_design_given_material,systemParameter);
initialDesign = mean([designSpaceLowerBound;designSpaceUpperBound],1);

objectiveFunction = @(x)x(:,2);
equalityConstraintFunction = @(x) x(:,1)-5.1;
inequalityConstraintFunction = @(x) [x(:,3)-95, x(:,4)-38, 26-x(:,5), x(:,6)+20, 5-x(:,7), 50-x(:,8)];

[optimalSpring,massOptimal] = design_optimize_quantities_of_interest(...
    bottomUpMapping,...
    initialDesign,...
    designSpaceLowerBound,...
    designSpaceUpperBound,...
    objectiveFunction,...
    'InequalityConstraintFunction',inequalityConstraintFunction,...
    'EqualityConstraintFunction',equalityConstraintFunction,...
    'OptimizationMethodFunction',@optimization_fmincon_wrapper,...
    'OptimizationMethodOptions',{'Display','iter-detailed'});

quantityOfInterest = bottomUpMapping.response(optimalSpring);
fprintf('\nDesign Variables:\n');
fprintf('                    Wire Diameter (d): %g mm\n',optimalSpring(1));
fprintf('                    Coil Diameter (D): %g mm\n',optimalSpring(2));
fprintf('               Number of Windings (n): %g\n',optimalSpring(3));
fprintf('Spacing of Wires in Relaxed State (s): %g mm\n',optimalSpring(4));
if(length(optimalSpring)>=5)
    fprintf('                        Density (rho): %g kg/mm^3\n',optimalSpring(5));
    fprintf('                  Shear Stiffness (G): %g MPa\n',optimalSpring(6));
    fprintf('                Yield Stress (sigmaY): %g MPa\n',optimalSpring(7));
end
fprintf('\nQuantities of Interest:\n');
fprintf('                         Stiffness (K): %g N/mm\n',quantityOfInterest(1));
fprintf('                              Mass (m): %g kg\n',quantityOfInterest(2));
fprintf('          Height of Operation Space(H): %g mm\n',quantityOfInterest(3));
fprintf('Outer Diameter of Operation Space (Wo): %g mm\n',quantityOfInterest(4));
fprintf('Inner Diameter of Operation Space (Wi): %g mm\n',quantityOfInterest(5));
fprintf('        Deformation to Compaction (uc): %g mm\n',quantityOfInterest(6));
fprintf('     Deformation to Plastic Yield (uY): %g mm\n',quantityOfInterest(7));
fprintf('  Deformation to End of Linearity (uL): %g mm\n',quantityOfInterest(8));


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

