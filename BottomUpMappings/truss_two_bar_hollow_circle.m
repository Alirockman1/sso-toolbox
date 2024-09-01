function performanceMeasure = truss_two_bar_hollow_circle(designSample,systemParameter)
%TWO_BAR_TRUSS_HOLLOW_CIRCLE Bottom-Up Mapping (Bar Selection Problem)
%   TWO_BAR_TRUSS_HOLLOW_CIRCLE is a bottom-up mapping for a two-bar
%   truss system (cross-sections for the bars: hollow circle) with a force   
%   applied at the tip.
%     |\
%     | \
%     |  \ 
%     |   \
%   h1|    \
%     |     \ bar1 (l1,A1,J1,m1)
%     | all  \
%    -|------- |
%     |alpha2/ | 
%     |     /  V F
%     |    / 
%   h2|   / bar2 (l2,A2,J2,m2)
%     |  /
%     | /
%     |/
%
%   PERFORMANCEMEASURE = TWO_BAR_TRUSS_HOLLOW_CIRCLE(DESIGNSAMPLE,
%   SYSTEMPARAMETER) receives the thickness, outer diameter and length of each 
%   of two truss bars in DESIGNSAMPLE, their properties and force applied in 
%   SYSTEMPARAMETER, and calculates the overall system performance 
%   PERFORMANCEMEASURE when said load is applied at the tip of the structure, 
%   including yielding and buckling factors.
%
%   Example
%       To test an example with a steel bar on top and a aluminum bar 
%       on the bottom:
%           designSample = [1e-2 1e-1 1 2e-2 1e-1 1.5]; % t1[m],D1[m],l1[m],t2[m],D2[m],l2[m] 
%           systemParameter = [210e9 7800  355e6... % (E,rho,sigmaY) bar 1 - steel
%                               70e9 2700  100e6... % (E,rho,sigmaY) bar 2 - aluminum
%                                  1 1000... % 1m distance; 1kN force
%                             ];
%           performanceMeasure = TWO_BAR_TRUSS_HOLLOW_CIRCLE(DV,PARAM);
%
%   Inputs:
%       - DESIGNSAMPLE : (nSample,6) double 
%           -- (1) thickness of bar 1 (top) (t1)
%           -- (2) outer diameter of bar 1 (top) (D1)
%           -- (3) length of bar 1 (top) (l1)
%           -- (4) thickness of bar 2 (bottom) (t2)
%           -- (5) outer diameter of bar 2 (bottom) (D2)
%           -- (6) length of bar 2 (bottom) (l2)
%       - SYSTEMPARAMETER : (1,8) double
%           -- (1) Young's modulus of bar 1 (top) (E1)
%           -- (2) density of bar 1 (top) (rho1)
%           -- (3) yield strength of bar 1 (top) (sigmaYield2)
%           -- (4) Young's modulus of bar 2 (bottom) (E2)
%           -- (5) density of bar 2 (bottom) (rho2)
%           -- (6) yield strength of bar 2 (bottom) (sigmaYield2)
%           -- (7) force applied at the tip of the structure (F)
%           -- (8) distance from wall where force is applied (tip) (l)
%
%   Output:
%       - PERFORMANCEMEASURE : (nSample,29) double
%           -- (1) reserve yield force of bar 1 (rY1)
%           -- (2) reserve yield force of bar 2 (rY2)
%           -- (3) reserve buckling force of bar 1 (rB1)
%           -- (4) reserve buckling force of bar 2 (rB2)
%           -- (5) total mass of structure (m)
%           -- (6) total tip displacement (u)
%           -- (7) (undeformed) angle of bar 1 (alpha1)
%           -- (8) (undeformed) angle of bar 2 (alpha2)
%           -- (9) height of wall-support of bar 1 (h1)
%           -- (10) height of wall-support of bar 2 (h2)
%           -- (11) inner diameter of bar 1 (d1)
%           -- (12) inner diameter of bar 2 (d2)
%           -- (13) cross-section area of bar 1 (A1)
%           -- (14) cross-section area of bar 2 (A2)
%           -- (15) moment of inertia of bar 1 (I1)
%           -- (16) moment of inertia of bar 2 (I2)
%           -- (17) angle-normalized force on bar 1 (S1)
%           -- (18) angle-normalized force on bar 2 (S2)
%           -- (19) tension force on bar 1 (F1)
%           -- (20) compression force on bar 2 (F2)
%           -- (21) yield force limit of bar 1 (FY1)
%           -- (22) yield force limit of bar 2 (FY2)
%           -- (23) buckling force limit of bar 1 (FB1)
%           -- (24) buckling force limit of bar 2 (FB2)
%           -- (25) mass of bar 1 (m1)
%           -- (26) mass of bar 2 (m2)
%           -- (27) normal-stiffness of bar 1 (k1)
%           -- (28) normal-stiffness of bar 2 (k2)
%           -- (29) compliance of total structure (NF)
%
%   See also truss_six_bar_2d_moving_node.
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

	%% Unwrap inputs
    % design variables
    wallThickness1 = designSample(:,1);
    outerDiameter1 = designSample(:,2);
    length1 = designSample(:,3);
    wallThickness2 = designSample(:,4);
    outerDiameter2 = designSample(:,5);
    length2 = designSample(:,6);
    
    % systemParametereters
    youngsModulus1 = systemParameter(1);
    desity1 = systemParameter(2);
    yieldStress1 = systemParameter(3);
    youngsModulus2 = systemParameter(4);
    density2 = systemParameter(5);
    yieldStress2 = systemParameter(6);
    appliedForce = systemParameter(7);
    tipDistance = systemParameter(8);
    
    
    %% Geometric system properties
    % angles
    alpha1 = acos(tipDistance./length1);
    alpha2 = acos(tipDistance./length2);
    
    % vertical distance to center
    height1 = sin(alpha1).*length1;
    height2 = sin(alpha2).*length2;
    
    % inner diameter
    innerDiameter1 = outerDiameter1-2*wallThickness1;
    innerDiameter2 = outerDiameter2-2*wallThickness2;
    
    % areas
    area1 = pi.*outerDiameter1.^2/4 - pi.*innerDiameter1.^2/4;
    area2 = pi.*outerDiameter2.^2/4 - pi.*innerDiameter2.^2/4;

    % moment of inertia
    momentOfInertia1 = pi/64.*(outerDiameter1.^4 - innerDiameter1.^4);
    momentOfInertia2 = pi/64.*(outerDiameter2.^4 - innerDiameter2.^4);
    
    
    %% Force applied to each bar 
    normalizedUnitForce1 = cos(alpha2)./sin(alpha1+alpha2);
    normalizedUnitForce2 = cos(alpha1)./sin(alpha1+alpha2);
    force1 = appliedForce .* normalizedUnitForce1;
    force2 = - appliedForce .* normalizedUnitForce2;
    
    
    %% Yielding reserve
    yieldForce1 = area1.*yieldStress1;
    yieldForce2 = area2.*yieldStress2;
    yieldReserve1 = (yieldForce1 - abs(force1))./yieldForce1;
    yieldReserve2 = (yieldForce2 - abs(force2))./yieldForce2;
    
    
    %% Buckling reserve
    bucklingForce1 = pi^2 * youngsModulus1 * momentOfInertia1./(length1.^2);
    bucklingForce2 = pi^2 * youngsModulus2 * momentOfInertia2./(length2.^2);
    reserveBuckling1 = (bucklingForce1 + force1)./bucklingForce1;
    reserveBuckling2 = (bucklingForce2 + force2)./bucklingForce2;
    
    
    %% Total mass
    mass1 = desity1.*area1.*length1;
    mass2 = density2.*area2.*length2;
    totalMass = mass1 + mass2;
    
    
    %% Total displacement
    stiffness1 = youngsModulus1.*area1./length1;
    stiffness2 = youngsModulus2.*area2./length2;
    totalCompliance = (normalizedUnitForce1.^2)./stiffness1 + (normalizedUnitForce2.^2)./stiffness2;
    tipDisplacement = appliedForce.*totalCompliance;
    
    
    %% Wrap Outputs
    performanceMeasure(:,1) = yieldReserve1;
    performanceMeasure(:,2) = yieldReserve2;
    performanceMeasure(:,3) = reserveBuckling1;
    performanceMeasure(:,4) = reserveBuckling2;
    performanceMeasure(:,5) = totalMass;
    performanceMeasure(:,6) = tipDisplacement;
    performanceMeasure(:,7) = alpha1;
    performanceMeasure(:,8) = alpha2;
    performanceMeasure(:,9) = height1;
    performanceMeasure(:,10) = height2;
    performanceMeasure(:,11) = innerDiameter1;
    performanceMeasure(:,12) = innerDiameter2;
    performanceMeasure(:,13) = area1;
    performanceMeasure(:,14) = area2;
    performanceMeasure(:,15) = momentOfInertia1;
    performanceMeasure(:,16) = momentOfInertia2;
    performanceMeasure(:,17) = normalizedUnitForce1;
    performanceMeasure(:,18) = normalizedUnitForce2;
    performanceMeasure(:,19) = force1;
    performanceMeasure(:,20) = force2;
    performanceMeasure(:,21) = yieldForce1;
    performanceMeasure(:,22) = yieldForce2;
    performanceMeasure(:,23) = bucklingForce1;
    performanceMeasure(:,24) = bucklingForce2;
    performanceMeasure(:,25) = mass1;
    performanceMeasure(:,26) = mass2;
    performanceMeasure(:,27) = stiffness1;
    performanceMeasure(:,28) = stiffness2;
    performanceMeasure(:,29) = totalCompliance;
end