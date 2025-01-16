function finalPointError = brachistochrone_solve_analytical(radiusAngle,distanceXY)
%BRACHISTOCHRONE_SOLVE_ANALYTICAL Find brachistochrone parameters for problem
%   BRACHISTOCHRONE_SOLVE_ANALYTICAL allows for the computation of the  
%   brachistochrone parameters connecting points A (start) and B (end) when 
%   used together with a non-linear system solving function, like 'fsolve'. 
%   Point A is assumed to be in the origin (0,0).
%   Point B is defined by the problem parameters input in the form of horizontal
%   and vertical (x,y) distances to point A.
%
%   FINALPOINTERROR = BRACHISTOCHRONE_SOLVE_ANALYTICAL(RADIUSANGLE,DISTANCEXY)
%   receives the candidate cycloid radius and angle in RADIUSANGLE and the
%   expected distances between points B (end) and A (start) DISTANCEXY, 
%   returning the resulting quadratic error of the final position when using
%   those candidate cycloid parameters in FINALPOINTERROR.
%
%   Example
%       The brachistochrone curve connecting the two points can be found using 
%       this function like so:
%       distanceX = 10; % horizontal distance
%       distanceY = 5; % vertical distance
%       radiusAngleOptimal = fsolve(@(x)brachistochrone_solve_analytical(x,[distanceX,distanceY]),[1 1]);
%
%   Inputs:
%       - RADIUSANGLE : (1,2) double
%           -- (1) : radius of cycloid (r)
%           -- (2) : angle through rotation of cycloid (theta)
%       - DISTANCEXY : (1,2) double
%           -- (1) : total horizontal distance (X)
%           -- (2) : total vertical distance (Y)
%   Output:
%       - FINALPOINTERROR : (1,2) double
%           -- (1) : quadradic errors in final horizontal (X) position
%           -- (1) : quadradic errors in final vertical (Y) position
%
%   See also fsolve, bead_slide_time.
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

    cycloidRadius = radiusAngle(1); % r
    cycloidAngle = radiusAngle(2); % theta
    
    finalPointError(1) = (cycloidRadius*(cycloidAngle-sin(cycloidAngle)) - distanceXY(1))^2;
    finalPointError(2) = (cycloidRadius*(-1+cos(cycloidAngle)) + distanceXY(2))^2;
end