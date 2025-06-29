function performanceMeasure = car_crash_2d(designSample,systemParameter)
%CAR_CRASH_2D Bottom-Up Mapping (Simplified Car Crash - 2 Dimensional)
%   CAR_CRASH_2D is a further simplified version of the car crash problem, where
%	only the compression forces are defined as design variables.
%   The vehicle has two structures to absorb the energy; each one can exert 
%	F1/F2 force through a critical displacement of maximally d1c/d2c, 
%	respectively. The vehicle starts with velocity v0 and has mass m.
%
%   PERFORMANCEMEASURE = CAR_CRASH_2D(DESIGNSAMPLE,SYSTEMPARAMETER) 
%   receives the forces in DESIGNSAMPLE and the vehicle mass, critical 
%   displacements and initial speed in SYSTEMPARAMETER and returns the remaining 
%   energy, maximum acceleration and order of deformation of the vehicle in 
% 	PERFORMANCEMEASURE.
%
%   Example
%       Taking a simple case with two tests:
%           vehicleMass = 2000; % [kg]
%           displacementCritical1 = 0.3; % [m]
%           displacementCritical2 = 0.3; % [m]
%           initialVelocity = 15.6; % [m/s]
%           force1 = [350e3;300e3]; % [N]
%           force2 = [350e3;400e3]; % [N]
%           designSample = [force1,force2];
%           systemParameter = [vehicleMass,displacementCritical1,displacementCritical2,initialVelocity];
%           performanceMeasure = car_crash_2d(designSample,systemParameter);
%       PERFORMANCEMEASURE will contain the quantities of interest (energy 
%       remaining, maximum acceleration, order of deformation) for the pairs of 
%       forces (F1,F2) = (350e3,350e3) and (300e3,400e3). 
%
%   Inputs: 
%       - DESIGNSAMPLE : (nSample,2) double
%           -- (1) : force of component 1 (F1)
%           -- (2) : force of component 2 (F2)
%       - SYSTEMPARAMETER : (1,4) double  
%           -- (1) : mass of vehicle (m)
%           -- (2) : critical displacement of component 1 (d1c)
%           -- (3) : critical displacement of component 2 (d2c)
%           -- (4) : initial vehicle speed (v0)
%
%   Output:
%       - PERFORMANCEMEASURE : (nSample,3) double
%           -- (1) energy remaining after crash (Erem)
%           -- (2) maximum acceleration during crash (amax)
%           -- (3) order of deformation of components (order)
%
%   See also car_crash_5d.
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

    % create wrapper to call 5d function
    nSample = size(designSample,1);
	designSample5d = nan(nSample,5);
	designSample5d(:,[1,2]) = designSample;
	designSample5d(:,[3,4,5]) = repmat(systemParameter(1:3),nSample,1);

	% call more complex function
	performanceMeasure = car_crash_5d(designSample5d,systemParameter(4));
end