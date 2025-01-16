function performanceMeasure = car_crash_5d(designSample,systemParameter)
%CAR_CRASH_5D Bottom-Up Mapping (Simplified Car Crash - 5 Dimensional)
%   CAR_CRASH_5D tests designs for the car crash problem, where we use a 
%   simplified model to ascertain if a vehicle would succeed in not harming 
%   its passenger during a frontal car crash against a rigid wall.
%   The vehicle has two structures to absorb the energy; each one can
%   exert F1/F2 force during a displacement of maximally d1c/d2c,
%   respectively. The vehicle starts with velocity v0 and has mass m.
%
%   PERFORMANCEMEASURE = CAR_CRASH_5D(DESIGNSAMPLE,SYSTEMPARAMETER) receives 
%   the forces, vehicle mass and critical displacements in DESIGNSAMPLE, the 
%   initial speed in SYSTEMPARAMETER, and returns the remaining energy, maximum 
%   acceleration and order of deformation of the vehicle in PERFORMANCEMEASURE.
%
%   Example
%       Taking a simple case with two tests:
%           force1 = [350e3;300e3]; % [N]
%           force2 = [350e3;400e3]; % [N]
%           vehicleMass = [2000;1500]; % [kg]
%           criticalDisplacement1 = [0.3;0.2]; % [m]
%           criticalDisplacement2 = [0.3;0.4]; % [m] 
%           initialVehicleSpeed = 15.6; % [m/s]
%           designSample = [force1,force2,vehicleMass,criticalDisplacement1,criticalDisplacement2];
%           systemParameter = initialVehicleSpeed;
%           performanceMeasure = car_crash_5d(designSample,systemParameter);
%       PERFORMANCEMEASURE will contain the performance measures (energy 
%       remaining, maximum acceleration, order of deformation) for the entries 
%       (F1,F2,m,d1c,d2c) = (350e3,350e3,2000,0.3,0.3) and 
%       (300e3,400e3,1500,0.2,0.4). 
%
%   Inputs: 
%       - DESIGNSAMPLE : (nSample,5) double
%           -- (1) : force of component 1 (F1)
%           -- (2) : force of component 2 (F2)
%           -- (3) : mass of vehicle (m)
%           -- (4) : critical displacement of component 1 (d1c)
%           -- (5) : critical displacement of component 2 (d2c)
%       - SYSTEMPARAMETER : double
%           -- (1) : initial vehicle speed (v0)
%
%   Output:
%       - PERFORMANCEMEASURE : (nSample,3) double
%           -- (1) : energy remaining after crash (Erem)
%           -- (2) : maximum acceleration during crash (amax)
%           -- (3) : order of deformation of components (order)
%
%   See also car_crash_2d.
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

	% Unwrap inputs
	force1 = designSample(:,1); % F1
	force2 = designSample(:,2); % F2
    vehicleMass = designSample(:,3); % m
	criticalDisplacement1 = designSample(:,4); % d1c
	criticalDisplacement2 = designSample(:,5); % d2c
    initialVehicleSpeed = systemParameter(1); % v0
	
	% Energy remaining = energy vehcile has - deformation energy
	% Erem = 1/2*m*v0^2 - (F1*d1c + F2*d2c)
	energyRemaining = (1/2).*vehicleMass.*(initialVehicleSpeed.^2) - (force1.*criticalDisplacement1 + force2.*criticalDisplacement2);
	
	% Maximum Acceleration = force of strongest collapse load / total mass
	% amax = F2/m
	maximumAcceleration = force2./vehicleMass;

	% Deformation order = force of component 1 (lesser) - force of component 2 (greater)
	% order = F1 - F2
	deformationOrder = force1 - force2;
	
	% Wrap outputs
	performanceMeasure(:,1) = energyRemaining;
	performanceMeasure(:,2) = maximumAcceleration;
	performanceMeasure(:,3) = deformationOrder;
end