function performanceMeasure = one_mass_oscillator_continuous_stiffness_only(designSample,systemParameter)
%ONE_MASS_OSCILLATOR_CONTINUOUS_STIFFNESS_ONLY Bottom-Up Mapping (Mass-Spring System)
%   ONE_MASS_OSCILLATOR_CONTINUOUS_STIFFNESS_ONLY simulates a single mass
%   oscillator with nonlinear spring characteristics and linear damping. This is
%   a simplified version of ONE_MASS_OSCILLATOR_CONTINUOUS where only the spring
%   stiffness is allowed to vary, while the damping coefficient remains constant.
%
%   PERFORMANCEMEASURE = ONE_MASS_OSCILLATOR_CONTINUOUS_STIFFNESS_ONLY(DESIGNSAMPLE,
%   SYSTEMPARAMETER) receives the spring stiffness function in DESIGNSAMPLE and
%   the system configuration (including damping) in SYSTEMPARAMETER, returning
%   the maximum displacement for each design in PERFORMANCEMEASURE.
%
%   Example
%       For a nonlinear spring with k(x)=1000*x^3:
%           mass = 1; % [kg]
%           initialPosition = 0.1; % [m]
%           initialVelocity = 0; % [m/s]
%           timeSpan = 10; % [s]
%           dampingCoefficient = 10; % [Ns/m]
%           springFunction = @(x) 1000*x.^3;
%           systemParameter = [mass,initialPosition,initialVelocity,timeSpan,dampingCoefficient];
%           performanceMeasure = one_mass_oscillator_continuous_stiffness_only(springFunction,systemParameter);
%
%   Inputs:
%       - DESIGNSAMPLE : (nSample,1) cell
%           -- spring stiffness function handle k(x)
%       - SYSTEMPARAMETER : (1,5) double
%           -- (1) : mass [kg]
%           -- (2) : initial position [m]
%           -- (3) : initial velocity [m/s]
%           -- (4) : simulation time span [s]
%           -- (5) : damping coefficient [Ns/m]
%
%   Outputs:
%       - PERFORMANCEMEASURE : (nSample,1) double
%           -- maximum absolute displacement [m]
%
%   See also one_mass_oscillator_continuous.
%
%   Copyright 2025 Eduardo Rodrigues Della Noce
%   SPDX-License-Identifier: Apache-2.0

    dampingCoefficient = systemParameter(5);
    dampingFunction = @(x) dampingCoefficient.*x;

    designSampleNew = cell(size(designSample,1),2);
    designSampleNew(:,1) = designSample;
    designSampleNew(:,2) = {dampingFunction};

    performanceMeasure = one_mass_oscillator_continuous(designSampleNew,systemParameter(1:4));
end
