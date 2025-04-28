function performanceMeasure = one_mass_oscillator_continuous(designSample,systemParameter)
%ONE_MASS_OSCILLATOR_CONTINUOUS Bottom-Up Mapping (Mass-Spring-Damper System)
%   ONE_MASS_OSCILLATOR_CONTINUOUS simulates a single mass oscillator with
%   nonlinear spring and damping characteristics defined by continuous functions.
%   The system consists of a mass connected to a spring and damper, where both
%   the spring stiffness and damping coefficient can vary with displacement and
%   velocity respectively.
%
%   PERFORMANCEMEASURE = ONE_MASS_OSCILLATOR_CONTINUOUS(DESIGNSAMPLE,SYSTEMPARAMETER)
%   receives the spring stiffness and damping coefficient functions in DESIGNSAMPLE
%   and the system configuration in SYSTEMPARAMETER, returning the maximum
%   displacement for each design in PERFORMANCEMEASURE.
%
%   Example
%       For a linear spring (k=1000 N/m) and damper (c=10 Ns/m):
%           mass = 1; % [kg]
%           initialPosition = 0.1; % [m]
%           initialVelocity = 0; % [m/s]
%           timeSpan = 10; % [s]
%           springFunction = @(x) 1000*x;
%           damperFunction = @(v) 10*v;
%           designSample = {springFunction,damperFunction};
%           systemParameter = [mass,initialPosition,initialVelocity,timeSpan];
%           performanceMeasure = one_mass_oscillator_continuous(designSample,systemParameter);
%
%   Inputs:
%       - DESIGNSAMPLE : (nSample,2) cell
%           -- (1) : spring stiffness function handle k(x)
%           -- (2) : damping coefficient function handle c(v)
%       - SYSTEMPARAMETER : (1,4) double
%           -- (1) : mass [kg]
%           -- (2) : initial position [m]
%           -- (3) : initial velocity [m/s]
%           -- (4) : simulation time span [s]
%
%   Outputs:
%       - PERFORMANCEMEASURE : (nSample,1) double
%           -- maximum absolute displacement [m]
%
%   See also one_mass_oscillator_continuous_stiffness_only, ode45.
%
%   Copyright 2025 Eduardo Rodrigues Della Noce
%   SPDX-License-Identifier: Apache-2.0

    mass = systemParameter(1);
    initialPosition = systemParameter(2);
    initialVelocity = systemParameter(3);
    timeSpan = systemParameter(4);
    
    nSample = size(designSample,1);
    performanceMeasure = nan(nSample,1);
    for i=1:nSample
        currentSpringStiffnessFunction = designSample{i,1};
        currentDampingCoefficientFunction = designSample{i,2};

        currentSystem = @(time,state) one_mass_oscillator_continuous_derivative(time,state,currentSpringStiffnessFunction,currentDampingCoefficientFunction,mass);
        [~,solution] = ode45(currentSystem,[0 timeSpan],[initialPosition;initialVelocity]);
        
        [~,largestDisplacementIndex] = max(abs(solution(:,1)));
        performanceMeasure(i) = solution(largestDisplacementIndex,1);
    end
end

function derivativeState = one_mass_oscillator_continuous_derivative(time,state,springStiffnessFunction,dampingCoefficientFunction,mass)
%ONE_MASS_OSCILLATOR_CONTINUOUS_DERIVATIVE State derivative for oscillator system
%   Helper function that computes the state derivatives for the oscillator system
%   given the current state and system parameters.
%
%   Inputs:
%       - TIME : double
%           -- current time [s]
%       - STATE : (2,1) double
%           -- (1) : position [m]
%           -- (2) : velocity [m/s]
%       - SPRINGSTIFFNESSFUNCTION : function_handle
%           -- spring force function k(x)
%       - DAMPINGCOEFFICIENTFUNCTION : function_handle
%           -- damping force function c(v)
%       - MASS : double
%           -- oscillator mass [kg]
%
%   Outputs:
%       - DERIVATIVESTATE : (2,1) double
%           -- (1) : velocity [m/s]
%           -- (2) : acceleration [m/s^2]

    position = state(1);
    velocity = state(2);

    currentSprintgForce = springStiffnessFunction(position);
    currentDampingForce = dampingCoefficientFunction(velocity);

    forceApplied = -currentSprintgForce-currentDampingForce;

    derivativeState = [...
        velocity;
        forceApplied./mass];
end
