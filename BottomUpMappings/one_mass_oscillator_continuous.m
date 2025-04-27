function performanceMeasure = one_mass_oscillator_continuous(designSample,systemParameter)
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
    position = state(1);
    velocity = state(2);

    currentSprintgForce = springStiffnessFunction(position);
    currentDampingForce = dampingCoefficientFunction(velocity);

    forceApplied = -currentSprintgForce-currentDampingForce;

    derivativeState = [...
        velocity;
        forceApplied./mass];
end
