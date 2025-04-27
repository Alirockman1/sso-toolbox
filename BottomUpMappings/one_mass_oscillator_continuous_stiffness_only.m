function performanceMeasure = one_mass_oscillator_continuous_stiffness_only(designSample,systemParameter)
    dampingCoefficient = systemParameter(5);
    dampingFunction = @(x) dampingCoefficient.*x;

    designSampleNew = cell(size(designSample,1),2);
    designSampleNew(:,1) = designSample;
    designSampleNew(:,2) = {dampingFunction};

    performanceMeasure = one_mass_oscillator_continuous(designSampleNew,systemParameter(1:4));
end
