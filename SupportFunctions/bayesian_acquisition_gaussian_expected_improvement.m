function expectedImprovement = bayesian_acquisition_gaussian_expected_improvement(gaussianRegressionModelObjective,designSample,objetiveOptimalCurrent,explorationFactor)
	% find predictions for sample points to be evaluated
	[predictionExpectedValue,predictionStandardDeviation] = gaussianRegressionModelObjective.predict(designSample);

	% compute predicted improvement based solely on expected values
	predictedImprovement = predictionExpectedValue - objetiveOptimalCurrent - explorationFactor;

	% compute variation information (based on how unknown each sample is)
	normalTransformation = predictedImprovement./predictionStandardDeviation;
	normalTransformation(predictionStandardDeviation==0) = 0;
	sampleCumulativeDistribution = normcdf(normalTransformation);
	sampleProbabilityDensity = normpdf(normalTransformation);

	% use both predicted expected values and variation information to estimate the actual expected improvement
	expectedImprovement = predictedImprovement.*sampleCumulativeDistribution + predictionStandardDeviation.*sampleProbabilityDensity;
	expectedImprovement(predictionStandardDeviation==0) = 0;
end