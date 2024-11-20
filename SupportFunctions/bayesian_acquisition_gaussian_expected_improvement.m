function expectedImprovement = bayesian_acquisition_gaussian_expected_improvement(designSample,gaussianRegressionModelObjective,objetiveOptimalCurrent,explorationFactor)
%BAYESIAN_ACQUISITION_GAUSSIAN_EXPECTED_IMPROVEMENT for optimization
%	BAYESIAN_ACQUISITION_GAUSSIAN_EXPECTED_IMPROVEMENT is an acquisition 
%	function for bayesian optimization using gaussian process regression for
%	the surrogate model, which uses the expected improvement for the 
%	acquisition.
%	NOTE: contrary to most literature, this assumes a bayesian optimization
%	being used to minimize the objective function instead of maximize, so an
%	improvement here is the predicted value being lesser than current optimum.
%	This value should be maximized inside the bayesian optimizer to find the 
%	best candidate for a new optimum.
%	Nice visualization of different acquisition functions and the impact of
%	the exploration factor: https://distill.pub/2020/bayesian-optimization/ .
%
%	EXPECTEDIMPROVEMENT = BAYESIAN_ACQUISITION_GAUSSIAN_EXPECTED_IMPROVEMENT(
%	DESIGNSAMPLE,GAUSSIANREGRESSIONMODELOBJECTIVE,OBJETIVEOPTIMALCURRENT,
%	EXPLORATIONFACTOR) uses the current regression model 
%	GAUSSIANREGRESSIONMODELOBJECTIVE, current best known value 
%	OBJETIVEOPTIMALCURRENT and the exploration factor EXPLORATIONFACTOR to 
%	compute the expected improvement EXPECTEDIMPROVEMENT for each sample point
%	DESIGNSAMPLE.
%
%	Input:
%		- DESIGNSAMPLE : (nSample,nDimension) double
%		- GAUSSIANREGRESSIONMODELOBJECTIVE : RegressionGP
%		- OBJETIVEOPTIMALCURRENT : double
%		- EXPLORATIONFACTOR : double
%	Output:
%		- EXPECTEDIMPROVEMENT : (nSample,1) double
%
%   See also optimization_bayesian, 
%	bayesian_acquisition_gaussian_probability_of_improvement.
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


	% find predictions for sample points to be evaluated
	[predictionExpectedValue,predictionStandardDeviation] = gaussianRegressionModelObjective.predict(designSample);

	% compute predicted improvement based solely on expected values
    %   note: objetiveOptimalCurrent+explorationFactor and predictionExpectedValue 
    %   are inverted compared to standard implementations, as we are attempting to  
    %   minimize the function instead of maximize, so our improvement is having 
    %   predictionExpectedValue be lesser instead of greater
	predictedImprovement = objetiveOptimalCurrent+explorationFactor - predictionExpectedValue;

	% normalize the predicted improvement to a gaussian curve
	normalizedPredictedImprovement = predictedImprovement./predictionStandardDeviation;
	normalizedPredictedImprovement(predictionStandardDeviation==0) = 0;

    % get probability of improvement and uncertainty of improvement
    %	--> model that with: cumulative distribution and probability density
    %	--> (1) the higher the predicted improvement, the higher the probability
    %		that is true (CDF)
    %	--> (2) the closer the predicted improvement is to 0, the larger the 
    %		uncertainty of whether this is an improvement or not (PDF)
	probabilityOfImprovement = normcdf(normalizedPredictedImprovement);
	uncertaintyOfImprovement = normpdf(normalizedPredictedImprovement);

	% use both information to determine an expected improvement
	%	--> first term: the (probabilistic) expected value of improvement
	%	--> second term: the uncertainty term 
	expectedImprovement = predictedImprovement.*probabilityOfImprovement + predictionStandardDeviation.*uncertaintyOfImprovement;
	expectedImprovement(predictionStandardDeviation==0) = 0;
end