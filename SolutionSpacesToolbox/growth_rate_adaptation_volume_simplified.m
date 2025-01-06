function growthAdaptationFactor = growth_rate_adaptation_volume_simplified(purity,targetPurity,nDimension,fractionAcceptableIncreaseMeasure,varargin)
	parser = inputParser;
	parser.addParameter('MinimumGrowthExponent',1);
	parser.parse(varargin{:});
	options = parser.Results;

	growthExponent = nDimension - (nDimension-options.MinimumGrowthExponent)*fractionAcceptableIncreaseMeasure;

	growthAdaptationFactor = ((1-targetPurity)./(1-purity)).^(1./growthExponent);
end