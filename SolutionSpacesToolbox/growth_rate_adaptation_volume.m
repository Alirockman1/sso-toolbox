function growthAdaptationFactor = growth_rate_adaptation_volume(purity,targetPurity,nDimension,fractionAcceptableIncreaseMeasure,varargin)
	parser = inputParser;
	parser.addParaemeter('MinimumGrowthExponent',1);
	parser.parse(varargin{:});
	options = parser.Results;

	growthExponent = nDimension - (nDimension-options.MinimumGrowthExponent)*fractionAcceptableIncreaseMeasure;

	growthAdaptationFactor = (((1-targetPurity)*purity)./((1-purity)*targetPurity)).^(1./growthExponent);
end