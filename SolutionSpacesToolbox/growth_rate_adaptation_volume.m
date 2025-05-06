function growthAdaptationFactor = growth_rate_adaptation_volume(purity,targetPurity,nDimension,fractionAcceptableIncreaseMeasure,varargin)
%GROWTH_RATE_ADAPTATION_VOLUME Adapt growth rate with purity using volume ratio
%	GROWTH_RATE_ADAPTATION_VOLUME gives the growth rate adaptation factor 
%	computed using an extrapolation based on the volume proportionality with
%	purity. This is based on the assumption that growth rate and purity are 
%	related in the form g^d ~ (1-a)/a.  
%
%	GROWTHADAPTATIONFACTOR = GROWTH_RATE_ADAPTATION_VOLUME(PURITY,TARGETPURITY,
%	NDIMENSION,FRACTIONACCEPTABLEINCREASEMEASURE) receives the current purity 
%	PURITY, the target purity TARGETPURITY, the problem dimension NDIMENSION, 
%	and the fraction of good volume added in the previous iteration 
%	FRACTIONACCEPTABLEINCREASEMEASURE, and returns the growth adaptation factor 
%	GROWTHADAPTATIONFACTOR.
%
%   See also sso_box_stochastic, sso_component_stochastic, 
%	growth_rate_adaptation_linear, growth_rate_adaptation_volume_conservative.
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

	parser = inputParser;
	parser.addParameter('MinimumGrowthExponent',1);
	parser.parse(varargin{:});
	options = parser.Results;

	growthExponent = nDimension - (nDimension-options.MinimumGrowthExponent)*fractionAcceptableIncreaseMeasure;

	growthAdaptationFactor = (((1-targetPurity)*purity)./((1-purity)*targetPurity)).^(1./growthExponent);
end