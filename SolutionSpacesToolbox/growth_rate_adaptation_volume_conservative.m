function growthAdaptationFactor = growth_rate_adaptation_volume_conservative(purity,targetPurity,nDimension,fractionAcceptableIncreaseMeasure,varargin)
%GROWTH_RATE_ADAPTATION_VOLUME_CONSERVATIVE Variant of volume ratio adaptation
%	GROWTH_RATE_ADAPTATION_VOLUME_CONSERVATIVE gives the growth rate adaptation  
%	factor computed using an extrapolation based on the volume proportionality 
%	with purity. This is based on the assumption that growth rate and purity are 
%	approximately related in the form g^d ~ 1-a.  
%
%	GROWTHADAPTATIONFACTOR = GROWTH_RATE_ADAPTATION_VOLUME_CONSERVATIVE(PURITY,
%	TARGETPURITY,NDIMENSION,FRACTIONACCEPTABLEINCREASEMEASURE) receives the 
%	current purity PURITY, the target purity TARGETPURITY, the problem dimension
%	NDIMENSION, and the fraction of good volume added in the previous iteration 
%	FRACTIONACCEPTABLEINCREASEMEASURE, and returns the growth adaptation factor 
%	GROWTHADAPTATIONFACTOR.
%
%   See also sso_box_stochastic, sso_component_stochastic, 
%	growth_rate_adaptation_linear, growth_rate_adaptation_volume.
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

	growthAdaptationFactor = ((1-targetPurity)./(1-purity)).^(1./growthExponent);
end