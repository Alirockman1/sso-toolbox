function growthAdaptationFactor = growth_rate_adaptation_linear(purity,targetPurity,varargin)
%GROWTH_RATE_ADAPTATION_LINEAR Adapt growth rate linearly with purity
%	GROWTH_RATE_ADAPTATION_LINEAR gives the growth rate adaptation factor 
%	computed from a ratio between the current purity and target purity.
%
%	GROWTHADAPTATIONFACTOR = GROWTH_RATE_ADAPTATION_LINEAR(PURITY,TARGETPURITY)
%	receives the current purity PURITY and the target purity TARGETPURITY
%	and returns the growth adaptation factor GROWTHADAPTATIONFACTOR, which is 
%	the ratio of purity and the target purity.
%
%   See also sso_box_stochastic, sso_component_stochastic, 
%	growth_rate_adaptation_volume, growth_rate_adaptation_volume_conservative.
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
	growthAdaptationFactor = purity/targetPurity;
end