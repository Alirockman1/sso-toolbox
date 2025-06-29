function isClosest = design_closest_point_different_label(designSample,label,varargin)
%DESIGN_CLOSEST_POINT_DIFFERENT_LABEL 
%	DESIGN_CLOSEST_POINT_DIFFERENT_LABEL finds the design sample points 
%	that are the closest to other points that have a different label than them.
%	This can be achieved with different distance metrics.
%	
%	ISCLOSEST = DESIGN_CLOSEST_POINT_DIFFERENT_LABEL(DESIGNSAMPLE,LABEL) finds
%	the design sample points in DESIGNSAMPLE that are the closest to other
%	design points with a different logical label LABEL, returning that as
%	a logical array ISCLOSEST. This uses all default options of 'knnsearch'.
%
%	ISCLOSEST = DESIGN_CLOSEST_POINT_DIFFERENT_LABEL(...NAME,VALUE,...) allows
%	for using different options for 'knnsearch', as those name-value pair 
%	arguments are passed directly to said function.
%
%   Input:
%		- DESIGNSAMPLE : (nSample,nDesignVariable) double
%		- LABEL : (nSample,1) logical
%		- Name-value pair arguments: passed directly to 'knnsearch'
%
%   Output:
%		- ISCLOSEST : (nSample,1) logical
%
%   See also knnsearch, design_find_boundary_samples.
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

	positiveSample = designSample(label,:);
	negativeSample = designSample(~label,:);

	iPositiveBoundary = knnsearch(positiveSample,negativeSample,varargin{:});
	iNegativeBoundary = knnsearch(negativeSample,positiveSample,varargin{:});

	nSample = size(designSample,1);
	
	isClosestPositive = false(nSample,1);
	isClosestPositive(convert_index_base(label,iPositiveBoundary,'backward')) = true;

	isClosestNegative = false(nSample,1);
	isClosestNegative(convert_index_base(~label,iNegativeBoundary,'backward')) = true;

	isClosest = isClosestPositive | isClosestNegative;
end