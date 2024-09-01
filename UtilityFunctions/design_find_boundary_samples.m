function isBoundary = design_find_boundary_samples(designSample,label,distanceMetric)
%DESIGN_FIND_BOUNDARY_SAMPLES Estimate boundary points by using varying metrics
%	DESIGN_FIND_BOUNDARY_SAMPLES uses different distance measures to compute 
%	the design sample points that are closest to other points with a different 
%	label, and collages all results into an estimation of the designs that would
%	define a boundary between the regions of different labels.
%
%	ISBOUNDARY = DESIGN_FIND_BOUNDARY_SAMPLES(DESIGNSAMPLE,LABEL) receives the
%	design sample points DESIGNSAMPLE and their respective labels LABEL, 
%	returning a flag as to whether each design is in the boundary or not 
%	ISBOUNDARY. 
%
%	ISBOUNDARY = DESIGN_FIND_BOUNDARY_SAMPLES(DESIGNSAMPLE,LABEL,DISTANCEMETRIC)
%	allows one to additionally choose which distance metrics are to be used
%	to find the closest designs. By default, the distance metrics used are 
%	'euclidean', 'seuclidean', 'cityblock', 'chebychev' and 'mahalanobis'.
%
%   Input:
%		- DESIGNSAMPLE : (nSample,nDesignVariable) double
%		- LABEL : (nSample,1) logical
%		- DISTANCEMETRIC : (1,nDistance) cell
%
%   Output:
%		- ISBOUNDARY : (nSample,1) logical
%
%   See also design_closest_point_different_label.
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
	
	if(nargin<3 || isempty(distanceMetric))
		distanceMetric = {...
		    'euclidean',...
		    'seuclidean',...
		    ...'fasteuclidean',... % already using euclidean
		    ...'fastseuclidean',... % already using seuclidean
		    'cityblock',...
		    'chebychev',...
		    ...'minkowski',... % equivalent to euclidean if no other parameter is given
		    'mahalanobis'...
		    ...'cosine',... % does not have wanted behavior
		    ...'correlation',... % includes 'interior' designs
		    ...'spearman',... % includes 'interior' designs
		    ...'hamming',... % includes 'interior' designs
		    ...'jaccard',... % includes 'interior' designs
		};
	end

    nSample = size(designSample,1);
    nDistance = length(distanceMetric);

    isClosest = false(nSample,nDistance);
    for i=1:nDistance
        isClosest(:,i) = design_closest_point_different_label(...
            designSample,...
            label,...
            'K',1,...
            'IncludeTies',false,...
            'Distance',distanceMetric{i});
    end
    isBoundary = any(isClosest,2);
end

