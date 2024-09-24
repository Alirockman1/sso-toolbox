function [convexHullIndex,measure] = compute_convex_hull(convexHullPoint,varargin)
%COMPUTE_CONVEX_HULL Compute the convex hull for any number of dimensions
%   COMPUTE_CONVEX_HULL computes the convex hull for the points given.
%
%   CONVEXHULLINDEX = COMPUTE_CONVEX_HULL(CONVEXHULLPOINT) receives the sample 
%   points CONVEXHULLPOINT and returns the indices of the vertices of the convex
%   hull that includes all points CONVEXHULLINDEX. 
%
%   [CONVEXHULLINDEX,MEASURE] = COMPUTE_CONVEX_HULL(CONVEXHULLPOINT) also  
%   returns the measure (1D: length, 2D: area, 3D: volume, ...) of the convex 
%   hull.
%
%   [...] = COMPUTE_CONVEX_HULL(...,NAME,VALUE,...) allows for the specification
%   of additional options. These get passed directly to 'convhull' or 
%   'convhulln' depending on the dimension of the problem. For 2D and 3D 
%   problems, the default value is {'Simplify',true}, while for other 
%   dimensions, it is empty.
%
%   Input:
%       - CONVEXHULLPOINT : (nSample,nDimension) double
%
%   Output:
%       - CONVEXHULLINDEX : (nFacet,nDimension) integer
%       - MEASURE : double
%
%   See also convhull, convhulln, find_facet_reference_point_normal, 
%   is_in_convex_hull_with_facet_normal.
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

	nDimension = size(convexHullPoint,2);
    if(nDimension==2 || nDimension==3)
        defaultConvhullOptions = {'Simplify',true};
    else
        defaultConvhullOptions = {};
    end
    [~,convhullOptions] = merge_name_value_pair_argument(defaultConvhullOptions,varargin{:});

    if(nDimension==1)
        % take maximum and minimum, then use as normals simply +1 for lower limit
        % and -1 for upper limit
        [pointLowerBoundary,iLowerBoundary] = min(convexHullPoint);
        [pointUpperBoundary,iUpperBoundary] = max(convexHullPoint);
        convexHullIndex = [iLowerBoundary;iUpperBoundary];
        measure = pointUpperBoundary - pointLowerBoundary;
    elseif(nDimension==2)
    	[convexHullIndex,measure] = convhull(convexHullPoint,convhullOptions{:});
		nEdge = size(convexHullIndex,1)-1;
        convexHullIndexNew = nan(nEdge,2);
        for i=1:nEdge
            convexHullIndexNew(i,:) = [convexHullIndex(i),convexHullIndex(i+1)];
        end
        convexHullIndex = convexHullIndexNew;
	elseif(nDimension==3)
    	[convexHullIndex,measure] = convhull(convexHullPoint,convhullOptions{:});
    else
    	[convexHullIndex,measure] = convhulln(convexHullPoint,convhullOptions{:});
    end
end