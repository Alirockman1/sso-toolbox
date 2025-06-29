function [minDistanceOrder,minPathDistance] = design_sort_min_distance(queryDesign,varargin)
%DESIGN_SORT_MIN_DISTANCE Sort designs by minimum distance in sequence
%	DESIGN_SORT_MIN_DISTANCE returns an ordered index of the designs in a 
%	query that have the minimum distance to the previous selected one. 
%	The order starts at the first design of the query, then goes to the design
%	closest to the first one, then to the one closest to the second (not 
%	counting the first), and so on.
%	This is a greedy / nearest-neighbor strategy to solve a problem similar to 
%	a hamiltonian weighted path problem, except instead of the total distance
%	being the quantity to be minimized, it is the maximum distance traveled
%	at any point in the path, while moving to the closest point available.
%
%	MINDISTANCEORDER = DESIGN_SORT_MIN_DISTANCE(QUERYDESIGN) orders the designs 
%	present in the query QUERYDESIGN in such a way that each design is the 
%	closest available to the previous one, returning the order in 
%	MINDISTANCEORDER. To compute distances, 'knnsearch' is used. All given 
%	points in the query are used as initial points from which the path is 
%	computed, and the one with the smallest maximum distance between two points
%	is chosen as best.
%
%	MINDISTANCEORDER = DESIGN_SORT_MIN_DISTANCE(QUERYDESIGN,NINITIAL) allows one 
%	to set the number of points with which it starts finding the path NINITIAL. 
%	This goes in order of the designs given in QUERYDESIGN. In this way, one
%	can turn a quadratic-scaling problem into a linear one, at the cost of the
%	quality of the result.
%
%	MINDISTANCEORDER = DESIGN_SORT_MIN_DISTANCE(...,NAME,VALUE,...) allows one 
%	to set additional options for distance computation; these name-value pairs 
%	are used in 'knnsearch' directly.
%
%	[MINDISTANCEORDER,MINPATHDISTANCE] = DESIGN_SORT_MIN_DISTANCE(...)  
%	additionally returns the computed distance between each design and the 
%	previous one. The total path length can be computed by summing all values
%	in this array.
%
%   Inputs:
%       - QUERYDESIGN : (nQuery,nDesignVariable) double
%		- NINITIAL : double
%
%   Outputs:
%       - MINDISTANCEORDER : (nQuery,1) double
%		- MINPATHDISTANCE : (nQuery,1) double
%
%   See also knnsearch, design_select_max_distance.
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

	parser = inputParser;
	parser.KeepUnmatched = true;
	parser.StructExpand = true;
	parser.addOptional('NumberInitialPointsTested',[],@(x)isnumeric(x)&&isscalar(x)||isempty(x));
	parser.parse(varargin{:});
	nInitial = parser.Results.NumberInitialPointsTested;
	knnsearchOptions = namedargs2cell(parser.Unmatched);

	% initialize recurring arrays
	nSample = size(queryDesign,1);
	designOrder = nan(nSample,1);
    pathDistance = nan(nSample,1);
    availableQuery = true(nSample,1);

	if(isempty(nInitial) || nInitial>=nSample)
		nInitial = nSample;
	end
	if(nInitial<=0)
		nInitial = 1;
	end

    % initialize outputs
    shortestMaxDistance = inf;
   	minDistanceOrder = [];
   	minPathDistance = [];

   	% consider all designs as possible starting points and find the one that 
   	% leads to the best result
    for i=1:nInitial
    	% reset
    	designOrder(:) = nan;
    	pathDistance(:) = nan;
    	availableQuery(:) = true;

		% start from index i
		designOrder(1) = i; % first design
	    pathDistance(1) = 0; % start with 0 path distance
	    availableQuery(i) = false;

		% find distances to that design
	    [~,distance] = knnsearch(queryDesign(i,:),queryDesign(availableQuery,:),knnsearchOptions{:});

		% select the design with greatest distance and update the distance array
		% with the distance relative to the chosen design
		for j=2:nSample
			% get the queryDesign with the greatest distance to obtained samples
			[pathDistance(j),iAvailableSelect] = min(distance);
			designOrder(j) = convert_index_base(availableQuery,iAvailableSelect,'backward');
			availableQuery(designOrder(j)) = false;

			% update distances with new entry, assuming it will be evaluated
			[~,distance] = knnsearch(queryDesign(designOrder(j),:),queryDesign(availableQuery,:),knnsearchOptions{:});
		end

		% save if best result so far
		maxDistance = max(pathDistance);
		if(maxDistance<shortestMaxDistance)
			shortestMaxDistance = maxDistance;
		   	minDistanceOrder = designOrder;
		   	minPathDistance = pathDistance;
		end
	end
end