function [isInside,score] = is_in_convex_hull_linear_programming(convexHullPoint,queryPoint,varargin)
%IS_IN_CONVEX_HULL_LINEAR_PROGRAMMING without computing the convex hull itself
%	IS_IN_CONVEX_HULL_LINEAR_PROGRAMMING allows one to verify if a set of query 
%	points is inside the ocnvex hull formed by a different set of points without
%	having to compute the convex hull itself. 
%	This is achieved by using a property of any convex hull that all points 
%	inside it are linear combinations of the points inside it with linear 
%	coefficients that sum up to 1. 
%	Therefore, a linear programming optimization problem can be solved to see
%	if such a representation is possible for the query points.
%	Reference:
%	- https://stackoverflow.com/questions/4901959/find-if-a-point-is-inside-a-convex-hull-for-a-set-of-points-without-computing-th#11731437
%	- "Convex Optimization", Stephen Boyd and Lieven Vandenberghe, Cambridge 
%	University Press - https://web.stanford.edu/~boyd/cvxbook/
%
%	ISINSIDE = IS_IN_CONVEX_HULL_LINEAR_PROGRAMMING(CONVEXHULLPOINT,QUERYPOINT)
%	checks for each point in QUERYPOINT if it is inside the convex hull 
%	encompassing the points CONVEXHULLPOINT, returning the label in ISINSIDE.
%
%	ISINSIDE = IS_IN_CONVEX_HULL_LINEAR_PROGRAMMING(...,NAME,VALUE,...) allows
%	for the specification of additional options for 'linprog'. These name-value
%	pair arguments get passed to 'optimoptions' set for 'linprog'.
%
%	[ISINSIDE,SCORE] = IS_IN_CONVEX_HULL_LINEAR_PROGRAMMING(...) also returns
%	a score, which is -1 for successful optimizations (is inside) and 1 for
%	failures (is outside).
%
%	Input:
%		- CONVEXHULLPOINT : (nConvexHullPoint,nDesignVariable) double
%		- QUERYPOINT : (nQueryPoint,nDesignVariable) double
%		- Name-value pair arguments: passed directly to optimoptions.
%
%	Output:
%		- ISINSIDE : (nQueryPoint,1) logical
%		- SCORE : (nQueryPoint,1) double
%
%   See also linprog, optimoptions, is_in_convex_hull_with_face.
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

	nConvexHullPoint = size(convexHullPoint,1);
	equalityConstraintMatrix = [convexHullPoint';ones(1,nConvexHullPoint)];
    objectiveCoefficientVector = ones(1,nConvexHullPoint);
    lowerBound = zeros(1,nConvexHullPoint);
    upperBound = ones(1,nConvexHullPoint);

    % linprog options
    defaultLinprogOptions = {'Display','off'};
    [~,mergedLinprogOptions] = merge_name_value_pair_argument(defaultLinprogOptions,varargin);
    linprogOptions = optimoptions('linprog',mergedLinprogOptions{:});

	nSample = size(queryPoint,1);
	isInside = false(nSample,1);
    score = nan(nSample,1);
	for i=1:nSample
		equalityConstraintRightHandeSide = [queryPoint(i,:)';1];
		[~,~,exitFlag] = linprog(objectiveCoefficientVector,...
			[],... % A
			[],... % b
			equalityConstraintMatrix,... % Aeq
			equalityConstraintRightHandeSide,... % beq
			lowerBound,... % lb
			upperBound,... % ub
			linprogOptions);

        isInside(i) = (exitFlag==1); % converged
	    score(i) = isInside(i)*1 + (~isInside(i))*(-1);
    end
end