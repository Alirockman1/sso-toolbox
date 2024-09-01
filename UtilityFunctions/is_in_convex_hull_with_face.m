function [isInside,score] = is_in_convex_hull_with_face(facePoint,faceNormalIn,queryPoint)
%IS_IN_CONVEX_HULL_WITH_FACE Check if point is inside convex hull given faces
%   IS_IN_CONVEX_HULL_WITH_FACE uses the information on the faces (edges,
%   planes, hyperplanes, ...) of a convex hull to check if a query point is 
%   inside or outside that hull. 
%
%   ISINSIDE = IS_IN_CONVEX_HULL_WITH_FACE(FACEPOINT,FACENORMALIN,QUERYPOINT)
%   receives points that are in each convex hull face FACEPOINT, vectors which 
%   are normal to said faces and point to the inside of the hull FACENORMALIN, 
%   and query points QUERYPOINT, and returns whether or not each point is 
%   inside the hull in ISINSIDE (true = inside, false = outside).
%
%   [ISINSIDE,SCORE] = IS_IN_CONVEX_HULL_WITH_FACE(...) also returns a score
%   for each point SCORE. This score is negative for points inside the hull and
%   positive for points outside, and the farther the result is from 0, the 
%   farther away the point is from any face.
%
%   Input:
%       - FACEPOINT : (nFace,nDimension) double
%       - FACENORMALIN : (nFace,nDimension) double
%       - QUERYPOINT : (nSample,nDimension) double
%
%   Output:
%       - ISINSIDE : (nSample,1) logical
%       - SCORE : (nSample,1) double
%
%   See also convex_hull_face.
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

    % see if dot product of each design with the normal vectors is negative
    % if so --> inside
    % if not --> outside
    nQuery = size(queryPoint,1);
    isInside = false(nQuery,1);
    score = nan(nQuery,1);
    for i=1:nQuery
        dotProduct = dot(facePoint - queryPoint(i,:),faceNormalIn,2);

        isInside(i) = all(dotProduct<=0);
        score(i) = max(dotProduct);
    end
end