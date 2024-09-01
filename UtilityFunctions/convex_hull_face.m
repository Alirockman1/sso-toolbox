function [convexHullIndex,measure,facePoint,faceNormalIn] = convex_hull_face(convexHullPoint,varargin)
%CONVEX_HULL_FACE Compute the convex hull and return information on the faces
%   CONVEX_HULL_FACE computes the convex hull for the points given if necessary
%   and gets information regarding the faces (edges, planes, hyperplanes, ...)
%   of said convex hull, computing a point that lies within each plane, and 
%   a normal vector that points inwards to the convex hull.
%
%   CONVEXHULLINDEX = CONVEX_HULL_FACE(CONVEXHULLPOINT) receives the sample 
%   points CONVEXHULLPOINT and returns the indices of the vertices of the convex
%   hull that includes all points CONVEXHULLINDEX. 
%
%   [CONVEXHULLINDEX,MEASURE] = CONVEX_HULL_FACE(CONVEXHULLPOINT) also returns 
%   the measure (1D: length, 2D: area, 3D: volume, ...) of the convex hull.
%
%   [CONVEXHULLINDEX,MEASURE] = CONVEX_HULL_FACE(CONVEXHULLPOINT,
%   CONVEXHULLINDEX) allows one to already specify the convex hull index of the
%   vertices. In this case, the convex hull is not recomputed, and MEASURE is
%   left empty. CONVEXHULLINDEX may be altered if it was previously computed
%   for a 2D problem using 'convhull' to suit the style of output from 
%   'convhulln'.
%
%   [CONVEXHULLINDEX,MEASURE,FACEPOINT,FACENORMALIN] = CONVEX_HULL_FACE(...) 
%   additionally returns both points that lie within each face of the convex
%   hull FACEPOINT and a normalized vector which is normal to the respective
%   face and points to the inside of the convex hull FACENORMALIN.
%
%   [...] = CONVEX_HULL_FACE(...,NAME,VALUE,...) allows for the specification
%   of additional options. These are:
%       - 'ConvhullOptions' : options for the computation of the convex hull,
%       either 'convhull' or ' convhulln' depending on problem dimension. For
%       2D and 3D problems, the default value is {'Simplify',true}, while for 
%       other dimensions, it is empty.
%       - 'FacePointPreference' : preference regarding which point to choose
%       as the reference FACEPOINT. There are two possibilities:
%           -- 'center' : center of the face.
%           -- 'vertice' : first vertice in the convex hull index.
%
%   Input:
%       - CONVEXHULLPOINT : (nSample,nDimension) double
%       - CONVEXHULLINDEX : (nFace,1) integer (2D) OR (nFace,nDimension) integer
%       - 'ConvhullOptions' : cell
%       - 'FacePointPreference' : char OR string
%
%   Output:
%       - CONVEXHULLINDEX : (nFace,nDimension) integer
%       - MEASURE : double
%       - FACEPOINT : (nFace,nDimension) double
%       - FACENORMALIN : (nFace,nDimension) double
%
%   See also convhull, convhulln, is_in_convex_hull_with_face.
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
	parser.addOptional('ConvexHullIndex',[]);
	parser.addParameter('ConvhullOptions',{});
    parser.addParameter('FacePointPreference','center');
	parser.parse(varargin{:});
	options = parser.Results;

    nDimension = size(convexHullPoint,2);
    if(nDimension==2 || nDimension==3)
        defaultConvhullOptions = {'Simplify',true};
    else
        defaultConvhullOptions = {};
    end
    [~,convhullOptions] = merge_name_value_pair_argument(defaultConvhullOptions,options.ConvhullOptions{:});

    measure = [];
    convexHullIndex = options.ConvexHullIndex;
    if(nDimension==1)
        % take maximum and minimum, then use as normals simply +1 for lower limit
        % and -1 for upper limit
        [pointLowerBoundary,iLowerBoundary] = min(convexHullPoint);
        [pointUpperBoundary,iUpperBoundary] = max(convexHullPoint);
        vectorNormalToFace = [1;-1];
        convexHullIndex = [iLowerBoundary;iUpperBoundary];
        measure = pointUpperBoundary - pointLowerBoundary;
        clear pointLowerBoundary iLowerBoundary pointUpperBoundary iUpperBoundary
    elseif(nDimension==2)
    	if(isempty(convexHullIndex))
    		[convexHullIndex,measure] = convhull(convexHullPoint,convhullOptions{:});
    	end

    	% find edges of convex hull
    	if(size(convexHullIndex,2)==1) % convex hull index was computed with convhull
    		nEdge = size(convexHullIndex,1)-1;
            convexHullIndexNew = nan(nEdge,2);
            for i=1:nEdge
                convexHullIndexNew(i,:) = [convexHullIndex(i),convexHullIndex(i+1)];
            end
            convexHullIndex = convexHullIndexNew;
            clear nEdge convexHullIndexNew
        end
        convexHullEdge = convexHullPoint(convexHullIndex(:,2),:) - convexHullPoint(convexHullIndex(:,1),:);

    	% vector with dot product = 0 --> normal to edge/face
	    vectorNormalToFace = [convexHullEdge(:,2),-convexHullEdge(:,1)];
        clear convexHullEdge
	elseif(nDimension==3)
		if(isempty(convexHullIndex))
    		[convexHullIndex,measure] = convhull(convexHullPoint,convhullOptions{:});
    	end

    	% find both edges for each triangle/plane, with reference being the first index
        verticeReference = convexHullPoint(convexHullIndex(:,1),:);
        vectice1 = convexHullPoint(convexHullIndex(:,2),:);
        vertice2 = convexHullPoint(convexHullIndex(:,3),:);
    	convexHullEdge1 = vectice1 - verticeReference;
		convexHullEdge2 = vertice2 - verticeReference;

		% for vector normal to plane, use the cross product of both edges
		vectorNormalToFace = cross(convexHullEdge1,convexHullEdge2,2);
        clear verticeReference vectice1 vertice2 convexHullEdge1 convexHullEdge2
    else
    	if(isempty(convexHullIndex))
    		[convexHullIndex,measure] = convhulln(convexHullPoint,convhullOptions{:});
    	end

    	% follow similar procedure to previous elements --> find a vector that is
    	% completely perpendicular (dot product = 0) to all the edges (starting from
    	% one reference point).
        nFace = size(convexHullIndex,1);
        vectorNormalToFace = nan(nFace,nDimension);
        for i=1:nFace
            verticeReference = convexHullPoint(convexHullIndex(i,1),:);
            verticeOther = convexHullPoint(convexHullIndex(i,2:end),:);
            convexHullEdge = verticeOther - verticeReference;

            % want to find x such that A*x=0, where A has the edges in each row
            %   --> that imples dot(edge1,x)=0, dot(edge2,x)=0, ...
            % to do so, find the nullspace of a matrix
    	    vectorNormalToFace(i,:) = null(convexHullEdge)';
        end
        clear nFace verticeReference verticeOther convexHullEdge
    end
    faceNormalIn = vectorNormalToFace./vecnorm(vectorNormalToFace,2,2);

    % get reference face point
    nFace = size(faceNormalIn,1);
    facePoint = nan(nFace,nDimension);
    if(strcmpi(options.FacePointPreference,'vertice'))
        facePoint = convexHullPoint(convexHullIndex(:,1),:);
    else % center
        for i=1:nFace
            facePoint(i,:) = mean(convexHullPoint(convexHullIndex(i,:),:),1);
        end
    end

    % make sure normal vectors are pointing inside
	convexHullCenter = mean(convexHullPoint,1);
	distanceCenterToVertice = facePoint - convexHullCenter;
	dotProduct = dot(distanceCenterToVertice,faceNormalIn,2);
	wrongOrientation = (dotProduct>0);
	faceNormalIn(wrongOrientation,:) = -faceNormalIn(wrongOrientation,:);
end