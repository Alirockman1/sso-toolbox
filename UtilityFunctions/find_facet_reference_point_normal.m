function [facetReferencePoint,facetNormal] = find_facet_reference_point_normal(facetVertex,facetIndex,varargin)
%FIND_FACET_REFERENCE_POINT_NORMAL Express planes using one point/normal vector
%   FIND_FACET_REFERENCE_POINT_NORMAL takes facets represented by given vertices
%   and indices and finds for each facet a point that lies within that 
%   (hyper)plane, and a vector which is normal to said (hyper)plane/facet.
%
%   FACETREFERENCEPOINT = FIND_FACET_REFERENCE_POINT_NORMAL(FACETVERTEX,
%   FACETINDEX) receives sample points FACETVERTEX and the index which 
%   determines which ones belong to which facet FACETINDEX, and returns points 
%   where each one belongs to the respective facet FACETREFERENCEPOINT.
%
%   [FACETREFERENCEPOINT,FACETNORMAL] = FIND_FACET_REFERENCE_POINT_NORMAL(...)
%   additionally returns for each facet a normalized vector which is normal to  
%   it FACETNORMAL.
%
%   [...] = FIND_FACET_REFERENCE_POINT_NORMAL(...,NAME,VALUE,...) allows for the
%   specification of additional options. These are:
%       - 'FacePointPreference' : preference regarding which point to choose
%       as the reference FACEPOINT. There are two possibilities:
%           -- 'center' : center of the face.
%           -- 'vertice' : first vertice in the convex hull index.
%
%   Input:
%       - FACETVERTEX : (nSample,nDimension) double
%       - FACETINDEX : (nFacet,nDimension) integer
%       - 'FacePointPreference' : char OR string
%
%   Output:
%       - FACETREFERENCEPOINT : (nFacet,nDimension) double
%       - FACETNORMAL : (nFacet,nDimension) double
%
%   See also compute_convex_hull, is_in_convex_hull_with_facet_normal, 
%   find_triangulation_facets.
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
    parser.addParameter('FacePointPreference','center');
	parser.parse(varargin{:});
	options = parser.Results;

    nDimension = size(facetVertex,2);
    nFacet = size(facetIndex,1);

    % get reference face point
    facetReferencePoint = nan(nFacet,nDimension);
    if(strcmpi(options.FacePointPreference,'vertex'))
        facetReferencePoint = facetVertex(facetIndex(:,1),:);
    else % center
        for i=1:nFacet
            facetReferencePoint(i,:) = mean(facetVertex(facetIndex(i,:),:),1);
        end
    end

    % get normals
    if(nargout>1)
        if(nDimension==1)
            vectorNormalToFacet = ones(nFacet,1);
        elseif(nDimension==2)
        	% find edges of shape
            verticeReference = facetVertex(facetIndex(:,1),:);
            verticeNext = facetVertex(facetIndex(:,2),:);
            edge = verticeNext - verticeReference;

        	% vector with dot product = 0 --> normal to edge/face
    	    vectorNormalToFacet = [edge(:,2),-edge(:,1)];
    	elseif(nDimension==3)
        	% find both edges for each triangle/plane, with reference being the first index
            verticeReference = facetVertex(facetIndex(:,1),:);
            vectice1 = facetVertex(facetIndex(:,2),:);
            vertice2 = facetVertex(facetIndex(:,3),:);
        	edge1 = vectice1 - verticeReference;
    		edge2 = vertice2 - verticeReference;

    		% for vector normal to plane, use the cross product of both edges
    		vectorNormalToFacet = cross(edge1,edge2,2);
        else
        	% follow similar procedure to previous elements --> find a vector that is
        	% completely perpendicular (dot product = 0) to all the edges (starting from
        	% one reference point).
            vectorNormalToFacet = nan(nFacet,nDimension);
            for i=1:nFacet
                verticeReference = facetVertex(facetIndex(i,1),:);
                verticeOther = facetVertex(facetIndex(i,2:end),:);
                edge = verticeOther - verticeReference;

                % want to find x such that A*x=0, where A has the edges in each row
                %   --> that imples dot(edge1,x)=0, dot(edge2,x)=0, ...
                % to do so, find the nullspace of a matrix
        	    vectorNormalToFacet(i,:) = null(edge)';
            end
        end
        facetNormal = vectorNormalToFacet./vecnorm(vectorNormalToFacet,2,2);
    end
end