function [freeBoundaryFacet,freeBoundaryVertex] = find_free_boundary_facets(delaunayIndex)
%FIND_FREE_BOUNDARY_FACETS Find boundary facets of delaunay triangulation
%   FIND_FREE_BOUNDARY_FACETS finds the facets of a triangulation which are 
%   boundary facets; this is determined by checking how many simpleces from the
%   triangulation contain each facet, and if it is only one, than that facet
%   is a free boundary facet. This is a generalization for any number of 
%   dimensions of 'freeBoundary'.
%
%   FREEBOUNDARYFACET = FIND_FREE_BOUNDARY_FACETS(DELAUNAYINDEX) receives the 
%   delaunay triangulation index DELAUNAYINDEX and returns the facet indices 
%   that are on the boundary FREEBOUNDARYFACET.
%
%   [FREEBOUNDARYFACET,FREEBOUNDARYVERTEX] = FIND_FREE_BOUNDARY_FACETS(...) 
%   additionally returns the indices of vertices that are on those said facet
%   boundaries FREEBOUNDARYVERTEX.
%
%   Inputs:
%       - DELAUNAYINDEX : (nSimplex,nSimplexVertex) integer
%
%   Output:
%       - FREEBOUNDARYFACET : (nFreeBoundaryFacet,nFacetVertex) integer
%       - FREEBOUNDARYVERTEX : (nVertexUnique,1) integer
%
%   See also: delaunay, delaunayn, CandidateSpaceDelaunay, freeBoundary.
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


    nSimplex = size(delaunayIndex,1);
    nVertexPerSimplex = size(delaunayIndex,2);
    nDimension = nVertexPerSimplex - 1;
    nFacetPerSimplex = nDimension + 1;
    nVertexPerFacet = nVertexPerSimplex - 1;

    % first, build indices of possible local combinations of vertices for each facet
    verticesLocalIndex = 1:nVertexPerSimplex;
    simplexLocalFacet = nan(nFacetPerSimplex,nVertexPerFacet);
    for i=1:nFacetPerSimplex
        facetIndex = verticesLocalIndex;
        facetIndex(i) = [];
        simplexLocalFacet(i,:) = facetIndex;
    end
    clear verticesLocalIndex facetIndex nVertexPerSimplex

    % second, build faces based on index
    verticesGlobalIndex = sort(delaunayIndex,2);
    triangulationFacet = nan(nSimplex*nFacetPerSimplex,nDimension);
    for i = 1:nSimplex
        simplexFacet = reshape(verticesGlobalIndex(i,simplexLocalFacet),nFacetPerSimplex,nVertexPerFacet);
        simplexIndex = 1 + (nDimension+1)*(i-1) + [0:nDimension];
        triangulationFacet(simplexIndex,:) = simplexFacet;
    end
    clear verticesGlobalIndex simplexFacet simplexIndex nDimension nSimplex nFacetPerSimplex nVertexPerFacet

    % third, find unique facets and in how many simpleces they appear
    [uniqueFacet, ~, iUniqueToTriangulation] = unique(triangulationFacet, 'rows');
    nSimplexWithFacet = accumarray(iUniqueToTriangulation,1);

    % finally, facets that appear only once are free-boundary facets
    freeBoundaryFacet = uniqueFacet(nSimplexWithFacet==1, :);

    if(nargout>1)
        % if necessary, also find all vertex indices that are part of these boundaries
        freeBoundaryVertex = unique(freeBoundaryFacet(:));
    end
end
