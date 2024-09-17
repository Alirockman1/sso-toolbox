function [freeBoundaryFacet,freeBoundaryVertex] = find_free_boundary_facets(delaunayIndex)

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
