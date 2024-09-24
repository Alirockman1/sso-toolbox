classdef CandidateSpaceDelaunay < CandidateSpaceBase
%CANDIDATESPACEDELAUNAY Candidate Space defined by delaunay triangulation
%   CANDIDATESPACEDELAUNAY can be used to define candidate spaces for  
%   component solution space computation. In this case, the candidate space is 
%   defined by the result of a Delaunay triangulation where the triangulation
%   is done with all the designs considered 'inside' and where the simpleces
%   that contain designs considered 'outside' are removed.
%
%   This class is derived from CandidateSpaceBase.
%
%   CANDIDATESPACEDELAUNAY properties:
%       - DesignSampleDefinition : design sample points used in the candidate 
%       space definition.
%       - IsInsideDefinition : logical labels of the design sample points used   
%       in the candidate space definition regarding whether they are inside or 
%       outside the candidate space (true = inside, false = outside).
%       - IsShapeDefinition : logical labels of the design sample points used in 
%       the candidate space definition regarding whether these designs are 
%       directly related to the shape of the candidate space (or in other words,
%       the design sample points that define the boundary itself and are used
%       to identify if a design is inside/outside this space).
%       - ActiveDesign : sample points that are labeled as inside.
%       - Measure : measure of the candidate space (usually area/volume/etc).
%       - SamplingBox : bounding box around the internal region of the candidate
%       space that can be used to better attempt to sample inside said region.
%       - SamplingBoxSlack : a value used to give the defined amount of slack
%       for the SamplingBox, varying between the strict bounding box around the 
%       internal region and a larger bounding box which may contain negative
%       designs in its edges. 
%       - DelaunayIndex : resulting indexing of the Delaunay triangulation after
%       removing the simpleces with 'outside' designs in them.
%       - DelaunaySimplex : all the simpleces which come as a result of the
%       triangulation, including their vertices, hull index, and measure.
%       - DelaunaynOptions : options for when using 'delaunayn' to compute
%       the triangulation.
%       - SimplexConvexHullOptions : options for when using 
%       'compute_convex_hull' to compute the properties of individual simpleces.
%
%   CANDIDATESPACEDELAUNAY methods:
%       - define_candidate_space : create a candidate space based on design 
%       samples that are labeled as inside/outside.
%       - grow_candidate_space : expand the candidate space by a given factor.
%       - is_in_candidate_space : verify if given design samples are inside 
%       the candidate space.
%       - plot_candidate_space : visualize 1D/2D/3D candidate spaces in given
%       figure.
%
%   See also CandidateSpaceBase.
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

	properties (SetAccess = protected)
        %DESIGNSAMPLEDEFINITION Design sample points used in candidate space definition
        %   DESIGNSAMPLEDEFINITION are the sample points used in the definition of the 
        %   current candidate space.
        %
        %   DESIGNSAMPLEDEFINITION : (nSample,nDesignVariable) double
        %
        %   See also IsInsideDefinition, IsShapeDefinition.
        DesignSampleDefinition
        
        %ISINSIDEDEFINITION Labels of sample points used in candidate space definition
        %   ISINSIDEDEFINITION are the labels of design samples used in the definition 
        %   of the current candidate space. A label of 'true' indicates the respective 
        %   design point is inside the candidate space.
        %
        %   ISINSIDEDEFINITION : (nSample,1) logical
        %
        %   See also DesignSampleDefinition, IsShapeDefinition.
        IsInsideDefinition

        %ISSHAPEDEFINITION Labels if sample points from definition contributes to shape
        %   ISSHAPEDEFINITION is a logical array where 'true' values indicate that that
        %   design point (in the respective row) from the definition sample actively 
        %   contributes to the shape of the candidate space. 
        %   In this case, these are the free boundary facets index points.
        %
        %   ISSHAPEDEFINITION : (nSample,1) logical
        %
        %   See also DesignSampleDefinition, IsInsideDefinition, 
        %   find_free_boundary_facets.
        IsShapeDefinition
        
        %DELAUNAYINDEX Simplex vertices index resulting from Delaunay triangulation
        %   DELAUNAYINDEX is the tesselation index given by 'delaunayn'.
        %
        %   DELAUNAYINDEX : (nSimplex,nDimension+1) integer
        %   
        %   See also delaunayn.
        DelaunayIndex

        %DELAUNAYSIMPLEX Simpleces information (each regarded as a convex hull)
        %   DELAUNAYSIMPLEX saves the relevant information of each simplex individually,
        %   where they are regarded as separate convex hulls.
        %
        %   DELAUNAYSIMPLEX : (nSimplex,1) structure
        %       - Vertices : (nVertexPerSimplex,nDimension+1) double
        %       - HullIndex : (nFacePerSimplex,nDimension) double
        %       - Measure : double
        %
        %   See also delaunayn, compute_convex_hull.
        DelaunaySimplex

        %DELAUNAYNOPTIONS Options when computing 'delaunayn'
        %   DELAUNAYNOPTIONS are extra options which may be used when computing the 
        %   delaunay tesselation using 'delaunayn'.
        %
        %   DELAUNAYNOPTIONS : (1,nOptions) cell
        %
        %   See also delaunayn.
        DelaunaynOptions
        
        %SIMPLEXCONVEXHULLOPTIONS Options when computing convex hull of each simplex
        %   SIMPLEXCONVEXHULLOPTIONS are extra options which may be used when computing  
        %   the convex hull of each simplex using 'compute_convex_hull'.
        %
        %   SIMPLEXCONVEXHULLOPTIONS : (1,nOptions) cell
        %   
        %   See also compute_convex_hull.
        SimplexConvexHullOptions
    end

    properties (SetAccess = protected, Dependent)
    	%MEASURE Size measure of the candidate space
        %   MEASURE is a value that works as the measure of the candidate space. This 
        %   may be its volume, or normalized volume relative to its design space, or
        %   some other metric.
        %   In this case, this is the sum of areas/volumes/... of the simplices as 
        %   computed.
        %
        %   MEASURE : double
        %
        %   See also DesignSampleDefinition, IsInsideDefinition, convhull, convhulln.   
        Measure
    end

    methods
    	function obj = CandidateSpaceDelaunay(designSpaceLowerBound,designSpaceUpperBound,varargin)
        %CANDIDATESPACEDELAUNAY Constructor
        %   CANDIDATESPACEDELAUNAY is a constructor initializes an object of this class.
        %   
        %   OBJ = CANDIDATESPACEDELAUNAY(DESIGNSPACELOWERBOUND,
        %   DESIGNSPACEUPPERBOUND) creates an object of the CandidateSpaceDelaunay
        %   class and sets its design space boundaries to its input values 
        %   DESIGNSPACELOWERBOUND and DESIGNSPACEUPPERBOUND. Other properties are set to
        %   empty.
        %
        %   OBJ = CANDIDATESPACEDELAUNAY(...,NAME,VALUE,...) also allows one to set
        %   specific options for the object. This can be 
        %       - 'SamplingBoxSlack' : where the boundaries of the sampling box 
        %       will be relative to the strictest bounding box and the most relaxed
        %       bounding box. A value of 0 means no slack and therefore the sampling
        %       box will be the most strict one possible, and 1 means the sampling
        %       box will be the most relaxed.
        %       - 'DelaunaynOptions' : options for the computation of the Delaunay 
        %       tesselation using 'delaunayn'. Default is empty.
        %       - 'SimplexConvexHullOptions' : options for the computation of the convex
        %       hull of each simplex using 'compute_convex_hull'. Default is empty.
        %
        %   Inputs:
        %       - DESIGNSPACELOWERBOUND : (1,nDesignVariable) double
        %       - DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double
        %       - 'SamplingBoxSlack' : double
        %       - 'DelaunaynOptions' : (1,nOptionsDelaunayn) cell
        %       - 'SimplexConvexHullOptions' : (1.nOptionsConvexHull) cell
        %
        %   Outputs:
        %       - OBJ : CandidateSpaceDelaunay
        %
        %   See also delaunayn, compute_convex_hull.

            % parse inputs
            parser = inputParser;
            parser.addRequired('designSpaceLowerBound',@(x)isnumeric(x)&&(size(x,1)==1));
            parser.addRequired('designSpaceUpperBound',@(x)isnumeric(x)&&(size(x,1)==1));
            parser.addParameter('SamplingBoxSlack',0.5,@(x)isnumeric(x)&&isscalar(x)&&(x>=0)&&(x<=1));
            parser.addParameter('DelaunaynOptions',{});
            parser.addParameter('SimplexConvexHullOptions',{});
            parser.parse(designSpaceLowerBound,designSpaceUpperBound,varargin{:});
            
            obj.DesignSpaceLowerBound = parser.Results.designSpaceLowerBound;
            obj.DesignSpaceUpperBound = parser.Results.designSpaceUpperBound;
            obj.SamplingBoxSlack = parser.Results.SamplingBoxSlack;
            obj.DelaunaynOptions = parser.Results.DelaunaynOptions;
            obj.SimplexConvexHullOptions = parser.Results.SimplexConvexHullOptions;

            obj.DelaunayIndex = [];
            obj.DelaunaySimplex = struct(...
                'Vertices',[],...
            	'HullIndex',[],...
            	'Measure',[]);
        end

        function obj = define_candidate_space(obj,designSample,isInside)
        %DEFINE_CANDIDATE_SPACE Initial definition of the candidate space
        %   DEFINE_CANDIDATE_SPACE uses labeled design samples to define the inside / 
        %   outside regions of the candidate space. For CandidateSpaceDelaunay, this
        %   creates a tesselation with all designs labeled as 'inside', and removes any
        %   simplices which contain designs labeled 'outside'.
        %
        %   OBJ = OBJ.DEFINE_CANDIDATE_SPACE(DESIGNSAMPLE) receives the design samle
        %   points in DESIGNSAMPLE and returns a candidate space object OBJ with the new
        %   definition, assuming all designs are inside the candidate space.
        %
        %   OBJ = OBJ.DEFINE_CANDIDATE_SPACE(DESIGNSAMPLE,ISINSIDE) additionally 
        %   receives the inside/outside (true/false) labels of each design point in 
        %   ISINSIDE.
        %
        %   Inputs:
        %       - OBJ : CandidateSpaceDelaunay
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %       - ISINSIDE : (nSample,1) logical
        %   
        %   Outputs:
        %       - OBJ : CandidateSpaceDelaunay
        %   
        %   See also delaunayn, compute_convex_hull, find_triangulation_facets.

        	if(nargin<3)
                isInside = true(size(designSample,1),1);
            end
            
            if(any(designSample<obj.DesignSpaceLowerBound,'all') || ...
                any(designSample>obj.DesignSpaceUpperBound,'all'))
                error('CandidateSpaceConvexHull:Definition:OutsideDesignSpace',...
                    'Specified designs are outside the design space');
            end

            obj.DesignSampleDefinition = designSample;
            obj.IsInsideDefinition = isInside;
            insideSample = designSample(isInside,:);

            % create base delauney triangulation
            obj.DelaunayIndex = delaunayn(insideSample,obj.DelaunaynOptions{:});

            % remove all simpleces with 'outside' designs inside them
            if(any(~isInside))
                % find which simplices have 'outside' designs inside
                outsideInsideSimplex = tsearchn(insideSample,obj.DelaunayIndex,designSample(~isInside,:));
                isInsideWrong = unique(outsideInsideSimplex(~isnan(outsideInsideSimplex)));
    
                % delete simplices where bad designs are inside
                obj.DelaunayIndex(isInsideWrong,:) = [];
            end

            % get free facets to also define shape
            obj.IsShapeDefinition = false(size(designSample,1),1);
            [~,~,freeBoundaryVertex] = find_triangulation_facets(obj.DelaunayIndex);
            obj.IsShapeDefinition(convert_index_base(isInside,freeBoundaryVertex,'backward')) = true;

            % treat each simplex as a convex hull
            nSimplex = size(obj.DelaunayIndex,1);
            for i=1:nSimplex
            	simplexPoint = insideSample(obj.DelaunayIndex(i,:),:);
            	[convexHullIndex,measure] = compute_convex_hull(simplexPoint,obj.SimplexConvexHullOptions{:});

            	obj.DelaunaySimplex(i) = struct(...
            		'Vertices',simplexPoint,...
            		'HullIndex',convexHullIndex,...
            		'Measure',measure);
            end
        end

        function obj = grow_candidate_space(obj,growthRate)
        %GROW_CANDIDATE_SPACE Expansion of candidate space by given factor
        %   GROW_CANDIDATE_SPACE will grow the region considered inside the current 
        %   candidate space by the factor given. Said growth is done in a fixed rate 
        %   defined by the input relative to the design space.
        %   This is done by finding the center of each simplex that forms the 
        %   tesselation and expanding its vertices opposite to that; this is similar to  
        %   same process with the convex hull, but each simplex is treated separately.
        %
        %   OBJ = OBJ.GROW_CANDIDATE_SPACE(GROWTHRATE) will growth the candidate space 
        %   defined in OBJ by a factor of GROWTHRATE. This is an isotropic expansion of 
        %   the candidate space by a factor of the growth rate times the size of the 
        %   design space.
        %
        %   Inputs:
        %       - OBJ : CandidateSpaceDelaunay
        %       - GROWTHRATE : double
        %   
        %   Outputs:
        %       - OBJ : CandidateSpaceDelaunay
        %   
        %   See also define_candidate_space, is_in_candidate_space.

            nSimplex = size(obj.DelaunayIndex,1);
            nDimension = size(obj.DesignSampleDefinition,2);
            designSampleNew = nan(nSimplex*(nDimension+1),nDimension);
            delauynayIndexNew = nan(nSimplex,nDimension+1);
            designSpaceFactor = obj.DesignSpaceUpperBound - obj.DesignSpaceLowerBound;
            designSpace = [obj.DesignSpaceLowerBound;obj.DesignSpaceUpperBound];
        
            for i=1:nSimplex
                simplexVertices = obj.DelaunaySimplex(i).Vertices;
                simplexCenter = mean(simplexVertices,1);
                distanceVertice = simplexVertices - simplexCenter;
                directionGrowth = distanceVertice./vecnorm(distanceVertice,2,2);

                % find maximum growth rate not to escape design space
                maxGrowthRate = region_limit_line_search([],simplexVertices,designSpaceFactor.*directionGrowth,designSpace);
                verticeGrowthRate = min(growthRate,maxGrowthRate);

                % grow
                verticesNew = simplexVertices + verticeGrowthRate.*designSpaceFactor.*directionGrowth;
                verticesNew = max(verticesNew, obj.DesignSpaceLowerBound); % lower bound limit
                verticesNew = min(verticesNew, obj.DesignSpaceUpperBound); % upper bound limit
                verticesNew = unique(verticesNew,'rows');
                
                [convexHullIndex,measure] = compute_convex_hull(verticesNew,obj.SimplexConvexHullOptions{:});

                obj.DelaunaySimplex(i) = struct(...
                    'Vertices',verticesNew,...
                    'HullIndex',convexHullIndex,...
                    'Measure',measure);

                simplexIndex = 1 + (nDimension+1)*(i-1) + [0:nDimension];
                designSampleNew(simplexIndex,:) = verticesNew;
                delauynayIndexNew(i,:) = simplexIndex;
            end
            obj.DesignSampleDefinition = designSampleNew;
            obj.DelaunayIndex = delauynayIndexNew;

            nSample = size(designSampleNew,1);
            obj.IsInsideDefinition = true(nSample,1);
            obj.IsShapeDefinition = true(nSample,1);
        end

        function [label, score] = is_in_candidate_space(obj,designSample)
        %IS_IN_CANDIDATE_SPACE Verification if given design samples are inside
        %   IS_IN_CANDIDATE_SPACE uses the currently defined candidate space to 
        %   determine if given design sample points are inside or outside the candidate 
        %   space.
        %
        %   ISINSIDE = OBJ.IS_IN_CANDIDATE_SPACE(DESIGNSAMPLE) receives the design
        %   samples in DESIGNSAMPLE and returns whether or not they are inside the 
        %   candidate space in ISINSIDE. For ISINSIDE values of 'true', it means the 
        %   respective design is inside the candidate space, while 'false' means it is 
        %   outside.
        %
        %   [ISINSIDE,SCORE] = OBJ.IS_IN_CANDIDATE_SPACE(...) also returns a SCORE value
        %   for each sample point; negative values of SCORE indicate the design sample 
        %   is inside the candidate space, and positive values indicate it is outside. 
        %   Designs with lower/higher SCORE are further from the boundary, with 0 
        %   representing that they are exactly at the boundary.
        %
        %   Inputs:
        %       - OBJ : CandidateSpaceDelaunay
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %
        %   Outputs:
        %       - ISINSIDE : (nSample,1) logical
        %       - SCORE : (nSample,1) double
        %   
        %   See also tsearchn.

            [insideSimplex,barycentricCoordinate] = tsearchn(obj.ActiveDesign,obj.DelaunayIndex,designSample);
            isInsideSpace = ~isnan(insideSimplex);
            isInBoundary = ismember(designSample,obj.DesignSampleDefinition(obj,IsShapeDefinition,:),'rows');
            label = isInsideSpace | isInBoundary;

            if(nargout>1)
        	    nSample = size(designSample,1);;
        	    score = nan(nSample,1);
                score(isInBoundary) = 0;
                score(isInsideSpace) = max(barycentricCoordinate(isInsideSpace,:)-1,[],2);
                %TODO: barycentric coordinate perhaps not most appropriate, as that only tells
                %how 'inside' that is of the particular simplex, not the tesselation as a whole
                score(~label) = 1;
            end
        end

        function plotHandle = plot_candidate_space(obj,figureHandle,varargin)
        %PLOT_CANDIDATE_SPACE Visualization of the boundary of the canidate space 2D/3D
        %   PLOT_CANDIDATE_SPACE allows for the visualization of the boundary of the
        %   candidate space in the given figure. 
        %
        %   PLOTHANDLE = OBJ.PLOT_CANDIDATE_SPACE(FIGUREHANDLE) plots the boundary of
        %   the candidate space in figure FIGUREHANDLE, returning the handle of the 
        %   object plot PLOTHANDLE.
        %
        %   PLOTHANDLE = OBJ.PLOT_CANDIDATE_SPACE(...,NAME,VALUE) allows the 
        %   specification for additional options in the process. These options should 
        %   refer to 'patch'.
        %
        %   Input:
        %       - OBJ : CandidateSpaceDelaunay
        %       - FIGUREHANDLE : Figure
        %
        %   Output:
        %       - PLOTHANDLE : patch object
        %
        %   See also patch.

        	figure(figureHandle);
            nDimension = size(obj.DesignSampleDefinition,2);
        	if((nDimension==2) || (nDimension==3))
        		plotHandle = patch(...
        			'Faces',obj.DelaunayIndex,...
        			'Vertices',obj.ActiveDesign,...
        			varargin{:});
        	end
        end

        function measure = get.Measure(obj)
        	measure = 0;
        	nSimplex = size(obj.DelaunayIndex,1);
        	for i=1:nSimplex
        		measure = measure + obj.DelaunaySimplex(i).Measure;
        	end
        end
    end
end