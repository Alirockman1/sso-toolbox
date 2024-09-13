classdef CandidateSpaceDelaunay < CandidateSpaceBase
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
        %   In this case, these are the convex hull index points.
        %
        %   ISSHAPEDEFINITION : (nSample,1) logical
        %
        %   See also DesignSampleDefinition, IsInsideDefinition.
        IsShapeDefinition
        
        DelaunayIndex

        DelaunaySimplex

        DelaunaynOptions

        SimplexConvexHullOptions
    end

    properties (SetAccess = protected, Dependent)
    	%MEASURE Size measure of the candidate space
        %   MEASURE is a value that works as the measure of the candidate space. This 
        %   may be its volume, or normalized volume relative to its design space, or
        %   some other metric.
        %   In this case, this is the area/volume of the convex hull as computed.
        %
        %   MEASURE : double
        %
        %   See also DesignSampleDefinition, IsInsideDefinition, convhull, convhulln.   
        Measure
    end

    methods
    	function obj = CandidateSpaceDelaunay(designSpaceLowerBound,designSpaceUpperBound,varargin)
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

            obj.GrowthDistanceOptions = {};
            obj.DelaunayIndex = [];
            obj.DelaunaySimplex = struct(...
                'Vertices',[],...
            	'HullIndex',[],...
            	'Measure',[],...
            	'FacePoint',[],...
            	'FaceNormal',[]);
        end

        function obj = define_candidate_space(obj,designSample,isInside)
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

            if(any(~isInside))
                % find which simplices have 'outside' designs inside
                outsideInsideSimplex = tsearchn(insideSample,obj.DelaunayIndex,designSample(~isInside,:));
                isInsideWrong = unique(outsideInsideSimplex(~isnan(outsideInsideSimplex)));
                boundarySimplex = unique(reshape(obj.DelaunayIndex(isInsideWrong,:),[],1));
                
                % shape is defined by these designs which force the exclusion / are excluded
                obj.IsShapeDefinition = false(size(designSample,1),1);
                obj.IsShapeDefinition(convert_index_base(~isInside,~isnan(outsideInsideSimplex),'backward')) = true;
                obj.IsShapeDefinition(convert_index_base(isInside,boundarySimplex,'backward')) = true;
                obj.IsShapeDefinition = obj.IsShapeDefinition | ...
                    design_find_boundary_samples(obj.DesignSampleDefinition,obj.IsInsideDefinition);
    
                % delete simplices where bad designs are inside
                obj.DelaunayIndex(isInsideWrong,:) = [];
            else
                obj.IsShapeDefinition = true(size(designSample,1),1);
            end

            % treat each simplex as a convex hull
            nSimplex = size(obj.DelaunayIndex,1);
            for i=1:nSimplex
            	simplexPoint = insideSample(obj.DelaunayIndex(i,:),:);
            	[convexHullIndex,measure,convexHullFacePoint,convexHullFaceNormal] = ...
            		convex_hull_face(simplexPoint,obj.SimplexConvexHullOptions{:});

            	obj.DelaunaySimplex(i) = struct(...
            		'Vertices',simplexPoint,...
            		'HullIndex',convexHullIndex,...
            		'Measure',measure,...
            		'FacePoint',convexHullFacePoint,...
            		'FaceNormal',convexHullFaceNormal);
            end
        end

        function obj = grow_candidate_space(obj,growthRate)
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
                
                [convexHullIndex,measure,convexHullFacePoint,convexHullFaceNormal] = ...
                    convex_hull_face(verticesNew,obj.SimplexConvexHullOptions{:});

                obj.DelaunaySimplex(i) = struct(...
                    'Vertices',verticesNew,...
                    'HullIndex',convexHullIndex,...
                    'Measure',measure,...
                    'FacePoint',convexHullFacePoint,...
                    'FaceNormal',convexHullFaceNormal);

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
            [insideSimplex,barycentricCoordinate] = tsearchn(obj.ActiveDesign,obj.DelaunayIndex,designSample);
            isInsideSpace = ~isnan(insideSimplex);
            isInBoundary = ismember(designSample,obj.ActiveDesign,'rows');
            label = isInsideSpace | isInBoundary;

            if(nargout>1)
        	    nSample = size(designSample,1);
        	    nSimplex = size(obj.DelaunayIndex,1);
        	    scoreSimplex = nan(nSample,nSimplex);
                scoreSimplex(isInBoundary,1) = 0;
                scoreSimplex(isInsideSpace,1) = max(barycentricCoordinate(isInsideSpace,:)-1,[],2);
        	    
                outsideSample = designSample(~label,:);
        	    for i=1:nSimplex
        		    [~,score] = is_in_convex_hull_with_face(...
        			    obj.DelaunaySimplex(i).FacePoint,...
        			    obj.DelaunaySimplex(i).FaceNormal,...
        			    outsideSample);
				    scoreSimplex(~label,i) = score;
                end
        	    score = min(scoreSimplex,[],2);
            end
        end

        function plotHandle = plot_candidate_space(obj,figureHandle,varargin)
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