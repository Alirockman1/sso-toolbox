classdef CandidateSpaceAlphaShape < CandidateSpaceBase
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
        
        CandidateAlphaShape

        AlphaShapeOptions

        CriticalAlphaType
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
    	function obj = CandidateSpaceAlphaShape(designSpaceLowerBound,designSpaceUpperBound,varargin)
    		parser = inputParser;
            parser.addRequired('designSpaceLowerBound',@(x)isnumeric(x)&&(size(x,1)==1));
            parser.addRequired('designSpaceUpperBound',@(x)isnumeric(x)&&(size(x,1)==1));
            parser.addParameter('SamplingBoxSlack',0.5,@(x)isnumeric(x)&&isscalar(x)&&(x>=0)&&(x<=1));
            parser.addParameter('GrowthDistanceOptions',{});
            parser.addParameter('AlphaShapeOptions',{});
            parser.addParameter('CriticalAlphaType','one-region');
            parser.addParameter('FillHoles',true);
            parser.parse(designSpaceLowerBound,designSpaceUpperBound,varargin{:});

            
            obj.DesignSpaceLowerBound = parser.Results.designSpaceLowerBound;;
            obj.DesignSpaceUpperBound = parser.Results.designSpaceUpperBound;
            obj.SamplingBoxSlack = parser.Results.SamplingBoxSlack;
            obj.GrowthDistanceOptions = parser.Results.GrowthDistanceOptions;
            obj.CriticalAlphaType = parser.Results.CriticalAlphaType;

            if(~parser.Results.FillHoles)
                defaultAlphaShapeOptions = {};
            else
                defaultAlphaShapeOptions = {'HoleThreshold',prod(designSpaceUpperBound-designSpaceLowerBound)};
            end
            [~,obj.AlphaShapeOptions] = merge_name_value_pair_argument(...
                defaultAlphaShapeOptions,parser.Results.AlphaShapeOptions);

            obj.CandidateAlphaShape = [];
            obj.DesignSampleDefinition = [];
            obj.IsInsideDefinition = [];
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

            % boundary Definition
            obj.DesignSampleDefinition = designSample;
            obj.IsInsideDefinition = isInside;

            % create base alpha shape with all points
            obj.CandidateAlphaShape = alphaShape(designSample(isInside,:),obj.AlphaShapeOptions{:});

            % check for the possibilities where results may change
            possibleAlphaRadius = obj.CandidateAlphaShape.alphaSpectrum;

            % eliminate possibilites where not all points are included; 
            % additionally, flip so lower index means smaller volume
            minimumAlphaRadius = obj.CandidateAlphaShape.criticalAlpha(obj.CriticalAlphaType);
            possibleAlphaRadius = flip(possibleAlphaRadius(possibleAlphaRadius>=minimumAlphaRadius));

            iStart = 1;
            iEnd = size(possibleAlphaRadius,1);
            converged = false;
            while(~converged)
            	iMiddle = iStart + floor((iEnd-iStart)/2);

            	% find boundary with this shrink factor
            	alphaRadiusTest = possibleAlphaRadius(iMiddle);
            	obj.CandidateAlphaShape = alphaShape(designSample(isInside,:),alphaRadiusTest,obj.AlphaShapeOptions{:});

            	% see if any points that should be inside aren't
            	isInsideTest = obj.is_in_candidate_space(designSample(~isInside,:));
            	if(any(isInsideTest))
            		iEnd = iMiddle; % move upper boundary --> more restrictive
            	else
            		iStart = iMiddle; % move lower boundary --> less restrictive
            	end

            	converged = (iEnd-iStart<=1);
            end
            
            % test least restrictive option; if not valid, use previous
            alphaRadiusTest = possibleAlphaRadius(iEnd);
            obj.CandidateAlphaShape = alphaShape(designSample(isInside,:),alphaRadiusTest,obj.AlphaShapeOptions{:});
            isInsideTest = obj.is_in_candidate_space(designSample(~isInside,:));
            if(any(isInsideTest))
                alphaRadiusTest = possibleAlphaRadius(iStart);
                obj.CandidateAlphaShape = alphaShape(designSample(isInside,:),alphaRadiusTest,obj.AlphaShapeOptions{:});
            end


            % find which indexes belong in the boundary
            nSample = size(designSample,1);
            boundaryIndex = unique(reshape(obj.CandidateAlphaShape.boundaryFacets,[],1));
            globalBoundaryIndex = convert_index_base(isInside,boundaryIndex,'backward');
            obj.IsShapeDefinition = ismember((1:nSample)',globalBoundaryIndex);
    	end

    	function [label, score] = is_in_candidate_space(obj,designSample)
    		label = inShape(obj.CandidateAlphaShape,designSample);
    		[~,score] = nearestNeighbor(obj.CandidateAlphaShape,designSample);
    		score(label) = -score(label);
    	end

    	function plotHandle = plot_candidate_space(obj,figureHandle,varargin)
    		figure(figureHandle);
    		plotHandle = plot(obj.CandidateAlphaShape,varargin{:});
    	end

    	function measure = get.Measure(obj)
    		if(size(obj.DesignSampleDefinition,2)==2)
    			measure = obj.CandidateAlphaShape.area;
    		else
    			measure = obj.CandidateAlphaShape.volume;
    		end
    	end
    end
end