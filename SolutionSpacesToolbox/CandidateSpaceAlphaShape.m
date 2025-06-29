classdef CandidateSpaceAlphaShape < CandidateSpaceBase
%CANDIDATESPACEALPHASHAPE Candidate Space defined by delaunay triangulation
%   CANDIDATESPACEALPHASHAPE can be used to define candidate spaces for  
%   component solution space computation. In this case, the candidate space is 
%   defined mainly through an alphaShape object.
%   As this uses alphaShape as a base, this type of candidate space can only
%   be used for 2D and 3D problems.
%
%   This class is derived from CandidateSpaceBase.
%
%   CANDIDATESPACEALPHASHAPE properties:
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
%       - GrowthDistanceOptions : options to be used when computing the distance 
%       that the candidate space will grow to. In its default base 
%       implementation, these are options to 'knnsearch'.
%       - AlphaShapeOptions : additional options when training the alphaShape 
%       object. 
%       - CriticalAlphaType : criterion for the minimum value of alpha radius
%       allowed.
%   CANDIDATESPACEALPHASHAPE methods:
%       - define_candidate_space : create a candidate space based on design 
%       samples that are labeled as inside/outside.
%       - grow_candidate_space : expand the candidate space by a given factor.
%       - is_in_candidate_space : verify if given design samples are inside 
%       the candidate space.
%       - plot_candidate_space : visualize 1D/2D/3D candidate spaces in given
%       figure.
%
%   See also CandidateSpaceBase, alphaShape.
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
        %   In this case, these are the convex hull index points.
        %
        %   ISSHAPEDEFINITION : (nSample,1) logical
        %
        %   See also DesignSampleDefinition, IsInsideDefinition.
        IsShapeDefinition

        %GROWTHDISTANCEOPTIONS Options when computing distances in the growth algorithm
        %   GROWTHDISTANCEOPTIONS are the extra options to be used when finding the
        %   distance between sample points in the growth algorithm. 
        %   By default in its base implementation, these would be options for the 
        %   'knnsearch' funciton.
        %
        %   GROWTHDISTANCEOPTIONS : (1,nOption) cell
        %
        %   See also grow_candidate_space, knnsearch.
        GrowthDistanceOptions
        
        %ALPHASHAPEOBJECT alphaShape object to be used for candidate space definition
        %   ALPHASHAPEOBJECT contains the alphaShape which can be used to determine
        %   if query points are inside or outside the region.
        %
        %   ALPHASHAPEOBJECT : alphaShape object
        %
        %   See also alphaShape.
        AlphaShapeObject

        %ALPHASHAPEOPTIONS Options for when defining alphaShape
        %   ALPHASHAPEOPTIONS allows the specification of additional options when 
        %   computing alphaShape from the given points.
        %
        %   ALPHASHAPEOPTIONS : (1,nOptions) cell
        %
        %   See also alphaShape.
        AlphaShapeOptions

        %CRITICALALPHATYPE How to determine the allowed possible alpha radius
        %   CRITICALALPHATYPE determines the criterion for the minimum alpha radius to
        %   be chosen during the definition of the candidate space.
        %   'one-region' should be used for simply-connected candidate spaces, while
        %   'all-points' can be used for more general representations.
        %
        %   CRITICALALPHATYPE : char OR string
        %
        %   See also criticalAlpha.
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
        %CANDIDATESPACEALPHASHAPE Constructor
        %   CANDIDATESPACEALPHASHAPE is a constructor initializes an object of this 
        %   class.
        %   
        %   OBJ = CANDIDATESPACEALPHASHAPE(DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND)
        %   creates an object of the CandidateSpaceAlphaShape class and sets its design 
        %   space boundaries to its input values DESIGNSPACELOWERBOUND and 
        %   DESIGNSPACEUPPERBOUND. Other properties are set to empty.
        %
        %   OBJ = CANDIDATESPACEALPHASHAPE(...,NAME,VALUE,...) also allows one to set
        %   specific options for the object. This can be 
        %       - 'SamplingBoxSlack' : where the boundaries of the sampling box 
        %       will be relative to the strictest bounding box and the most relaxed
        %       bounding box. A value of 0 means no slack and therefore the sampling
        %       box will be the most strict one possible, and 1 means the sampling
        %       box will be the most relaxed.
        %       - 'AlphaShapeOptions' : additional options when training the alphaShape 
        %       object. By default, it is either empty, or 'HoleThreshold' is set to 
        %       fill all holes in the object, depending on the 'FillHoles' input.
        %       - 'CriticalAlphaType' : criterion for the minimum value of alpha radius
        %       allowed, as seen in 'alphaShape.criticalAlpha'. Default: 'one-region'.
        %       - 'FillHoles' : logical flag to choose whether holes should always be
        %       filled or not. This should be set to true for simply-connected candidate
        %       spaces, and false if that is not a necessary restriction. Default: true.
        %
        %   Inputs:
        %       - DESIGNSPACELOWERBOUND : (1,nDesignVariable) double
        %       - DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double
        %       - 'SamplingBoxSlack' : double
        %       - 'DelaunaynOptions' : (1,nOptionsDelaunayn) cell
        %       - 'SimplexConvexHullOptions' : (1.nOptionsConvexHull) cell
        %
        %   Outputs:
        %       - OBJ : CandidateSpaceAlphaShape
        %
        %   See also alphaShape, criticalAlpha.

    		parser = inputParser;
            parser.addRequired('designSpaceLowerBound',@(x)isnumeric(x)&&(size(x,1)==1));
            parser.addRequired('designSpaceUpperBound',@(x)isnumeric(x)&&(size(x,1)==1));
            parser.addParameter('SamplingBoxSlack',0.5,@(x)isnumeric(x)&&isscalar(x)&&(x>=0)&&(x<=1));
            parser.addParameter('AlphaShapeOptions',{});
            parser.addParameter('CriticalAlphaType','one-region');
            parser.addParameter('FillHoles',true);
            parser.addParameter('GrowthDistanceOptions',{});
            parser.parse(designSpaceLowerBound,designSpaceUpperBound,varargin{:});
            
            obj.DesignSpaceLowerBound = parser.Results.designSpaceLowerBound;;
            obj.DesignSpaceUpperBound = parser.Results.designSpaceUpperBound;
            obj.SamplingBoxSlack = parser.Results.SamplingBoxSlack;
            obj.CriticalAlphaType = parser.Results.CriticalAlphaType;
            obj.GrowthDistanceOptions = parser.Results.GrowthDistanceOptions;

            if(~parser.Results.FillHoles)
                defaultAlphaShapeOptions = {};
            else
                defaultAlphaShapeOptions = {'HoleThreshold',prod(designSpaceUpperBound-designSpaceLowerBound)};
            end
            [~,obj.AlphaShapeOptions] = merge_name_value_pair_argument(...
                defaultAlphaShapeOptions,parser.Results.AlphaShapeOptions);

            obj.AlphaShapeObject = [];
            obj.DesignSampleDefinition = [];
            obj.IsInsideDefinition = [];
    	end

    	function obj = define_candidate_space(obj,designSample,isInside)
        %DEFINE_CANDIDATE_SPACE Initial definition of the candidate space
        %   DEFINE_CANDIDATE_SPACE uses labeled design samples to define the inside / 
        %   outside regions of the candidate space. For CandidateSpaceAlphaShape, this
        %   creates an alphaShape object containing all the 'inside' points, while 
        %   attempting not to include any of the 'outside' ones. This is achieved by
        %   performing binary search on the possible values of alpha radius until the
        %   largest possible one fitting the criteria is found.
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
        %       - OBJ : CandidateSpaceAlphaShape
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %       - ISINSIDE : (nSample,1) logical
        %   
        %   Outputs:
        %       - OBJ : CandidateSpaceAlphaShape
        %   
        %   See also alphaShape, alphaSpectrum.

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

            % create base alpha shape with all points
            obj.AlphaShapeObject = alphaShape(designSample(isInside,:),obj.AlphaShapeOptions{:});

            % check for the possibilities where results may change
            possibleAlphaRadius = obj.AlphaShapeObject.alphaSpectrum;

            % eliminate possibilites where not all points are included; 
            % additionally, flip so lower index means smaller volume
            minimumAlphaRadius = obj.AlphaShapeObject.criticalAlpha(obj.CriticalAlphaType);
            possibleAlphaRadius = flip(possibleAlphaRadius(possibleAlphaRadius>=minimumAlphaRadius));

            % use binary search to find biggest alpha radius that still does not include
            % 'outside' designs inside it
            iStart = 1;
            iEnd = size(possibleAlphaRadius,1);
            converged = false;
            while(~converged)
            	iMiddle = iStart + floor((iEnd-iStart)/2);

            	% find boundary with this alpha radius
            	alphaRadiusTest = possibleAlphaRadius(iMiddle);
            	obj.AlphaShapeObject = alphaShape(designSample(isInside,:),alphaRadiusTest,obj.AlphaShapeOptions{:});

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
            obj.AlphaShapeObject = alphaShape(designSample(isInside,:),alphaRadiusTest,obj.AlphaShapeOptions{:});
            isInsideTest = obj.is_in_candidate_space(designSample(~isInside,:));
            if(any(isInsideTest))
                alphaRadiusTest = possibleAlphaRadius(iStart);
                obj.AlphaShapeObject = alphaShape(designSample(isInside,:),alphaRadiusTest,obj.AlphaShapeOptions{:});
            end

            % find which indexes belong in the boundary
            nSample = size(designSample,1);
            boundaryIndex = unique(reshape(obj.AlphaShapeObject.boundaryFacets,[],1));
            globalBoundaryIndex = convert_index_base(isInside,boundaryIndex,'backward');
            obj.IsShapeDefinition = ismember((1:nSample)',globalBoundaryIndex);
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
        %       - OBJ : CandidateSpaceAlphaShape
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %   
        %   Outputs:
        %       - ISINSIDE : (nSample,1) logical
        %       - SCORE : (nSample,1) double
        %   
        %   See also inShape, nearestNeighbor.

    		label = obj.AlphaShapeObject.inShape(designSample);
    		[~,score] = obj.AlphaShapeObject.nearestNeighbor(designSample);
    		score(label) = -score(label);
    	end

        function obj = grow_candidate_space(obj,growthRate)
        %GROW_CANDIDATE_SPACE Expansion of candidate space by given factor
        %   GROW_CANDIDATE_SPACE will grow the region considered inside the current 
        %   candidate space by the factor given. Said growth is done in a fixed rate 
        %   defined by the input relative to the design space.
        %
        %   OBJ = OBJ.GROW_CANDIDATE_SPACE(GROWTHRATE) will growth the candidate space 
        %   defined in OBJ by a factor of GROWTHRATE. This is an isotropic expansion of 
        %   the candidate space by a factor of the growth rate times the size of the 
        %   design space.
        %
        %   Inputs:
        %       - OBJ : CandidateSpaceBase
        %       - GROWTHRATE : double
        %   
        %   Outputs:
        %       - OBJ : CandidateSpaceBase
        %   
        %   See also is_in_candidate_space.

            [designSampleExpanded,labelExpanded] = grow_sample_region_positive_label(...
                obj.DesignSampleDefinition,...
                obj.IsInsideDefinition,...
                obj.DesignSpaceLowerBound,...
                obj.DesignSpaceUpperBound,...
                growthRate,...
                obj.GrowthDistanceOptions{:});
            
            % train new candidate space
            obj = obj.define_candidate_space(designSampleExpanded,labelExpanded);
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
        %   refer to 'plot' of alphaShape objects.
        %
        %   Input:
        %       - OBJ : CandidateSpaceAlphaShape
        %       - FIGUREHANDLE : Figure
        %
        %   Output:
        %       - PLOTHANDLE : patch object
        %
        %   See also alphaShape.plot.

    		figure(figureHandle);
    		plotHandle = obj.AlphaShapeObject.plot(varargin{:});
    	end

    	function measure = get.Measure(obj)
    		if(size(obj.DesignSampleDefinition,2)==2)
    			measure = obj.AlphaShapeObject.area;
    		else
    			measure = obj.AlphaShapeObject.volume;
    		end
    	end
    end
end