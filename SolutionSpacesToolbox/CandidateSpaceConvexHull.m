classdef CandidateSpaceConvexHull < CandidateSpaceBase
%CANDIDATESPACECONVEXHULL Candidate Space defined by a convex hull
%   CANDIDATESPACECONVEXHULL can be used to define candidate spaces for  
%   component solution space computation. In this case, the candidate space is 
%   defined by a convex hull that encompasses all inside designs.
%   Note that if the trimming operation does not result in convex-shaped 
%   solutions, this candidate space definition may also include "outside" 
%   designs inside of it.
%
%   This class is derived from CandidateSpaceBase.
%
%   CANDIDATESPACECONVEXHULL properties:
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
%       - ConvexHullIndex : indexing for the definition of the convex hull.
%
%   CANDIDATESPACECONVEXHULL methods:
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
        %   In this case, these are the convex hull index points.
        %
        %   ISSHAPEDEFINITION : (nSample,1) logical
        %
        %   See also DesignSampleDefinition, IsInsideDefinition.
        IsShapeDefinition

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

        %CONVEXHULLINDEX Indexing for the definition of the convex hull
        %   CONVEXHULLINDEX is the result of the process to create a convex
        %   hull, indexing the vertices that define all the simplexes.
        %
        %   CONVEXHULLINDEX : (nHullFace,nDesignVariable) integer
        %
        %   See also convhull, convhulln.   
        ConvexHullIndex

        %CONVEXHULLFACEPOINT Reference point of each face in the convex hull
        %   CONVEXHULLFACEPOINT contains reference points in each face of the convex 
        %   hull, which can be used (together with the face normals) to determine if 
        %   given points are inside the convex hull or not.
        %
        %   CONVEXHULLFACEPOINT : (nHullFace,nDesignVariable) double
        %
        %   See also convex_hull_plane, is_in_convex_hull_with_plane.   
        ConvexHullFacePoint

        %CONVEXHULLFACENORMAL Face-normal vectors pointing to the interior of the hull
        %   CONVEXHULLFACENORMAL contains vectors which are normal to each face of the 
        %   convex hull and point inwards, which can be used (together with the face 
        %   points) to determine if given points are inside the convex hull or not.
        %
        %   CONVEXHULLFACENORMAL : (nHullFace,nDesignVariable) double
        %
        %   See also convex_hull_plane, is_in_convex_hull_with_plane.   
        ConvexHullFaceNormal
    end
    
    methods
        function obj = CandidateSpaceConvexHull(designSpaceLowerBound,designSpaceUpperBound,varargin)
        %CANDIDATESPACECONVEXHULL Constructor
        %   CANDIDATESPACECONVEXHULL is a constructor initializes an object of
        %   this class.
        %
        %   OBJ = CANDIDATESPACECONVEXHULL(DESIGNSPACELOWERBOUND,
        %   DESIGNSPACEUPPERBOUND) creates an object of the 
        %   CandidateSpaceConvexHull class and sets its design space boundaries
        %   to its input values DESIGNSPACELOWERBOUND and DESIGNSPACEUPPERBOUND.
        %   Other properties are set to empty.
        %
        %   OBJ = CANDIDATESPACECONVEXHULL(...,NAME,VALUE,...) also allows one to set
        %   specific options for the object. This can be 
        %       - 'SamplingBoxSlack' : where the boundaries of the sampling box 
        %       will be relative to the strictest bounding box and the most relaxed
        %       bounding box. A value of 0 means no slack and therefore the sampling
        %       box will be the most strict one possible, and 1 means the sampling
        %       box will be the most relaxed.
        %
        %   Inputs:
        %       - DESIGNSPACELOWERBOUND : (1,nDesignVariable) double
        %       - DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double
        %       - 'SamplingBoxSlack' : double
        %
        %   Outputs:
        %       - OBJ : CandidateSpaceConvexHull
        %   
        %   See also convex_hull_plane.

            % parse inputs
            parser = inputParser;
            parser.addRequired('designSpaceLowerBound',@(x)isnumeric(x)&&(size(x,1)==1));
            parser.addRequired('designSpaceUpperBound',@(x)isnumeric(x)&&(size(x,1)==1));
            parser.addParameter('SamplingBoxSlack',0.5,@(x)isnumeric(x)&&isscalar(x)&&(x>=0)&&(x<=1));
            parser.parse(designSpaceLowerBound,designSpaceUpperBound,varargin{:});

            
            obj.DesignSpaceLowerBound = parser.Results.designSpaceLowerBound;;
            obj.DesignSpaceUpperBound = parser.Results.designSpaceUpperBound;
            obj.SamplingBoxSlack = parser.Results.SamplingBoxSlack;

            obj.ConvexHullIndex = [];
            obj.ConvexHullFacePoint = [];
            obj.ConvexHullFaceNormal = [];
            obj.Measure = [];
            obj.DesignSampleDefinition = [];
            obj.IsInsideDefinition = [];
        end
        
        function obj = define_candidate_space(obj,designSample,isInside)
        %DEFINE_CANDIDATE_SPACE Initial definition of the candidate space
        %   DEFINE_CANDIDATE_SPACE uses labeled design samples to define the inside / 
        %   outside regions of the candidate space. For CandidateSpaceConvexHull, this
        %   means a convex hull is created around the inside designs.
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
        %       - OBJ : CandidateSpaceConvexHull
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %       - ISINSIDE : (nSample,1) logical
        %   
        %   Outputs:
        %       - OBJ : CandidateSpaceConvexHull
        %   
        %   See also convex_hull_plane, is_in_candidate_space.

            if(nargin<3)
                isInside = true(size(designSample,1),1);
            end
            
            if(any(designSample<obj.DesignSpaceLowerBound,'all') || ...
                any(designSample>obj.DesignSpaceUpperBound,'all'))
                error('CandidateSpaceConvexHull:Definition:OutsideDesignSpace',...
                    'Specified designs are outside the design space');
            end
            
            % Convex Hull Definition
            obj.DesignSampleDefinition = designSample;
            obj.IsInsideDefinition = isInside;
            convexHullPoint = designSample(isInside,:);
            
            % compute convex hull and find reference points / normals to each facet
            [obj.ConvexHullIndex,obj.Measure] = compute_convex_hull(convexHullPoint);
            [obj.ConvexHullFacePoint,obj.ConvexHullFaceNormal] = find_facet_reference_point_normal(convexHullPoint,obj.ConvexHullIndex);

            % make sure normal vectors are pointing inside
            convexHullCenter = mean(convexHullPoint,1);
            distanceCenterToVertice = obj.ConvexHullFacePoint - convexHullCenter;
            dotProduct = dot(distanceCenterToVertice,obj.ConvexHullFaceNormal,2);
            wrongOrientation = (dotProduct>0);
            obj.ConvexHullFaceNormal(wrongOrientation,:) = -obj.ConvexHullFaceNormal(wrongOrientation,:);

            % label designs which contribute to shape definition
            nSample = size(designSample,1);
            globalConvexHullIndex = convert_index_base(isInside,obj.ConvexHullIndex(:),'backward');
            obj.IsShapeDefinition = ismember((1:nSample)',globalConvexHullIndex);
        end
        
        function obj = grow_candidate_space(obj,growthRate)
        %GROW_CANDIDATE_SPACE Expansion of candidate space by given factor
        %   GROW_CANDIDATE_SPACE will grow the region considered inside the current 
        %   candidate space by the factor given. Said growth is done in a fixed rate 
        %   defined by the input relative to the design space.
        %   This is done by finding the center of the convex hull and then making all 
        %   inside designs move opposite to that direction. 
        %
        %   OBJ = OBJ.GROW_CANDIDATE_SPACE(GROWTHRATE) will growth the candidate space 
        %   defined in OBJ by a factor of GROWTHRATE. This is an isotropic expansion of 
        %   the candidate space by a factor of the growth rate times the size of the 
        %   design space.
        %
        %   Inputs:
        %       - OBJ : CandidateSpaceConvexHull
        %       - GROWTHRATE : double
        %   
        %   Outputs:
        %       - OBJ : CandidateSpaceConvexHull
        %   
        %   See also define_candidate_space, is_in_candidate_space.

            center = mean(obj.ActiveDesign,1);
            
            distances = obj.ActiveDesign - center;
            directionGrowth = distances./vecnorm(distances,2,2);
            designSpaceFactor = obj.DesignSpaceUpperBound - obj.DesignSpaceLowerBound;
            designSpace = [obj.DesignSpaceLowerBound;obj.DesignSpaceUpperBound];

            % find maximum growth rate not to escape design space
            maxGrowthRate = region_limit_line_search([],obj.ActiveDesign,designSpaceFactor.*directionGrowth,designSpace);
            sampleGrowthRate = min(growthRate,maxGrowthRate);

            newSamples = obj.ActiveDesign + sampleGrowthRate.*designSpaceFactor.*directionGrowth;
            newSamples = max(newSamples, obj.DesignSpaceLowerBound); % lower bound limit
            newSamples = min(newSamples, obj.DesignSpaceUpperBound); % upper bound limit
            newSamples = unique([obj.ActiveDesign;newSamples],'rows');
            obj = obj.define_candidate_space(newSamples);
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
        %       - OBJ : CandidateSpaceConvexHull
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %   
        %   Outputs:
        %       - ISINSIDE : (nSample,1) logical
        %       - SCORE : (nSample,1) double
        %   
        %   See also is_in_convex_hull_with_plane.

            [label,score] = is_in_convex_hull_with_facet_normal(obj.ConvexHullFacePoint,obj.ConvexHullFaceNormal,designSample);
            isInBoundary = ismember(designSample,obj.ActiveDesign,'rows');
            label = label | isInBoundary;
            score(isInBoundary) = -abs(score(isInBoundary));
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
        %   specification for additional options in the process. For 2D candidate 
        %   spaces, these options should refer to 'plot', and for 3D spaces, they
        %   should refer to 'trisurf'.
        %
        %   Input:
        %       - OBJ : CandidateSpaceConvexHull
        %       - FIGUREHANDLE : Figure
        %
        %   Output:
        %       - PLOTHANDLE : line OR trisurf-object
        %
        %   See also plot_convex_hull_2d, plot_convex_hull_3d.

            nDimension = size(obj.ActiveDesign,2);
            if(nDimension==2)
                plotHandle = plot_convex_hull_2d(figureHandle,obj.ActiveDesign,obj.ConvexHullIndex,varargin{:});
            elseif(nDimension==3)
                plotHandle = plot_convex_hull_3d(figureHandle,obj.ActiveDesign,obj.ConvexHullIndex,varargin{:});
            else
                plotHandle = [];
            end
        end
    end
end