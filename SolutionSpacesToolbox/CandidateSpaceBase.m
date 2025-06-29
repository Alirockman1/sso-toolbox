classdef (Abstract) CandidateSpaceBase
%CANDIDATESPACEBASE Candidate Space abstract class for component solution spaces
%   CANDIDATESPACEBASE provides the necessary interfaces for particular 
%   implementations that may serve as candidate spaces when performing the
%   computation of component solution spaces.
%
%   As an abstract class, particular implementations have to be created using 
%   this as base for one to be able to create objects.
%
%   CANDIDATESPACEBASE properties:
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
%
%   CANDIDATESPACEBASE methods:
%       - define_candidate_space : create a candidate space based on design 
%       samples that are labeled as inside/outside.
%       - grow_candidate_space : expand the candidate space by a given factor.
%       - is_in_candidate_space : verify if given design samples are inside 
%       the candidate space.
%       - plot_candidate_space : visualize 1D/2D/3D candidate spaces in given
%       figure.
%
%   See also CandidateSpaceConvexHull, CandidateSpaceSvm.
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

    properties (SetAccess = protected, Abstract)
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
        %   contributes to the shape of the candidate space. For example, for a
        %   convex hull, this would be the convex hull index points.
        %
        %   ISSHAPEDEFINITION : (nSample,1) logical
        %
        %   See also DesignSampleDefinition, IsInsideDefinition.
        IsShapeDefinition
        
        %MEASURE Size measure of the candidate space
        %   MEASURE is a value that works as the measure of the candidate space.
        %   This may be its volume, or normalized volume relative to its design
        %   space, or some other metric.
        %
        %   MEASURE : double
        %
        %   See also DesignSampleDefinition, IsInsideDefinition.   
        Measure
    end

    properties (Dependent)
        %ACTIVEDESIGN All sample points that are labeled as inside candidate space
        %   ACTIVEDESIGN are all the sample points labeled as inside the candidate 
        %   space.
        %
        %   ACTIVEDESIGN : (nInside,nDesignVariable) double
        %
        %   See also DesignSampleDefinition, IsInsideDefinition.
        ActiveDesign

        %SAMPLINGBOX Bounding box of inside region used to help with sampling
        %   SAMPLINGBOX is a bounding box formed around the internal region of the
        %   candidate space. It can be used to facilitate trying to sample inside said
        %   space.
        %
        %   SAMPLINGBOX : (2,nDesignVariable) double
        %       - (1) : lower boundary of the design box
        %       - (2) : upper boundary of the design box
        %
        %   See also SamplingBoxSlack.
        SamplingBox
    end

    properties (SetAccess = protected)
        %DESIGNSPACELOWERBOUND Lower boundary of the design space
        %   DESIGNSPACELOWERBOUND is the lower boundary of the design space. Design 
        %   samples which are given outside the design space are not considered in its 
        %   definition / are considered outside of its region.
        %
        %   DESIGNSPACELOWERBOUND : (1,nDesignVariable) double
        %
        %   See also DesignSpaceUpperBound.
        DesignSpaceLowerBound

        %DESIGNSPACEUPPERBOUND Upper boundary of the design space
        %   DESIGNSPACEUPPERBOUND is the upper boundary of the design space. Design 
        %   samples which are given outside the design space are not considered in its 
        %   definition / are considered outside region.
        %
        %   DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double
        %
        %   See also DesignSpaceLowerBound.
        DesignSpaceUpperBound

        %SAMPLINGBOXSLACK Slack allowed for for the sampling box
        %   SAMPLINGBOXSLACK defines where the boundaries of the sampling box will be 
        %   relative to the strictest bounding box and the most relaxed bounding box. A 
        %   value of 0 means no slack and therefore the sampling box will be the most 
        %   strict one possible, and 1 means the sampling box will be the most relaxed.
        %
        %   SAMPLINGBOXSLACK : double
        %
        %   See also SamplingBox, design_bounding_box.
        SamplingBoxSlack
    end

    methods (Abstract)
        %DEFINE_CANDIDATE_SPACE Initial definition of the candidate space
        %   DEFINE_CANDIDATE_SPACE uses labeled design samples to define the inside / 
        %   outside regions of the candidate space.
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
        %       - OBJ : CandidateSpaceBase
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %       - ISINSIDE : (nSample,1) logical
        %   
        %   Outputs:
        %       - OBJ : CandidateSpaceBase
        %   
        %   See also is_in_candidate_space.
        obj = define_candidate_space(obj,designSample,isInside)

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
        %       - OBJ : CandidateSpaceBase
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %   
        %   Outputs:
        %       - ISINSIDE : (nSample,1) logical
        %       - SCORE : (nSample,1) double
        %   
        %   See also is_in_design_box, is_in_convex_hull_with_plane.
        [isInside,score] = is_in_candidate_space(obj,designSample)


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
        obj = grow_candidate_space(obj,growthRate)
    end

    methods
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
        %   spaces, these options should refer to 'contour', and for 3D spaces, they
        %   should refer to 'patch'.
        %
        %   Input:
        %       - OBJ : CandidateSpaceBase
        %       - FIGUREHANDLE : Figure
        %
        %   Output:
        %       - PLOTHANDLE : contour object OR patch object
        %
        %   See also plot_candidate_space_2d, plot_candidate_space_3d.

            nDimension = size(obj.DesignSampleDefinition,2);
            if(nDimension==2)
                plotHandle = plot_candidate_space_2d(figureHandle,obj,varargin{:});
            elseif(nDimension==3)
                plotHandle = plot_candidate_space_3d(figureHandle,obj,varargin{:});
            else
                plotHandle = [];
            end
        end

        function activeDesign = get.ActiveDesign(obj)
            activeDesign = obj.DesignSampleDefinition(obj.IsInsideDefinition,:);
        end

        function samplingBox = get.SamplingBox(obj)
            [boundingBoxStrict,boundingBoxRelaxed] = design_bounding_box(...
                obj.DesignSampleDefinition,obj.IsInsideDefinition);
            samplingBox = (1-obj.SamplingBoxSlack).*boundingBoxStrict + ...
                obj.SamplingBoxSlack.*boundingBoxRelaxed;
        end
    end
end