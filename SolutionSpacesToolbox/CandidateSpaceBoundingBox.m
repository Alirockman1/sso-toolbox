classdef CandidateSpaceBoundingBox < CandidateSpaceBase
%CANDIDATESPACEBOUNDINGBOX Candidate Space defined by a bounding box
%   CANDIDATESPACEBOUNDINGBOX can be used to define candidate spaces for  
%   component solution space computation. In this case, the candidate space is 
%   defined by intervals for each design variable that create a bounding 
%   box that encompasses all "inside" designs.
%   Note that if the trimming operation does not result in box-shaped solutions,
%   this candidate space definition may also include "outside" designs inside
%   it.
%
%   This class is derived from CandidateSpaceBase.
%
%   CANDIDATESPACEBOUNDINGBOX properties:
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
%
%   CANDIDATESPACEBOUNDINGBOX methods:
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

    properties (SetAccess = protected, Dependent)    
        %MEASURE Size measure of the candidate space
        %   MEASURE is a value that works as the measure of the candidate space. In this
        %   case, the area/volume/hypervolume of the bounding box is used.
        %
        %   MEASURE : double
        %
        %   See also DesignSampleDefinition, IsInsideDefinition.  
        Measure
    end

    properties (SetAccess = protected)
        %DESIGNSAMPLEDEFINITION Design sample points used in candidate space definition
        %   DESIGNSAMPLEDEFINITION are the sample points used in the definition of the 
        %   current candidate space.
        %
        %   DESIGNSAMPLEDEFINITION : (nSample,nDesignVariable) double
        %
        %   See also IsInsideDefinition, IsIsShapeDefinition.
        DesignSampleDefinition
        
        %ISINSIDEDEFINITION Labels of sample points used in candidate space definition
        %   ISINSIDEDEFINITION are the labels of design samples used in the definition 
        %   of the current candidate space. A label of 'true' indicates the respective 
        %   design point is inside the candidate space.
        %
        %   ISINSIDEDEFINITION : (nSample,1) logical
        %
        %   See also DesignSampleDefinition, IsIsShapeDefinition.
        IsInsideDefinition

        %ISSHAPEDEFINITION Labels if sample points from definition contributes to shape
        %   ISSHAPEDEFINITION is a logical array where 'true' values indicate that that
        %   design point (in the respective row) from the definition sample actively 
        %   contributes to the shape of the candidate space. 
        %   In this case, these are the closest inside/outside design sample points for
        %   each boundary.
        %
        %   ISSHAPEDEFINITION : (nSample,1) logical
        %
        %   See also DesignSampleDefinition, IsInsideDefinition.
        IsShapeDefinition

        %BOUNDINGBOX Bounding box which defines the candidate space
        %   BOUNDINGBOX is the bounding box which is created when the candidate
        %   space is defined. It contains all the designs labeled as inside.
        %
        %   BOUNDINGBOX : (2,nDesignVariable) double
        %       - (1) : lower boundary of the design box
        %       - (2) : upper boundary of the design box
        %
        %   See also DesignSampleDefinition, IsInsideDefinition.
        BoundingBox
    end
    
    methods
        function obj = CandidateSpaceBoundingBox(designSpaceLowerBound,designSpaceUpperBound,varargin)
        %CANDIDATESPACEBOUNDINGBOX Constructor
        %   CANDIDATESPACEBOUNDINGBOX is a constructor initializes an object of
        %   this class.
        %
        %   OBJ = CANDIDATESPACEBOUNDINGBOX(DESIGNSPACELOWERBOUND,
        %   DESIGNSPACEUPPERBOUND) creates an object of the CandidateSpaceBoundingBox
        %   class and sets its design space boundaries to its input values 
        %   DESIGNSPACELOWERBOUND and DESIGNSPACEUPPERBOUND. Other properties are set to
        %   empty.
        %
        %   OBJ = CANDIDATESPACEBOUNDINGBOX(...,NAME,VALUE,...) allows for the 
        %   specification of additional optional values. These are:
        %       - 'SamplingBoxSlack' : slack used when computing the bounding box used
        %       for sampling. Default value is 0.
        %
        %   Inputs:
        %       - DESIGNSPACELOWERBOUND : (1,nDesignVariable) double
        %       - DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double
        %       - 'SamplingBoxSlack' : double
        %
        %   Outputs:
        %       - OBJ : CandidateSpaceBoundingBox object
        %   
        %   See also .

            parser = inputParser;
            parser.addParameter('SamplingBoxSlack',0);
            parser.parse(varargin{:});

            obj.BoundingBox = [];
            obj.DesignSampleDefinition = [];
            obj.IsInsideDefinition = [];
            obj.IsShapeDefinition = [];
            obj.SamplingBoxSlack = parser.Results.SamplingBoxSlack;
            obj.DesignSpaceLowerBound = designSpaceLowerBound;
            obj.DesignSpaceUpperBound = designSpaceUpperBound;
        end
        
        function obj = define_candidate_space(obj,designSample,isInside)
        %DEFINE_CANDIDATE_SPACE Initial definition of the candidate space
        %   DEFINE_CANDIDATE_SPACE uses labeled design samples to define the inside / 
        %   outside regions of the candidate space. For CandidateSpaceBoundingBox, this 
        %   means a bounding box is created around the inside designs, and everything 
        %   inside that box is then labeled as being inside the space as well.
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
        %       - OBJ : CandidateSpaceBoundingBox
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %       - ISINSIDE : (nSample,1) logical
        %   
        %   Outputs:
        %       - OBJ : CandidateSpaceBoundingBox
        %   
        %   See also is_in_candidate_space.

            if(nargin<3)
                isInside = true(size(designSample,1),1);
            end

            obj.DesignSampleDefinition = designSample;
            obj.IsInsideDefinition = isInside;

            obj.BoundingBox = design_bounding_box(obj.DesignSampleDefinition,obj.IsInsideDefinition);

            % find designs that produce the result
            nDimension = size(obj.BoundingBox,2);
            isShapeDefinitionIndex = [];
            for i=1:nDimension
                otherDimension = true(1,nDimension);
                otherDimension(i) = false;

                % find which designs are in the expanded region of the current edges
                insideBoxOtherDimension = is_in_design_box(designSample(:,otherDimension),obj.BoundingBox(:,otherDimension));

                % find designs outside/inside closest to the boundary / edge
                distanceToLowerBoundary = designSample(insideBoxOtherDimension,i) - obj.BoundingBox(1,i);
                distanceToUpperBoundary = obj.BoundingBox(2,i) - designSample(insideBoxOtherDimension,i);
                [~,orderLower] = sort(abs(distanceToLowerBoundary));
                [~,orderUpper] = sort(abs(distanceToUpperBoundary));

                % closest to lower boundary/edge, outside box
                boundaryDesignLowerOutside = find(distanceToLowerBoundary(orderLower)<0,1,'first');
                % closest to lower boundary/edge, inside box
                boundaryDesignLowerInside = find(distanceToLowerBoundary(orderLower)>=0,1,'first');
                % closest to upper boundary/edge, outside box
                boundaryDesignUpperOutside = find(distanceToUpperBoundary(orderUpper)<0,1,'first');
                % closest to upper boundary/edge, inside box
                boundaryDesignUpperInside = find(distanceToUpperBoundary(orderUpper)>=0,1,'first');

                % add those to shape definition index
                isShapeDefinitionIndex = [...
                    isShapeDefinitionIndex;...
                    convert_index_base(insideBoxOtherDimension,orderLower(boundaryDesignLowerOutside),'backward');...
                    convert_index_base(insideBoxOtherDimension,orderLower(boundaryDesignLowerInside),'backward');...
                    convert_index_base(insideBoxOtherDimension,orderUpper(boundaryDesignUpperOutside),'backward');...
                    convert_index_base(insideBoxOtherDimension,orderUpper(boundaryDesignUpperInside),'backward')];
            end

            nSample = size(designSample,1);
            obj.IsShapeDefinition = ismember((1:nSample)',isShapeDefinitionIndex);
        end
        
        function obj = grow_candidate_space(obj,growthRate)
        %GROW_CANDIDATE_SPACE Expansion of candidate space by given factor
        %   GROW_CANDIDATE_SPACE will grow the region considered inside the current 
        %   candidate space by the factor given. Said growth is done in a fixed rate 
        %   defined by the input relative to the design space.
        %   This is done by extending the bounding box equally in all dimensions.
        %
        %   OBJ = OBJ.GROW_CANDIDATE_SPACE(GROWTHRATE) will growth the candidate space 
        %   defined in OBJ by a factor of GROWTHRATE. This is an isotropic expansion of 
        %   the candidate space by a factor of the growth rate times the size of the 
        %   design space.
        %
        %   Inputs:
        %       - OBJ : CandidateSpaceBoundingBox
        %       - GROWTHRATE : double
        %   
        %   Outputs:
        %       - OBJ : CandidateSpaceBoundingBox
        %   
        %   See also define_candidate_space, design_box_grow_fixed.

            % Expand candidate solution box in both sides of each interval isotroply
            obj.BoundingBox = design_box_grow_fixed(obj.BoundingBox,...
                obj.DesignSpaceLowerBound,...
                obj.DesignSpaceUpperBound,...
                growthRate);
        end
        
        function [isInside,score] = is_in_candidate_space(obj,designSample)
        %IS_IN_CANDIDATE_SPACE Verification if given design samples are inside space
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
        %   See also is_in_design_box.

            [isInside,score] = is_in_design_box(designSample,obj.BoundingBox);
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
        %   specification for additional options in the process. 
        %   These options refer to 'patch'.
        %
        %   Input:
        %       - OBJ : CandidateSpaceBase
        %       - FIGUREHANDLE : Figure
        %
        %   Output:
        %       - PLOTHANDLE : patch object
        %
        %   See also plot_design_box_2d, plot_design_box_3d, patch.

            nDimension = size(obj.DesignSampleDefinition,2);
            if(nDimension==2)
                plotHandle = plot_design_box_2d(figureHandle,obj.BoundingBox,varargin{:});
            elseif(nDimension==3)
                plotHandle = plot_design_box_3d(figureHandle,obj.BoundingBox,varargin{:});
            else
                plotHandle = [];
            end
        end

        function isShapeDefinition = get.IsShapeDefinition(obj)
            isShapeDefinition = obj.IsShapeDefinition;
        end
        
        function measure = get.Measure(obj)
            measure = prod(obj.BoundingBox(2,:) - obj.BoundingBox(1,:));
        end
    end
end