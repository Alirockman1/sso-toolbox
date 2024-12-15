classdef CandidateSpaceCornerBoxRemoval < CandidateSpaceBase
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
        
        %
        AnchorPoint

        % 
        CornerDirection
    end

    properties (SetAccess = protected, Dependent)
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
        %   In this case, this is the sum of areas/volumes/... of the simplices as 
        %   computed.
        %
        %   MEASURE : double
        %
        %   See also DesignSampleDefinition, IsInsideDefinition, convhull, convhulln.   
        Measure
    end
    
    methods
        function obj = CandidateSpaceCornerBoxRemoval(designSpaceLowerBound,designSpaceUpperBound,varargin)
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

            obj.DesignSampleDefinition = [];
            obj.AnchorPoint = [];
            obj.CornerDirection = [];
        end
        
        function obj = define_candidate_space(obj,designSample,trimmingInformation)
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
            obj.DesignSampleDefinition = designSample;

            if(isempty(trimmingInformation))
                return;
            end

            obj.AnchorPoint = vertcat(trimmingInformation.Anchor);
            obj.CornerDirection = vertcat(trimmingInformation.CornerDirection);
        end

        function obj = update_candidate_space(obj,designSample,isInside,trimmingInformation)
            if(isempty(obj.DesignSampleDefinition) || isempty(obj.AnchorPoint))
                obj = obj.define_candidate_space(designSample,trimmingInformation);
                return;
            end

            if(isempty(trimmingInformation))
                return;
            end

            anchorPointNew = vertcat(trimmingInformation.Anchor);
            cornerDirectionNew = vertcat(trimmingInformation.CornerDirection);

            obj.AnchorPoint = [obj.AnchorPoint;anchorPointNew];
            obj.CornerDirection = [obj.CornerDirection;cornerDirectionNew];
            
            % keep samples in inside/outside
            [~,iLowerBoundaryAll] = min(obj.DesignSampleDefinition,[],1);
            [~,iUpperBoundaryAll] = max(obj.DesignSampleDefinition,[],1);

            isInsideDefinition = obj.IsInsideDefinition;
            [~,iLowerBoundaryInside] = min(obj.DesignSampleDefinition(isInsideDefinition,:),[],1);
            [~,iUpperBoundaryInside] = max(obj.DesignSampleDefinition(isInsideDefinition,:),[],1);
            iBoundaryInside = convert_index_base(isInsideDefinition,[iLowerBoundaryInside,iUpperBoundaryInside]','backward');

            obj.DesignSampleDefinition = unique(...
                [obj.DesignSampleDefinition([iLowerBoundaryAll,iUpperBoundaryAll,iBoundaryInside'],:);...
                designSample;...
                obj.AnchorPoint],'rows');
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
            designSpaceFactor = obj.DesignSpaceUpperBound - obj.DesignSpaceLowerBound;
            designSpace = [obj.DesignSpaceLowerBound;obj.DesignSpaceUpperBound];

            center = mean(obj.DesignSampleDefinition(obj.IsInsideDefinition,:),1);
            distanceToCenter = obj.DesignSampleDefinition - center;
            directionGrowth = distanceToCenter./vecnorm(distanceToCenter,2,2);

            maxGrowthRate = region_limit_line_search([],obj.DesignSampleDefinition,designSpaceFactor.*directionGrowth,designSpace);
            sampleGrowthRate = min(growthRate,maxGrowthRate);
            designSampleNew = obj.DesignSampleDefinition + sampleGrowthRate.*designSpaceFactor.*directionGrowth;
            designSampleNew = min(max(designSampleNew,obj.DesignSpaceLowerBound),obj.DesignSpaceUpperBound);
            obj.DesignSampleDefinition = unique([obj.DesignSampleDefinition;designSampleNew],'rows');

            if(~isempty(obj.AnchorPoint))
                % connect each anchor to its respective corner
                anchorCorner = nan(size(obj.AnchorPoint,1),size(obj.DesignSpaceLowerBound,2));
                for i=1:size(obj.DesignSpaceLowerBound,2)
                    anchorCorner(~obj.CornerDirection(:,i),i) = obj.DesignSpaceLowerBound(i);
                    anchorCorner(obj.CornerDirection(:,i),i) = obj.DesignSpaceUpperBound(i);
                end
                directionGrowth = anchorCorner - obj.AnchorPoint;
                
                maxGrowthRate = region_limit_line_search([],obj.AnchorPoint,designSpaceFactor.*directionGrowth,designSpace);
                anchorGrowthRate = min(growthRate,maxGrowthRate);
                anchorPointNew = obj.AnchorPoint + anchorGrowthRate.*designSpaceFactor.*directionGrowth;
                obj.AnchorPoint = min(max(anchorPointNew,obj.DesignSpaceLowerBound),obj.DesignSpaceUpperBound);

                % update definition
                obj.DesignSampleDefinition = unique([obj.DesignSampleDefinition;obj.AnchorPoint],'rows');
            end
        end
        
        function [isInside, score] = is_in_candidate_space(obj,designSample)
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

            nSample = size(designSample,1);
            if(isempty(obj.DesignSampleDefinition))
                isInside = true(nSample,1);
                score = zeros(nSample,1);
                return;
            elseif(isempty(obj.AnchorPoint))
                boundingBox = design_bounding_box(obj.DesignSampleDefinition);
                [isInside,score] = is_in_design_box(designSample,boundingBox);
                return;
            else
                isInside = true(nSample,1);
                score = nan(nSample,1);
            end
            
            nAnchor = size(obj.AnchorPoint,1);
            for i=1:nAnchor
                designLesser = (designSample-obj.AnchorPoint(i,:)<=0);
                designGreater = (designSample-obj.AnchorPoint(i,:)>=0);

                combinationLesser = all(designLesser(:,~obj.CornerDirection(i,:)),2);
                combinationGreater = all(designGreater(:,obj.CornerDirection(i,:)),2);
                isInside(combinationLesser & combinationGreater) = false;
            end
            score(isInside) = -1;
            score(~isInside) = 1;
        end

        function isInside = get.IsInsideDefinition(obj)
            isInside = obj.is_in_candidate_space(obj.DesignSampleDefinition);
        end

        function isShapeDefinition = get.IsShapeDefinition(obj)
            if(isempty(obj.DesignSampleDefinition))
                isShapeDefinition = [];
            else
                [~,iLowerBoundaryAll] = min(obj.DesignSampleDefinition,[],1);
                [~,iUpperBoundaryAll] = max(obj.DesignSampleDefinition,[],1);

                isInsideDefinition = obj.IsInsideDefinition;
                [~,iLowerBoundaryInside] = min(obj.DesignSampleDefinition(isInsideDefinition,:),[],1);
                [~,iUpperBoundaryInside] = max(obj.DesignSampleDefinition(isInsideDefinition,:),[],1);
                iBoundaryInside = convert_index_base(isInsideDefinition,[iLowerBoundaryInside,iUpperBoundaryInside]','backward');

                isShapeDefinition = false(size(obj.DesignSampleDefinition,1),1);
                isShapeDefinition([iLowerBoundaryAll,iUpperBoundaryAll,iBoundaryInside']) = true;

                isInBoundary = design_find_boundary_samples(obj.DesignSampleDefinition,isInsideDefinition);
                isShapeDefinition = isShapeDefinition | isInBoundary;
                
                if(~isempty(obj.AnchorPoint))
                    isShapeDefinition(ismember(obj.DesignSampleDefinition,obj.AnchorPoint,'rows')) = true;
                end
            end
        end

        function volume = get.Measure(obj)
            nSample = size(obj.DesignSampleDefinition,1);
            samplingBox = obj.SamplingBox;
            
            volumeSample = sampling_latin_hypercube(samplingBox,nSample);
            isInside = obj.is_in_candidate_space(volumeSample);
            volumeFactor = sum(isInside) / size(isInside,1);
            volume = volumeFactor * prod(samplingBox(2,:) - samplingBox(1,:));
        end
    end
end