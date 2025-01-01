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

        %ISINSIDEDEFINITION Labels of sample points used in candidate space definition
        %   ISINSIDEDEFINITION are the labels of design samples used in the definition 
        %   of the current candidate space. A label of 'true' indicates the respective 
        %   design point is inside the candidate space.
        %
        %   ISINSIDEDEFINITION : (nSample,1) logical
        %
        %   See also DesignSampleDefinition, IsShapeDefinition.
        IsInsideDefinition
        
        %
        AnchorPoint

        % 
        CornerDirection

        %
        DetachTolerance
    end

    properties (SetAccess = protected, Dependent)
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
            parser.addRequired('DesignSpaceLowerBound',@(x)isnumeric(x)&&(size(x,1)==1));
            parser.addRequired('DesignSpaceUpperBound',@(x)isnumeric(x)&&(size(x,1)==1));
            parser.addParameter('DetachTolerance',0.00);
            parser.parse(designSpaceLowerBound,designSpaceUpperBound,varargin{:});

            
            obj.DesignSpaceLowerBound = parser.Results.DesignSpaceLowerBound;
            obj.DesignSpaceUpperBound = parser.Results.DesignSpaceUpperBound;
            obj.DetachTolerance = parser.Results.DetachTolerance;

            obj.DesignSampleDefinition = [];
            obj.IsInsideDefinition = [];
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

            if(~isempty(trimmingInformation))
                obj.AnchorPoint = vertcat(trimmingInformation.Anchor);
                obj.CornerDirection = vertcat(trimmingInformation.CornerDirection);
            end
            obj.IsInsideDefinition = obj.is_in_candidate_space(designSample,false);
        end

        function obj = update_candidate_space(obj,designSample,isInside,labelViable,trimmingInformation)
            if(isempty(obj.DesignSampleDefinition) || isempty(obj.AnchorPoint))
                obj = obj.define_candidate_space(designSample,trimmingInformation);
                return;
            end

            if(isempty(trimmingInformation))
                return;
            end
            anchorPointNew = vertcat(trimmingInformation.Anchor);
            cornerDirectionNew = vertcat(trimmingInformation.CornerDirection);
            
            isNewAnchor = [false(size(obj.AnchorPoint,1),1);true(size(anchorPointNew,1),1)];
            obj.AnchorPoint = [obj.AnchorPoint;anchorPointNew];
            obj.CornerDirection = [obj.CornerDirection;cornerDirectionNew];

            % verify if there are any redundant anchor points
            % -> region removed includes a different anchor with same corner direction
            nAnchor = size(obj.AnchorPoint,1);
            isRedundantAnchor = false(nAnchor,1);
            for i = 1:nAnchor
                currentAnchor = obj.AnchorPoint(i,:);
                currentCornerDirection = obj.CornerDirection(i,:);

                isDesignLesser = (obj.AnchorPoint - currentAnchor<0);
                isDesignGreater = (obj.AnchorPoint - currentAnchor>0);
                combinationLesser = all(isDesignLesser(:,~currentCornerDirection),2);
                combinationGreater = all(isDesignGreater(:,currentCornerDirection),2);
                isSameDirection = all(obj.CornerDirection==currentCornerDirection,2);
                isRedundantAnchor(combinationLesser & combinationGreater & isSameDirection) = true;
            end
            obj.AnchorPoint(isRedundantAnchor,:) = [];
            obj.CornerDirection(isRedundantAnchor,:) = [];
            isNewAnchor(isRedundantAnchor) = [];

            % check for detachments
            if(obj.DetachTolerance>0)
                nAnchor = size(obj.AnchorPoint,1);
                
                lowerBound = min(obj.DesignSampleDefinition(obj.IsInsideDefinition,:),[],1);
                upperBound = max(obj.DesignSampleDefinition(obj.IsInsideDefinition,:),[],1);
                allowedSlack = obj.DetachTolerance.*(upperBound - lowerBound);

                for i=1:nAnchor
                    distanceToAnchor = obj.AnchorPoint - obj.AnchorPoint(i,:);
                    distanceToAnchor(:,obj.CornerDirection(i,:)) = -distanceToAnchor(:,obj.CornerDirection(i,:));

                    [maximumSlack,iDimension] = max(distanceToAnchor,[],2);
                    
                    shouldCollapse = (abs(maximumSlack)<=allowedSlack(:,iDimension)');
                    dimensionsToCollapse = (distanceToAnchor>=0) & (distanceToAnchor<=allowedSlack) & (obj.CornerDirection~=obj.CornerDirection(i,:));
                    currentAnchor = repmat(obj.AnchorPoint(i,:),sum(shouldCollapse),1);
                    obj.AnchorPoint(shouldCollapse & dimensionsToCollapse) = currentAnchor(dimensionsToCollapse(shouldCollapse,:));
                end
            end

            % keep samples in inside/outside
            [~,iLowerBoundaryAll] = min(obj.DesignSampleDefinition,[],1);
            [~,iUpperBoundaryAll] = max(obj.DesignSampleDefinition,[],1);

            insideSample = obj.DesignSampleDefinition(obj.IsInsideDefinition,:);
            [~,iLowerBoundaryInside] = min(insideSample,[],1);
            [~,iUpperBoundaryInside] = max(insideSample,[],1);
            iBoundaryInside = convert_index_base(obj.IsInsideDefinition,[iLowerBoundaryInside,iUpperBoundaryInside]','backward');

            obj.DesignSampleDefinition = unique(...
                [obj.DesignSampleDefinition([iLowerBoundaryAll,iUpperBoundaryAll,iBoundaryInside'],:);...
                designSample;...
                obj.AnchorPoint],'rows');
            obj.IsInsideDefinition = obj.is_in_candidate_space(obj.DesignSampleDefinition,false);
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

            %maxGrowthRate = region_limit_line_search([],obj.DesignSampleDefinition,designSpaceFactor.*directionGrowth,designSpace);
            %sampleGrowthRate = min(growthRate,maxGrowthRate);
            sampleGrowthRate = growthRate;
            designSampleNew = obj.DesignSampleDefinition + sampleGrowthRate.*designSpaceFactor.*directionGrowth;
            designSampleNew = min(max(designSampleNew,obj.DesignSpaceLowerBound),obj.DesignSpaceUpperBound);
            obj.DesignSampleDefinition = unique([obj.DesignSampleDefinition;designSampleNew],'rows');

            if(~isempty(obj.AnchorPoint))
                % connect each anchor to its respective corner
                %anchorCorner = nan(size(obj.AnchorPoint,1),size(obj.DesignSpaceLowerBound,2));
                %for i=1:size(obj.DesignSpaceLowerBound,2)
                %    anchorCorner(~obj.CornerDirection(:,i),i) = obj.DesignSpaceLowerBound(i);
                %    anchorCorner(obj.CornerDirection(:,i),i) = obj.DesignSpaceUpperBound(i);
                %end
                %distanceAnchorToCorner = anchorCorner - obj.AnchorPoint;
                %directionGrowthCorner = distanceAnchorToCorner./vecnorm(distanceAnchorToCorner,2,2);

                % grow away from center
                %distanceCenterToAnchor = obj.AnchorPoint - center;
                %directionGrowthCenter = distanceCenterToAnchor./vecnorm(distanceCenterToAnchor,2,2);
                
                %directionGrowthCorner = directionGrowthCenter;
                %directionGrowthCenter = directionGrowthCorner;
                %directionGrowth = (directionGrowthCorner + directionGrowthCenter)/2;

                %directionGrowth = nan(size(obj.AnchorPoint,1),size(obj.DesignSpaceLowerBound,2));
                %for i=1:size(obj.DesignSpaceLowerBound,2)
                %    directionGrowth(~obj.CornerDirection(:,i),i) = -1;
                %    directionGrowth(obj.CornerDirection(:,i),i) = 1;
                %end

                isInsideDefinition = obj.is_in_candidate_space(obj.DesignSampleDefinition,false);
                insideSample = obj.DesignSampleDefinition(isInsideDefinition,:);
                nAnchor = size(obj.AnchorPoint,1);
                nDimension = size(obj.DesignSpaceLowerBound,2);
                directionGrowth = zeros(nAnchor,nDimension);
                for i=1:nAnchor
                    distanceToAnchor = insideSample - obj.AnchorPoint(i,:);
                    distanceToAnchor(:,obj.CornerDirection(i,:)) = -distanceToAnchor(:,obj.CornerDirection(i,:));
                    [~,iDimension] = max(distanceToAnchor,[],2);
                    for j=1:nDimension
                        directionGrowth(i,j) = sum(iDimension==j);
                    end
                end
                %[~,preferentialDirection] = max(directionGrowth,[],2);
                %directionGrowth = zeros(size(directionGrowth));
                %directionGrowth(:,preferentialDirection) = 1;
                directionGrowth(~obj.CornerDirection) = -directionGrowth(~obj.CornerDirection);

                directionGrowth = directionGrowth./vecnorm(directionGrowth,2,2);
                anchorPointNew = obj.AnchorPoint + growthRate.*designSpaceFactor.*directionGrowth;
                anchorPointNew = min(max(anchorPointNew,obj.DesignSpaceLowerBound),obj.DesignSpaceUpperBound);
                
                % don't include anchors that were moved to the boundary corners
                isAnchorInLowerBoundary = (anchorPointNew==obj.DesignSpaceLowerBound);
                isAnchorInUpperBoundary = (anchorPointNew==obj.DesignSpaceUpperBound);
                isAnchorInCorner = all(isAnchorInLowerBoundary|isAnchorInUpperBoundary,2);

                % only move anchors that were active before
                obj.AnchorPoint = anchorPointNew(~isAnchorInCorner,:);
                obj.CornerDirection = obj.CornerDirection(~isAnchorInCorner,:);

                % check for detachments
                if(obj.DetachTolerance>0)
                    nAnchor = size(obj.AnchorPoint,1);
                    
                    lowerBound = min(obj.DesignSampleDefinition(obj.IsInsideDefinition,:),[],1);
                    upperBound = max(obj.DesignSampleDefinition(obj.IsInsideDefinition,:),[],1);
                    allowedSlack = obj.DetachTolerance.*(upperBound - lowerBound);

                    for i=1:nAnchor
                        distanceToAnchor = obj.AnchorPoint - obj.AnchorPoint(i,:);
                        distanceToAnchor(:,obj.CornerDirection(i,:)) = -distanceToAnchor(:,obj.CornerDirection(i,:));

                        [maximumSlack,iDimension] = max(distanceToAnchor,[],2);
                        
                        shouldCollapse = (abs(maximumSlack)<=allowedSlack(:,iDimension)');
                        dimensionsToCollapse = (distanceToAnchor>=0) & (distanceToAnchor<=allowedSlack) & (obj.CornerDirection~=obj.CornerDirection(i,:));
                        currentAnchor = repmat(obj.AnchorPoint(i,:),sum(shouldCollapse),1);
                        obj.AnchorPoint(shouldCollapse & dimensionsToCollapse) = currentAnchor(dimensionsToCollapse(shouldCollapse,:));
                    end
                end

                % update definition
                obj.DesignSampleDefinition = unique([obj.DesignSampleDefinition;obj.AnchorPoint],'rows');
            end
            obj.IsInsideDefinition = obj.is_in_candidate_space(obj.DesignSampleDefinition,false);
        end
        
        function [isInside, score] = is_in_candidate_space(obj,designSample,includeBoundingBox)
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
            nDimension = size(designSample,2);
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
            
            % check if it is in any of the corners
            nAnchor = size(obj.AnchorPoint,1);
            if(nAnchor<=nSample)
                for i=1:nAnchor
                    distanceToAnchor = designSample-obj.AnchorPoint(i,:);
                    distanceToAnchor(:,obj.CornerDirection(i,:)) = -distanceToAnchor(:,obj.CornerDirection(i,:));
                    
                    isInside(all(distanceToAnchor<0,2)) = false;
                    
                    anchorScore = -max(distanceToAnchor,[],2);
                    score = max(score,anchorScore);
                end
            else
                for i=1:nSample
                    distanceToAnchor = designSample(i,:) - obj.AnchorPoint;
                    distanceToAnchor(obj.CornerDirection) = -distanceToAnchor(obj.CornerDirection);
                    
                    isInside(i) = ~any(all(distanceToAnchor<0,2),1);
                    score(i) = max(-max(distanceToAnchor,[],2),[],1);
                end
            end

            % check if it's inside the bounding box
            if(nargin<3 || includeBoundingBox)
                boundingBox = design_bounding_box(obj.DesignSampleDefinition,obj.IsInsideDefinition);
                [isInsideBounding,scoreBounding] = is_in_design_box(designSample,boundingBox);
                isInside(~isInsideBounding) = false;
                score = max(score,scoreBounding);
            end
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
                
                if(~isempty(obj.AnchorPoint))
                    isShapeDefinition(ismember(obj.DesignSampleDefinition,obj.AnchorPoint,'rows')) = true;
                end
            end
        end

        function samplingBox = get.SamplingBox(obj)
            samplingBox = design_bounding_box(...
                obj.DesignSampleDefinition,obj.IsInsideDefinition);
        end

        function volume = get.Measure(obj)
            nSample = size(obj.DesignSampleDefinition,1);
            samplingBox = obj.SamplingBox;
            
            volumeSample = sampling_random(samplingBox,nSample);
            isInside = obj.is_in_candidate_space(volumeSample);
            volumeFactor = sum(isInside) / size(isInside,1);
            volume = volumeFactor * prod(samplingBox(2,:) - samplingBox(1,:));
        end
    end
end