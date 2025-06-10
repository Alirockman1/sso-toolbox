classdef CandidateSpaceCornerBoxRemoval < CandidateSpaceBase
%CANDIDATESPACECORNERBOXREMOVAL Candidate Space using corner-box removal
%   CANDIDATESPACECORNERBOXREMOVAL defines a candidate space by removing
%   certain "corner boxes" from the feasible region. These corner boxes are
%   identified via anchor points and a corner-direction logic, and designs
%   that lie inside these boxes are excluded from the candidate space. 
%
%   CANDIDATESPACECORNERBOXREMOVAL is derived from CandidateSpaceBase.
%
%   CANDIDATESPACECORNERBOXREMOVAL properties:
%       - DesignSampleDefinition : sample points used to define the candidate 
%       space.
%       - IsInsideDefinition : logical flags specifying whether each sample 
%       point is inside or outside the space (true = inside).
%       - AnchorPoint : anchor points that determine the corners to remove.
%       - CornerDirection : directions specifying how each anchor trims the 
%       space.
%       - DetachTolerance : tolerance for merging or discarding anchors.
%       - IsShapeDefinition : logical flags indicating which points directly 
%       shape the boundary of the candidate space.
%       - Measure : approximate measure (area/volume) of the candidate space.
%       - SamplingBox : bounding box around the current inside region for 
%       sampling.
%
%   CANDIDATESPACECORNERBOXREMOVAL methods:
%       - generate_candidate_space : initialize the candidate space using sample 
%       points and trimming information (anchors/corners).
%       - update_candidate_space : update the definition with new samples, 
%       possibly adding or removing corner anchors.
%       - expand_candidate_space : grow the space by a factor within the design
%       space.
%       - is_in_candidate_space : determine whether new points are inside the 
%       space.
%
%   See also: CandidateSpaceBase.
%
%   Copyright 2025 Eduardo Rodrigues Della Noce
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
        
        %ANCHORPOINT Anchor points used to define corner regions to remove
        %   ANCHORPOINT represents specific reference points that, combined with 
        %   CornerDirection, identify which corner boxes (i.e., rectangular subregions)
        %   in the design space should be removed.
        %
        %   ANCHORPOINT : (nAnchor, nDesignVariable) double
        %
        %   See also CornerDirection, generate_candidate_space.
        AnchorPoint

        %CORNERDIRECTION Corner-direction flags specifying the removal region
        %   CORNERDIRECTION is an array of booleans indicating for each dimension 
        %   whether the corner is "upper" or "lower." If true in dimension j, the 
        %   removal corner extends toward the upper boundary in that dimension, else 
        %   it extends toward the lower boundary.
        %
        %   CORNERDIRECTION : (nAnchor, nDesignVariable) logical
        %
        %   See also AnchorPoint.
        CornerDirection

        %DETACHTOLERANCE Tolerance for merging or discarding nearby anchors
        %   DETACHTOLERANCE defines a threshold for how close anchors must be 
        %   before they are treated as the same anchor point. A value of zero 
        %   means no merging. Larger values cause anchors that lie within a 
        %   fraction of the bounding box range to be collapsed or removed.
        %
        %   DETACHTOLERANCE : double
        %
        %   See also update_candidate_space.
        DetachTolerance

        %NORMALIZEGROWTHDIRECTION Determine if growth direction should be normalized
        %   When true, the growth direction is normalized to the design space.
        %
        %   NORMALIZEGROWTHDIRECTION : logical
        NormalizeGrowthDirection

        %CHECKREDUNDANTRIMMINGGROWTH Determine if redundant trimming is checked
        %   When true, the redundant trimming is checked.
        %
        %   CHECKREDUNDANTRIMMINGGROWTH : logical
        CheckRedundantTrimmingGrowth

        %CHECKREDUNDANTRIMMINGUPDATE Determine if redundant trimming is checked
        %   When true, the redundant trimming is checked.
        %
        %   CHECKREDUNDANTRIMMINGUPDATE : logical
        CheckRedundantTrimmingUpdate

        %CHECKDUPLICATEPOINTSGROWTH Determine if duplicate points are checked
        %   When true, the duplicate points are checked.
        %
        %   CHECKDUPLICATEPOINTSGROWTH : logical
        CheckDuplicatePointsGrowth

        %CHECKDUPLICATEPOINTSUPDATE Determine if duplicate points are checked
        %   When true, the duplicate points are checked.
        %
        %   CHECKDUPLICATEPOINTSUPDATE : logical
        CheckDuplicatePointsUpdate

        %MEASUREESTIMATIONFACTOR Factor to estimate the measure of the candidate space
        %   MEASUREESTIMATIONFACTOR is a factor that is used to estimate the measure of the 
        %   candidate space. This is used to determine the number of samples to use when 
        %   generating the candidate space.
        %
        %   MEASUREESTIMATIONFACTOR : double
        MeasureEstimationFactor
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
            parser.addParameter('DetachTolerance',0);
            parser.addParameter('NormalizeGrowthDirection',false,@islogical);
            parser.addParameter('CheckRedundantTrimmingGrowth',true,@islogical);
            parser.addParameter('CheckRedundantTrimmingUpdate',true,@islogical);
            parser.addParameter('CheckDuplicatePointsGrowth',true,@islogical);
            parser.addParameter('CheckDuplicatePointsUpdate',true,@islogical);
            parser.addParameter('MeasureEstimationFactor',10);
            parser.parse(designSpaceLowerBound,designSpaceUpperBound,varargin{:});

            obj.DesignSpaceLowerBound = parser.Results.DesignSpaceLowerBound;
            obj.DesignSpaceUpperBound = parser.Results.DesignSpaceUpperBound;
            obj.DetachTolerance = parser.Results.DetachTolerance;
            obj.NormalizeGrowthDirection = parser.Results.NormalizeGrowthDirection;
            obj.CheckRedundantTrimmingGrowth = parser.Results.CheckRedundantTrimmingGrowth;
            obj.CheckRedundantTrimmingUpdate = parser.Results.CheckRedundantTrimmingUpdate;
            obj.CheckDuplicatePointsGrowth = parser.Results.CheckDuplicatePointsGrowth;
            obj.CheckDuplicatePointsUpdate = parser.Results.CheckDuplicatePointsUpdate;
            obj.MeasureEstimationFactor = parser.Results.MeasureEstimationFactor;

            obj.DesignSampleDefinition = [];
            obj.IsInsideDefinition = [];
            obj.AnchorPoint = [];
            obj.CornerDirection = [];
        end
        
        function obj = generate_candidate_space(obj,designSample,trimmingInformation)
        %GENERATE_CANDIDATE_SPACE Define the candidate space from sample points
        %   GENERATE_CANDIDATE_SPACE uses sample points and trimming information
        %   to determine which corners should be removed. The inside definition 
        %   is then updated accordingly.
        %
        %   OBJ = OBJ.GENERATE_CANDIDATE_SPACE(DESIGNSAMPLE,TRIMMINGINFORMATION)
        %   sets DesignSampleDefinition to DESIGNSAMPLE, and if 
        %   trimmingInformation is provided, populates AnchorPoint and
        %   CornerDirection. 
        %
        %   Inputs:
        %       - OBJ : CandidateSpaceCornerBoxRemoval
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %       - TRIMMINGINFORMATION : (nCorner,1) struct
        %           -- Anchor : (1,nDimension) double
        %           -- CornerDirection : (1,nDimension) logical
        %
        %   Outputs:
        %       - OBJ : CandidateSpaceCornerBoxRemoval
        %
        %   See also AnchorPoint, CornerDirection, is_in_candidate_space.

            obj.DesignSampleDefinition = designSample;

            if(~isempty(trimmingInformation))
                obj.AnchorPoint = vertcat(trimmingInformation.Anchor);
                obj.CornerDirection = vertcat(trimmingInformation.CornerDirection);
            end
            obj.IsInsideDefinition = obj.is_in_candidate_space(designSample,false);
        end

        function obj = update_candidate_space(obj,designSample,isInside,trimmingInformation)
        %UPDATE_CANDIDATE_SPACE Incorporate new corner anchors or samples
        %   UPDATE_CANDIDATE_SPACE merges existing anchors with new ones from
        %   trimmingInformation, checks for redundant anchors, and merges new 
        %   sample points. The inside definition is recalculated.
        %
        %   OBJ = OBJ.UPDATE_CANDIDATE_SPACE(DESIGNSAMPLE,ISINSIDE,TRIMMINGINFORMATION)
        %   updates the internal state of the object by appending anchors and 
        %   corner directions from TRIMMINGINFORMATION. It also checks for
        %   redundant anchors and merges them if needed. 
        %
        %   Inputs:
        %       - OBJ : CandidateSpaceCornerBoxRemoval
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %       - ISINSIDE : (nSample,1) logical
        %       - TRIMMINGINFORMATION : struct
        %
        %   Outputs:
        %       - OBJ : CandidateSpaceCornerBoxRemoval
        %
        %   See also generate_candidate_space, AnchorPoint, CornerDirection.

            if(isempty(obj.DesignSampleDefinition) || isempty(obj.AnchorPoint))
                obj = obj.generate_candidate_space(designSample,trimmingInformation);
                return;
            end

            if(~isempty(trimmingInformation))
                anchorPointNew = vertcat(trimmingInformation.Anchor);
                cornerDirectionNew = vertcat(trimmingInformation.CornerDirection);
                
                obj.AnchorPoint = [obj.AnchorPoint;anchorPointNew];
                obj.CornerDirection = [obj.CornerDirection;cornerDirectionNew];
    
                % verify if there are any redundant anchor points
                % -> region removed includes a different anchor with same corner direction
                if(obj.CheckRedundantTrimmingUpdate)
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
                end

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
            end

            % keep samples in inside/outside
            [~,iLowerBoundaryAll] = min(obj.DesignSampleDefinition,[],1);
            [~,iUpperBoundaryAll] = max(obj.DesignSampleDefinition,[],1);

            insideSample = obj.DesignSampleDefinition(obj.IsInsideDefinition,:);
            [~,iLowerBoundaryInside] = min(insideSample,[],1);
            [~,iUpperBoundaryInside] = max(insideSample,[],1);
            iBoundaryInside = convert_index_base(obj.IsInsideDefinition,[iLowerBoundaryInside,iUpperBoundaryInside]','backward');

            obj.DesignSampleDefinition = [...
                obj.DesignSampleDefinition([iLowerBoundaryAll,iUpperBoundaryAll,iBoundaryInside'],:);...
                designSample;...
                obj.AnchorPoint];
            if(obj.CheckDuplicatePointsUpdate)
                obj.DesignSampleDefinition = unique(obj.DesignSampleDefinition,'rows');
            end
            obj.IsInsideDefinition = obj.is_in_candidate_space(obj.DesignSampleDefinition,false);
        end
        
        function obj = expand_candidate_space(obj,growthRate)
        %EXPAND_CANDIDATE_SPACE Expansion of candidate space by given factor
        %   EXPAND_CANDIDATE_SPACE will grow the region considered inside the current 
        %   candidate space by the factor given. Said growth is done in a fixed rate 
        %   defined by the input relative to the design space.
        %   This is done by finding the center of the convex hull and then making all 
        %   inside designs move opposite to that direction. 
        %
        %   OBJ = OBJ.EXPAND_CANDIDATE_SPACE(GROWTHRATE) will growth the candidate space 
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
        %   See also generate_candidate_space, is_in_candidate_space.
            designSpaceFactor = obj.DesignSpaceUpperBound - obj.DesignSpaceLowerBound;
            designSpace = [obj.DesignSpaceLowerBound;obj.DesignSpaceUpperBound];
            
            if(obj.NormalizeGrowthDirection)
                designSpaceNormalization = designSpaceFactor;
            else
                designSpaceNormalization = 1;
            end

            center = mean(obj.DesignSampleDefinition(obj.IsInsideDefinition,:),1);
            distanceToCenter = (obj.DesignSampleDefinition - center)./designSpaceNormalization;
            normalizedDirectionGrowth = distanceToCenter./vecnorm(distanceToCenter,2,2);

            % positive direction - away from center
            directionGrowth = designSpaceFactor.*normalizedDirectionGrowth;
            maxGrowthRate = region_limit_line_search([],obj.DesignSampleDefinition,directionGrowth,designSpace);
            sampleGrowthRate = min(growthRate,maxGrowthRate);
            designSampleNewPositive = obj.DesignSampleDefinition + sampleGrowthRate.*directionGrowth;

            % negative direction - towards center
            directionGrowth = -designSpaceFactor.*normalizedDirectionGrowth;
            maxGrowthRate = region_limit_line_search([],obj.DesignSampleDefinition,directionGrowth,designSpace);
            sampleGrowthRate = min(growthRate,maxGrowthRate);
            designSampleNewNegative = obj.DesignSampleDefinition + sampleGrowthRate.*directionGrowth;

            % aggregate
            designSampleNew = min(max([designSampleNewPositive;designSampleNewNegative],obj.DesignSpaceLowerBound),obj.DesignSpaceUpperBound);
            obj.DesignSampleDefinition = [obj.DesignSampleDefinition;designSampleNew];

            if(~isempty(obj.AnchorPoint))
                isInsideDefinition = obj.is_in_candidate_space(obj.DesignSampleDefinition,false);
                insideSample = obj.DesignSampleDefinition(isInsideDefinition,:);
                nAnchor = size(obj.AnchorPoint,1);
                nDimension = size(obj.DesignSpaceLowerBound,2);
                directionGrowth = zeros(nAnchor,nDimension);
                for i=1:nAnchor
                    % check distances to anchor in each dimension (positive -> inside)
                    distanceToAnchor = (insideSample - obj.AnchorPoint(i,:))./designSpaceNormalization;
                    distanceToAnchor(:,obj.CornerDirection(i,:)) = -distanceToAnchor(:,obj.CornerDirection(i,:));
                    [~,iDimension] = max(distanceToAnchor,[],2);
                    for j=1:nDimension
                        directionGrowth(i,j) = sum(iDimension==j);
                    end
                end
                directionGrowth(~obj.CornerDirection) = -directionGrowth(~obj.CornerDirection);

                directionGrowth = directionGrowth./vecnorm(directionGrowth,2,2);
                anchorPointNew = obj.AnchorPoint + growthRate.*designSpaceFactor.*directionGrowth;
                anchorPointNew = min(max(anchorPointNew,obj.DesignSpaceLowerBound),obj.DesignSpaceUpperBound);
                
                % don't include anchors that were moved to the boundary corners
                isAnchorInLowerBoundary = (anchorPointNew<=obj.DesignSpaceLowerBound);
                isAnchorInUpperBoundary = (anchorPointNew>=obj.DesignSpaceUpperBound);
                isAnchorInCorner = all(isAnchorInLowerBoundary|isAnchorInUpperBoundary,2);

                obj.AnchorPoint = anchorPointNew(~isAnchorInCorner,:);
                obj.CornerDirection = obj.CornerDirection(~isAnchorInCorner,:);

                if(obj.CheckRedundantTrimmingGrowth)
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
                end

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
            end
            obj.DesignSampleDefinition = [obj.DesignSampleDefinition;obj.AnchorPoint];
            if(obj.CheckDuplicatePointsGrowth)
                obj.DesignSampleDefinition = unique(obj.DesignSampleDefinition,'rows');
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
            nSample = obj.MeasureEstimationFactor*size(obj.DesignSampleDefinition,1);
            samplingBox = obj.SamplingBox;
            
            volumeSample = sampling_random(samplingBox,nSample);
            isInside = obj.is_in_candidate_space(volumeSample);
            volumeFactor = sum(isInside) / size(isInside,1);
            volume = volumeFactor * prod(samplingBox(2,:) - samplingBox(1,:));
        end
    end
end