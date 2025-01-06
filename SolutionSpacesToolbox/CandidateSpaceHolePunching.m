classdef CandidateSpaceHolePunching < CandidateSpaceBase
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
        HoleSize

        %
        TrimmingApplicationSlack
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
        function obj = CandidateSpaceHolePunching(designSpaceLowerBound,designSpaceUpperBound,varargin)
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
            parser.addParameter('SamplingBoxSlack',0.5,@(x)isnumeric(x)&&isscalar(x)&&(x>=0)&&(x<=1));
            parser.addParameter('TrimmingApplicationSlack',0.01,@(x)isnumeric(x)&&isscalar(x)&&(x>0)&&(x<=1));
            parser.parse(designSpaceLowerBound,designSpaceUpperBound,varargin{:});

            
            obj.DesignSpaceLowerBound = parser.Results.DesignSpaceLowerBound;;
            obj.DesignSpaceUpperBound = parser.Results.DesignSpaceUpperBound;
            obj.SamplingBoxSlack = parser.Results.SamplingBoxSlack;
            obj.TrimmingApplicationSlack = parser.Results.TrimmingApplicationSlack;

            obj.DesignSampleDefinition = [];
            obj.AnchorPoint = [];
            obj.HoleSize = [];
        end
        
        function obj = generate_candidate_space(obj,designSample,trimmingInformation)
        %GENERATE_CANDIDATE_SPACE Initial definition of the candidate space
        %   GENERATE_CANDIDATE_SPACE uses labeled design samples to define the inside / 
        %   outside regions of the candidate space. For CandidateSpaceConvexHull, this
        %   means a convex hull is created around the inside designs.
        %
        %   OBJ = OBJ.GENERATE_CANDIDATE_SPACE(DESIGNSAMPLE) receives the design samle
        %   points in DESIGNSAMPLE and returns a candidate space object OBJ with the new
        %   definition, assuming all designs are inside the candidate space.
        %
        %   OBJ = OBJ.GENERATE_CANDIDATE_SPACE(DESIGNSAMPLE,ISINSIDE) additionally 
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
            obj.HoleSize = vertcat(trimmingInformation.HoleSize);
        end

        function obj = update_candidate_space(obj,designSample,isInside,trimmingInformation)
            if(isempty(obj.DesignSampleDefinition) || isempty(obj.AnchorPoint))
                obj = obj.generate_candidate_space(designSample,trimmingInformation);
                return;
            end

            if(isempty(trimmingInformation))
                return;
            end

            anchorPointNew = vertcat(trimmingInformation.Anchor);
            holeSizeNew = vertcat(trimmingInformation.HoleSize);
            
            isNewAnchor = [false(size(obj.AnchorPoint,1),1);true(size(anchorPointNew,1),1)];
            obj.AnchorPoint = [obj.AnchorPoint;anchorPointNew];
            obj.HoleSize = [obj.HoleSize;holeSizeNew];

            % keep samples in inside/outside
            [~,iLowerBoundaryAll] = min(obj.DesignSampleDefinition,[],1);
            [~,iUpperBoundaryAll] = max(obj.DesignSampleDefinition,[],1);

            isInsideDefinition = obj.IsInsideDefinition;
            insideSample = obj.DesignSampleDefinition(isInsideDefinition,:);
            [~,iLowerBoundaryInside] = min(insideSample,[],1);
            [~,iUpperBoundaryInside] = max(insideSample,[],1);
            iBoundaryInside = convert_index_base(isInsideDefinition,[iLowerBoundaryInside,iUpperBoundaryInside]','backward');

            % determine slack
            if(obj.TrimmingApplicationSlack<1 && any(isNewAnchor))
                % update to also include the new samples
                insideSample = [insideSample;designSample(isInside,:)];
                holeLowerBound = obj.AnchorPoint - obj.HoleSize./2;
                holeUpperBound = obj.AnchorPoint + obj.HoleSize./2;

                nDimension = size(insideSample,2);
                nAnchor = size(obj.AnchorPoint,1);
                maximumSlackUpperBound = nan(1,nDimension);
                maximumSlackLowerBound = nan(1,nDimension);
                for i=1:nAnchor
                    if(~isNewAnchor(i))
                        continue;
                    end

                    % move upper bound (if applicable)
                    distanceToUpperCorner = insideSample - holeUpperBound(i,:);
                    [maximumDistanceInside,iDimension] = max(distanceToUpperCorner,[],2);
                    for j=1:nDimension
                        distanceToBoundary = maximumDistanceInside(iDimension==j);
                        if(~isempty(distanceToBoundary))
                            maximumSlackUpperBound(j) = min(distanceToBoundary);
                        else
                            maximumSlackUpperBound(j) = obj.DesignSpaceUpperBound(j) - holeUpperBound(i,j);
                        end
                    end

                    % move lower bound (if applicable)
                    distanceToUpperCorner = holeLowerBound(i,:) - insideSample;
                    [maximumDistanceInside,iDimension] = max(distanceToUpperCorner,[],2);
                    for j=1:nDimension
                        distanceToBoundary = maximumDistanceInside(iDimension==j);
                        if(~isempty(distanceToBoundary))
                            maximumSlackLowerBound(j) = min(distanceToBoundary);
                        else
                            maximumSlackLowerBound(j) = holeLowerBound(i,j) - obj.DesignSpaceLowerBound(j);
                        end
                    end

                    % move boundaries
                    holeLowerBound(i,:) = holeLowerBound(i,:) - obj.TrimmingApplicationSlack*maximumSlackLowerBound;
                    holeUpperBound(i,:) = holeUpperBound(i,:) + obj.TrimmingApplicationSlack*maximumSlackUpperBound;
                end
                obj.AnchorPoint = (holeLowerBound+holeUpperBound)./2;
                obj.HoleSize = holeUpperBound - holeLowerBound;
            end

            obj.DesignSampleDefinition = unique(...
                [obj.DesignSampleDefinition([iLowerBoundaryAll,iUpperBoundaryAll,iBoundaryInside'],:);...
                designSample;...
                obj.AnchorPoint],'rows');
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

            isInsideDefinition = obj.IsInsideDefinition;
            center = mean(obj.DesignSampleDefinition(isInsideDefinition,:),1);
            distanceToCenter = obj.DesignSampleDefinition - center;
            directionGrowth = distanceToCenter./vecnorm(distanceToCenter,2,2);

            if(isempty(obj.AnchorPoint))
                insideFunction = [];
            else
                insideFunction = @(x)obj.is_in_candidate_space(x);
            end

            maxGrowthRate = region_limit_line_search(insideFunction,obj.DesignSampleDefinition,designSpaceFactor.*directionGrowth,designSpace);
            sampleGrowthRate = min(growthRate,maxGrowthRate);
            designSampleGrown = obj.DesignSampleDefinition + sampleGrowthRate.*designSpaceFactor.*directionGrowth;
            designSampleGrown = min(max(designSampleGrown,obj.DesignSpaceLowerBound),obj.DesignSpaceUpperBound);
            isInsideGrown = isInsideDefinition;

            maxGrowthRate = region_limit_line_search(insideFunction,obj.DesignSampleDefinition,-designSpaceFactor.*directionGrowth,designSpace);
            sampleGrowthRate = min(growthRate,maxGrowthRate);
            designSampleShrunk = obj.DesignSampleDefinition - sampleGrowthRate.*designSpaceFactor.*directionGrowth;
            designSampleShrunk = min(max(designSampleShrunk,obj.DesignSpaceLowerBound),obj.DesignSpaceUpperBound);
            isInsideShrunk = isInsideDefinition;

            isInsideNew = [isInsideDefinition;isInsideGrown;isInsideShrunk];
            [obj.DesignSampleDefinition,uniqueIndex] = unique([obj.DesignSampleDefinition;designSampleGrown;designSampleShrunk],'rows');
            isInsideNew = isInsideNew(uniqueIndex);

            if(~isempty(obj.AnchorPoint))
                holeLowerBound = obj.AnchorPoint - obj.HoleSize./2;
                holeUpperBound = obj.AnchorPoint + obj.HoleSize./2;

                isInside = obj.IsInsideDefinition;
                insideSample = obj.DesignSampleDefinition(isInside,:);
                nAnchor = size(obj.AnchorPoint,1);
                nDimension = size(obj.AnchorPoint,2);
                
                upperCornerMove = zeros(nAnchor,nDimension);
                lowerCornerMove = zeros(nAnchor,nDimension);
                for i=1:nAnchor
                    % move upper bound (if applicable)
                    distanceToUpperCorner = insideSample - holeUpperBound(i,:);
                    [maximumDistanceInside,iDimension] = max(distanceToUpperCorner,[],2);
                    for j=1:nDimension
                        distanceToBoundary = maximumDistanceInside(iDimension==j);
                        if(~isempty(distanceToBoundary))
                            upperCornerMove(i,j) = max(designSpaceFactor(j)*growthRate - min(distanceToBoundary),0);
                        end
                    end

                    % move lower bound (if applicable)
                    distanceToUpperCorner = holeLowerBound(i,:) - insideSample;
                    [maximumDistanceInside,iDimension] = max(distanceToUpperCorner,[],2);
                    for j=1:nDimension
                        distanceToBoundary = maximumDistanceInside(iDimension==j);
                        if(~isempty(distanceToBoundary))
                            lowerCornerMove(i,j) = max(designSpaceFactor(j)*growthRate - min(distanceToBoundary),0);
                        end
                    end
                end
                lowerCornerMove = max(designSpaceFactor.*growthRate.*lowerCornerMove./vecnorm(lowerCornerMove,2,2),0);
                upperCornerMove = max(designSpaceFactor.*growthRate.*upperCornerMove./vecnorm(upperCornerMove,2,2),0);
                holeLowerBound = holeLowerBound + lowerCornerMove;
                holeUpperBound = holeUpperBound - upperCornerMove;
                
                % remove holes where the boundaries overlap
                removedHole = any(holeLowerBound>=holeUpperBound,2);
                holeLowerBound = holeLowerBound(~removedHole,:);
                holeUpperBound = holeUpperBound(~removedHole,:);

                % remove any hole where designs that should be inside aren't
                designSampleInside = obj.DesignSampleDefinition(isInsideNew,:);
                nAnchorNew = size(holeLowerBound,1);
                removedHole = false(nAnchorNew,1);
                for i=1:nAnchorNew
                    anchorBox = [holeLowerBound(i,:);holeUpperBound(i,:)];
                    isInside = is_in_design_box(designSampleInside,anchorBox);
                    if(any(isInside))
                        removedHole(i) = true;
                    end
                end
                holeLowerBound = holeLowerBound(~removedHole,:);
                holeUpperBound = holeUpperBound(~removedHole,:);

                % update anchor points and hole sizes
                obj.AnchorPoint = (holeLowerBound+holeUpperBound)./2;
                obj.HoleSize = holeUpperBound - holeLowerBound;

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
            
            nAnchor = size(obj.AnchorPoint,1);
            holeLowerBound = obj.AnchorPoint - obj.HoleSize./2;
            holeUpperBound = obj.AnchorPoint + obj.HoleSize./2;
            for i=1:nAnchor
                anchorBox = [holeLowerBound(i,:);holeUpperBound(i,:)];
                [isInsideAnchor,scoreAnchor] = is_in_design_box(designSample,anchorBox);

                isInside(isInsideAnchor & scoreAnchor~=0) = false;
                score = max(score,-scoreAnchor);
            end
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

                %if(any(isInsideDefinition(1)~=isInsideDefinition))
                %    isInBoundary = design_find_boundary_samples(obj.DesignSampleDefinition,isInsideDefinition);
                %    isShapeDefinition = isShapeDefinition | isInBoundary;
                %end
                
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