classdef CandidateSpacePlanarTrimming < CandidateSpaceBase
%CANDIDATESPACEPLANARTRIMMING Candidate Space defined by planar trimming
%   CANDIDATESPACEPLANARTRIMMING can be used to define candidate spaces for  
%   component solution space computation. In this case, the candidate space is 
%   defined by a set of planes that successively trim the design space, labeling
%   regions as inside or outside. 
%
%   This class is derived from CandidateSpaceBase.
%
%   CANDIDATESPACEPLANARTRIMMING properties:
%       - DesignSampleDefinition : design sample points used in the candidate 
%       space definition.
%       - IsInsideDefinition : logical labels of the design sample points used 
%       in the candidate space definition regarding whether they are inside or 
%       outside the candidate space (true = inside, false = outside).
%       - AnchorPoint : reference anchor points (one per plane) used for 
%       trimming.
%       - PlaneOrientationAnchor : orientation vectors of each plane that define
%       the direction pointing "inside".
%       - IsShapeDefinition : logical labels of the design sample points that 
%       actively define the space boundary (e.g., extreme samples, anchors).
%       - Measure : measure of the candidate space (e.g., approximate area or 
%       volume).
%       - SamplingBox : bounding box around the internal region of the candidate
%       space that can be used to help with sampling.
%
%   CANDIDATESPACEPLANARTRIMMING methods:
%       - generate_candidate_space : create a candidate space based on design 
%       samples and plane definitions (via trimmingInformation).
%       - update_candidate_space : re-generate or refine it with new data.
%       - expand_candidate_space : expand the candidate space by a given factor.
%       - is_in_candidate_space : verify if given design sample points are  
%       inside the candidate space.
%       - plot_candidate_space : visualize 1D/2D/3D candidate spaces in the 
%       specified figure.
%
%   See also CandidateSpaceBase.
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

        %ANCHORPOINT Reference anchor points used for trimming planes
        %   ANCHORPOINT collects reference anchor points that define each plane.
        %   One plane is defined per anchor point, and the plane orientation 
        %   in PlaneOrientationAnchor indicates the region considered "inside."
        %
        %   ANCHORPOINT : (nPlane, nDesignVariable) double
        %
        %   See also PlaneOrientationAnchor, is_in_candidate_space.
        AnchorPoint

        %PLANEORIENTATIONANCHOR Orientation vectors of each trimming plane
        %   PLANEORIENTATIONANCHOR contains orientation vectors for each plane 
        %   that define which side of the plane is considered inside. Each row
        %   corresponds to the orientation vector of a plane. 
        %
        %   PLANEORIENTATIONANCHOR : (nPlane, nDesignVariable) double
        %
        %   See also AnchorPoint, is_in_candidate_space.
        PlaneOrientationAnchor
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
        function obj = CandidateSpacePlanarTrimming(designSpaceLowerBound,designSpaceUpperBound,varargin)
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
            parser.parse(designSpaceLowerBound,designSpaceUpperBound,varargin{:});

            
            obj.DesignSpaceLowerBound = parser.Results.designSpaceLowerBound;;
            obj.DesignSpaceUpperBound = parser.Results.designSpaceUpperBound;

            obj.DesignSampleDefinition = [];
            obj.AnchorPoint = [];
            obj.PlaneOrientationAnchor = [];
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
        %   OBJ = OBJ.GENERATE_CANDIDATE_SPACE(DESIGNSAMPLE,TRIMMINGINFORMATION) uses
        %   information from the trimming operation TRIMMINGINFORMATION to gather the 
        %   considered anchor points and plane orientations, using that to determine 
        %   the inside region.
        %
        %   Inputs:
        %       - OBJ : CandidateSpaceConvexHull
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %       - TRIMMINGINFORMATION : struct
        %   
        %   Outputs:
        %       - OBJ : CandidateSpaceConvexHull
        %   
        %   See also convex_hull_plane, is_in_candidate_space.
            obj.DesignSampleDefinition = designSample;

            if(~isempty(trimmingInformation))
                obj.AnchorPoint = vertcat(trimmingInformation.Anchor);
                obj.PlaneOrientationAnchor = vertcat(trimmingInformation.PlaneOrientationInside);
            end
            obj.IsInsideDefinition = obj.is_in_candidate_space(designSample,false);
        end

        function obj = update_candidate_space(obj,designSample,isInside,trimmingInformation)
        %UPDATE_CANDIDATE_SPACE Re-generate or refine candidate space with new data
        %   UPDATE_CANDIDATE_SPACE can receive new samples and new trimming
        %   planes, then adjust the candidate space definition accordingly.
        %
        %   OBJ = OBJ.UPDATE_CANDIDATE_SPACE(DESIGNSAMPLE,ISINSIDE,TRIMMINGINFORMATION)
        %   merges previous definitions with the new data. If no anchor info or
        %   plane orientation is given, the trimming is not changed. 
        %
        %   Inputs:
        %       - OBJ : CandidateSpacePlanarTrimming
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %       - ISINSIDE : (nSample,1) logical
        %       - TRIMMINGINFORMATION : struct
        %
        %   Outputs:
        %       - OBJ : CandidateSpacePlanarTrimming
        %
        %   See also generate_candidate_space, AnchorPoint, PlaneOrientationAnchor.

            if(isempty(obj.DesignSampleDefinition) || isempty(obj.AnchorPoint))
                obj = obj.generate_candidate_space(designSample,trimmingInformation);
                return;
            end

            if(isempty(trimmingInformation))
                return;
            end
            anchorPointNew = vertcat(trimmingInformation.Anchor);
            planeOrientationNew = vertcat(trimmingInformation.PlaneOrientationInside);

            obj.AnchorPoint = [obj.AnchorPoint;anchorPointNew];
            obj.PlaneOrientationAnchor = [obj.PlaneOrientationAnchor;planeOrientationNew];

            % project points to planes - if no point is inside, plane is redundant and can be removed
            nAnchor = size(obj.AnchorPoint,1);
            hullPointMinMax = [];
            isRedundantAnchor = false(nAnchor,1);
            for i=1:nAnchor
                dotProduct = sum((obj.AnchorPoint(i,:) - obj.DesignSampleDefinition).*obj.PlaneOrientationAnchor(i,:),2);
                distanceToPlane = obj.PlaneOrientationAnchor(i,:).*dotProduct;
                candidateHullPoint = obj.DesignSampleDefinition + distanceToPlane;

                isInside = obj.is_in_candidate_space(candidateHullPoint,false);
                if(~any(isInside))
                    isRedundantAnchor(i) = true;
                end

                hullPoint = [hullPointMinMax;candidateHullPoint(isInside,:)];
                [~,iHullPointMin] = min(hullPoint,[],1);
                [~,iHullPointMax] = max(hullPoint,[],1);
                hullPointMinMax = hullPoint([iHullPointMin,iHullPointMax],:);
            end
            obj.AnchorPoint(isRedundantAnchor,:) = [];
            obj.PlaneOrientationAnchor(isRedundantAnchor,:) = [];
            
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
                hullPointMinMax;...
                obj.AnchorPoint],'rows');
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

            center = mean(obj.DesignSampleDefinition(obj.IsInsideDefinition,:),1);
            distanceToCenter = obj.DesignSampleDefinition - center;
            directionGrowth = distanceToCenter./vecnorm(distanceToCenter,2,2);

            maxGrowthRate = region_limit_line_search([],obj.DesignSampleDefinition,designSpaceFactor.*directionGrowth,designSpace);
            sampleGrowthRate = min(growthRate,maxGrowthRate);
            designSampleNew = obj.DesignSampleDefinition + sampleGrowthRate.*designSpaceFactor.*directionGrowth;
            designSampleNew = min(max(designSampleNew,obj.DesignSpaceLowerBound),obj.DesignSpaceUpperBound);
            obj.DesignSampleDefinition = [obj.DesignSampleDefinition;designSampleNew];
            
            hullPointMinMax = [];
            if(~isempty(obj.AnchorPoint))
                %directionGrowth = obj.AnchorPoint - center;
                directionGrowth = -obj.PlaneOrientationAnchor;
                
                directionGrowth = directionGrowth./vecnorm(directionGrowth,2,2);
                anchorPointNew = obj.AnchorPoint + growthRate.*designSpaceFactor.*directionGrowth;
                anchorPointNew = min(max(anchorPointNew,obj.DesignSpaceLowerBound),obj.DesignSpaceUpperBound);

                isAnchorInLowerBoundary = (anchorPointNew==obj.DesignSpaceLowerBound);
                isAnchorInUpperBoundary = (anchorPointNew==obj.DesignSpaceUpperBound);
                isAnchorInCorner = all(isAnchorInLowerBoundary|isAnchorInUpperBoundary,2);
                obj.AnchorPoint = anchorPointNew(~isAnchorInCorner,:);
                obj.PlaneOrientationAnchor = obj.PlaneOrientationAnchor(~isAnchorInCorner,:);

                % project points to planes - if no point is inside, plane is redundant and can be removed
                nAnchor = size(obj.AnchorPoint,1);
                isRedundantAnchor = false(nAnchor,1);
                hullPointMinMax = [];
                for i=1:nAnchor
                    dotProduct = sum((obj.AnchorPoint(i,:) - obj.DesignSampleDefinition).*obj.PlaneOrientationAnchor(i,:),2);
                    distanceToPlane = obj.PlaneOrientationAnchor(i,:).*dotProduct;
                    candidateHullPoint = obj.DesignSampleDefinition + distanceToPlane;

                    isInside = obj.is_in_candidate_space(candidateHullPoint,false);
                    if(~any(isInside))
                        isRedundantAnchor(i) = true;
                    end
                    
                    hullPoint = [hullPointMinMax;candidateHullPoint(isInside,:)];
                    [~,iHullPointMin] = min(hullPoint,[],1);
                    [~,iHullPointMax] = max(hullPoint,[],1);
                    hullPointMinMax = hullPoint([iHullPointMin,iHullPointMax],:);
                end
                obj.AnchorPoint(isRedundantAnchor,:) = [];
                obj.PlaneOrientationAnchor(isRedundantAnchor,:) = [];
            end
            obj.DesignSampleDefinition = unique([obj.DesignSampleDefinition;hullPointMinMax;obj.AnchorPoint],'rows');
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
            if(nAnchor<=nSample)
                for i=1:nAnchor
                    dotProduct = sum((obj.AnchorPoint(i,:) - designSample).*obj.PlaneOrientationAnchor(i,:),2);

                    isInside(dotProduct>0) = false;
                    score = max(score,dotProduct);
                end
            else
                for i=1:nSample
                    dotProduct = dot(obj.AnchorPoint - designSample(i,:),obj.PlaneOrientationAnchor,2);

                    isInside(i) = all(dotProduct<=0);
                    score(i) = max(dotProduct);
                end
            end

            % check if it's inside the bounding box
            if(nargin<3 || includeBoundingBox)
                boundingBox = design_bounding_box(obj.DesignSampleDefinition,obj.IsInsideDefinition);
                boundingBox = min(max(boundingBox,obj.DesignSpaceLowerBound),obj.DesignSpaceUpperBound);
                [isInsideBounding,scoreBounding] = is_in_design_box(designSample,boundingBox);
                isInside(~isInsideBounding) = false;
                score = max(score,scoreBounding);
            end
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

            nDimension = size(obj.DesignSampleDefinition,2);
            if(nDimension>3)
                if(nargout>0)
                    plotHandle = [];
                end
                return;
            end

            hullPoint = [];
            nAnchor = size(obj.AnchorPoint,1);
            % project the points into each plane
            for i=1:nAnchor
                dotProduct = sum((obj.AnchorPoint(i,:) - obj.DesignSampleDefinition).*obj.PlaneOrientationAnchor(i,:),2);
                distanceToPlane = obj.PlaneOrientationAnchor(i,:).*dotProduct;
                candidateHullPoint = obj.DesignSampleDefinition + distanceToPlane;
                candidateHullPoint = min(max(candidateHullPoint,obj.DesignSpaceLowerBound),obj.DesignSpaceUpperBound);

                isInside = obj.is_in_candidate_space(candidateHullPoint);
                hullPoint = unique([hullPoint;candidateHullPoint(isInside,:)],'rows');
            end

            % create a convex hull based on that
            convexHullIndex = compute_convex_hull(hullPoint);

            % plot convex hull
            if(nDimension==1)
                plotHandle = plot_convex_hull_1d(figureHandle,hullPoint,convexHullIndex,varargin{:});
            elseif(nDimension==2)
                plotHandle = plot_convex_hull_2d(figureHandle,hullPoint,convexHullIndex,varargin{:});
            elseif(nDimension==3)
                plotHandle = plot_convex_hull_3d(figureHandle,hullPoint,convexHullIndex,varargin{:});
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
            nSample = 10*size(obj.DesignSampleDefinition,1);
            samplingBox = obj.SamplingBox;
            
            volumeSample = sampling_random(samplingBox,nSample);
            isInside = obj.is_in_candidate_space(volumeSample);
            volumeFactor = sum(isInside) / size(isInside,1);
            volume = volumeFactor * prod(samplingBox(2,:) - samplingBox(1,:));
        end
    end
end