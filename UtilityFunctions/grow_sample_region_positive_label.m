function [designSampleExpanded,labelExpanded] = grow_sample_region_positive_label(designSample,label,designSpaceLowerBound,designSpaceUpperBound,growthRate,varargin)
%GROW_SAMPLE_REGION_POSITIVE_LABEL Expansion of labeled region by given factor
%   GROW_SAMPLE_REGION_POSITIVE_LABEL will expand the region labeled as positive 
%   by the factor given. Said growth is done in a fixed rate defined by the 
%   input relative to the design space. This is an isotropic expansion.
%
%   [DESIGNSAMPLEEXPANDED,LABELEXPANDED] = GROW_SAMPLE_REGION_POSITIVE_LABEL(
%   DESIGNSAMPLE,LABEL,DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND,GROWTHRATE)
%   takes the design sample points DESIGNSAMPLE which are labeled LABEL,
%   the design space DESIGNSPACELOWERBOUND and DESIGNSPACEUPPERBOUND, and the
%   growth rate factor GROWTHRATE, and produces new sample points 
%   DESIGNSAMPLEEXPANDED which are labeled LABELEXPANDED depending on their
%   distances to the original positive-label points. The expansion is done both
%   outwards and inwards.
%
%   [...] = GROW_SAMPLE_REGION_POSITIVE_LABEL(...,NAME,VALUE,...) allows the 
%   specification of additional arguments; here, these arguments are used when
%   computing distances for defining the growth with 'knnsearch'.
%
%   Inputs:
%       - DESIGNSAMPLE : (nSample,nDesignVariable) double
%       - LABEL : (nSample,1) logical
%       - DESIGNSPACELOWERBOUND : (1,nDesignVariable) double
%       - DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double
%       - GROWTHRATE : double
%       - Name-value pair arguments: passed to 'knnsearch'
%
%   Output:
%       - DESIGNSAMPLEEXPANDED : (3*nSample,nDesignVariable) double
%       - LABELEXPANDED : (3*nSample,1) logical
%
%   See also: knnsearch, CandidateSpaceSvm, CandidateSpaceDecisionTree.
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

	% grow designs in training data
	positiveSample = designSample(label,:);
    centerPositive = mean(positiveSample,1);
    
    % find direction fo growth based on center of inside designs
    distanceToCenter = designSample - centerPositive;
    directionGrowth = distanceToCenter./vecnorm(distanceToCenter,1,2);
    
    % prepare for growth operation
    designSpaceFactor = designSpaceUpperBound - designSpaceLowerBound;
    designSpace = [designSpaceLowerBound;designSpaceUpperBound];

    % grow the current samples
    maxGrowthRate = region_limit_line_search([],designSample,designSpaceFactor.*directionGrowth,designSpace);
    sampleGrowthRate = min(growthRate,maxGrowthRate);
    grownSample = designSample + sampleGrowthRate.*designSpaceFactor.*directionGrowth;

    % shrink the current samples
    maxGrowthRate = region_limit_line_search([],designSample,designSpaceFactor.*(-directionGrowth),designSpace);
    sampleGrowthRate = min(growthRate,maxGrowthRate);
    shrunkSample = designSample + sampleGrowthRate.*designSpaceFactor.*(-directionGrowth);
    
    % limit new samples to design space
    designSampleExpanded = [designSample;grownSample;shrunkSample];
    designSampleExpanded = max(designSampleExpanded, designSpaceLowerBound); % lower bound limit
    designSampleExpanded = min(designSampleExpanded, designSpaceUpperBound); % upper bound limit
    designSampleExpanded = unique(designSampleExpanded,'rows');
    
    % see which samples are within a distance of growth rate of the previous inside region
    referenceSample = positiveSample./designSpaceFactor;
    querySample = designSampleExpanded./designSpaceFactor;
    % [~,distanceInside] = dsearchn(referenceSample,querySample);
    [~,distanceInside] = knnsearch(referenceSample,querySample,...
        'K',1,'IncludeTies',false,varargin{:});
    labelExpanded = (distanceInside <= growthRate);
end