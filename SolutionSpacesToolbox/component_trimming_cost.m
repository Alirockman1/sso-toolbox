function [removalCost,iChoice] = component_trimming_cost(designSample,isKeep,isRemove,isRemain,removalCandidate,ineligibleCandidate,varargin)
%COMPONENT_TRIMMING_COST Cost of removing selected designs from component space
%   COMPONENT_TRIMMING_COST computes the cost associated with performing the
%   candidate trimming operations on the component space. 
%
%   REMOVALCOST = COMPONENT_TRIMMING_COST(DESIGNSAMPLE,ACTIVEKEEP,
%   REMOVALCANDIDATE) receives the component design sample points in 
%   DESIGNSAMPLE, the label for which points should be kept on ACTIVEKEEP, and
%   all the trimming operation candidates for this component in 
%   REMOVALCANDIDATE, and returns the cost for each of those removal candidates
%   in REMOVALCOST. This cost is defined by the amount of points that are 
%   supposed to be kept that are being removed.
%
%   REMOVALCOST = COMPONENT_TRIMMING_COST(...NAME,VALUE,...) allows the 
%   specification of name-value pair arguments. These can be:
%       - 'CostType' : how to compute the cost of each removal candidate option. 
%       This can be one of the following:
%           -- 'NumberKeep' : number of design points labeled as 'keep' which 
%           are being removed.
%           -- 'VolumeKeep' : estimated volume of the desired region being 
%           removed.
%       Alternatively, a function handle can be used to define custom weighting.
%       This function must have the form 
%       'cost = f(designSample,isKeep,removalCandidate)'. 
%
%   Input:
%       - DESIGNSAMPLE : (nSample,nDesignVariable) double
%       - ACTIVEKEEP : (nSample,1) logical
%       - REMOVALCANDIDATE : (nSample,nRemovalCandidate) logical
%       - 'CostType' : char OR string OR function_handle
%
%   Output:
%       - REMOVALCOST : (1,nRemovalCandidate) double
%
%   See also component_trimming_choice, component_trimming_operation.
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

    parser = inputParser;
    parser.addParameter('CostType','NumberKeep');
    parser.addParameter('TieBreakerType','NumberRemoved');
    parser.parse(varargin{:});
    options = parser.Results;
    
    nRemovalCandidate = size(removalCandidate,2);
    removalCost = nan(1,nRemovalCandidate);
    if(isa(options.CostType,'function_handle'))
        removalCost = options.CostType(designSample,isKeep,removalCandidate);
    elseif(strcmpi(options.CostType,'RemovedKeepVolume'))
        nSample = size(designSample,1);
        keepFraction = sum(removalCandidate & isKeep,1)./nSample;

        for i=1:nRemovalCandidate
            boundingBoxRemoval = design_bounding_box(designSample,removalCandidate(:,i));
            totalVolumeRemoved = prod(boundingBoxRemoval(2,:) - boundingBoxRemoval(1,:));
            removalCost(i) = keepFraction(i) * totalVolumeRemoved;
        end
    elseif(strcmpi(options.CostType,'NumberKeep'))
        % number of designs being kept
        removalCost = -sum(~removalCandidate & isKeep & isRemain & ~isRemove,1);
    else
        % estimated volume being kept
        removalCost = nan(1,nRemovalCandidate);
        for i=1:nRemovalCandidate
            isKeepMaintain = isKeep & isRemain & ~removalCandidate(:,i);
            if(sum(isKeepMaintain)>0)
                boundingBoxRemainKeep = design_bounding_box(designSample,isKeep & isRemain & ~removalCandidate(:,i));
                pointInsideBox = is_in_design_box(designSample,boundingBoxRemainKeep);
                volumeFraction = sum(pointInsideBox & isKeep & isRemain & ~removalCandidate(:,i))/sum(pointInsideBox);
                removalCost(i) = -volumeFraction*prod(boundingBoxRemainKeep(2,:) - boundingBoxRemainKeep(1,:));
            else
                removalCost(i) = 0;
            end
        end
    end

    removalCost(ineligibleCandidate) = inf;
    [sortedCost,iCost] = sort(removalCost);
    iTieBreaker = (removalCost==sortedCost(1));

    nTieBreaker = sum(iTieBreaker);
    if(nTieBreaker>1)
        if(strcmpi(options.TieBreakerType,'NumberRemain'))
            % tie-breaker: total amount of points eliminated
            removalCostTieBreaker = -sum(isRemain & ~isRemove & ~removalCandidate(:,iTieBreaker),1);
        elseif(strcmpi(options.TieBreakerType,'NumberRemoved'))
            removalCostTieBreaker = -sum(isRemain & isRemove & removalCandidate(:,iTieBreaker),1);
        else%if(strcmpi(options.TieBreakerType,'VolumeRemain'))
            % tie-breaker: volume that remains
            removalCostTieBreaker = nan(1,nTieBreaker);
            for i=1:nTieBreaker
                iCurrent = convert_index_base(iTieBreaker',i,'backward');
                boundingBoxRemain = design_bounding_box(designSample,isRemain & ~removalCandidate(:,iCurrent));
                pointInsideBox = is_in_design_box(designSample,boundingBoxRemain);
                volumeFraction = sum(pointInsideBox & isRemain & ~removalCandidate(:,iCurrent))/sum(pointInsideBox);
                removalCostTieBreaker(i) = -volumeFraction*prod(boundingBoxRemain(2,:) - boundingBoxRemain(1,:));
            end
        end

        [~,iCostTieBreaker] = sort(removalCostTieBreaker);
        iChoice = convert_index_base(iTieBreaker',iCostTieBreaker(1),'backward');
    else
        iChoice = iCost(1);
    end
    removalCost = removalCost(iChoice);
end

