function iChoice = component_trimming_cost(designSample,isViable,isExclude,isInsideComponent,isInsideAll,removalCandidate,ineligibleCandidate,varargin)
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
    parser.addParameter('SelectionCriteria',{'NumberInsideViable','NumberInsideExclude','NumberExclude','NumberViable','NumberInsideNotExclude'});
    parser.parse(varargin{:});
    options = parser.Results;
    
    nCandidate = size(removalCandidate,2);
    removalCost = nan(1,nCandidate);
    nCriteria = length(options.SelectionCriteria);

    if(all(ineligibleCandidate))
        isTie = true(1,nCandidate);
    else
        isTie = ~ineligibleCandidate;
    end
    nTie = sum(isTie);

    i = 1;
    while(nTie>1 && i<=nCriteria)
        criterionCurrent = options.SelectionCriteria{i};
        removalCandidateCurrent = removalCandidate(:,isTie);

        if(isa(criterionCurrent,'function_handle'))
            removalCost = options.CostType(designSample,isViable,isExclude,isInsideComponent,isInsideAll,removalCandidateCurrent);
        elseif(strcmpi(criterionCurrent,'NumberInsideViable'))
            removalCost = -sum(isViable & isInsideAll & ~removalCandidateCurrent,1);
        elseif(strcmpi(criterionCurrent,'NumberViable'))
            removalCost = -sum(isViable & isInsideComponent & ~removalCandidateCurrent,1);
        elseif(strcmpi(criterionCurrent,'NumberInsideExclude'))
            removalCost = -sum(isExclude & isInsideAll & removalCandidateCurrent,1);
        elseif(strcmpi(criterionCurrent,'NumberExclude'))
            removalCost = -sum(isExclude & isInsideComponent & removalCandidateCurrent,1);
        elseif(strcmpi(criterionCurrent,'NumberInsideNotExclude'))
            removalCost = -sum(isInsideComponent & ~isExclude & ~removalCandidateCurrent,1);
        end

        isTie(isTie) = (removalCost==min(removalCost));
        nTie = sum(isTie);
        i = i+1;
    end
    iChoice = find(isTie,1);
end

