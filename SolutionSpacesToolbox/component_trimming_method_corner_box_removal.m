function [removalCandidate,removalInformation] = component_trimming_method_corner_box_removal(designSampleComponent,iRemove,isKeep,varargin)
%COMPONENT_TRIMMING_METHOD_CORNER_BOX_REMOVAL Component SSO Trimming
%   COMPONENT_TRIMMING_METHOD_CORNER_BOX_REMOVAL uses the corner box removal
%   method to find the sample points for candidate removal during the trimming 
%   operation of the component SSO procedure.
%   In corner box removal, for a component of n-dimensiones, there will be 2^n 
%   trimming possibilities, and each one will be a combination where the sample 
%   points that have design variable values of less tahn / greater than the 
%   removal anchor are selected.
%
%   REMOVALCANDIDATE = COMPONENT_TRIMMING_METHOD_CORNER_BOX_REMOVAL(
%   DESIGNSAMPLECOMPONENT,IREMOVE) receives the design sample points 
%   of the component space in DESIGNSAMPLECOMPONENT and the index of the design
%   to be removed IREMOVE, and returns all trimming possibilities in 
%   REMOVALCANDIDATE. 
%
%   Input:
%       - DESIGNSAMPLECOMPONENT : (nSample,nComponentDesignVariable) double
%       - IREMOVE : integer
%
%   Output:
%       - REMOVALCANDIDATE : (nSample,nCandidate) logical
%
%   See also component_trimming_operation, component_trimming_leanness, 
%   component_trimming_method_planar_trimming.
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

    parser = inputParser;
    parser.addParameter('CornersToTest','all');
    parser.addParameter('TrimmingSlack',0.5);
    parser.addParameter('ConsiderOnlyKeepInSlack',true);
    parser.parse(varargin{:});
    options = parser.Results;

    nDesignVariable = size(designSampleComponent,2);

    % create all combinations of (lesser than / greater than) for each design variable
    if(strcmpi(options.CornersToTest,'all'))
        combination = logical(round(fullfact(2*ones(nDesignVariable,1)) - 1));
    elseif(strcmpi(options.CornersToTest,'away'))
        center = mean(designSampleComponent(isKeep,:),1);
        combination = ((designSampleComponent(iRemove,:) - center)>=0);
    end

    distanceToAnchorBase = designSampleComponent-designSampleComponent(iRemove,:);

    nCombination = size(combination,1);
    nSample = size(designSampleComponent,1);
    removalCandidate = false(nSample,nCombination);
    anchorPoint = repmat(designSampleComponent(iRemove,:),nCombination,1);
    maximumSlack = nan(1,nDesignVariable);
    anchorSlack = nan(1,nDesignVariable);

    for i=1:nCombination
        % See which points fullfill criteria
        distanceToAnchor = distanceToAnchorBase;
        distanceToAnchor(:,combination(i,:)) = -distanceToAnchor(:,combination(i,:));
        removalCandidateCurrent = all(distanceToAnchor<0,2);
        anchorPoint = designSampleComponent(iRemove,:);

        if(options.TrimmingSlack<1)
            lowerBound = min(designSampleComponent,[],1);
            upperBound = max(designSampleComponent,[],1);

            isInsideSlack = ~removalCandidateCurrent;
            if(options.ConsiderOnlyKeepInSlack)
                isInsideSlack = isInsideSlack & isKeep;
            end
            insideToAnchorDistance = distanceToAnchor(isInsideSlack,:);

            [maximumSlackInside,iDimension] = max(insideToAnchorDistance,[],2);
            for j=1:nDesignVariable
                allowedSlack = maximumSlackInside(iDimension==j);
                if(isempty(allowedSlack))
                    if(~combination(i,j))
                        maximumSlack(j) = upperBound(j) - anchorPoint(j);
                    else
                        maximumSlack(j) = anchorPoint(j) - lowerBound(j);
                    end
                else
                    maximumSlack(j) = min(allowedSlack);
                end
            end
            anchorSlack(~combination(i,:)) = maximumSlack(~combination(i,:));
            anchorSlack(combination(i,:)) = -maximumSlack(combination(i,:));
            anchorPoint = anchorPoint + (1-options.TrimmingSlack)*anchorSlack;

            distanceToAnchor = designSampleComponent - anchorPoint;
            distanceToAnchor(:,combination(i,:)) = -distanceToAnchor(:,combination(i,:));
        end
        
        removalCandidate(:,i) = all(distanceToAnchor<0,2);
        if(nargout>1)
            removalInformation(i).Anchor = anchorPoint;
            removalInformation(i).CornerDirection = combination(i,:);
        end
    end
end