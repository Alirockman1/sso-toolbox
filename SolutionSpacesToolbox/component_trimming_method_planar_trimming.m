function removalCandidate = component_trimming_method_planar_trimming(designSampleComponent,iRemove,iKeep,varargin)
%COMPONENT_TRIMMING_METHOD_PLANAR_TRIMMING Component SSO Convex Trimming
%   COMPONENT_TRIMMING_METHOD_PLANAR_TRIMMING uses the planar trimming method to
%   find the sample points for candidate removal during the trimming operation 
%   of the component SSO procedure.
%   In planar trimming, the removal operations are done by planes passing 
%   through the removal point and which are perpendicular to all design points
%   that are desired to be kept. This means there are always a number of
%   removal candidates equal to the number of points that should be kept. Using
%   this trimming method for the component SSO trimming operation will always 
%   result in convex spaces.
%
%   REMOVALCANDIDATE = COMPONENT_TRIMMING_METHOD_PLANAR_TRIMMING(
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
%   component_trimming_method_corner_box_removal.
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
    
    % normalize values for better numerical behavior
    designSampleComponent = designSampleComponent./(max(designSampleComponent,[],1)-min(designSampleComponent,[],1));

    % find all distances
    distanceAll = designSampleComponent(iRemove,:) - designSampleComponent;
    normalizedDistanceAll = distanceAll./vecnorm(distanceAll,2,2);
    normalizedDistanceKeep = normalizedDistanceAll(iKeep,:);

    % for each plane, points being removed are all whose dot product between the 
    % distance to the anchor and each normal is non-positive
    % note: full product could be written as a matrix multiplication, but in my
    % experience, that actually makes the performance worse
    nSample = size(designSampleComponent,1);
    nRemovalCandidate = size(normalizedDistanceKeep,1);
    removalCandidate = false(nSample,nRemovalCandidate);
    for i=1:nRemovalCandidate
        dotProduct = sum(normalizedDistanceKeep(i,:).*normalizedDistanceAll,2);
        removalCandidate(:,i) = (dotProduct<=0);
    end
end