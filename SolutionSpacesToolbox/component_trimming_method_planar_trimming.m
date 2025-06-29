function removalCandidate = component_trimming_method_planar_trimming(designSampleComponent,iRemove,isKeep,varargin)
%COMPONENT_TRIMMING_METHOD_PLANAR_TRIMMING Component SSO Convex Trimming
%   COMPONENT_TRIMMING_METHOD_PLANAR_TRIMMING uses the planar trimming method to
%   find the sample points for candidate removal during the trimming operation 
%   of the component SSO procedure.
%   In planar trimming, the removal operations are done by planes passing 
%   through the removal point and which are perpendicular to all reference 
%   design points (which can be, for example, those that are desired to be 
%   kept). Using this trimming method for the component SSO trimming operation  
%   will always result in convex spaces.
%
%   REMOVALCANDIDATE = COMPONENT_TRIMMING_METHOD_PLANAR_TRIMMING(
%   DESIGNSAMPLECOMPONENT,IREMOVE,ISKEEP) receives the design sample points 
%   of the component space in DESIGNSAMPLECOMPONENT, the index of the design
%   to be removed IREMOVE, the label for designs that one wants to keep ISKEEP,
%   and returns all trimming possibilities in REMOVALCANDIDATE. In this case,
%   the reference designs are all designs to be kept.
%
%   REMOVALCANDIDATE = COMPONENT_TRIMMING_METHOD_PLANAR_TRIMMING(...NAME,VALUE,
%   ...) allows the specification of additional parameters. In this case:
%       - 'ReferenceDesigns' : which designs should be used for finding the
%       planes used for trimming. The following options are available:
%           -- 'all' : all design sample points are used (even ones which
%           one does not necessarily want to keep).
%           -- 'keep' : designs labeled as being desired to be kept.
%           -- 'boundary-center' : designs labeled as keep which are at the 
%           limits for each design variable (maximum/minimum values) and the
%           computed center of designs to keep.
%           -- 'center' : only the computed center of designs to keep.
%
%   Input:
%       - DESIGNSAMPLECOMPONENT : (nSample,nComponentDesignVariable) double
%       - IREMOVE : integer
%       - ISKEEP : (nSample,1) logical
%       - 'ReferenceDesigns' : char OR string
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

    parser = inputParser;
    parser.addParameter('ReferenceDesigns','keep');
    parser.parse(varargin{:});
    options = parser.Results;

    [upperValueDesignVariable,upperBoundaryDesignIndexKeep] = max(designSampleComponent(isKeep,:),[],1);
    [lowerValueDesignVariable,lowerBoundaryDesignIndexKeep] = min(designSampleComponent(isKeep,:),[],1);
    
    % normalize values for better numerical behavior
    designSampleComponent = designSampleComponent./(upperValueDesignVariable-lowerValueDesignVariable);

    if(strcmpi(options.ReferenceDesigns,'all'))
        designReference = designSampleComponent;
        designReference(iRemove,:) = [];
    elseif(strcmpi(options.ReferenceDesigns,'keep'))
        designReference = designSampleComponent(isKeep,:);
    elseif(strcmpi(options.ReferenceDesigns,'boundary-center'))
        boundaryIndexKeep = unique([upperBoundaryDesignIndexKeep,lowerBoundaryDesignIndexKeep])';
        boundaryIndex = convert_index_base(isKeep,boundaryIndexKeep,'backward');
        designReference = designSampleComponent(boundaryIndex,:);
        designReference = [designReference;mean(designSampleComponent(isKeep,:),1)];
    elseif(strcmpi(options.ReferenceDesigns,'center'))
        designReference = mean(designSampleComponent(isKeep,:),1);
    end

    % find all distances
    distanceAll = designSampleComponent(iRemove,:) - designSampleComponent;
    normalizedDistanceAll = distanceAll./vecnorm(distanceAll,2,2);

    distanceReference = designSampleComponent(iRemove,:) - designReference;
    normalizedDistanceReference = distanceReference./vecnorm(distanceReference,2,2);

    % for each plane, points being removed are all whose dot product between the 
    % distance to the anchor and each normal is non-positive
    % note: full product could be written as a matrix multiplication, but in my
    % experience, that actually makes the performance worse
    nSample = size(designSampleComponent,1);
    nRemovalCandidate = size(normalizedDistanceReference,1);
    removalCandidate = false(nSample,nRemovalCandidate);
    for i=1:nRemovalCandidate
        dotProduct = sum(normalizedDistanceReference(i,:).*normalizedDistanceAll,2);
        removalCandidate(:,i) = (dotProduct<=0);
    end
end