function removalCandidate = component_trimming_method_corner_box_removal(designSampleComponent,iRemove,varargin)
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

    nDesignVariable = size(designSampleComponent,2);
    combination = logical(round(fullfact(2*ones(nDesignVariable,1)) - 1));

    designLesser = (designSampleComponent-designSampleComponent(iRemove,:)<=0);
    designGreater = (designSampleComponent-designSampleComponent(iRemove,:)>=0);

    nCombination = size(combination,1);
    nSample = size(designSampleComponent,1);
    removalCandidate = false(nSample,nCombination);
    for i=1:nCombination
        % See which points fullfill criteria
        combinationCandidateLesser = all(designLesser(:,~combination(i,:)),2);
        combinationCandidateGreater = all(designGreater(:,combination(i,:)),2);
        removalCandidate(:,i) = combinationCandidateLesser & combinationCandidateGreater;
    end
end