function removalCandidate = component_trimming_method_hole_punching(designSampleComponent,iRemove,varargin)
%COMPONENT_TRIMMING_METHOD_HOLE_PUNCHING Component SSO Trimming
%   COMPONENT_TRIMMING_METHOD_HOLE_PUNCHING uses the hole punching method to 
%   find the sample points for candidate removal during the trimming operation 
%   of the component SSO procedure.
%   In hole punching, all design points within a certain distance of the 
%   removal point are trimmed.
%
%   REMOVALCANDIDATE = COMPONENT_TRIMMING_METHOD_HOLE_PUNCHING(
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
%       - REMOVALCANDIDATE : (nSample,1) logical
%
%   See also component_trimming_operation, component_trimming_leanness, 
%   component_trimming_method_planar_trimming, 
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

    nSample = size(designSampleComponent,1);
    nDesignVariable = size(designSampleComponent,2);
    designSampleMax = max(designSampleComponent,[],1);
    designSampleMin = min(designSampleComponent,[],1);
    holeSize = (designSampleMax-designSampleMin)*((1/sqrt(nSample))^(1/nDesignVariable));
    
    removalCandidate = all(abs(designSampleComponent-designSampleComponent(iRemove,:))<=holeSize./2,2);
end