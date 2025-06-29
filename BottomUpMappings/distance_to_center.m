function performanceMeasure = distance_to_center(designSample,systemParameter)
%DISTANCE_TO_CENTER Bottom-up Mapping (n-Dimensional L2 Distance Test)
%   DISTANCE_TO_CENTER is a bottom-up mapping that returns the L2-distance of 
%   the given designs to a center defined as a parameter. When such measure is
%   given a (lower and) upper limit, the resulting regions of good designs is a
%   (hollow) sphere. Multiple centers can be defined, with the distance being
%   returned for each one.
%
%   PERFORMANCEMEASURE = DISTANCE_TO_CENTER(DESIGNSAMPLE) 
%   receives the spatial coordinates of the design points in DESIGNSAMPLE
%   and returns in PERFORMANCEMEASURE the L2 distance to the origin (0,0,...).
%
%   PERFORMANCEMEASURE = DISTANCE_TO_CENTER(DESIGNSAMPLE,SYSTEMPARAMETER) 
%   allows the specification of center(s) other than the origin in 
%   SYSTEMPARAMETER. Each center must be given as a new row in the array.
%
%   Input: 
%       - DESIGNSAMPLE : (nSample,d) double 
%           -- (i,j) : spatial coordinate j of each sample point i
%       - SYSTEMPARAMETER : (nCenter,d) double 
%           -- (i,j) : spatial coordinate j of reference center i
%
%   Output:
%       - PERFORMANCEMEASURE : (nSample,nCenter) double 
%           -- (i,j) : L2-norm distance i to reference center j
%
%   See also sphere_nd.
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

    % if center was not given, assume origin
    if(nargin<2 || isempty(systemParameter))
        systemParameter = zeros(1,size(designSample,2));
    end

    % calculate distance to each center
    nCenter = size(systemParameter,1);
    nSample = size(designSample,1);
    performanceMeasure = nan(nSample,nCenter);
    for i=1:nCenter
        performanceMeasure(:,i) = vecnorm(designSample-systemParameter(i,:),2,2);
    end
end