function performanceMeasure = springs_serial_parallel_total_stiffness(designSample,systemParameter)
%SPRINGS_SERIAL_PARALLEL_TOTAL_STIFFNESS Bottom-up Mapping (Spring System)
%   SPRINGS_SERIAL_PARALLEL_TOTAL_STIFFNESS calculates the total stiffness of
%   a system containing multiple parallel spring entries, where each parallel
%   entry consists of springs connected in series.
%   
%   Example configuration:
%   o-spring1-+-spring2-o
%   |
%   o-----spring3-------o
%   For this system, systemParameter would be {[1,2], [3]} indicating two
%   parallel entries (first with springs 1+2 in series, second with spring 3)
%
%   PERFORMANCEMEASURE = SPRINGS_SERIAL_PARALLEL_TOTAL_STIFFNESS(DESIGNSAMPLE,
%   SYSTEMPARAMETER) calculates the total stiffness of the spring system by:
%   1. Calculating equivalent stiffness for each parallel entry (series springs)
%   2. Summing stiffness of all parallel entries
%
%   Input:
%       - designSample : (nSample,nSpring) double
%       - systemParameter : cell array
%
%   Output:
%       - performanceMeasure : (nSample,1) double
%
%   See also SPRINGS_SERIAL_STIFFNESS, SPRINGS_PARALLEL_STIFFNESS.
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

    nParallelEntry = length(systemParameter);
    nSample = size(designSample,1);
    performanceMeasure = zeros(nSample,1);
    for i=1:nParallelEntry
        entryInverseStiffness = sum(1./designSample(:,systemParameter{i}),2);
        performanceMeasure = performanceMeasure + 1./entryInverseStiffness;
    end
end