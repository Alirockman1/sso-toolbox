function designBoxGrown = design_box_grow_fixed(designBox,designSpaceLowerBound,designSpaceUpperBound,growthRate)
%DESIGN_BOX_GROW_FIXED Grow design box with fixed rate relative to design space
%   DESIGN_BOX_GROW_FIXED extends the lengths of a design box to each side
%   by a given fixed rate which is proportional to the size of the given design 
%   space. The edges/limits of the box, however, always stay within the 
%   boundaries of the design space.
%
%   DESIGNBOXGROWN = DESIGN_BOX_GROW_FIXED(DESIGNBOX,DESIGNSPACELOWERBOUND,
%   DESIGNSPACEUPPERBOUND,GROWTHRATE) receives a design box DESIGNBOX, the 
%   lower/upper boundaries of the design space DESIGNSPACELOWERBOUND/
%   DESIGNSPACEUPPERBOUND, and grows that input design box by the given growth
%   rate GROWTHRATE, returning the new box DESIGNBOXGROWN. Said growth takes the
%   the form, for each dimension: 
%       - NewLowerLimit = OldLowerLimit - GrowthRate * LengthOfDesignSpace
%       - NewUpperLimit = OldUpperLimit + GrowthRate * LengthOfDesignSpace
%   The limits of the new box are also limited by the design space.
%
%   Input:
%       - DESIGNBOX : (2,nDesignVariable) double
%       - DESIGNSPACELOWERBOUND : (1,nDesignVariable) double
%       - DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double
%       - GROWTHRATE : double
%   Output:
%       - DESIGNBOXGROWN : (2,nDesignVariable) double
%
%   See also sso_box_stochastic, candidatespaceboundingbox.
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

    % grow isotropically based on design space size
    designBoxGrown = designBox + growthRate*(designSpaceUpperBound-designSpaceLowerBound).*[-1;1];

    % Expand maximal to Design Space bounds
    designBoxGrown(1,:) = max(designBoxGrown(1,:), designSpaceLowerBound); % lower bound limit
    designBoxGrown(2,:) = min(designBoxGrown(2,:), designSpaceUpperBound); % upper bound limit
end