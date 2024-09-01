function [isInside,score] = is_in_design_box(designSample,designBox)
%IS_IN_DESIGN_BOX Test if design sample points are inside a design box
%   IS_IN_DESIGN_BOX tests if the given design sample points are inside a design
%   box (lower and upper boundaries for each design variable).
%
%   ISINSIDE = IS_IN_DESIGN_BOX(DESIGNSAMPLE,DESIGNBOX) checks if the design 
%   sample points DESIGNSAMPLE are inside the design box DESIGNBOX, checking for  
%   each sample if it is inside ('true') or outside ('false') of said box, 
%   aggregating results in ISINSIDE.
%
%   [ISINSIDE,SCORE] = IS_IN_DESIGN_BOX(...) also returns a score for each 
%   sample point; said score can be computed as follows:
%       - lower limit score: (box lower limit - design variable value) 
%                            / interval size
%       - upper limit score: (design variable value - box upper limit)
%                            / interval size
%       - score : maximum between lower limit and upper limit score
%
%   Input:
%       - DESIGNSAMPLE : (nSample,nDesignVariable)
%       - DESIGNBOX : (2,nDesignVariable)
%           -- (1) : lower boundary of the design box
%           -- (2) : upper boundary of the design box
%
%   Output:
%       - INSIDE : (nSample,1) logical
%       - SCORE : (nSample,1) double
%
%   See also is_in_convex_hull_with_face.
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

    % samples are inside box if all their individual entries are within 
    % the defined lower and upper boundaries of the box
    isInside = all(designSample>=designBox(1,:),2) & all(designSample<=designBox(2,:),2);

    intervalFactor = designBox(2,:) - designBox(1,:);
    score = max([designBox(1,:)-designSample,designSample-designBox(2,:)]./[intervalFactor,intervalFactor],[],2);
end