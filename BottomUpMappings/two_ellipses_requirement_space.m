function [performanceMeasure, physicalFeasibilityMeasure] = two_ellipses_requirement_space(designSample,systemParameter)
%TWO_ELLIPSES_REQUIREMENT_SPACE Bottom-Up Mapping (Requirement Spaces Testing)
%   TWO_ELLIPSES_REQUIREMENT_SPACE is a bottom-up mapping characterized by two 
%   ellipsoids:
%       - Designs inside the first ellipse meet the requirements of the problem
%       - Designs inside the second ellipse are physically feasible
%   With this, one can test any methods they wish for requirement spaces.
%   Designs are inside each sphere if the returned measure has a value less than
%   or equal to 1.
%
%   The characteristics of the ellipses can be set using system parameters. 
%   However, they also have default values for quick analysis.
%   Unless specified, the first ellipsoid (performance) has the following 
%   characteristics:
%       - Center: [v1,v2] = [11,5]
%       - Principal Semiaxes: [a,b] = [7,2]
%       - Rotation Angle: theta = 1.8 rad
%   Unless specified, the second ellipsoid (physical feasibility) has the  
%   following characteristics:
%       - Center: [v1,v2] = [5,6]
%       - Principal Semiaxes: [a,b] = [7,2]
%       - Rotation Angle: theta = 0.2 rad
%
%   PERFORMANCEMEASURE = TWO_ELLIPSES_REQUIREMENT_SPACE(DESIGNSAMPLE) receives 
%   the spatial coordinates of each design sample in DESIGNSAMPLE and returns 
%   the ellipse-norm regarding the performance of each in PERFORMANCEMEASURE.
%
%   PERFORMANCEMEASURE = TWO_ELLIPSES_REQUIREMENT_SPACE(DESIGNSAMPLE,
%   SYSTEMPARAMETER) allows one to specify the characteristics of each ellipse.
%   The first row corresponds to the performance ellipse, and the second row 
%   to the physical feasibility ellipse.
%
%   [PERFORMANCEMEASURE,PHYSICALFEASIBILITYMEASURE] = 
%   TWO_ELLIPSES_REQUIREMENT_SPACE(...)  
%   also returns the ellipse-norm of each design regarding their physical 
%   feasibility in PHYSICALFEASIBILITYMEASURE;
%
%   Input: 
%       - DESIGNSAMPLE : (nSample,2) double
%           -- (1) : spatial coordinate 1
%           -- (2) : spatial coordinate 2
%       - SYSTEMPARAMETER [optional] : (2,5) double
%           -- (1) : ellipse center in 1st dimension (vx)
%           -- (2) : ellipse center in 2nd dimension (vy)
%           -- (3) : semi-major axis in 1st dimension (a)
%           -- (4) : semi-major axis in 2nd dimension (b)
%           -- [optional] (5) : rotation angle around center (theta)
%
%   Output:
%       - PERFORMANCEMEASURE : (nSample,1) double
%       - PHYSICALFEASIBILITYMEASURE : (nSample,1) double
%
%   See also ellipse_2d.
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

    if(nargin<2 || isempty(systemParameter))
        systemParameter = [...
            11, 5, 7, 2, 1.8; ... % performance: vx, vy, a, b, theta
            5, 6, 7, 2, 0.2]; % physical feasibility: vx, vy, a, b, theta
        end
    
    performanceMeasure = ellipse_2d(designSample,systemParameter(1,:));
    physicalFeasibilityMeasure = ellipse_2d(designSample,systemParameter(2,:));
end

