function performanceMeasure = sphere_nd(designSample)
%SPHERE_ND Bottom-Up Mapping (n-Dimensional Benchmark Test)
%   SPHERE_ND is a bottom-up mapping that returns the radial distance of a 
%   n-dimensional design to a sphere of a particular radius; the radius is 
%   constructed such that the optimal box-shaped solution space for this problem
%   is fixed and known. For a design to be inside the sphere, the returned
%   value must be lesser than or equal to 0.
%   SPHERE_ND can be used to benchmark the performance/accuracy of an algorithm 
%   to compute optimal solution boxes. With the way this is constructed, the 
%   analytical solution for each design variable is always the interval 
%   [-0.5,0.5].
% 
%   Work with a solution space that is a n-dimensional sphere with center 
%   at [0,0,0,...,0]; the biggest box that fits in a sphere is (diagonal):
%       - 1-dimensional: a = (1/2) * R
%       - 2-dimensional: a^2 + a^2 = (2R)^2 
%           -> a = sqrt(2) * R
%       - 3-dimensional: a^2 + a^2 + a^2 = (2R)^2 
%           -> a = sqrt(4/3) * R
%       - In general: n*a^2 = (2R)^2 
%           -> a = 2/sqrt(n) * R
%   Choosing a radius of solution space sphere R = sqrt(n)/2 will
%   always leads to a maximum solution box of side length a=1: 
%   [-0.5, 0.5]^n.
%
%   PERFORMANCEMEASURE = SPHERE_ND(DESIGNSAMPLE) receives the
%   the spatial coordinates [x1,x2,x3,...] of the design samples in DESIGNSAMPLE
%   and returns the radial distance to the shell of the n-dimensional sphere
%   in PERFORMANCEMEASURE. For a design to be inside the sphere, the 
%   PERFORMANCEMEASURE value should be lesser than or equal to zero.
%
%   Example: 
%       For a three-dimensional problem:
%           designSample = [0.5,0.3,0.1 ; 2,3,1]; % x1,x2,x3
%           performanceMeasure = sphere_nd(designSample);
%       PERFORMANCEMEASURE will return the distance of these points to the shell
%       of the sphere of radius sqrt(3)/2; negative numbers means it is inside
%       said sphere. In this example, the first design is inside the sphere,
%       but the second one is not.
%
%   Input: 
%       - DESIGNSAMPLE : (nSample,nDimension) double
%
%   Output:
%       - PERFORMANCEMEASURE : (nSample,1) double
%
%   See also distance_to_center, ellipse_2d, ellipsoid_3d.
%
%   First version by Julian Stumpf.
%   Re-written and documented by Eduardo Rodrigues Della Noce.
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

    % dimension of design space d
    nDimension = size(designSample,2);

    % Check which designs are inside sphere of radius sqrt(d)/2
    radius = sqrt(nDimension)/2;
    distanceRadial = vecnorm(designSample,2,2);
    performanceMeasure = distanceRadial - radius;
end