function performanceMeasure = ellipse_2d(designSample,systemParameter)
%ELLIPSE_2D Bottom-up Mapping (Designs in 2-Dimensional Ellipse)
%   ELLIPSE_2D is a bottom-up mapping that returns the ellipse-norm of each
%	design given the ellipse's quadric expression: (x-v)^T * A * (x-v), where 
%	'x' is the tested design, 'v' is the center of the ellipse, and 'A' is the 
%	ellipse matrix which encodes its information (major semi-axes, rotation).
%	For a design to be inside the ellipse, the returned value must be less 
%	than or equal to 1.
%
%   PERFORMANCEMEASURE = ELLIPSE_2D(DESIGNSAMPLE,SYSTEMPARAMETER) returns
%	for each DESIGNSAMPLE the result of the quadric expression in 
%	PERFORMANCEMEASURE, with SYSTEMPARAMETER containing: (1,2) center of the 
%	ellipse, (3,4) semi-major axes, and (5) rotation angle.
%
%   Input: 
%       - DESIGNSAMPLE : (nSample,2) double 
%       - SYSTEMPARAMETER : (1,5) double 
%           -- (1) : ellipse center in 1st dimension (vx)
%           -- (2) : ellipse center in 2nd dimension (vy)
%           -- (3) : semi-major axis in 1st dimension (a)
%           -- (4) : semi-major axis in 2nd dimension (b)
%           -- [optional] (5) : rotation angle around center (theta)
%
%   Output:
%       - PERFORMANCEMEASURE : (nSample,1) double 
%
%   See also ellipsoid_3d, distance_to_center, sphere_nd.
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

	ellipseCenter = systemParameter(1:2);
	ellipseMajorSemiAxis = systemParameter(3:4);

	if(length(systemParameter)<5)
		ellipseRotationAngle = 0;
	else
		ellipseRotationAngle = systemParameter(5);
	end	

	% compute rotation matrix assuming an angle in direction z+
	rotationMatrix = [...
		cos(ellipseRotationAngle),-sin(ellipseRotationAngle);...
		sin(ellipseRotationAngle),cos(ellipseRotationAngle)];

	% compute matrix of quadric form of general ellipse
	semiAxesMatrix = diag(1./(ellipseMajorSemiAxis).^2);
	ellipseMatrix = rotationMatrix * semiAxesMatrix * rotationMatrix';

	nSample = size(designSample,1);
	performanceMeasure = nan(nSample,1);
	for i=1:nSample
		distanceToEllipseCenter = (designSample(i,:) - ellipseCenter)';
		performanceMeasure(i) = distanceToEllipseCenter' * ellipseMatrix * distanceToEllipseCenter;
	end
end