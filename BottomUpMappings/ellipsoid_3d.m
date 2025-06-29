function performanceMeasure = ellipsoid_3d(designSample,systemParameter)
%ELLIPSOID_3D Bottom-up Mapping (Designs in 3-Dimensional Ellipsoid)
%   ELLIPSOID_3D is a bottom-up mapping that returns the ellipsoid-norm of each  
%	design given the ellipsoid's quadric expression: (x-v)^T * A * (x-v), where
%	'x' is the tested design, 'v' is the center of the ellipse, and 'A' is the 
%	ellipsoid matrix which encodes its information (major semi-axes, rotations).
%	For a design to be inside the ellipsoid, the returned value must be less 
%	than or equal to 1.
%
%   PERFORMANCEMEASURE = ELLIPSOID_3D(DESIGNSAMPLE,SYSTEMPARAMETER) returns
%	for each DESIGNSAMPLE the result of the quadric expression in 
%	PERFORMANCEMEASURE, with SYSTEMPARAMETER containing: (1,2,3) center of the 
%	ellipse, (4,5,6) semi-major axes, and (7,8,9) rotation angles.
%
%   Input: 
%       - DESIGNSAMPLE : (nSample,2) double 
%       - SYSTEMPARAMETER : (1,9) double 
%           -- (1) : ellipse center in 1st dimension (vx)
%           -- (2) : ellipse center in 2nd dimension (vy)
%           -- (3) : ellipse center in 3rd dimension (vz)
%           -- (4) : semi-major axis in 1st dimension (a)
%           -- (5) : semi-major axis in 2nd dimension (b)
%           -- (6) : semi-major axis in 3rd dimension (c)
%           -- [optional] (7) : rotation angle around center about 1st dimension 
%			(Euler-angle alpha)
%           -- [optional] (8) : rotation angle around center about 2nd dimension 
%			(Euler-angle beta)
%           -- [optional] (9) : rotation angle around center about 3rd dimension 
%			(Euler-angle gamma)
%
%   Output:
%       - PERFORMANCEMEASURE : (nSample,1) double 
%
%   See also ellipse_2d, distance_to_center, sphere_nd.
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

	ellipsoidCenter = systemParameter(1:3);
	ellipsoidMajorSemiAxis = systemParameter(4:6);

	if(length(systemParameter)<7)
		ellipsoidRotationAngle = [0 0 0];
	else
		ellipsoidRotationAngle = systemParameter(7:9);
	end	

	% compute rotation matrix from the matrices of each rotation
	rotationMatrixX = [...
		1 0 0;...
		0 cos(ellipsoidRotationAngle(1)),-sin(ellipsoidRotationAngle(1));...
		0 sin(ellipsoidRotationAngle(1)),cos(ellipsoidRotationAngle(1))];
	rotationMatrixY = [...
		cos(ellipsoidRotationAngle(2)),0,sin(ellipsoidRotationAngle(2));...
		0 1 0;...
		-sin(ellipsoidRotationAngle(2)),0,cos(ellipsoidRotationAngle(2))];
	rotationMatrixZ = [...
		cos(ellipsoidRotationAngle(3)),-sin(ellipsoidRotationAngle(3)),0;...
		sin(ellipsoidRotationAngle(3)),cos(ellipsoidRotationAngle(3)),0;...
		0 0 1];
	rotationMatrixTotal = rotationMatrixZ * rotationMatrixY * rotationMatrixX;

	% compute matrix of quadric form of general ellipsoid
	semiAxesMatrix = diag(1./(ellipsoidMajorSemiAxis).^2);
	ellipsoidMatrix = rotationMatrixTotal * semiAxesMatrix * rotationMatrixTotal';

	nSample = size(designSample,1);
	performanceMeasure = nan(nSample,1);
	for i=1:nSample
		distanceToEllipsoidCenter = (designSample(i,:) - ellipsoidCenter)';
		performanceMeasure(i) = distanceToEllipsoidCenter' * ellipsoidMatrix * distanceToEllipsoidCenter;
	end
end