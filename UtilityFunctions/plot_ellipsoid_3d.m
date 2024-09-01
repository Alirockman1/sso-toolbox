function handleSurf = plot_ellipsoid_3d(figureHandle,ellipsoidCenter,ellipsoidMajorSemiAxis,varargin)
%PLOT_ELLIPSOID_3D Visualize ellipsoid (or sphere) in figure
%	PLOT_ELLIPSOID_3D plots an ellipsoid in the given figure given the 
%	ellipsoid's center, major semi-axes, and rotation angle. If only one major 
%	semi-axis length is given, it is assumed the figure is a sphere instead.
%
%	PLOT_ELLIPSOID_3D(FIGUREHANDLE,ELLIPSOIDCENTER,ELLIPSOIDMAJORSEMIAXIS) plots
%	in figure FIGUREHANDLE an ellipsoid with center on ELLIPSOIDCENTER and with  
%	major semi-axes ELLIPSOIDMAJORSEMIAXIS. If ELLIPSOIDMAJORSEMIAXIS is a 
%	single value, a sphere is plotted instead, with that value being the assumed
%	radius.
%
%	PLOT_ELLIPSOID_3D(FIGUREHANDLE,ELLIPSOIDCENTER,ELLIPSOIDMAJORSEMIAXIS,
%	ROTATIONANGLE) allows one to rotate the ellipsoid by the given angles in 
%	radians; (orientation: [+X,+Y,+Z]). When applying the rotation matrices,
%	the order of application is first rotation in X, then Y, and finally Z.
%
%	PLOT_ELLIPSOID_3D(...NAME,VALUE,...) also allows for the choice of 
%	additional options. These are:
%		- 'NumberPoints' : number of points used for each angle parameter 
%		(polar, azimuth), so the ellipsoid will have a total of points which
%		is equal to the square of this. Default: 10.
%		- 'SurfOptions' : options to be used in the 'surf' function. Default
%		is empty.
%
%   Input:
%		- FIGUREHANDLE : Figure
%		- ELLIPSOIDCENTER : (1,3) double
%		- ELLIPSOIDMAJORSEMIAXIS : (1,3) double OR double
%		- ROTATIONANGLE : (1,3) double
%		- 'NumberPoints' : double
%		- 'SurfOptions' : (1,nOption) cell
%
%   Output:
%		- HANDLEPATCH : Line
%
%   See also surf, meshgrid.
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

	parser = inputParser;
	parser.addOptional('RotationAngle',[0 0 0]);
	parser.addParameter('NumberPoints',10);
	parser.addParameter('SurfOptions',{});
	parser.parse(varargin{:});
	options = parser.Results;

	% if only one semi-axis is given - assume sphere
	if(length(ellipsoidMajorSemiAxis)==1)
		ellipsoidMajorSemiAxis([2,3]) = ellipsoidMajorSemiAxis(1);
	end

	% generate points based on parametrized description
	polarAngleInterval = linspace(0,pi,options.NumberPoints); % theta
	azimuthAngleInterval = linspace(0,2*pi,options.NumberPoints); % phi
	[polarAngleGrid, azimuthAngleGrid] = meshgrid(polarAngleInterval,azimuthAngleInterval);
	polarAnglePoint = polarAngleGrid(:);
	azimuthAnglePoint = azimuthAngleGrid(:);
	xPointLocal = ellipsoidMajorSemiAxis(1).*sin(polarAnglePoint).*cos(azimuthAnglePoint);
	yPointLocal = ellipsoidMajorSemiAxis(2).*sin(polarAnglePoint).*sin(azimuthAnglePoint);
	zPointLocal = ellipsoidMajorSemiAxis(3).*cos(polarAnglePoint);

	% generate rotation matrix
	rotationMatrixX = [...
		1 0 0;...
		0 cos(options.RotationAngle(1)),-sin(options.RotationAngle(1));...
		0 sin(options.RotationAngle(1)),cos(options.RotationAngle(1))];
	rotationMatrixY = [...
		cos(options.RotationAngle(2)),0,sin(options.RotationAngle(2));...
		0 1 0;...
		-sin(options.RotationAngle(2)),0,cos(options.RotationAngle(2))];
	rotationMatrixZ = [...
		cos(options.RotationAngle(3)),-sin(options.RotationAngle(3)),0;...
		sin(options.RotationAngle(3)),cos(options.RotationAngle(3)),0;...
		0 0 1];
	rotationMatrixTotal = rotationMatrixZ * rotationMatrixY * rotationMatrixX;

	% translate/rotate local points
	xPoint = ellipsoidCenter(1) + ...
		xPointLocal*rotationMatrixTotal(1,1) + yPointLocal*rotationMatrixTotal(1,2) + zPointLocal*rotationMatrixTotal(1,3);
	yPoint = ellipsoidCenter(2) + ...
		xPointLocal*rotationMatrixTotal(2,1) + yPointLocal*rotationMatrixTotal(2,2) + zPointLocal*rotationMatrixTotal(2,3);
	zPoint = ellipsoidCenter(3) + ...
		xPointLocal*rotationMatrixTotal(3,1) + yPointLocal*rotationMatrixTotal(3,2) + zPointLocal*rotationMatrixTotal(3,3);

	xPointGrid = reshape(xPoint, size(polarAngleGrid));
	yPointGrid = reshape(yPoint, size(polarAngleGrid));
	zPointGrid = reshape(zPoint, size(polarAngleGrid));

	figure(figureHandle);
	handleSurf = surf(xPointGrid,yPointGrid,zPointGrid,options.SurfOptions{:});

	if(nargout<1)
        clear handleSurf
    end
end