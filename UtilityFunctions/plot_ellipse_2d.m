function handlePatch = plot_ellipse_2d(figureHandle,ellipseCenter,ellipseMajorSemiAxis,varargin)
%PLOT_ELLIPSE_2D Visualize ellipse (or circle) in figure
%	PLOT_ELLIPSE_2D plots an ellipse in the given figure given the ellipse's
%	center, major semi-axes, and rotation angle. If only one major semi-axis 
%	length is given, it is assumed the figure is a circle instead.
%
%	PLOT_ELLIPSE_2D(FIGUREHANDLE,ELLIPSECENTER,ELLIPSEMAJORSEMIAXIS) plots in 
%	figure FIGUREHANDLE an ellipse with center on ELLIPSECENTER and with major 
%	semi-axes ELLIPSEMAJORSEMIAXIS. If ELLIPSEMAJORSEMIAXIS is a single value,
%	a circle is plotted instead, with that value being the assumed radius.
%
%	PLOT_ELLIPSE_2D(FIGUREHANDLE,ELLIPSECENTER,ELLIPSEMAJORSEMIAXIS,
%	ROTATIONANGLE) allows one to rotate the ellipse by the given angle
%	in radians (orientation: +Z).
%
%	PLOT_ELLIPSE_2D(...NAME,VALUE,...) also allows for the choice of additional
%	options. These are:
%		- 'NumberPoints' : number of points used to build the ellipse. 
%		Default: 100.
%		- 'PatchOptions' : options to be used in the 'patch' function. Default
%		is empty.
%
%   Input:
%		- FIGUREHANDLE : Figure
%		- ELLIPSECENTER : (1,2) double
%		- ELLIPSEMAJORSEMIAXIS : (1,2) double OR double
%		- ROTATIONANGLE : double
%		- 'NumberPoints' : double
%		- 'PatchOptions' : (1,nOption) cell
%
%   Output:
%		- HANDLEPATCH : Line
%
%   See also patch.
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
	parser.addOptional('RotationAngle',0);
	parser.addParameter('NumberPoints',100);
	parser.addParameter('PatchOptions',{});
	parser.parse(varargin{:});
	options = parser.Results;

	% if only one semi-axis is given - assume circle
	if(length(ellipseMajorSemiAxis)==1)
		ellipseMajorSemiAxis(2) = ellipseMajorSemiAxis(1);
	end

	% generate points based on parametrized description
	anglePoint = linspace(0,2*pi,options.NumberPoints);
	xPointLocal = ellipseMajorSemiAxis(1).*cos(anglePoint);
	yPointLocal = ellipseMajorSemiAxis(2).*sin(anglePoint);

	% translate/rotate local points
	xPoint = ellipseCenter(1) + ...
		xPointLocal*cos(options.RotationAngle) - yPointLocal*sin(options.RotationAngle);
	yPoint = ellipseCenter(2) + ...
		xPointLocal*sin(options.RotationAngle) + yPointLocal*cos(options.RotationAngle);

	figure(figureHandle);
	handlePatch = patch('XData',xPoint,'YData',yPoint,options.PatchOptions{:});

	if(nargout<1)
        clear handlePatch
    end
end