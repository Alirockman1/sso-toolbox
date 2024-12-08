function color = colormap_cmocean(varargin)
%COLORMAP_CMOCEAN Using cmocean colormaps in MATLAB (Color-blindness Friendly)
%	COLORMAP_CMOCEAN creates a MATLAB [R,G,B] triplet matrix for colormaps using  
%	one of the color schemes from 'cmocean'.
%	Link to project description: https://matplotlib.org/cmocean/
%	Link to repository with data used: https://github.com/matplotlib/cmocean
%	Citation:
%		Thyng, K.M., C.A. Greene, R.D. Hetland, H.M. Zimmerle, and S.F. DiMarco. 
%		2016. True colors of oceanography: Guidelines for effective and accurate 
%		colormap selection. Oceanography 29(3):9–13. 
%		http://dx.doi.org/10.5670/oceanog.2016.66
%
%	COLOR = COLORMAP_CMOCEAN returns the colormap as a three-column matrix of 
%	RGB triplets. 
%
%	COLOR = COLORMAP_CMOCEAN(MAP) allows one to choose one of the color maps
%	available. The options are:
%		- 'algae'
%		- 'amp'
%		- 'balance'
%		- 'curl'
%		- 'deep'
%		- 'delta'
%		- 'dense'
%		- 'diff'
%		- 'gray'
%		- 'haline'
%		- 'ice'
%		- 'matter'
%		- 'oxy'
%		- 'phase'
%		- 'rain'
%		- 'solar'
%		- 'speed'
%		- 'tarn'
%		- 'tempo'
%		- 'thermal'
%		- 'topo'
%		- 'turbid'
%	More about each choice should be seen on the official documentation:
%	https://matplotlib.org/cmocean/ . Default: 'algae'.
%
%	COLOR = COLORMAP_CMOCEAN(MAP,NLEVEL) allows one to choose how many levels
%	should be used. This is interpolated from the levels available for each map.
%	Default is empty (no interpolation, base map is used).
%
%	... = COLORMAP_CMOCEAN(...,NAME,VALUE,...) allows one to specify the 
%	following arguments:
%		- 'InterpolationOptions' : options for the interpolation done in case
%		a specific number of levels is given. Default is empty.
%		- 'Invert' : boolean flag to indicate whether to invert the color scale.
%
%	Input:
%		- MAP : char OR string
%		- NLEVEL : integer
%		- 'InterpolationOptions' : cell
%		- 'Invert' : logical
%		
%	Output:
%		- COLOR : (nColor,3)
%
%   See also colormap.
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
	parser.addOptional('Map','algae',@(x)ischar(x)||isstring(x));
	parser.addOptional('NumberOfLevels',[],@(x)isscalar(x)||isempty(x));
	parser.addParameter('InterpolationOptions',{});
	parser.addParameter('Invert',false,@(x)islogical(x));
	parser.parse(varargin{:});
	options = parser.Results;

	if(~options.Invert)
		chosenFormat = '-rgb.txt';
	else
		chosenFormat = '_i-rgb.txt';
	end
	filename = sprintf('cmocean-rgb/%s%s',lower(options.Map),chosenFormat);
	colormapBase = readmatrix(filename);
	nLevelBase = size(colormapBase,1);

	if(isempty(options.NumberOfLevels))
		color = colormapBase;
	else
		levelBase = 1:nLevelBase;
		levelDesired = linspace(1,nLevelBase,options.NumberOfLevels);
		color = interp1(levelBase,colormapBase,levelDesired,options.InterpolationOptions{:});
	end
end
