function color = colormap_viridis(varargin)
%COLORMAP_VIRIDIS Using viridis colormaps in MATLAB (Color-blindness Friendly)
%	COLORMAP_VIRIDIS creates a MATLAB [R,G,B] triplet matrix for colormaps using  
%	one of the color schemes from color scheme from 'viridis'.
%	Link to repository with data used: https://github.com/sjmgarnier/viridisLite
%
%	COLOR = COLORMAP_VIRIDIS returns the colormap as a three-column matrix of 
%	RGB triplets. 
%
%	COLOR = COLORMAP_VIRIDIS(MAP) allows one to choose one of the color maps
%	available. The options are:
%		- 'magma' (Option A)
%		- 'inferno' (Option B)
%		- 'plasma' (Option C)
%		- 'viridis' (Option D)
%		- 'cividis' (Option E)
%		- 'rocket' (Option F)
%		- 'mako' (Option G)
%		- 'turbo' (Option H)
%	More about each choice should be seen on the official documentation:
%	https://github.com/sjmgarnier/viridisLite/ . Default: 'viridis'.
%
%	COLOR = COLORMAP_VIRIDIS(MAP,NLEVEL) allows one to choose how many levels
%	should be used. This is interpolated from the levels available for each map.
%	Default is empty (no interpolation, base map is used).
%
%	... = COLORMAP_VIRIDIS(...,NAME,VALUE,...) allows one to specify the 
%	following arguments:
%		- 'InterpolationOptions' : options for the interpolation done in case
%		a specific number of levels is given. Default is empty.
%
%	Input:
%		- MAP : char OR string
%		- NLEVEL : integer
%		- 'InterpolationOptions' : cell
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
	parser.addOptional('Map','viridis',@(x)ischar(x)||isstring(x));
	parser.addOptional('NumberOfLevels',[],@(x)isscalar(x)||isempty(x));
	parser.addParameter('InterpolationOptions',{});
	parser.parse(varargin{:});
	options = parser.Results;

	filename = sprintf('viridis-rgb/%s-rgb.csv',lower(options.Map));
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
