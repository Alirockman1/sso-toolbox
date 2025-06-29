function color = color_palette_okabe_ito(index)
%COLOR_PALETTE_OKABE_ITO Colorblind-friendly color palette from Okabe, Ito
%	COLOR_PALETTE_OKABE_ITO returns a color palette developed by Masataka Okabe
%	and Kei Ito which is more accessible to colorblind people.
%	Source of palette: https://jfly.uni-koeln.de/color/
%	<a href="https://davidmathlogic.com/colorblind/#%23000000-%23E69F00-%2356B4E9-%23009E73-%23F0E442-%230072B2-%23D55E00-%23CC79A7">Testing the colorblind appearance</a>
%
%	COLOR = COLOR_PALETTE_OKABE_ITO returns an array with all the colors in the 
%	palette. These are expressed as [R,G,B] triplets with values from 0 to 1,
%	with each row being a new color.
%
%	COLOR = COLOR_PALETTE_OKABE_ITO(INDEX) allows the choice of colors specified
%	by their index INDEX. Can be left empty to choose all colors. Default is
%	empty. Colors can also be specified by name:
%		- (1) 'black' : #000000
%		- (2) 'orange' : #E69F00
%		- (3) 'sky-blue' : #56B4E9
%		- (4) 'bluish-green' : #009E73
%		- (5) 'yellow' : #F0E442
%		- (6) 'blue' : #0072B2
%		- (7) 'vermillion' : #D55E00
%		- (8) 'reddish-purple' : #CC79A7
%
%	Input:
%		- INDEX : (nChoice) integer OR (nChoice) cell
%		
%	Output:
%		- COLOR : (8,3) OR (nChoice,3) double 
%
%   See also color_palette_ibm, color_palette_tol, orderedcolors, 
%	colororder.
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

	colorPalette = [...
		  0   0   0;... % black #000000
		230 159   0;... % orange #E69F00
		 86 180 233;... % sky blue #56B4E9
		  0 158 115;... % bluish green #009E73
		240 228  66;... % yellow #F0E442
		  0 114 178;... % blue #0072B2
		213  94   0;... % vermillion #D55E00
		204 121 167;... % reddish purple #CC79A7
		]./255;
	colorName = {...
		'black';...
		'orange';...
		'sky-blue';...
		'bluish-green';...
		'yellow';...
		'blue';...
		'vermillion';...
		'reddish-purple';...
		};

	if(nargin<1 || isempty(index))
		index = 1:size(colorPalette,1);
	end

	if(ischar(index)||isstring(index))
		index = {index};
	end
	if(iscell(index))
		nEntry = length(index);
		numericalIndex = nan(nEntry,1);

		for i=1:nEntry
			choice = strcmpi(colorName,index{i});

			if(~any(choice))
				error('ColorPaletteOkabeIto:NameNotFound',['Color ',index{i},...
					' not found for palette Okabe-Ito.']);
			end

			numericalIndex(i) = find(choice);
        end

        index = numericalIndex;
	end

	color = colorPalette(index,:);
end