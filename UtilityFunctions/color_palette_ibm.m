function color = color_palette_ibm(index)
%COLOR_PALETTE_IBM Color Blind Accessible color palette from IBM Design Library
%	COLOR_PALETTE_IBM returns a color palette developed by IBM which is supposed
%	which is more accessible to colorblind people.
%	<a href="https://davidmathlogic.com/colorblind/#%23648FFF-%23785EF0-%23DC267F-%23FE6100-%23FFB000">Testing the colorblind appearance</a>
%
%	COLOR = COLOR_PALETTE_IBM returns an array with all the colors in the 
%	palette. These are expressed as [R,G,B] triplets with values from 0 to 1,
%	with each row being a new color.
%
%	COLOR = COLOR_PALETTE_IBM(INDEX) allows the choice of colors specified by
%	their index INDEX. Can be left empty to choose all colors. Default is
%	empty. Colors can also be specified by name:
%		- (1) 'blue' : #648fff
%		- (2) 'purple' : #785EF0
%		- (3) 'pink' : #DC267F
%		- (4) 'orange' : #FE6100
%		- (5) 'yellow' : #FFB000
%
%	Input:
%		- INDEX : (nChoice) integer OR (nChoice) cell
%		
%	Output:
%		- COLOR : (5,3) OR (nChoice,3) double 
%
%   See also color_palette_okabe_ito, color_palette_tol, orderedcolors, 
%	colororder.
%
%   Copyright 2025 Eduardo Rodrigues Della Noce
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
		100 143 255;... % blue #648fff
		120  94 240;... % purple #785EF0
		220  38 127;... % pink #DC267F
		254  97   0;... % orange #FE6100
		255 176   0;... % yellow #FFB000
		]./255;
	colorName = {...
		'blue';...
		'purple';...
		'pink';...
		'orange';...
		'yellow';...
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
				error('ColorPaletteIbm:NameNotFound',['Color ',index{i},...
					' not found for palette IBM.']);
			end

			numericalIndex(i) = find(choice);
        end

        index = numericalIndex;
	end

	color = colorPalette(index,:);
end