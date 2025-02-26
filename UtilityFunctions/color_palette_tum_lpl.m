function color = color_palette_tum_lpl(index)
%COLOR_PALETTE_TUM_LPL Color palette from TUM LPL
%	COLOR_PALETTE_TUM_LPL returns a color palette based on the Technical 
%	University of Munich (TUM) corporate design colors.
%
%	COLOR = COLOR_PALETTE_TUM_LPL returns an array with all the colors in the 
%	palette. These are expressed as [R,G,B] triplets with values from 0 to 1,
%	with each row being a new color.
%
%	COLOR = COLOR_PALETTE_TUM_LPL(INDEX) allows the choice of colors specified
%	by their index INDEX. Can be left empty to choose all colors. Default is
%	empty. Colors can also be specified by name:
%		- (1) 'TUM-blue' : #0065BD
%		- (2) 'black' : #000000
%		- (3) 'white' : #FFFFFF
%		- (4) 'TUM-light-blue' : #005293
%		- (5) 'TUM-dark-blue' : #003359
%		- (6) 'blue' : #64A0C8
%		- (7) 'light-blue' : #98C6EA
%		- (8) 'green' : #A2AD00
%		- (9) 'vermillion' : #E37222
%		- (10) 'grey' : #DAD7CB
%		- (11) 'TUM-light-blue-80' : #3375A9
%		- (12) 'TUM-light-blue-50' : #7FA8C9
%		- (13) 'TUM-light-blue-20' : #CCDCE9
%		- (14) 'TUM-dark-blue-80' : #335C7A
%		- (15) 'TUM-dark-blue-50' : #7F99AC
%		- (16) 'TUM-dark-blue-20' : #CCD6DE
%		- (17) 'grey-80' : #404040
%		- (18) 'grey-50' : #7F7F7F
%		- (19) 'grey-20' : #D9D9D9
%
%	Input:
%		- INDEX : (nChoice) integer OR (nChoice) cell
%		
%	Output:
%		- COLOR : (19,3) OR (nChoice,3) double 
%
%   See also color_palette_ibm, color_palette_okabe_ito, color_palette_tol, 
%	orderedcolors, colororder.
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
        0   101 189;... % primary color: TUM blue #0065BD
        0     0   0;... % primary color: black #000000
        255   255 255;... % primary color: white #FFFFFF
        0    82 147;... % secondary color: TUM light blue #005293
        0    51  89;... % secondary color: TUM dark blue #003359
        100   160 200;... % accent color: blue #64A0C8
        152   198 234;... % accent color: light blue #98C6EA
        162   173   0;... % accent color: green #A2AD00
        227   114  34;... % accent color: vermillion #E37222
        218   215 203;... % accent color: grey #DAD7CB
        51   117 169;... % TUM-light-blue variant at 80% #3375A9
        127   168 201;... % TUM-light-blue variant at 50% #7FA8C9
        204   220 233;... % TUM-light-blue variant at 20% #CCDCE9
        51    92 122;... % TUM-dark-blue variant at 80% #335C7A
        127   153 172;... % TUM-dark-blue variant at 50% #7F99AC
        204   214 222;... % TUM-dark-blue variant at 20% #CCD6DE
        64    64  64;... % grey variant: grey-80 #404040
        127   127 127;... % grey variant: grey-50 #7F7F7F
        217   217 217;... % grey variant: grey-20 #D9D9D9
    ]./255;
    colorName = {...
        'TUM-blue';...
        'black';...
        'white';...
        'TUM-light-blue';...
        'TUM-dark-blue';...
        'blue';...
        'light-blue';...
        'green';...
        'vermillion';...
        'grey';...
        'TUM-light-blue-80';...
        'TUM-light-blue-50';...
        'TUM-light-blue-20';...
        'TUM-dark-blue-80';...
        'TUM-dark-blue-50';...
        'TUM-dark-blue-20';...
        'grey-80';...
        'grey-50';...
        'grey-20';...
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
                error('ColorPaletteTumLpl:NameNotFound',['Color ',index{i},...
                    ' not found for palette TUM-LPL.']);
            end

            numericalIndex(i) = find(choice);
        end

        index = numericalIndex;
    end

    color = colorPalette(index,:);
end
