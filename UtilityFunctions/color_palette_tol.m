function color = color_palette_tol(index,paletteType)
%COLOR_PALETTE_TOL Colorblind-friendly color palette from Dr. Paul Tol
%	COLOR_PALETTE_TOL returns a color palette developed by Dr. Paul Tol which is 
%	more accessible to colorblind people. Multiple palettes are available, and 
%	the choice of which one to get can be specified.
%	Source of palettes: https://personal.sron.nl/~pault/
%	Testing the colorblind appearance: https://davidmathlogic.com/colorblind/
%
%	COLOR = COLOR_PALETTE_TOL returns an array with all the colors in the 
%	palette. These are expressed as [R,G,B] triplets with values from 0 to 1,
%	with each row being a new color.
%
%	COLOR = COLOR_PALETTE_TOL(INDEX) allows the choice of colors specified
%	by their index INDEX. Can be left empty to choose all colors. Default is
%	empty.
%
%	COLOR = COLOR_PALETTE_TOL(INDEX,PALETTETYPE) also allows the 
%	specfication of which color palette is desired PALETTETYPE. The following 
%	options are available: 
%		- 'bright' : 7 colors, discrete use.
%		- 'high-contrast' : 5 colors, discrete use.
%		- 'vibrant' : 7 colors, discrete use.
%		- 'muted' : 10 colors, discrete use, last color for bad data.
%		- 'medium-contrast' : 8 colors, discrete use.
%		- 'pale' : 6 colors, discrete use.
%		- 'dark' : 6 colors, discrete use.
%		- 'light' : 9 colors, discrete use.
%		- 'sunset' : 12 colors, can be interpolated, last color for bad data.
%		- 'nightfall' : 18 colors, can be interpolated, last color for bad data.
%		- 'burd' : 10 colors, can be interpolated, last color for bad data.
%		- 'prgn' : 10 colors, can be interpolated, last color for bad data.
%		- 'ylorbr' : 10 colors, can be interpolated, last color for bad data.
%		- 'iridescent' : 24 colors, can be interpolated, last color for bad 
%		data.
%		- 'incandescent' : 12 colors, can be interpolated, last color for bad 
%		data.
%		- 'discrete-rainbow' : 15 colors, can be interpolated, last color for 
%		bad data.
%		- 'discrete-rainbow-23' : 24 colors, can be interpolated, last color for 
%		bad data.
%		- 'smooth-rainbow' : 35 colors, can be interpolated, last color for bad 
%		data.
%	More about each palette can be found on https://personal.sron.nl/~pault/ .
%	Default choice is 'vibrant'.
%
%	Input:
%		- INDEX : (nChoice) integer
%		- PALETTETYPE : string
%		
%	Output:
%		- COLOR : (nColor,3) OR (nChoice,3) double 
%
%   See also color_palette_ibm, color_palette_okabe_ito, orderedcolors, 
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


	if(nargin<2 || isempty(paletteType))
		paletteType = 'vibrant';
	end

	switch lower(paletteType)
		case 'bright'
			colorPalette = [... % default order
				 68 119 170;... % blue #4477AA
				238 102 119;... % red #EE6677
				 34 136  51;... % green #228833
				204 187  68;... % yellow #CCBB44
				102 204 238;... % cyan #66CCEE
				170  51 119;... % purple #AA3377
				187 187 187;... % grey #BBBBBB
				]./255;
		case 'high-contrast'
			colorPalette = [... % default order
				255 255 255;... % white #FFFFFF
				  0  68 136;... % blue #004488
				221 170  51;... % yellow #DDAA33
				187  85 102;... % red #BB5566
				  0   0   0;... % black #000000
				]./255;
		case 'vibrant'
			colorPalette = [... % default order
				238 119  51;... % orange #EE7733
				  0 119 187;... % blue #0077BB
				 51 187 238;... % cyan #33BBEE
				 238  51 119;... % magenta #EE3377
				204  51  17;... % red #CC3311
				  0 153 136;... % teal #009988
				187 187 187;... % grey #BBBBBB
				]./255;
		case 'muted'
			colorPalette = [... % default order
				204 102 119;... % rose #CC6677
				 51  34 136;... % indigo #332288
				221 204 119;... % sand #DDCC77
				 17 119  51;... % green #117733
				136 204 238;... % cyan #88CCEE
				136  34  85;... % wine #882255
				170  68 153;... % purple #44AA99
				153 153  51;... % olive #999933
				 68 170 153;... % teal #AA4499
				221 221 221;... % pale grey #DDDDDD - bad data
				]./255;
		case 'medium-contrast'
			colorPalette = [... % default order
				255 255 255;... % white #FFFFFF
				102 153 204;... % light blue #6699CC
				  0  68 136;... % dark blue #004488
				238 204 102;... % light yellow #EECC66
				153  68  85;... % dark red #994455
				153 119   0;... % dark yellow #997700
				238 153 170;... % light red #EE99AA
				  0   0   0;... % black #000000
				]./255;
		case 'pale'
			colorPalette = [...
				187 204 238;... % pale blue #BBCCEE
				204 238 255;... % pale cyan #CCEEFF
				204 221 170;... % pale green #CCDDAA
				238 238 187;... % pale yellow #EEEEBB
				255 204 204;... % pale red #FFCCCC
				221 221 221;... % pale grey #DDDDDD
				]./255;
		case 'dark'
			colorPalette = [...
				 34  34  85;... % dark blue #222255
				 34  85  85;... % dark cyan #225555
				 34  85  34;... % dark green #225522
				102 102  51;... % dark yellow #666633
				102  51  51;... % dark red #663333
				 85  85  85;... % dark grey #555555
				]./255;
		case 'light'
			colorPalette = [... % default order
				119 170 221;... % light blue #77AADD
				238 136 102;... % orange #EE8866
				238 221 136;... % light yellow #EEDD88
				255 170 187;... % pink #FFAABB
				153 221 255;... % light cyan #99DDFF
				 68 187 153;... % mint #44BB99
				187 204  51;... % pear #BBCC33
				170 170   0;... % olive #AAAA00
				221 221 221;... % pale grey #DDDDDD
				]./255;
		case 'sunset'
			colorPalette = [...
				 54  75 154;... #364B9A
				 74 123 183;... #4A7BB7
				110 166 205;... #6EA6CD
				152 202 225;... #98CAE1
				194 228 239;... #C2E4EF
				234 236 204;... #EAECCC
				254 218 139;... #FEDA8B
				253 179 102;... #FDB366
				246 126  75;... #F67E4B
				221  61  45;... #DD3D2D
				165   0  38;... #A50026
				255 255 255;... #FFFFFF - bad data
				]./255;
		case 'nightfall'
			colorPalette = [...
				 18  90  86;... #125A56
				  0 118 123;... #00767B
				 35 143 157;... #238F9D
				 66 157 198;... #42A7C6
				 96 188 233;... #60BCE9
				157 204 239;... #9DCCEF
				198 219 237;... #C6DBED
				222 230 231;... #DEE6E7
				236 234 218;... #ECEADA
				240 230 189;... #F0E6B2
				249 213 118;... #F9D576
				255 185  84;... #FFB954
				253 154  86;... #FD9A44
				245 118  52;... #F57634
				233  76  31;... #E94C1F
				209  24   7;... #D11807
				160  24  19;... #A01813
				255 255 255;... #FFFFFF - bad data
				]./255;
		case 'burd'
			colorPalette = [...
				 33 102 172;... #2166AC
				 67 147 195;... #4393C3
				146 197 222;... #92C5DE
				209 229 240;... #D1E5F0
				247 247 247;... #F7F7F7
				253 219 199;... #FDDBC7
				244 165 130;... #F4A582
				214  96  77;... #D6604D
				178 24   43;... #B2182B
				255 238 153;... # FFEE99 - bad data
				]./255;
		case 'prgn'
			colorPalette = [...
				118  42 131;... #762A83
				153 112 171;... #9970AB
				194 165 207;... #C2A5CF
				231 212 232;... #E7D4E8
				247 247 247;... #F7F7F7
				217 240 211;... #D9F0D3
				172 211 158;... #ACD39E
				 90 174  97;... #5AAE61
				 27 120  55;... #1B7837
				255 238 153;... #FFEE99 - bad data
				]./255;
		case 'ylorbr'
			colorPalette = [...
				255 255 229;... #FFFFE5
				255 247 188;... #FFF7BC
				254 227 145;... #FEE391
				254 196  79;... #FEC44F
				251 154  41;... #FB9A29
				236 112  20;... #EC7014
				204  76   2;... #CC4C02
				153  52   4;... #993404
				102  37   6;... #662506
				136 136 136;... #888888 - bad data
				]./255;
		case 'iridescent'
			colorPalette = [...
				254 251 233;... #FEFBE9
				252 247 213;... #FCF7D5
				245 243 193;... #F5F3C1
				234 240 181;... #EAF0B5
				221 236 191;... #DDECBF
				208 231 202;... #D0E7CA
				194 227 210;... #C2E3D2
				181 221 216;... #B5DDD8
				168 216 220;... #A8D8DC
				155 210 255;... #9BD2E1
				141 203 228;... #8DCBE4
				129 196 231;... #81C4E7
				123 188 231;... #7BBCE7
				126 178 228;... #7EB2E4
				136 156 221;... #88A5DD
				147 152 210;... #9398D2
				155 138 196;... #9B8AC4
				157 125 178;... #9D7DB2
				154 112 158;... #9A709E
				144  99 136;... #906388
				128  87 112;... #805770
				104  73  87;... #684957
				 70  53  58;... #46353A
				153 153 153;... #999999 - bad data
				]./255;
		case 'incandescent'
			colorPalette = [...
				206 255 255;... #CEFFFF
				198 247 214;... #C6F7D6
				162 244 155;... #A2F49B
				187 228  83;... #BBE453
				213 206   4;... #D5CE04
				231 181   3;... #E7B503
				241 153   3;... #F19903
				246 121  11;... #F6790B
				249  73   2;... #F94902
				228   5  21;... #E40515
				168   0   3;... #A80003
				136 136 136;... #888888 - bad data
				]./255;
		case 'discrete-rainbow'
			colorPalette = [...
				209 187 215;... #D1BBD7
				174 118 163;... #AE76A3
				136  46 114;... #882E72
				 25 101 176;... #1965B0
				 82 137 199;... #5289C7
				123 175 222;... #7BAFDE
				 78 178 101;... #4EB265
				144 201 135;... #90C987
				202 224 171;... #CAE0AB
				247 240  86;... #F7F056
				246 193  65;... #F6C141
				241 147  45;... #F1932D
				232  96  28;... #E8601C
				220   5  12;... #DC050C
				119 119 119;... #777777 - bad data
				]./255;
		case 'discrete-rainbow-23'
			colorPalette = [...
				232 236 251;... #E8ECFB
				217 204 227;... #D9CCE3
				202 172 203;... #CAACCB
				186 141 180;... #BA8DB4
				170 111 158;... #AA6F9E
				153  79 136;... #994F88
				136  46 114;... #882E72
				 25 101 176;... #1965B0
				 67 125 191;... #437DBF
				 97 149 207;... #6195CF
				123 175 222;... #7BAFDE
				 78 178 101;... #4EB265
				144 201 135;... #90C987
				202 224 171;... #CAE0AB
				247 240  86;... #F7F056
				247 203  69;... #F7CB45
				244 167  54;... #F4A736
				238 128  38;... #EE8026
				230  85  24;... #E65518
				220   5  12;... #DC050C
				165  23  14;... #A5170E
				114  25  14;... #72190E
				 66  21  10;... #42150A
				119 119 119;... #777777 - bad data
				]./255;
		case 'smooth-rainbow'
			colorPalette = [...
				232 236 251;... #E8ECFB
				221 216 239;... #DDD8EF
				209 193 225;... #D1C1E1
				195 168 209;... #C3A8D1
				181 143 194;... #B58FC2
				167 120 180;... #A778B4
				155  98 167;... #9B62A7
				140  78 153;... #8C4E99
				111  76 155;... #6F4C9B
				 96  89 169;... #6059A9
				 85 104 184;... #5568B8
				 78 121 197;... #4E79C5
				 77 138 198;... #4D8AC6
				 78 150 188;... #4E96BC
				 84 158 179;... #549EB3
				 89 165 169;... #59A5A9
				 96 171 158;... #60AB9E
				105 177 144;... #69B190
				119 183 125;... #77B77D
				140 188 104;... #8CBC68
				166 190  84;... #A6BE54
				190 188  72;... #BEBC48
				209 181  65;... #D1B541
				221 170  60;... #DDAA3C
				228 156  57;... #E49C39
				231 140  53;... #E78C35
				230 121  50;... #E67932
				228  99  45;... #E4632D
				223  72  40;... #DF4828
				218  34  34;... #DA2222
				184  34  30;... #B8221E
				149  33  27;... #95211B
				114  30  23;... #721E17
				 82  26  19;... #521A13
				102 102 102;... #666666 - bad data
				]./255;
	end

	if(nargin<1 || isempty(index))
		index = 1:size(colorPalette,1);
	end

	color = colorPalette(index,:);
end