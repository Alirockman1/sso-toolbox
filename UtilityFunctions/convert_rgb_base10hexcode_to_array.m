function arrayRgb = convert_rgb_base10hexcode_to_array(hex10Rgb,nDigit)
%CONVERT_RGB_BASE10HEXCODE_TO_ARRAY Convert hex code in base 10 to RGB array
%   CONVERT_RGB_BASE10HEXCODE_TO_ARRAY converts an RGB value that is a hex code
%   given as an integer in base 10 to a [R,G,B] array. 
%   As an example of base 10 hex code, for RGB using two 'digits' (16^2 = 256   
%   shades for each color), the hex codes would convert to base 10 as follows:
%       - Pure red: #0000ff --> 255 (256^1-1)
%       - Pure green: #00ff00 --> 65280 (256^2-1 - 255)
%       - Pure blue: #ff0000 --> 16711680 (256^3-1 - 65280 - 255)
%   CONVERT_RGB_BASE10HEXCODE_TO_ARRAY takes integers like the ones on the 
%   right-hand side and converts them to the MATLAB standard [R,G,B] array,
%   where each row is one color and each column is a number between 0 and 1
%   indicating intensity.
%
%   ARRAYRGB = CONVERT_RGB_BASE10HEXCODE_TO_ARRAY(HEX10RGB) converts the hex
%   code written in base 10 HEX10RGB and converts it to the [R,G,B] array 
%   ARRAYRGB. It is assumed the hex code had 2 digits (so the number of
%   shades for each color is 256).
%
%   ARRAYRGB = CONVERT_RGB_BASE10HEXCODE_TO_ARRAY(HEX10RGB,NDIGIT) allows the
%   specification of the number of hex code digits NDIGIT.
%
%   Input:
%       - HEX10RGB : (nEntry,1) integer
%       - NDIGIT : integer OR (nEntry,1) integer
%
%   Output:
%       - ARRAYRGB : (nEntry,3) double
%
%   See also hex2rgb, rgb2hex.
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

    %
    if(nargin<2 || isempty(nDigit))
        nDigit = 2;
    end  

    baseDivider = round(16.^nDigit); % for 2 digits, 256

    % R
    quotientAfterRed = floor(hex10Rgb./baseDivider);
    red = round(hex10Rgb - quotientAfterRed.*baseDivider);
    
    % G
    quotientAfterGreen = floor(quotientAfterRed./baseDivider);
    green = round(quotientAfterRed - quotientAfterGreen.*baseDivider);
    
    % B
    quotientAfterBlue = floor(quotientAfterGreen./baseDivider); % should be 0
    blue = round(quotientAfterGreen - quotientAfterBlue.*baseDivider);

    % wrap
    arrayRgb = [red(:),green(:),blue(:)]/255;
end