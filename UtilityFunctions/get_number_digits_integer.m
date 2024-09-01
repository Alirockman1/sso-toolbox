function nDigit = get_number_digits_integer(integerArray)
%GET_NUMBER_DIGITS_INTEGER Count how many digits each integer has
%   GET_NUMBER_DIGITS_INTEGER computes how many digits each integer has in its
%   base 10 representation. '9' has 1, '42' has 2, '255' has 3, etc.
%   
%   NDIGIT = GET_NUMBER_DIGITS_INTEGER(INTEGERARRAY) counts the number of digits
%   NDIGIT for each integer in INTEGERARRAY.
%
%   Input:
%       - INTEGERARRAY : (1,n) integer
%   Output:
%       - NDIGIT : (1,n) integer
%
%   See also log10.
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

    nDigit = floor(log10(abs(integerArray))) + 1;
    nDigit(isnan(nDigit) | isinf(nDigit)) = 1;
end