function output = ternary_array(conditionArray,outputTrue,outputFalse)
%TERNARY_ARRAY Element-wise conditional ternary operator (if-else)
%   TERNARY_ARRAY returns its second argument when the first evaluates as true,
%	and the third when not. This is done individually for each element of the
%	condition array.
%
%	OUTPUT = TERNARY_ARRAY(CONDITIONARRAY,OUTPUTTRUE,OUTPUTFALSE) assigns to 
%	each element of OUTPUT the value OUTPUTTRUE if the respective condition  
%	CONDITIONARRAY is 'true' for that element, and assigns the value OUTPUTFALSE
%	if it is not.
%
%   Input:
%       - CONDITION : (n,1) logical
%		- OUTPUTTRUE : (n,1) any
%		- OUTPUTFALSE : (n,1) any
%   Output:
%       - OUTPUT : (n,1) any
%
%   See also ternary, conditional_default_value_assignment.
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

	output = outputFalse;
	output(conditionArray) = outputTrue(conditionArray);
end