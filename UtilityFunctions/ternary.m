function output = ternary(condition,outputTrue,outputFalse)
%TERNARY Conditional ternary operator (if-else)
%   TERNARY returns its second argument when the first evaluates as true, and
%	the third when not.
%
%	OUTPUT = TERNARY(CONDITION,OUTPUTTRUE,OUTPUTFALSE) assigns to OUTPUT the
%	value OUTPUTTRUE if the condition CONDITION is 'true', and assigns the
%	value OUTPUTFALSE if it is not.
%
%   Input:
%       - CONDITION : logical
%		- OUTPUTTRUE : any
%		- OUTPUTFALSE : any
%   Output:
%       - OUTPUT : any
%
%   See also ternary_array, conditional_default_value_assignment.
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

	if(condition)
		output = outputTrue;
	else
		output = outputFalse;
	end
end