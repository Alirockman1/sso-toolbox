function newValue = conditional_default_value_assignment(value,defaultValue,condition)
%CONDITIONAL_DEFAULT_VALUE_ASSIGNMENT Assign value based on given condition
%	CONDITIONAL_DEFAULT_VALUE_ASSIGNMENT will return a default value if the 
%	condition applied on its input evaluates as 'true'; otherwise, the output 
%	retains the value of the input.
%
%	NEWVALUE = CONDITIONAL_DEFAULT_VALUE_ASSIGNMENT(VALUE,DEFAULTVALUE) returns  
%	as NEWVALUE the default output DEFAULTVALUE if the input VALUE is empty;
%	if not, it returns the input VALUE itself instead.
%
%	NEWVALUE = CONDITIONAL_DEFAULT_VALUE_ASSIGNMENT(VALUE,DEFAULTVALUE,
%	CONDITION) allows for the specification of the condition CONDITION to be 
%	used to determine whether the default output should be used or not. If said  
%	condition applied to the input VALUE is evaluated as 'true', the default 
%	output DEFAULTVALUE is returned; otherwise, like before, the input VALUE is 
%	returned instead.
%
%   Inputs:
%       - VALUE : any
%       - DEFAULTVALUE : any
%       - CONDITION : function_handle
%
%   Output:
%       - NEWVALUE : any
%
%   See also: isempty, ternary, ternary_array.
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

	if(nargin<3)
		condition = @isempty;
	end
	newValue = ternary(condition(value),defaultValue,value);
end