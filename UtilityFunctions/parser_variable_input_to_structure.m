function outputStructure = parser_variable_input_to_structure(varargin)
%PARSER_VARIABLE_INPUT_TO_STRUCTURE Turn name-value pair arguments to structure
%	PARSER_VARIABLE_INPUT_TO_STRUCTURE takes name-value pair arguments and
%	structures and unifies them into a single structure. 
%	Names are assumed to be case insensitive.
%
%	OUTPUTSTRUCTURE = PARSER_VARIABLE_INPUT_TO_STRUCTURE(NAME1,VALUE1,...)
%	builds a structure OUTPUTSTRUCTURE with fieldnames as specified in NAME
%	inputs and associated values VALUE. 
%
%	OUTPUTSTRUCTURE = PARSER_VARIABLE_INPUT_TO_STRUCTURE(...,STRUCTURE,...)
%	expands the given structure and takes its own fieldnames and values 
%	as name-value pair arguments. 
%
%	Input:
%		- VARARGIN: name-value pair arguments OR structures
%
%	Output:
%		- OUTPUTSTRUCTURE : structure
%
%   See also inputParser.
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
	parser.CaseSensitive = false;
	parser.KeepUnmatched = true;
	parser.StructExpand = true;
	parser.parse(varargin{:});
	outputStructure = parser.Unmatched;
end