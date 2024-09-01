function [mergedStructure,mergedCell] = merge_name_value_pair_argument(varargin)
%MERGE_NAME_VALUE_PAIR_ARGUMENT merges its entries as name-value pair arguments
%	MERGE_NAME_VALUE_PAIR_ARGUMENT takes an arbitrary number of inputs (which 
%	may be structures or cells with name-value pair arguments) and merges them 
%	into a single structure/cell with said name-value pairs.
%
%	MERGEDSTRUCTURE = MERGE_NAME_VALUE_PAIR_ARGUMENT(INPUT1,INPUT2,...) crawls
%	through all inputs and analyses them as name-pair arguments, merging them
%	into the same structure; if two inputs have properties with the same name, 
%	the inputs that come later in the input have higher priority and their 
%	values will overwrite the previous ones.
%	These inputs can be:
%		- Structures.
%		- Cells with the format {'Name1',Value1,'Name2',Value2,...}.
%		- Maps 
%		- Dictionaries.
%	The result of the merging is returned as a structure in MERGEDSTRUCTURE.
%
%	[MERGEDSTRUCTURE,MERGEDCELL] = MERGE_NAME_VALUE_PAIR_ARGUMENT(...) does the
%	same, but also returns the output as a cell with name-value pair entries 
%	in MERGEDCELL (same format as specified above).
%
%	See also struct, namedargs2cell, parser_variable_input_to_structure.
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

	% merge all entries into a structure
	currentStructure = {};
	for i=1:size(varargin,2)
		currentEntry = varargin{i};
		currentNames = {};
		currentValues = {};
		currentStructure{i} = struct;

		if(isa(currentEntry,'cell'))
			currentNames = {currentEntry{1:2:end-1}}'; % first entry of each pair
			currentValues = {currentEntry{2:2:end}}'; % second entry of each pair
	    elseif(isa(currentEntry,'struct'))
	    	currentNames = fieldnames(varargin{i});
	    	for j=1:size(currentNames,1)
	    		currentValues{j} = currentEntry.(currentNames{j});
	    	end
	    elseif(isa(currentEntry,'containers.Map') || isa(currentEntry,'dictionary'))
	    	currentNames = currentEntry.keys';
	    	currentValues = currentEntry.values';
	    else
	    	error('MergeNameValuePair:TypeNotIdentified','Entry %d not a type identified.',i);
	    end

	    for j=1:size(currentNames,1)
	        currentStructure{i}.(currentNames{j}) = currentValues{j};
	    end
	end

	% unwrap and reorganize into cell if necessary
	mergedStructure = parser_variable_input_to_structure(currentStructure{:});
	if(nargout>1)
		mergedCell = namedargs2cell(mergedStructure);
	end
end