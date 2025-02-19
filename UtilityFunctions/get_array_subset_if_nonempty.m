function newArray = get_array_subset_if_nonempty(array, rowIndices, columnIndices)
%GET_ARRAY_SUBSET_IF_NONEMPTY Return array subset or empty based on array emptiness
%	GET_ARRAY_SUBSET_IF_NONEMPTY returns the subset of the input array specified by
%	row and column indices if the array is not empty. If the array is empty,
%	returns an empty array. Empty row or column indices will return all rows/columns.
%
%	NEWARRAY = GET_ARRAY_SUBSET_IF_NONEMPTY(ARRAY, ROWINDICES, COLUMNINDICES) returns:
%	   - ARRAY(ROWINDICES, COLUMNINDICES) when ARRAY is not empty
%	   - Empty array [] when ARRAY is empty
%	If ROWINDICES is empty, all rows are selected. If COLUMNINDICES is empty or
%	omitted, all columns are selected.
%
%   Inputs:
%       - ARRAY : any
%           Array to check for emptiness
%       - ROWINDICES : vector | []
%           Row indices to select (empty for all rows)
%       - COLUMNINDICES : vector | [] (optional)
%           Column indices to select (empty or omitted for all columns)
%
%   Output:
%       - NEWARRAY : same type as ARRAY | []
%           Subset of array or empty based on array emptiness
%
%   See also: isempty, conditional_default_value_assignment.
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

    if isempty(array)
        newArray = [];
        return;
    end
    
    % Handle column indices
    if nargin < 3 || isempty(columnIndices)
        columnIndices = ':';
    end
    % Convert empty indices to "all" selector
    if nargin < 2 || isempty(rowIndices)
        rowIndices = ':';
    end
    
    newArray = array(rowIndices, columnIndices);
end 