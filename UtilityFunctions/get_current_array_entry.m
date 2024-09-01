function output = get_current_array_entry(array,index)
%GET_CURRENT_ARRAY_ENTRY Return value with given index or the last available one
%   GET_CURRENT_ARRAY_ENTRY either returns the value with the specified index
%   from the given array, or the final array value in case the index is out-of-
%   bounds.
%
%   OUTPUT = GET_CURRENT_ARRAY_ENTRY(ARRAY,INDEX) assigns to OUTPUT the value of
%   ARRAY in index INDEX. If said INDEX is larger than the size of ARRAY, the
%   last value in ARRAY is returned instead. If array is multi-dimensional,
%   index is counted in a column-major fashion.
%
%   Input:
%       - ARRAY : (m,n) array
%       - INDEX : integer
%   Output:
%       - OUTPUT : array value
%
%   See also end.
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

    if(index>=numel(array))
        output = array(end);
    else
        output = array(index);
    end
end