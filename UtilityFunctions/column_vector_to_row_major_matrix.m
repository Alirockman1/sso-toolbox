function rowMajorMatrix = column_vector_to_row_major_matrix(columnVector,nColumn)
%COLUMN_VECTOR_TO_ROW_MAJOR_MATRIX Convert representation of entries
%	COLUMN_VECTOR_TO_ROW_MAJOR_MATRIX takes a column vector and changes its
%	shape to a row-major matrix.
%
%	For example, the column vector:
%		[1]
%		[2]
%		[3]
%		[4]
%		[5]
%		[6]
%	Would return, with 2 columns:
%		[1,2]
%		[3,4]
%		[5,6]
%	Or, with 3 columns:
%		[1,2,3]
%		[4,5,6]
%
%	ROWMAJORMATRIX = COLUMN_VECTOR_TO_ROW_MAJOR_MATRIX(COLUMNVECTOR,NCOLUMN) 
%	takes the column vector COLUMNVECTOR and returns its row-major matrix  
%	ROWMAJORMATRIX with the number of columns NCOLUMN.
%
%   Input:
%		- COLUMNVECTOR : (nEntry,1) double
%		- NCOLUMN : integer
%
%   Output:
%		- ROWMAJORMATRIX : (nEntry/nColumn,nColumn) double
%
%   See also row_major_matrix_to_column_vector.
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

	rowMajorMatrix = reshape(columnVector,nColumn,[])';
end