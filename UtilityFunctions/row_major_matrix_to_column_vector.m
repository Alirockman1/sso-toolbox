function columnVector = row_major_matrix_to_column_vector(rowMajorMatrix)
%ROW_MAJOR_MATRIX_TO_COLUMN_VECTOR Convert representation of entries
%	ROW_MAJOR_MATRIX_TO_COLUMN_VECTOR takes a row-major matrix and changes its
%	shape to a column vector.
%
%	For example, when using ROW_MAJOR_MATRIX_TO_COLUMN_VECTOR on the following 
%	matrices with two columns and three columns:
%		[1,2]	[1,2,3]
%		[3,4]	[4,5,6]
%		[5,6]
%	They will both return the column vector:
%		[1]
%		[2]
%		[3]
%		[4]
%		[5]
%		[6]
%	For these same examples, the traditional colon operator 'matrix(:)' would 
%	result in:
%		[1]   	[1]
%		[3]		[4]
%		[5]		[2]
%		[2]		[5]
%		[4]		[3]
%		[6]		[6]
%	For the first and second examples, respectively. 
%
%	COLUMNVECTOR = ROW_MAJOR_MATRIX_TO_COLUMN_VECTOR(ROWMAJORMATRIX) takes the
%	row-major matrix ROWMAJORMATRIX and returns its column vector form 
%	COLUMNVECTOR.
%
%   Input:
%		- ROWMAJORMATRIX : (nRow,nColumn) double
%
%   Output:
%		- COLUMNVECTOR : (nRow*nColumn,1) double
%
%   See also column_vector_to_row_major_matrix.
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

	columnVector = reshape(rowMajorMatrix',[],1);
end