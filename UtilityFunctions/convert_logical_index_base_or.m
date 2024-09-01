function varargout = convert_logical_index_base_or(varargin)
%CONVERT_LOGICAL_INDEX_BASE_OR Perform multiple 'or' and change base of arrays 
%	CONVERT_LOGICAL_INDEX_BASE_OR takes many boolean arrays of equal size and 
%	produces the result of their logical 'or' operation, as well as the 
%	equivalent input arrays but expressed in the new size/base.
%	
%	OROUTPUT = CONVERT_LOGICAL_INDEX_BASE_OR(ARRAY1,ARRAY2,...) performs an
%	logical 'or' operation between all input arrays ARRAY1,ARRAY2,... and 
%	returns the output of said operation applied to all inputs in OROUTPUT.
%
%	[OROUTPUT,NEWARRAY1,NEWARRAY2,...] = CONVERT_LOGICAL_INDEX_BASE_OR(...)
%	also returns the input arrays' equivalent logical indexing in the output
%	arrays NEWARRAY1,NEWARRAY2,.... This means the following: given a generic
%	array X and its new base newX = X(OROUTPUT), the following holds true:
%	newX(NEWARRAY1)==X(ARRAY1), newX(NEWARRAY2)==X(ARRAY2), etc.
%
%	Input:
%		- ARRAY : (n,1) logical
%	Output:
%		- OROUTPUT : (m,1) logical
%		- NEWARRAY : (m,1) logical
%
%   See also or, varargin, varargout.
%
%   Copyright 2024 Eduardo Rodrigues Della Noce
%   SPDX-License-Identifier: Apache-2.0

	% find final output logical array
	jFinalOutput = varargin{1};
	nEntry = size(varargin,2); 
	for i=2:nEntry
		jFinalOutput = (jFinalOutput|varargin{i});
	end
	
	% create a linear index array with the conversion of indeces from 
	% old base to new base; for example, if you want to know what 
	% index 3 in the old base is in the new base, you can look for
	% iFinalOutput(3), or iFinalOutput([false false true])
	sizeInput = size(jFinalOutput,1);
	iConvertIndex = nan(sizeInput,1);
	sizeFinalOutput = sum(jFinalOutput);
	iConvertIndex(jFinalOutput) = 1:sizeFinalOutput;

	% for each entry, find the conversion between old and new;
	% use the converter built above for that. meaning here, for each entry,
	% get all the indices where value was true in the old base, and
	% set them to true in the new base
	varargout{1} = jFinalOutput;
	for i=1:nEntry
		varargout{i+1} = false(sizeFinalOutput,1);
		varargout{i+1}(iConvertIndex(varargin{i})) = true;
	end
end