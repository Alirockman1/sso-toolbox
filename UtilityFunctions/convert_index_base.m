function [linearConvertedQuery,logicalConvertedQuery] = convert_index_base(baseConversion,indexQuery,conversionType)
%CONVERT_INDEX_BASE Turn query indices to a different index base
%	CONVERT_INDEX_BASE takes a query array expressed in a base of 
%   elements and converts the indices expressed in that query array to
%	a different base, which is handled by a logical index array of conversion. 
%   The conversion can be either:
%       - A 'backward' conversion, where the query array is in the reduced base
%       and must be converted back to the full base.
%       - A 'forward' conversion, where the query array is in the full base
%       and must be converted forward to the reduced base.
%
%	LINEARCONVERTEDQUERY = CONVERT_INDEX_BASE(BASECONVERSION,INDEXQUERY,
%   CONVERSIONTYPE) takes the indices expressed in INDEXQUERY and transforms 
%   them into the equivalent indices in the respective full/reduced base 
%   depending on CONVERSIONTYPE. In either case, the base conversion is 
%   specified by BASECONVERSION, with 'true' values indicating where the 
%   reduced base has an equivalent index in the full base. Each column of  
%   INDEXQUERY is converted independently.
%       - In a 'backward' conversion : it is assumed INDEXQUERY is expressed
%       in the reduced base, and therefore, LINEARCONVERTEDQUERY will be the 
%       indices in the full base that correspond to the indices in INDEXQUERY.
%       - In a 'forward' conversion : it is assumed INDEXQUERY is expressed
%       in the full base, and therefore, LINEARCONVERTEDQUERY will be the 
%       indices in the reduced base that correspond to the indices in 
%       INDEXQUERY.
%   If INDEXQUERY is given as a many-column logical array, LINEARCONVERTEDQUERY
%   returns an array with 'NaN' values filling columns which had less query
%   elements then others.
%
%	[LINEARCONVERTEDQUERY,LOGICALCONVERTEDQUERY] = CONVERT_INDEX_BASE(...) also 
%	returns the logical indices in the new base for the query 
%	LOGICALCONVERTEDQUERY; if the query was a linear index array with repeated  
%	elements, that will not be reflected on LOGICALCONVERTEDQUERY.
%
%	Example 'backward': assume you have an array of numbers, and you want to 
%	know the indices of the ones that are negative, from the largest in absolute
%	value to the smallest. That can be done as follows:
%		manyNumbers = [-4 1 -2 6 -5 2 8 -9 -1 2 3];
%		jNegativeNumbers = (manyNumbers<0)';
%		[~,orderNegativeOnly] = sort(manyNumbers(jNegativeNumbers));
%		orderNegative = CONVERT_INDEX_BASE(jNegativeNumbers,orderNegativeOnly,'backward');
%		manyNumbers(orderNegative)
%	This process will return indices [8;5;1;3;9], which correspond to numbers
%	[-9,-5,-4,-2,-1] as displayed in the console, as we wanted.
%
%   Example 'forward': assume you want to know what the index of the '-2' entry 
%   of the array [-4 1 -2 6 -5 2 8 -9 -1 2 3] will be after you remove all 
%   non-negative numbers. This can be done as follows:
%       manyNumbers = [-4 1 -2 6 -5 2 8 -9 -1 2 3];
%       jNegativeNumbers = (manyNumbers<0)';
%       negativeNumbers = manyNumbers(jNegativeNumbers);
%       iMinusTwo = find(manyNumbers==-2,1,'first'); % result: 3
%       iMinusTwoNew = CONVERT_INDEX_BASE(jNegativeNumbers,iMinusTwo,'forward');
%   This process will return index 2, because in the new array of only negative
%   numbers [-4 -2 -5 -9 -1], '-2' is the second entry. 
%
%   Input:
%       - BASECONVERSION : (n,1) logical
%       - INDEXQUERY : (m,p) integer OR logical 
%       - CONVERSIONTYPE : char OR string
%   Output:
%       - LINEARCONVERTEDQUERY : (m,p) integer
%       - LOGICALCONVERTEDQUERY : (n,p) ('backward') OR (nT,p) ('forward') logical
%
%   See also convert_logical_index_base_or.
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
        conversionType = [];
    end

    if(strcmpi(conversionType,'backward'))
        % new arrays will have as many elements as the original base
        sizeNewBase = size(baseConversion,1);
        % connect the logical conversion to full-size linear indices
        iBaseConversion = [1:sizeNewBase]';
        iBaseConversion = iBaseConversion(baseConversion);
    else%if(strcmpi(conversionType,'forward'))
        % new arrays will have as many elements as there were 'true's in conversion
        sizeNewBase = sum(baseConversion);
        % new linear indices are only increased for each 'true' in the base
        iBaseConversion = cumsum(baseConversion);
        iBaseConversion(~baseConversion) = nan;
    end

    nConversionEntry = size(indexQuery,2);
    if(islogical(indexQuery))
        nElementEntry = sum(indexQuery,1);
    else
        nElementEntry = size(indexQuery,1);
    end

    % based on the conversion linear index, get the respective base linear indices
    linearConvertedQuery = nan(max(nElementEntry),nConversionEntry);
    logicalConvertedQuery = false(sizeNewBase,nConversionEntry);
    for i=1:size(indexQuery,2)
        nCurrentEntry = nElementEntry(i);
        linearConvertedQuery(1:nCurrentEntry,i) = iBaseConversion(indexQuery(:,i));

        % create the logical index from the linear one
        validIndex = (~isnan(linearConvertedQuery(:,i)));
        logicalConvertedQuery(linearConvertedQuery(validIndex,i),i) = true;
    end
end