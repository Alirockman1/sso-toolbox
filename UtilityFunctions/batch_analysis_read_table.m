function batchOptions = batch_analysis_read_table(filename)
%BATCH_ANALYSIS_READ_TABLE Process batch analysis table to structure
%   BATCH_ANALYSIS_READ_TABLE can read tables with the expected batch analysis
%   format and repackage the information in the form of a structure array.
%   The expected format for the table is:
%
%   ID  | ReferenceID | OPTION1NAME   | OPTION2NAME   | ...
%   id1 | refid1      | option1value1 | option2value1 | ...
%   id2 | refid2      | option1value2 | option2value2 | ...
%   .......................................................
%
%   This can be specified as an .csv file, an Excel file, or others. 
%   BATCH_ANALYSIS_READ_TABLE creates a structure array where each entry is
%   correspondent to a row in the original table, separating 'ID' and 
%   'ReferenceID' and merging the other options into the generic field 
%   'Options', which itself is a structure containing all the options.
%
%   BATCHOPTIONS = BATCH_ANALYSIS_READ_TABLE(FILENAME) read the file FILENAME 
%   and returns the structure with all options BATCHOPTIONS.
%
%   Input:
%       - FILENAME : char OR string
%
%   Output:
%       - BATCHPOPTIONS : (nTest,1) structure
%           - ID : char OR string
%           - ReferenceID : char OR string
%           - Options : structure
%
%   See also readtable.
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

    %% Read Excel Spreadsheets Information
    warning('OFF','MATLAB:table:ModifiedAndSavedVarnames'); % yes, I know, sheet names will be changed when converted to variables...
    % Design Variables
    opts = detectImportOptions(filename);
    opts = setvartype(opts, 'string');
    tableBatch = readtable(filename,opts);
    warning('ON','MATLAB:table:ModifiedAndSavedVarnames'); % ok, reset it now

    columnName = erase(tableBatch.Properties.VariableNames,{'ID','ReferenceID'});
    isEmpty = cellfun(@isempty,columnName);
    optionName = {columnName{~isEmpty}};

    %
    batchOptions = struct(...
        'ID',[],...
        'ReferenceID',[],...
        'Options',[]);

    % create structure
    nEntry = size(tableBatch,1);
    nOption = length(optionName);
    for i=1:nEntry
        entryCurrent.ID = char(tableBatch.ID{i});
        entryCurrent.ReferenceID = char(tableBatch.ReferenceID{i});

        for j=1:nOption
            optionCurrent = tableBatch.(optionName{j}){i};
            try % see if it is a numeric value, function handle, cell, ...
                entryCurrent.Options.(optionName{j}) = eval(optionCurrent);
            catch % if not, treat as string
                entryCurrent.Options.(optionName{j}) = char(optionCurrent);
            end
        end

        batchOptions(i) = entryCurrent;
    end
end


function d = readNumericEntries(s)
    commas = (s==',');
    s(commas) = '.';
    d = str2double(s);
end

