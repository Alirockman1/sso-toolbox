function structureExtracted = structure_extract_fields(structure,desiredFields,fieldMarkers,nestedFields)
%STRUCTURE_EXTRACT_FIELDS Extract specific fields / patterns from a structure
%   STRUCTURE_EXTRACT_FIELDS can be used to create a new structure with 
%   specified fields from another structure. It also contains ways of
%   copying fields that contain certain patterns, and fields with are nested
%   structures themselves.
%   
%   STRUCTUREEXTRACTED = STRUCTURE_EXTRACT_FIELDS(STRUCTURE) copies all fields
%   from the initial structure STRUCTURE to the new structure 
%   STRUCTUREEXTRACTED.
%
%   STRUCTUREEXTRACTED = STRUCTURE_EXTRACT_FIELDS(STRUCTURE,DESIREDFIELDS) only
%   copies the fields DESIREDFIELDS from the initial STRUCTURE to the new
%   STRUCTUREEXTRACTED. DESIREDFIELDS may be specified as a single-column cell,
%   in which case the field names remain the same, or it may be specified as
%   a two-column cell, in which case the names will be changed from the name
%   specified in the first column to the name from the second column when
%   copied to STRUCTUREEXTRACTED.
%
%   STRUCTUREEXTRACTED = STRUCTURE_EXTRACT_FIELDS(STRUCTURE,DESIREDFIELDS,
%   FIELDMARKERS) additionally copies all fields that have markers at the start
%   of their names as specified in FIELDMARKERS; when copied, the marker itself
%   gets erased. For example, if there is a field in STRUCTURE called 
%   'Trial_Case', and 'Trial_' is set as a marker in FIELDMARKERS, the value
%   of that field gets copied into STRUCTUREEXTRACTED with the name 'Case'.
%
%   STRUCTUREEXTRACTED = STRUCTURE_EXTRACT_FIELDS(STRUCTURE,DESIREDFIELDS,
%   FIELDMARKERS,NESTEDFIELDS) additionally copies all the nested fields from
%   specified fields NESTEDFIELDS into the new structure. For example, if
%   'Trial' is specified as a nested field in NESTEDFIELDS and STRUCTURE has
%   a field which is itself a structure like 'Trial.Case','Trial.Date', then
%   those values get copied into STRUCTUREEXTRACTED as 'Case' and 'Date' 
%   respectively.
%
%   Input:
%       - STRUCTURE  : structure
%       - DESIREDFIELDS : (d,1) OR (d,2) cell
%       - FIELDMARKERS : (m,1) cell
%       - NESTEDFIELDS : (n,1) cell
%   Output:
%       - STRUCTUREEXTRACTED : structure
%
%   See also struct, fieldnames.
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

    if(nargin<4)
        nestedFields = [];
        if(nargin<3)
            fieldMarkers = [];
            if(nargin<2)
                desiredFields = fieldnames(structure);
            end
        end
    end
    
    % if only one input and not in cell, change to cell type
    if(isstring(desiredFields) || ischar(desiredFields))
        desiredFields = {desiredFields};
    end
    if(isstring(fieldMarkers) || ischar(fieldMarkers))
        fieldMarkers = {fieldMarkers};
    end
    if(isstring(nestedFields) || ischar(nestedFields))
        nestedFields = {nestedFields};
    end
    
    % assume either one or two column inputs
    if(size(desiredFields,2)>2)
        desiredFields = desiredFields';
    end
    
    % assume single-column input
    if(size(fieldMarkers,2)>1)
        fieldMarkers = fieldMarkers';
    end
    if(size(nestedFields,2)>1)
        nestedFields = nestedFields';
    end
    
    % get names of fields
    if(size(desiredFields,2)==1)
        oldNames = desiredFields;
        newNames = desiredFields;
    elseif(size(desiredFields,2)==2)
        oldNames = desiredFields(:,1);
        newNames = desiredFields(:,2);
    end
    
    % initialize
    structureExtracted = struct;
    
    % copy main values required
    for i=1:size(oldNames,1)
        if(isfield(structure,oldNames{i}))
            structureExtracted.(newNames{i}) = structure.(oldNames{i});
        end
    end
    
    % copy values with marker attached at the start
    originalFields = fieldnames(structure);
    for i=1:size(fieldMarkers,1)
        casematch = startsWith(originalFields,fieldMarkers{i},'IgnoreCase',true);
        indCaseMatch = find(casematch);
        
        for j=1:size(indCaseMatch,1)
            oldName = originalFields{indCaseMatch(j)};
            newName = erase(oldName,fieldMarkers{i});
            
            if(isempty(newName))
                structureExtracted.(oldName) = structure.(oldName);
            else
                structureExtracted.(newName) = structure.(oldName);
            end
        end
    end
    
    % search for nested fields
    for i=1:size(nestedFields,1)
        casematch = strcmpi(originalFields,nestedFields{i});
        
        if(any(casematch))
            currentField = originalFields{casematch};
            if(isstruct(structure.(currentField)))
                nested = fieldnames(structure.(currentField));
                for j=1:size(nested,1)
                    structureExtracted.(nested{j}) = structure.(currentField).(nested{j});
                end
            end
        end
    end
end