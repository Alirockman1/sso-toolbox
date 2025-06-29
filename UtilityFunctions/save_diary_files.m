function saveFolder = save_diary_files(folderName,archiveFiles,additionalSearchFolders)
%SAVE_DIARY_FILES Create a folder where console is logged and files can be saved
%   SAVE_DIARY_FILES creates a folder where all data generated from the 
%   currently running script may be saved. Also initializes 'diary', so all 
%   console input/output is also saved in a .txt file.
%
%   SAVEFOLDER = SAVE_DIARY_FILES(FOLDERNAME) creates a folder named FOLDERNAME
%   with the following structure, starting from the primary working directory: 
%   '/RESULTS/<FOLDERNAME>/Data-<Timestamp>/'. It initializes 'diary' as a text 
%   file in said folder as well, writing the console outputs to this file. The 
%   path to that folder is returned as SAVEFOLDER. For FOLDERNAME, a suggestion 
%   is to use the name of the script, which can be obtained with 'mfilename'.
%
%   SAVEFOLDER = SAVE_DIARY_FILES(FOLDERNAME,ARCHIVEFILES) also creates an 
%   additional folder '/Files/' inside SAVEFOLDER and creates copies of files 
%   ARCHIVEFILES and puts said copies there. This can be used for archival 
%   purposes. Files are searched in all folders currently in the 'path' and
%   the primary working directory.
%
%   SAVEFOLDER = SAVE_DIARY_FILES(FOLDERNAME,ARCHIVEFILES,
%   ADDITIONALSEARCHFOLDERS) allows for the specification of additional folders
%   where the files of ARCHIVEFILES may be found.
%
%   Input:
%       - FOLDERNAME : string
%       - ARCHIVEFILES : (nFile,1) cell
%       - ADDITIONALSEARCHFOLDERS : (nFolder,1) cell
%
%   Output:
%       - SAVEFOLDER : string
%
%   See also diary, mfilename, mkdir, path, isfile, copyfile, get_date_time.
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

    %% Process Inputs
    if(nargin<3)
        additionalSearchFolders = {};
        if(nargin<2)
            archiveFiles = {};
        end
    end
    
    
    %% Create Folder where Information will be stored
    saveFolder = ['./RESULTS/',folderName,'/Data-',get_date_time,'/'];
    mkdir(saveFolder);
    

    %% Start diary (Console Transcription)
    diary([saveFolder,'Console Output.txt']);
    

    %% Find Files and Copy Them to Folder
    if(~isempty(archiveFiles))
        visibleFolders = [pwd;split(path,';')];
        visibleFolders = [additionalSearchFolders;visibleFolders];

        archiveFolder = [saveFolder,'Files/'];
        mkdir(archiveFolder);
        for i=1:length(archiveFiles)
            % Test each combination of possible path
            foundFile = false;
            for j=1:length(visibleFolders)
                % Check in which folder the file is located
                candidatePath = [visibleFolders{j},archiveFiles{i}];
                if(isfile(candidatePath))
                    archiveFiles{i} = candidatePath;
                    foundFile = true;
                    break;
                end
            end
            
            % File not found -> error
            if(~foundFile)
                warning('File %s not found, please make sure it exists, is in a visible folder and is named correctly.',archiveFiles{i})
            else
                copyfile(archiveFiles{i},archiveFolder);
            end
        end
    end
end