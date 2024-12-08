%setup_sso_toolbox Initial setup for the use of the SSO toolbox
%   setup_sso_toolbox adds to the MATLAB path all the folders of the SSO 
%   Toolbox, allowing the use of the functions in them without any extra steps.
%   It also adds the Python folder for bottom-up mappings to the Py-specific
%   path if a compatible Python environment is found.
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


%% Cleanup
close all;
fclose all;
clear all;
clc;
more off;
diary off;


%% get folder of setup
setupFilePath = mfilename('fullpath');
if contains(setupFilePath,'LiveEditorEvaluationHelper')
    setupFilePath = matlab.desktop.editor.getActiveFilename;
end
toolboxFolderPath = replace(erase(setupFilePath,mfilename),'\','/');
addpath(toolboxFolderPath);
clear setupFilePath


%% Add Folders (and potentially Subfolders) to Path
% 'foldername',includeSubfolder
addDirectory = {...
    ... % public folders
    'BottomUpMappings/',true;...
    'Data/',true;...
    'SolutionSpacesToolbox/',true;...
    'SupportFunctions/',true;...
    'SurrogateModelingToolbox/',true;...
    'UtilityFunctions/',true;...
    ... % private folders
    'PrivateData/',false;...
    'PrivateBottomUpMappings/',true};

for i=1:size(addDirectory,1)
    currentPath = [toolboxFolderPath,addDirectory{i,1}];

    if(addDirectory{i,2})
        % also add subfolders
        totalCurrentPath = genpath(currentPath);
    else
        % add only main folder
        totalCurrentPath = currentPath;
    end

    if(exist(currentPath, 'dir'))
       addpath(totalCurrentPath);
    end
end
clear addDirectory i currentPath totalCurrentPath


%% Python Setup
pythonEnviroment = pyenv; % check if environment is present
if(~isempty(pythonEnviroment.Version) && ~strcmpi(pythonEnviroment.Version,''))
    % add bottom-up mapping folder to MATLAB path - unnecessary for functionality
    % addpath('./BottomUpMappingsPython/');
    
    % add bottom-up mapping folder to Python path
    % see if folder with files is in python path already; if not, add
    pythonBottomUpMappingFolder = [toolboxFolderPath,'BottomUpMappingsPython/'];
    compatibleBottomUpMappingFolder = replace(pythonBottomUpMappingFolder,'/','\\'); % python extracts paths with two backslashes
    pythonPath = split(char(py.sys.path),','); % get py path in cell matlab format
    if(~any(contains(pythonPath,compatibleBottomUpMappingFolder))) 
        % insert it at the start of path list
        insert(py.sys.path,int32(0),pythonBottomUpMappingFolder);
    end
    
    clear pythonBottomUpMappingFolder compatibleBottomUpMappingFolder pythonPath
end
clear pythonEnviroment


%% clear remaining variables
clear toolboxFolderPath  

