classdef ConsoleLogging
%CONSOLELOGGING Logging utility for displaying information on the console
%	CONSOLELOGGING allows for logging functionality in MATLAB similar to that
%	of functions like log4j and python.logging, done directly in the MATLAB 
%	console.
%	
%	When the constructor is called, CONSOLELOGGING gets a given severity level 
%	of logging associated with it. After that, only calls with a logging level 
%	higher than the specified one are actually written to console. The severity 
%	levels available are:
%		- (1) 'ALL'
%		- (2) 'TRACE'
%		- (3) 'DEBUG'
%		- (4) 'INFO'
%		- (5) 'WARN'
%		- (6) 'CRITICAL'
%		- (7) 'OFF'
%	For example, if CONSOLELOGGING was initialized with a severity logging level
%	of 'INFO', only calls of the 'critical', 'warn' and 'info' methods will 
%	produce an output in the console.
%	
%	CONSOLELOGGING properties:
%		- SeverityLevelName : names of all logging levels, in severity order
%		- SeverityLogLevel : set level of logging severity for the object
%
%	CONSOLELOGGING methods: 
%		- critical : writes a message with level CRITICAL on the console
%		- warn : writes a message with level WARN on the console
%		- info : writes a message with level INFO on the console
%		- debug : writes a message with level DEBUG on the console
%		- trace : writes a message with level TRACE on the console
%
%	For similar functionality that writes logs to a file instead of the console,
%	check 'log4m' by Luke Winslow:
%	https://www.mathworks.com/matlabcentral/fileexchange/37701-log4m
%
%   See also error, warning, fprintf.
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

	properties (Constant)
		%SEVERITYLEVELNAME Logging severity level names
		%	SEVERITYLEVELNAME contains the names for the levels of logging available,
		%	in order of highest to lowest severity.
		%
		%	SEVERITYLEVELNAME : (1,8) cell
		%
		%	See also SeverityLogLevel.
		SeverityLevelName = {'ALL', 'TRACE', 'DEBUG', 'INFO', 'WARN', 'CRITICAL', 'OFF'};
	end

	properties
		%SEVERITYLOGLEVEL Object logging severity level
		%	SEVERITYLOGLEVEL sets the minimum level of severity for messages to be
		%	displayed on the console.
		%
		%	SEVERITYLOGLEVEL : integer
		%
		%	See also SeverityLevelName.
		SeverityLogLevel
	end

	methods
		function obj = ConsoleLogging(SeverityLogLevel)
		%CONSOLELOGGING Constructor
		%   CONSOLELOGGING initializes a logger with the given severity level.
        %
        %   OBJ = CONSOLELOGGING(LOGLEVEL) receives a severity level in LOGLEVEL
        %	and returns a logger in OBJ with which only messages with severity
        %	level equal or higher than LOGLEVEL are displayed. LOGLEVEL can be
        %	given as a case-insensitive string (for example, 'info') or as an
        %	integer (for example, 4).
        %
        %   Inputs:
        %       - LOGLEVEL : string OR integer
        %
        %   Outputs:
        %       - OBJ : ConsoleLogging
        %   
        %   See also trace, debug, info, warn, critical.
			if(isnumeric(SeverityLogLevel))
            	obj.SeverityLogLevel = SeverityLogLevel;
            else
            	iLevelChoice = find(strcmpi(SeverityLogLevel,ConsoleLogging.SeverityLevelName));
            	if(isempty(iLevelChoice))
            		error('ConsoleLogging:Constructor:LevelNotFound','Logging level not accepted.');
            	else
            		obj.SeverityLogLevel = iLevelChoice;
            	end
            end
		end

		function trace(obj,varargin)
		%TRACE Trace-level messages
		%   TRACE displays the given message on the console if the severity log
		%	level is set to show messages with said level of severity or less.
		%	The message is shown via the 'fprintf' function, and as such, the 
		%	given arguments should be compatible with that.
        %
        %   TRACE(MESSAGE) displays the MESSAGE in the console if log severity
        %	level is set accordingly.
        %
        %   Inputs:
        %       - MESSAGE : 'fprintf' input argument
        %   
        %   See also fprintf, debug, info, warn, critical.
			currentLevel = ConsoleLogging.get_method_level(dbstack);
			if(obj.SeverityLogLevel<=currentLevel)
				fprintf(varargin{:});
			end
		end

		function debug(obj,varargin)
		%DEBUG Debug-level messages
		%   DEBUG displays the given message on the console if the severity log
		%	level is set to show messages with said level of severity or less.
		%	The message is shown via the 'fprintf' function, and as such, the 
		%	given arguments should be compatible with that.
        %
        %   DEBUG(MESSAGE) displays the MESSAGE in the console if log severity
        %	level is set accordingly.
        %
        %   Inputs:
        %       - MESSAGE : 'fprintf' input argument
        %   
        %   See also fprintf, trace, info, warn, critical.
			currentLevel = ConsoleLogging.get_method_level(dbstack);
			if(obj.SeverityLogLevel<=currentLevel)
				fprintf(varargin{:});
			end
		end

		function info(obj,varargin)
		%INFO Info-level messages
		%   INFO displays the given message on the console if the severity log
		%	level is set to show messages with said level of severity or less.
		%	The message is shown via the 'fprintf' function, and as such, the 
		%	given arguments should be compatible with that.
        %
        %   INFO(MESSAGE) displays the MESSAGE in the console if log severity
        %	level is set accordingly.
        %
        %   Inputs:
        %       - MESSAGE : 'fprintf' input argument
        %   
        %   See also fprintf, trace, debug, warn, critical.
			currentLevel = ConsoleLogging.get_method_level(dbstack);
			if(obj.SeverityLogLevel<=currentLevel)
				fprintf(varargin{:});
			end
		end

		function warn(obj,varargin)
		%WARN Warn-level messages
		%   WARN displays the given message on the console if the severity log
		%	level is set to show messages with said level of severity or less.
		%	The message is shown via the 'warning' function, and as such, the 
		%	given arguments should be compatible with that.
        %
        %   WARN(MESSAGE) displays the MESSAGE in the console if log severity
        %	level is set accordingly.
        %
        %   Inputs:
        %       - MESSAGE : 'warning' input argument
        %   
        %   See also warning, trace, debug, info, critical.
			currentLevel = ConsoleLogging.get_method_level(dbstack);
			if(obj.SeverityLogLevel<=currentLevel)
				warning(varargin{:});
			end
		end

		function critical(obj,varargin)
		%CRITICAL Critical-level messages
		%   CRITICAL displays the given message on the console if the severity 
		%	log level is set to show messages with said level of severity or 
		%	less.
		%	The message is shown via the 'error' function, and as such, the 
		%	given arguments should be compatible with that.
        %
        %   CRITICAL(MESSAGE) displays the MESSAGE in the console if log 
        %	severity level is set accordingly.
        %
        %   Inputs:
        %       - MESSAGE : 'error' input argument
        %   
        %   See also error, trace, debug, info, warn.
            currentLevel = ConsoleLogging.get_method_level(dbstack);
			if(obj.SeverityLogLevel<=currentLevel)
				error(varargin{:});
			end
        end
    end
    
    methods (Static, Access = private)
        function currentLevel = get_method_level(functionStack)
        %GET_METHOD_LEVEL Return severity level based on method name
		%   GET_METHOD_LEVEL returns the severiy level of a method based on its
		%	name; said name is compared against the order set in 
		%	'SeverityLevelName'.
        %
        %   CURRENTLEVEL = CONSOLELOGGING.GET_METHOD_LEVEL(FUNCTIONSTACK) uses 
        %	the function information extracted from its stack, in particular its
        %	name, to determine its severity level CURRENTLEVEL based on the  
        %	other definitions of the class.
        %
        %   Input:
        %       - FUNCTIONSTACK : struct
        %
        %	Output:
        %		- CURRENTLEVEL : integer
        %   
        %   See also dbstack, SeverityLevelName.
			currentMethod = erase(functionStack(1).name,'ConsoleLogging.');
			currentLevel = find(strcmpi(currentMethod,ConsoleLogging.SeverityLevelName));
        end
    end
end