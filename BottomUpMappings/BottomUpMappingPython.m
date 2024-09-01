classdef BottomUpMappingPython < BottomUpMappingBase
%BOTTOMUPMAPPINGPYTHON Bottom-Up Mapping from given Python function
%	BOTTOMUPMAPPINGPYTHON takes an appropriate function and respective 
%	parameters and uses that as bottom-up mapping for systems design.
%	To use this, you must ensure you have a compatible version of Python 
%	installed. You can see which versions are compatible with your MATLAB 
%	installation at:
%	- https://www.mathworks.com/support/requirements/python-compatibility.html
%	You must also have whatever modules you want to use installed for that 
%	version.
%
%   This class is derived from BottomUpMappingBase.
%
%	BOTTOMUPMAPPINGPYTHON methods:
%		- response : a function that receives the design samples and returns 
%       the system response for those variables in terms of performance and 
%       physical feasibility.
%
%   See also BottomUpMappingBase.
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

	properties(SetAccess = protected)
		%PYMODSYSTEM Python Module with Imported Library and Main Response Function
		%	PYMODSYSTEM is a Python module with the imported function given in 
		%	SystemFunctionName.
		%
		%	PYMODSYSTEM : py.module
		%
		%	See also SystemFunctionName.
		PyModSystem

		%PYMODPHYSICALFEASIBILITY Python Module with Physical Feasibility Function
		%	PYMODPHYSICALFEASIBILITY is a Python module with the imported function given
		%	in PhysicalFeasibilityFunction.
		%
		%	PYMOD : py.module
		%
		%	See also PhysicalFeasibilityFunction.
		PyModPhysicalFeasibility
	end

	properties
		%SYSTEMFUNCTIONNAME Python function name for system performance / physical feasibility
		%	SYSTEMFUNCTIONNAME is the name of the function to be used to get the system 
		%	response, given during initialization with the constructor. It must be given
		%	as a string, with the format '<name>.py'.
		%	
		%	This function must accept the following inputs (in order):
		%		- designSample : an array of size (nSample,nDesignVariable).
		%		- systemParameter : constant system parameters (type: 
		%		function-dependent).
		%	And must have the following outputs (in order):
		%		- performanceMeasure : an array of size (nSample,nPerformanceIndicator).
		%		- [optional] physicalFeasibilityMeasure : an array of size (nSample,
		%		nComponent).
		%
		%	SYSTEMFUNCTIONNAME : string
		%
		%	See also PyModSystem.
		SystemFunctionName

		%SYSTEMPARAMETER Constant system parameters
		%	SYSTEMPARAMETER is a generic variable that is used to store constant system 
		%	parameters, or in other words, values that do not change within the analysis
		%	of one system (but may change between different systems).
		%	SYSTEMPARAMETER can have any variable type. It is given during initilization 
		%	with the constructor.
		%
		%	SYSTEMPARAMETER : function-dependent
		%
		%	See also SystemFunctionName.
		SystemParameter

		%PHYSICALFEASIBILITYDEFAULTVALUE Measure used in the absence of specific output
		%	PHYSICALFEASIBILITYDEFAULTVALUE is used when the function used for system 
		%	performance does not include an output for the physical feasibility measure,
		%	nor was any specific function for physical feasibility been given.
		%
		%	PHYSICALFEASIBILITYDEFAULTVALUE : double OR (1,nComponent) double
		%
		%	See also SystemFunctionName, PhysicalFeasibilityFunction.
		PhysicalFeasibilityDefaultValue

		%PHYSICALFEASIBILITYFUNCTION Python function name
		%	PHYSICALFEASIBILITYFUNCTION is the name of the function used to get the
		%	physical feasibility measures of the system. This is an optional input that
		%	can be used in case the main response function does not include physical 
		%	feasibility information, but that is still relevant to the system. It must 
		%	be given as a string, with the format '<name>.py'.
		%	
		%	This function must accept the following inputs (in order):
		%		- designSample : an array of size (nSample,nDesignVariable).
		%		- physicalFeasibilityParameter : constant system parameters (type: function-dependent)
		%	And must have the following outputs (in order):
		%		- physicalFeasibilityMeasure : an array of size (nSample,nComponent).
		%
		%	PHYSICALFEASIBILITYFUNCTION : function_handle
		%
		%	See also PhysicalFeasibilityParameter.
		PhysicalFeasibilityFunction

		%PHYSICALFEASIBILITYPARAMETER Constant system parameters
		%	PHYSICALFEASIBILITYPARAMETER is a generic variable that is used to store 
		%	constant system parameters relevant specifically to computation of physical
		%	feasibility measures.
		%	PHYSICALFEASIBILITYPARAMETER can have any variable type. It is optionally 
		%	given during initilization with the constructor. If not given, it assumes
		%	the same value as the main Parameter.
		%
		%	PHYSICALFEASIBILITYPARAMETER : function-dependent
		%
		%	See also Parameter, PhysicalFeasibilityFunction.
		PhysicalFeasibilityParameter
	end

	methods
        function obj = BottomUpMappingPython(SystemfunctionName,varargin)
		%BOTTOMUPMAPPINGPYTHON Constructor
		%   BOTTOMUPMAPPINGPYTHON uses a given Python function (specified with its name
		%	as a string) and other parameters to create a bottom-up mapping. 
        %
        %   OBJ = BOTTOMUPMAPPINGPYTHON(SystemFUNCTIONNAME) receives the defined function 
        %	in SystemFUNCTIONNAME and returns a BottomUpMappingPython object in OBJ. In 
        %	said case, all other input arguments are assumed to have their default 
        %	values.
        %
        %   OBJ = BOTTOMUPMAPPINGPYTHON(SystemFUNCTIONNAME,PARAMETER) also allows the use of 
        %	constant system parameters when computing with the main function, passed in 
        %	PARAMETER. Default value is empty.
        %
        %   OBJ = BOTTOMUPMAPPINGPYTHON(...,NAME1,VALUE1,...) also allows
        %   for setting custom options for the computation procedure, passed as 
        %   name-value pairs. These are:
        %       - 'PhysicalFeasibilityDefaultValue' : default value used for the measure
        %		when no information regarding physical feasibility is available. 
        %		Default: -1.
        %       - 'PhysicalFeasibilityFunction' : function handle to be used to compute
        %		physical feasibility measures when that is not a part of the main 
        %		response function. Default: [].
        %       - 'PhysicalFeasibilityParameter' : constant system parameters used
        %		exclusively when computing physical feasibility measures with the
        %		correspondent specific function. By default it is given the same value 
        %		as PARAMETER.
        %
        %   Inputs:
        %       - FUNCTIONHANDLE : function_handle
        %       - SYSTEMPARAMETER : function-dependent
        %       - 'PhysicalFeasibilityDefaultValue' : double OR (1,nPhysicalFeasibility)
        %		double
        %       - 'PhysicalFeasibilityFunction' : function_handle
        %       - 'PhysicalFeasibilityParameter' : function-dependent
        %   
        %   Outputs:
        %       - OBJ : BottomUpMappingPython
        %
        %	See also response.
        
			parser = inputParser;
			parser.addOptional('SystemParameter',[]);
			parser.addParameter('PhysicalFeasibilityDefaultValue',-1);
			parser.addParameter('PhysicalFeasibilityFunction',[]);
			parser.addParameter('PhysicalFeasibilityParameter',[]);
			parser.parse(varargin{:});

			obj.SystemFunctionName = erase(SystemfunctionName,'.py');
    		pyMod = py.importlib.import_module(obj.SystemFunctionName);
    		obj.PyModSystem = py.importlib.reload(pyMod);

    		obj.SystemParameter = parser.Results.SystemParameter;

    		obj.PhysicalFeasibilityDefaultValue = parser.Results.PhysicalFeasibilityDefaultValue;


    		if(~isempty(parser.Results.PhysicalFeasibilityFunction))
    			obj.PhysicalFeasibilityFunction = erase(parser.Results.PhysicalFeasibilityFunction,'.py');
	    		pyModPhysicalFeasibility = py.importlib.import_module(obj.PhysicalFeasibilityFunction);
	    		obj.PyModPhysicalFeasibility = py.importlib.reload(pyModPhysicalFeasibility);
    		end

    		obj.PhysicalFeasibilityParameter = conditional_default_value_assignment(...
				parser.Results.PhysicalFeasibilityParameter,parser.Results.SystemParameter);
		end

		function [performanceMeasure,physicalFeasibilityMeasure,output] = response(obj,designSample)
		%RESPONSE Getting system response based on the design variables
		%	RESPONSE is the main function of any bottom-up mapping, where design samples 
		%	are given and the system response in terms of performance and physical 
		%	feasibility is returned.
        %
        %   PERFORMANCEMEASURE = OBJ.RESPONSE(DESIGNSAMPLE) returns the performance
        % 	measure of each given design in PERFORMANCEMEASURE. This is given as an 
        %	array with each row being correspondent to the respective design.
        %
        %   [PERFORMANCEMEASURE,PHYSICALFEASIBILITYMEASURE] = OBJ.RESPONSE(...) also 
        %	returns the physical feasibility measure of each given design in 
        %	PHYSICALFEASIBILITYMEASURE; similar to PERFORMANCEMEASURE, this is given as 
        %	an array with each row being correspondent to the respective design.
        %
        %   [PERFORMANCEMEASURE,PHYSICALFEASIBILITYMEASURE,SYSTEMOUTPUT] = 
        %	OBJ.RESPONSE(...) also returns any more information from the system response 
        %	computation in SYSTEMOUTPUT; for this class, it is empty.
        %
        %   Inputs:
        %       - OBJ : BottomUpMappingPython
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %   
        %   Outputs:
        %       - PERFORMANCEMEASURE : (nSample,nPerformance) double
        %       - PHYSICALFEASIBILITYDEFICIT : (nSample,nPhysicalFeasibility) double
        %       - SYSTEMOUTPUT : empty
        %
        %   See also CarCrashPythonBottomUpMapping.

			% reload to be safe - avoids errors that can occur if you change the function
		    pyMod = py.importlib.reload(obj.PyModSystem);
		    
		    % Call the function, passing design variables and parameters as inputs
		    results = pyMod.(obj.SystemFunctionName)(designSample,obj.SystemParameter);

		    if(isa(results,'py.tuple'))
		        performanceMeasure = double(results{1});
		        physicalFeasibilityMeasure = double(results{2});
		    else
		        performanceMeasure = double(results);
                physicalFeasibilityMeasure = [];
		    end

		    % physical feasibility analysis
            if(isempty(physicalFeasibilityMeasure))
                if(isempty(obj.PhysicalFeasibilityFunction))
                    physicalFeasibilityMeasure = repmat(obj.PhysicalFeasibilityDefaultValue,size(designSample,1),1);
                else
                    pyModPhysicalFeasibility = py.importlib.reload(obj.PyModPhysicalFeasibility); 
                    physicalFeasibilityMeasure = pyModPhysicalFeasibility.(...
                    	obj.PhysicalFeasibilityFunction)(designSample,obj.PhysicalFeasibilityParameter);
                end
            end

		    % no special output
			output = [];
		end
	end
end