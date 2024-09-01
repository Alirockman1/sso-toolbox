classdef BottomUpMappingFunction < BottomUpMappingBase
%BOTTOMUPMAPPINGFUNCTION Bottom-Up Mapping from given MATLAB system function
%	BOTTOMUPMAPPINGFUNCTION takes an appropriate MATLAB function and respective 
%	parameters and uses that as bottom-up mapping for systems design.
%
%   This class is derived from BottomUpMappingBase.
%
%	BOTTOMUPMAPPINGFUNCTION methods:
%		- response : a function that receives the design sample points and  
%       returns the system response for those variables in terms of performance 
%       and physical feasibility.
%
%   See also BottomUpMappingBase, BottomUpMappingPython.
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

	properties
		%SYSTEMFUNCTION MATLAB function handle for system response
		%	SYSTEMFUNCTION is the handle of the function to be used to get the system 
		%	response, given during initialization to the constructor.
		%	
		%	This function must accept the following inputs (in order):
		%		- designSample : an array of size (nSample,nDesignVariable).
		%		- [optional] systemParameter : constant system parameters (type: 
		%		function-dependent)
		%	And must have the following outputs (in order):
		%		- performanceMeasure : an array of size (nSample,nPerformanceMeasure).
		%		- [optional] physicalFeasibilityMeasure : an array of size (nSample,
		%		nComponent).
		%		- [optional] systemOutput : any other output which may be desired 
		%		(normally to be saved).
		%
		%	SYSTEMFUNCTION : function_handle
		%
		%	See also Parameter.
		SystemFunction

		%SYSTEMPARAMETER Constant system parameters
		%	SYSTEMPARAMETER is a generic variable that is used to store constant system 
		%	parameters, or in other words, values that do not change within the analysis
		%	of one system.
		%	SYSTEMPARAMETER can have any variable type. It is given during initilization
		%	to the constructor.
		%
		%	SYSTEMPARAMETER : function-dependent
		%
		%	See also SystemFunction.
		SystemParameter

		%PHYSICALFEASIBILITYDEFAULTVALUE Measure used in the absence of specific output
		%	PHYSICALFEASIBILITYDEFAULTVALUE is used when the function used for system 
		%	performance does not include an output for the physical feasibility measure,
		%	nor was any specific function for physical feasibility been given.
		%
		%	PHYSICALFEASIBILITYDEFAULTVALUE : double OR (1,nComponent) double
		%
		%	See also SystemFunction, PhysicalFeasibilityFunction.
		PhysicalFeasibilityDefaultValue

		%PHYSICALFEASIBILITYFUNCTION MATLAB function handle for physical feasibility
		%	PHYSICALFEASIBILITYFUNCTION is the handle of the function used to get the
		%	physical feasibility measures of the system. This is an optional input that
		%	can be used in case the main response function does not include physical 
		%	feasibility information, but that is still relevant to the system.
		%	This function must accept the following inputs (in order):
		%		- designSample : an array of size (nSample,nDesignVariable).
		%		- [optional] physicalFeasibilityParameter : constant system parameters 
		%		(type: function-dependent)
		%	And must have the following outputs (in order):
		%		- physicalFeasibilityMeasure : an array of size (nSample,nComponent).
		%		- physicalFeasibilityOutput : function-dependent.
		%
		%	PHYSICALFEASIBILITYFUNCTION : function_handle
		%
		%	See also SystemFunction.
		PhysicalFeasibilityFunction

		%PHYSICALFEASIBILITYPARAMETER Constant system parameters of physical feasibility
		%	PHYSICALFEASIBILITYPARAMETER is a generic variable that is used to store 
		%	constant system parameters relevant specifically to computation of physical
		%	feasibility measures.
		%	PHYSICALFEASIBILITYPARAMETER can have any variable type. It is optionally 
		%	given during initilization with the constructor. 
		%
		%	PHYSICALFEASIBILITYPARAMETER : function-dependent
		%
		%	See also SystemParameter.
		PhysicalFeasibilityParameter
	end

	methods
		function obj = BottomUpMappingFunction(SystemFunction,varargin)
		%BOTTOMUPMAPPINGFUNCTION Constructor
		%   BOTTOMUPMAPPINGFUNCTION uses a given MATLAB function handle which computes 
		%	system response and other parameters to create a bottom-up mapping. 
        %
        %   OBJ = BOTTOMUPMAPPINGFUNCTION(SYSTEMFUNCTION) receives the defined function 
        %	in SYSTEMFUNCTION and returns a BottomUpMappingFunction object in OBJ.
        %
        %   OBJ = BOTTOMUPMAPPINGFUNCTION(SYSTEMFUNCTION,SYSTEMPARAMETER) also allows
        %   the use of constant system parameters when computing with the main function, 
        %	passed in SYSTEMPARAMETER. Default value is empty.
        %
        %   OBJ = BOTTOMUPMAPPINGFUNCTION(...,NAME1,VALUE1,...) also allows
        %   for setting custom options for the computation procedure, passed as 
        %   name-value pair arguments. These are:
        %       - 'PhysicalFeasibilityDefaultValue' : default value used for the measure
        %		when no information regarding physical feasibility is available.
        %		Default: -1.
        %       - 'PhysicalFeasibilityFunction' : function handle to be used to compute
        %		physical feasibility measures when that is not a part of the main 
        %		response function. Default: [].
        %       - 'PhysicalFeasibilityParameter' : constant system parameters used
        %		exclusively when computing physical feasibility measures with the
        %		correspondent specific function. Default value is the same as 
        %		SYSTEMPARAMETER.
        %
        %   Inputs:
        %       - SYSTEMFUNCTION : function_handle
        %       - SYSTEMPARAMETER : function-dependent
        %       - 'PhysicalFeasibilityDefaultValue' : double OR (1,nPhysicalFeasibility)
        %		double
        %       - 'PhysicalFeasibilityFunction' : function_handle
        %       - 'PhysicalFeasibilityParameter' : function-dependent
        %   
        %   Outputs:
        %       - OBJ : BottomUpMappingFunction
        %
        %	See also response.

			parser = inputParser;
			parser.addOptional('SystemParameter',[]);
			parser.addParameter('PhysicalFeasibilityDefaultValue',-1);
			parser.addParameter('PhysicalFeasibilityFunction',[]);
			parser.addParameter('PhysicalFeasibilityParameter',[]);
			parser.parse(varargin{:});

			obj.SystemFunction = SystemFunction;
			obj.SystemParameter = parser.Results.SystemParameter;

			obj.PhysicalFeasibilityDefaultValue = parser.Results.PhysicalFeasibilityDefaultValue;
			obj.PhysicalFeasibilityFunction = parser.Results.PhysicalFeasibilityFunction;
			obj.PhysicalFeasibilityParameter = conditional_default_value_assignment(...
				parser.Results.PhysicalFeasibilityParameter,parser.Results.SystemParameter);
		end

		function [performanceMeasure,physicalFeasibilityMeasure,systemOutput] = response(obj,designSample)
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
        %	computation in SYSTEMOUTPUT; for this class, it directly passes the extra 
        %	'output' from the main system MATLAB function and physical feasibility 
        %	function, if it has any. 
        %
        %   Inputs:
        %       - OBJ : DesignBottomUpMappingFunction
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %   
        %   Outputs:
        %       - PERFORMANCEMEASURE : (nSample,nPerformance) double
        %       - PHYSICALFEASIBILITYDEFICIT : (nSample,nPhysicalFeasibility) double
        %       - SYSTEMOUTPUT : structure
        %			-- 'SystemFunctionOutput' : function-dependent
        %			-- 'PhysicalFeasibilityFunctionOutput' : function-dependent
        %
        %   See also car_crash_2d, distance_to_center, two_ellipsoids.

			if(nargin(obj.SystemFunction)==1)
				systemFunctionInputArgument = {designSample};
			else
				systemFunctionInputArgument = {designSample,obj.SystemParameter};
			end

			if(~isempty(obj.PhysicalFeasibilityFunction))
				if(nargin(obj.PhysicalFeasibilityFunction)==1)
					physicalFeasibilityFunctionInputArgument = {designSample};
				else
					physicalFeasibilityFunctionInputArgument = {designSample,obj.PhysicalFeasibilityParameter};
				end
			end

			% performance analysis
			systemFunctionOutput = [];
			if(nargout(obj.SystemFunction)==1) % physical feasibility is not considered
			    performanceMeasure = obj.SystemFunction(systemFunctionInputArgument{:});
                physicalFeasibilityMeasure = [];
			elseif(nargout(obj.SystemFunction)==2)
			    [performanceMeasure,physicalFeasibilityMeasure] = obj.SystemFunction(systemFunctionInputArgument{:});
			else
				[performanceMeasure,physicalFeasibilityMeasure,systemFunctionOutput] = obj.SystemFunction(systemFunctionInputArgument{:});
			end

			% physical feasibility analysis
			physicalFeasibilityFunctionOutput = [];
            if(isempty(physicalFeasibilityMeasure))
			    if(isempty(obj.PhysicalFeasibilityFunction))
				    physicalFeasibilityMeasure = repmat(obj.PhysicalFeasibilityDefaultValue,size(designSample,1),1);
			    elseif(nargout(obj.PhysicalFeasibilityFunction)==1)
				    physicalFeasibilityMeasure = obj.PhysicalFeasibilityFunction(physicalFeasibilityFunctionInputArgument{:});
			    else
				    [physicalFeasibilityMeasure,physicalFeasibilityFunctionOutput] = ...
				    	obj.PhysicalFeasibilityFunction(designSample,obj.PhysicalFeasibilityParameter);
                end
            end

			% special output from function if applicable
			systemOutput.SystemFunctionOutput = systemFunctionOutput;
			systemOutput.PhysicalFeasibilityFunctionOutput = physicalFeasibilityFunctionOutput;
		end
	end
end