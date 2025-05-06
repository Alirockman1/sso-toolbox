classdef BottomUpMappingFixedVariables < BottomUpMappingBase
%BOTTOMUPMAPPINGFIXEDVARIABLES Bottom-Up Mapping with fixed design variables
%	BOTTOMUPMAPPINGFIXEDVARIABLES takes another bottom-up mapping and fixes certain
%   design variables to predefined values. When the response method is called, it
%   completes the design variable array with the fixed values before passing it to
%   the base mapping.
%
%   This class is derived from BottomUpMappingBase.
%
%	BOTTOMUPMAPPINGFIXEDVARIABLES methods:
%		- response : a function that receives the design sample points, completes them
%       with fixed values for specified variables, and returns the system response
%       for those variables in terms of performance and physical feasibility.
%
%   See also BottomUpMappingBase, BottomUpMappingFunction, BottomUpMappingPython.
%
%   Copyright 2025 Eduardo Rodrigues Della Noce
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
		%BASEBOTTOMUPMAPPING Base bottom-up mapping to be used
		%	BASEBOTTOMUPMAPPING is the underlying bottom-up mapping that will be called after
		%	completing the design variables with fixed values.
		%
		%	BASEBOTTOMUPMAPPING : BottomUpMappingBase
		%
		%	See also IsFixedVariable, FixedVariableValues.
		BaseBottomUpMapping

		%ISFIXEDVARIABLE Logical index of fixed design variables
		%	ISFIXEDVARIABLE is the logical index which indicates which design
		%	variables are fixed (true) as opposed to free (false).
		%	This is used during the response method to reconstruct the complete design
		%	variable input as expected by the base mapping.
		%
		%	ISFIXEDVARIABLE : (1,nDesignVariable) logical
		%
		%	See also BaseBottomUpMapping, FixedVariableValues.
		IsFixedVariable

		%FIXEDVARIABLEVALUES Values for the fixed design variables
		%	FIXEDVARIABLEVALUES contains the values to be used for the fixed design
		%	variables when completing the design sample.
		%
		%	FIXEDVARIABLEVALUES : (1,nFixedDesignVariable) double
		%
		%	See also IsFixedVariable.
		FixedVariableValues
	end

	methods
		function [obj,freeLowerBound,freeUpperBound,freeInitialDesign,freeComponentIndex] = BottomUpMappingFixedVariables(baseBottomUpMapping, isFixedVariable, fixedVariableValues, varargin)
		%BOTTOMUPMAPPINGFIXEDVARIABLES Constructor
		%   BOTTOMUPMAPPINGFIXEDVARIABLES creates an object that wraps another bottom-up
		%   mapping and fixes certain design variables to predefined values.
        %
        %   OBJ = BOTTOMUPMAPPINGFIXEDVARIABLES(BASEBOTTOMUPMAPPING, ISFIXEDVARIABLE, 
        %   FIXEDVARIABLEVALUES) creates a new BottomUpMappingFixedVariables object
        %   that uses BASEBOTTOMUPMAPPING as the underlying mapping. ISFIXEDVARIABLE is a
        %   logical array indicating which design variables are fixed, and 
        %   FIXEDVARIABLEVALUES contains the values for those fixed variables.
        %
        %   OBJ = BOTTOMUPMAPPINGFIXEDVARIABLES(BASEBOTTOMUPMAPPING, ISFIXEDVARIABLE, 
        %   FIXEDVARIABLEVALUES, DESIGNSPACELOWERBOUND, DESIGNSPACEUPPERBOUND) also 
        %   adapts the design space boundaries DESIGNSPACELOWERBOUND and 
        %   DESIGNSPACEUPPERBOUND by removing the fixed variables.
        %
        %   OBJ = BOTTOMUPMAPPINGFIXEDVARIABLES(BASEBOTTOMUPMAPPING, ISFIXEDVARIABLE, 
        %   FIXEDVARIABLEVALUES, DESIGNSPACELOWERBOUND, DESIGNSPACEUPPERBOUND, 
        %   INITIALDESIGN) also adapts the initial design INITIALDESIGN by removing 
        %   the fixed variables.
        %
        %   OBJ = BOTTOMUPMAPPINGFIXEDVARIABLES(BASEBOTTOMUPMAPPING, ISFIXEDVARIABLE, 
        %   FIXEDVARIABLEVALUES, DESIGNSPACELOWERBOUND, DESIGNSPACEUPPERBOUND, 
        %   INITIALDESIGN, COMPONENTINDEX) also adapts the component index 
        %   COMPONENTINDEX by converting indices to account for the removal of fixed 
        %   variables.
        %
        %   [OBJ,FREELOWERBOUND] = BOTTOMUPMAPPINGFIXEDVARIABLES(...) additionally 
        %   returns the adapted lower bound of the design space with fixed variables 
        %   removed.
        %
        %   [OBJ,FREELOWERBOUND,FREEUPPERBOUND] = BOTTOMUPMAPPINGFIXEDVARIABLES(...) 
        %   additionally returns the adapted upper bound of the design space with fixed 
        %   variables removed.
        %
        %   [OBJ,FREELOWERBOUND,FREEUPPERBOUND,FREEINITIALDESIGN] = 
        %   BOTTOMUPMAPPINGFIXEDVARIABLES(...) additionally returns the adapted initial 
        %   design with fixed variables removed.
        %
        %   [OBJ,FREELOWERBOUND,FREEUPPERBOUND,FREEINITIALDESIGN,FREECOMPONENTINDEX] = 
        %   BOTTOMUPMAPPINGFIXEDVARIABLES(...) additionally returns the adapted 
        %   component index with indices converted to account for the removal of fixed 
        %   variables.
        %
        %   Inputs:
        %       - BASEBOTTOMUPMAPPING : BottomUpMappingBase
        %       - ISFIXEDVARIABLE : (1,nDesignVariable) logical
        %       - FIXEDVARIABLEVALUES : (1,nFixedVariable) double
        %       - DESIGNSPACELOWERBOUND : (1,nDesignVariable) double
        %       - DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double
        %       - INITIALDESIGN : (1,nDesignVariable) double
        %       - COMPONENTINDEX : (1,nComponent) cell
        %   
        %   Outputs:
        %       - OBJ : BottomUpMappingFixedVariables
        %       - FREELOWERBOUND : (1,nFreeDesignVariable) double
        %       - FREEUPPERBOUND : (1,nFreeDesignVariable) double
        %       - FREEINITIALDESIGN : (1,nFreeDesignVariable) double
        %       - FREECOMPONENTINDEX : (1,nComponent) cell
        %
        %	See also response, convert_index_base.

			% Parse input arguments
			parser = inputParser;
			parser.addOptional('DesignSpaceLowerBound', []);
			parser.addOptional('DesignSpaceUpperBound', []);
			parser.addOptional('InitialDesign', []);
			parser.addOptional('ComponentIndex', {});
			parser.parse(varargin{:});
			options = parser.Results;

			% Validate inputs
			if(~isa(baseBottomUpMapping, 'BottomUpMappingBase'))
				error('BottomUpMappingFixedVariables:InvalidBaseMapping', ...
					'Base mapping must be a BottomUpMappingBase object');
			end
			
			if(~islogical(isFixedVariable))
				error('BottomUpMappingFixedVariables:InvalidIsFixedVariable', ...
					'Fixed variable index must be a logical array');
			end
			
			if(length(fixedVariableValues) ~= sum(isFixedVariable))
				error('BottomUpMappingFixedVariables:InvalidFixedVariableValues', ...
					'Number of fixed variable values must match the number of fixed variables');
			end
			
			% Assign properties
			obj.BaseBottomUpMapping = baseBottomUpMapping;
			obj.IsFixedVariable = isFixedVariable;
			obj.FixedVariableValues = fixedVariableValues;

			% Adapt design space lower bound
			if ~isempty(options.DesignSpaceLowerBound)
				freeLowerBound = options.DesignSpaceLowerBound(~isFixedVariable);
			else
				freeLowerBound = [];
			end

			% Adapt design space upper bound
			if ~isempty(options.DesignSpaceUpperBound)
				freeUpperBound = options.DesignSpaceUpperBound(~isFixedVariable);
			else
				freeUpperBound = [];
			end

			% Adapt initial design
			if ~isempty(options.InitialDesign)
				freeInitialDesign = options.InitialDesign(~isFixedVariable);
			else
				freeInitialDesign = [];
			end

			% Adapt component index
			freeComponentIndex = cell(size(options.ComponentIndex));
			for i = 1:length(freeComponentIndex)
				if(~isempty(options.ComponentIndex{i}))
					convertedIndex = convert_index_base(~isFixedVariable(:), options.ComponentIndex{i}(:), 'forward');

					freeComponentIndex{i} = convertedIndex(~isnan(convertedIndex));
				end
			end
		end

		function [performanceMeasure, physicalFeasibilityMeasure, systemOutput] = response(obj, designSample)
		%RESPONSE Getting system response based on the design sample points
		%	RESPONSE completes the design sample with fixed values for specified variables
		%	before passing it to the base mapping's response method.
        %
        %   PERFORMANCEMEASURE = OBJ.RESPONSE(DESIGNSAMPLE) returns the performance
        % 	measure of each given design in PERFORMANCEMEASURE. This is returned as an 
        %	array with each row being correspondent to the respective design.
        %
        %   [PERFORMANCEMEASURE,PHYSICALFEASIBILITYMEASURE] = OBJ.RESPONSE(...) also 
        %	returns the physical feasibility measure of each given design in 
        %	PHYSICALFEASIBILITYMEASURE; similar to PERFORMANCEMEASURE, this is given as 
        %	an array with each row being correspondent to the respective design.
        %
        %   [PERFORMANCEMEASURE,PHYSICALFEASIBILITYMEASURE,SYSTEMOUTPUT] = 
        %	OBJ.RESPONSE(...) also returns any more information from the system response 
        %	computation in SYSTEMOUTPUT; this is directly passed from the base mapping's
        %   response method.
        %
        %   Inputs:
        %       - OBJ : BottomUpMappingFixedVariables
        %       - DESIGNSAMPLE : (nSample,nFreeDesignVariable) double
        %   
        %   Outputs:
        %       - PERFORMANCEMEASURE : (nSample,nPerformance) double
        %       - PHYSICALFEASIBILITYDEFICIT : (nSample,nPhysicalFeasibility) double
        %       - SYSTEMOUTPUT : class-defined
        %
        %   See also BottomUpMappingBase.

			% Get dimensions
			nSample = size(designSample, 1);
			nDesignVariable = length(obj.IsFixedVariable);
			nFreeVariables = sum(~obj.IsFixedVariable);
			
			% Validate input dimensions
			if size(designSample, 2) ~= nFreeVariables
				error('BottomUpMappingFixedVariables:InvalidDesignSampleSize', ...
					'Design sample must have %d columns (number of free variables)', nFreeVariables);
			end
			
			% Create complete design sample with fixed values
			completeDesignSample = zeros(nSample, nDesignVariable);
			
			% Fill in the free variables
			completeDesignSample(:, ~obj.IsFixedVariable) = designSample;
			
			% Fill in the fixed variables (replicate for each sample)
			completeDesignSample(:, obj.IsFixedVariable) = repmat(obj.FixedVariableValues, nSample, 1);
			
			% Call the base mapping's response method
			[performanceMeasure, physicalFeasibilityMeasure, systemOutput] = obj.BaseBottomUpMapping.response(completeDesignSample);
		end
	end
end 