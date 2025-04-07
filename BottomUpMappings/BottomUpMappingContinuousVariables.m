classdef BottomUpMappingContinuousVariables < BottomUpMappingBase
    %BOTTOMUPMAPPINGCONTINUOUSVARIABLES Bottom-Up Mapping with continuous design variables
    %	BOTTOMUPMAPPINGCONTINUOUSVARIABLES takes another bottom-up mapping and fixes certain
    %   design variables to predefined values. When the response method is called, it
    %   completes the design variable array with the fixed values before passing it to
    %   the base mapping.
    %
    %   This class is derived from BottomUpMappingBase.
    %
    %	BOTTOMUPMAPPINGCONTINUOUSVARIABLES methods:
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
            InterpolationMethodFunction
    
            %FIXEDVARIABLEVALUES Values for the fixed design variables
            %	FIXEDVARIABLEVALUES contains the values to be used for the fixed design
            %	variables when completing the design sample.
            %
            %	FIXEDVARIABLEVALUES : (1,nFixedDesignVariable) double
            %
            %	See also IsFixedVariable.
            InterpolationMethodOptions

            %DESIGNVARIABLEINDEXMAPPING Mapping of design variables to indices
            %	DESIGNVARIABLEINDEXMAPPING is the mapping of design variables to indices.
            %
            %	DESIGNVARIABLEINDEXMAPPING : (1,nDesignVariable) double
            %
            %	See also DiscreteLowerBound, DiscreteUpperBound, DiscreteInitialDesign, DiscreteComponentIndex.
            DesignVariableIndexMapping

            % BaseVariablesDomain Domain of the base variables
            %   BASEVARIABLESDOMAIN is the domain of the base variables.
            %
            %   BASEVARIABLESDOMAIN : (1,nBaseVariable) double
            %
            %   See also.
            BaseVariablesStencils
        end
    
        methods
            function [obj,functionalDiscreteBound,discreteInitialDesign,discreteComponentIndex] = BottomUpMappingContinuousVariables(baseBottomUpMapping, baseVariableDesignSpace, functionalDesignSpace, initialDesign, varargin)
            %BOTTOMUPMAPPINGCONTINUOUSVARIABLES Constructor
            %   BOTTOMUPMAPPINGCONTINUOUSVARIABLES creates an object that wraps another bottom-up
            %   mapping and fixes certain design variables to predefined values.
            %
            %   OBJ = BOTTOMUPMAPPINGCONTINUOUSVARIABLES(BASEBOTTOMUPMAPPING, ISFIXEDVARIABLE, 
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
                parser.addOptional('ComponentIndex', {});
                parser.addParameter('InterpolationMethodFunction',@interp1);
                parser.addParameter('InterpolationMethodOptions',{'pchip','extrap'});
                parser.addParameter('NumberOfStencils',1);
                parser.parse(varargin{:});
                options = parser.Results;
    
                % Validate inputs
                if(~isa(baseBottomUpMapping, 'BottomUpMappingBase'))
                    error('BottomUpMappingFixedVariables:InvalidBaseMapping', ...
                        'Base mapping must be a BottomUpMappingBase object');
                end
                
                % Assign properties
                obj.BaseBottomUpMapping = baseBottomUpMapping;
                obj.InterpolationMethodFunction = options.InterpolationMethodFunction;
                obj.InterpolationMethodOptions = options.InterpolationMethodOptions;

                % for the number of stencils given, discretize design space,
                %   initial design and components
                nDesignVariable = size(functionalDesignSpace,2);

                nStencil = options.NumberOfStencils;
                if(isscalar(nStencil))
                    nStencil = repmat(nStencil,1,nDesignVariable);
                end
    
                % Adapt design space / initial design
                functionalDiscreteBound = [];
                discreteInitialDesign = [];
                obj.DesignVariableIndexMapping = [];
                obj.BaseVariablesStencils = cell(1,nDesignVariable);
                for i=1:nDesignVariable
                    functionalDiscreteBound = [functionalDiscreteBound,functionalDesignSpace(:,i).*ones(2,nStencil(i))];
                    discreteInitialDesign = [discreteInitialDesign,initialDesign(i)*ones(1,nStencil(i))];

                    obj.DesignVariableIndexMapping = [obj.DesignVariableIndexMapping,i*ones(1,nStencil(i))];
                    obj.BaseVariablesStencils{i} = linspace(baseVariableDesignSpace(1,i),baseVariableDesignSpace(2,i),nStencil(i));
                end
    
                % Adapt component index
                discreteComponentIndex = cell(size(options.ComponentIndex));
                for i = 1:length(discreteComponentIndex)
                    discreteComponentIndex{i} = [];
                    if(~isempty(options.ComponentIndex{i}))
                        for j=1:length(options.ComponentIndex{i})
                            discreteComponentIndex{i} = [discreteComponentIndex{i};find(obj.DesignVariableIndexMapping == options.ComponentIndex{i}(j))'];
                        end
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
                nDesignVariable = max(obj.DesignVariableIndexMapping);
                
                % Create complete design sample with fixed values
                functionalDesignSample = cell(nSample, nDesignVariable);
                for i=1:nSample
                    % For each design variable, create a function handle that interpolates
                    % between the stencil points provided in designSample
                    for j = 1:nDesignVariable
                        % Find all indices in the mapping that correspond to this design variable
                        variableIndices = find(obj.DesignVariableIndexMapping == j);
                        
                        if ~isempty(variableIndices)
                            % Extract the stencil values for this design variable
                            stencilValues = designSample(i, variableIndices);
                            
                            % Create a function handle that uses the interpolation method
                            % to compute values at any parameter point
                            functionalDesignSample{i, j} = @(xInterp) obj.InterpolationMethodFunction(obj.BaseVariablesStencils{j}, stencilValues, xInterp, obj.InterpolationMethodOptions{:});
                        end
                    end
                end
                
                % Call the base mapping's response method
                [performanceMeasure, physicalFeasibilityMeasure, systemOutput] = obj.BaseBottomUpMapping.response(functionalDesignSample);
            end
        end
    end 