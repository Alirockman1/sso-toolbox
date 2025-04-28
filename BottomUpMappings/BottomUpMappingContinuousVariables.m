classdef BottomUpMappingContinuousVariables < BottomUpMappingBase
    %BOTTOMUPMAPPINGCONTINUOUSVARIABLES Bottom-Up Mapping with continuous design variables
    %   BOTTOMUPMAPPINGCONTINUOUSVARIABLES takes a base bottom-up mapping and converts
    %   discrete design variables into continuous functions through interpolation.
    %   It manages the conversion between discrete points and continuous functions,
    %   handles fixed points in the design space, and provides a consistent interface
    %   for optimization with continuous variables.
    %
    %   This class is derived from BottomUpMappingBase.
    %
    %   BOTTOMUPMAPPINGCONTINUOUSVARIABLES methods:
    %       - response : receives design sample points as function handles and returns
    %       the system response in terms of performance and physical feasibility
    %
    %   See also BottomUpMappingBase, BottomUpMappingFunction.
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
            %   BASEBOTTOMUPMAPPING is the underlying bottom-up mapping that will be 
            %   called after converting discrete design points into continuous functions
            %   through interpolation.
            %
            %   BASEBOTTOMUPMAPPING : BottomUpMappingBase
            %
            %   See also InterpolationMethodFunction, InterpolationMethodOptions.
            BaseBottomUpMapping
    
            %INTERPOLATIONMETHODFUNCTION Function handle for interpolation
            %   INTERPOLATIONMETHODFUNCTION is the function used to interpolate between
            %   discrete design points to create continuous functions. By default, this
            %   is MATLAB's interp1 function.
            %
            %   INTERPOLATIONMETHODFUNCTION : function_handle
            %
            %   See also InterpolationMethodOptions.
            InterpolationMethodFunction
    
            %INTERPOLATIONMETHODOPTIONS Options for interpolation method
            %   INTERPOLATIONMETHODOPTIONS contains the options passed to the
            %   interpolation function. By default, these are {'pchip','extrap'} for
            %   piecewise cubic hermite interpolation with extrapolation.
            %
            %   INTERPOLATIONMETHODOPTIONS : cell
            %
            %   See also InterpolationMethodFunction.
            InterpolationMethodOptions

            %DESIGNVARIABLEINDEXMAPPING Mapping between design variables and indices
            %   DESIGNVARIABLEINDEXMAPPING maps each discrete point to its corresponding
            %   design variable index. This is used to reconstruct continuous functions
            %   from discrete points during interpolation.
            %
            %   DESIGNVARIABLEINDEXMAPPING : (1,nDiscretizedPoints) double
            %
            %   See also BaseVariablesStencils.
            DesignVariableIndexMapping

            %BASEVARIABLESSTENCILS Stencil points for base variables
            %   BASEVARIABLESSTENCILS contains the discretization points (stencils) for
            %   each base variable where the interpolation is performed.
            %
            %   BASEVARIABLESSTENCILS : (1,nDesignVariable) cell
            %
            %   See also DesignVariableIndexMapping.
            BaseVariablesStencils

            %FIXEDPOINTINDEX Indices of fixed points in design space
            %   FIXEDPOINTINDEX contains the indices of points in the design space that
            %   are fixed to specific values and should not be interpolated.
            %
            %   FIXEDPOINTINDEX : (1,nDesignVariable) cell
            %
            %   See also FixedPointValue.
            FixedPointIndex 

            %FIXEDPOINTVALUE Values at fixed points
            %   FIXEDPOINTVALUE contains the values that should be used at fixed points
            %   in the design space instead of interpolated values.
            %
            %   FIXEDPOINTVALUE : (1,nDesignVariable) cell
            %
            %   See also FixedPointIndex.
            FixedPointValue
        end
    
        methods
            function [obj,functionalDiscreteBound,discreteInitialDesign,discreteComponentIndex] = BottomUpMappingContinuousVariables(baseBottomUpMapping, baseVariableDesignSpace, functionalDesignSpace, initialDesign, varargin)
            %BOTTOMUPMAPPINGCONTINUOUSVARIABLES Constructor
            %   BOTTOMUPMAPPINGCONTINUOUSVARIABLES creates an object that converts discrete
            %   design variables into continuous functions through interpolation.
            %
            %   OBJ = BOTTOMUPMAPPINGCONTINUOUSVARIABLES(BASEBOTTOMUPMAPPING,
            %   BASEVARIABLEDESIGNSPACE, FUNCTIONALDESIGNSPACE, INITIALDESIGN) creates a
            %   new object that uses BASEBOTTOMUPMAPPING as the underlying mapping.
            %   BASEVARIABLEDESIGNSPACE defines the domain for interpolation,
            %   FUNCTIONALDESIGNSPACE defines the range, and INITIALDESIGN provides
            %   initial values.
            %
            %   [OBJ,FUNCTIONALDISCRETEBOUND,DISCRETEINITIALDESIGN,DISCRETECOMPONENTINDEX] =
            %   BOTTOMUPMAPPINGCONTINUOUSVARIABLES(...) additionally returns the discretized
            %   bounds, initial design, and component indices adapted for the continuous
            %   representation.
            %
            %   [...] = BOTTOMUPMAPPINGCONTINUOUSVARIABLES(...,NAME,VALUE) specifies
            %   additional options using name-value pairs:
            %       - 'InterpolationMethodFunction' : function handle for interpolation
            %         (default: @interp1)
            %       - 'InterpolationMethodOptions' : cell array of options for interpolation
            %         (default: {'pchip','extrap'})
            %       - 'NumberOfStencils' : number of discretization points
            %         (default: 1)
            %       - 'FixedPointsBaseVariables' : points with fixed values
            %         (default: [])
            %       - 'FixedPointsFunctionalValues' : values at fixed points
            %         (default: [])
            %       - 'ComponentIndex' : indices for component grouping
            %         (default: {})
            %
            %   Inputs:
            %       - BASEBOTTOMUPMAPPING : BottomUpMappingBase
            %       - BASEVARIABLEDESIGNSPACE : (2,nDesignVariable) double
            %       - FUNCTIONALDESIGNSPACE : (2,nDesignVariable) double
            %       - INITIALDESIGN : (1,nDesignVariable) double
            %   
            %   Outputs:
            %       - OBJ : BottomUpMappingContinuousVariables
            %       - FUNCTIONALDISCRETEBOUND : (2,nDiscretizedPoints) double
            %       - DISCRETEINITIALDESIGN : (1,nDiscretizedPoints) double
            %       - DISCRETECOMPONENTINDEX : (1,nComponent) cell
            %
            %   See also response, interp1.
    
                % Parse input arguments
                parser = inputParser;
                parser.addOptional('ComponentIndex', {});
                parser.addParameter('InterpolationMethodFunction',@interp1);
                parser.addParameter('InterpolationMethodOptions',{'pchip','extrap'});
                parser.addParameter('NumberOfStencils',1);
                parser.addParameter('FixedPointsBaseVariables',[]);
                parser.addParameter('FixedPointsFunctionalValues',[]);
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

                % Adapt fixed points
                if(~isempty(options.FixedPointsBaseVariables))
                    obj.FixedPointIndex = cell(1,nDesignVariable);
                    obj.FixedPointValue = cell(1,nDesignVariable);
                    for i=1:nDesignVariable
                        isFixed = ~isnan(options.FixedPointsBaseVariables(:,i));
                        nFixed = sum(isFixed);

                        [obj.BaseVariablesStencils{i},sortOrder] = unique([obj.BaseVariablesStencils{i},options.FixedPointsBaseVariables(isFixed,i)'],'sorted');
                        originalOrder = 1:length(obj.BaseVariablesStencils{i});

                        [~,obj.FixedPointIndex{i}] = ismember(originalOrder(originalOrder(end-nFixed+1:end)),sortOrder(:)');
                        obj.FixedPointValue{i} = options.FixedPointsFunctionalValues(isFixed,i)';
                    end
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
            %   RESPONSE converts discrete design points into continuous functions through
            %   interpolation before passing them to the base mapping's response method.
            %
            %   PERFORMANCEMEASURE = OBJ.RESPONSE(DESIGNSAMPLE) returns the performance
            %   measure for each design in PERFORMANCEMEASURE. Each row corresponds to a
            %   design sample point.
            %
            %   [PERFORMANCEMEASURE,PHYSICALFEASIBILITYMEASURE] = OBJ.RESPONSE(...)
            %   additionally returns the physical feasibility measure for each design in
            %   PHYSICALFEASIBILITYMEASURE.
            %
            %   [PERFORMANCEMEASURE,PHYSICALFEASIBILITYMEASURE,SYSTEMOUTPUT] =
            %   OBJ.RESPONSE(...) additionally returns any extra information from the
            %   system response computation in SYSTEMOUTPUT.
            %
            %   Inputs:
            %       - OBJ : BottomUpMappingContinuousVariables
            %       - DESIGNSAMPLE : (nSample,nDiscretizedPoints) double
            %   
            %   Outputs:
            %       - PERFORMANCEMEASURE : (nSample,nPerformance) double
            %       - PHYSICALFEASIBILITYMEASURE : (nSample,nPhysicalFeasibility) double
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
                            stencilValuesSample = designSample(i, variableIndices);

                            % If the design variable is a fixed point, use the fixed point value
                            if(~isempty(obj.FixedPointIndex) && ~isempty(obj.FixedPointIndex{j}))
                                stencilValues = nan(1,length(obj.BaseVariablesStencils{j}));
                                stencilValues(obj.FixedPointIndex{j}) = obj.FixedPointValue{j};
                                stencilValues(isnan(stencilValues)) = stencilValuesSample;
                            else
                                stencilValues = stencilValuesSample;
                            end
                            
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