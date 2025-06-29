function [designOptimal,objectiveOptimal,optimizationOutput,systemResponseData] = design_optimize_quantities_of_interest(bottomUpMapping,initialDesign,designSpaceLowerBound,designSpaceUpperBound,objectiveFunction,varargin)
%DESIGN_OPTIMIZE_QUANTITIES_OF_INTEREST Optimization with bottom-up mapping
%   DESIGN_OPTIMIZE_QUANTITIES_OF_INTEREST allows for easier optimization of
%   the results of a bottom-up mapping response. The objective and constraint 
%   functions (both equality and inequality) can be specified and derived 
%   directly from the performance (and physical feasibility) measures.
%   To avoid having to perform the same computation of the response twice, 
%   nested functions are used. Reference:
%	- https://www.mathworks.com/help/optim/ug/objective-and-nonlinear-constraints-in-the-same-function.html
%
%   DESIGNOPTIMAL = DESIGN_OPTIMIZE_QUANTITIES_OF_INTEREST(BOTTOMUPMAPPING,
%   INITIALDESIGN,DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND,OBJECTIVEFUNCTION)
%   minimizes the value of OBJECTIVEFUNCTION starting from the initial design
%   INITIALDESIGN, in the design space with boundaries DESIGNSPACELOWERBOUND and
%   DESIGNSPACEUPPERBOUND, returning the optimal design DESIGNOPTIMAL. The 
%   values for the objective function are computed from the response of 
%   BOTTOMUPMAPPING. The objective function must have one of the following two
%   interfaces:
%       - objectiveValue = f(performanceMeasure)
%       - objectiveValue = f(performanceMeasure,physicalFeasibilityMeasure,...,
%       NAME,VALUE,...)
%
%   DESIGNOPTIMAL = DESIGN_OPTIMIZE_QUANTITIES_OF_INTEREST(...,NAME,VALUE,...) 
%   allows for the specification of additional options. These are:
%       - 'OptimizationMethodFunction' : optimization method to be used; see
%       specific implementations. Default: @optimization_fmincon_wrapper.
%       - 'OptimizationMethodOptions' : options for that optimization method.
%       Default: empty.
%       - 'ObjectiveOptions' : options to be passed on to the objective
%       function. Default: empty.
%       - 'InequalityConstraintFunction' : function to determine the inequality 
%       constraint values. Has same interface as OBJECTIVEFUNCTION. Default is
%       empty.
%       - 'InequalityConstraintOptions' : options to be passed on to the 
%       inequality constraint function. Default: empty.
%       - 'EqualityConstraintFunction' : function to determine the equality 
%       constraint values. Has same interface as OBJECTIVEFUNCTION. Default is
%       empty.
%       - 'EqualityConstraintOptions' : options to be passed on to the 
%       equality constraint function. Default: empty.
%       - 'IsFixedDesignVariable' : logical index for each design variable
%       to indicate if it should remain constant during the optimization. 
%       Variables with 'true' label here are not varied and are kept at their 
%       value from the initial design. Default is false for all variables (all 
%       can be varied to find the optimum).
%
%   [DESIGNOPTIMAL,OBJECTIVEOPTIMAL] = 
%   DESIGN_OPTIMIZE_QUANTITIES_OF_INTEREST(...) also returns the objective 
%   function value of the optimal design OBJECTIVEOPTIMAL.
%
%   [DESIGNOPTIMAL,OBJECTIVEOPTIMAL,OPTIMIZATIONOUTPUT] = 
%   DESIGN_OPTIMIZE_QUANTITIES_OF_INTEREST(...) also returns any extra outputs  
%   from the optimization methoo used OPTIMIZATIONOUTPUT.
%
%   [DESIGNOPTIMAL,OBJECTIVEOPTIMAL,OPTIMIZATIONOUTPUT,SYSTEMRESPONSEDATA] = 
%   DESIGN_OPTIMIZE_QUANTITIES_OF_INTEREST(...) additionally returns information
%   regarding the points that were evaluated and the response results in 
%   SYSTEMRESPONSEDATA. 
%
%   See also design_optimize_performance_score.
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

    parser = inputParser;
    parser.addParameter('OptimizationMethodFunction',@optimization_fmincon_wrapper);
    parser.addParameter('OptimizationMethodOptions',{});
    parser.addParameter('ObjectiveOptions',{});
    parser.addParameter('InequalityConstraintFunction',[]);
    parser.addParameter('InequalityConstraintOptions',{});
    parser.addParameter('EqualityConstraintFunction',[]);
    parser.addParameter('EqualityConstraintOptions',{});
    parser.addParameter('IsFixedDesignVariable',[]);

    parser.parse(varargin{:});
    options = parser.Results;

    nDesignVariable = size(initialDesign,2);
    if(isempty(options.IsFixedDesignVariable))
        options.IsFixedDesignVariable = false(1,nDesignVariable);
    end

    isOutputSystemResponseData = (nargout>3);
    if(isOutputSystemResponseData)
        systemResponseData = struct(...
            'DesignSample',[],...
            'PerformanceMeasure',[],...
            'PhysicalFeasibilityMeasure',[],...
            'BottomUpMappingOutput',[]);
        iEvaluation = 1;
    end

	designLastVariable = [];
	performanceMeasure = [];
	physicalFeasibilityMeasure = [];

    initialDesignFixed = initialDesign(options.IsFixedDesignVariable);
    initialDesignVariable = initialDesign(~options.IsFixedDesignVariable);

    if(size(designSpaceLowerBound,2)==nDesignVariable)
        designSpaceLowerBoundVariable = designSpaceLowerBound(~options.IsFixedDesignVariable);
    else
        designSpaceLowerBoundVariable = designSpaceLowerBound;
    end
    if(size(designSpaceUpperBound,2)==nDesignVariable)
        deisgnSpaceUpperBoundVariable = designSpaceUpperBound(~options.IsFixedDesignVariable);
    else
        deisgnSpaceUpperBoundVariable = designSpaceUpperBound;
    end

	[designOptimalVariable,objectiveOptimal,optimizationOutput] = options.OptimizationMethodFunction(...
        @objfun,...
        initialDesignVariable,...
        designSpaceLowerBoundVariable,...
        deisgnSpaceUpperBoundVariable,...
        @constr,...
        options.OptimizationMethodOptions{:});

    designOptimal = nan(1,nDesignVariable);
    designOptimal(options.IsFixedDesignVariable) = initialDesignFixed;
    designOptimal(~options.IsFixedDesignVariable) = designOptimalVariable;
    return;

	% nested objective function
	function objective = objfun(designCurrentVariable)
        if ~isequal(designCurrentVariable,designLastVariable) % Check if computation is necessary
            nSample = size(designCurrentVariable,1);
            designCurrent = nan(nSample,nDesignVariable);
            designCurrent(:,options.IsFixedDesignVariable) = initialDesignFixed;
            designCurrent(:,~options.IsFixedDesignVariable) = designCurrentVariable;

            [performanceMeasure,physicalFeasibilityMeasure,bottomUpMappingOutput] = bottomUpMapping.response(designCurrent);
            designLastVariable = designCurrentVariable;

            if(isOutputSystemResponseData)
                systemResponseData(iEvaluation) = struct(...
                    'DesignSample',designCurrent,...
                    'PerformanceMeasure',performanceMeasure,...
                    'PhysicalFeasibilityMeasure',physicalFeasibilityMeasure,...
                    'BottomUpMappingOutput',bottomUpMappingOutput);
                iEvaluation = iEvaluation+1;
            end
        end
        % Now compute objective function
        if(nargin(objectiveFunction)==1)
            objective = objectiveFunction(performanceMeasure);
        else
            objective = objectiveFunction(performanceMeasure,physicalFeasibilityMeasure,options.ObjectiveOptions{:});
        end
    end

    % nested constraint function
    function [constraintInequality,constraintEquality] = constr(designCurrentVariable)
        if ~isequal(designCurrentVariable,designLastVariable) % Check if computation is necessary
            nSample = size(designCurrentVariable,1);
            designCurrent = nan(nSample,nDesignVariable);
            designCurrent(:,options.IsFixedDesignVariable) = initialDesignFixed;
            designCurrent(:,~options.IsFixedDesignVariable) = designCurrentVariable;

            [performanceMeasure,physicalFeasibilityMeasure,bottomUpMappingOutput] = bottomUpMapping.response(designCurrent);
            designLastVariable = designCurrentVariable;

            if(isOutputSystemResponseData)
                systemResponseData(iEvaluation) = struct(...
                    'DesignSample',designCurrent,...
                    'PerformanceMeasure',performanceMeasure,...
                    'PhysicalFeasibilityMeasure',physicalFeasibilityMeasure,...
                    'BottomUpMappingOutput',bottomUpMappingOutput);
                iEvaluation = iEvaluation+1;
            end
        end
        % Now compute constraint function
        if(~isempty(options.InequalityConstraintFunction))
            if(nargin(options.InequalityConstraintFunction)==1)
                constraintInequality = options.InequalityConstraintFunction(performanceMeasure);
            else
                constraintInequality = options.InequalityConstraintFunction(performanceMeasure,...
                    physicalFeasibilityMeasure,options.InequalityConstraintOptions{:});
            end
        else
            constraintInequality = [];
        end
        if(~isempty(options.EqualityConstraintFunction))
            if(nargin(options.InequalityConstraintFunction)==1)
                constraintEquality = options.EqualityConstraintFunction(performanceMeasure);
            else
                constraintEquality = options.EqualityConstraintFunction(performanceMeasure,...
                    physicalFeasibilityMeasure,options.EqualityConstraintOptions{:});
            end
        else
            constraintEquality = [];
        end
    end
end