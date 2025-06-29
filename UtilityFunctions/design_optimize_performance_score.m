function [designOptimal,scoreOptimal,optimizationOutput,evaluationData] = design_optimize_performance_score(designEvaluator,initialDesign,designSpaceLowerBound,designSpaceUpperBound,varargin)
%DESIGN_OPTIMIZE_PERFORMANCE_SCORE Minimize performance deficit while physical
%   DESIGN_OPTIMIZE_PERFORMANCE_SCORE minimizes the performance score for
%   a given design evaluator subject to the constraint that the result must also
%   be physically feasible.
%   To avoid having to perform the same computation of the evaluation twice,
%   nested functions are used. Reference:
%	- https://www.mathworks.com/help/optim/ug/objective-and-nonlinear-constraints-in-the-same-function.html
%
%   DESIGNOPTIMAL = DESIGN_OPTIMIZE_PERFORMANCE_SCORE(DESIGNEVALUATOR,
%   INITIALDESIGN,DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND) optimizes for 
%   DESIGNEVALUATOR starting from the initial design INITIALDESIGN, in the 
%   design space with boundaries DESIGNSPACELOWERBOUND and
%   DESIGNSPACEUPPERBOUND, returning the optimal design DESIGNOPTIMAL.
%
%   DESIGNOPTIMAL = DESIGN_OPTIMIZE_PERFORMANCE_SCORE(...,NAME,VALUE,...) allows
%   for the specification of additional options. These are:
%       - 'OptimizationMethodFunction' : optimization method to be used; see
%       specific implementations. Default: @optimization_fmincon_wrapper.
%       - 'OptimizationMethodOptions' : options for that optimization method.
%       Default: empty.
%       - 'PerformanceDeficitWeight' : weight for each performance deficit 
%       used when computing the score (which serves as objective). Default: 1.
%       - 'PhysicalFeasibilityDeficitWeight' : weight for each physical 
%       feasibility deficit, which just multiplies each value (all are used
%       as cosntraints). Default: 1.
%       - 'CompensationAspaceIndex' : logical index for each design variable
%       to indicate if it belongs to the A-space (in case of a compensation-type
%       problem). Variables with 'true' label here are not varied and are 
%       kept at their value from the initial design. Default is false for all
%       variables (all can be varied to find the optimum).
%
%   [DESIGNOPTIMAL,SCOREOPTIMAL] = DESIGN_OPTIMIZE_PERFORMANCE_SCORE(...) also
%   returns the score of the optimal design SCOREOPTIMAL.
%
%   [DESIGNOPTIMAL,SCOREOPTIMAL,OPTIMIZATIONOUTPUT] = 
%   DESIGN_OPTIMIZE_PERFORMANCE_SCORE(...) also returns any extra outputs from 
%   the optimization method used OPTIMIZATIONOUTPUT.
%
%   [DESIGNOPTIMAL,SCOREOPTIMAL,OPTIMIZATIONOUTPUT,EVALUATIONDATA] = 
%   DESIGN_OPTIMIZE_PERFORMANCE_SCORE(...) additionally returns information
%   regarding the points that were evaluated and their evaluation results in 
%   EVALUATIONDATA. 
%
%   Input:
%       - DESIGNEVALUATOR : DesignEvaluatorBase
%       - INITIALDESIGN : (1,nDesignVariable) double
%       - DESIGNSPACELOWERBOUND : (1,nDesignVariable) double
%       - DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double
%       - 'OptimizationMethodFunction' : function_handle
%       - 'OptimizationMethodOptions' : (1,nOptionOptimization) cell
%       - 'PerformanceDeficitWeight' : (1,nPerformance) double
%       - 'PhysicalFeasibilityDeficitWeight' : (1,nPhysicalFeasibility) double
%       - 'CompensationAspaceIndex' : (1,nDesignVariable) logical
%
%   Output:
%       - DESIGNOPTIMAL : (1,nDesignVariable) double
%       - SCOREOPTIMAL : double
%       - OPTIMIZATIONOUTPUT : structure
%       - EVALUATIONDATA : structure
%           -- DesignSample : (nSample,nDesignVariable) double
%           -- PerformanceDeficit : (nSample,nPerformance) double
%           -- PhysicalFeasibilityDeficit : (nSample,nPhysicalFeasibility) 
%           double
%           -- EvaluatorOutput : class-dependent
%
%   See also design_optimize_quantities_of_interest.
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
    parser.addParameter('PerformanceDeficitWeight',1);
    parser.addParameter('PhysicalFeasibilityDeficitWeight',1);
    parser.addParameter('CompensationAspaceIndex',[]);

    parser.parse(varargin{:});
    options = parser.Results;

    nDesignVariable = size(initialDesign,2);
    if(isempty(options.CompensationAspaceIndex))
        options.CompensationAspaceIndex = false(1,nDesignVariable);
    end

    isOutputEvaluationData = (nargout>3);
    if(isOutputEvaluationData)
        evaluationData = struct(...
            'DesignSample',[],...
            'PerformanceDeficit',[],...
            'PhysicalFeasibilityDeficit',[],...
            'EvaluatorOutput',[]);
        iEvaluation = 1;
    end

	designLastBspace = [];
	performanceDeficit = [];
	physicalFeasibilityDeficit = [];

    initialDesignAspace = initialDesign(options.CompensationAspaceIndex);
    initialDesignBspace = initialDesign(~options.CompensationAspaceIndex);

    if(size(designSpaceLowerBound,2)==nDesignVariable)
        bspaceLowerBound = designSpaceLowerBound(~options.CompensationAspaceIndex);
    else
        bspaceLowerBound = designSpaceLowerBound;
    end
    if(size(designSpaceUpperBound,2)==nDesignVariable)
        bspaceUpperBound = designSpaceUpperBound(~options.CompensationAspaceIndex);
    else
        bspaceUpperBound = designSpaceUpperBound;
    end

	objectiveFunction = @objfun;
	constraintFunction = @constr;

	[designOptimalBspace,scoreOptimal,optimizationOutput] = options.OptimizationMethodFunction(...
        objectiveFunction,...
        initialDesignBspace,...
        bspaceLowerBound,...
        bspaceUpperBound,...
        constraintFunction,...
        options.OptimizationMethodOptions{:});

    designOptimal = nan(1,nDesignVariable);
    designOptimal(options.CompensationAspaceIndex) = initialDesignAspace;
    designOptimal(~options.CompensationAspaceIndex) = designOptimalBspace;

	% nested objective function
	function objective = objfun(designCurrentBspace)
        if ~isequal(designCurrentBspace,designLastBspace) % Check if computation is necessary
            nSample = size(designCurrentBspace,1);
            designCurrent = nan(nSample,nDesignVariable);
            designCurrent(:,options.CompensationAspaceIndex) = initialDesignAspace;
            designCurrent(:,~options.CompensationAspaceIndex) = designCurrentBspace;

            [performanceDeficit,physicalFeasibilityDeficit,evaluatorOutput] = designEvaluator.evaluate(designCurrent);
            designLastBspace = designCurrentBspace;

            if(isOutputEvaluationData)
                evaluationData(iEvaluation) = struct(...
                    'DesignSample',designCurrent,...
                    'PerformanceDeficit',performanceDeficit,...
                    'PhysicalFeasibilityDeficit',physicalFeasibilityDeficit,...
                    'EvaluatorOutput',evaluatorOutput);
                iEvaluation = iEvaluation+1;
            end
        end
        % Now compute objective function
        [~,objective] = design_deficit_to_label_score(performanceDeficit,options.PerformanceDeficitWeight);
    end

    % nested constraint function
    function [constraintInequality,constraintEquality] = constr(designCurrentBspace)
        if ~isequal(designCurrentBspace,designLastBspace) % Check if computation is necessary
            nSample = size(designCurrentBspace,1);
            designCurrent = nan(nSample,nDesignVariable);
            designCurrent(:,options.CompensationAspaceIndex) = initialDesignAspace;
            designCurrent(:,~options.CompensationAspaceIndex) = designCurrentBspace;

            [performanceDeficit,physicalFeasibilityDeficit,evaluatorOutput] = designEvaluator.evaluate(designCurrent);
            designLastBspace = designCurrentBspace;

            if(isOutputEvaluationData)
                evaluationData(iEvaluation) = struct(...
                    'DesignSample',designCurrent,...
                    'PerformanceDeficit',performanceDeficit,...
                    'PhysicalFeasibilityDeficit',physicalFeasibilityDeficit,...
                    'EvaluatorOutput',evaluatorOutput);
                iEvaluation = iEvaluation+1;
            end
        end
        % Now compute constraint function
        constraintInequality = options.PhysicalFeasibilityDeficitWeight.*physicalFeasibilityDeficit;
        constraintEquality = [];
    end
end