function [designOptimal,objectiveOptimal,optimizationOutput] = optimization_patternsearch_wrapper(objectiveFunction,initialDesign,designSpaceLowerBound,designSpaceUpperBound,constraintFunction,varargin)
%OPTIMIZATION_PATTERNSEARCH_WRAPPER for bottom-up mapping / evaluator optimization
%	OPTIMIZATION_PATTERNSEARCH_WRAPPER acts as a wrapper for 'patternsearch' when linear
%	equality and inequality constraints are not being considered. It also
%	allows direct specification of options to be used in 'optimoptions' for
%	'patternsearch'.
%
%	DESIGNOPTIMAL = OPTIMIZATION_PATTERNSEARCH_WRAPPER(OBJECTIVEFUNCTION,
%	INITIALDESIGN,DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND,
%	CONSTRAINTFUNCTION) minimizes OBJECTIVEFUNCTION starting on INITIALDESIGN,
%	on the design space defined by DESIGNSPACELOWERBOUND and 
%	DESIGNSPACEUPPERBOUND, subject to the constraint functions 
%	CONSTRAINTFUNCTION being less than or equal to zero, returning the optimal
%	design DESIGNOPTIMAL. 
%
%	DESIGNOPTIMAL = OPTIMIZATION_PATTERNSEARCH_WRAPPER(...NAME,VALUE,...) allows 
%	for the specification of additional options for 'patternsearch'. These  
%	values get passed directly to the options function 'optimoptions'.
%
%	[DESIGNOPTIMAL,OBJECTIVEOPTIMAL] = OPTIMIZATION_PATTERNSEARCH_WRAPPER(...) 
%	also returns the value of the objective function OBJECTIVEOPTIMAL for the 
%	optimal design.
%
%	[DESIGNOPTIMAL,OBJECTIVEOPTIMAL,OPTIMIZATIONOUTPUT] = 
%	OPTIMIZATION_PATTERNSEARCH_WRAPPER(...) returns additional information 
%	regarding the optimization procedure. In particular:
%			-- ExitCondition : exit flag of 'patternsearch'
%			-- InformationOptimizationProcess : output information of 
%			'patternsearch'
%
%   Input:
%		- OBJECTIVEFUNCTION : function_handle
%		- INITIALDESIGN : (1,nDesignVariable) double
%		- DESIGNSPACELOWERBOUND : (1,nDesignVariable) double
%		- DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double
%		- CONSTRAINTFUNCTION : function_handle
%		- Name-value pair arguments : directly passed to 'optimoptions'.
%
%   Output:
%		- DESIGNOPTIMAL : (1,nDesignVariable) double
%		- OBJECTIVEOPTIMAL : double
%		- OPTIMIZATIONOUTPUT : structure
%			-- ExitCondition : integer
%			-- InformationOptimizationProcess : structure
%
%   See also patternsearch, optimoptions, 
%	design_optimize_quantities_of_interest,  design_optimize_performance_score.
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

	options = optimoptions('patternsearch',varargin{:});

	[designOptimal,objectiveOptimal,exitflag,output] = patternsearch(...
		objectiveFunction,...
		initialDesign,...
		[],[],[],[],... % A, b, Aeq, beq
		designSpaceLowerBound,...
		designSpaceUpperBound,...
		constraintFunction,...
		options);

	optimizationOutput.ExitCondition = exitflag;
	optimizationOutput.InformationOptimizationProcess = output;
end