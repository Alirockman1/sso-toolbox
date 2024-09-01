function [designOptimal,objectiveOptimal,optimizationOutput] = optimization_ga_wrapper(objectiveFunction,initialDesign,designSpaceLowerBound,designSpaceUpperBound,constraintFunction,varargin)
%OPTIMIZATION_GA_WRAPPER for bottom-up mapping / evaluator optimization
%	OPTIMIZATION_GA_WRAPPER acts as a wrapper for 'ga' (genetic algorithm) when 
%	linear equality and inequality constraints are not being considered. It also
%	allows direct specification of options to be used in 'optimoptions' for
%	'ga'.
%
%	DESIGNOPTIMAL = OPTIMIZATION_GA_WRAPPER(OBJECTIVEFUNCTION,
%	INITIALDESIGN,DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND,
%	CONSTRAINTFUNCTION) minimizes OBJECTIVEFUNCTION starting on INITIALDESIGN,
%	on the design space defined by DESIGNSPACELOWERBOUND and 
%	DESIGNSPACEUPPERBOUND, subject to the constraint functions 
%	CONSTRAINTFUNCTION being less than or equal to zero, returning the optimal
%	design DESIGNOPTIMAL. 
%
%	DESIGNOPTIMAL = OPTIMIZATION_GA_WRAPPER(...NAME,VALUE,...) allows for
%	the specification of additional options for 'ga'. These values get 
%	passed directly to the options function 'optimoptions'.
%
%	[DESIGNOPTIMAL,OBJECTIVEOPTIMAL] = OPTIMIZATION_GA_WRAPPER(...) also
%	returns the value of the objective function OBJECTIVEOPTIMAL for the optimal
%	design.
%
%	[DESIGNOPTIMAL,OBJECTIVEOPTIMAL,OPTIMIZATIONOUTPUT] = 
%	OPTIMIZATION_GA_WRAPPER(...) returns additional information regarding
%	the optimization procedure. In particular:
%			-- ExitCondition : exit flag of 'ga'
%			-- InformationOptimizationProcess : output information of 'ga'
%			-- Population : final population
%			-- GAScores : scores as calculated by 'ga'
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
%			-- Population : (nPopulation,nDesignVariable) double
%			-- GAScores : (nPopulation,1) double
%
%   See also ga, optimoptions, design_optimize_quantities_of_interest, 
%	design_optimize_performance_score.
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

	options = optimoptions('ga',varargin{:});

	nDesignVariable = size(initialDesign,2);
	[designOptimal,objectiveOptimal,exitflag,output,population,scores] = ga(...
		objectiveFunction,...
		nDesignVariable,...
		[],[],[],[],... % A, b, Aeq, beq
		designSpaceLowerBound,...
		designSpaceUpperBound,...
		constraintFunction,...
		[],... % intcon
		options);

	optimizationOutput.ExitCondition = exitflag;
	optimizationOutput.InformationOptimizationProcess = output;
	optimizationOutput.Population = population;
	optimizationOutput.GAScores = scores;
end