classdef (Abstract) DesignEvaluatorBase
%DESIGNEVALUATORBASE Design Evaluator abstract class for top-down design
%   DESIGNEVALUATORBASE provides the interface for particular implementations  
%   that may serve as evaluators when doing top-down development/design.
%   The goal of this class and its derivations is to receive as input possible
%   design sample points and return if those designs meet the performance 
%   requirements or not.
%   Additionally, if requirement spaces are being considered, also if those
%   designs are physically feasible or not.
%
%   As an abstract class, particular implementations have to be created using 
%   this as base for one to be able to create objects.
%
%   DESIGNEVALUATORBASE methods:
%       - evaluate : a function that receives the design sample points and  
%       returns their deficit in terms of performance and physical feasibility 
%       relative to the requirements.
%
%   See also DesignEvaluatorBottomUpMapping, DesignEvaluatorCompensation.
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

    methods (Abstract)
        %EVALUATE Evaluation of designs' performance relative to requirements
        %   EVALUATE is the main function of any design evaluator, where design sample
        %   points are given and judgement must be made as to whether or not that design
        %   meets the system requirements.
        %
        %   PERFORMANCEDEFICIT = OBJ.EVALUATE(DESIGNSAMPLE) returns the performance
        %   deficit PERFORMANCEDEFICIT of each given design sample point DESIGNSAMPLE.
        %   This is given as an array where negative values mean that performance metric
        %   meets the requirements, and positive values mean that performance metric 
        %   does not meet the requirements.
        %
        %   [PERFORMANCEDEFICIT,PHYSICALFEASIBILITYDEFICIT] = OBJ.EVALUATE(...) also
        %   returns the physical feasibility deficit of each given design in
        %   PHYSICALFEASIBILITYDEFICIT; similar to PERFORMANCEDEFICIT, negative values 
        %   means the specific component is physically feasible (within the specified 
        %   physical feasibility requirements), and positive values means that component
        %   is physically infeasible.
        %
        %   [PERFORMANCEDEFICIT,PHYSICALFEASIBILITYDEFICIT,EVALUATIONOUTPUT] = 
        %   OBJ.EVALUATE(...) also returns any more information from the evaluation
        %   process in EVALUATIONOUTPUT; what that is and its type is are entirely 
        %   dependent on the specific class implementation.
        %
        %   Inputs:
        %       - OBJ : DesignEvaluatorBase
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %   
        %   Outputs:
        %       - PERFORMANCEDEFICIT : (nSample,nPerformance) double
        %       - PHYSICALFEASIBILITYDEFICIT : (nSample,nPhysicalFeasibility) double
        %       - EVALUATIONOUTPUT : class-defined
        %
        %   See also design_measure_to_deficit.
        [performanceDeficit,physicalFeasibilityDeficit,evaluationOutput] = evaluate(obj,designSample)
    end
end