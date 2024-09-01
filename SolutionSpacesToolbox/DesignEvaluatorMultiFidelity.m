classdef DesignEvaluatorMultiFidelity < DesignEvaluatorBase
%DESIGNEVALUATORMULTIFIDELITY Design Evaluator that mixes high- and low-fidelity
%   DESIGNEVALUATORMULTIFIDELITY allows one to use two different design 
%   evaluators, one which is high-fidelity and one which is low-fidelity. In 
%   general, the low-fidelity one is used first on all requested sample points, 
%   and for the designs where it is deemed that there is uncertainty, the high-
%   fidelity model is then used. This is done to minimize the function calls
%   to the high-fidelity model (which is assumed to be more expensive).
%
%   This class is derived from DesignEvaluatorBase.
%
%   DESIGNEVALUATORMULTIFIDELITY methods:
%       - evaluate : a function that receives the design samples and returns 
%       their deficit in terms of performance and physical feasibility relative
%       to the requirements.
%
%   See also DesignEvaluatorBottomUpMapping, DesignEvaluatorFastForward.
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
        %DESIGNEVALUATORLOWFIDELITY Low-fidelity model
        %   DESIGNEVALUATORLOWFIDELITY is the low-fidelity design evaluator. It is used
        %   every time for every sample point that has to be evaluated. It is based on
        %   its output (and how uncertain it is) that the high-fidelity model may or
        %   may not be called.
        %
        %   DESIGNEVALUATORLOWFIDELITY : DesignEvaluatorBase
        %   
        %   See also DesignEvaluatorHighFidelity.
        DesignEvaluatorLowFidelity

        %DESIGNEVALUATORHIGHFIDELITY High-fidelity model
        %   DESIGNEVALUATORHIGHFIDELITY is the high-fidelity design evaluator. It is 
        %   only used when the output from the low-fidelity model is deemed uncertain.
        %
        %   DESIGNEVALUATORHIGHFIDELITY : DesignEvaluatorBase
        %   
        %   See also DesignEvaluatorLowFidelity.
        DesignEvaluatorHighFidelity

        %DECISIONFUNCTION Decide if the output from low-fidelity model is uncertain
        %   DECISIONFUNCTION is a function handle for the function used to determine 
        %   which sample points that were evaluated using the low-fidelity model have
        %   an uncertain output (meaning, it's not so clear whether the design is 
        %   actually good/bad or physically feasible/infeasible).
        %
        %   DECISIONFUNCTION : function_handle
        %   
        %   See also DecisionOptions.
        DecisionFunction

        %DECISIONOPTIONS Extra options for the uncertainty decision function
        %   DECISIONOPTIONS are any extra options necessary for the uncertainty decision
        %   function to be used.
        %
        %   DECISIONOPTIONS : (1,nOption) cell
        %   
        %   See also DecisionFunction.
        DecisionOptions
    end
    
    methods
        function obj = DesignEvaluatorMultiFidelity(designEvaluatorLowFidelity,designEvaluatorHighFidelity,decisionFunction,varargin)
        %DESIGNEVALUATORMULTIFIDELITY Constructor
        %   DESIGNEVALUATORMULTIFIDELITY uses two different design evaluators to setup
        %   a DesignEvaluatorMultiFidelity object.
        %
        %   OBJ = DESIGNEVALUATORMULTIFIDELITY(DESIGNEVALUATORLOWFIDELITY,
        %   DESIGNEVALUATORHIGHFIDELITY,DECISIONFUNCTION) receives a low-fidelity 
        %   design evaluator in DESIGNEVALUATORLOWFIDELITY, a high-fidelity design 
        %   evaluator DESIGNEVALUATORHIGHFIDELITY, and an uncertainty decision function
        %   in DECISIONFUNCTION, which are used to create a DesignEvaluatorMultiFidelity
        %   object OBJ. The outputs of the low-fidelity design evaluator are fed into
        %   the decision function, which must determine which ones are uncertain. For
        %   those selected as uncertain, the high-fidelity model is then used.
        %
        %   OBJ = DESIGNEVALUATORMULTIFIDELITY(DESIGNEVALUATORLOWFIDELITY,
        %   DESIGNEVALUATORHIGHFIDELITY,DECISIONFUNCTION,DECISIONOPTIONS) also allows 
        %   the specification of additional options for the decision function. Default 
        %   value is empty.
        %
        %   Input:
        %       - DESIGNEVALUATORLOWFIDELITY : DesignEvaluatorBase
        %       - DESIGNEVALUATORHIGHFIDELITY : DesignEvaluatorBase
        %       - DECISIONFUNCTION : function_handle
        %       - DECISIONOPTIONS : (1,nOption) cell
        %   
        %   Output:
        %       - OBJ : DesignEvaluatorMultiFidelity
        %   
        %   See also evaluate.

            parser = inputParser;
            parser.addOptional('DecisionOptions',{});
            parser.parse(varargin{:});

            obj.DesignEvaluatorLowFidelity = designEvaluatorLowFidelity;
            obj.DesignEvaluatorHighFidelity = designEvaluatorHighFidelity;
            obj.DecisionFunction = decisionFunction;
            obj.DecisionOptions = parser.Results.DecisionOptions;
        end
        
        function [performanceDeficit,physicalFeasibilityDeficit,evaluationOutput] = evaluate(obj,designSample)
        %EVALUATE Evaluation of designs' performance relative to requirements
        %   EVALUATE is the main function of any design evaluator, where design sample
        %   points are given and judgement must be made as to whether or not that design
        %   meets the system requirements.
        %   Here, the low-fidelity model is used first on all samples, and for those 
        %   where the result is considered uncertain, the high-fidelity model is then
        %   used.
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
        %   process in EVALUATIONOUTPUT; in this case, that is the extra output of both 
        %   the low- and high-fidelity models, as well as the logical array indicating
        %   which samples were uncertain.
        %
        %   Input:
        %       - OBJ : DesignEvaluatorBase
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %   
        %   Output:
        %       - PERFORMANCEDEFICIT : (nSample,nPerformance) double
        %       - PHYSICALFEASIBILITYDEFICIT : (nSample,nPhysicalFeasibility) double
        %       - EVALUATIONOUTPUT : structure
        %           -- IsUncertain : (nSample,1) logical
        %           -- EvaluatorOutputLowFidelity : class-dependent
        %           -- EvaluatorOutputHighFidelity : class-dependent
        %
        %   See also DesignEvaluatorBase, DesignEvaluatorMultiFidelity.
            
            nSample = size(designSample,1);

            [performanceDeficitLowFidelity,physicalFeasibilityDeficitLowFidelity,evaluatorOutputLowFidelity] = ...
                obj.DesignEvaluatorLowFidelity.evaluate(designSample);

            isUncertain = obj.DecisionFunction(...
                designSample,...
                performanceDeficitLowFidelity,...
                physicalFeasibilityDeficitLowFidelity,...
                evaluatorOutputLowFidelity,...
                obj.DecisionOptions{:});

            if(any(isUncertain))
                [performanceDeficitHighFidelity,physicalFeasibilityDeficitHighFidelity,evaluatorOutputHighFidelity] = ...
                    obj.DesignEvaluatorHighFidelity.evaluate(designSample(isUncertain,:));
            else
                performanceDeficitHighFidelity = [];
                physicalFeasibilityDeficitHighFidelity = [];
                evaluatorOutputHighFidelity = [];
            end

            % performance
            nPerformanceDeficitLowFidelity = size(performanceDeficitLowFidelity,2);
            nPerformanceDeficitHighFidelity = size(performanceDeficitHighFidelity,2);
            nPerformanceDeficit = max(nPerformanceDeficitLowFidelity,nPerformanceDeficitHighFidelity);

            performanceDeficit = nan(nSample,nPerformanceDeficit);
            performanceDeficit(~isUncertain,1:nPerformanceDeficitLowFidelity) = performanceDeficitLowFidelity(~isUncertain,:);
            performanceDeficit(isUncertain,1:nPerformanceDeficitHighFidelity) = performanceDeficitHighFidelity;

            % physical feasibility
            nPhysicalFeasibilityDeficitLowFidelity = size(physicalFeasibilityDeficitLowFidelity,2);
            nPhysicalFeasibilityDeficitHighFidelity = size(physicalFeasibilityDeficitHighFidelity,2);
            nPhysicalFeasibilityDeficit = max(nPhysicalFeasibilityDeficitLowFidelity,nPhysicalFeasibilityDeficitHighFidelity);

            physicalFeasibilityDeficit = nan(nSample,nPhysicalFeasibilityDeficit);
            physicalFeasibilityDeficit(~isUncertain,1:nPhysicalFeasibilityDeficitLowFidelity) = physicalFeasibilityDeficitLowFidelity(~isUncertain,:);
            physicalFeasibilityDeficit(isUncertain,1:nPhysicalFeasibilityDeficitHighFidelity) = physicalFeasibilityDeficitHighFidelity;

            % wrap extra outputs
            evaluationOutput.IsUncertain = isUncertain;
            evaluationOutput.EvaluatorOutputLowFidelity = evaluatorOutputLowFidelity;
            evaluationOutput.EvaluatorOutputHighFidelity = evaluatorOutputHighFidelity;
        end
    end
end