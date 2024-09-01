classdef DesignEvaluatorFastForward < DesignEvaluatorBase
%DESIGNEVALUATORFASTFORWARD Design Evaluator that uses Fast Forward models
%   DESIGNEVALUATORFASTFORWARD takes design fast-forward models and adapts the 
%   interface calls so they may also work as design evaluators.
%
%   This class is derived from DesignEvaluatorBase.
%
%   DESIGNEVALUATORFASTFORWARD methods:
%       - evaluate : a function that receives the design samples and returns 
%       their deficit in terms of performance and physical feasibility relative
%       to the requirements.
%
%   See also DesignFastForwardBase.
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
        %MODELPERFORMANCE Design Fast-forward Model used to predict performance
        %   MODELPERFORMANCE is the design fast-forward model used to predict
        %   whether a given design meets the requirements or not.
        %
        %   MODELPERFORMANCE : DesignFastForward object
        %   
        %   See also DesignFastForward, evaluate.
        modelPerformance
        
        %MODELPHYSICALFEASIBILITY Model used to predict physical feasibility
        %   MODELPHYSICALFEASIBILITY is the design fast-forward model used to 
        %   predict whether a given design is physically feasible or not.
        %
        %   MODELPHYSICALFEASIBILITY : DesignFastForward object
        %   
        %   See also DesignFastForward, evaluate.
        modelPhysicalFeasibility
    end
    
    methods
        function obj = DesignEvaluatorFastForward(modelPerformance,modelPhysicalFeasibility)
        %DESIGNEVALUATORFASTFORWARD Constructor
        %   DESIGNEVALUATORFASTFORWARD uses design fast-forward models to 
        %   create a design evaluator.
        %
        %   OBJ = DESIGNEVALUATORFASTFORWARD(MODELPERFORMANCE) uses the design
        %   fast-forward model specified in MODELPERFORMANCE to create a new object of
        %   the DesignEvaluatorFastForward class in OBJ. Said object uses the 
        %   fast-forward model to determine if a given design is good or not when the
        %   'evaluate' method is called.
        %
        %   OBJ = DESIGNEVALUATORFASTFORWARD(MODELPERFORMANCE,MODELPHYSICALFEASIBILITY) 
        %   also takes a second fast-forward model given in MODELPHYSICALFEASIBILITY and
        %   uses it to determine the physical feasibility of designs when the 'evaluate'
        %   method is called. 
        %
        %   Input:
        %       - MODELPERFORMANCE : DesignFastForwardBase
        %       - MODELPHYSICALFEASIBILITY : DesignFastForwardBase
        %   
        %   Output:
        %       - OBJ : DesignEvaluatorFastForward
        %   
        %   See also evaluate.

            if(nargin<2)
                modelPhysicalFeasibility = [];
                if(nargin<1)
                    modelPerformance = [];
                end
            end

            obj.modelPerformance = modelPerformance;
            obj.modelPhysicalFeasibility = modelPhysicalFeasibility;
        end
        
        function [performanceDeficit,physicalFeasibilityDeficit,evaluationOutput] = evaluate(obj,designSample)
        %EVALUATE Evaluation of designs' performance relative to requirements
        %   EVALUATE is the main function of any design evaluator, where design sample
        %   points are given and judgement must be made as to whether or not that design
        %   meets the system requirements.
        %   Here, the given design fast-forward models are used to make such a decision.
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
        %   process in EVALUATIONOUTPUT; in this case, that is the outputs
        %   of the predictions of the design fast-forward models.
        %
        %   Inputs:
        %       - OBJ : DesignEvaluator
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %   
        %   Outputs:
        %       - PERFORMANCEDEFECIT : (nSample,nPerformance) double
        %       - PHYSICALFEASIBILITYDEFICIT : (nSample,nPhysicalFeasibility) double
        %       - EVALUATIONOUTPUT : structure
        %           -- PerformanceOutput : class-defined
        %           -- PhysicalFeasibilityOutput : class-defined
        %   
        %   See also DesignFastForward.

            % performance prediction
            if(~isempty(obj.modelPerformance))
                [performanceDeficit,performanceOutput] = obj.modelPerformance.predict_designs(designSample);
                % predict designs outside training region as bad
                boundingBox = obj.modelPerformance.BoundingBoxPositive;
                insidePredictionRange = is_in_design_box(designSample,boundingBox);
                performanceDeficit(~insidePredictionRange,:) = abs(performanceDeficit(~insidePredictionRange,:));
            else
                % no model, assume all good
                performanceDeficit = -ones(size(designSample,1),1);
                performanceOutput = [];
            end

            % physical feasibility prediction
            if(~isempty(obj.modelPhysicalFeasibility))
                [physicalFeasibilityDeficit,physicalFeasibilityOutput] = obj.modelPhysicalFeasibility.predict_designs(designSample);
                % predict designs outside training region as physically infeasible
                boundingBox = obj.modelPhysicalFeasibility.BoundingBoxPositive;
                insidePredictionRange = is_in_design_box(designSample,boundingBox);
                physicalFeasibilityDeficit(~insidePredictionRange,:) = abs(physicalFeasibilityDeficit(~insidePredictionRange,:));
            else
                % no model, assume all physically feasible
                physicalFeasibilityDeficit = -ones(size(designSample,1),1);
                physicalFeasibilityOutput = [];
            end

            % wrap outputs
            evaluationOutput.PerformanceOutput = performanceOutput;
            evaluationOutput.PhysicalFeasibilityOutput = physicalFeasibilityOutput;
        end
    end
end