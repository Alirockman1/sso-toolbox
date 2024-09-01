classdef DesignEvaluatorBottomUpMapping < DesignEvaluatorBase
%DESIGNEVALUATORBOTTOMUPMAPPING Design Evaluator that uses bottom-up mappings
%   DESIGNEVALUATORBOTTOMUPMAPPING takes a bottom-up mapping object and uses it 
%   to determine if given designs are good/bad (and physically feasible/
%   infeasible).
%
%   This class is derived from DesignEvaluatorBase.
%
%   DESIGNEVALUATORBOTTOMUPMAPPING methods:
%       - evaluate : a function that receives the design sample points and  
%       returns their deficit in terms of performance and physical feasibility 
%       relative to the requirements.
%
%   See also DesignEvaluatorBase, DesignEvaluatorFastForward.
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
        %BOTTOMUPMAPPING Bottom-up Mapping Object
        %   BOTTOMUPMAPPING is the bottom-up mapping object that returns the performance 
        %   measures (and potentially physical feasibility measures) for a given design.
        %
        %   BOTTOMUPMAPPING : BottomUpMappingBase
        %   
        %   See also evaluate, BottomUpMappingBase.
        BottomUpMapping

        %PERFORMANCELOWERLIMIT Performance Measure Critical Lower Limit
        %   PERFORMANCELOWERLIMIT contains the system requirements for the performance 
        %   measures, in particular the requirements that specify the minimum values 
        %   those measures can assume.
        %
        %   PERFORMANCELOWERLIMIT : (1,nPerformance) double
        %   
        %   See also design_measure_to_deficit, evaluate, PerformanceUpperLimit.
        PerformanceLowerLimit

        %PERFORMANCEUPPERLIMIT Performance Measure Critical Upper Limit
        %   PERFORMANCEUPPERLIMIT contains the system requirements for the performance 
        %   measures, in particular the requirements that specify the maximum values 
        %   those measures can assume.
        %
        %   PERFORMANCEUPPERLIMIT : (1,nPerformance) double
        %   
        %   See also design_measure_to_deficit, evaluate, PerformanceLowerLimit.
        PerformanceUpperLimit

        %PERFORMANCENORMALIZATIONFACTOR Performance Measure Normalization Factor
        %   PERFORMANCENORMALIZATIONFACTOR contains normilization factors for the 
        %   performance measures that are used when computing the individual performance
        %   deficits for a given design.
        %   If any entry assumes the value 'nan', the normalization factor is computed 
        %   on-the-fly based on the lower/upper limits and measure values.
        %   If empty, all entries will assume the value 'nan'. 
        %
        %   PERFORMANCENORMALIZATIONFACTOR : (1,nPerformance) double
        %   
        %   See also design_measure_to_deficit, evaluate.
        PerformanceDefaultNormalizationFactor

        %PHYSICALFEASIBILITYLOWERLIMIT Physical Feasibility Critical Lower Limit
        %   PHYSICALFEASIBILITYLOWERLIMIT contains the system requirements for the 
        %   physical feasibility measures, in particular the requirements that specify 
        %   the minimum values those measures can assume.
        %
        %   PHYSICALFEASIBILITYLOWERLIMIT : (1,nPhysicalFeasibility) double
        %   
        %   See also design_measure_to_deficit, evaluate.
        PhysicalFeasibilityLowerLimit

        %PHYSICALFEASIBILITYUPPERLIMIT Physical Feasibility Critical Upper Limit
        %   PHYSICALFEASIBILITYUPPERLIMIT contains the system requirements for the 
        %   physical feasibility measures, in particular the requirements that specify 
        %   the maximum values those measures can assume.
        %
        %   PHYSICALFEASIBILITYUPPERLIMIT : (1,nPhysicalFeasibility) double
        %   
        %   See also design_measure_to_deficit, evaluate.
        PhysicalFeasibilityUpperLimit

        %PHYSICALFEASIBILITYNORMALIZATIONFACTOR Physical Feasibility Factor
        %   PHYSICALFEASIBILITYNORMALIZATIONFACTOR contains normalization factors for 
        %   the physical feasibility measures that are used when computing the 
        %   individual physical feasibility deficits for a given design.
        %   If any entry assumes the value 'nan', the normalization factor is computed 
        %   on-the-fly based on the lower/upper limits and measure values.
        %   If empty, all entries will assume the value 'nan'. 
        %
        %   PHYSICALFEASIBILITYNORMALIZATIONFACTOR : (1,nPhysicalFeasibility) double
        %   
        %   See also design_measure_to_deficit, evaluate.
        PhysicalFeasibilityDefaultNormalizationFactor
    end
    
    methods
        function obj = DesignEvaluatorBottomUpMapping(bottomUpMapping,varargin)
        %DESIGNEVALUATORBOTTOMUPMAPPING Constructor
        %   DESIGNEVALUATORBOTTOMUPMAPPING uses a given bottom-up mapping object and 
        %   other parameters to create a design evaluator. 
        %
        %   OBJ = DESIGNEVALUATORBOTTOMUPMAPPING(BOTTOMUPMAPPING) receives the defined
        %   bottom-up mapping object in BOTTOMUPMAPPING and returns a 
        %   DESIGNEVALUATORBOTTOMUPMAPPING object in OBJ. In said case, all other 
        %   input arguments are assumed to have their default values.
        %
        %   OBJ = DESIGNEVALUATORBOTTOMUPMAPPING(BOTTOMUPMAPPING,PERFORMANCELOWERLIMIT) 
        %   allows setting the minimum values acceptable for the performance measures, 
        %   passed in PERFORMANCELOWERLIMIT. Default value is '-inf'.
        %
        %   OBJ = DESIGNEVALUATORBOTTOMUPMAPPING(BOTTOMUPMAPPING,PERFORMANCELOWERLIMIT,
        %   PERFORMANCEUPPERLIMIT) allows setting the maximum values acceptable for the 
        %   performance measures, passed in PERFORMANCEUPPERLIMIT. Default value is '0'.
        %
        %   OBJ = DESIGNEVALUATORBOTTOMUPMAPPING(...,NAME1,VALUE1,...) also allows for 
        %   setting custom options for the evaluation procedure, passed as 
        %   name-value pairs. These are:
        %       - 'PerformanceNormalizationFactor' : fixed normalization factors
        %       to use when computing the performance deficit. Default value is
        %       empty, meaning said factor is computed on-the-fly.
        %       - 'PhysicalFeasibilityLowerLimit' : minimum values acceptable 
        %       for the physical feasibility measures. Default value is '-inf'.
        %       - 'PhysicalFeasibilityUpperLimit' : maximum values acceptable 
        %       for the physical feasibility measures. Default value is '0'.
        %       - 'PhysicalFeasibilityNormalizationFactor' : fixed normalization 
        %       factors to use when computing the performance deficit. Default 
        %       value is empty, meaning said factor is computed on-the-fly.
        %
        %   Inputs:
        %       - BOTTOMUPMAPPING : BottomUpMappingBase
        %       - PARAMETER : function-dependent
        %       - PERFORMANCELOWERLIMIT : (1,nPerformance) double
        %       - PERFORMANCEUPPERLIMIT : (1,nPerformance) double
        %       - 'PerformanceNormalizationFactor' : (1,nPerformance) double
        %       - 'PhysicalFeasibilityLowerLimit' : (1,nPhysicalFeasibility) double
        %       - 'PhysicalFeasibilityUpperLimit' : (1,nPhysicalFeasibility) double
        %       - 'PhysicalFeasibilityNormalizationFactor' : (1,nPhysicalFeasibility) 
        %       double
        %   
        %   Outputs:
        %       - OBJ : DesignEvaluatorBottomUpMapping
        %   
        %   See also BottomUpMappingBase, evaluate.

            parser = inputParser;
            parser.addRequired('BottomUpMapping', @(x)isa(x,'BottomUpMappingBase'));
            parser.addOptional('PerformanceLowerLimit',-inf,@(x)size(x,1)==1);
            parser.addOptional('PerformanceUpperLimit',0,@(x)size(x,1)==1);
            parser.addParameter('PerformanceNormalizationFactor',[],@(x)isnumeric(x)&&(size(x,1)==1 || isempty(x)));
            parser.addParameter('PhysicalFeasibilityLowerLimit',-inf,@(x)isnumeric(x)&&(size(x,1)==1 || isempty(x)));
            parser.addParameter('PhysicalFeasibilityUpperLimit',0,@(x)isnumeric(x)&&(size(x,1)==1 || isempty(x)));
            parser.addParameter('PhysicalFeasibilityNormalizationFactor',[],@(x)isnumeric(x)&&(size(x,1)==1 || isempty(x)));
            parser.parse(bottomUpMapping,varargin{:});

            obj.BottomUpMapping = parser.Results.BottomUpMapping;
            obj.PerformanceLowerLimit = parser.Results.PerformanceLowerLimit;
            obj.PerformanceUpperLimit = parser.Results.PerformanceUpperLimit;
            obj.PerformanceDefaultNormalizationFactor = parser.Results.PerformanceNormalizationFactor;
            obj.PhysicalFeasibilityLowerLimit = parser.Results.PhysicalFeasibilityLowerLimit;
            obj.PhysicalFeasibilityUpperLimit = parser.Results.PhysicalFeasibilityUpperLimit;
            obj.PhysicalFeasibilityDefaultNormalizationFactor = parser.Results.PhysicalFeasibilityNormalizationFactor;
        end
        
        function [performanceDeficit,physicalFeasibilityDeficit,evaluationOutput] = evaluate(obj,designSample)
        %EVALUATE Evaluation of designs' performance relative to requirements
        %   EVALUATE is the main function of any design evaluator, where design sample
        %   points are given and judgement must be made as to whether or not that design
        %   meets the system requirements.
        %   Here, the given bottom-up mapping and limits are used to make said decision.
        %   The bottom-up mapping is used to acquire the performance measures (and 
        %   physical feasibility measures, if applicable), and those are used to compute
        %   their deficits relative to the given requirements.
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
        %   process in EVALUATIONOUTPUT; in this case, that is the performance and
        %   physical feasibility measures which are obtained from the response of the 
        %   bottom-up mapping.
        %
        %   Input:
        %       - OBJ : DesignEvaluatorBase
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %   
        %   Output:
        %       - PERFORMANCEDEFICIT : (nSample,nPerformance) double
        %       - PHYSICALFEASIBILITYDEFICIT : (nSample,nPhysicalFeasibility) double
        %       - EVALUATIONOUTPUT : structure
        %           -- PerformanceMeasure : (nSample,nPerformance) double
        %           -- PhysicalFeasibilityMeasure : (nSample,nPhysicalFeasibility) 
        %           double
        %           -- BottomUpMappingOutput : class-defined
        %           -- PerformanceNormalizationFactor : (1,nPerformance) double
        %           -- PhysicalFeasibilityNormalizationFactor : (1,nPhysicalFeasibility) 
        %           double
        %
        %   See also design_measure_to_deficit.

            % call function
            [performanceMeasure,physicalFeasibilityMeasure,bottomUpMappingOutput] = obj.BottomUpMapping.response(designSample);
            evaluationOutput.BottomUpMappingOutput = bottomUpMappingOutput;

            % wrap performance
            [performanceDeficit,performanceNormalizationFactor] = ...
                design_measure_to_deficit(...
                    performanceMeasure,...
                    obj.PerformanceLowerLimit,...
                    obj.PerformanceUpperLimit,...
                    obj.PerformanceDefaultNormalizationFactor);
            evaluationOutput.PerformanceMeasure = performanceMeasure;
            evaluationOutput.PerformanceNormalizationFactor = performanceNormalizationFactor;

            % wrap physical feasibility
            [physicalFeasibilityDeficit,physicalFeasibilityNormalizationFactor] = ...
                design_measure_to_deficit(...
                    physicalFeasibilityMeasure,...
                    obj.PhysicalFeasibilityLowerLimit,...
                    obj.PhysicalFeasibilityUpperLimit,...
                    obj.PhysicalFeasibilityDefaultNormalizationFactor);
            evaluationOutput.PhysicalFeasibilityMeasure = physicalFeasibilityMeasure;
            evaluationOutput.PhysicalFeasibilityNormalizationFactor = physicalFeasibilityNormalizationFactor;
        end
    end
end