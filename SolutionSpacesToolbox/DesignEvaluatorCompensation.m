classdef DesignEvaluatorCompensation < DesignEvaluatorBase
%DESIGNEVALUATORCOMPENSATION Evaluator for Solution-compensation Spaces
%   DESIGNEVALUATORCOMPENSATION separates the design variables in two spaces:
%   an A-space (for early decision variables) and a B-space (for late decision
%   variables). The evaluator receives design samples in A-space, and uses
%   optimization to find the design variables in B-space that give the best
%   result for that A-space design. 
%   This can be used for the computation of solution-compensation spaces.
%
%   This class is derived from DesignEvaluatorBase.
%
%   DESIGNEVALUATORCOMPENSATION methods:
%       - evaluate : a function that receives the design sample points and  
%       returns their deficit in terms of performance and physical feasibility 
%       relative to the requirements.
%
%   See also DesignEvaluatorBase, DesignEvaluatorBottomUpMapping.
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
        %BASEDESIGNEVALUATOR Design evaluator whose performance score will be optimized
        %   BASEDESIGNEVALUATOR is the basic design evaluator of the problem, which will
        %   be used for individual evaluations and whose score will be optimized during
        %   the A-space evaluation process.
        %
        %   BASEDESIGNEVALUATOR : DesignEvaluatorBase
        %   
        %   See also evaluate.
        BaseDesignEvaluator

        %COMPENSATIONASPACEINDEX Logical index of A-space design variables
        %   COMPENSATIONASPACEINDEX is the logical index which indicates which design
        %   variables belong in the A-space as opposed to the B-space. 
        %   This is used during the construction of the object to separate the design 
        %   space, components, and initial design if necessary, as well as during 
        %   evaluation to reconstruct the design variable input as expected by the base  
        %   design evaluator.
        %
        %   COMPENSATIONASPACEINDEX : (1,nDesignVariable)
        %   
        %   See also BaseDesignEvaluator.
        CompensationAspaceIndex

        %BSPACELOWERBOUND Design space lower boundary for the variables in B-space
        %   BSPACELOWERBOUND is the lower boundary of the design space for the design
        %   variables in B-space, which serves as lower boundary for the optimization
        %   procedure for the evaluation of every sample.
        %
        %   BSPACELOWERBOUND : (1,nBspaceDesignVariable)
        %
        %   See also BspaceUpperBound, CompensationAspaceIndex.
        BspaceLowerBound

        %BSPACEUPPERBOUND Design space upper boundary for the variables in B-space
        %   BSPACEUPPERBOUND is the upper boundary of the design space for the design
        %   variables in B-space, which serves as upper boundary for the optimization
        %   procedure for the evaluation of every sample.
        %
        %   BSPACEUPPERBOUND : (1,nBspaceDesignVariable)
        %   
        %   See also BspaceLowerBound, CompensationAspaceIndex.
        BspaceUpperBound

        %BSPACEINITIALDESIGN Starting design point for variables in B-space
        %   BSPACEINITIALDESIGN is the initial starting point for the variables in 
        %   B-space when the evaluation process begins for a new sample. For each 
        %   subsequent sample point, however, the optimized values of the previous
        %   sample point is used as new starting point for the optimization.
        %
        %   BSPACEINITIALDESIGN : (1,nBspaceDesignVariable)
        %   
        %   See also CompensationAspaceIndex.
        BspaceInitialDesign

        %OPTIMIZATIONOPTIONS Options for the optimization of the score
        %   OPTIMIZATIONOPTIONS are the options for the score optimization procedure.
        %   See the documentation of 'design_optimize_performance_score' for more 
        %   details.
        %
        %   OPTIMIZATIONOPTIONS : (1,nOptionOptimization) cell
        %   
        %   See also design_optimize_performance_score.
        OptimizationOptions

        %ASPACEORDERPASSES Number of passes in the sorting process of A-space samples
        %   ASPACEORDERPASSES is the number of passes done when sorting the A-space 
        %   samples by distance from each other (to minimize said distance). 
        %   This sorting is done before the evaluation procedure, to then use the 
        %   previous optimized B-space result as the new initial design for each 
        %   subsequent sample point and hopefully reduce the number of optimization 
        %   iterations necessary.
        %   
        %   ASPACEORDERPASSES : integer
        %   
        %   See also design_sort_min_distance.
        AspaceOrderPasses

        %ASPACEORDERKNNSEARCHOPTIONS Options for the sorting process of A-space samples
        %   ASPACEORDERKNNSEARCHOPTIONS are the options used by 'knnsearch' when 
        %   sorting the A-space samples.
        %   This sorting is done before the evaluation procedure, to then use the 
        %   previous optimized B-space result as the new initial design for each 
        %   subsequent sample point and hopefully reduce the number of optimization 
        %   iterations necessary.
        %
        %   ASPACEORDERKNNSEARCHOPTIONS : (1,nOptionDistance) cell
        %   
        %   See also design_sort_min_distance, knnsearch.
        AspaceOrderKnnsearchOptions
    end

    properties
        %ACCEPTFULLSPACESAMPLE Use B-space part of sample points as initial design
        %   ACCEPTFULLSPACESAMPLE determines whether the input of the evaluate method
        %   can be given in the full space (including both A-space and B-space design
        %   variable values). If enabled, the B-space part of the sample points is used
        %   as initial guess for the optimization; in this case, the A-space sample 
        %   points are not ordered, nor is the previous result used.
        %   ACCEPTFULLSPACESAMPLE is introduced only to stop the use of such 
        %   functionality in cases where it wouldn't be appropriate, such as the 
        %   computation of solution spaces.
        %
        %   ACCEPTFULLSPACESAMPLE : logical
        %   
        %   See also evaluate.
        AcceptFullSpaceSample
    end
    
    methods
        function [obj,aspaceLowerBound,aspaceUpperBound,aspaceInitialDesign,aspaceComponentIndex] = DesignEvaluatorCompensation(...
            baseDesignEvaluator,compensationAspaceIndex,designSpaceLowerBound,designSpaceUpperBound,varargin)
        %DESIGNEVALUATORCOMPENSATION Constructor
        %   DESIGNEVALUATORCOMPENSATION creates an object of the 
        %   DesignEvaluatorCompensation class, while also adapting other important 
        %   variables (such as design space, component index) to A-space as well.
        %
        %   OBJ = DESIGNEVALUATORCOMPENSATION(BASEDESIGNEVALUATOR,
        %   COMPENSATIONASPACEINDEX,DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND) uses 
        %   the basic design evaluator BASEDESIGNEVALUATOR, the compensation index 
        %   COMPENSATIONASPACEINDEX and the design space boundaries 
        %   DESIGNSPACELOWERBOUND and DESIGNSPACEUPPERBOUND to create a 
        %   DesignEvaluatorCompensation object OBJ. The basic design evaluator is
        %   used in the evaluation procedure and its score for each sample point
        %   is minimized. The compensation index is used to determine which design 
        %   variables belong in the A-space (true) and which ones belong in the B-space
        %   (false). The boundaries of the design space can be given either in full or
        %   exclusively for the B-space.
        %
        %   OBJ = DESIGNEVALUATORCOMPENSATION(BASEDESIGNEVALUATOR,
        %   COMPENSATIONASPACEINDEX,DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND,
        %   INITIALDESIGN) also specifies the initial design INITIALDESIGN, where its
        %   B-space part is then used as starting point for the first optimization 
        %   during the evaluation procedure. It can be given either in full or
        %   exclusively for the B-space.
        %
        %   OBJ = DESIGNEVALUATORCOMPENSATION(BASEDESIGNEVALUATOR,
        %   COMPENSATIONASPACEINDEX,DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND,
        %   INITIALDESIGN,COMPONENTINDEX) also specifies the component index 
        %   COMPONENTINDEX; this is only used for conversion of that to the A-space.
        %
        %   OBJ = DESIGNEVALUATORCOMPENSATION(...NAME,VALUE,...) also allows for setting
        %   custom options for the evaluation procedure, passed as name-value pair arguments
        %   arguments. These are:
        %       - 'OptimizationOptions' : optimization options used by 
        %       'design_optimize_performance_score'. Default: 
        %       {'OptimizationMethodOptions',{'Display','none'},
        %       'CompensationAspaceIndex',compensationAspaceIndex}.
        %       - 'AspaceOrderPasses' : number of passes when ordering the A-space 
        %       sample points by distance before evaluation. Default: 1.
        %       - 'AspaceOrderKnnsearchOptions' : options used when ordering the A-space
        %       sample points by distance before evaluation. Default is empty.
        %        - 'AcceptFullSpaceSample' : flag to know if the input sample points for
        %       the evaluate method can be given in the complete design space instead
        %       of just the A-space part. For situations like computing solution-
        %       compensation spaces, such inputs would be undesirable and would suggest
        %       an error in setup. Default: false.
        %
        %   [OBJ,ASPACELOWERBOUND] = DESIGNEVALUATORCOMPENSATION(...) additionally 
        %   returns the part of the lower bound of the design space that is exclusive to
        %   the A-space.
        %
        %   [OBJ,ASPACELOWERBOUND,ASPACEUPPERBOUND] = DESIGNEVALUATORCOMPENSATION(...)
        %   additionally returns the part of the upper bound of the design space that is
        %   exclusive to the A-space.
        %
        %   [OBJ,ASPACELOWERBOUND,ASPACEUPPERBOUND,ASPACEINITIALDESIGN] = 
        %   DESIGNEVALUATORCOMPENSATION(...) additionally returns the part of the 
        %   initial design that is exclusive to the A-space.
        %
        %   [OBJ,ASPACELOWERBOUND,ASPACEUPPERBOUND,ASPACEINITIALDESIGN,
        %   ASPACECOMPONENTINDEX] = DESIGNEVALUATORCOMPENSATION(...) additionally  
        %   returns the component index with indices adapted to the A-space.
        %
        %   Input:
        %       - BASEDESIGNEVALUATOR : DesignEvaluatorBase
        %       - COMPENSATIONASPACEINDEX : (1,nDesignVariable) logical
        %       - DESIGNSPACELOWERBOUND : (1,nDesignVariable) double OR 
        %       (1,nBspaceDesignVariable) double
        %       - DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double OR 
        %       (1,nBspaceDesignVariable) double
        %       - INITIALDESIGN : (1,nDesignVariable) double OR 
        %       (1,nBspaceDesignVariable) double
        %       - COMPONENTINDEX : (1,nComponent) cell
        %       - 'OptimizationOptions' : (1,nOptionsOptimization) cell
        %       - 'AspaceOrderPasses' : integer
        %       - 'AspaceOrderKnnsearchOptions' : (1,nOptionsDistance) cell
        %   
        %   Output:
        %       - OBJ : DesignEvaluatorCompensation
        %       - ASPACELOWERBOUND : (1,nAspaceDesignVariable) double
        %       - ASPACEUPPERBOUND : (1,nAspaceDesignVariable) double
        %       - ASPACEINITIALDESIGN : (1,nAspaceDesignVariable) double
        %       - ASPACECOMPONENTINDEX : (1,nComponent) cell
        %   
        %   See also evaluate, design_optimize_performance_score, 
        %   design_sort_min_distance.

            parser = inputParser;
            parser.addOptional('InitialDesign',[]);
            parser.addOptional('ComponentIndex',{});
            parser.addParameter('OptimizationOptions',{});
            parser.addParameter('AspaceOrderPasses',1);
            parser.addParameter('AspaceOrderKnnsearchOptions',{});
            parser.addParameter('AcceptFullSpaceSample',false,@islogical);
            parser.parse(varargin{:});
            options = parser.Results;

            nDesignVariable = size(compensationAspaceIndex,2);
            if(size(designSpaceLowerBound,2)==nDesignVariable)
                bspaceLowerBound = designSpaceLowerBound(~compensationAspaceIndex);
                aspaceLowerBound = designSpaceLowerBound(compensationAspaceIndex);
            else
                bspaceLowerBound = designSpaceLowerBound;
                aspaceLowerBound = [];
            end
            if(size(designSpaceUpperBound,2)==nDesignVariable)
                bspaceUpperBound = designSpaceUpperBound(~compensationAspaceIndex);
                aspaceUpperBound = designSpaceUpperBound(compensationAspaceIndex);
            else
                bspaceUpperBound = designSpaceUpperBound;
                aspaceUpperBound = [];
            end
            if(size(options.InitialDesign,2)==nDesignVariable)
                bspaceInitialDesign = options.InitialDesign(~compensationAspaceIndex);
                aspaceInitialDesign = options.InitialDesign(compensationAspaceIndex);
            else
                bspaceInitialDesign = options.InitialDesign;
                aspaceInitialDesign = [];
            end

            aspaceComponentIndex = cell(size(options.ComponentIndex));
            for i=1:length(aspaceComponentIndex)
                aspaceComponentIndex{i} = convert_index_base(compensationAspaceIndex,options.ComponentIndex{i}','forward')';
            end

            obj.BaseDesignEvaluator = baseDesignEvaluator;
            obj.CompensationAspaceIndex = compensationAspaceIndex;
            obj.BspaceLowerBound = bspaceLowerBound;
            obj.BspaceUpperBound = bspaceUpperBound;
            obj.BspaceInitialDesign = bspaceInitialDesign;
            [~,obj.OptimizationOptions] = merge_name_value_pair_argument(...
                    {'OptimizationMethodOptions',{'Display','none'}},...
                    parser.Results.OptimizationOptions,...
                    {'CompensationAspaceIndex',compensationAspaceIndex});
            obj.AspaceOrderPasses = options.AspaceOrderPasses;
            obj.AspaceOrderKnnsearchOptions = options.AspaceOrderKnnsearchOptions;
            obj.AcceptFullSpaceSample = options.AcceptFullSpaceSample;
        end
        
        function [performanceDeficit,physicalFeasibilityDeficit,evaluationOutput] = evaluate(obj,designSample)
        %EVALUATE Evaluation of designs' performance relative to requirements
        %   EVALUATE is the main function of any design evaluator, where design sample
        %   points are given and judgement must be made as to whether or not that design
        %   meets the system requirements.
        %   Here, a compensation factor is also taken into consideration, where given
        %   a sample of the early decision variables ("A-space"), an optimization 
        %   sub-problem is solved for each sample point where late-decision variables 
        %   ("B-space") are varied to find an optimal design for those early-decision 
        %   variables. 
        %   Furthermore, if the design sample points are given in A-space, they are
        %   ordered to minimize the distance from one sample to the other, and the
        %   optimized B-space design from the previous A-space sample is used as the 
        %   new initial design in the optimization process. If the design sample
        %   points given are in the full space, however, the B-space design variables
        %   are used as initial designs for each sample point instead, and no 
        %   ordering is used.
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
        %   process in EVALUATIONOUTPUT. In this case, this includes information on the
        %   ordering of A-space sample points, the optimal designs in B-space for each
        %   A-space sample, and the outputs of both the optimization method and the
        %   basic design evaluator.
        %
        %   Input:
        %       - OBJ : DesignEvaluatorBase
        %       - DESIGNSAMPLE : (nSample,nDesignVariable) double
        %   
        %   Output:
        %       - PERFORMANCEDEFICIT : (nSample,nPerformance) double
        %       - PHYSICALFEASIBILITYDEFICIT : (nSample,nPhysicalFeasibility) double
        %       - EVALUATIONOUTPUT : structure
        %           -- AspaceMinPathOrder : (nSample,1) integer
        %           -- AspaceMinPathDistance : (nSample,1) double
        %           -- BspaceOptimal : (nSample,nBspaceDesignVariable) double
        %           -- PerformanceScoreOptimal : (nSample,1) double
        %           -- OptimizationMethodOutput : function-dependent
        %           -- BaseEvaluatorOutputOptimization : class-dependent
        %           -- BaseEvaluatorOutputOptimal : class-dependent
        %   
        %   See also design_sort_min_distance, design_optimize_performance_score.
            
            nSample = size(designSample,1);
            nDesignVariable = size(obj.CompensationAspaceIndex,2);
            nDesignVariableAspace = sum(obj.CompensationAspaceIndex);
            nDesignVariableBspace = size(obj.BspaceLowerBound,2);

            isFullSpaceSample = (size(designSample,2)==nDesignVariable);
            if(isFullSpaceSample)
                if(obj.AcceptFullSpaceSample)
                    % use the given B-space value as initial designs for each optimization
                    bspaceInitialSample = designSample(:,~obj.CompensationAspaceIndex);
                    designSample = designSample(:,obj.CompensationAspaceIndex);

                    % do not order
                    minDistanceOrder = 1:nSample;
                    minPathDistance = [];
                else
                    error('DesignEvaluatorCompensation:FullSpaceSampleNotAccepted',...
                        ['Sample points for compensation given in full space form.\n',...
                        'If you are computing solution-compensation spaces, this suggests an error in setup.\n',...
                        'If not, and this is expected, remember to change property ''AcceptFullSpaceSample'' to true.']);
                end
            elseif(size(designSample,2)==nDesignVariableAspace)
                % start B-space sample in the middle if value was not given
                bspaceInitial = obj.BspaceInitialDesign;
                if(isempty(bspaceInitial))
                    bspaceInitial = mean([obj.BspaceLowerBound;obj.BspaceUpperBound],1);
                end
    
                % order A-space designs to assist with convergence
                [minDistanceOrder,minPathDistance] = design_sort_min_distance(designSample,...
                    obj.AspaceOrderPasses,obj.AspaceOrderKnnsearchOptions{:});
            else
                error('DesignEvaluatorCompensation:SampleNotCompatibleSize',...
                    'Sample given to DesignEvaluatorCompensation is not of a compatible size.');
            end

            % initialize outputs and other information to be logged
            scoreOptimal = nan(nSample,1);
            optimizationOutput = cell(nSample,1);
            evaluatorOutput = cell(nSample,1);
            bspaceOptimal = nan(nSample,nDesignVariableBspace);

            % find the optimal B-space values for each design
            for i=1:nSample
                iCurrent = minDistanceOrder(i);

                % if full sample was given, use that for the initial design of B-space
                if(isFullSpaceSample)
                    bspaceInitial = bspaceInitialSample(iCurrent,:);
                end

                initialDesign = nan(1,size(obj.CompensationAspaceIndex,2));
                initialDesign(obj.CompensationAspaceIndex) = designSample(iCurrent,:);
                initialDesign(~obj.CompensationAspaceIndex) = bspaceInitial;

                [optimalDesign,scoreOptimalCurrent,optimizationOutputCurrent,evaluatorOutputCurrent] = ...
                    design_optimize_performance_score(...
                        obj.BaseDesignEvaluator,...
                        initialDesign,...
                        obj.BspaceLowerBound,...
                        obj.BspaceUpperBound,...
                        obj.OptimizationOptions{:});

                % log extra information
                bspaceOptimal(iCurrent,:) = optimalDesign(~obj.CompensationAspaceIndex);
                scoreOptimal(iCurrent) = scoreOptimalCurrent;
                optimizationOutput{iCurrent} = optimizationOutputCurrent;
                evaluatorOutput{iCurrent} = evaluatorOutputCurrent;

                % start next iteration already where you are
                % --> with A-space ordered, this should improve convergence time
                if(~isFullSpaceSample)
                    bspaceInitial = bspaceOptimal(iCurrent,:);
                end
            end

            % evaluate optimal designs
            optimalDesignSample = nan(nSample,nDesignVariable);
            optimalDesignSample(:,obj.CompensationAspaceIndex) = designSample;
            optimalDesignSample(:,~obj.CompensationAspaceIndex) = bspaceOptimal;

            [performanceDeficit,physicalFeasibilityDeficit,evaluatorOutputFinal] = ...
                obj.BaseDesignEvaluator.evaluate(optimalDesignSample);

            % save extra information
            evaluationOutput.AspaceMinPathOrder = minDistanceOrder;
            evaluationOutput.AspaceMinPathDistance = minPathDistance;
            evaluationOutput.BspaceOptimal = bspaceOptimal;
            evaluationOutput.PerformanceScoreOptimal = scoreOptimal;
            evaluationOutput.OptimizationMethodOutput = optimizationOutput;
            evaluationOutput.BaseEvaluatorOutputOptimization = evaluatorOutput;
            evaluationOutput.BaseEvaluatorOutputOptimal = evaluatorOutputFinal;
        end
    end
end