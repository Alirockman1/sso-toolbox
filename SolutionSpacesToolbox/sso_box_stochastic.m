function [candidateBox,optimizationData] = sso_box_stochastic(designEvaluator,initialBox,designSpaceLowerBound,designSpaceUpperBound,varargin)
%SSO_BOX_STOCHASTIC Box-shaped solution spaces optimization (Stochastic method)
%   SSO_BOX_STOCHASTIC uses a modified version of the stochastic method to 
%   compute optimal soluton (or requirement) spaces. 
%
%   CANDIDATEBOX = SSO_BOX_STOCHASTIC(DESIGNEVALUATOR,INITIALBOX,
%   DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND) starts the algorithm in 
%   INITIALBOX and finds the optimum box-shaped solution (or requirement) spaces
%   within the design space defined by DESIGNSPACELOWERBOUND and 
%   DESIGNSPACEUPPERBOUND, evaluating design sample points with DESIGNEVALUATOR
%   and returning the optimal box CANDIDATEBOX.
%
%   CANDIDATEBOX = SSO_BOX_STOCHASTIC(DESIGNEVALUATOR,INITIALBOX,
%   DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND,OPTIONS) also allows one to
%   change the options being used through the output of 'sso_stochastic_options'
%   OPTIONS.
%
%   [CANDIDATEBOX,OPTIMIZATIONDATA] = SSO_BOX_STOCHASTIC(...) additionally returns
%   the fixed problem and iteration data in OPTIMIZATIONDATA.
%
%
%   Input:
%       - DESIGNEVALUATOR : DesignEvaluatorBase
%       - INITIALBOX : (1,nDesignVariable) double OR (2,nDesignVariable) double
%       - DESIGNSPACELOWERBOUND : (1,nDesignVariable) double
%       - DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double
%       - OPTIONS : structure
%
%   Output:
%       - CANDIDATEBOX  : (2,nDesignVariable) double
%       - OPTIMIZATIONDATA : struct
%           - DesignEvaluator : DesignEvaluatorBase 
%           - InitialBox : (2,nDesignVariable) double
%           - DesignSpaceLowerBound : (1,nDesignVariable) double
%           - DesignSpaceUpperBound : (1,nDesignVariable) double
%           - Options : struct
%           - InitialRNGState : struct
%           - IterationData : (nIteration,1) struct
%               - EvaluatedDesignSamples : (nSample,nDesignVariable) double
%               - EvaluationOutput : class-dependent
%               - Phase : integer
%               - GrowthRate : double
%               - DesignScore : (nSample,1) double
%               - IsGoodPerformance : (nSample,1) logical
%               - IsPhysicallyFeasible : (nSample,1) logical
%               - IsAcceptable : (nSample,1) logical
%               - IsUseful : (nSample,1) logical
%               - CandidateBoxBeforeTrim : (2,nDesignVariable) double
%               - CandidateBoxAfterTrim : (2,nDesignVariable) double
%
%   See also sso_stochastic_options.
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
    
    % Options
    defaultOptions = sso_stochastic_options('box');
    inputOptions = parser_variable_input_to_structure(varargin{:});
    options = merge_name_value_pair_argument(defaultOptions,inputOptions);
    
    %% extract options as necessary
    % requirement spaces
    requirementSpacesType = options.RequirementSpacesType;
    applyLeannessEachTrim = strcmpi(options.ApplyLeanness,'always');
    applyLeannessFinalTrim = any(strcmpi(options.ApplyLeanness,{'always','end-only'}));

    % trimming
    [~,trimmingOperationOptions] = merge_name_value_pair_argument(...
        {'MeasureFunction',options.MeasureFunction,...
            'MeasureOptions',options.MeasureOptions},...
        options.TrimmingOperationOptions);

    % logging verbosity
    console = ConsoleLogging(options.LoggingLevel);
    
    % Initial Candidate Box
    if(size(initialBox,1)==1)
        candidateBox = [initialBox;initialBox]; % single point
    elseif(size(initialBox,1)==2)
        candidateBox = initialBox; % candidate box
    else
        console.error('SSOBoxOptStochastic:InitialGuessWrong','Error. Initial guess incompatible in ''sso_box_stochastic''.');
    end
    
    % Initial Measure
    measureTrimmed = options.MeasureFunction(candidateBox, [], options.MeasureOptions{:});
    if(isinf(measureTrimmed) || isnan(measureTrimmed))
        measureTrimmed = 0;
    end
    nDimension = size(designSpaceLowerBound,2);
    
    %% Log Initialization
    isOutputOptimizationData = (nargout>=2);
    if(isOutputOptimizationData)
        optimizationData = struct(...
            'DesignEvaluator',designEvaluator,...
            'InitialBox',candidateBox,...
            'DesignSpaceLowerBound',designSpaceLowerBound,...
            'DesignSpaceUpperBound',designSpaceUpperBound,...
            'Options',options,...
            'InitialRNGState',rng);
        iLog = 1;
    end

    %% Phase I - Exploration
    iExploration = 1;
    growthRate = options.GrowthRate;
    hasConvergedExploration = false;
    while((~hasConvergedExploration) && (iExploration<=options.MaxIterExploration))
        %% Modification Step B - Growth: Extend Candidate Box
        console.info([repelem('=',120),'\n']);
        console.info('Initiating Phase I - Exploration: Iteration %d\n',iExploration);
        
        % Change growth rate depending on previous result
        if(iExploration>1 && options.UseAdaptiveGrowthRate)
            console.info('Adapting growth rate... ');
            tic

            purity = max(min(purity,options.MaximumGrowthPurity),options.MinimumGrowthPurity);

            increaseMeasure = measureGrown - measurePrevious;
            increaseMeasureAcceptable = max(measureGrown*purity - measurePrevious,0);
            fractionAcceptableIncreaseMeasure = increaseMeasureAcceptable/increaseMeasure;

            % Change step size to a bigger or smaller value depending on whether
            % the achieved purity is smaller or larger than the desired one
            growthAdaptationFactor = options.GrowthAdaptationFactorFunction(...
                purity,options.TargetAcceptedRatioExploration,nDimension,fractionAcceptableIncreaseMeasure,...
                options.GrowthAdaptationFactorOptions{:});
            growthAdaptationFactor = max(min(growthAdaptationFactor,options.MaximumGrowthAdaptationFactor),options.MinimumGrowthAdaptationFactor);

            growthRate = growthAdaptationFactor * growthRate;
            growthRate = max(min(growthRate,options.MaximumGrowthRate),options.MinimumGrowthRate);

            console.info('Elapsed time is %g seconds.\n',toc);
        end
        measurePrevious = measureTrimmed;
        
        % Where design variables aren't fixed, expand candidate solution box 
        %   in both sides of each interval isotroply
        console.info('Growing candidate box... ');
        tic

        candidateBoxGrown = design_box_grow_fixed(candidateBox,designSpaceLowerBound,designSpaceUpperBound,growthRate);
        
        console.info('Elapsed time is %g seconds.\n',toc);
        console.debug('- Current Growth Rate: %g\n',growthRate);
        
        %% Sample inside the current candidate box
        % get current number of samples
        nSample = get_current_array_entry(options.NumberSamplesPerIterationExploration,iExploration);
        
        % Generate samples that are to be evaluated
        designSample = sso_box_sub_generate_new_sample_points(...
            candidateBoxGrown,...
            nSample,...
            options.SamplingMethodFunction,...
            options.SamplingMethodOptions,...
            console);
        
        % Evaluate the samples
        [isGoodPerformance,isPhysicallyFeasible,score,outputEvaluation] = sso_box_sub_evaluate_sample_points(...
            designEvaluator,...
            designSample,...
            console);
        
        % Label samples according to desired requirement spaces problem type
        [isAcceptable,isUseful] = sso_box_sub_label_samples_requirement_spaces(...
            requirementSpacesType,...
            isGoodPerformance,...
            isPhysicallyFeasible,...
            console);
        
        % Count number of labels
        [nAcceptable,nUseful,nAcceptableUseful] = sso_box_sub_count_label_acceptable_useful(...
            isAcceptable,...
            isUseful,...
            console);
        
        % Compute candidate box measure
        measureGrown = sso_box_sub_compute_candidate_box_measure(...
            candidateBoxGrown,...
            nSample,...
            nUseful,...
            options.MeasureFunction,...
            options.MeasureOptions,...
            console);
        
        % Box may have grown too much + "unlucky sampling" w/ no good
        % points, go back in this case
        purity = nAcceptable/nSample;
        if(nAcceptable==0 || nUseful==0 || nAcceptableUseful==0 || purity<options.MinimumPurityReset)
            console.warn('SSOOptBox:BadSampling',['Not enough good/useful points found, ',...
                'rolling back and reducing growth rate...']);

            if(isOutputOptimizationData)
                console.info('Logging relevant information... ');
                tic
                
                optimizationData.IterationData(iLog) = struct(... 
                    ... % system data
                    'EvaluatedDesignSamples',designSample,...
                    'EvaluationOutput',outputEvaluation,...
                    ... % algorithm data
                    'Phase',1,...
                    'GrowthRate',growthRate,...
                    ... % problem data
                    'DesignScore',score,...
                    'IsGoodPerformance',isGoodPerformance,...
                    'IsPhysicallyFeasible',isPhysicallyFeasible,...
                    'IsAcceptable',isAcceptable,...
                    'IsUseful',isUseful,...
                    ... % trimming data
                    'CandidateBoxBeforeTrim',candidateBoxGrown,...
                    'CandidateBoxAfterTrim',candidateBox);
                
                % Finalizing
                iLog = iLog + 1;
                console.info('Elapsed time is %g seconds.\n',toc);
            end

            % limit growth rate to this current value
            %options.MaximumGrowthRate = growthRate;

            iExploration = iExploration + 1;
            continue;
        end
        

        %% Modification Step A - Trimming: Remove Bad Points
        % find trimming order
        orderTrim = options.TrimmingOrderFunction(~isAcceptable,score,options.TrimmingOrderOptions{:});
        [candidateBoxTrimmed,measureTrimmed] = sso_box_sub_trimming_operation(...
            candidateBoxGrown,...
            measureGrown,...
            designSample,...
            isAcceptable,...
            isUseful,...
            orderTrim,...
            options.TrimmingOperationFunction,...
            trimmingOperationOptions,...
            console);

        if(applyLeannessEachTrim) 
            trimmingOrder = trimming_order(~isUseful,score,'OrderPreference','score-low-to-high');
            candidateBoxTrimmed = box_trimming_leanness(designSample,isUseful,trimmingOrder,candidateBoxTrimmed);
            % get estimate of new measure
            insideBoxNew = is_in_design_box(designSample,candidateBoxTrimmed);
            measureTrimmed = sso_box_sub_compute_candidate_box_measure(...
                candidateBox,...
                sum(insideBoxNew),...
                sum(insideBoxNew & isUseful),...
                options.MeasureFunction,...
                options.MeasureOptions,...
                console);
        end
        
        %% Convergence Criteria
        console.info('Checking convergence... ');
        tic
        
        % Stop phase I if measure doesn't change significantly from step to step
        if(iExploration >= options.MaxIterExploration)
            hasConvergedExploration = true;
        elseif (~options.FixIterNumberExploration && abs(measureTrimmed-measurePrevious)/measureTrimmed < options.ToleranceMeasureChangeExploration)
            hasConvergedExploration = true;
        end
        
        console.info('Elapsed time is %g seconds.\n',toc);  
        
        if(isOutputOptimizationData)
            %% Save Data
            console.info('Logging relevant information... ');
            tic

            optimizationData.IterationData(iLog) = struct(...
                ... % system data
                'EvaluatedDesignSamples',designSample,...
                'EvaluationOutput',outputEvaluation,...
                ... % algorithm data
                'Phase',1,...
                'GrowthRate',growthRate,...
                ... % problem data
                'DesignScore',score,...
                'IsGoodPerformance',isGoodPerformance,...
                'IsPhysicallyFeasible',isPhysicallyFeasible,...
                'IsAcceptable',isAcceptable,...
                'IsUseful',isUseful,...
                ... % trimming data
                'CandidateBoxBeforeTrim',candidateBoxGrown,...
                'CandidateBoxAfterTrim',candidateBoxTrimmed);
            
            % Finalizing
            iLog = iLog + 1;
            
            console.info('Elapsed time is %g seconds.\n',toc);
        end
        
        
        %% Prepare for next iteration
        console.info('Done with iteration %d!\n\n',iExploration);
    
        candidateBox = candidateBoxTrimmed;
        iExploration = iExploration + 1;
    end
    console.info('\nDone with Phase I - Exploration in iteration %d!\n\n',iExploration-1);

    
    
    %% Intermediary: Sample inside current (End of Exporation) candidate box
    console.info([repelem('=',120),'\n']);
    console.info('Initiating transition to Phase II - Consolidation...\n');
    
    
    %% Phase II - Consolidation
    % Iteration start
    if(options.FixIterNumberConsolidation && (options.MaxIterConsolidation==0))
        convergenceConsolidation = true;
    else
        convergenceConsolidation = false;
    end
    iConsolidation = 1;
    while((~convergenceConsolidation) && (iConsolidation<=options.MaxIterConsolidation))
        console.info([repelem('=',120),'\n']);
        console.info('Initiating Phase II - Consolidation: Iteration %d...\n',iConsolidation);
        
        % get current number of samples
        nSample = get_current_array_entry(options.NumberSamplesPerIterationConsolidation,iConsolidation);
        
        % Generate samples that are to be evaluated
        designSample = sso_box_sub_generate_new_sample_points(...
            candidateBox,...
            nSample,...
            options.SamplingMethodFunction,...
            options.SamplingMethodOptions,...
            console);

        % Evaluate the samples
        [isGoodPerformance,isPhysicallyFeasible,score,outputEvaluation] = sso_box_sub_evaluate_sample_points(...
            designEvaluator,...
            designSample,...
            console);

        % Label samples according to desired requirement spaces problem type
        [isAcceptable,isUseful] = sso_box_sub_label_samples_requirement_spaces(...
            options.RequirementSpacesType,...
            isGoodPerformance,...
            isPhysicallyFeasible,...
            console);

        % Count number of labels
        [nAcceptable,nUseful,nAcceptableUseful] = sso_box_sub_count_label_acceptable_useful(isAcceptable,isUseful,console);

        % No viable design found; throw error
        if(nAcceptable==0 || nUseful==0 || nAcceptableUseful==0)
            console.critical('SSOOptBox:BadSampling',['No good/useful points found, ',...
                'please retry process with different parameters / looser requirements.']);
        end

        % Compute candidate box measure
        measure = sso_box_sub_compute_candidate_box_measure(...
            candidateBox,...
            nSample,...
            nUseful,...
            options.MeasureFunction,...
            options.MeasureOptions,...
            console);
        
        
        %% Convergence Check - Purity
        if(~options.FixIterNumberConsolidation && nAcceptable/nSample>=tolerancePurityConsolidation)
            convergenceConsolidation = true;
        end
        
        
        %% Modification Step A (Trimming): Remove Bad Points
        if(~convergenceConsolidation)
            orderTrim = options.TrimmingOrderFunction(~isAcceptable,score,options.TrimmingOrderOptions{:});
            [candidateBoxTrimmed,measureTrimmed] = sso_box_sub_trimming_operation(...
                candidateBox,...
                measure,...
                designSample,...
                isAcceptable,...
                isUseful,...
                orderTrim,...
                options.TrimmingOperationFunction,...
                options.TrimmingOperationOptions,...
                console);
        end

        if(applyLeannessEachTrim) 
            trimmingOrder = trimming_order(~isUseful,score,'OrderPreference','score-low-to-high');
            candidateBoxTrimmed = box_trimming_leanness(designSample,isUseful,trimmingOrder,candidateBoxTrimmed);
            % get estimate of new measure
            insideBoxNew = is_in_design_box(designSample,candidateBoxTrimmed);
            measureTrimmed = sso_box_sub_compute_candidate_box_measure(...
                candidateBox,...
                sum(insideBoxNew),...
                sum(insideBoxNew & isUseful),...
                options.MeasureFunction,...
                options.MeasureOptions,...
                console);
        end
        
                                                                
        %% Convergence check - Number of Iterations
        if(iConsolidation >= options.MaxIterConsolidation)
            convergenceConsolidation = true;
        end
        
        if(isOutputOptimizationData)
            %% Save Data
            console.info('Logging relevant information... ');
            tic

            optimizationData.IterationData(iLog) = struct(...
                ... % system data
                'EvaluatedDesignSamples',designSample,...
                'EvaluationOutput',outputEvaluation,...
                ... % algorithm data
                'Phase',2,...
                'GrowthRate',[],...
                ... % problem data
                'DesignScore',score,...
                'IsGoodPerformance',isGoodPerformance,...
                'IsPhysicallyFeasible',isPhysicallyFeasible,...
                'IsAcceptable',isAcceptable,...
                'IsUseful',isUseful,...
                ... % trimming data
                'CandidateBoxBeforeTrim',candidateBox,...
                'CandidateBoxAfterTrim',candidateBoxTrimmed);
            
            % Finalizing
            iLog = iLog + 1;
            console.info('Elapsed time is %g seconds.\n',toc);
        end
        
        
        %% Prepare for next iteration
        console.info('Done with iteration %d!\n\n',iConsolidation);
        candidateBox = candidateBoxTrimmed;
        iConsolidation = iConsolidation + 1;
    end
    console.info('\nDone with Phase II - Consolidation in iteration %d!\n\n',iConsolidation-1);
    
    
    %% Check for the leanness condition
    if(applyLeannessFinalTrim)
        trimmingOrder = trimming_order(~isUseful,score,'OrderPreference','score-low-to-high');
        candidateBox = box_trimming_leanness(designSample,isUseful,trimmingOrder,candidateBox);
    end
end


%% Generate New Samples
function designSample = sso_box_sub_generate_new_sample_points(candidateBox,nSample,samplingFunction,samplingOptions,console)
    console.info('Generating new sample points in candidate box... ');
    tic
    
    designSample = samplingFunction(candidateBox,nSample,samplingOptions{:});
    
    console.info('Elapsed time is %g seconds.\n',toc);
    console.debug('- Number of samples generated: %g\n',nSample);
end


%% Evaluate Samples
function [isGoodPerformance,isPhysicallyFeasible,score,outputEvaluation] = sso_box_sub_evaluate_sample_points(designEvaluator,designSample,console)
    console.info('Evaluating sample points... ');
    tic
    
    [performanceDeficit,physicalFeasibilityDeficit,outputEvaluation] = designEvaluator.evaluate(designSample);
    
    [isGoodPerformance,score] = design_deficit_to_label_score(performanceDeficit);
    
    if(~isempty(physicalFeasibilityDeficit))
        isPhysicallyFeasible = design_deficit_to_label_score(physicalFeasibilityDeficit);
    else
        isPhysicallyFeasible = true(size(designSample,1),1);
    end

    console.info('Elapsed time is %g seconds.\n',toc);
    nSamples = size(designSample,1);
    console.debug('- Number of good samples: %g (%g%%)\n',sum(isGoodPerformance),100*sum(isGoodPerformance)/nSamples);
    console.debug('- Number of bad samples: %g (%g%%)\n',sum(~isGoodPerformance),100*sum(~isGoodPerformance)/nSamples);
    console.debug('- Number of physically feasible samples: %g (%g%%)\n',sum(isPhysicallyFeasible),100*sum(isPhysicallyFeasible)/nSamples);
    console.debug('- Number of physically infeasible samples: %g (%g%%)\n',sum(~isPhysicallyFeasible),100*sum(~isPhysicallyFeasible)/nSamples);
end


%% Label Samples
function [isAcceptable,isUseful] = sso_box_sub_label_samples_requirement_spaces(requirementSpacesType,isGoodPerformance,isPhysicallyFeasible,console)
    console.info('Creating labels for each design... ');
    tic
    
    [isAcceptable,isUseful] = design_requirement_spaces_label(requirementSpacesType,isGoodPerformance,isPhysicallyFeasible);
    
    console.info('Elapsed time is %g seconds.\n',toc);
    nSample = size(isGoodPerformance,1);
    console.debug('- Number of accepted samples: %g (%g%%)\n',sum(isAcceptable),100*sum(isAcceptable)/nSample);
    console.debug('- Number of rejected samples: %g (%g%%)\n',sum(~isAcceptable),100*sum(~isAcceptable)/nSample);
    console.debug('- Number of useful samples: %g (%g%%)\n',sum(isUseful),100*sum(isUseful)/nSample);
    console.debug('- Number of useless samples: %g (%g%%)\n',sum(~isUseful),100*sum(~isUseful)/nSample);
end


%% Count Labels
function [nAcceptable,nUseful,nAcceptableUseful] = sso_box_sub_count_label_acceptable_useful(isAcceptable,isUseful,console)
    nAcceptable = sum(isAcceptable);
    nUseful = sum(isUseful);
    nAcceptableUseful = sum(isAcceptable & isUseful);
    
    nSamples = size(isAcceptable,1);
    console.debug('- Number of accepted samples: %g (%g%%)\n',nAcceptable,100*nAcceptable/nSamples);
    console.debug('- Number of useful samples: %g (%g%%)\n',nUseful,100*nUseful/nSamples);
    console.debug('- Number of accepted and useful samples: %g (%g%%)\n',nAcceptableUseful,100*nAcceptableUseful/nSamples);
end


%% Compute Measure
function measure = sso_box_sub_compute_candidate_box_measure(candidateBox,nSample,nUseful,measureFunction,measureOptions,console)
    console.info('Computing candidate box measure... ');
    tic

    measure = measureFunction(candidateBox,nUseful/nSample,measureOptions{:});
    
    console.info('Elapsed time is %g seconds.\n',toc);
    console.debug('- Current candidate box measure: %g\n',measure);
end

%% Trimming Operation
function [candidateBoxTrimmed,measureTrimmed] = sso_box_sub_trimming_operation(candidateBox,measure,designSample,isAcceptable,isUseful,orderTrim,trimmingMethodFunction,trimmingOperationOptions,console)
    console.info('Performing box trimming operation... ');
    tic
    
    if(sum(isAcceptable)~=size(isAcceptable,1)) 
        % if there are bad designs, perform trimming
        labelViable = isAcceptable & isUseful;
        [candidateBoxTrimmed,measureTrimmed] = trimmingMethodFunction(designSample,labelViable,orderTrim,candidateBox,trimmingOperationOptions{:});
    else
        % no trimming necessary
        candidateBoxTrimmed = candidateBox;
        measureTrimmed = measure;
    end
    
    console.info('Elapsed time is %g seconds.\n',toc);
    insideBox = is_in_design_box(designSample,candidateBoxTrimmed);
    acceptedInsideBox = insideBox & isAcceptable;
    usefulInsideBox = insideBox & isUseful;
    nSample = size(designSample,1);
    console.debug('- Number of samples removed from candidate box: %g (%g%%)\n',sum(~insideBox),100*sum(~insideBox)/nSample);
    console.debug('- Number of acceptable samples lost: %g\n',sum(~insideBox & isAcceptable));
    console.debug('- Number of useful samples lost: %g\n',sum(~insideBox & isUseful));
    console.debug('- Number of acceptable samples inside trimmed candidate box: %g\n',sum(acceptedInsideBox));
    console.debug('- Number of useful samples inside trimmed candidate box: %g\n',sum(usefulInsideBox));
    console.debug('- Trimmed candidate box measure: %g (Relative shrinkage: %g%%)\n',measureTrimmed,-100*(measureTrimmed-measure)/measure);
end

