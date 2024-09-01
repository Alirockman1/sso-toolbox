function [componentSolutionSpace,problemData,iterationData] = sso_component_stochastic(designEvaluator,initialDesign,designSpaceLowerBound,designSpaceUpperBound,componentIndex,varargin)
%SSO_COMPONENT_STOCHASTIC Component solution spaces optimization (stochastic)
%   SSO_COMPONENT_STOCHASTIC computes the optimal component solution spaces 
%   using a stochastic approach, where samples are generated inside the 
%   candidate spaces and are used with the trimming operation to determine what
%   will be the new inside / outside of each space.
%
%   COMPONENTSOLUTIONSPACE = SSO_COMPONENT_STOCHASTIC(DESIGNEVALUATOR,
%   INITIALDESIGN,DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND,COMPONENTINDEX) 
%   starts the algorithm in INITIALDESIGN and finds the optimum component 
%   solution (or requirement) spaces within the design space defined by 
%   DESIGNSPACELOWERBOUND and DESIGNSPACEUPPERBOUND, evaluating design sample 
%   points with DESIGNEVALUATOR and returning the optimal component solution (or 
%   requirement) spaces COMPONENTSOLUTIONSPACE.
%
%   COMPONENTSOLUTIONSPACE = SSO_COMPONENT_STOCHASTIC(DESIGNEVALUATOR,
%   INITIALDESIGN,DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND,COMPONENTINDEX,
%   OPTIONS) also allows one to change the options being used through the output
%   of 'sso_stochastic_options' OPTIONS.
%
%   [COMPONENTSOLUTIONSPACE,PROBLEMDATA] = SSO_COMPONENT_STOCHASTIC(...)  
%   additionally returns the fixed problem data in PROBLEMDATA.
%
%   [COMPONENTSOLUTIONSPACE,PROBLEMDATA,ITERATIONDATA] = 
%   SSO_COMPONENT_STOCHASTIC(...) additionally returns data generated at each 
%   iteration of the process ITERATIONDATA.
%
%
%   Input:
%       - DESIGNEVALUATOR : DesignEvaluatorBase
%       - INITIALDESIGN : (1,nDesignVariable) double OR (2,nDesignVariable) 
%       double
%       - DESIGNSPACELOWERBOUND : (1,nDesignVariable) double
%       - DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double
%       - COMPONENTINDEX : (1,nComponent) cell
%       - OPTIONS : structure
%
%   Output:
%       - CANDIDATEBOX  : (2,nDesignVariable) double
%       - PROBLEMDATA : struct
%           - DesignEvaluator : DesignEvaluatorBase 
%           - InitialBox : (2,nDesignVariable) double
%           - DesignSpaceLowerBound : (1,nDesignVariable) double
%           - DesignSpaceUpperBound : (1,nDesignVariable) double
%           - Options : struct
%           - InitialRNGState : struct
%       - ITERATIONDATA : (nIteration,1) struct
%           - EvaluatedDesignSamples : (nSample,nDesignVariable) double
%           - EvaluationOutput : class-dependent
%           - Phase : integer
%           - GrowthRate : double
%           - NumberPaddingSamplesGenerated : double
%           - PaddingSamplesUsed : (nPaddingSample,nDesignVariable) double
%           - DesignScore : (nSample,1) double
%           - IsGoodPerformance : (nSample,1) logical
%           - IsPhysicallyFeasible : (nSample,1) logical
%           - IsAcceptable : (nSample,1) logical
%           - isUseful : (nSample,1) logical
%           - SamplingBoxBeforeTrim : (2,nDesignVariable) double
%           - SamplingBoxAfterTrim : (2,nDesignVariable) double
%           - CandidateSpacesBeforeTrim : (1,nComponent) CandidateSpaceBase
%           - CandidateSpacesAfterTrim : (1,nComponent) CandidateSpaceBase
%
%   See also sso_stochastic_options.
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
    
    % Options
    defaultOptions = sso_stochastic_options('component');
    inputOptions = parser_variable_input_to_structure(varargin{:});
    options = merge_name_value_pair_argument(defaultOptions,inputOptions);
    
    
    %% extract options as necessary
    % requirement spaces
    applyLeannessEachTrim = strcmpi(options.ApplyLeanness,'always');
    applyLeannessFinalTrim = any(strcmpi(options.ApplyLeanness,{'always','end-only'}));

    % sampling
    [~,candidateSpaceSamplingOptions] = merge_name_value_pair_argument(...
        structure_extract_fields(options,...
            {'SamplingMethodFunction';...
                'SamplingMethodOptions'}),...
        options.CandidateSpaceSamplingOptions);
    
    % trimming
    [~,trimmingOperationOptions] = merge_name_value_pair_argument(...
        structure_extract_fields(options,...
            {'TrimmingCostFunction';...
                'TrimmingCostOptions';...
                'TrimmingComponentChoiceFunction';...
                'TrimmingComponentChoiceOptions';...
                'TrimmingMethodFunction';...
                'TrimmingMethodOptions'}),...
        options.TrimmingOperationOptions);
    
    % verbosity
    console = ConsoleLogging(options.LoggingLevel);
    
    
    %% Initial Setup
    % Initial Candidate Space
    if(size(initialDesign,1)==1)
        samplingBox = [initialDesign;initialDesign]; % single point
    elseif(size(initialDesign,1)==2)
        samplingBox = initialDesign; % candidate box
    else
        console.error('SSOBoxOptStochastic:InitialDesignWrong','Error. Initial guess/region incompatible in ''sso_component_stochastic''.');
    end
    
                             
    %% Log Initialization
    if(nargout>=2)
        problemData = struct(...
            'DesignClassifier',designEvaluator,...
            'InitialDesign',samplingBox,...
            'DesignSpaceLowerBound',designSpaceLowerBound,...
            'DesignSpaceUpperBound',designSpaceUpperBound,...
            'ComponentIndex',{componentIndex},...
            'Options',options,...
            'InitialRNGState',rng);
    end

    isOutputIterationData = (nargout>=3);
    if(isOutputIterationData)
        iterationData = struct(... 
            ... % system data
        	'EvaluatedDesignSamples',[],...
            'EvaluationOutput',[],...
            ... % algorithm data
            'Phase',[],...
            'GrowthRate',[],...
            'NumberPaddingSamplesGenerated',[],...
            'PaddingSamplesUsed',[],...
            ... % problem data
            'DesignScore',[],...
            'IsGoodPerformance',[],...
            'IsPhysicallyFeasible',[],...
            'IsAcceptable',[],...
            'IsUseful',[],...
            ... % trimming data
            'SamplingBoxBeforeTrim',[],...
            'SamplingBoxAfterTrim',[],... 
            'CandidateSpacesBeforeTrim',[],...
            'CandidateSpacesAfterTrim',[]);
        iLog = 1;
    end
    
    
    %% Exploration: While mu(Omega) is changing, grow box
    % Create first instantiation of the candidate spaces for each component
    for i=1:size(componentIndex,2)
        candidateSpace(i) = options.CandidateSpaceConstructorExploration(...
            designSpaceLowerBound(componentIndex{i}),...
            designSpaceUpperBound(componentIndex{i}),...
            options.CandidateSpaceOptionsExploration{:});
        candidateSpaceDefined = false;
    end
    
    % initialize iterators and loop
    iExploration = 1;
    growthRate = options.GrowthRate;
    convergedExploration = false;
    while((~convergedExploration) && (iExploration<=options.MaxIterExploration))
        console.info([repelem('=',120),'\n']);
        console.info('Initiating Phase I - Exploration: Iteration %d\n',iExploration);

        if(options.UseAdaptiveGrowthRate && iExploration>1)
            console.info('Adapting growth rate... ');
            tic

            % Change step size to a bigger or smaller value depending on whether
            % there were or weren't enough good design points in the step
            growthRate = (1/options.TargetAcceptedRatioExploration)*(nAcceptable/nSample)*growthRate;
            growthRate = max(growthRate,options.MinimumGrowthRate);

            console.info('Elapsed time is %g seconds.\n',toc);
        end

        %% Modification Step B - Growth: Extend Candidate Space
        console.info('Growing candidate space... ');
        tic

        if(~candidateSpaceDefined)
            % Expand candidate solution box in both sides of each interval isotroply
            samplingBoxGrown = design_box_grow_fixed(samplingBox,designSpaceLowerBound,designSpaceUpperBound,growthRate);
            candidateSpaceGrown = candidateSpace;
        else
            for i=1:size(componentIndex,2)
                candidateSpaceGrown(i) = candidateSpace(i).grow_candidate_space(growthRate);
                samplingBoxGrown(:,componentIndex{i}) = candidateSpaceGrown(i).SamplingBox;
            end
        end
        
        console.info('Elapsed time is %g seconds.\n',toc);
        console.debug('- Current Growth Rate: %g\n',growthRate);
        
        
        %% Sample inside the current candidate box
        % get current number of samples
        nSample = get_current_array_entry(options.NumberSamplesPerIterationExploration,iExploration);
        
        % Generate samples that are to be evaluated
        console.info('Generating new sample points in candidate space... ');
        tic
        
        if(~candidateSpaceDefined)
            designSample = component_sso_generate_new_sample_points_sampling_box(...
                samplingBoxGrown,...
                nSample,...
                options.SamplingMethodFunction,...
                options.SamplingMethodOptions,...
                console);
            paddingSample = [];
        else
            [designSample,paddingSample] = options.CandidateSpaceSamplingFunction(...
                candidateSpaceGrown,...
                componentIndex,...
                nSample,...
                candidateSpaceSamplingOptions{:});
        end
        
        % truncate padding samples if there are too many, or generate more
        nPaddingGenerated = size(paddingSample,1);
        if(nPaddingGenerated>=options.NumberPaddingSamples)
            paddingSample = paddingSample(1:options.NumberPaddingSamples,:);
        else
            nPaddingMissing = options.NumberPaddingSamples - nPaddingGenerated;
            extraPadding = options.SamplingMethodFunction(samplingBoxGrown,nPaddingMissing,options.SamplingMethodOptions{:});
            paddingSample = [paddingSample;extraPadding];
        end


        console.info('Elapsed time is %g seconds.\n',toc);
        console.debug('- Number of samples generated: %g\n',size(designSample,1));
        console.debug('- Number of padding samples: %g\n',size(paddingSample,1));
            
        % Evaluate the samples
        [isGoodPerformance,isPhysicallyFeasible,score,outputEvaluation] = ...
            component_sso_evaluate_sample_points(designEvaluator,designSample,console);
        
        % Label samples according to desired requirement spaces problem type
        [isAcceptable,isUseful] = component_sso_label_requirement_spaces(...
            options.RequirementSpacesType,...
            isGoodPerformance,...
            isPhysicallyFeasible,...
            console);
        
        % Count number of labels
        [nAcceptable,nUseful,nAccetableUseful] = component_sso_count_acceptable_useful(isAcceptable,isUseful,console);
            
        % Box may have grown too much + "unlucky sampling" w/ no good
        % points, go back in this case
        if(nAcceptable==0 || nUseful==0 || nAccetableUseful==0)
            console.warn('SSOOptBox:BadSampling',['No good points found, ',...
                'rolling back and reducing growth rate to minimum...']);
            
            if(isOutputIterationData)
                console.info('Logging relevant information... ');
                tic
                
                iterationData(iLog) = struct(... 
                    ... % system data
                    'EvaluatedDesignSamples',designSample,...
                    'EvaluationOutput',outputEvaluation,...
                    ... % algorithm data
                    'Phase',1,...
                    'GrowthRate',growthRate,...
                    'NumberPaddingSamplesGenerated',nPaddingGenerated,...
                    'PaddingSamplesUsed',paddingSample,...
                    ... % problem data
                    'DesignScore',score,...
                    'IsGoodPerformance',isGoodPerformance,...
                    'IsPhysicallyFeasible',isPhysicallyFeasible,...
                    'IsAcceptable',isAcceptable,...
                    'IsUseful',isUseful,...
                    ... % trimming data
                    'SamplingBoxBeforeTrim',samplingBoxGrown,...
                    'SamplingBoxAfterTrim',[],... 
                    'CandidateSpacesBeforeTrim',candidateSpaceGrown,...
                    'CandidateSpacesAfterTrim',[]);
                
                % Finalizing
                iLog = iLog + 1;
                console.info('Elapsed time is %g seconds.\n',toc);
            end

            % Rollback Settings
            growthRate = options.MinimumGrowthRate;

            % proceed to next iteration
            iExploration = iExploration + 1;
            continue;
        end
        

        %% Modification Step A - Trimming: Remove Bad Points
        console.info('Performing component trimming operation... ');
        tic
        
        isPadding = false(size(designSample,1),1);
        [trimmingSample,trimmingIsAcceptable,trimmingIsUseful,trimmingScore,isPadding] = component_sso_prepare_trimming_samples(...
            designSample,...
            isAcceptable,...
            isUseful,...
            score,...
            isPadding,...
            paddingSample,...
            options.UsePaddingSamplesInTrimming);
        
        % define which samples are inside the candidate space
        activeComponent = true(size(trimmingSample,1),size(componentIndex,2));
        if(candidateSpaceDefined)
            for i=1:size(componentIndex,2)
                activeComponent(:,i) = candidateSpaceGrown(i).is_in_candidate_space(trimmingSample(:,componentIndex{i}));
            end
        end

        % define order of trimming operation for samples that must be excluded
        trimmingOrder = options.TrimmingOrderFunction(~trimmingIsAcceptable,trimmingScore,options.TrimmingOrderOptions{:});

        % trim
        trimmingLabelViable = trimmingIsAcceptable & trimmingIsUseful;
        activeComponent = component_trimming_operation(trimmingSample,trimmingLabelViable,trimmingOrder,componentIndex,activeComponent,trimmingOperationOptions{:});
        active = all(activeComponent,2);
        console.info('Elapsed time is %g seconds.\n',toc);

        if(applyLeannessEachTrim)
            console.info('Applying the leanness condition... ');

            labelRemoveLeanness = ~trimmingIsUseful & ~isPadding;
            labelKeep = trimmingIsAcceptable & trimmingIsUseful;
            trimmingOrder = trimming_order(labelRemoveLeanness,trimmingScore,'OrderPreference','score-low-to-high');
            activeComponent = component_trimming_leanness(trimmingSample,labelKeep,trimmingOrder,componentIndex,activeComponent,trimmingOperationOptions{:});
            active = all(activeComponent,2);

            console.info('Elapsed time is %g seconds.\n',toc);
        end
        
        acceptedInsideSpace = active & trimmingIsAcceptable;
        usefulInsideSpace = active & trimmingIsUseful;
        % console.debug('- Number of samples removed from candidate space: %g (%g%%)\n',sum(~active),100*sum(~active)/size(active,2));
        % console.debug('- Number of acceptable samples lost: %g\n',sum(~active & trimmingIsAcceptable));
        % console.debug('- Number of useful samples lost: %g\n',sum(~active & trimmingIsUseful));
        % console.debug('- Number of acceptable samples inside trimmed candidate space: %g\n',sum(acceptedInsideSpace));
        % console.debug('- Number of useful samples inside trimmed candidate space: %g\n',sum(usefulInsideSpace));
        
        
        %% Retraining Candidate Spaces
        console.info('Updating candidate spaces... ');
        tic
        [candidateSpaceTrimmed,samplingBoxTrimmed] = component_sso_update_candidate_spaces(trimmingSample,activeComponent,componentIndex,candidateSpace);
        candidateSpaceDefined = true;
        console.info('Elapsed time is %g seconds.\n',toc);
        
        
        %% Convergence Criteria
        console.info('Checking convergence... ');
        tic
        convergenceCriterion = all(abs((samplingBoxTrimmed-samplingBox)./samplingBox)<=options.ToleranceMeasureChangeExploration,'all');
        
        % Stop phase I if mu doesn't change significantly from step to step
        if(~options.FixIterNumberExploration && convergenceCriterion)
            convergedExploration = true;
        end
        console.info('Elapsed time is %g seconds.\n',toc);
        
        
        %% Save Data
        if(isOutputIterationData)
            console.info('Logging relevant information... ');
            tic
            
            iterationData(iLog) = struct(... 
                ... % system data
                'EvaluatedDesignSamples',designSample,...
                'EvaluationOutput',outputEvaluation,...
                ... % algorithm data
                'Phase',1,...
                'GrowthRate',growthRate,...
                'NumberPaddingSamplesGenerated',nPaddingGenerated,...
                'PaddingSamplesUsed',paddingSample,...
                ... % problem data
                'DesignScore',score,...
                'IsGoodPerformance',isGoodPerformance,...
                'IsPhysicallyFeasible',isPhysicallyFeasible,...
                'IsAcceptable',isAcceptable,...
                'IsUseful',isUseful,...
                ... % trimming data
                'SamplingBoxBeforeTrim',samplingBoxGrown,...
                'SamplingBoxAfterTrim',samplingBoxTrimmed,... 
                'CandidateSpacesBeforeTrim',candidateSpaceGrown,...
                'CandidateSpacesAfterTrim',candidateSpaceTrimmed);
            
            % Finalizing
            iLog = iLog + 1;
            console.info('Elapsed time is %g seconds.\n',toc);
        end
        
        
        %% Prepare for next iteration
        console.info('Done with iteration %d!\n\n',iExploration);
        samplingBox = samplingBoxTrimmed;
        candidateSpace = candidateSpaceTrimmed;
        iExploration = iExploration + 1;
    end
    console.info('\nDone with Phase I - Exploration in iteration %d!\n\n',iExploration-1);
    
    
    %% Phase II - Consolidation
    console.info([repelem('=',120),'\n']);
    console.info('Initiating Phase II - Consolidation...\n');
    
    % Create first instantiation of the candidate spaces for each component
    % in addition, find which designs help define the shape of each space.
    isShapeDefinition = false(size(activeComponent));
    for i=1:size(componentIndex,2)
        candidateSpaceConsolidation(i) = options.CandidateSpaceConstructorConsolidation(...
            designSpaceLowerBound(componentIndex{i}),...
            designSpaceUpperBound(componentIndex{i}),...
            options.CandidateSpaceOptionsConsolidation{:});
        
        candidateSpaceConsolidation(i) = candidateSpaceConsolidation(i).define_candidate_space(...
            candidateSpace(i).DesignSampleDefinition,candidateSpace(i).IsInsideDefinition);

        samplingBox(:,componentIndex{i}) = candidateSpaceConsolidation(i).SamplingBox;
        isShapeDefinition(:,i) = candidateSpaceConsolidation(i).IsShapeDefinition;
    end
    candidateSpace = candidateSpaceConsolidation;
    clear candidateSpaceConsolidation
    
    % create new arrays for designs to be used/cconsidered
    isKept = false(size(activeComponent,1),3);
    isKept(:,1) = any(isShapeDefinition,2);
    isKept(:,2) = (~isPadding) & options.UsePreviousEvaluatedSamplesConsolidation;
    isKept(:,3) = (isPadding) & options.UsePreviousPaddingSamplesConsolidation;
    keepSample = any(isKept,2);

    previousTrimmingSample = trimmingSample(keepSample,:);
    previousTrimmingisAcceptable = trimmingIsAcceptable(keepSample,:);
    previousTrimmingisUseful = trimmingIsUseful(keepSample,:);
    previousTrimmingScore = trimmingScore(keepSample,:);
    previousIsPadding = isPadding(keepSample,:);
    iConsolidation = 1;
    convergedConsolidation = false;
    while((~convergedConsolidation) && (iConsolidation<=options.MaxIterConsolidation))
        console.info([repelem('=',120),'\n']);
        console.info('Initiating Phase II - Consolidation: Iteration %d\n',iConsolidation);


        %% Sample inside the current candidate space
        % get current number of samples
        nSample = get_current_array_entry(options.NumberSamplesPerIterationConsolidation,iConsolidation);
        
        console.info('Generating new sample points in candidate space... ');
        tic

        [designSample,paddingSample] = options.CandidateSpaceSamplingFunction(...
                candidateSpace,...
                componentIndex,...
                nSample,...
                candidateSpaceSamplingOptions{:});

        % truncate padding samples if there are too many, or generate more
        nPaddingGenerated = size(paddingSample,1);
        if(nPaddingGenerated>=options.NumberPaddingSamples)
            paddingSample = paddingSample(1:options.NumberPaddingSamples,:);
        else
            nPaddingMissing = options.NumberPaddingSamples - nPaddingGenerated;
            extraPadding = options.SamplingMethodFunction(samplingBox,nPaddingMissing,options.SamplingMethodOptions{:});
            paddingSample = [paddingSample;extraPadding];
        end
        
        console.info('Elapsed time is %g seconds.\n',toc);
        console.debug('- Number of samples generated: %g\n',size(designSample,1));
        console.debug('- Number of padding samples: %g\n',size(paddingSample,1));
        
        % Evaluate the samples
        [isGoodPerformance,isPhysicallyFeasible,score,outputEvaluation] = ...
            component_sso_evaluate_sample_points(designEvaluator,designSample,console);
        
        % Label samples according to desired requirement spaces problem type
        [isAcceptable,isUseful] = component_sso_label_requirement_spaces(...
            options.RequirementSpacesType,isGoodPerformance,isPhysicallyFeasible,console);
        
        % Count number of labels
        [nAcceptable,nUseful,nAcceptableUseful] = component_sso_count_acceptable_useful(isAcceptable,isUseful,console);
        
        
        %% Convergence Criteria - Purity
        samplePurity = nAcceptable/nSample;
        puritySatisfied = (samplePurity>=options.TolerancePurityConsolidation);
        if(~options.MaxIterConsolidation && puritySatisfied)
            convergedConsolidation = true;
        end
        

        %% Modification Step A - Trimming: Remove Bad Points
        performTrim = (~convergedConsolidation && nAcceptable<nSample && ~puritySatisfied);
        if(performTrim)
            console.info('Performing component trimming operation... ');
            tic
            consideredSample = [previousTrimmingSample;designSample];
            consideredisAcceptable = [previousTrimmingisAcceptable;isAcceptable];
            consideredisUseful = [previousTrimmingisUseful;isUseful];
            consideredScore = [previousTrimmingScore;score];
            consideredIsPadding = [previousIsPadding;false(size(designSample,1),1)];

            [trimmingSample,trimmingIsAcceptable,trimmingIsUseful,trimmingScore,trimmingIsPadding] = ...
                component_sso_prepare_trimming_samples(...
                consideredSample,...
                consideredisAcceptable,...
                consideredisUseful,...
                consideredScore,...
                consideredIsPadding,...
                paddingSample,...
                options.UsePaddingSamplesInTrimming);

            
            % define which samples are inside the candidate spaces
            activeComponent = true(size(trimmingSample,1),size(componentIndex,2));
            for i=1:size(componentIndex,2)
                activeComponent(:,i) = candidateSpace(i).is_in_candidate_space(trimmingSample(:,componentIndex{i}));
            end

            % define order of trimming operation for samples that must be excluded
            trimmingOrder = options.TrimmingOrderFunction(~trimmingIsAcceptable,trimmingScore,options.TrimmingOrderOptions{:});

            % trim
            trimmingLabelViable = trimmingIsAcceptable & trimmingIsUseful;
            activeComponent = component_trimming_operation(trimmingSample,trimmingLabelViable,trimmingOrder,componentIndex,activeComponent,trimmingOperationOptions{:});
            active = all(activeComponent,2);
            console.info('Elapsed time is %g seconds.\n',toc);

            if(applyLeannessEachTrim)
                console.info('Applying the leanness condition... ');

                labelRemoveLeanness = ~trimmingIsUseful & ~trimmingIsPadding;
                labelKeep = trimmingIsAcceptable & trimmingIsUseful;
                trimmingOrder = trimming_order(labelRemoveLeanness,trimmingScore,'OrderPreference','score-low-to-high');
                activeComponent = component_trimming_leanness(trimmingSample,labelKeep,trimmingOrder,componentIndex,activeComponent,trimmingOperationOptions{:});
                active = all(activeComponent,2);

                console.info('Elapsed time is %g seconds.\n',toc);
            end
            
            acceptedInsideSpace = active & trimmingIsAcceptable;
            usefulInsideSpace = active & trimmingIsUseful;
            % console.debug('- Number of samples removed from candidate space: %g (%g%%)\n',sum(~active),100*sum(~active)/size(active,2));
            % console.debug('- Number of acceptable samples lost: %g\n',sum(~active & trimmingIsAcceptable));
            % console.debug('- Number of useful samples lost: %g\n',sum(~active & trimmingIsUseful));
            % console.debug('- Number of acceptable samples inside trimmed candidate space: %g\n',sum(acceptedInsideSpace));
            % console.debug('- Number of useful samples inside trimmed candidate space: %g\n',sum(usefulInsideSpace));
            
            %% Retraining Candidate Spaces
            console.info('Updating candidate spaces... ');
            tic
            
            [candidateSpaceTrimmed,samplingBoxTrimmed,isShapeDefinition] = ...
                component_sso_update_candidate_spaces(trimmingSample,activeComponent,componentIndex,candidateSpace);

            %% update samples being kept
            isKept = false(size(activeComponent,1),3);
            isKept(:,1) = any(isShapeDefinition,2);
            isKept(:,2) = (~trimmingIsPadding) & options.UsePreviousEvaluatedSamplesConsolidation;
            isKept(:,3) = (trimmingIsPadding) & options.UsePreviousPaddingSamplesConsolidation;
            keepSample = any(isKept,2);

            previousTrimmingSample = trimmingSample(keepSample,:);
            previousTrimmingisAcceptable = trimmingIsAcceptable(keepSample,:);
            previousTrimmingisUseful = trimmingIsUseful(keepSample,:);
            previousTrimmingScore = trimmingScore(keepSample,:);
            previousIsPadding = trimmingIsPadding(keepSample,:);
            
            console.info('Elapsed time is %g seconds.\n',toc);
        else
            console.info('Purity criterium satisfied; no trimming is necessary.\n');
            candidateSpaceTrimmed = candidateSpace;
            samplingBoxTrimmed = samplingBox;
        end
        
        %% Convergence check - Number of Iterations
        console.info('Checking convergence... ');
        tic
        
        if(iConsolidation >= options.MaxIterConsolidation)
            convergedConsolidation = true;
        end
        
        console.info('Elapsed time is %g seconds.\n',toc);
        
        
        %% Save Data
        if(isOutputIterationData)
            console.info('Logging relevant information... ');
            tic
            
            iterationData(iLog) = struct(... 
                ... % system data
                'EvaluatedDesignSamples',designSample,...
                'EvaluationOutput',outputEvaluation,...
                ... % algorithm data
                'Phase',2,...
                'GrowthRate',0,...
                'NumberPaddingSamplesGenerated',nPaddingGenerated,...
                'PaddingSamplesUsed',paddingSample,...
                ... % problem data
                'DesignScore',score,...
                'IsGoodPerformance',isGoodPerformance,...
                'IsPhysicallyFeasible',isPhysicallyFeasible,...
                'IsAcceptable',isAcceptable,...
                'IsUseful',isUseful,...
                ... % trimming data
                'SamplingBoxBeforeTrim',samplingBox,...
                'SamplingBoxAfterTrim',samplingBoxTrimmed,... 
                'CandidateSpacesBeforeTrim',candidateSpace,...
                'CandidateSpacesAfterTrim',candidateSpaceTrimmed);
            
            % Finalizing
            iLog = iLog + 1;
            console.info('Elapsed time is %g seconds.\n',toc);
        end
        
        
        %% Prepare for next iteration
        console.info('Done with iteration %d!\n\n',iConsolidation);
        samplingBox = samplingBoxTrimmed;
        candidateSpace = candidateSpaceTrimmed;
        iConsolidation = iConsolidation + 1;
    end
    console.info('\nDone with Phase II - Consolidation!\n\n');
    
    
    %% leanness condition
    if(applyLeannessFinalTrim)
        console.info('Applying the leanness condition... ');
        labelRemoveLeanness = ~trimmingIsUseful & ~trimmingIsPadding;
        labelKeep = trimmingIsAcceptable & trimmingIsUseful;
        trimmingOrder = trimming_order(labelRemoveLeanness,trimmingScore,'OrderPreference','score-low-to-high');
        activeComponent = component_trimming_leanness(trimmingSample,labelKeep,trimmingOrder,componentIndex,activeComponent,trimmingOperationOptions{:});
        console.info('Elapsed time is %g seconds.\n',toc);

        console.info('Updating candidate spaces... ');
        candidateSpace = component_sso_update_candidate_spaces(trimmingSample,activeComponent,componentIndex,candidateSpace);
        console.info('Elapsed time is %g seconds.\n\n',toc);
    end
    
    %% 
    componentSolutionSpace = candidateSpace;
end

%% Generate New Samples
function dv = component_sso_generate_new_sample_points_sampling_box(candidateBox,nSample,samplingMethodFunction,samplingMethodOptions,console)
    console.info('Generating new sample points in candidate box... ');
    tic
    
    dv = samplingMethodFunction(candidateBox,nSample,samplingMethodOptions{:});
    
    toc
    console.debug('- Number of samples generated: %g\n',nSample);
end


%% Evaluate Samples
function [labelReq,labelPf,score,outputEvaluation] = component_sso_evaluate_sample_points(designEvaluator,dv,console)
    console.info('Evaluating sample points... ');
    tic
    
    [performanceDeficit,physicalFeasibilityDeficit,outputEvaluation] = designEvaluator.evaluate(dv);
    [labelReq,score] = design_deficit_to_label_score(performanceDeficit);
    labelPf = design_deficit_to_label_score(physicalFeasibilityDeficit);
    
    toc
    console.debug('- Number of good samples: %g (%g%%)\n',sum(labelReq),100*sum(labelReq)/size(labelReq,1));
    console.debug('- Number of bad samples: %g (%g%%)\n',sum(~labelReq),100*sum(~labelReq)/size(labelReq,1));
    console.debug('- Number of physically feasible samples: %g (%g%%)\n',sum(labelPf),100*sum(labelPf)/size(labelPf,1));
    console.debug('- Number of physically infeasible samples: %g (%g%%)\n',sum(~labelPf),100*sum(~labelPf)/size(labelPf,1));
end


%% Label Samples
function [labelAcc_sample,labelUse_sample] = component_sso_label_requirement_spaces(RequirementSpacesType,labelReq_sample,labelPF_sample,console)
    console.info('Creating labels for each design... ');
    tic
    
    [labelAcc_sample,labelUse_sample] = design_requirement_spaces_label(RequirementSpacesType,labelReq_sample,labelPF_sample);
    
    toc
end


%% Count Labels
function [mAcc,mUse,mAccUse] = component_sso_count_acceptable_useful(labelAcc_sample,labelUse_sample,console)
    mAcc = sum(labelAcc_sample);
    mUse = sum(labelUse_sample);
    mAccUse = sum(labelAcc_sample&labelUse_sample);
    
    console.debug('- Number of accepted samples: %g (%g%%)\n',mAcc,100*mAcc/size(labelAcc_sample,1));
    console.debug('- Number of useful samples: %g (%g%%)\n',mUse,100*mUse/size(labelUse_sample,1));
    console.debug('- Number of accepted and useful samples: %g (%g%%)\n',mAccUse,100*mAccUse/size(labelAcc_sample,1));
end


%%
function [candidateSpaces,samplingBox,isShapeDefinition] = component_sso_update_candidate_spaces(dvTrim,activeComponent,componentIndex,candidateSpaces)
    samplingBox = nan(2,size(dvTrim,2));
    isShapeDefinition = false(size(activeComponent));
    for i=1:size(componentIndex,2)
        % Eliminate from sampling designs that were taken out by other components
        Xtrain = dvTrim(:,componentIndex{i});
        Ytrain = activeComponent(:,i);

        candidateSpaces(i) = candidateSpaces(i).define_candidate_space(Xtrain,Ytrain);
        samplingBox(:,componentIndex{i}) = candidateSpaces(i).SamplingBox;

        isShapeDefinition(:,i) = candidateSpaces(i).IsShapeDefinition;
    end
end


%% 
function [dvTrim,labelAccTrim,labelUseTrim,scoreTrim,isPadding] = component_sso_prepare_trimming_samples(dv,labelAcc,labelUse,score,isPadding,dvPad,usePadInTrimming)
    dvTrim = dv;
    labelAccTrim = labelAcc;
    labelUseTrim = labelUse;
    scoreTrim = score;
    
    if(usePadInTrimming)
        % for trimming purposes, may consider the designs for padding unnaceptable and useless
        dvTrim = [dvTrim;dvPad];
        labelAccTrim = [labelAccTrim;true(size(dvPad,1),1)];
        labelUseTrim = [labelUseTrim;false(size(dvPad,1),1)];
        scoreTrim = [scoreTrim;-ones(size(dvPad,1),1)];
        isPadding = [isPadding;true(size(dvPad,1),1)];
    end
end