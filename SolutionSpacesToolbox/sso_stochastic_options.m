function options = sso_stochastic_options(solutionSpaceType,varargin)
%SSO_STOCHASTIC_OPTIONS Options to be used for stochastic SSO (box/component)
%   SSO_STOCHASTIC_OPTIONS returns a structure that contains all options 
%   required by either 'sso_box_stochastic' or 'sso_component_stochastic' 
%   depending on the choice. 
%
%   OPTIONS = SSO_STOCHASTIC_OPTIONS(SOLUTIONSPACETYPE)
%   creates such a structure OPTIONS with the options necessary for either
%   'sso_box_stochastic' if SOLUTIONSPACETYPE is 'box' or 
%   'sso_component_stochastic' if SOLUTIONSPACETYPE is 'component'.
%
%   OPTIONS = SSO_STOCHASTIC_OPTIONS(...,NAME,VALUE,...) allows one to specify
%   those options by name-value pairs. 
%
%   Common options (available for both box and component SSO):
%       - 'RequirementSpacesType' : type of requirement spaces to be used 
%         during trimming; available options are 'Omega0','Omega1','Omega2', 
%         and 'Omega3'. Default: 'Omega2'.
%
%       - 'ApplyLeanness' : determine when/if the leanness condition of 
%         requirement spaces should be applied. 
%         Available: 'always','end-only','never'. Default: 'always'.
%
%       - 'NumberSamplesPerIteration' : integer, number of samples to be 
%         evaluated per iteration (both in exploration and consolidation 
%         phases) if the phase-specific options are not set. Default: 100.
%
%       - 'SamplingMethodFunction' : function handle used for sampling
%         stochastically inside candidate spaces. Default: 
%         @sampling_latin_hypercube.
%
%       - 'SamplingMethodOptions' : cell array of additional options for 
%         the sampling function. Default: {}.
%
%       - 'NumberSamplesPerIterationExploration' : integer, number of samples 
%         used during the exploration phase. If empty, defaults to 
%         'NumberSamplesPerIteration'.
%
%       - 'NumberSamplesPerIterationConsolidation' : integer, number of samples 
%         used during the consolidation phase. If empty, defaults to 
%         'NumberSamplesPerIteration'.
%
%       - 'GrowthRate' : double, fixed growth rate of the box/candidate spaces 
%         during exploration if adaptive growth is not used. Default: 0.1.
%
%       - 'UseAdaptiveGrowthRate' : logical, whether the adaptive growth 
%         rate method should be used. Default: false.
%
%       - 'MinimumGrowthRate' : double, minimum growth rate if using 
%         adaptive growth. Default: 0.
%
%       - 'MaximumGrowthRate' : double, maximum growth rate if using 
%         adaptive growth. Default: 0.2.
%
%       - 'MinimumGrowthPurity' : double, lower bound for the considered purity
%         if using an adaptive approach for trimming/expansion. Default: 0.001.
%
%       - 'MaximumGrowthPurity' : double, upper bound for the considered purity
%         if using an adaptive approach. Default: 0.999.
%
%       - 'GrowthAdaptationFactorFunction' : function handle used to adapt 
%         the growth rate based on measured volume/purity/etc. 
%         Default: @growth_rate_adaptation_volume.
%
%       - 'GrowthAdaptationFactorOptions' : cell array of additional options 
%         for the growth adaptation factor function. Default: {}.
%
%       - 'MinimumGrowthAdaptationFactor' : double, minimum factor the
%         growth rate can be multiplied by when adapted. Default: 0.2.
%
%       - 'MaximumGrowthAdaptationFactor' : double, maximum factor the
%         growth rate can be multiplied by when adapted. Default: 1.5.
%
%       - 'TargetAcceptedRatioExploration' : double, target ratio of accepted 
%         designs in exploration to drive adaptive growth. Default: 0.7.
%
%       - 'ToleranceMeasureChangeExploration' : double, how much the measure 
%         (box volume or candidate space volume) can change between iterations 
%         before being considered converged. Default: 1e-2.
%
%       - 'TolerancePurityConsolidation' : double, required purity level in 
%         consolidation for stopping. Default: 1.0.
%
%       - 'MaxIter' : integer, maximum total number of iterations for the entire
%         algorithm if the exploration and consolidation phase limits are 
%         not directly specified. Default: 40.
%
%       - 'MaxFunctionEvaluations' : integer, if specified, the maximum total 
%         number of function evaluations. This may be used to derive 
%         'MaxIterExploration' and 'MaxIterConsolidation' if those are empty.
%         Default: [].
%
%       - 'MaxIterExploration' : integer, maximum number of iterations 
%         allowed in the exploration phase. If empty, it may be derived from 
%         'MaxFunctionEvaluations' or 'MaxIter'. Default: [].
%
%       - 'MaxIterConsolidation' : integer, maximum number of iterations 
%         allowed in the consolidation phase. If empty, it may be derived from 
%         'MaxFunctionEvaluations' or 'MaxIter'. Default: [].
%
%       - 'FixIterNumberExploration' : logical, whether to fix the 
%         exploration phase to 'MaxIterExploration' or allow early stopping 
%         upon convergence. Default: true.
%
%       - 'FixIterNumberConsolidation' : logical, whether to fix the 
%         consolidation phase to 'MaxIterConsolidation' or allow early 
%         stopping upon convergence. Default: true.
%
%       - 'TrimmingOrderFunction' : function handle used to determine the 
%         order of trimming non-acceptable points or regions. Default: 
%         @trimming_order.
%
%       - 'TrimmingOrderOptions' : cell array of additional options for the 
%         trimming order function. Default: {}.
%
%       - 'TrimmingOperationFunction' : function handle used to perform 
%         trimming operations for the solution space (box or component). 
%         If 'box' SSO, defaults to @box_trimming_operation; if 'component' 
%         SSO, defaults to @component_trimming_operation.
%
%       - 'TrimmingOperationOptions' : cell array of additional options 
%         for the trimming operation. Default: {}.
%
%       - 'LoggingLevel' : controls the logging verbosity. 
%         Default: 'info' (see ConsoleLogging.SeverityLevelName).
%
%   Options exclusive to box SSO (SOLUTIONSPACETYPE == 'box'):
%       - 'MeasureFunction' : function handle to compute the "size" 
%         (measure) of the candidate box. Default: @box_measure_volume.
%
%       - 'MeasureOptions' : cell array of additional options for the measure 
%         function. Default: {}.
%
%   Options exclusive to component SSO (SOLUTIONSPACETYPE == 'component'):
%       - 'MinimumNumberPaddingSamples' : integer, minimum number of 
%         padding sample points to generate for each iteration. Default: 1000.
%
%       - 'MaximumNumberPaddingSamples' : integer or empty, maximum number 
%         of padding samples. If empty, no explicit maximum is enforced. 
%         Default: [].
%
%       - 'CandidateSpaceSamplingFunction' : function handle to sample 
%         candidate spaces. Default: 
%         @candidate_space_sampling_individual_feasible.
%
%       - 'CandidateSpaceSamplingOptions' : cell array of additional options 
%         for candidate space sampling. Default: {}.
%
%       - 'TrimmingMethodFunction' : function handle describing which 
%         "method" of trimming is performed. Default: 
%         @component_trimming_method_planar_trimming.
%
%       - 'TrimmingMethodOptions' : cell array of additional options for 
%         the trimming method function. Default: {}.
%
%       - 'TrimmingCostFunction' : function handle to compute the cost(s) 
%         of removing certain portions of the design space. Default: 
%         @component_trimming_cost.
%
%       - 'TrimmingCostOptions' : cell array of additional options for the 
%         trimming cost function. Default: {}.
%
%       - 'TrimmingComponentChoiceFunction' : function handle to choose 
%         among multiple possible trimming components based on cost or 
%         other criteria. Default: @component_trimming_component_choice.
%
%       - 'TrimmingComponentChoiceOptions' : cell array of additional 
%         options for the trimming component choice function. Default: {}.
%
%       - 'TrimmingOptimalChoiceFunction' : function handle to pick a 
%         possibly optimal trimming strategy after potential partial choices. 
%         Default: @component_trimming_optimal_choice.
%
%       - 'TrimmingOptimalChoiceOptions' : cell array of additional options 
%         for the trimming optimal choice function. Default: {}.
%
%       - 'UsePaddingSamplesInTrimming' : logical, whether the padding 
%         samples are also considered in trimming. Default: true.
%
%       - 'UseShapeSamplesExploration' : logical, whether shape 
%         (component-based) sample points should be used for exploration steps. 
%         Default: true.
%
%       - 'ShapeSamplesUsefulExploration' : logical, whether shape samples 
%         from the exploration phase are considered useful for subsequent 
%         steps. Default: false.
%
%       - 'ShapeSamplesUsefulConsolidation' : logical, whether shape 
%         samples are considered useful during consolidation. 
%         Default: true.
%
%       - 'UsePreviousEvaluatedSamplesConsolidation' : logical, whether 
%         previously evaluated samples from exploration are retained and 
%         reused during consolidation. Default: false.
%
%       - 'CandidateSpaceConstructor' : function handle for constructing 
%         the underlying candidate space representation. Default: 
%         @CandidateSpaceConvexHull.
%
%       - 'CandidateSpaceOptions' : cell array of options passed to the 
%         candidate space constructor. Default: {}.
%
%   Input:
%       - SOLUTIONSPACETYPE : char or string, 'box' or 'component'
%       - Various NAME,VALUE pairs as described above
%
%   Output:
%       - OPTIONS : struct of parsed options
%
%   See also sso_box_stochastic, design_requirement_spaces_label.
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

    parser = inputParser;

    %% common parameters
    parser.addParameter('RequirementSpacesType','Omega2',@(x)any(strcmpi(x,{'Omega0','Omega1','Omega2','Omega3'})));
    parser.addParameter('ApplyLeanness','always',@(x)any(strcmpi(x,{'always','end-only','never'})));
    % sampling
    parser.addParameter('NumberSamplesPerIteration',100,@(x)isnumeric(x)&&all(x>0));
    parser.addParameter('SamplingMethodFunction',@sampling_latin_hypercube,@(x)isa(x,'function_handle'));
    parser.addParameter('SamplingMethodOptions',{},@(x)iscell(x));
    parser.addParameter('NumberSamplesPerIterationExploration',[],@(x)isnumeric(x)&&all(x>0));
    parser.addParameter('NumberSamplesPerIterationConsolidation',[],@(x)isnumeric(x)&&all(x>0));
    % growth
    parser.addParameter('GrowthRate',0.1,@(x)isnumeric(x)&&isscalar(x)&&(x>0));
    parser.addParameter('UseAdaptiveGrowthRate',false,@(x)islogical(x)&&isscalar(x));
    parser.addParameter('MinimumGrowthRate',0,@(x)isnumeric(x)&&isscalar(x)&&(x>=0));
    parser.addParameter('MaximumGrowthRate',0.2,@(x)isnumeric(x)&&isscalar(x)&&(x>0));
    parser.addParameter('MinimumGrowthPurity',0.001,@(x)isnumeric(x)&&isscalar(x)&&(x>=0));
    parser.addParameter('MaximumGrowthPurity',0.999,@(x)isnumeric(x)&&isscalar(x)&&(x>0));
    parser.addParameter('GrowthAdaptationFactorFunction',@growth_rate_adaptation_volume,@(x)isa(x,'function_handle'));
    parser.addParameter('GrowthAdaptationFactorOptions',{},@(x)iscell(x));
    parser.addParameter('MinimumGrowthAdaptationFactor',0.2,@(x)isnumeric(x)&&isscalar(x)&&(x>=0));
    parser.addParameter('MaximumGrowthAdaptationFactor',1.5,@(x)isnumeric(x)&&isscalar(x)&&(x>0));
    parser.addParameter('TargetAcceptedRatioExploration',0.7,@(x)isnumeric(x)&&isscalar(x)&&(x>0));
    parser.addParameter('MinimumPurityReset',0.15,@(x)isnumeric(x)&&isscalar(x)&&(x>=0));
    % loop control
    parser.addParameter('ToleranceMeasureChangeExploration',1e-2,@(x)isnumeric(x)&&isscalar(x)&&(x>0));
    parser.addParameter('TolerancePurityConsolidation',1.0,@(x)isnumeric(x)&&isscalar(x)&&(x>=0));
    parser.addParameter('MaxIter',40,@(x)isnumeric(x)&&isscalar(x)&&(x>=0));
    parser.addParameter('MaxFunctionEvaluations',[],@(x)isnumeric(x)&&isscalar(x)&&(x>=0));
    parser.addParameter('MaxIterExploration',[],@(x)isnumeric(x)&&isscalar(x)&&(x>=0));
    parser.addParameter('MaxIterConsolidation',[],@(x)isnumeric(x)&&isscalar(x)&&(x>=0));
    parser.addParameter('FixIterNumberExploration',true,@(x)islogical(x)&&isscalar(x));
    parser.addParameter('FixIterNumberConsolidation',true,@(x)islogical(x)&&isscalar(x));
    % trimming
    parser.addParameter('TrimmingOrderFunction',@trimming_order,@(x)isa(x,'function_handle'));
    parser.addParameter('TrimmingOrderOptions',{},@(x)iscell(x));
    parser.addParameter('TrimmingOperationOptions',{},@(x)iscell(x));
    % logging verbosity
    parser.addParameter('LoggingLevel','info',@(x)any(strcmpi(x,ConsoleLogging.SeverityLevelName)));


    %% box parameters
    if(strcmpi(solutionSpaceType,'box'))
        % trimming / box measure
        parser.addParameter('TrimmingOperationFunction',@box_trimming_operation,@(x)isa(x,'function_handle'));
        parser.addParameter('MeasureFunction',@box_measure_volume,@(x)isa(x,'function_handle'));
        parser.addParameter('MeasureOptions',{},@(x)iscell(x));
    end


    %% component parameters
    if(strcmpi(solutionSpaceType,'component'))
        % sampling
        parser.addParameter('MinimumNumberPaddingSamples',1000,@(x)isnumeric(x)&&all(x>=0));
        parser.addParameter('MaximumNumberPaddingSamples',[],@(x)isnumeric(x)&&all(x>=0));
        parser.addParameter('CandidateSpaceSamplingFunction',@candidate_space_sampling_individual_feasible,@(x)isa(x,'function_handle'));
        parser.addParameter('CandidateSpaceSamplingOptions',{},@(x)iscell(x));
        % trimming
        parser.addParameter('TrimmingOperationFunction',@component_trimming_operation,@(x)isa(x,'function_handle'));
        parser.addParameter('TrimmingMethodFunction',@component_trimming_method_planar_trimming,@(x)isa(x,'function_handle'));
        parser.addParameter('TrimmingMethodOptions',{},@(x)iscell(x));
        parser.addParameter('TrimmingCostFunction',@component_trimming_cost,@(x)isa(x,'function_handle'));
        parser.addParameter('TrimmingCostOptions',{},@(x)iscell(x));
        parser.addParameter('TrimmingComponentChoiceFunction',@component_trimming_component_choice,@(x)isa(x,'function_handle'));
        parser.addParameter('TrimmingComponentChoiceOptions',{},@(x)iscell(x));
        parser.addParameter('TrimmingOptimalChoiceFunction',@component_trimming_optimal_choice,@(x)isa(x,'function_handle'));
        parser.addParameter('TrimmingOptimalChoiceOptions',{},@(x)iscell(x));
        parser.addParameter('UsePaddingSamplesInTrimming',true,@(x)islogical(x)&&isscalar(x));
        parser.addParameter('UseShapeSamplesExploration',true,@(x)islogical(x)&&isscalar(x));
        parser.addParameter('ShapeSamplesUsefulExploration',false,@(x)islogical(x)&&isscalar(x));
        parser.addParameter('ShapeSamplesUsefulConsolidation',true,@(x)islogical(x)&&isscalar(x));
        parser.addParameter('UsePreviousEvaluatedSamplesConsolidation',false,@(x)islogical(x)&&isscalar(x));
        % candidate spaces
        parser.addParameter('CandidateSpaceConstructor',@CandidateSpaceConvexHull,@(x)isa(x,'function_handle'));
        parser.addParameter('CandidateSpaceOptions',{},@(x)iscell(x));
    end    

    parser.parse(varargin{:});
    options = parser.Results;
    
    %% Options Processing
    % number of samples per iteration
    if(isempty(options.NumberSamplesPerIterationExploration))
        options.NumberSamplesPerIterationExploration = options.NumberSamplesPerIteration;
    end
    if(isempty(options.NumberSamplesPerIterationConsolidation))
        options.NumberSamplesPerIterationConsolidation = options.NumberSamplesPerIteration;
    end

    % maximum number of iterations per phase by number of function evaluations
    if(~isempty(options.MaxFunctionEvaluations))
        if(isempty(options.MaxIterExploration) && isempty(options.MaxIterConsolidation))
            % distribute function evaluations equally
            maxFunctionEvaluationsExploration = round(options.MaxFunctionEvaluations/2);
            maxFunctionEvaluationsConsolitation = round(options.MaxFunctionEvaluations/2);
        elseif(isempty(options.MaxIterConsolidation))
            maxFunctionEvaluationsConsolitation = options.MaxFunctionEvaluations - ...
                options.MaxIterExploration*options.NumberSamplesPerIterationExploration;
        elseif(isempty(options.MaxIterExploration))
            maxFunctionEvaluationsExploration = options.MaxFunctionEvaluations - ...
                options.MaxIterConsolidation*options.NumberSamplesPerIterationConsolidation;
        end

        options.MaxIterExploration = round(maxFunctionEvaluationsExploration / ...
            options.NumberSamplesPerIterationExploration);
        options.MaxIterConsolidation = round(maxFunctionEvaluationsConsolitation / ...
            options.NumberSamplesPerIterationConsolidation);
    end

    % maximum number of iterations per phase by total number of iterations
    if(~isempty(options.MaxIter))
        if(isempty(options.MaxIterExploration) && isempty(options.MaxIterConsolidation))
            % assume equal distribution of iterations
            options.MaxIterExploration = options.MaxIter/2;
            options.MaxIterConsolidation = options.MaxIter - options.MaxIterExploration;
        elseif(isempty(options.MaxIterConsolidation))
            options.MaxIterConsolidation = options.MaxIter - options.MaxIterExploration;
        elseif(isempty(options.MaxIterExploration))
            options.MaxIterExploration = options.MaxIter - options.MaxIterConsolidation;
        end
    end
end

