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
%   those options. Available for both box SSO and component SSO are:
%       - 'RequirementSpacesType' : type of requirement spaces to be used during
%       trimming; avialable options are 'Omega0','Omega1','Omega2' and 'Omega3'.
%       See more about this on 'design_requirement_spaces_label'. Default value
%       is 'Omega2'.
%
%       - 'ApplyLeanness' : determine when/if the leanness condition of 
%       requirement spaces should be applied. 'always' means it is applied at 
%       the end of each trimming step, 'end-only' only applies it after the 
%       last consolidation step, and 'never' the condition is not applied.
%
%       - 'NumberSamplesPerIteration' : number of samples to be evaluated in each
%       iteration of the algorithm, in both exploration and consolidation 
%       phases. Default value is 100.
%
%       - 'SamplingMethodFunction' : function used to stochastically sample 
%       inside a given design box, or base sampling method for sampling the 
%       candidate spaces in component SSO. Default: @sampling_latin_hypercube.
%
%       - 'SamplingMethodOptions' : options for the (base) sampling function. 
%       Default is empty.
%
%       - 'NumberSamplesPerIterationExploration' : number of samples to be used
%       during the exploration phase; set to 'NumberSamplesPerIteration' if 
%       empty.
%
%       - 'NumberSamplesPerIterationConsolidation' : number of samples to be 
%       used during the consolidation phase; set to 'NumberSamplesPerIteration'  
%       if empty.
%
%       - 'GrowthRate' : fixed growth rate of the box / candidate spaces during 
%       exploration. Default: 0.2.
%
%       - 'UseAdaptiveGrowthRate' : flag to determine whether the adaptive 
%       growth rate method should be used to try to keep a given ratio of 
%       accepted / not accepted designs. Default: false.
%
%       - 'MinimumGrowthRate' : minimum value the growth rate can assume if the
%       adaptive growth method is used. Default value is '0.01'.
%
%       - 'TargetAcceptedRatioExploration' : target value for the adaptive 
%       growth method. Default value is '0.8'.
%
%       - 'ToleranceMeasureChangeExploration' : tolerance for how much the 
%       overall measure of the box / candidate space can vary between iterations
%       before it is considered the algorithm has converged. Default: 0.01.
%
%       - 'TolerancePurityConsolidation' : tolerance for the level of purity 
%       required for the consolidation to be considered done. Default: 1.0.
%
%       - 'MaxIter' : maximum number of iterations for the whole algorithm. 
%       Default: 40.
%
%       - 'MaxFunctionEvaluations' : maximum number of function evaluations. If
%       specified, the maximum number of iterations for the exploration and 
%       consolidation phases will be derived from this (if they have not been
%       specified themselves). Default is empty.
%
%       - 'MaxIterExploration' : maximum number of iterations allowed during the
%       exploration phase. If empty, it first sees if that value can be
%       obtained from the maximum number of function evaluations; if not, 
%       it uses maximum number of iterations. Default is empty.
%
%       - 'MaxIterConsolidation' : maximum number of iterations allowed during the
%       consolidation phase. If empty, it first sees if that value can be
%       obtained from the maximum number of function evaluations; if not, 
%       it uses maximum number of iterations. Default is empty.
%
%       - 'FixIterNumberExploration' : flag to determine if the number of 
%       iterations performed during exploration should be fixed to its maximum
%       value ('true'), or if convergence criteria should be checked to stop
%       the algorithm ('false'). Default: 'true'.
%
%       - 'FixIterNumberConsolidation' : flag to determine if the number of 
%       iterations performed during consolidation should be fixed to its maximum
%       value ('true'), or if convergence criteria should be checked to stop
%       the algorithm ('false'). Default: 'true'.
%
%       - 'TrimmingOrderFunction' : function to be used when determining
%       the order that non-acceptable designs should be trimmed out of the
%       candidate space. Default: @trimming_order.
%
%       - 'TrimmingOrderOptions' : options to be used for for the trimming
%       order function. Default: empty.
%
%       - 'TrimmingOperationFunction' : function to be used when performing the 
%       trimming operation for each component. Default is 
%       'box_trimming_operation' for box SSO and 'component_trimming_operation'
%       for component SSO.
%
%       - 'TrimmingOperationOptions' : options to be used for for the trimming
%       operation function. Default is empty. 
%
%       - 'LoggingLevel' : setting for what should be logged to console. See 
%       more on 'ConsoleLogging'. Default: 'info'. 
%   
%   Available exclusively for box SSO are:
%       - 'MeasureFunction' : function to get the measure of a candidate box.
%       Default: @box_measure_volume.
%
%       - 'MeasureOptions' : options for the measure function. Default: empty.
%
%   Available exclusively for component SSO are:
%       - 'NumberPaddingSamples' : number of padding sample points to generate.
%       These are points that used during trimming and candidate space 
%       definition, but don't get evaluated. Default: 1000.
%
%       - 'CandidateSpaceSamplingFunction' : function to sample the candidate 
%       spaces. Default: @candidate_space_sampling_individual_feasible.
%       
%       - 'CandidateSpaceSamplingOptions' : options for the sampling of the 
%       candidate space. 'SamplingMethodFunction' and 'SamplingMethodOptions'
%       get added to this when the algorithm is actually running. 
%       Default is empty.
%
%       - 'TrimmingMethodFunction' : method of trimming used, as in how the
%       unnaceptable designs are removed. Default: 
%       @component_trimming_method_planar_trimming.
%
%       - 'TrimmingMethodOptions' : options for the trimming method. Default is
%       empty.
%
%       - 'TrimmingCostFunction' : function to compute the costs of the removal 
%       of different trimming candidates. Default: @component_trimming_cost.
%
%       - 'TrimmingCostOptions' : options for the trimming cost function.
%       Default is empty.
%
%       - 'TrimmingComponentChoiceFunction' : function to choose from the 
%       available component trimmings which one should be done based on the 
%       cost. Default: @component_trimming_choice.
%
%       - 'TrimmingComponentChoiceOptions' : options for the trimming choice 
%       function. Default is empty.
%
%       - 'UsePaddingSamplesInTrimming' : flag to decide if padding samples 
%       should be considered in the trimming operation (and, as a consequence,
%       on candidate space definition/update). Default: true.
%
%       - 'UsePreviousEvaluatedSamplesConsolidation' : flag to decide whether
%       all previously evaluated samples should be kept and remain in the 
%       trimming process / candidate space update during consolidation. Default:
%       false.
%
%       - 'UsePreviousPaddingSamplesConsolidation' : flag to decide whether
%       all previous padding samples should be kept and remain in the trimming
%       process / candidate space update during consolidation. Default: false.
%
%       - 'CandidateSpaceConstructorExploration' : constructor function handle
%       for the candidate spaces to be used during exploration. Default:
%       @CandidateSpaceConvexHull.
%
%       - 'CandidateSpaceConstructorConsolidation' : constructor function handle
%       for the candidate spaces to be used during consolidation. Default:
%       @CandidateSpaceConvexHull.
%
%       - 'CandidateSpaceOptionsExploration' : options to be used when first
%       initializing the candidate spaces with the constructors during 
%       exploration. Default is empty.
%
%       - 'CandidateSpaceOptionsConsolidation' : options to be used when first
%       initializing the candidate spaces with the constructors during 
%       exploration. Default is empty.
%
%   
%   Input:
%       - SOLUTIONSPACETYPE : char OR string
%       - 'RequirementSpacesType' : char OR string OR function_handle
%       - 'ApplyLeanness' : logical
%       - 'NumberSamplesPerIteration' : integer 
%       - 'SamplingMethodFunction' : function_handle
%       - 'SamplingMethodOptions' : (1,nOptionSamplingMethod) cell
%       - 'NumberSamplesPerIterationExploration' : integer
%       - 'NumberSamplesPerIterationConsolidation' : integer
%       - 'GrowthRate' : double
%       - 'AdaptiveGrowthRate' : logical
%       - 'MinimumGrowthRate' : double
%       - 'TargetAcceptedRatioExploration' : double
%       - 'ToleranceMeasureChangeExploration' : double
%       - 'TolerancePurityConsolidation' : double
%       - 'MaxIter' : integer
%       - 'MaxFunctionEvaluations' : integer
%       - 'MaxIterExploration' : integer
%       - 'MaxIterConsolidation' : integer
%       - 'FixIterNumberExploration' : logical
%       - 'FixIterNumberConsolidation' : logical
%       - 'TrimmingOrderFunction' : function_handle
%       - 'TrimmingOrderOptions' : (1,nOptionOrder) cell
%       - 'TrimmingOperationFunction' : function_handle
%       - 'TrimmingOperationOptions' : (1,nOptionOperation) cell
%       - 'LoggingLevel' : char OR string OR integer
%       - 'MeasureFunction' : function_handle
%       - 'MeasureOptions' : (1,nOptionMeasure) cell
%       - 'NumberPaddingSamples' : integer
%       - 'CandidateSpaceSamplingFunction' : function_handle
%       - 'CandidateSpaceSamplingOptions' : (1,nOptionCandidateSampling) cell
%       - 'TrimmingMethodFunction' : function_handle
%       - 'TrimmingMethodOptions' : (1,nOptionTrimmingMethod) cell
%       - 'TrimmingCostFunction' : function_handle
%       - 'TrimmingCostOptions' : (1,nOptionCost) cell
%       - 'TrimmingComponentChoiceFunction' : function_handle
%       - 'TrimmingComponentChoiceOptions' : (1,nOptionChoice) cell
%       - 'UsePaddingSamplesInTrimming' : logical
%       - 'UsePreviousEvaluatedSamplesConsolidation' : logical
%       - 'UsePreviousPaddingSamplesConsolidation' : logical
%       - 'CandidateSpaceConstructorExploration' : funciton_handle
%       - 'CandidateSpaceConstructorConsolidation' : funciton_handle
%       - 'CandidateSpaceOptionsExploration' : (1,nOptionCandidate) cell
%       - 'CandidateSpaceOptionsConsolidation' : (1,nOptionCandidate) cell
%
%   Output:
%       - OPTIONS : struct
%
%   See also sso_box_stochastic, design_requirement_spaces_label.
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
    parser.addParameter('GrowthRate',0.2,@(x)isnumeric(x)&&isscalar(x)&&(x>0));
    parser.addParameter('UseAdaptiveGrowthRate',false,@(x)islogical(x)&&isscalar(x));
    parser.addParameter('MinimumGrowthRate',0.01,@(x)isnumeric(x)&&isscalar(x)&&(x>0));
    parser.addParameter('TargetAcceptedRatioExploration',0.8,@(x)isnumeric(x)&&isscalar(x)&&(x>0));
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
        parser.addParameter('NumberPaddingSamples',1000,@(x)isnumeric(x)&&all(x>0));
        parser.addParameter('CandidateSpaceSamplingFunction',@candidate_space_sampling_individual_feasible,@(x)isa(x,'function_handle'));
        parser.addParameter('CandidateSpaceSamplingOptions',{},@(x)iscell(x));
        % trimming
        parser.addParameter('TrimmingOperationFunction',@component_trimming_operation,@(x)isa(x,'function_handle'));
        parser.addParameter('TrimmingMethodFunction',@component_trimming_method_planar_trimming,@(x)isa(x,'function_handle'));
        parser.addParameter('TrimmingMethodOptions',{},@(x)iscell(x));
        parser.addParameter('TrimmingCostFunction',@component_trimming_cost,@(x)isa(x,'function_handle'));
        parser.addParameter('TrimmingCostOptions',{},@(x)iscell(x));
        parser.addParameter('TrimmingComponentChoiceFunction',@component_trimming_choice,@(x)isa(x,'function_handle'));
        parser.addParameter('TrimmingComponentChoiceOptions',{},@(x)iscell(x));
        parser.addParameter('UsePaddingSamplesInTrimming',true,@(x)islogical(x)&&isscalar(x));
        parser.addParameter('UsePreviousEvaluatedSamplesConsolidation',false,@(x)islogical(x)&&isscalar(x));
        parser.addParameter('UsePreviousPaddingSamplesConsolidation',false,@(x)islogical(x)&&isscalar(x));
        % candidate spaces
        parser.addParameter('CandidateSpaceConstructorExploration',@CandidateSpaceConvexHull,@(x)isa(x,'function_handle'));
        parser.addParameter('CandidateSpaceConstructorConsolidation',@CandidateSpaceConvexHull,@(x)isa(x,'function_handle'));
        parser.addParameter('CandidateSpaceOptionsExploration',{},@(x)iscell(x));
        parser.addParameter('CandidateSpaceOptionsConsolidation',{},@(x)iscell(x));
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

