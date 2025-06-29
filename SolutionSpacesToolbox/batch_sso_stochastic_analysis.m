function [solutionSpace,problemData,iterData,algoData,batchOptions] = batch_sso_stochastic_analysis(batchOptions,designEvaluator,initialPoint,designSpaceLowerBound,designSpaceUpperBound,varargin)
%BATCH_SSO_STOCHASTIC_ANALYSIS Solving same SSO problem with different options
%   BATCH_SSO_STOCHASTIC_ANALYSIS can be used to procedurally find optimal
%   solution spaces for a given problem with different options. The results
%   can then be used to better understand either how different algorithm 
%   parameters affect the convergence / final solution.
%   
%   SOLUTIONSPACE = BATCH_SSO_STOCHASTIC_ANALYSIS(BATCHOPTIONS,DESIGNEVALUATOR,
%   INITIALPOINT,DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND) takes the options
%   specified in BATCHOPTIONS (as either a structure or a table file) and solves
%   the box SSO problem with those options, using DESIGNEVALUATOR and design 
%   space boundaries DESIGNSPACELOWERBOUND and DESIGNSPACEUPPERBOUND, starting 
%   every test at INITIALPOINT. It returns the resulting solution boxes in 
%   SOLUTIONSPACE.
%
%   SOLUTIONSPACE = BATCH_SSO_STOCHASTIC_ANALYSIS(BATCHOPTIONS,DESIGNEVALUATOR,
%   INITIALPOINT,DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND,COMPONENTINDEX) 
%   instead solves the component SSO problem, with components defined in 
%   COMPONENTINDEX. In this case, SOLUTIONSPACE is a cell with all the component
%   solution spaces.
%
%   SOLUTIONSPACE = BATCH_SSO_STOCHASTIC_ANALYSIS(...,NAME,VALUE,...) also
%   allows the specification of the additional name-value pair arguments.
%   These are:
%       - 'FixedOptions' : options that are the same through all tests. Assumed
%       to be empty if not specified.
%
%   [SOLUTIONSPACE,PROBLEMDATA] = BATCH_SSO_STOCHASTIC_ANALYSIS(...) 
%   additionally returns the PROBLEMDATA structure for each test, containing
%   all initial information.
%
%   [SOLUTIONSPACE,PROBLEMDATA,ITERDATA] = BATCH_SSO_STOCHASTIC_ANALYSIS(...) 
%   additionally returns the ITERDATA structure for each test, containing
%   all information generated at each iteration of the stochastic procedure.
%
%   [SOLUTIONSPACE,PROBLEMDATA,ITERDATA,ALGODATA] = 
%   BATCH_SSO_STOCHASTIC_ANALYSIS(...) additionally returns the ALGODATA 
%   structure for each test, containing the algorithm performance data, mainly
%   used to generate metrics plots.
%
%   [SOLUTIONSPACE,PROBLEMDATA,ITERDATA,ALGODATA,BATCHOPTIONS] = 
%   BATCH_SSO_STOCHASTIC_ANALYSIS(...) finally also returns the (processed) 
%   options used for each test of the batch analysis in BATCHOPTIONS;
%   this output is the same as the input with the same name if the options were
%   already processed then, but will differ if a table file was given instead.
%
%   Input:
%       - BATCHOPTIONS : cell OR string OR (nTest,1) structure
%       - DESIGNEVALUATOR : DesignEvaluatorBase
%       - INITIALPOINT : (1,nDesignVariable) double
%       - DESIGNSPACELOWERBOUND : (1,nDesignVariable) double
%       - DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double
%       - COMPONENTINDEX : (1,nComponent) cell
%       - 'FixedOptions' : (1,nExtraOption) cell
%
%   Output:
%       - SOLUTIONSPACE : (nTest,1) cell
%       - PROBLEMDATA : (nTest,1) cell
%       - ITERDATA : (nTest,1) cell
%       - ALGODATA : (nTest,1) structure
%       - BATCHOPTIONS : (nTest,1) structure
%
%   See also sso_box_stochastic, sso_box_stochastic_options, 
%   sso_component_stochastic, sso_component_stochastic_options,
%   batch_analysis_read_table, 
%   plot_batch_sso_stochastic_analysis_component_metrics.
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
    % requirement spaces
    parser.addOptional('ComponentIndex',{});
    parser.addParameter('FixedOptions',{});
    parser.parse(varargin{:});
    options = parser.Results;

    % process batch information
    if(ischar(batchOptions) || isstring(batchOptions))
        batchOptions = batch_analysis_read_table(batchOptions);
    end

    % save rng state
    rngState = rng;
    
    %% run batch tests for component SSO
    nTest = length(batchOptions);
    solutionSpace = cell(nTest,1);
    problemData = cell(nTest,1);
    iterData = cell(nTest,1);
    for i=1:nTest
        % create folder to save data    
        testOptions = namedargs2cell(batchOptions(i).Options);
        rng(rngState);

        if(isempty(options.ComponentIndex))
            ssoOptions = sso_stochastic_options('box',options.FixedOptions{:},testOptions{:});
            [solutionSpace{i},problemData{i},iterData{i}] = sso_box_stochastic(...
                designEvaluator,...
                initialPoint,...
                designSpaceLowerBound,...
                designSpaceUpperBound,...
                ssoOptions);
            algoData(i) = postprocess_sso_box_stochastic(problemData{i},iterData{i});
        else
            ssoOptions = sso_stochastic_options('component',options.FixedOptions{:},testOptions{:});
            [solutionSpace{i},problemData{i},iterData{i}] = sso_component_stochastic(...
                designEvaluator,...
                initialPoint,...
                designSpaceLowerBound,...
                designSpaceUpperBound,...
                options.ComponentIndex,...
                ssoOptions);
            algoData(i) = postprocess_sso_component_stochastic(problemData{i},iterData{i});
        end
    end
end