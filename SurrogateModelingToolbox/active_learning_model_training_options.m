function options = active_learning_model_training_options(varargin)
%ACTIVE_LEARNING_MODEL_TRAINING_OPTIONS Options for the active learning training
%   ACTIVE_LEARNING_MODEL_TRAINING_OPTIONS returns a structure that contains 
%   all options required by 'active_learning_model_training'.
%
%   OPTIONS = ACTIVE_LEARNING_MODEL_TRAINING_OPTIONS(...NAME,VALUE,...) creates 
%   such a structure OPTIONS where . All others are set to their default.
%       - FastForwardModelInitial : initial fast-forward model, which can be 
%       given as a DesignFastForward object, or a constructor. Default: 
%       @DesignFastForwardAnn.
%
%       - FastForwardModelOptions : options for the design fast-forward model if
%       it is given as a constructor. Default is empty.
%
%       - NewModelInitialRegion : if the model was given with a constructor,
%       an initial region is required to start sampling/training there. Default 
%       is the middle of the design space.
%
%       - NewModelInitialRegionGrowth : this can set a fixed factor for the
%       growth of the region above before an initial sampling. Default: 0.1.
%
%       - ActiveLearningSamplingFunction : function used to generate candidate  
%       samples for probing and testing. Default: 
%       @active_learning_sampling_random.
%
%       - ActiveLearningSamplingOptions : options for the sampling function. 
%       Default is empty.
%
%       - MeasureChoice : choice of measure to compute the surrogate model for,
%       either 'performance' or 'physical-feasibility'. Default is 
%       'performance'.
%
%       - NumberSamplesEvaluation : number of new samples to be evaluated for 
%       each iteration of training. Default value is '100'.
%
%       - NumberCandidatesInitial : the initial number of candidates that can
%       be chosen to be evaluated (becoming evaluated design samples). 
%       
%       - FunctionUpdateNumberCandidates : a function to update the number of
%       candidates for each iteration.
%
%       - MaxIter : maximum number of iterations for the training algorithm.
%       Default value is '20'.
%
%       - FixIter : flag to flag to determine if the number of iterations 
%       performed be fixed to its maximum value ('true'), or if convergence 
%       criteria should be checked to stop the algorithm ('false'). Default 
%       value is 'true'.
%
%       - ExplorationErrorTolerance : the tolerance of error rate for 
%       exploratory samples; if the number of iterations isn't fixed,
%       then one criterion to determine if the model is good enough is that the
%       error rate is smaller than this tolerance after a given number of times.
%       Default value is '0.01' (1% error rate).
%
%       - CyclesToConverge : number of times it takes for the tolerance to be 
%       met for it to be considered the method has converged and the model is
%       good enough. Default value is '3'.
%
%       - ExplorationToBoundaryRatio : ratio of design samples chosen to as
%       exploratory, meaning they are not necessarily close to the boundary
%       between positive and negative samples. Set to '1' for only exploratory
%       designs to be used during training, and '0' for only boundary designs.
%       Default value is '0.5'.
%
%       - LoggingLevel : setting for what should be logged to console. 
%       Default value is 'info'. See more on 'ConsoleLogging'.
%
%   
%   Inputs:
%       - 'FastForwardModelInitial' : DesignFastForwardBase OR function_handle
%       - 'FastForwardModelOptions' : cell
%       - 'NewModelInitialRegion' : (2,nDesignVariable) double
%       - 'NewModelInitialRegionGrowth' : double
%       - 'SamplingMethod' : function handle
%       - 'SamplingOptions' : cell
%       - 'VectorNormType' : double
%       - 'DistanceOptions' : cell
%       - 'LineSearchMaxIter' : double
%       - 'MeasureChoice' : char OR string
%       - 'NumberSamplesEvaluation' : double
%       - 'NumberCandidatesInitial' : double
%       - 'FunctionUpdateNumberCandidates' : function_handle
%       - 'MaxIter' : double
%       - 'FixIter' : logical
%       - 'ExplorationErrorTolerance' : double
%       - 'CyclesToConverge' : double
%       - 'ExplorationToBoundaryRatio' : double
%       - 'LoggingLevel' : char OR string
%
%   Outputs:
%       - OPTIONS : struct
%
%   See also active_learning_model_training.
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
    parser.addParameter('FastForwardModelInitial',@DesignFastForwardAnn,@(x)isa(x,'function_handle')||isa(x,'DesignFastForward'));
    parser.addParameter('FastForwardModelOptions',{},@(x)iscell(x));
    parser.addParameter('ActiveLearningSamplingFunction',@active_learning_sampling_random,@(x)isa(x,'function_handle'));
    parser.addParameter('ActiveLearningSamplingOptions',{},@(x)iscell(x));
    parser.addParameter('NewModelInitialRegion',[],@(x)isnumeric(x));
    parser.addParameter('NewModelInitialRegionGrowth',0.1,@(x)isnumeric(x)&&isscalar(x));
    parser.addParameter('MeasureChoice','solution',@(x)ischar(x)||isstring(x));
    parser.addParameter('NumberSamplesEvaluation',100,@(x)isnumeric(x));
    parser.addParameter('NumberCandidatesInitial',1000,@(x)isnumeric(x)&&isscalar(x));
    parser.addParameter('FunctionUpdateNumberCandidates',@(nCandidate0, nCandidate, iter)[nCandidate0*iter],@(x)isa(x,'function_handle'));
    parser.addParameter('MaxIter',20,@(x)isnumeric(x)&&isscalar(x));
    parser.addParameter('FixIter',false,@(x)islogical(x)&&isscalar(x));
    parser.addParameter('ExplorationErrorTolerance',0.01,@(x)isnumeric(x)&&isscalar(x));
    parser.addParameter('CyclesToConverge',3,@(x)isnumeric(x)&&isscalar(x));
    parser.addParameter('ExplorationToBoundaryRatio',0.50,@(x)isnumeric(x)&&isscalar(x));
    parser.addParameter('LoggingLevel','info',@(x)any(strcmpi(x,ConsoleLogging.SeverityLevelName)));
    parser.parse(varargin{:});
    options = parser.Results;
end

