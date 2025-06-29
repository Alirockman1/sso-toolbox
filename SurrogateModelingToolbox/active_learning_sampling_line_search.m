function [designSample,scorePrediction,labelBoundary] = active_learning_sampling_line_search(fastForwardModel,designSpace,nDesiredSample,nCandidate,varargin)
%ACTIVE_LEARNING_SAMPLING_LINE_SEARCH Sampling for training
%   ACTIVE_LEARNING_SAMPLING_LINE_SEARCH uses the positive samples
%   in the current design fast-forward model to create new samples close to
%   the boundary and across the positive/negative regions. 
%   In the first step, line-search is performed to look for the 
%   positive/negative boundary  starting from the positive samples and moving 
%   in randomly-determined  directions. Once the step sizes to hit the boundary 
%   have been found, those step-sizes are used to determine the boundary 
%   candidates, and the positive/negative region candidates are determined by 
%   randomly choosing steps between the limits. 
%   From all candidates, the ones with the biggest distance to the already
%   evaluated samples used in the training of the fast-forward model are 
%   selected.
%
%   DESIGNSAMPLE = ACTIVE_LEARNING_SAMPLING_LINE_SEARCH(
%   FASTFORWARDMODEL,DESIGNSPACE,NDESIREDSAMPLE,NCANDIDATE) uses the current
%   fast-forward model FASTFORWARDMODEL and produces design samples DESIGNSAMPLE
%   with distribution/properties as requested in NDESIREDSAMPLE; in the process,
%   a number of candidates NCANDIDATE are used. NDESIREDSAMPLE is composed of 
%   four fields specifying how many samples are needed in regards to region/
%   proximity to the boundary: 'PositiveBoundary', 'NegativeBoundary',
%   'PositiveExploratory' and 'NegativeExploratory'.
%
%   DESIGNSAMPLE = ACTIVE_LEARNING_SAMPLING_LINE_SEARCH(
%   FASTFORWARDMODEL,DESIGNSPACE,NDESIREDSAMPLE,NCANDIDATE,OPTIONS) also allows
%   for the specification of options for the method in OPTIONS. This can be 
%   given as a structure or as name-value pair arguments. Said options are:
%       - VectorNormType : the type of vector norm being used ('2''for L2, etc.)
%       - MaxIterSearch : how many search iterations are done during bracketing/
%       bissectioning when performing line-search.
%       - KnnSearchOptions : options for the 'knnsearch' call when looking for
%       designs with the maximum distance to the ones already evaluated.
%
%   [DESIGNSAMPLE,SCOREPREDICTION] = 
%   ACTIVE_LEARNING_SAMPLING_LINE_SEARCH(...) also returns the 
%   predicted score SCOREPREDICTION for each sample in DESIGNSAMPLE. This may
%   allow for further analysis in terms of the error rate of the current
%   fast-forward model.
%
%   [DESIGNSAMPLE,SCOREPREDICTION,LABELBOUNDARY] = 
%   ACTIVE_LEARNING_SAMPLING_LINE_SEARCH(...) also returns the labels
%   of each sample regarding if they were selected due to being close to the
%   positive/negative boundary ('true') or not ('false'). Those marked as 'true'
%   are therefore boundary samples, and the ones marked as 'false' are 
%   exploratory samples.
%
%   Inputs:
%       - FASTFORWARDMODEL : DesignFastForward
%       - DESIGNSPACE : (2,nDesignVariable) double
%       - NDESIREDSAMPLE : struct
%           -- PositiveBoundary : integer
%           -- NegativeBoundary : integer
%           -- PositiveExploratory : integer
%           -- NegativeExploratory : integer
%       - NCANDIDATE : integer
%       - OPTIONS : struct OR name-value pair arguments
%           -- VectorNormType : double
%           -- MaxIterSearch : double
%           -- KnnSearchOptions : cell
%
%   Outputs:
%       - DESIGNSAMPLE : (nSample,nDesignVariable) double
%       - SCOREPREDICTION : (nSample,1) double
%       - LABELBOUNDARY : (nSample,1) logical
%
%   See also active_learning_model_training, design_select_max_distance.
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
    parser.addParameter('VectorNormType',2,@(x)isnumeric(x));
    parser.addParameter('MaxLineSearchIterations',[],@(x)isnumeric(x));
    parser.addParameter('knnsearchOptions',{},@(x)iscell(x));
    parser.addParameter('InitialPointSelectionFunction',@randsample,@(x)isa(x,'function_handle'));
    parser.addParameter('InitialPointSelectionOptions',{true},@(x)iscell(x));
    parser.addParameter('DirectionSamplingFunction',@sampling_random,@(x)isa(x,'function_handle'));
    parser.addParameter('DirectionSamplingOptions',{},@(x)iscell(x));
    parser.addParameter('BoundarySamplingFunction',@sample_line_search_boundary,@(x)isa(x,'function_handle'));
    parser.addParameter('BoundarySamplingOptions',{},@(x)iscell(x));
    parser.addParameter('ExploratorySamplingFunction',@sample_line_search_exploratory,@(x)isa(x,'function_handle'));
    parser.addParameter('ExploratorySamplingOptions',{},@(x)iscell(x));
    parser.parse(varargin{:});
    options = parser.Results;

	% get the designs used in training
	positiveTraining = fastForwardModel.DesignSampleTrain(fastForwardModel.LabelTrain,:);
    negativeTraining = fastForwardModel.DesignSampleTrain(~fastForwardModel.LabelTrain,:);

    % start from the positive samples
    nInitialSample = size(positiveTraining,1);
    iInitial = options.InitialPointSelectionFunction(nInitialSample,nCandidate,options.InitialPointSelectionOptions{:});
    initialPoint = positiveTraining(iInitial,:);
    nDimension = size(designSpace,2);

    % directions where to perform linesearch
    direction = options.DirectionSamplingFunction([-ones(1,nDimension);ones(1,nDimension)],nCandidate,options.DirectionSamplingOptions{:});
    normalizedDirection = direction./vecnorm(direction,options.VectorNormType,2);

    % look for designs inside the positive region
	regionCriterium = @(x)[design_deficit_to_label_score(fastForwardModel.predict_designs(x))];
    [stepSizePositiveMax,stepSizeNegativeMin] = region_limit_line_search(regionCriterium,initialPoint,normalizedDirection,designSpace,options.MaxLineSearchIterations);
    
    % look for designs inside the negative region
    initialPointNegative = initialPoint + stepSizeNegativeMin.*normalizedDirection;
    [stepSizeNegativeMax,~] = region_limit_line_search(regionCriterium,initialPointNegative,normalizedDirection,designSpace,options.MaxLineSearchIterations);

    % processing possible positive boundary candidate
    candidateBoundaryPositive = options.BoundarySamplingFunction(initialPoint,stepSizePositiveMax,normalizedDirection,options.BoundarySamplingOptions{:});
    [candidateBoundaryPositive,scoreBoundaryPositive] = select_correct_label_candidates(fastForwardModel,candidateBoundaryPositive,true);

    % processing possible negative boundary candidate
    candidateBoundaryNegative = options.BoundarySamplingFunction(initialPoint,stepSizeNegativeMin,normalizedDirection,options.BoundarySamplingOptions{:});
    [candidateBoundaryNegative,scoreBoundaryNegative] = select_correct_label_candidates(fastForwardModel,candidateBoundaryNegative,false);

    % processing possible positive exploratory candidates
    candidateExploratoryPositive = options.ExploratorySamplingFunction(initialPoint,0,stepSizePositiveMax,normalizedDirection,options.ExploratorySamplingOptions{:});
    [candidateExploratoryPositive,scoreExploratoryPositive] = select_correct_label_candidates(fastForwardModel,candidateExploratoryPositive,true);

    % processing possible negative exploratory candidates
    candidateExploratoryNegative = options.ExploratorySamplingFunction(initialPoint,stepSizeNegativeMin,stepSizeNegativeMax,normalizedDirection,options.ExploratorySamplingOptions{:});
    [candidateExploratoryNegative,scoreExploratoryNegative] = select_correct_label_candidates(fastForwardModel,candidateExploratoryNegative,false);

    % process if not enough negative samples - add to positive boundary
    availableNegativeBoundary = size(candidateBoundaryNegative,1);
    availableNegativeExploratory = size(candidateExploratoryNegative,1);
    remainingNegativeBoundary = max(nDesiredSample.NegativeBoundary - availableNegativeBoundary,0);
    remainingNegativeExploratory = max(nDesiredSample.NegativeExploratory - availableNegativeExploratory,0);

    nSelect.PositiveBoundary = nDesiredSample.PositiveBoundary + remainingNegativeBoundary + remainingNegativeExploratory;
    nSelect.PositiveExploratory = nDesiredSample.PositiveExploratory;
    nSelect.NegativeBoundary = nDesiredSample.NegativeBoundary - remainingNegativeBoundary;
    nSelect.NegativeExploratory = nDesiredSample.NegativeExploratory - remainingNegativeExploratory;

    % build array with all samples already evaluated
    alreadyEvaluated = [positiveTraining;negativeTraining];

    % select positive boundary designs
    candidateBoundaryPositive = select_appropriate_boundary_candidates(...
        candidateBoundaryPositive,scoreBoundaryPositive,nSelect.PositiveBoundary); % remove candidates on design space boundary
    iSelected = design_select_max_distance(...
    	alreadyEvaluated,candidateBoundaryPositive,nSelect.PositiveBoundary,options.knnsearchOptions{:});
    candidateBoundaryPositive = candidateBoundaryPositive(iSelected,:);
    scoreBoundaryPositive = scoreBoundaryPositive(iSelected);
    alreadyEvaluated = [alreadyEvaluated;candidateBoundaryPositive];

    % select negative boundary designs
    candidateBoundaryNegative = select_appropriate_boundary_candidates(...
        candidateBoundaryNegative,scoreBoundaryNegative,nSelect.NegativeBoundary); % remove candidates on design space boundary
    iSelected = design_select_max_distance(...
    	alreadyEvaluated,candidateBoundaryNegative,nSelect.NegativeBoundary,options.knnsearchOptions{:});
    candidateBoundaryNegative = candidateBoundaryNegative(iSelected,:);
    scoreBoundaryNegative = scoreBoundaryNegative(iSelected);
    alreadyEvaluated = [alreadyEvaluated;candidateBoundaryNegative];

    % select positive exploratory designs
    iSelected = design_select_max_distance(...
    	alreadyEvaluated,candidateExploratoryPositive,nSelect.PositiveExploratory,options.knnsearchOptions{:});
    candidateExploratoryPositive = candidateExploratoryPositive(iSelected,:);
    scoreExploratoryPositive = scoreExploratoryPositive(iSelected);
    alreadyEvaluated = [alreadyEvaluated;candidateExploratoryPositive];

    % select negative exploratory designs
    iSelected = design_select_max_distance(...
    	alreadyEvaluated,candidateExploratoryNegative,nSelect.NegativeExploratory,options.knnsearchOptions{:});
    candidateExploratoryNegative = candidateExploratoryNegative(iSelected,:);
    scoreExploratoryNegative = scoreExploratoryNegative(iSelected);

    % concatenate final result
    designSample = [...
        candidateBoundaryPositive;...
        candidateBoundaryNegative;...
    	candidateExploratoryPositive;...
        candidateExploratoryNegative];
    scorePrediction = [...
        scoreBoundaryPositive;...
        scoreBoundaryNegative;...
        scoreExploratoryPositive;...
        scoreExploratoryNegative];
    labelBoundary = [...
        true(size(candidateBoundaryPositive,1),1);...
        true(size(candidateBoundaryNegative,1),1);...
        false(size(candidateExploratoryPositive,1),1);...
        false(size(candidateExploratoryNegative,1),1)];
end

function [candidate,score] = select_correct_label_candidates(fastForwardModel,candidate,wantsPositive)
    [label,score] = design_deficit_to_label_score(fastForwardModel.predict_designs(candidate));
    candidate = candidate(label==wantsPositive,:);
    score = score(label==wantsPositive);
    [~,order] = sort(abs(score));
    candidate = candidate(order,:);
    score = score(order);
end

function candidateSelected = select_appropriate_boundary_candidates(candidate,score,nMinSelect)
    % get average score, and remove the candidates that have (absolute) score higher than that
    iSelected = (abs(score)<=mean(abs(score)));
    nSelected = sum(iSelected);
    if(nSelected<nMinSelect)
        missingSelect = nMinSelect-nSelected;
        firstScoreAboveMean = find(~iSelected,1,'first');
        iSelected(firstScoreAboveMean:firstScoreAboveMean+missingSelect-1) = true;
    end

    candidateSelected = candidate(iSelected,:);
end

function candidateBoundary = sample_line_search_boundary(initialPoint,stepSize,direction,varargin)
    candidateBoundary = initialPoint + stepSize.*direction;
end

function candidateExploratory = sample_line_search_exploratory(initialPoint,stepSizeMin,stepSizeMax,direction,varargin)
    nCandidate = size(initialPoint,1);
    stepSize = stepSizeMin + rand(nCandidate,1,varargin{:}).*(stepSizeMax-stepSizeMin);
    candidateExploratory = initialPoint + stepSize.*direction;
end