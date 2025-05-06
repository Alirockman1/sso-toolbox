function [designSample,scorePrediction,labelBoundary] = active_learning_sampling_random(fastForwardModel,designSpace,nDesiredSample,nCandidate,varargin)
%ACTIVE_LEARNING_SAMPLING_RANDOM Sampling for training
%   ACTIVE_LEARNING_SAMPLING_RANDOM uses pure random sampling to generate 
%	candidate samples to be used in the active learning process. While not 
%	enough samples of one label have been produced, the process repeats the
%	sampling process until the desired number is achieved.
%
%   DESIGNSAMPLE = ACTIVE_LEARNING_SAMPLING_RANDOM(FASTFORWARDMODEL,DESIGNSPACE,
%   NDESIREDSAMPLE,NCANDIDATE) uses the current fast-forward model 
%	FASTFORWARDMODEL and produces design samples DESIGNSAMPLE inside the design 
%	space DESIGNSPACE with distribution/properties as requested in 
%	NDESIREDSAMPLE; in the process, a number of candidates NCANDIDATE are used. 
%	NDESIREDSAMPLE is composed of four fields specifying how many samples are 
%	needed in regards to region/proximity to the boundary: 'PositiveBoundary', 
%	'NegativeBoundary', 'PositiveExploratory' and 'NegativeExploratory'.
%
%   DESIGNSAMPLE = ACTIVE_LEARNING_SAMPLING_RANDOM(...NAME,VALUE,...) also 
%   allows for the specification of options for the method in OPTIONS. This can 
%   be given as a structure or as name-value pair arguments. Said options are:
%       - 'SamplingMethodFunction' : base sampling method. Default: 
%		@sampling_random.
%       - 'SamplingMethodOptions' : options for the base sampling method. 
%		Default is empty.
%       - 'BoundingBoxGrowth' : the sampling box is created by growing the 
%		strict bounding box by this fixed rate. Default: 0.1.
%		- 'MaxCandidateSample' : maximum number of candidates that can be 
%		generated. Default: 100 times NCANDIDATE.
%
%   [DESIGNSAMPLE,SCOREPREDICTION] = 
%   ACTIVE_LEARNING_SAMPLING_RANDOM(...) also returns the 
%   predicted score SCOREPREDICTION for each sample in DESIGNSAMPLE. This may
%   allow for further analysis in terms of the error rate of the current
%   fast-forward model.
%
%   [DESIGNSAMPLE,SCOREPREDICTION,LABELBOUNDARY] = 
%   ACTIVE_LEARNING_SAMPLING_RANDOM(...) also returns the labels
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
%   See also active_learning_model_training.
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
	parser.addParameter('SamplingMethodFunction',@sampling_random,@(x)isa(x,'function_handle'));
    parser.addParameter('SamplingMethodOptions',{},@(x)iscell(x));
    parser.addParameter('BoundingBoxGrowth',0.1,@(x)isnumeric(x));
    parser.addParameter('MaxCandidateSample',100*nCandidate,@(x)isnumeric(x));
    parser.parse(varargin{:});
    options = parser.Results;

    % find box to sample from around current good region
    sampleTraining = fastForwardModel.DesignSampleTrain;
    labelTraining = fastForwardModel.LabelTrain;
    samplingBoxGrowthRate = options.BoundingBoxGrowth;
    samplingBox = design_bounding_box(sampleTraining,labelTraining);

    % initialize sampling process
    candidateSample = [];
    label = false(0,1); % empty logical array
    score = [];
    nDesiredPositive = nDesiredSample.Positive+nDesiredSample.Boundary;
    nDesiredNegative = nDesiredSample.Negative+nDesiredSample.Boundary;
    notEnoughPositive = true;
    notEnoughNegative = true;
    nGeneratedSample = 0;
    while((notEnoughPositive || notEnoughNegative) && (nGeneratedSample<=options.MaxCandidateSample))
    	if(notEnoughNegative)
    		samplingBoxGrowthRate = 2*samplingBoxGrowthRate;
    	end
    	if(notEnoughPositive)
    		samplingBoxGrowthRate = samplingBoxGrowthRate/2;
    	end

    	% grow sampling box
    	samplingBoxGrown = design_box_grow_fixed(samplingBox,designSpace(1,:),designSpace(2,:),samplingBoxGrowthRate);

    	% sample inside this box and evaluate
		newCandidateSample = options.SamplingMethodFunction(samplingBoxGrown,nCandidate,options.SamplingMethodOptions{:});
		[newLabel,newScore] = design_deficit_to_label_score(fastForwardModel.predict_designs(newCandidateSample));

		% concatenate full information
		candidateSample = [candidateSample;newCandidateSample];
		label = [label;newLabel];
		score = [score;newScore];

		% check convergence
		availablePositive = sum(label);
		notEnoughPositive = (availablePositive<nDesiredPositive);

		availableNegative = sum(~label);
		notEnoughNegative = (availableNegative<nDesiredNegative);

		nGeneratedSample = size(candidateSample,1);
    end

    % process if not enough samples 
    nAvailablePositive = sum(label);
    nAvailableNegative = sum(~label);
    nSelectPositive = min(nDesiredSample.Positive,nAvailablePositive);
    nSelectNegative = min(nDesiredSample.Negative,nAvailableNegative);

    % add points not available to the boundary selection
    remainingPositive = max(nDesiredSample.Positive - nAvailablePositive,0);
    remainingNegative = max(nDesiredSample.Negative - nAvailableNegative,0);
    nSelectBoundary = nDesiredSample.Boundary + remainingPositive + remainingNegative;

    % select boundary candidates
    [candidateBoundary,scoreBoundary,iSelectBoundary] = select_boundary_designs(candidateSample,score,nSelectBoundary);

    % remove boundary from available selection
    candidateSample(iSelectBoundary,:) = [];
    score(iSelectBoundary) = [];
    label(iSelectBoundary) = [];

    % select positive candidates
	[candidateExploratoryPositive,scoreExploratoryPositive] = ...
		select_exploratory_designs(candidateSample,score,label,nSelectPositive);

	% select negative candidates
	[candidateExploratoryNegative,scoreExploratoryNegative] = ...
		select_exploratory_designs(candidateSample,score,~label,nSelectNegative);

	% concatenate final result
    designSample = [...
        candidateBoundary;...
    	candidateExploratoryPositive;...
        candidateExploratoryNegative];
    scorePrediction = [...
        scoreBoundary;...
        scoreExploratoryPositive;...
        scoreExploratoryNegative];
    labelBoundary = [...
        true(size(candidateBoundary,1),1);...
        false(size(candidateExploratoryPositive,1),1);...
        false(size(candidateExploratoryNegative,1),1)];
end

function [candidateBoundary,scoreBoundary,iSelect] = select_boundary_designs(candidateSample,score,nSelect)
	[~,order] = sort(abs(score));
	score = score(order);
	candidateSample = candidateSample(order,:);

	candidateBoundary = candidateSample(1:nSelect,:);
	scoreBoundary = score(1:nSelect);
	iSelect = order(1:nSelect);
end

function [candidateExploratory,scoreExploratory] = select_exploratory_designs(candidateSample,score,label,nSelect)
    availableSample = candidateSample(label,:);
    availableScore = score(label);

	permutationAvailable = randperm(size(availableSample,1));
    iSelect = permutationAvailable(1:nSelect);
	candidateExploratory = availableSample(iSelect,:);
	scoreExploratory = availableScore(iSelect);
end